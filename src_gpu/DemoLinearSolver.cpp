#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX_GpuAsyncArray.H>
#include <AMReX_GpuReduce.H>

#include "DemoLinearSolver.h"
using namespace amrex;

DemoLinearSolver::DemoLinearSolver(const BoxArray& ba,
                                   const DistributionMapping& dm, int nc,
                                   int restart, const Geometry& geom,
                                   bool batchedIn)
    : nComp(nc),
      restartLength(restart),
      batched(batchedIn),
      residual(ba, dm, nc, 1),
      work(ba, dm, nc, 1) {
  owner = residual.OwnerMask(geom.periodicity());
  if (batched) {
    packedBasis = std::make_unique<MultiFab>(ba, dm, nc * (restart + 1), 0);
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    deviceProjection.resize(restart);
    hostProjection.resize(restart);
#endif
  }
  for (int n = 0; n <= restart; ++n)
    basis.push_back(std::make_unique<MultiFab>(ba, dm, nc, 1));
}
Real DemoLinearSolver::dot(const MultiFab& a, const MultiFab& b) const {
  return MultiFab::Dot(*owner, a, 0, b, 0, nComp, 0);
}
// Coefficients are copied once; each output element traverses the packed basis
// in one kernel instead of launching a SAXPY for every Arnoldi vector.
namespace {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
// One scalar per Arnoldi vector stays on device across all boxes. Every
// block reduces a contiguous cell chunk, then contributes one atomic partial.
// All launches and the sole device-to-host copy use stream zero, avoiding a
// host-returning reduction and synchronization for each basis vector.
void project_basis(const MultiFab& work, const MultiFab& packedBasis,
                   const iMultiFab& owner, int nc, int count,
                   Gpu::DeviceVector<Real>& deviceProjection,
                   Gpu::PinnedVector<Real>& hostProjection) {
  constexpr int nthread = 256;
  Real* projection = deviceProjection.data();
  Gpu::Device::setStreamIndex(0);
  ParallelFor(count,
              [=] AMREX_GPU_DEVICE(int n) noexcept { projection[n] = 0; });
  for (MFIter mfi(work, MFItInfo().SetNumStreams(1).DisableDeviceSync());
       mfi.isValid(); ++mfi) {
    const auto box = mfi.validbox();
    const auto small = box.smallEnd();
    const auto length = box.length();
    const Long cells = box.numPts();
    const Long chunks = (cells + nthread - 1) / nthread;
    AMREX_ALWAYS_ASSERT(chunks <= std::numeric_limits<int>::max() / count);
    const auto w = work.const_array(mfi);
    const auto vectors = packedBasis.const_array(mfi);
    const auto mask = owner.const_array(mfi);
    amrex::launch<nthread>(
        static_cast<int>(chunks * count), Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE() {
          const int n = blockIdx.x / chunks;
          const Long cell = (blockIdx.x % chunks) * nthread + threadIdx.x;
          Real value = 0;
          if (cell < cells) {
            const int i = small[0] + cell % length[0];
            const int j = small[1] + (cell / length[0]) % length[1];
            const int k = small[2] + cell / (Long(length[0]) * length[1]);
            if (mask(i, j, k))
              for (int c = 0; c < nc; ++c)
                value += w(i, j, k, c) * vectors(i, j, k, n * nc + c);
          }
          const Real sum = Gpu::blockReduceSum<nthread>(value);
          if (threadIdx.x == 0)
            Gpu::Atomic::AddNoRet(projection + n, sum);
        });
  }
  Gpu::copyAsync(Gpu::deviceToHost, deviceProjection.begin(),
                 deviceProjection.begin() + count, hostProjection.begin());
  Gpu::streamSynchronize();
}
#endif

void add_basis(MultiFab& destination, const MultiFab& packedBasis, int nc,
               const std::vector<Real>& coefficients, int count) {
  for (MFIter mfi(destination); mfi.isValid(); ++mfi) {
    // MFIter can select a different stream for each box. Copy and release
    // coefficients on that same stream as the consumer kernel.
    Gpu::AsyncArray<Real> deviceCoefficients(coefficients.data(), count);
    const Real* coeff = deviceCoefficients.data();
    const auto dst = destination.array(mfi);
    const auto vectors = packedBasis.const_array(mfi);
    ParallelFor(mfi.validbox(), nc,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  Real value = dst(i, j, k, c);
                  for (int n = 0; n < count; ++n)
                    value += coeff[n] * vectors(i, j, k, n * nc + c);
                  dst(i, j, k, c) = value;
                });
  }
}
} // namespace
DemoSolveResult DemoLinearSolver::solve(MultiFab& x, const MultiFab& rhs,
                                        const Operator& apply, Real tolerance,
                                        int maxIter) {
  DemoSolveResult result;
  const Real rhsNorm = std::sqrt(dot(rhs, rhs));
  const Real target = tolerance * std::max(rhsNorm, Real(1e-14));
  auto true_residual = [&]() {
    apply(x, work);
    MultiFab::LinComb(residual, 1, rhs, 0, -1, work, 0, 0, nComp, 0);
    return std::sqrt(dot(residual, residual));
  };
  result.initialResidual = true_residual();
  result.residual = result.initialResidual;
  while (result.iterations < maxIter && result.residual > target) {
    const int m = std::min(restartLength, maxIter - result.iterations);
    std::vector<Real> h((m + 1) * m, 0), cs(m), sn(m), g(m + 1, 0), y(m, 0);
    auto H = [&](int i, int j) -> Real& { return h[j * (m + 1) + i]; };
    MultiFab::LinComb(*basis[0], 1 / result.residual, residual, 0, 0, residual,
                      0, 0, nComp, 0);
    if (batched)
      MultiFab::Copy(*packedBasis, *basis[0], 0, 0, nComp, 0);
    g[0] = result.residual;
    int used = 0;
    for (int j = 0; j < m; ++j) {
      apply(*basis[j], work);
      if (batched) {
        // Twice classical Gram-Schmidt: all projections in a pass use the
        // same work vector, then one vector update applies them together.
        // GPU projections share one device reduction buffer and one transfer
        // per pass; MPI also reduces the entire projection vector at once.
        std::vector<Real> projection(j + 1);
        for (int pass = 0; pass < 2; ++pass) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
          if (Gpu::inLaunchRegion()) {
            project_basis(work, *packedBasis, *owner, nComp, j + 1,
                          deviceProjection, hostProjection);
            std::copy_n(hostProjection.begin(), j + 1, projection.begin());
          } else
#endif
          {
            for (int i = 0; i <= j; ++i)
              projection[i] =
                  MultiFab::Dot(*owner, work, 0, *basis[i], 0, nComp, 0, true);
          }
          ParallelDescriptor::ReduceRealSum(projection.data(), j + 1);
          for (int i = 0; i <= j; ++i) {
            H(i, j) += projection[i];
            projection[i] = -projection[i];
          }
          add_basis(work, *packedBasis, nComp, projection, j + 1);
        }
      } else {
        // Retain twice modified Gram-Schmidt as the reference algorithm.
        for (int pass = 0; pass < 2; ++pass)
          for (int i = 0; i <= j; ++i) {
            const Real projection = dot(work, *basis[i]);
            H(i, j) += projection;
            MultiFab::Saxpy(work, -projection, *basis[i], 0, 0, nComp, 0);
          }
      }
      H(j + 1, j) = std::sqrt(dot(work, work));
      const Real next = H(j + 1, j);
      if (next > 0) {
        MultiFab::LinComb(*basis[j + 1], 1 / next, work, 0, 0, work, 0, 0,
                          nComp, 0);
        if (batched)
          MultiFab::Copy(*packedBasis, *basis[j + 1], 0, (j + 1) * nComp, nComp,
                         0);
      }
      for (int i = 0; i < j; ++i) {
        const Real t = cs[i] * H(i, j) + sn[i] * H(i + 1, j);
        H(i + 1, j) = -sn[i] * H(i, j) + cs[i] * H(i + 1, j);
        H(i, j) = t;
      }
      const Real norm = std::hypot(H(j, j), H(j + 1, j));
      cs[j] = norm > 0 ? H(j, j) / norm : 1;
      sn[j] = norm > 0 ? H(j + 1, j) / norm : 0;
      H(j, j) = norm;
      H(j + 1, j) = 0;
      g[j + 1] = -sn[j] * g[j];
      g[j] *= cs[j];
      used = j + 1;
      ++result.iterations;
      if (std::abs(g[j + 1]) <= target ||
          next <= std::numeric_limits<Real>::epsilon() * norm)
        break;
    }
    bool singular = false;
    for (int i = used - 1; i >= 0; --i) {
      Real r = g[i];
      for (int j = i + 1; j < used; ++j)
        r -= H(i, j) * y[j];
      if (std::abs(H(i, i)) <= std::numeric_limits<Real>::min()) {
        singular = true;
        break;
      }
      y[i] = r / H(i, i);
    }
    if (singular)
      break;
    if (batched)
      add_basis(x, *packedBasis, nComp, y, used);
    else
      for (int j = 0; j < used; ++j)
        MultiFab::Saxpy(x, y[j], *basis[j], 0, 0, nComp, 0);
    result.residual = true_residual();
    if (!std::isfinite(result.residual))
      break;
  }
  result.converged =
      std::isfinite(result.residual) && result.residual <= target;
  return result;
}
