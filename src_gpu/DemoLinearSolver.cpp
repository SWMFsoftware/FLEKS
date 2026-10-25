#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX.H>
#include <AMReX_GpuReduce.H>
#include <AMReX_ParallelDescriptor.H>

#include "DemoLinearSolver.h"
using namespace amrex;

// On the GPU the batched solver avoids host waits: AMReX vector operations,
// MFIter loops, and local FillBoundary copies all synchronize the stream
// unless they run in a Gpu::NoSyncRegion, and every MultiFab::Dot returns to
// the host. Inside NoSyncRegion this file uses only single-launch
// whole-MultiFab kernels on the current stream, never MFIter loops, which
// could otherwise spread boxes over several streams without a closing sync.
// The caller's operator runs outside these regions and keeps its own
// synchronization.

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
    // project() keeps one register-array partial per projection plus norm.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(restart + 2 <= 64,
                                     "gmres_batched supports gmres_restart<=62");
    packedBasis = std::make_unique<MultiFab>(ba, dm, nc * (restart + 1), 1);
    for (int n = 0; n <= restart; ++n)
      basis.push_back(std::make_unique<MultiFab>(*packedBasis, amrex::make_alias,
                                                 n * nc, nc));
    // Two passes, each with up to restart+1 projections and one norm.
    deviceProjection.resize(2 * (restart + 2));
    hostProjection.resize(2 * (restart + 2));
    deviceCoefficients.resize(restart + 1);
    hostCoefficients.resize(restart + 1);
    std::vector<BoxRange> ranges;
    for (int local = 0; local < work.local_size(); ++local) {
      const Box box = ba[work.IndexArray()[local]];
      BoxRange range;
      for (int d = 0; d < 3; ++d) {
        range.lo[d] = box.smallEnd(d);
        range.length[d] = box.length(d);
      }
      range.offset = localCells;
      ranges.push_back(range);
      localCells += box.numPts();
    }
    boxRanges.resize(ranges.size());
    Gpu::copyAsync(Gpu::hostToDevice, ranges.begin(), ranges.end(),
                   boxRanges.begin());
    Gpu::streamSynchronize();
  } else {
    for (int n = 0; n <= restart; ++n)
      basis.push_back(std::make_unique<MultiFab>(ba, dm, nc, 1));
  }
}

Real DemoLinearSolver::dot(const MultiFab& a, const MultiFab& b) const {
  return MultiFab::Dot(*owner, a, 0, b, 0, nComp, 0);
}

Real DemoLinearSolver::waited_dot(const MultiFab& a, const MultiFab& b) {
  ++statistics.hostWaits;
  return dot(a, b);
}

template <class F>
void DemoLinearSolver::section(double DemoSolverStats::*field, F&& operation) {
  if (!timing) {
    operation();
    return;
  }
  Gpu::synchronize();
  const double start = ParallelDescriptor::second();
  operation();
  Gpu::synchronize();
  statistics.*field += ParallelDescriptor::second() - start;
}

void DemoLinearSolver::project(int count, int pass, bool withNorm) {
  Real* projection = deviceProjection.data() + pass * (restartLength + 2);
  const int nValue = count + (withNorm ? 1 : 0);
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
  if (Gpu::inLaunchRegion()) {
    // One launch over every local cell of every box. Each block reduces its
    // cells for every requested vector, then adds one atomic partial per
    // value; the projections stay on the device for subtract_projection.
    constexpr int nthread = 256;
    const int nc = nComp;
    const Long cells = localCells;
    const int nbox = static_cast<int>(boxRanges.size());
    const BoxRange* ranges = boxRanges.data();
    const auto w = work.const_arrays();
    const auto vectors = packedBasis->const_arrays();
    const auto mask = owner->const_arrays();
    ParallelFor(nValue,
                [=] AMREX_GPU_DEVICE(int n) noexcept { projection[n] = 0; });
    const int nblocks = static_cast<int>(
        std::min<Long>((cells + nthread - 1) / nthread, 1024));
    amrex::launch<nthread>(nblocks, Gpu::gpuStream(), [=] AMREX_GPU_DEVICE() {
      Real partial[64];
      for (int n = 0; n < nValue; ++n)
        partial[n] = 0;
      for (Long cell = Long(blockIdx.x) * nthread + threadIdx.x; cell < cells;
           cell += Long(gridDim.x) * nthread) {
        int box = 0;
        while (box + 1 < nbox && ranges[box + 1].offset <= cell)
          ++box;
        const BoxRange& range = ranges[box];
        const Long local = cell - range.offset;
        const int i = range.lo[0] + local % range.length[0];
        const int j = range.lo[1] + (local / range.length[0]) % range.length[1];
        const int k =
            range.lo[2] + local / (Long(range.length[0]) * range.length[1]);
        if (!mask[box](i, j, k))
          continue;
        for (int c = 0; c < nc; ++c) {
          const Real value = w[box](i, j, k, c);
          for (int n = 0; n < count; ++n)
            partial[n] += value * vectors[box](i, j, k, n * nc + c);
          if (withNorm)
            partial[count] += value * value;
        }
      }
      for (int n = 0; n < nValue; ++n) {
        const Real sum = Gpu::blockReduceSum<nthread>(partial[n]);
        if (threadIdx.x == 0)
          Gpu::Atomic::AddNoRet(projection + n, sum);
      }
    });
  } else
#endif
  {
    // Host path: deviceProjection is host memory in a CPU build.
    for (int n = 0; n < count; ++n)
      projection[n] =
          MultiFab::Dot(*owner, work, 0, *basis[n], 0, nComp, 0, true);
    if (withNorm)
      projection[count] = MultiFab::Dot(*owner, work, 0, work, 0, nComp, 0, true);
  }
  if (ParallelDescriptor::NProcs() > 1) {
    // The vector update needs global projections: one wait and one MPI
    // reduction per pass.
    Real* host = hostProjection.data() + pass * (restartLength + 2);
    Gpu::copyAsync(Gpu::deviceToHost, projection, projection + nValue, host);
    Gpu::streamSynchronize();
    ++statistics.hostWaits;
    ParallelDescriptor::ReduceRealSum(host, nValue);
    Gpu::copyAsync(Gpu::hostToDevice, host, host + nValue, projection);
  }
}

void DemoLinearSolver::subtract_projection(int count, int pass) {
  const Real* projection = deviceProjection.data() + pass * (restartLength + 2);
  const int nc = nComp;
  const auto w = work.arrays();
  const auto vectors = packedBasis->const_arrays();
  ParallelFor(work, IntVect(0), nc,
              [=] AMREX_GPU_DEVICE(int box, int i, int j, int k, int c) noexcept {
                Real value = w[box](i, j, k, c);
                for (int n = 0; n < count; ++n)
                  value -= projection[n] * vectors[box](i, j, k, n * nc + c);
                w[box](i, j, k, c) = value;
              });
}

void DemoLinearSolver::scale_into(MultiFab& destination, const MultiFab& source,
                                  Real scale) {
  const auto dst = destination.arrays();
  const auto src = source.const_arrays();
  ParallelFor(destination, IntVect(0), nComp,
              [=] AMREX_GPU_DEVICE(int box, int i, int j, int k, int c) noexcept {
                dst[box](i, j, k, c) = scale * src[box](i, j, k, c);
              });
}

void DemoLinearSolver::residual_from(const MultiFab& rhs) {
  const auto r = residual.arrays();
  const auto b = rhs.const_arrays();
  const auto ax = work.const_arrays();
  ParallelFor(residual, IntVect(0), nComp,
              [=] AMREX_GPU_DEVICE(int box, int i, int j, int k, int c) noexcept {
                r[box](i, j, k, c) = b[box](i, j, k, c) - ax[box](i, j, k, c);
              });
}

void DemoLinearSolver::add_solution(MultiFab& x, const std::vector<Real>& y,
                                    int used) {
  // The pinned staging buffer is reused only after the next host wait (the
  // true-residual dot product), by which time this copy has completed.
  std::copy_n(y.begin(), used, hostCoefficients.begin());
  Gpu::copyAsync(Gpu::hostToDevice, hostCoefficients.begin(),
                 hostCoefficients.begin() + used, deviceCoefficients.begin());
  const Real* coefficient = deviceCoefficients.data();
  const int nc = nComp;
  const auto dst = x.arrays();
  const auto vectors = packedBasis->const_arrays();
  ParallelFor(x, IntVect(0), nc,
              [=] AMREX_GPU_DEVICE(int box, int i, int j, int k, int c) noexcept {
                Real value = dst[box](i, j, k, c);
                for (int n = 0; n < used; ++n)
                  value += coefficient[n] * vectors[box](i, j, k, n * nc + c);
                dst[box](i, j, k, c) = value;
              });
}

// Twice classical Gram-Schmidt for Arnoldi column j. Both passes run on the
// device back to back; on one rank the host waits once, after the second
// pass, for the Hessenberg column and the norm. The second pass also
// reduces w.w for the vector it projects, so the norm of the final vector is
// sqrt(w.w - sum p_i^2) for the orthonormal basis. If that difference loses
// too many digits, the norm is recomputed with an explicit dot product.
void DemoLinearSolver::solve_batched_column(int j, std::vector<Real>& column,
                                            Real& next) {
  const int count = j + 1, stride = restartLength + 2;
  {
    Gpu::NoSyncRegion noSync;
    project(count, 0, false);
    subtract_projection(count, 0);
    project(count, 1, true);
    subtract_projection(count, 1);
    Gpu::copyAsync(Gpu::deviceToHost, deviceProjection.begin(),
                   deviceProjection.begin() + 2 * stride,
                   hostProjection.begin());
  }
  Gpu::streamSynchronize();
  ++statistics.hostWaits;
  const Real* first = hostProjection.data();
  const Real* second = hostProjection.data() + stride;
  Real removed = 0;
  for (int i = 0; i < count; ++i) {
    column[i] = first[i] + second[i];
    removed += second[i] * second[i];
  }
  const Real before = second[count];
  const Real after = before - removed;
  if (after > Real(1.e-2) * before) {
    next = std::sqrt(after);
  } else {
    next = std::sqrt(waited_dot(work, work));
    ++statistics.normFallbacks;
  }
}

DemoSolveResult DemoLinearSolver::solve(MultiFab& x, const MultiFab& rhs,
                                        const Operator& apply, Real tolerance,
                                        int maxIter, bool zeroInitialGuess) {
  DemoSolveResult result;
  ++statistics.solves;
  const Real rhsNorm = std::sqrt(waited_dot(rhs, rhs));
  const Real target = tolerance * std::max(rhsNorm, Real(1e-14));
  auto true_residual = [&]() {
    Real norm = 0;
    section(&DemoSolverStats::residualSeconds, [&] {
      apply(x, work);
      if (batched) {
        Gpu::NoSyncRegion noSync;
        residual_from(rhs);
      } else {
        MultiFab::LinComb(residual, 1, rhs, 0, -1, work, 0, 0, nComp, 0);
        ++statistics.hostWaits;
      }
      norm = std::sqrt(waited_dot(residual, residual));
    });
    return norm;
  };
  if (zeroInitialGuess) {
    // x == 0: the residual is rhs and its norm is already known.
    if (batched) {
      Gpu::NoSyncRegion noSync;
      scale_into(residual, rhs, 1);
    } else {
      MultiFab::Copy(residual, rhs, 0, 0, nComp, 0);
      ++statistics.hostWaits;
    }
    result.initialResidual = rhsNorm;
  } else {
    result.initialResidual = true_residual();
  }
  result.residual = result.initialResidual;
  while (result.iterations < maxIter && result.residual > target) {
    const int m = std::min(restartLength, maxIter - result.iterations);
    std::vector<Real> h((m + 1) * m, 0), cs(m), sn(m), g(m + 1, 0), y(m, 0);
    auto H = [&](int i, int j) -> Real& { return h[j * (m + 1) + i]; };
    if (batched) {
      Gpu::NoSyncRegion noSync;
      scale_into(*basis[0], residual, 1 / result.residual);
    } else {
      MultiFab::LinComb(*basis[0], 1 / result.residual, residual, 0, 0,
                        residual, 0, 0, nComp, 0);
      ++statistics.hostWaits;
    }
    g[0] = result.residual;
    int used = 0;
    std::vector<Real> column(m + 1, 0);
    for (int j = 0; j < m; ++j) {
      section(&DemoSolverStats::operatorSeconds,
              [&] { apply(*basis[j], work); });
      Real next = 0;
      section(&DemoSolverStats::orthogonalizeSeconds, [&] {
        if (batched) {
          solve_batched_column(j, column, next);
          for (int i = 0; i <= j; ++i)
            H(i, j) = column[i];
        } else {
          // Retain twice modified Gram-Schmidt as the reference algorithm.
          for (int pass = 0; pass < 2; ++pass)
            for (int i = 0; i <= j; ++i) {
              const Real projection = waited_dot(work, *basis[i]);
              H(i, j) += projection;
              MultiFab::Saxpy(work, -projection, *basis[i], 0, 0, nComp, 0);
              ++statistics.hostWaits;
            }
          next = std::sqrt(waited_dot(work, work));
        }
      });
      H(j + 1, j) = next;
      section(&DemoSolverStats::updateSeconds, [&] {
        if (next > 0) {
          if (batched) {
            Gpu::NoSyncRegion noSync;
            scale_into(*basis[j + 1], work, 1 / next);
          } else {
            MultiFab::LinComb(*basis[j + 1], 1 / next, work, 0, 0, work, 0, 0,
                              nComp, 0);
            ++statistics.hostWaits;
          }
        }
      });
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
      ++statistics.iterations;
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
    section(&DemoSolverStats::updateSeconds, [&] {
      if (batched) {
        Gpu::NoSyncRegion noSync;
        add_solution(x, y, used);
      } else {
        for (int j = 0; j < used; ++j) {
          MultiFab::Saxpy(x, y[j], *basis[j], 0, 0, nComp, 0);
          ++statistics.hostWaits;
        }
      }
    });
    result.residual = true_residual();
    if (!std::isfinite(result.residual))
      break;
  }
  result.converged =
      std::isfinite(result.residual) && result.residual <= target;
  return result;
}
