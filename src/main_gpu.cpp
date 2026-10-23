// Standalone test of AMReX FabArray/MultiFab operations on the GPU. It is the
// first step towards running FLEKS on GPUs and does not use any FLEKS class.
//
// Build: ./Config.pl -install -compiler=nvfortran,nvcc; make GPUTEST
// Run:   mpirun -np 2 bin/main_gpu.exe [n_cell=128] [max_grid_size=32]
//                                      [n_iter=100]
//
// Every check compares the GPU result with an exact or host-computed value and
// prints PASS/FAIL. The exit code is non-zero if any check fails.

#include <cmath>
#include <string>

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace {

constexpr Real twoPi = 2 * 3.14159265358979323846;

// Periodic test function on the unit cube, evaluated at the cell center of
// (i,j,k). It is also valid for ghost cells outside the domain.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
exact_value(int i, int j, int k, GpuArray<Real, 3> const& dx) noexcept {
  const Real x = (i + 0.5) * dx[0];
  const Real y = (j + 0.5) * dx[1];
  const Real z = (k + 0.5) * dx[2];
  return 1 + std::sin(twoPi * x) * std::cos(twoPi * y) +
         0.5 * std::sin(twoPi * z);
}

// Exact 7-point discrete Laplacian of exact_value() for a uniform grid.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
exact_laplacian(int i, int j, int k, GpuArray<Real, 3> const& dx) noexcept {
  const Real x = (i + 0.5) * dx[0];
  const Real y = (j + 0.5) * dx[1];
  const Real z = (k + 0.5) * dx[2];
  const Real dx2 = dx[0] * dx[0];
  const Real factor = (2 * std::cos(twoPi * dx[0]) - 2) / dx2;
  return 2 * factor * std::sin(twoPi * x) * std::cos(twoPi * y) +
         0.5 * factor * std::sin(twoPi * z);
}

int nFail = 0;

void report(const std::string& name, Real err, Real tol) {
  const bool isPass = err <= tol;
  if (!isPass)
    ++nFail;
  Print() << (isPass ? "  PASS  " : "  FAIL  ") << name << ": error = " << err
          << " (tol = " << tol << ")\n";
}

// Global max |mf - expected(i,j,k)| over the valid cells and nGhost ghosts.
template <typename F>
Real max_error(const MultiFab& mf, int nGhost, F const& expected) {
  auto const& arrs = mf.const_arrays();
  Real err = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, mf,
                       IntVect(nGhost),
                       [=] AMREX_GPU_DEVICE(int iBox, int i, int j, int k)
                           -> GpuTuple<Real> {
                         return { std::abs(arrs[iBox](i, j, k) -
                                           expected(i, j, k)) };
                       });
  ParallelDescriptor::ReduceRealMax(err);
  return err;
}

void apply_laplacian(MultiFab& lap, const MultiFab& mf,
                     GpuArray<Real, 3> const& dx) {
  auto const& dst = lap.arrays();
  auto const& src = mf.const_arrays();
  const Real dx2 = dx[0] * dx[0];
  ParallelFor(lap, [=] AMREX_GPU_DEVICE(int iBox, int i, int j, int k) {
    auto const& a = src[iBox];
    dst[iBox](i, j, k) = (a(i - 1, j, k) + a(i + 1, j, k) + a(i, j - 1, k) +
                          a(i, j + 1, k) + a(i, j, k - 1) + a(i, j, k + 1) -
                          6 * a(i, j, k)) /
                         dx2;
  });
}

void print_device_info() {
  Print() << "AMReX " << Version() << ", MPI ranks = "
          << ParallelDescriptor::NProcs() << "\n";
#ifdef AMREX_USE_GPU
  AllPrint() << "  rank " << ParallelDescriptor::MyProc() << " -> GPU "
             << Gpu::Device::deviceId() << " (" << Gpu::Device::deviceName()
             << ", " << Gpu::Device::totalGlobalMem() / (1024 * 1024)
             << " MiB)\n";
  Print() << "  The_Arena is managed memory: "
          << (The_Arena()->isManaged() ? "yes" : "no") << "\n";
#else
  Print() << "  AMReX is compiled without GPU support\n";
#endif
}

void run_tests() {
  int nCell = 128;
  int maxGridSize = 32;
  int nIter = 100;
  {
    ParmParse pp;
    pp.query("n_cell", nCell);
    pp.query("max_grid_size", maxGridSize);
    pp.query("n_iter", nIter);
  }

  print_device_info();

  // Periodic unit cube
  const Box domain(IntVect(0), IntVect(nCell - 1));
  const RealBox realBox({ 0., 0., 0. }, { 1., 1., 1. });
  const Array<int, 3> isPeriodic{ 1, 1, 1 };
  const Geometry geom(domain, realBox, CoordSys::cartesian, isPeriodic);
  const auto dx = geom.CellSizeArray();

  BoxArray ba(domain);
  ba.maxSize(maxGridSize);
  const DistributionMapping dm(ba);
  Print() << "Domain " << nCell << "^3, " << ba.size() << " boxes of "
          << maxGridSize << "^3\n\n";

  const int nGhost = 1;
  MultiFab mf(ba, dm, 1, nGhost);
  const auto exact = [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
    return exact_value(i, j, k, dx);
  };

  // 1. Fill the valid cells on the device
  {
    auto const& arrs = mf.arrays();
    ParallelFor(mf, [=] AMREX_GPU_DEVICE(int iBox, int i, int j, int k) {
      arrs[iBox](i, j, k) = exact(i, j, k);
    });
    Gpu::streamSynchronize();
    report("ParallelFor fill", max_error(mf, 0, exact), 1e-14);
  }

  // 2. Built-in reductions compared with a brute-force host reference
  {
    Real sumRef = 0, maxRef = 0, sum2Ref = 0, minRef = 1e30;
    for (int k = 0; k < nCell; ++k)
      for (int j = 0; j < nCell; ++j)
        for (int i = 0; i < nCell; ++i) {
          const Real v = exact_value(i, j, k, dx);
          sumRef += v;
          sum2Ref += v * v;
          maxRef = std::max(maxRef, std::abs(v));
          minRef = std::min(minRef, v);
        }
    // The summation order differs from the host loop, so the sums agree only
    // to round-off accumulated over all cells.
    report("MultiFab::sum", std::abs(mf.sum(0) - sumRef) / sumRef, 1e-12);
    report("MultiFab::norm0", std::abs(mf.norm0(0) - maxRef), 1e-14);
    report("MultiFab::norm2", std::abs(mf.norm2(0) - std::sqrt(sum2Ref)),
           1e-10 * std::sqrt(sum2Ref));

    // Custom fused reduction: sum of squares and minimum in one kernel
    auto const& arrs = mf.const_arrays();
    ReduceOps<ReduceOpSum, ReduceOpMin> reduceOp;
    ReduceData<Real, Real> reduceData(reduceOp);
    reduceOp.eval(mf, IntVect(0), reduceData,
                  [=] AMREX_GPU_DEVICE(int iBox, int i, int j, int k)
                      -> GpuTuple<Real, Real> {
                    const Real v = arrs[iBox](i, j, k);
                    return { v * v, v };
                  });
    auto const result = reduceData.value(reduceOp);
    Real sum2 = get<0>(result);
    Real minVal = get<1>(result);
    ParallelDescriptor::ReduceRealSum(sum2);
    ParallelDescriptor::ReduceRealMin(minVal);
    report("ReduceOps sum(v^2)", std::abs(sum2 - sum2Ref) / sum2Ref, 1e-12);
    report("ReduceOps min(v)", std::abs(minVal - minRef), 1e-14);
  }

  // 3. Periodic ghost cell exchange across boxes and MPI ranks
  {
    mf.setBndry(-999.);
    mf.FillBoundary(geom.periodicity());
    report("FillBoundary (periodic)", max_error(mf, nGhost, exact), 1e-13);
  }

  // 4. MultiFab arithmetic
  {
    MultiFab mf2(ba, dm, 1, 0);
    mf2.setVal(3.);
    MultiFab::Saxpy(mf2, 2., mf, 0, 0, 1, 0); // mf2 = 2*mf + 3
    report("MultiFab::Saxpy",
           max_error(mf2, 0,
                     [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                       return 2 * exact(i, j, k) + 3;
                     }),
           1e-13);

    MultiFab mf3(ba, dm, 1, 0);
    MultiFab::LinComb(mf3, 0.5, mf, 0, -1., mf2, 0, 0, 1, 0);
    report("MultiFab::LinComb",
           max_error(mf3, 0,
                     [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                       return -1.5 * exact(i, j, k) - 3;
                     }),
           1e-13);
  }

  // 5. Stencil kernel that reads the ghost cells
  MultiFab lap(ba, dm, 1, 0);
  {
    apply_laplacian(lap, mf, dx);
    report("7-point Laplacian",
           max_error(lap, 0,
                     [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                       return exact_laplacian(i, j, k, dx);
                     }),
           1e-7);
  }

  // 6. Copy to a different BoxArray and DistributionMapping
  {
    BoxArray baNew(domain);
    baNew.maxSize(maxGridSize / 2);
    const DistributionMapping dmNew(baNew);
    MultiFab mfNew(baNew, dmNew, 1, nGhost);
    mfNew.setVal(-999.);
    mfNew.ParallelCopy(mf, 0, 0, 1, 0, nGhost, geom.periodicity());
    report("ParallelCopy to " + std::to_string(baNew.size()) + " boxes",
           max_error(mfNew, nGhost, exact), 1e-13);
  }

  // 7. Copy to host memory and check on the CPU, independent of any GPU
  //    reduction
  {
    Real err = 0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const auto& fab = mf[mfi];
      FArrayBox hostFab(fab.box(), 1, The_Pinned_Arena());
      Gpu::dtoh_memcpy_async(hostFab.dataPtr(), fab.dataPtr(), fab.nBytes());
      Gpu::streamSynchronize();
      auto const& a = hostFab.const_array();
      amrex::LoopOnCpu(fab.box(), [&](int i, int j, int k) {
        err = std::max(err, std::abs(a(i, j, k) - exact_value(i, j, k, dx)));
      });
    }
    ParallelDescriptor::ReduceRealMax(err);
    report("Device to host copy", err, 1e-13);
  }

  // 8. Timing
  {
    const double nBytes = 2. * sizeof(Real) * nCell * nCell * nCell;

    Gpu::synchronize();
    ParallelDescriptor::Barrier();
    double t0 = second();
    for (int iter = 0; iter < nIter; ++iter)
      apply_laplacian(lap, mf, dx);
    Gpu::synchronize();
    ParallelDescriptor::Barrier();
    const double tLap = (second() - t0) / nIter;

    t0 = second();
    for (int iter = 0; iter < nIter; ++iter)
      mf.FillBoundary(geom.periodicity());
    Gpu::synchronize();
    ParallelDescriptor::Barrier();
    const double tFill = (second() - t0) / nIter;

    Print() << "\nTiming (average of " << nIter << " iterations)\n"
            << "  7-point Laplacian: " << tLap * 1e3 << " ms, "
            << nBytes / tLap / 1e9 << " GB/s effective\n"
            << "  FillBoundary:      " << tFill * 1e3 << " ms\n";
  }
}

// Stage MPI messages through host memory unless amrex.use_gpu_aware_mpi is
// given. The shared-memory transport (OMPI_MCA_btl=self,vader) cannot read
// device buffers and crashes when AMReX detects a CUDA-aware MPI.
void set_default_parameters() {
  ParmParse pp("amrex");
  if (!pp.contains("use_gpu_aware_mpi"))
    pp.add("use_gpu_aware_mpi", 0);
}

} // namespace

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, set_default_parameters);
  run_tests();
  const int nFailTotal = nFail;
  if (nFailTotal == 0)
    Print() << "\nALL AMREX GPU TESTS PASSED\n";
  else
    Print() << "\n" << nFailTotal << " AMREX GPU TEST(S) FAILED\n";
  amrex::Finalize();
  return nFailTotal == 0 ? 0 : 1;
}
