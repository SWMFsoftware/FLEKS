#include <AMReX_Print.H>

#include "DemoChecks.h"

using namespace amrex;

//==========================================================
DemoChecks::DemoChecks(const DemoGrid& grid, const DemoConfig& config)
    : nCell(config.nCell),
      maxGridSize(config.maxGridSize),
      nIter(config.nIter),
      domain(grid.geom.Domain()),
      geom(grid.geom),
      ba(grid.ba),
      dm(grid.dm) {
  mf.define(ba, dm, 1, nGhost);
  lap.define(ba, dm, 1, 0);
}

//==========================================================
int DemoChecks::run() {
  print_device_info();
  Print() << "Domain " << nCell << "^3, " << ba.size() << " boxes of "
          << maxGridSize << "^3\n\n";

  test_fill();
  test_reductions();
  test_fill_boundary();
  test_arithmetic();
  test_laplacian();
  test_parallel_copy();
  test_device_to_host();
  time_kernels();

  return nFail;
}

//==========================================================
// Fill the valid cells on the device
void DemoChecks::test_fill() {
  const auto dx = geom.CellSizeArray();
  const auto exact = [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k)
                         noexcept { return exact_value(i, j, k, dx); };

  auto const& arrs = mf.arrays();
  ParallelFor(mf, [=] AMREX_GPU_DEVICE(int iBox, int i, int j, int k) {
    arrs[iBox](i, j, k) = exact(i, j, k);
  });
  Gpu::streamSynchronize();
  report("ParallelFor fill", max_error(mf, 0, exact), 1e-14);
}

//==========================================================
// Built-in reductions compared with a brute-force host reference
void DemoChecks::test_reductions() {
  const auto dx = geom.CellSizeArray();

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
                [=] AMREX_GPU_DEVICE(int iBox, int i, int j,
                                     int k) -> GpuTuple<Real, Real> {
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

//==========================================================
// Periodic ghost cell exchange across boxes and MPI ranks
void DemoChecks::test_fill_boundary() {
  const auto dx = geom.CellSizeArray();
  const auto exact = [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k)
                         noexcept { return exact_value(i, j, k, dx); };

  mf.setBndry(-999.);
  mf.FillBoundary(geom.periodicity());
  report("FillBoundary (periodic)", max_error(mf, nGhost, exact), 1e-13);
}

//==========================================================
// MultiFab arithmetic
void DemoChecks::test_arithmetic() {
  const auto dx = geom.CellSizeArray();

  MultiFab mf2(ba, dm, 1, 0);
  mf2.setVal(3.);
  MultiFab::Saxpy(mf2, 2., mf, 0, 0, 1, 0); // mf2 = 2*mf + 3
  report("MultiFab::Saxpy",
         max_error(mf2, 0,
                   [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k)
                       noexcept { return 2 * exact_value(i, j, k, dx) + 3; }),
         1e-13);

  MultiFab mf3(ba, dm, 1, 0);
  MultiFab::LinComb(mf3, 0.5, mf, 0, -1., mf2, 0, 0, 1, 0);
  report("MultiFab::LinComb",
         max_error(mf3, 0,
                   [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept {
                     return -1.5 * exact_value(i, j, k, dx) - 3;
                   }),
         1e-13);
}

//==========================================================
// Stencil kernel that reads the ghost cells
void DemoChecks::test_laplacian() {
  const auto dx = geom.CellSizeArray();

  apply_laplacian();
  report("7-point Laplacian",
         max_error(lap, 0,
                   [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k)
                       noexcept { return exact_laplacian(i, j, k, dx); }),
         1e-7);
}

//==========================================================
// Copy to a different BoxArray and DistributionMapping
void DemoChecks::test_parallel_copy() {
  const auto dx = geom.CellSizeArray();
  const auto exact = [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k)
                         noexcept { return exact_value(i, j, k, dx); };

  BoxArray baNew(domain);
  baNew.maxSize(maxGridSize / 2);
  const DistributionMapping dmNew(baNew);
  MultiFab mfNew(baNew, dmNew, 1, nGhost);
  mfNew.setVal(-999.);
  mfNew.ParallelCopy(mf, 0, 0, 1, 0, nGhost, geom.periodicity());
  report("ParallelCopy to " + std::to_string(baNew.size()) + " boxes",
         max_error(mfNew, nGhost, exact), 1e-13);
}

//==========================================================
// Copy to host memory and check on the CPU, independent of any GPU reduction
void DemoChecks::test_device_to_host() {
  const auto dx = geom.CellSizeArray();

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

//==========================================================
void DemoChecks::time_kernels() {
  const double nBytes = 2. * sizeof(Real) * nCell * nCell * nCell;

  Gpu::synchronize();
  ParallelDescriptor::Barrier();
  double t0 = second();
  for (int iter = 0; iter < nIter; ++iter)
    apply_laplacian();
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

//==========================================================
// lap = 7-point Laplacian of mf. The ghost cells of mf must be filled.
void DemoChecks::apply_laplacian() {
  const auto dx = geom.CellSizeArray();
  const Real dx2 = dx[0] * dx[0];

  auto const& dst = lap.arrays();
  auto const& src = mf.const_arrays();
  ParallelFor(lap, [=] AMREX_GPU_DEVICE(int iBox, int i, int j, int k) {
    auto const& a = src[iBox];
    dst[iBox](i, j, k) =
        (a(i - 1, j, k) + a(i + 1, j, k) + a(i, j - 1, k) + a(i, j + 1, k) +
         a(i, j, k - 1) + a(i, j, k + 1) - 6 * a(i, j, k)) /
        dx2;
  });
}

//==========================================================
void DemoChecks::report(const std::string& name, Real err, Real tol) {
  const bool isPass = err <= tol;
  if (!isPass)
    ++nFail;
  Print() << (isPass ? "  PASS  " : "  FAIL  ") << name << ": error = " << err
          << " (tol = " << tol << ")\n";
}

//==========================================================
void DemoChecks::print_device_info() const {
  Print() << "AMReX " << Version()
          << ", MPI ranks = " << ParallelDescriptor::NProcs() << "\n";
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

//==========================================================
