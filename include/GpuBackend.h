#ifndef _GPUBACKEND_H_
#define _GPUBACKEND_H_

// Helpers that let the same FLEKS kernels run on CPUs and GPUs.
//
// FLEKS keeps one code path: kernels are AMREX_GPU_HOST_DEVICE lambdas
// launched by amrex::ParallelFor, which compiles to a plain loop in a CPU
// build. Code that must differ between the two builds lives in this file
// only (tools/check_gpu_kernels.py enforces it), so the physics never forks.

#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include <AMReX_Extension.H>
#include <AMReX_FabArray.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuAsyncArray.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Loop.H>
#include <AMReX_MultiFab.H>
#include <AMReX_OpenMP.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>
#include <AMReX_iMultiFab.H>

#ifdef AMREX_USE_CUDA
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_reduce.cuh>
#include <cuda/functional>
#include <cuda/iterator>
#endif

// True in a GPU build. Use it for run-time decisions (e.g. rejecting a
// feature that is not available on the GPU) instead of #ifdef.
#ifdef AMREX_USE_GPU
inline constexpr bool fleksUseGpu = true;
#else
inline constexpr bool fleksUseGpu = false;
#endif

// True in a CPU OpenMP build.
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
inline constexpr bool fleksUseOmp = true;
#else
inline constexpr bool fleksUseOmp = false;
#endif

// MFIter options for a loop that only launches kernels whose results are
// read by later kernels: everything runs on the default stream, so the loop
// does not synchronize the device when it ends (AMReX otherwise spreads the
// boxes over several streams and waits for all of them).
//
// FLEKS loops use it by default: all FLEKS host copies are synchronous (they
// wait for the stream), so no host code depends on the end-of-loop sync, and
// a later AMReX loop over several streams synchronizes the default stream
// before it starts. Small boxes make FLEKS launch-bound, and this sync, one
// per loop, was a large part of the run time.
inline amrex::MFItInfo stream_ordered_mfi(bool doTiling = false) {
  amrex::MFItInfo info;
  if (doTiling)
    info.EnableTiling();
  return info.UseDefaultStream().DisableDeviceSync();
}

// The same for a particle iterator (PIter), which takes the options by
// reference. A particle container with tiles keeps its tiling.
inline amrex::MFItInfo& stream_ordered_pti() {
  static thread_local amrex::MFItInfo info;
  info = stream_ordered_mfi();
  return info;
}

// A scope in which all GPU work, AMReX's own loops included, runs on one
// stream, in launch order. Freed device temporaries then need no wait (see
// sync_before_free()), and DeferredSums can queue its sums. The scope ends
// with one stream synchronization, so the caller sees finished results.
//
// This does not use amrex::Gpu::NoSyncRegion: AMReX frees the pinned staging
// buffer of a FabArray's arrays() table without waiting for the copy that
// reads it, so with NoSyncRegion a temporary MultiFab can hand its table to
// another one before the copy runs (an illegal-address error in AMR runs).
//
// Use it around a whole FLEKS step (Domain::update). Nested scopes are
// allowed; only the outermost one synchronizes.
class StreamOrderedRegion {
#ifdef AMREX_USE_GPU
public:
  StreamOrderedRegion() : outermost(!amrex::Gpu::inSingleStreamRegion()) {}
  ~StreamOrderedRegion() {
    if (outermost)
      amrex::Gpu::streamSynchronize();
  }
  StreamOrderedRegion(const StreamOrderedRegion&) = delete;
  StreamOrderedRegion& operator=(const StreamOrderedRegion&) = delete;

private:
  bool outermost;
  amrex::Gpu::SingleStreamRegion singleStream;
#endif
};

// Call before device temporaries go out of scope while kernels launched on
// them may still run. The arena hands a freed block out again at once. With
// one stream (inside a StreamOrderedRegion) the next kernel that uses the
// block runs after the ones that read it, so there is nothing to wait for.
// With several streams (MFIter spreads boxes over them), wait.
inline void sync_before_free() {
#ifdef AMREX_USE_GPU
  if (amrex::Gpu::numGpuStreams() > 1)
    amrex::Gpu::streamSynchronizeAll();
#endif
}

// Particle tile layout. The CPU keeps one particle tile per cell: the
// particles of a cell are contiguous (cache locality) and the per-cell
// algorithms (resampling, sources) work on a tile. The GPU uses one tile per
// box, because every tile costs a kernel launch.
inline constexpr bool useCellTiles = !fleksUseGpu;

// Deposition strategy. The CPU adds every particle in order (deposit_add is a
// plain +=). The GPU sums the particles of each cell in a fixed order and
// gathers per node (CellBins, deposit_to_nodes), which is reproducible run to
// run without floating-point atomics.
inline constexpr bool useCellBinDeposits = fleksUseGpu;

// Capacity of a KernelArray. FLEKS runs use a handful of species or fluids.
inline constexpr int nKernelArrayMax = 16;

// A short array of POD values passed to a kernel by value (as a kernel
// argument, so no device allocation or copy is needed). Use it for per-species
// or per-fluid tables that live in host containers such as amrex::Vector,
// which must not be dereferenced on the device.
template <class T> struct KernelArray {
  T data[nKernelArrayMax];
  int n = 0;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& operator[](int i) {
    return data[i];
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& operator[](int i) const {
    return data[i];
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int size() const { return n; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool empty() const { return n == 0; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T* begin() { return data; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T* end() { return data + n; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T* begin() const {
    return data;
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T* end() const {
    return data + n;
  }

  // Host-side container operations, so a KernelArray can replace a short
  // amrex::Vector member. New elements are zero.
  void resize(int m) {
    if (m > nKernelArrayMax) {
      amrex::Abort("KernelArray: increase nKernelArrayMax in GpuBackend.h");
    }
    for (int i = n; i < m; ++i)
      data[i] = T{};
    n = m;
  }
  void clear() { n = 0; }
};

// Copy a host container (anything with size() and operator[]) into a
// KernelArray.
template <class T, class Container>
KernelArray<T> make_kernel_array(const Container& c) {
  KernelArray<T> a;
  a.n = static_cast<int>(c.size());
  if (a.n > nKernelArrayMax) {
    amrex::Abort("make_kernel_array: increase nKernelArrayMax in GpuBackend.h");
  }
  for (int i = 0; i < a.n; ++i)
    a.data[i] = static_cast<T>(c[i]);
  return a;
}

// Add `value` to `*target` from a kernel whose threads may hit the same
// target (particle-to-mesh deposition). It is an atomic addition on the GPU
// and a plain one on the CPU, where ParallelFor visits the particles in order,
// so a CPU deposit adds the contributions in the same order as a serial loop.
//
// Floating-point atomics make the order of the additions, and so the round-
// off, vary from run to run on the GPU. Deposits that feed the solution should
// eventually use a deterministic per-cell reduction instead.
template <class T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void deposit_add(T* target, T value) {
#if AMREX_DEVICE_COMPILE
  amrex::Gpu::Atomic::AddNoRet(target, value);
#else
  *target += value;
#endif
}

// Sum the copies of the nodes shared by boxes (SumBoundary). A GPU build adds
// them in a fixed order, then makes all the copies of a node equal: the
// copies are sums in different orders, and AMReX's GPU FillBoundary copies a
// ghost node covered by several boxes with racing stores, so unequal copies
// would make the result change from run to run. A CPU build is unchanged.
template <class FAB>
void sum_boundary(amrex::FabArray<FAB>& mf, const amrex::Periodicity& period) {
  mf.SumBoundary(period, /*deterministic=*/true);
  if (fleksUseGpu && mf.is_nodal())
    mf.OverrideSync(period);
}

// Sum of f(i), i = 0..n-1 (f a device lambda returning Real), returned on the
// host. A CPU build adds the values in index order, like a serial loop. A GPU
// build stores the values and reduces them with cub::DeviceReduce, whose fixed
// decomposition gives the same result from run to run (amrex::ReduceOps does
// not guarantee that).
#ifdef AMREX_USE_CUDA
// A pinned host Real that receives the result of deterministic_sum().
inline amrex::Real* deterministic_sum_result() {
  static amrex::Real* p = static_cast<amrex::Real*>(
      amrex::The_Pinned_Arena()->alloc(sizeof(amrex::Real)));
  return p;
}
#endif

#ifdef AMREX_USE_CUDA
// Device scratch for the CUB reductions, kept for the whole run. Grows as
// needed.
inline void* reduction_scratch(std::size_t nBytes) {
  static void* p = nullptr;
  static std::size_t size = 0;
  if (nBytes > size) {
    if (p) {
      // The arena may hand the block out again at once: wait for the kernels
      // that still use it.
      amrex::Gpu::streamSynchronize();
      amrex::The_Arena()->free(p);
    }
    size = std::max(nBytes, 2 * size);
    p = amrex::The_Arena()->alloc(size);
  }
  return p;
}
#endif

// Fixed block size for CPU deterministic reductions and minimum work thresholds
// before spawning OpenMP threads (so small 1D/coarse loops stay on one thread
// with SIMD and avoid fork-join overhead).
inline constexpr int kDetSumBlockSize = 256;
inline constexpr int kMinOmpBoxCells = 2048;
inline constexpr int kMinOmp1DSize = 8192;
inline constexpr int kMinOmpTiles = 32;
inline constexpr int kMinOmpParticles = 4096;

// Performance-portable box kernel launch with size-gated OpenMP threading.
template <class L>
void parallel_for_omp(amrex::Box const& box, L const& f) noexcept {
#ifdef AMREX_USE_GPU
  amrex::ParallelFor(box, f);
#elif defined(AMREX_USE_OMP)
  if (omp_get_max_threads() <= 1 || amrex::OpenMP::in_parallel() ||
      box.numPts() < kMinOmpBoxCells) {
    amrex::ParallelFor(box, f);
    return;
  }
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);
#if (AMREX_SPACEDIM == 1)
#pragma omp parallel for schedule(static)
  for (int i = lo.x; i <= hi.x; ++i) {
    f(i, 0, 0);
  }
#elif (AMREX_SPACEDIM == 2)
  const int nj = hi.y - lo.y + 1;
  if (nj >= 2 * omp_get_max_threads()) {
#pragma omp parallel for schedule(static)
    for (int j = lo.y; j <= hi.y; ++j) {
      AMREX_PRAGMA_SIMD
      for (int i = lo.x; i <= hi.x; ++i) {
        f(i, j, 0);
      }
    }
  } else {
    const int ni = hi.x - lo.x + 1;
    const int nTot = nj * ni;
#pragma omp parallel for schedule(static)
    for (int idx = 0; idx < nTot; ++idx) {
      f(lo.x + (idx % ni), lo.y + (idx / ni), 0);
    }
  }
#else
  const int nk = hi.z - lo.z + 1;
  const int nj = hi.y - lo.y + 1;
  if (nk * nj >= 2 * omp_get_max_threads()) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
          f(i, j, k);
        }
      }
    }
  } else {
    const int ni = hi.x - lo.x + 1;
    const int nTot = nk * nj * ni;
#pragma omp parallel for schedule(static)
    for (int idx = 0; idx < nTot; ++idx) {
      const int i = lo.x + (idx % ni);
      const int t = idx / ni;
      f(i, lo.y + (t % nj), lo.z + (t / nj));
    }
  }
#endif
#else
  amrex::ParallelFor(box, f);
#endif
}

template <class T, class L, class M = std::enable_if_t<std::is_integral_v<T> > >
void parallel_for_omp(amrex::Box const& box, T ncomp, L const& f) noexcept {
#ifdef AMREX_USE_GPU
  amrex::ParallelFor(box, ncomp, f);
#elif defined(AMREX_USE_OMP)
  if (omp_get_max_threads() <= 1 || amrex::OpenMP::in_parallel() ||
      box.numPts() * static_cast<long>(ncomp) < kMinOmpBoxCells) {
    amrex::ParallelFor(box, ncomp, f);
    return;
  }
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);
#if (AMREX_SPACEDIM == 1)
#pragma omp parallel for collapse(2) schedule(static)
  for (T n = 0; n < ncomp; ++n) {
    for (int i = lo.x; i <= hi.x; ++i) {
      f(i, 0, 0, n);
    }
  }
#elif (AMREX_SPACEDIM == 2)
  const int nj = hi.y - lo.y + 1;
  if (static_cast<long>(ncomp) * nj >= 2 * omp_get_max_threads()) {
#pragma omp parallel for collapse(2) schedule(static)
    for (T n = 0; n < ncomp; ++n) {
      for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
          f(i, j, 0, n);
        }
      }
    }
  } else {
    const int ni = hi.x - lo.x + 1;
    const long nPlane = static_cast<long>(nj) * ni;
    const long nTot = static_cast<long>(ncomp) * nPlane;
#pragma omp parallel for schedule(static)
    for (long idx = 0; idx < nTot; ++idx) {
      const T n = static_cast<T>(idx / nPlane);
      const long rem = idx % nPlane;
      f(lo.x + static_cast<int>(rem % ni), lo.y + static_cast<int>(rem / ni), 0,
        n);
    }
  }
#else
  const int nk = hi.z - lo.z + 1;
  const int nj = hi.y - lo.y + 1;
  if (static_cast<long>(ncomp) * nk * nj >= 2 * omp_get_max_threads()) {
#pragma omp parallel for collapse(3) schedule(static)
    for (T n = 0; n < ncomp; ++n) {
      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
          AMREX_PRAGMA_SIMD
          for (int i = lo.x; i <= hi.x; ++i) {
            f(i, j, k, n);
          }
        }
      }
    }
  } else {
    const int ni = hi.x - lo.x + 1;
    const long nBox = static_cast<long>(nk) * nj * ni;
    const long nTot = static_cast<long>(ncomp) * nBox;
#pragma omp parallel for schedule(static)
    for (long idx = 0; idx < nTot; ++idx) {
      const T n = static_cast<T>(idx / nBox);
      const long rem = idx % nBox;
      const int i = lo.x + static_cast<int>(rem % ni);
      const long t = rem / ni;
      f(i, lo.y + static_cast<int>(t % nj), lo.z + static_cast<int>(t / nj), n);
    }
  }
#endif
#else
  amrex::ParallelFor(box, ncomp, f);
#endif
}

// Launch a box-loop across all local boxes of a MultiFab/FabArray in a single
// OpenMP parallel region on CPU (cutting per-box fork-join overhead by nBoxes),
// or per-box GPU kernels when AMREX_USE_GPU is enabled.
// `makeBox(mfi)` returns the Box for `mfi`; `makeKernel(mfi)` returns the
// `(int i, int j, int k)` functor.
template <class MF, class MakeBox, class MakeKernel>
void parallel_for_mf_omp(const MF& mf, MakeBox&& makeBox,
                         MakeKernel&& makeKernel) {
#ifdef AMREX_USE_GPU
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    amrex::ParallelFor(makeBox(mfi), makeKernel(mfi));
  }
#elif defined(AMREX_USE_OMP)
  if (omp_get_max_threads() <= 1 || amrex::OpenMP::in_parallel()) {
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
      amrex::ParallelFor(makeBox(mfi), makeKernel(mfi));
    }
    return;
  }
  using KernelT = decltype(makeKernel(std::declval<const amrex::MFIter&>()));
  std::vector<amrex::Box> boxes;
  std::vector<KernelT> kernels;
  long totalPts = 0;
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    amrex::Box bx = makeBox(mfi);
    if (bx.ok() && !bx.isEmpty()) {
      totalPts += bx.numPts();
      boxes.push_back(bx);
      kernels.push_back(makeKernel(mfi));
    }
  }
  const int nBoxes = static_cast<int>(boxes.size());
  if (nBoxes == 0)
    return;

  if (amrex::OpenMP::in_parallel() || totalPts < kMinOmpBoxCells) {
    for (int b = 0; b < nBoxes; ++b) {
      const auto lo = amrex::lbound(boxes[b]);
      const auto hi = amrex::ubound(boxes[b]);
      const auto& kern = kernels[b];
      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
          AMREX_PRAGMA_SIMD
          for (int i = lo.x; i <= hi.x; ++i) {
            kern(i, j, k);
          }
        }
      }
    }
  } else {
    std::vector<int> rowOffsets(nBoxes + 1, 0);
    for (int b = 0; b < nBoxes; ++b) {
      const auto lo = amrex::lbound(boxes[b]);
      const auto hi = amrex::ubound(boxes[b]);
      rowOffsets[b + 1] =
          rowOffsets[b] + (hi.y - lo.y + 1) * (hi.z - lo.z + 1);
    }
    const int totalRows = rowOffsets[nBoxes];
#pragma omp parallel for schedule(static)
    for (int globalRow = 0; globalRow < totalRows; ++globalRow) {
      int b = 0;
      while (b + 1 < nBoxes && globalRow >= rowOffsets[b + 1])
        ++b;
      const int localRow = globalRow - rowOffsets[b];
      const auto lo = amrex::lbound(boxes[b]);
      const auto hi = amrex::ubound(boxes[b]);
      const int nj = hi.y - lo.y + 1;
      const int j = lo.y + (localRow % nj);
      const int k = lo.z + (localRow / nj);
      const auto& kern = kernels[b];
      AMREX_PRAGMA_SIMD
      for (int i = lo.x; i <= hi.x; ++i) {
        kern(i, j, k);
      }
    }
  }
#else
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    amrex::ParallelFor(makeBox(mfi), makeKernel(mfi));
  }
#endif
}

template <class MF, class MakeBox, class MakeKernel>
void parallel_for_mf_omp(const MF& mf, int ncomp, MakeBox&& makeBox,
                         MakeKernel&& makeKernel) {
  if (ncomp <= 0)
    return;
#ifdef AMREX_USE_GPU
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    amrex::ParallelFor(makeBox(mfi), ncomp, makeKernel(mfi));
  }
#elif defined(AMREX_USE_OMP)
  if (omp_get_max_threads() <= 1 || amrex::OpenMP::in_parallel()) {
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
      amrex::ParallelFor(makeBox(mfi), ncomp, makeKernel(mfi));
    }
    return;
  }
  using KernelT = decltype(makeKernel(std::declval<const amrex::MFIter&>()));
  std::vector<amrex::Box> boxes;
  std::vector<KernelT> kernels;
  long totalPts = 0;
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    amrex::Box bx = makeBox(mfi);
    if (bx.ok() && !bx.isEmpty()) {
      totalPts += bx.numPts();
      boxes.push_back(bx);
      kernels.push_back(makeKernel(mfi));
    }
  }
  const int nBoxes = static_cast<int>(boxes.size());
  if (nBoxes == 0)
    return;

  if (amrex::OpenMP::in_parallel() ||
      totalPts * static_cast<long>(ncomp) < kMinOmpBoxCells) {
    for (int b = 0; b < nBoxes; ++b) {
      const auto lo = amrex::lbound(boxes[b]);
      const auto hi = amrex::ubound(boxes[b]);
      const auto& kern = kernels[b];
      for (int c = 0; c < ncomp; ++c) {
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
              kern(i, j, k, c);
            }
          }
        }
      }
    }
  } else {
    std::vector<int> rowOffsets(nBoxes + 1, 0);
    for (int b = 0; b < nBoxes; ++b) {
      const auto lo = amrex::lbound(boxes[b]);
      const auto hi = amrex::ubound(boxes[b]);
      rowOffsets[b + 1] =
          rowOffsets[b] + ncomp * (hi.y - lo.y + 1) * (hi.z - lo.z + 1);
    }
    const int totalRows = rowOffsets[nBoxes];
#pragma omp parallel for schedule(static)
    for (int globalRow = 0; globalRow < totalRows; ++globalRow) {
      int b = 0;
      while (b + 1 < nBoxes && globalRow >= rowOffsets[b + 1])
        ++b;
      const int localRow = globalRow - rowOffsets[b];
      const auto lo = amrex::lbound(boxes[b]);
      const auto hi = amrex::ubound(boxes[b]);
      const int nj = hi.y - lo.y + 1;
      const int nk = hi.z - lo.z + 1;
      const int njk = nj * nk;
      const int c = localRow / njk;
      const int rem = localRow % njk;
      const int j = lo.y + (rem % nj);
      const int k = lo.z + (rem / nj);
      const auto& kern = kernels[b];
      AMREX_PRAGMA_SIMD
      for (int i = lo.x; i <= hi.x; ++i) {
        kern(i, j, k, c);
      }
    }
  }
#else
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    amrex::ParallelFor(makeBox(mfi), ncomp, makeKernel(mfi));
  }
#endif
}

// OpenMP-accelerated MultiFab arithmetic helpers on CPU (since AMReX
// MultiFab::Add/Saxpy/mult are serial when AMReX USE_OMP=FALSE).
inline void mf_add(amrex::MultiFab& dst, const amrex::MultiFab& src, int scomp,
                   int dcomp, int ncomp, int ngrow) {
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
  parallel_for_mf_omp(
      dst, ncomp,
      [ngrow](const amrex::MFIter& mfi) {
        return amrex::grow(mfi.validbox(), ngrow);
      },
      [&](const amrex::MFIter& mfi) {
        const auto d = dst[mfi].array();
        const auto s = src[mfi].const_array();
        return [=](int i, int j, int k, int c) noexcept {
          d(i, j, k, dcomp + c) += s(i, j, k, scomp + c);
        };
      });
#else
  amrex::MultiFab::Add(dst, src, scomp, dcomp, ncomp, ngrow);
#endif
}

inline void mf_saxpy(amrex::MultiFab& dst, amrex::Real a,
                     const amrex::MultiFab& src, int scomp, int dcomp,
                     int ncomp, int ngrow) {
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
  parallel_for_mf_omp(
      dst, ncomp,
      [ngrow](const amrex::MFIter& mfi) {
        return amrex::grow(mfi.validbox(), ngrow);
      },
      [&](const amrex::MFIter& mfi) {
        const auto d = dst[mfi].array();
        const auto s = src[mfi].const_array();
        return [=](int i, int j, int k, int c) noexcept {
          d(i, j, k, dcomp + c) += a * s(i, j, k, scomp + c);
        };
      });
#else
  amrex::MultiFab::Saxpy(dst, a, src, scomp, dcomp, ncomp, ngrow);
#endif
}

inline void mf_mult(amrex::MultiFab& dst, amrex::Real val, int comp, int ncomp,
                    int ngrow) {
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
  parallel_for_mf_omp(
      dst, ncomp,
      [ngrow](const amrex::MFIter& mfi) {
        return amrex::grow(mfi.validbox(), ngrow);
      },
      [&](const amrex::MFIter& mfi) {
        const auto d = dst[mfi].array();
        return [=](int i, int j, int k, int c) noexcept {
          d(i, j, k, comp + c) *= val;
        };
      });
#else
  dst.mult(val, comp, ncomp, ngrow);
#endif
}

inline void mf_copy(amrex::MultiFab& dst, const amrex::MultiFab& src, int scomp,
                    int dcomp, int ncomp, int ngrow) {
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
  parallel_for_mf_omp(
      dst, ncomp,
      [ngrow](const amrex::MFIter& mfi) {
        return amrex::grow(mfi.validbox(), ngrow);
      },
      [&](const amrex::MFIter& mfi) {
        const auto d = dst[mfi].array();
        const auto s = src[mfi].const_array();
        return [=](int i, int j, int k, int c) noexcept {
          d(i, j, k, dcomp + c) = s(i, j, k, scomp + c);
        };
      });
#else
  amrex::MultiFab::Copy(dst, src, scomp, dcomp, ncomp, ngrow);
#endif
}

// Launch a 1D loop f(i) for i = 0..n-1 on the GPU, across OpenMP threads on
// the CPU when n is large enough and not already inside an OpenMP parallel
// region, or as a SIMD loop on the CPU.
template <class F> void parallel_for_1d(int n, F&& f) {
  if (n <= 0)
    return;
#ifdef AMREX_USE_GPU
  amrex::ParallelFor(n, std::forward<F>(f));
#elif defined(AMREX_USE_OMP)
  if (omp_get_max_threads() > 1 && !amrex::OpenMP::in_parallel() &&
      n >= kMinOmp1DSize) {
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
      f(i);
  } else {
    AMREX_PRAGMA_SIMD
    for (int i = 0; i < n; ++i)
      f(i);
  }
#else
  amrex::ParallelFor(n, std::forward<F>(f));
#endif
}

// Run host work items f(i), i = 0..n-1, across OpenMP threads (or serially
// without OpenMP). Each item must write a disjoint set of data, so the result
// does not depend on the thread count. Use `dynamic` when the items have very
// different costs. This is a host-only helper (f may use host data
// structures); it never launches a GPU kernel.
template <class F> void parallel_for_host_items(int n, F&& f, bool dynamic) {
  if (n <= 0)
    return;
#ifdef AMREX_USE_OMP
  if (omp_get_max_threads() > 1 && !amrex::OpenMP::in_parallel() &&
      n >= kMinOmpTiles) {
    if (dynamic) {
#pragma omp parallel for schedule(dynamic, 1)
      for (int i = 0; i < n; ++i)
        f(i);
    } else {
#pragma omp parallel for schedule(static)
      for (int i = 0; i < n; ++i)
        f(i);
    }
    return;
  }
#endif
  amrex::ignore_unused(dynamic);
  for (int i = 0; i < n; ++i)
    f(i);
}

// Run f() with a single OpenMP thread if cond is true (no-op otherwise, and
// without OpenMP). Used where an AMReX CPU routine is only deterministic per
// thread count: e.g. ParticleContainer::Redistribute gathers the particles
// that change tile into per-thread buffers level by level, so with several
// AMR levels the order of the particles in a destination tile (and hence the
// summation order of every later deposit) depends on OMP_NUM_THREADS.
template <class F> void run_single_threaded_if(bool cond, F&& f) {
  if (!cond || fleksUseGpu || amrex::OpenMP::in_parallel()) {
    f();
    return;
  }
  const int nThreads = amrex::OpenMP::get_max_threads();
  amrex::OpenMP::set_num_threads(1);
  f();
  amrex::OpenMP::set_num_threads(nThreads);
}

// Deterministic sum of N per-particle quantities over the local particles of
// level iLev, for the CPU (host particle memory). Each particle tile is
// summed sequentially into its own partial sum, the tiles are processed in
// parallel, and the partial sums are added in the fixed (grid, tile) order of
// the container, so the result does not depend on OMP_NUM_THREADS.
// f(const ParticleType& p, std::array<Real, N>& s) adds the particle to s.
template <int N, class PC, class F>
std::array<amrex::Real, N> deterministic_particle_sum(const PC& pc, int iLev,
                                                      F&& f) {
  using TileType = typename PC::ParticleTileType;
  std::array<amrex::Real, N> total{};
  if (iLev >= pc.numLevels())
    return total;
  const auto& plev = pc.GetParticles(iLev);
  std::vector<const TileType*> tiles;
  tiles.reserve(plev.size());
  for (const auto& kv : plev)
    tiles.push_back(&kv.second);

  const int nTile = static_cast<int>(tiles.size());
  std::vector<std::array<amrex::Real, N> > partial(nTile);
  parallel_for_host_items(
      nTile,
      [&](int it) {
        std::array<amrex::Real, N> s{};
        const TileType& ptile = *tiles[it];
        const auto* pstruct = ptile.GetArrayOfStructs()().data();
        const int np = static_cast<int>(ptile.numParticles());
        for (int ip = 0; ip < np; ++ip)
          f(pstruct[ip], s);
        partial[it] = s;
      },
      true);

  for (int it = 0; it < nTile; ++it)
    for (int c = 0; c < N; ++c)
      total[c] += partial[it][c];
  return total;
}

// Descriptor for a non-empty cell-sized particle tile (tile_size = {1, 1, 1})
// on the CPU. `grid` is the FAB index in ParticleBoxArray(iLev), `tile` is the
// local tile index, `cell` is the cell IntVect of the 1x1x1 tile box, and
// `ptile` points to the ParticleTile.
template <class TileType> struct CellTileItem {
  int grid = 0;
  int tile = 0;
  amrex::IntVect cell;
  TileType* ptile = nullptr;
};

inline bool is_home_nodal(const amrex::IntVect& loIdx,
                          const amrex::IntVect& cell) {
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (loIdx[d] < cell[d] - 1 || loIdx[d] > cell[d])
      return false;
  }
  return true;
}

// Collect all non-empty particle tiles at level iLev in deterministic
// (grid, tile) order, computing each tile's cell IntVect without red-black
// tree lookups per cell.
template <class PC> auto collect_cell_tiles(PC& pc, int iLev) {
  using TileType = std::conditional_t<std::is_const_v<PC>,
                                      const typename PC::ParticleTileType,
                                      typename PC::ParticleTileType>;
  std::vector<CellTileItem<TileType> > items;
  if (iLev >= pc.numLevels())
    return items;
  auto& plev = pc.GetParticles(iLev);
  if (plev.empty())
    return items;
  const amrex::BoxArray& ba = pc.ParticleBoxArray(iLev);
  items.reserve(plev.size());
  int lastGrid = -1;
  amrex::IntVect sml;
  amrex::IntVect len;
  for (auto& kv : plev) {
    if (kv.second.numParticles() == 0)
      continue;
    const int grid = kv.first.first;
    const int tile = kv.first.second;
    if (grid != lastGrid) {
      lastGrid = grid;
      const amrex::Box& bx = ba[grid];
      sml = bx.smallEnd();
      len = bx.length();
    }
    const amrex::IntVect cell(AMREX_D_DECL(
        sml[0] + (tile % len[0]), sml[1] + ((tile / len[0]) % len[1]),
        sml[2] + (tile / (len[0] * len[1]))));
    items.push_back(CellTileItem<TileType>{ grid, tile, cell, &kv.second });
  }
  return items;
}

template <class TileType>
inline long count_tile_particles(
    const std::vector<CellTileItem<TileType> >& items) {
  long total = 0;
  for (const auto& it : items)
    total += it.ptile->numParticles();
  return total;
}

// Run f(const CellTileItem<TileType>& item, int ord) across OpenMP threads for
// all non-empty cell tiles in `items` (in deterministic `ord = 0..nTiles-1`).
template <class TileType, class F>
void parallel_for_cell_tiles(const std::vector<CellTileItem<TileType> >& items,
                             F&& f) {
  const int nTiles = static_cast<int>(items.size());
  if (nTiles == 0)
    return;
#ifdef AMREX_USE_OMP
  if (!amrex::OpenMP::in_parallel() &&
      (nTiles >= kMinOmpTiles * 4 ||
       (nTiles >= 2 && count_tile_particles(items) >= kMinOmpParticles))) {
#pragma omp parallel for schedule(guided)
    for (int ord = 0; ord < nTiles; ++ord)
      f(items[ord], ord);
    return;
  }
#endif
  for (int ord = 0; ord < nTiles; ++ord)
    f(items[ord], ord);
}

template <class PC, class F>
void parallel_for_cell_tiles(PC& pc, int iLev, F&& f) {
  const auto items = collect_cell_tiles(pc, iLev);
  parallel_for_cell_tiles(items, std::forward<F>(f));
}

// Deterministic colored cell-tile deposition over `items` with stencil span
// `stencilSpan` per dimension (e.g. 2 for 2-node CIC, 3 for 3-cell CIC, 4 for
// Esirkepov current).
//
// Pass 1 partitions tiles into `stencilSpan^D` colors by `cell[d] mod
// stencilSpan` and runs each color in parallel across OpenMP threads: two
// tiles of the same color are separated by at least `stencilSpan` cells in at
// least one dimension and never touch the same mesh target. `f(item, ord,
// /*strayPass=*/false)` processes all particles whose stencil is confined to
// `item.cell` and returns true if the tile also holds any out-of-cell ("stray")
// particles (such as boundary ghost particles or freshly injected face
// particles before Redistribute).
//
// Pass 2 runs `f(item, ord, /*strayPass=*/true)` sequentially in deterministic
// `ord` order for only those tiles that returned true in Pass 1.
template <class TileType, class F>
void deposit_cell_tiles_colored(
    const std::vector<CellTileItem<TileType> >& items, int stencilSpan, F&& f) {
  const int nTiles = static_cast<int>(items.size());
  if (nTiles == 0)
    return;

#ifdef AMREX_USE_OMP
  if (omp_get_max_threads() > 1 && !amrex::OpenMP::in_parallel() &&
      (nTiles >= kMinOmpTiles * 4 ||
       (nTiles >= 4 && count_tile_particles(items) >= kMinOmpParticles))) {
    std::vector<unsigned char> hasStray(nTiles, 0);

    // Fast path: group cell tiles by (grid, highest-dim plane) and color
    // planes by `planeCoord % stencilSpan`. This reduces OpenMP barriers from
    // stencilSpan^nDim (27 in 3D) down to stencilSpan (3), eliminates false
    // sharing across threads, and gives each thread contiguous in-plane cells.
    bool usedPlaneStriping = false;
    if constexpr (AMREX_SPACEDIM >= 2) {
      constexpr int hiDim = AMREX_SPACEDIM - 1;
      int minGrid = items[0].grid, maxGrid = items[0].grid;
      int minK = items[0].cell[hiDim], maxK = items[0].cell[hiDim];
      for (int ord = 1; ord < nTiles; ++ord) {
        minGrid = std::min(minGrid, items[ord].grid);
        maxGrid = std::max(maxGrid, items[ord].grid);
        minK = std::min(minK, items[ord].cell[hiDim]);
        maxK = std::max(maxK, items[ord].cell[hiDim]);
      }
      const long gridSpan = static_cast<long>(maxGrid - minGrid + 1);
      const long kSpan = static_cast<long>(maxK - minK + 1);
      if (gridSpan * kSpan <= 65536) {
        const int tableSize = static_cast<int>(gridSpan * kSpan);
        std::vector<int> keyToPlane(tableSize, -1);
        std::vector<int> planeColor;
        std::vector<int> planeCounts;
        std::vector<int> itemPlane(nTiles);
        for (int ord = 0; ord < nTiles; ++ord) {
          const int g = items[ord].grid - minGrid;
          const int k = items[ord].cell[hiDim];
          const int key = g * static_cast<int>(kSpan) + (k - minK);
          int p = keyToPlane[key];
          if (p < 0) {
            p = static_cast<int>(planeColor.size());
            keyToPlane[key] = p;
            planeColor.push_back(((k % stencilSpan) + stencilSpan) %
                                 stencilSpan);
            planeCounts.push_back(0);
          }
          itemPlane[ord] = p;
          ++planeCounts[p];
        }
        const int nPlanes = static_cast<int>(planeColor.size());
        const int nThreads = omp_get_max_threads();
        if (nPlanes >= stencilSpan * nThreads) {
          usedPlaneStriping = true;
          std::vector<int> planeStarts(nPlanes + 1, 0);
          for (int p = 0; p < nPlanes; ++p)
            planeStarts[p + 1] = planeStarts[p] + planeCounts[p];
          std::vector<int> planeCursor(planeStarts.begin(),
                                       planeStarts.begin() + nPlanes);
          std::vector<int> planeOrds(nTiles);
          for (int ord = 0; ord < nTiles; ++ord)
            planeOrds[planeCursor[itemPlane[ord]]++] = ord;

          std::vector<int> colorCounts(stencilSpan, 0);
          for (int p = 0; p < nPlanes; ++p)
            ++colorCounts[planeColor[p]];
          std::vector<int> colorStarts(stencilSpan + 1, 0);
          for (int c = 0; c < stencilSpan; ++c)
            colorStarts[c + 1] = colorStarts[c] + colorCounts[c];
          std::vector<int> colorCursor(colorStarts.begin(),
                                       colorStarts.begin() + stencilSpan);
          std::vector<int> colorPlanes(nPlanes);
          for (int p = 0; p < nPlanes; ++p)
            colorPlanes[colorCursor[planeColor[p]]++] = p;

#pragma omp parallel
          {
            for (int color = 0; color < stencilSpan; ++color) {
              const int beg = colorStarts[color];
              const int end = colorStarts[color + 1];
              if (beg == end)
                continue;
#pragma omp for schedule(dynamic, 1)
              for (int pIdx = beg; pIdx < end; ++pIdx) {
                const int p = colorPlanes[pIdx];
                for (int idx = planeStarts[p]; idx < planeStarts[p + 1];
                     ++idx) {
                  const int ord = planeOrds[idx];
                  if (f(items[ord], ord, false))
                    hasStray[ord] = 1;
                }
              }
            }
          }
        }
      }
    }

    if (!usedPlaneStriping) {
      int nColors = 1;
      for (int d = 0; d < AMREX_SPACEDIM; ++d)
        nColors *= stencilSpan;

      std::vector<int> bucketCounts(nColors, 0);
      std::vector<int> tileColors(nTiles);
      for (int ord = 0; ord < nTiles; ++ord) {
        const auto& c = items[ord].cell;
        int color = 0;
        int stride = 1;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          const int rem = ((c[d] % stencilSpan) + stencilSpan) % stencilSpan;
          color += rem * stride;
          stride *= stencilSpan;
        }
        tileColors[ord] = color;
        ++bucketCounts[color];
      }

      std::vector<int> bucketStarts(nColors + 1, 0);
      for (int c = 0; c < nColors; ++c)
        bucketStarts[c + 1] = bucketStarts[c] + bucketCounts[c];

      std::vector<int> bucketCursor(bucketStarts.begin(),
                                    bucketStarts.begin() + nColors);
      std::vector<int> bucketOrds(nTiles);
      for (int ord = 0; ord < nTiles; ++ord)
        bucketOrds[bucketCursor[tileColors[ord]]++] = ord;

#pragma omp parallel
      {
        for (int color = 0; color < nColors; ++color) {
          const int beg = bucketStarts[color];
          const int end = bucketStarts[color + 1];
          if (beg == end)
            continue;
#pragma omp for schedule(guided)
          for (int idx = beg; idx < end; ++idx) {
            const int ord = bucketOrds[idx];
            if (f(items[ord], ord, false))
              hasStray[ord] = 1;
          }
        }
      }
    }

    for (int ord = 0; ord < nTiles; ++ord) {
      if (hasStray[ord])
        f(items[ord], ord, true);
    }
    return;
  }
#endif

  for (int ord = 0; ord < nTiles; ++ord) {
    if (f(items[ord], ord, false))
      f(items[ord], ord, true);
  }
}

template <class PC, class F>
void deposit_cell_tiles_colored(PC& pc, int iLev, int stencilSpan, F&& f) {
  const auto items = collect_cell_tiles(pc, iLev);
  deposit_cell_tiles_colored(items, stencilSpan, std::forward<F>(f));
}

template <class TileType, class F>
bool deposit_cell_tiles_thread_private(
    const std::vector<CellTileItem<TileType> >& items, amrex::MultiFab& mf,
    int nComp, F&& depositTileFn) {
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
  const int nThreads = omp_get_max_threads();
  if (amrex::OpenMP::in_parallel() || nThreads <= 1 ||
      static_cast<int>(items.size()) < kMinOmpTiles * 4) {
    return false;
  }

  const int nGrids = mf.boxArray().size();
  std::vector<int> gridToLocal(nGrids, -1);
  std::vector<amrex::FArrayBox*> dstFabs;
  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    gridToLocal[mfi.index()] = static_cast<int>(dstFabs.size());
    dstFabs.push_back(&mf[mfi]);
  }
  const int nLocal = static_cast<int>(dstFabs.size());
  std::vector<std::vector<int> > boxItems(nLocal);
  for (int ord = 0; ord < static_cast<int>(items.size()); ++ord) {
    const int lb = gridToLocal[items[ord].grid];
    if (lb >= 0)
      boxItems[lb].push_back(ord);
  }

  std::vector<amrex::FArrayBox> threadFabs(nThreads);
#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    amrex::FArrayBox& localFab = threadFabs[tid];

    for (int b = 0; b < nLocal; ++b) {
      const auto& bItems = boxItems[b];
      if (bItems.empty())
        continue;
      amrex::FArrayBox* dstFab = dstFabs[b];
      const amrex::Box fabBox = dstFab->box();
      localFab.resize(fabBox, nComp);
      localFab.setVal<amrex::RunOn::Host>(0.0);
      const amrex::Array4<amrex::Real> localArr = localFab.array();
      const int nBItems = static_cast<int>(bItems.size());

#pragma omp for schedule(guided)
      for (int idx = 0; idx < nBItems; ++idx) {
        const int ord = bItems[idx];
        depositTileFn(items[ord], ord, localArr);
      }

      const amrex::Array4<amrex::Real> dstArr = dstFab->array();
      const auto blo = amrex::lbound(fabBox);
      const auto bhi = amrex::ubound(fabBox);
      const int nj = bhi.y - blo.y + 1;
      const int nk = bhi.z - blo.z + 1;
      const int njk = nj * nk;
      const int totalRows = nComp * njk;
      const int ni = bhi.x - blo.x + 1;
#pragma omp for schedule(static)
      for (int row = 0; row < totalRows; ++row) {
        const int c = row / njk;
        const int rem = row % njk;
        const int j = blo.y + (rem % nj);
        const int k = blo.z + (rem / nj);
        amrex::Real* AMREX_RESTRICT dPtr = dstArr.ptr(blo.x, j, k, c);
        for (int t = 0; t < nThreads; ++t) {
          const amrex::Real* AMREX_RESTRICT sPtr =
              threadFabs[t].const_array().ptr(blo.x, j, k, c);
          AMREX_PRAGMA_SIMD
          for (int ii = 0; ii < ni; ++ii) {
            dPtr[ii] += sPtr[ii];
          }
        }
      }
    }
  }
  return true;
#else
  amrex::ignore_unused(items, mf, nComp, depositTileFn);
  return false;
#endif
}

template <class F> amrex::Real deterministic_sum(int n, F&& f) {
#ifdef AMREX_USE_CUDA
  if (n <= 0)
    return 0;
  // CUB reads the values from f directly (no temporary array); the result
  // goes to pinned memory, so one stream synchronization returns it.
  auto values = cuda::make_transform_iterator(
      cuda::counting_iterator<int>(0),
      cuda::proclaim_return_type<amrex::Real>(std::decay_t<F>(f)));
  amrex::Real* const result = deterministic_sum_result();
  std::size_t tmpBytes = 0;
  cub::DeviceReduce::Sum(nullptr, tmpBytes, values, result, n,
                         amrex::Gpu::gpuStream());
  cub::DeviceReduce::Sum(reduction_scratch(tmpBytes), tmpBytes, values, result,
                         n, amrex::Gpu::gpuStream());
  amrex::Gpu::streamSynchronize();
  return *result;
#else
  if (n <= 0)
    return 0;
  if (n <= kDetSumBlockSize) {
    amrex::Real sum = 0;
    for (int i = 0; i < n; ++i)
      sum += f(i);
    return sum;
  }
  const int nBlocks = (n + kDetSumBlockSize - 1) / kDetSumBlockSize;
  constexpr int kMaxStackBlocks = 256;
  amrex::Real stackBlocks[kMaxStackBlocks];
  std::vector<amrex::Real> heapBlocks;
  if (nBlocks > kMaxStackBlocks)
    heapBlocks.resize(nBlocks);
  amrex::Real* const blockSums =
      (nBlocks <= kMaxStackBlocks) ? stackBlocks : heapBlocks.data();
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(static) if (!amrex::OpenMP::in_parallel() && \
                                                  n >= kMinOmp1DSize)
#endif
  for (int b = 0; b < nBlocks; ++b) {
    const int iStart = b * kDetSumBlockSize;
    const int iEnd = std::min(n, iStart + kDetSumBlockSize);
    amrex::Real localSum = 0;
    for (int i = iStart; i < iEnd; ++i)
      localSum += f(i);
    blockSums[b] = localSum;
  }
  amrex::Real sum = 0;
  for (int b = 0; b < nBlocks; ++b)
    sum += blockSums[b];
  return sum;
#endif
}

// Sum of f(i), i = 0..n-1, written to *result, which is device memory in a
// GPU build: unlike deterministic_sum() this does not wait for the GPU, so
// later kernels can read the result without a round trip to the host. The
// summation order is that of deterministic_sum().
template <class F>
void deterministic_sum_to(int n, F&& f, amrex::Real* result) {
#ifdef AMREX_USE_CUDA
  if (n <= 0) {
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) { *result = 0; });
    return;
  }
  auto values = cuda::make_transform_iterator(
      cuda::counting_iterator<int>(0),
      cuda::proclaim_return_type<amrex::Real>(std::decay_t<F>(f)));
  std::size_t tmpBytes = 0;
  cub::DeviceReduce::Sum(nullptr, tmpBytes, values, result, n,
                         amrex::Gpu::gpuStream());
  cub::DeviceReduce::Sum(reduction_scratch(tmpBytes), tmpBytes, values, result,
                         n, amrex::Gpu::gpuStream());
#else
  *result = deterministic_sum(n, std::forward<F>(f));
#endif
}

// Deterministic replacement for amrex::MultiFab::Dot (whose OpenMP / GPU
// reduction order depends on the thread count): sum over the local FABs, in
// local index order, of sum_n x(i,j,k,xcomp+n) * y(i,j,k,ycomp+n) on the valid
// box grown by nghost, using deterministic_sum() inside each FAB. Cells where
// mask is 0 are skipped (mask may be nullptr). Unless local, the result is
// summed over all MPI ranks.
inline amrex::Real mf_dot(const amrex::MultiFab& x, int xcomp,
                          const amrex::MultiFab& y, int ycomp, int numcomp,
                          int nghost, const amrex::iMultiFab* mask = nullptr,
                          bool local = false) {
  amrex::Real sum = 0;
  for (amrex::MFIter mfi(x, false); mfi.isValid(); ++mfi) {
    const amrex::Box bx = amrex::grow(mfi.validbox(), nghost);
    const auto xa = x.const_array(mfi);
    const auto ya = y.const_array(mfi);
    const bool hasMask = (mask != nullptr);
    const amrex::Array4<const int> ma =
        hasMask ? mask->const_array(mfi) : amrex::Array4<const int>();
    const auto lo = amrex::lbound(bx);
    const int nx = bx.length(0);
    const int ny = (AMREX_SPACEDIM > 1) ? bx.length(1) : 1;
    const int nPts = static_cast<int>(bx.numPts());
    sum += deterministic_sum(nPts, [=] AMREX_GPU_HOST_DEVICE(int idx) {
      const int i = lo.x + idx % nx;
      const int j = lo.y + (idx / nx) % ny;
      const int k = lo.z + idx / (nx * ny);
      if (hasMask && !ma(i, j, k))
        return amrex::Real(0);
      amrex::Real t = 0;
      for (int n = 0; n < numcomp; ++n)
        t += xa(i, j, k, xcomp + n) * ya(i, j, k, ycomp + n);
      return t;
    });
  }
  if (!local)
    amrex::ParallelAllReduce::Sum(sum,
                                  amrex::ParallelContext::CommunicatorSub());
  return sum;
}

// Several deterministic_sum() results returned to the host together: add()
// starts the sum of f(i), i = 0..n-1, and values() returns the sums in the
// order of the add() calls (and forgets them). The sums are the same as
// those of deterministic_sum(), but a GPU build waits for the device once in
// values() instead of once per sum, e.g. one partial sum per box or tile.
//
// The sums share the CUB scratch, so they run on one stream: outside a
// StreamOrderedRegion (several streams) add() waits for each sum, as
// deterministic_sum() does.
class DeferredSums {
public:
  DeferredSums() = default;
  DeferredSums(const DeferredSums&) = delete;
  DeferredSums& operator=(const DeferredSums&) = delete;
  ~DeferredSums() {
#ifdef AMREX_USE_CUDA
    if (buf) {
      if (nPending > 0)
        amrex::Gpu::streamSynchronize();
      amrex::The_Pinned_Arena()->free(buf);
    }
#endif
  }

  template <class F> void add(int n, F&& f) {
#ifdef AMREX_USE_CUDA
    if (amrex::Gpu::numGpuStreams() > 1) {
      flush();
      done.push_back(deterministic_sum(n, std::forward<F>(f)));
      return;
    }
    if (!buf)
      buf = static_cast<amrex::Real*>(
          amrex::The_Pinned_Arena()->alloc(sizeof(amrex::Real) * capacity));
    if (nPending == capacity)
      flush();
    // The result goes straight to pinned host memory.
    deterministic_sum_to(n, std::forward<F>(f), buf + nPending++);
#else
    done.push_back(deterministic_sum(n, std::forward<F>(f)));
#endif
  }

  std::vector<amrex::Real> values() {
#ifdef AMREX_USE_CUDA
    flush();
#endif
    std::vector<amrex::Real> v;
    v.swap(done);
    return v;
  }

private:
#ifdef AMREX_USE_CUDA
  void flush() {
    if (nPending == 0)
      return;
    amrex::Gpu::streamSynchronize();
    done.insert(done.end(), buf, buf + nPending);
    nPending = 0;
  }

  static constexpr int capacity = 256;
  amrex::Real* buf = nullptr;
  int nPending = 0;
#endif
  std::vector<amrex::Real> done;
};

// Sum of the values f(i), i = 0..n-1, returned by the threads of one or more
// kernel launches (add() can be called many times before value()).
//
// A CPU build adds the values one by one, in launch and index order, which is
// exactly what a serial loop with `sum += ...` does; amrex::ReduceOps would
// add a partial sum per launch instead, which changes the round-off when the
// launches are small (e.g. one particle tile per cell). A GPU build uses
// deterministic_sum() for every launch, so the result
// is the same from run to run.
class KernelSum {
public:
  template <class F> void add(int n, F&& f) {
#ifdef AMREX_USE_GPU
    // One deterministic partial sum per launch, added in launch order.
    partials.add(n,
                 [=] AMREX_GPU_DEVICE(int i) { return amrex::get<0>(f(i)); });
#else
    for (int i = 0; i < n; ++i)
      sum += amrex::get<0>(f(i));
#endif
  }

  amrex::Real value() {
#ifdef AMREX_USE_GPU
    for (const amrex::Real v : partials.values())
      sum += v;
#endif
    return sum;
  }

private:
  amrex::Real sum = 0;
#ifdef AMREX_USE_GPU
  DeferredSums partials;
#endif
};

// Host access to FAB data for output code that reads it point by point. A
// CPU build returns the FAB's own data. A GPU build copies a FAB to pinned
// host memory the first time it is asked for and reuses the copy, so call
// clear() whenever the data may have changed (e.g. before every output).
class HostFabCache {
public:
  template <class FAB>
  amrex::Array4<typename FAB::value_type const> array(
      const amrex::FabArray<FAB>& fa, const amrex::MFIter& mfi) {
#ifdef AMREX_USE_GPU
    using T = typename FAB::value_type;
    const Key key{ static_cast<const void*>(&fa), mfi.index() };
    auto it = fabs.find(key);
    if (it == fabs.end()) {
      const FAB& src = fa[mfi];
      const std::size_t nBytes = src.nBytes();
      auto copy = std::make_shared<amrex::BaseFab<T> >(
          src.box(), src.nComp(), amrex::The_Pinned_Arena());
      amrex::Gpu::dtoh_memcpy(copy->dataPtr(), src.dataPtr(), nBytes);
      it = fabs.emplace(key, std::static_pointer_cast<void>(copy)).first;
    }
    return static_cast<const amrex::BaseFab<T>*>(it->second.get())
        ->const_array();
#else
    return fa[mfi].const_array();
#endif
  }

  void clear() { fabs.clear(); }

private:
  using Key = std::pair<const void*, int>;
  std::map<Key, std::shared_ptr<void> > fabs;
};

// Evaluate f(mfi, i, j, k, arr) on the host for every point of `box` of the
// FAB `mfi` of `fa`, where arr is an Array4 of the FAB, and store the result
// in the FAB.
//
// This is for code that needs host-only functions (virtual functions,
// Regions, std::function, ...). In a GPU build the FAB goes through a pinned
// host buffer, so keep it out of the time loop where possible.
template <class FAB, class F>
void host_fill_fab(amrex::FabArray<FAB>& fa, amrex::MFIter& mfi,
                   const amrex::Box& box, F&& f) {
#ifdef AMREX_USE_GPU
  using T = typename FAB::value_type;
  FAB& fab = fa[mfi];
  amrex::BaseFab<T> host(fab.box(), fab.nComp(), amrex::The_Pinned_Arena());
  // Start from the current values: f may leave some of them untouched.
  amrex::Gpu::dtoh_memcpy(host.dataPtr(), fab.dataPtr(), fab.nBytes());
  const auto arr = host.array();
  amrex::LoopOnCpu(box, [&](int i, int j, int k) { f(mfi, i, j, k, arr); });
  amrex::Gpu::htod_memcpy(fab.dataPtr(), host.dataPtr(), fab.nBytes());
#else
  const auto arr = fa[mfi].array();
  amrex::LoopOnCpu(box, [&](int i, int j, int k) { f(mfi, i, j, k, arr); });
#endif
}

// Make the FABs of `fa` host accessible for the lifetime of the object (for
// host-only code such as the user sources): a GPU build swaps the data of
// every FAB with a pinned host copy, and swaps it back (copying the possibly
// modified host values to the device) in the destructor or in restore(). A
// CPU build does nothing.
template <class FAB> class HostStaging {
public:
  explicit HostStaging(amrex::FabArray<FAB>& fa, bool writeBack = true)
      : fa_(fa), writeBack_(writeBack) {
#ifdef AMREX_USE_GPU
    for (amrex::MFIter mfi(fa_); mfi.isValid(); ++mfi) {
      FAB& fab = fa_[mfi];
      auto host = std::make_unique<FAB>(fab.box(), fab.nComp(),
                                        amrex::The_Pinned_Arena());
      amrex::Gpu::dtoh_memcpy(host->dataPtr(), fab.dataPtr(), fab.nBytes());
      swap_fabs(fab, *host); // fa_ now holds the host copy
      saved_.push_back(std::move(host));
    }
#endif
  }
  ~HostStaging() { restore(); }
  HostStaging(const HostStaging&) = delete;
  HostStaging& operator=(const HostStaging&) = delete;

  void restore() {
#ifdef AMREX_USE_GPU
    int i = 0;
    for (amrex::MFIter mfi(fa_); mfi.isValid() && !saved_.empty(); ++mfi) {
      FAB& fab = fa_[mfi];
      FAB& device = *saved_[i++];
      if (writeBack_)
        amrex::Gpu::htod_memcpy(device.dataPtr(), fab.dataPtr(), fab.nBytes());
      swap_fabs(fab, device);
    }
    saved_.clear();
#endif
  }

private:
  static void swap_fabs(FAB& a, FAB& b) {
    FAB tmp(std::move(a));
    a = std::move(b);
    b = std::move(tmp);
  }

  amrex::FabArray<FAB>& fa_;
  bool writeBack_;
  std::vector<std::unique_ptr<FAB> > saved_;
};

// Run f(particles, n) on the host for the n particles at `data` (device
// memory in a GPU build), e.g. for I/O or MPI. A GPU build copies the
// particles to a host buffer and, if writeBack, copies them back afterwards;
// a CPU build passes `data` itself.
template <class P, class F>
void host_particles(P* data, int n, bool writeBack, F&& f) {
  if (n <= 0)
    return;
#ifdef AMREX_USE_GPU
  amrex::Gpu::PinnedVector<P> host(n);
  amrex::Gpu::dtoh_memcpy(host.data(), data, sizeof(P) * n);
  f(host.data(), n);
  if (writeBack)
    amrex::Gpu::htod_memcpy(data, host.data(), sizeof(P) * n);
#else
  amrex::ignore_unused(writeBack);
  f(data, n);
#endif
}

// host_fill_fab() for every FAB of `fa`, on the valid box or, if withGhost,
// on the whole FAB including its ghost cells.
template <class FAB, class F>
void host_fill(amrex::FabArray<FAB>& fa, bool withGhost, F&& f) {
  for (amrex::MFIter mfi(fa); mfi.isValid(); ++mfi) {
    host_fill_fab(fa, mfi, withGhost ? mfi.fabbox() : mfi.validbox(), f);
  }
}

#endif
