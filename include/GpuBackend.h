#ifndef _GPUBACKEND_H_
#define _GPUBACKEND_H_

// Helpers that let the same FLEKS kernels run on CPUs and GPUs.
//
// FLEKS keeps one code path: kernels are AMREX_GPU_HOST_DEVICE lambdas
// launched by amrex::ParallelFor, which compiles to a plain loop in a CPU
// build. Code that must differ between the two builds lives in this file
// only (tools/check_gpu_kernels.py enforces it), so the physics never forks.

#include <algorithm>
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
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>

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
    amrex::Gpu::streamSynchronize();
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
// A pinned host double that receives the result of deterministic_sum().
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
  amrex::Real sum = 0;
  for (int i = 0; i < n; ++i)
    sum += f(i);
  *result = sum;
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
  amrex::Real sum = 0;
  for (int i = 0; i < n; ++i)
    sum += f(i);
  return sum;
#endif
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
    amrex::Real sum = 0;
    for (int i = 0; i < n; ++i)
      sum += f(i);
    done.push_back(sum);
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
