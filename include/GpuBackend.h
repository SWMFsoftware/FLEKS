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
constexpr bool fleksUseGpu = true;
#else
constexpr bool fleksUseGpu = false;
#endif

// Particle tile layout. The CPU keeps one particle tile per cell: the
// particles of a cell are contiguous (cache locality) and the per-cell
// algorithms (resampling, sources) work on a tile. The GPU uses one tile per
// box, because every tile costs a kernel launch.
constexpr bool useCellTiles = !fleksUseGpu;

// Deposition strategy. The CPU adds every particle in order (deposit_add is a
// plain +=). The GPU sums the particles of each cell in a fixed order and
// gathers per node (CellBins, deposit_to_nodes), which is reproducible run to
// run without floating-point atomics.
constexpr bool useCellBinDeposits = fleksUseGpu;

// Capacity of a KernelArray. FLEKS runs use a handful of species or fluids.
constexpr int nKernelArrayMax = 16;

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

// The particles of a tile grouped by cell ("bin"). perm()[i] is the i-th
// particle in bin order, and the particles of bin b are
// perm()[offsets()[b]] ... perm()[offsets()[b + 1] - 1]. A particle whose bin
// is nBins (or larger) is left out.
//
// Inside a bin the particles keep their storage order (build(np, nBins,
// binOf)), or are ordered by a unique particle key (build(np, nBins, binOf,
// idOf), e.g. built from the particle cpu and id). The storage order is the
// same from run to run except right after AMReX's GPU redistribution, which
// Particles::sort_particles() follows with the key order. So kernels that loop
// over the particles of a bin in this order give bitwise-reproducible results
// without floating-point atomics.
class CellBins {
public:
  // binOf(ip) returns the bin of particle ip (a device lambda); a bin keeps
  // the storage order.
  template <class F> void build(int np, int nBins, F&& binOf) {
    resize(np, nBins);
    int* const v = index.data();
#ifdef AMREX_USE_CUDA
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int i) { v[i] = i; });
#endif
    sort_by_bin(np, nBins, binOf, v);
  }

  // binOf(ip) returns the bin of particle ip and idOf(ip) its unique key (an
  // unsigned 64-bit integer); both must be device lambdas. A bin is ordered
  // by the key.
  template <class F, class G>
  void build(int np, int nBins, F&& binOf, G&& idOf) {
    resize(np, nBins);
#ifdef AMREX_USE_CUDA
    int* const v = index.data();
    // 1. Order the particles by their unique key.
    amrex::Gpu::DeviceVector<unsigned long long> ids(np), sortedIds(np);
    amrex::Gpu::DeviceVector<int> byIdVec(np);
    unsigned long long* const id = ids.data();
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int i) {
      id[i] = idOf(i);
      v[i] = i;
    });
    std::size_t tmpBytes = 0;
    cub::DeviceRadixSort::SortPairs(nullptr, tmpBytes, id, sortedIds.data(), v,
                                    byIdVec.data(), np, 0, 64,
                                    amrex::Gpu::gpuStream());
    amrex::Gpu::DeviceVector<char> tmp(tmpBytes);
    cub::DeviceRadixSort::SortPairs(tmp.data(), tmpBytes, id, sortedIds.data(),
                                    v, byIdVec.data(), np, 0, 64,
                                    amrex::Gpu::gpuStream());
    // 2. A stable sort by bin keeps the key order in a bin.
    sort_by_bin(np, nBins, binOf, byIdVec.data());
#else
    amrex::ignore_unused(idOf);
    sort_by_bin(np, nBins, binOf, index.data());
#endif
  }

  // One bin per cell of cellBox (bin = x-fastest index of the cell);
  // cellOf(ip, cell) gives the cell of particle ip, or false to leave the
  // particle out, as does a cell outside cellBox. A bin keeps the storage
  // order.
  template <class CellOf>
  void build_cells(int np, const amrex::Box& cellBox, CellOf&& cellOf) {
    const amrex::Dim3 lo = amrex::lbound(cellBox);
    const amrex::Dim3 len = amrex::length(cellBox);
    const int nBins = static_cast<int>(cellBox.numPts());
    build(np, nBins, [=] AMREX_GPU_DEVICE(int ip) {
      amrex::IntVect c;
      if (!cellOf(ip, c))
        return nBins;
      const int ijk[3] = { AMREX_D_DECL(c[0] - lo.x, c[1] - lo.y,
                                        c[2] - lo.z) };
      const int lenD[3] = { len.x, len.y, len.z };
      for (int d = 0; d < AMREX_SPACEDIM; ++d)
        if (ijk[d] < 0 || ijk[d] >= lenD[d])
          return nBins;
      return ijk[0] + len.x * ((AMREX_SPACEDIM > 1 ? ijk[1] : 0) +
                               len.y * (AMREX_SPACEDIM > 2 ? ijk[2] : 0));
    });
  }

  const int* perm() const { return sorted.data(); }
  const int* offsets() const { return offset.data(); }
  int n_bins() const { return nBin; }

  // Implementation (public only because nvcc does not allow device lambdas in
  // private member functions).
  void resize(int np, int nBins) {
    nBin = nBins;
    keys.resize(np);
    sortedKeys.resize(np);
    index.resize(np);
    sorted.resize(np);
    offset.resize(nBins + 2);
  }

  // Stable sort of the particles order[0..np) by bin into sorted, and the
  // start of every bin.
  template <class F>
  void sort_by_bin(int np, int nBins, F&& binOf, const int* order) {
    int* const k = keys.data();
    int* const v = index.data();
#ifdef AMREX_USE_CUDA
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int i) {
      const int p = order[i];
      const int b = binOf(p);
      k[i] = (b < 0 || b > nBins) ? nBins : b;
      v[i] = p;
    });
    int nBits = 1;
    while ((1L << nBits) <= nBins)
      ++nBits;
    std::size_t tmpBytes = 0;
    cub::DeviceRadixSort::SortPairs(nullptr, tmpBytes, k, sortedKeys.data(), v,
                                    sorted.data(), np, 0, nBits,
                                    amrex::Gpu::gpuStream());
    amrex::Gpu::DeviceVector<char> tmp(tmpBytes);
    cub::DeviceRadixSort::SortPairs(tmp.data(), tmpBytes, k, sortedKeys.data(),
                                    v, sorted.data(), np, 0, nBits,
                                    amrex::Gpu::gpuStream());
    // Start of every bin: binary search in the sorted keys.
    const int* const sk = sortedKeys.data();
    int* const off = offset.data();
    amrex::ParallelFor(nBins + 2, [=] AMREX_GPU_DEVICE(int b) {
      int lo = 0, hi = np;
      while (lo < hi) {
        const int mid = (lo + hi) / 2;
        if (sk[mid] < b)
          lo = mid + 1;
        else
          hi = mid;
      }
      off[b] = lo;
    });
    amrex::Gpu::streamSynchronize();
#else
    // Stable counting sort in storage order: a CPU run keeps its order.
    amrex::ignore_unused(order, v);
    for (int i = 0; i < np; ++i) {
      const int b = binOf(i);
      k[i] = (b < 0 || b > nBins) ? nBins : b;
    }
    std::vector<int> count(nBins + 2, 0);
    for (int i = 0; i < np; ++i)
      count[k[i] + 1]++;
    for (int b = 0; b <= nBins; ++b)
      count[b + 1] += count[b];
    for (int b = 0; b < nBins + 2; ++b)
      offset[b] = count[b];
    for (int i = 0; i < np; ++i)
      sorted[count[k[i]]++] = i;
#endif
  }

private:
  int nBin = 0;
  amrex::Gpu::DeviceVector<int> keys, sortedKeys, index, sorted, offset;
};

// Deposit particle quantities onto the S^D nodes (S = 2 by default) starting
// at the cell of each particle:
//   cellOf(ip, cell) -> bool   the cell (its lower node) of particle ip, or
//                              false to skip the particle;
//   contrib(ip, acc)           add the contribution of particle ip to
//                              acc[ii][jj][kk][n], n < NComp, the nodes
//                              cell + (ii, jj, kk), ii, jj, kk < S;
// dst(node, comp0 + n) receives the sums. cellBox must contain the cells of
// all the particles, and dst the nodes of cellBox.
//
// A CPU build adds every particle in order, exactly like a serial loop. A GPU
// build sums the particles of each cell in storage order (see CellBins) and
// then every node sums its 2^D cells in a fixed order, so the result is the
// same from run to run, without floating-point atomics.
//
// deposit_to_nodes_with() is the general form: add(i, j, k, n, value) stores
// the sums (e.g. into a FAB of a non-scalar type).
template <int NComp, int S = 2, class Add, class CellOf, class Contrib>
void deposit_to_nodes_with(int np, const amrex::Box& cellBox, Add&& add,
                           CellOf&& cellOf, Contrib&& contrib) {
  using amrex::Real;
  constexpr int nK = AMREX_SPACEDIM > 2 ? S : 1;
  if (np == 0)
    return;
#ifdef AMREX_USE_GPU
  const amrex::Dim3 lo = amrex::lbound(cellBox);
  const amrex::Dim3 len = amrex::length(cellBox);
  const int nBins = static_cast<int>(cellBox.numPts());
  CellBins bins;
  bins.build_cells(np, cellBox, cellOf);

  // 1. The sum of every cell, particle by particle in the bin order.
  constexpr int nAcc = S * S * nK * NComp;
  amrex::Gpu::DeviceVector<Real> scratch(static_cast<std::size_t>(nBins) *
                                         nAcc);
  Real* const sc = scratch.data();
  const int* const perm = bins.perm();
  const int* const off = bins.offsets();
  amrex::ParallelFor(nBins, [=] AMREX_GPU_DEVICE(int b) {
    Real acc[S][S][S][NComp];
    for (int kk = 0; kk < S; ++kk)
      for (int jj = 0; jj < S; ++jj)
        for (int ii = 0; ii < S; ++ii)
          for (int n = 0; n < NComp; ++n)
            acc[ii][jj][kk][n] = 0;
    for (int i = off[b]; i < off[b + 1]; ++i)
      contrib(perm[i], acc);
    Real* const out = sc + static_cast<std::size_t>(b) * nAcc;
    for (int kk = 0; kk < nK; ++kk)
      for (int jj = 0; jj < S; ++jj)
        for (int ii = 0; ii < S; ++ii)
          for (int n = 0; n < NComp; ++n)
            out[((kk * S + jj) * S + ii) * NComp + n] = acc[ii][jj][kk][n];
  });

  // 2. Every node sums its cells.
  amrex::Box nodeBox = amrex::surroundingNodes(cellBox);
  for (int d = 0; d < AMREX_SPACEDIM; ++d)
    nodeBox.growHi(d, S - 2);
  amrex::ParallelFor(nodeBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    const int ijk[3] = { i, j, k };
    Real sum[NComp];
    for (int n = 0; n < NComp; ++n)
      sum[n] = 0;
    for (int kk = 0; kk < nK; ++kk)
      for (int jj = 0; jj < S; ++jj)
        for (int ii = 0; ii < S; ++ii) {
          const int corner[3] = { ii, jj, kk };
          int c[3] = { 0, 0, 0 };
          bool inside = true;
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int l = d == 0 ? len.x : (d == 1 ? len.y : len.z);
            const int lo_d = d == 0 ? lo.x : (d == 1 ? lo.y : lo.z);
            c[d] = ijk[d] - corner[d] - lo_d;
            inside = inside && c[d] >= 0 && c[d] < l;
          }
          if (!inside)
            continue;
          const int b = c[0] + len.x * (c[1] + len.y * c[2]);
          const Real* const in = sc + static_cast<std::size_t>(b) * nAcc +
                                 ((kk * S + jj) * S + ii) * NComp;
          for (int n = 0; n < NComp; ++n)
            sum[n] += in[n];
        }
    for (int n = 0; n < NComp; ++n)
      add(i, j, k, n, sum[n]);
  });
  // The kernels read the bins and the scratch, which are freed on return.
  amrex::Gpu::streamSynchronize();
#else
  amrex::ignore_unused(cellBox);
  for (int ip = 0; ip < np; ++ip) {
    amrex::IntVect c;
    if (!cellOf(ip, c))
      continue;
    Real acc[S][S][S][NComp];
    for (int kk = 0; kk < S; ++kk)
      for (int jj = 0; jj < S; ++jj)
        for (int ii = 0; ii < S; ++ii)
          for (int n = 0; n < NComp; ++n)
            acc[ii][jj][kk][n] = 0;
    contrib(ip, acc);
    for (int n = 0; n < NComp; ++n)
      for (int kk = 0; kk < nK; ++kk)
        for (int jj = 0; jj < S; ++jj)
          for (int ii = 0; ii < S; ++ii)
            add(c[0] + ii, AMREX_SPACEDIM > 1 ? c[1] + jj : 0,
                AMREX_SPACEDIM > 2 ? c[AMREX_SPACEDIM > 2 ? 2 : 0] + kk : 0, n,
                acc[ii][jj][kk][n]);
  }
#endif
}

template <int NComp, class CellOf, class Contrib>
void deposit_to_nodes(int np, const amrex::Box& cellBox,
                      const amrex::Array4<amrex::Real>& dst, int comp0,
                      CellOf&& cellOf, Contrib&& contrib) {
  deposit_to_nodes_with<NComp>(
      np, cellBox,
      [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k, int n, amrex::Real v) {
        dst(i, j, k, comp0 + n) += v;
      },
      cellOf, contrib);
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
  amrex::Gpu::DeviceVector<char> tmp(tmpBytes);
  cub::DeviceReduce::Sum(tmp.data(), tmpBytes, values, result, n,
                         amrex::Gpu::gpuStream());
  amrex::Gpu::streamSynchronize();
  return *result;
#else
  amrex::Real sum = 0;
  for (int i = 0; i < n; ++i)
    sum += f(i);
  return sum;
#endif
}

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
    sum += deterministic_sum(
        n, [=] AMREX_GPU_DEVICE(int i) { return amrex::get<0>(f(i)); });
#else
    for (int i = 0; i < n; ++i)
      sum += amrex::get<0>(f(i));
#endif
  }

  amrex::Real value() { return sum; }

private:
  amrex::Real sum = 0;
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
