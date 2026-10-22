#ifndef _GPUBACKEND_H_
#define _GPUBACKEND_H_

// Helpers that let the same FLEKS kernels run on CPUs and GPUs.
//
// FLEKS keeps one code path: kernels are AMREX_GPU_HOST_DEVICE lambdas
// launched by amrex::ParallelFor, which compiles to a plain loop in a CPU
// build. Code that must differ between the two builds lives in this file
// only (tools/check_gpu_kernels.py enforces it), so the physics never forks.

#include <algorithm>

#include <AMReX_Extension.H>
#include <AMReX_FabArray.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Loop.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>

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

// Sum of the values f(i), i = 0..n-1, returned by the threads of one or more
// kernel launches (add() can be called many times before value()).
//
// A CPU build adds the values one by one, in launch and index order, which is
// exactly what a serial loop with `sum += ...` does; amrex::ReduceOps would
// add a partial sum per launch instead, which changes the round-off when the
// launches are small (e.g. one particle tile per cell). A GPU build uses
// amrex::ReduceOps, which combines the values in a fixed order, so the result
// is the same from run to run.
class KernelSum {
public:
  template <class F> void add(int n, F&& f) {
#ifdef AMREX_USE_GPU
    reduceOp.eval(n, reduceData, std::forward<F>(f));
#else
    for (int i = 0; i < n; ++i)
      sum += amrex::get<0>(f(i));
#endif
  }

  amrex::Real value() {
#ifdef AMREX_USE_GPU
    return amrex::get<0>(reduceData.value(reduceOp));
#else
    return sum;
#endif
  }

private:
#ifdef AMREX_USE_GPU
  amrex::ReduceOps<amrex::ReduceOpSum> reduceOp;
  amrex::ReduceData<amrex::Real> reduceData{ reduceOp };
#else
  amrex::Real sum = 0;
#endif
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
  FAB host(fa[mfi].box(), fa.nComp(), amrex::The_Pinned_Arena());
  // Start from the current values: f may leave some of them untouched.
  host.template copy<amrex::RunOn::Device>(fa[mfi]);
  amrex::Gpu::streamSynchronize();
  const auto arr = host.array();
  amrex::LoopOnCpu(box, [&](int i, int j, int k) { f(mfi, i, j, k, arr); });
  fa[mfi].template copy<amrex::RunOn::Device>(host);
  amrex::Gpu::streamSynchronize();
#else
  const auto arr = fa[mfi].array();
  amrex::LoopOnCpu(box, [&](int i, int j, int k) { f(mfi, i, j, k, arr); });
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
