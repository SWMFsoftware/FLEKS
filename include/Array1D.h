#ifndef _Arr1D_H_
#define _Arr1D_H_

#include <type_traits>

#include <AMReX_BaseFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_REAL.H>

// The default constructor leaves the values uninitialized, so that the type
// is trivial: AMReX then allocates FABs of it without a constructor kernel
// per FAB (e.g. the temporary FABs of every SumBoundary). Set the values
// before use (FabArray::setVal, or Arr1D(b)).
template <class T, const int n> struct Arr1D {
  T data[n];

  Arr1D() = default;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Arr1D(const T& b) {
    for (int i = 0; i < n; ++i)
      data[i] = b;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Arr1D& operator=(const T& b) {
    for (int i = 0; i < n; ++i)
      data[i] = b;
    return *this;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& operator[](const int i) {
    return data[i];
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& operator[](
      const int i) const {
    return data[i];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Arr1D& operator+=(const Arr1D& b) {
    for (int i = 0; i < n; ++i)
      data[i] += b.data[i];
    return *this;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Arr1D& operator*=(const Arr1D& b) {
    for (int i = 0; i < n; ++i)
      data[i] *= b.data[i];
    return *this;
  }

  template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U> > >
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Arr1D& operator*=(U b) {
    for (int i = 0; i < n; ++i)
      data[i] *= static_cast<T>(b);
    return *this;
  }
};

template <class T, const int n>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Arr1D<T, n> operator+(
    Arr1D<T, n> a, const Arr1D<T, n>& b) {
  a += b;
  return a;
}

template <class T, const int n, typename U,
          typename = std::enable_if_t<std::is_arithmetic_v<U> > >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Arr1D<T, n> operator*(Arr1D<T, n> a,
                                                               U b) {
  a *= b;
  return a;
}

template <class T, const int n, typename U,
          typename = std::enable_if_t<std::is_arithmetic_v<U> > >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Arr1D<T, n> operator*(U b,
                                                               Arr1D<T, n> a) {
  a *= b;
  return a;
}

#if (AMREX_SPACEDIM == 2)
constexpr int nMMComponents = 81;
constexpr int nCMMComponents = 9;
#else
constexpr int nMMComponents = 243;
constexpr int nCMMComponents = 27;
#endif

using RealMM = Arr1D<amrex::Real, nMMComponents>;
using RealCMM = Arr1D<amrex::Real, nCMMComponents>;

using NodeMMFab = amrex::MultiFab;
using CenterMMFab = amrex::FabArray<amrex::BaseFab<RealCMM> >;

static_assert(std::is_trivially_default_constructible_v<RealMM>&&
                  std::is_trivially_copyable_v<RealMM>,
              "RealMM must stay trivial (see Arr1D)");

#endif
