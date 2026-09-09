#ifndef _Arr1D_H_
#define _Arr1D_H_

#include <type_traits>

#include <AMReX_BaseFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_REAL.H>

template <class T, const int n> struct Arr1D {
  T data[n];

  Arr1D(const T& b = T(0)) {
    for (int i = 0; i < n; ++i)
      data[i] = b;
  }

  Arr1D& operator=(const T& b) {
    for (int i = 0; i < n; ++i)
      data[i] = b;
    return *this;
  }

  T& operator[](const int i) { return data[i]; }
  const T& operator[](const int i) const { return data[i]; }

  Arr1D& operator+=(const Arr1D& b) {
    for (int i = 0; i < n; ++i)
      data[i] += b.data[i];
    return *this;
  }

  Arr1D& operator*=(const Arr1D& b) {
    for (int i = 0; i < n; ++i)
      data[i] *= b.data[i];
    return *this;
  }

  template <typename U, typename = std::enable_if_t<std::is_arithmetic_v<U> > >
  Arr1D& operator*=(U b) {
    for (int i = 0; i < n; ++i)
      data[i] *= static_cast<T>(b);
    return *this;
  }
};

template <class T, const int n>
inline Arr1D<T, n> operator+(Arr1D<T, n> a, const Arr1D<T, n>& b) {
  a += b;
  return a;
}

template <class T, const int n, typename U,
          typename = std::enable_if_t<std::is_arithmetic_v<U> > >
inline Arr1D<T, n> operator*(Arr1D<T, n> a, U b) {
  a *= b;
  return a;
}

template <class T, const int n, typename U,
          typename = std::enable_if_t<std::is_arithmetic_v<U> > >
inline Arr1D<T, n> operator*(U b, Arr1D<T, n> a) {
  a *= b;
  return a;
}

using RealMM = Arr1D<amrex::Real, 243>;
using RealCMM = Arr1D<amrex::Real, 27>;

using NodeMMFab = amrex::FabArray<amrex::BaseFab<RealMM> >;
using CenterMMFab = amrex::FabArray<amrex::BaseFab<RealCMM> >;

#endif
