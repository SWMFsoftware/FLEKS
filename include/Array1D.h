#ifndef _ARRAY1D_H_
#define _ARRAY1D_H_

#include <iostream>

#include <AMReX_REAL.H>

template <class T, const int n> class Array1D {
public:
  T data[n];

  Array1D(){};
  Array1D(const T& b) {
    for (int i = 0; i < n; i++)
      data[i] = b;
  }

  Array1D(const Array1D<T, n>& b) {
    for (int i = 0; i < n; i++)
      data[i] = b.data[i];
  }

  Array1D<T, n>& operator=(const Array1D<T, n>& b) {
    for (int i = 0; i < n; i++)
      data[i] = b.data[i];
    return *this;
  }

  Array1D<T, n>& operator=(const T& b) {
    for (int i = 0; i < n; i++)
      data[i] = b;
    return *this;
  }

  Array1D<T, n>& operator+=(const Array1D<T, n>& b) {
    for (int i = 0; i < n; i++)
      data[i] += b.data[i];
    return *this;
  }

  Array1D<T, n>& operator<<(const Array1D<T, n>& b) {
    for (int i = 0; i < n; i++)
      data[i] += b.data[i];
    return *this;
  }
};

template <class T, const int n>
std::ostream& operator<<(std::ostream& out, const Array1D<T, n>& b) {
  for (int i = 0; i < n; i++)
    out << " i = " << i << " data[i] = " << b.data[i] << "\t";
  out << std::endl;
  return out;
}

using RealMM = Array1D<amrex::Real, 243>;
using RealCMM = Array1D<amrex::Real, 27>;

#endif
