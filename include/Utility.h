#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <unistd.h>

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_iMultiFab.H>

#include "Constants.h"
#include "Timer.h"

// Only works for x>-8;
inline int fastfloor(amrex::Real x) { return (int)(x + 8) - 8; }

inline int product(const amrex::IntVect& vect) {
  return AMREX_D_TERM(vect[0], *vect[1], *vect[2]);
}

inline amrex::Dim3 init_dim3(const int i) {
  amrex::Dim3 dim;
  dim.x = i;
  dim.y = i;
  dim.z = nDim > 2 ? i : 0;
  return dim;
}

// Return the firs integer in the string.
// Example: "abc123def345b" -> 123
inline int extract_int(std::string s) {
  std::size_t i0 = s.find_first_of("0123456789");
  std::size_t i1 = s.find_last_of("0123456789");
  return std::stoi(s.substr(i0, i1 - i0 + 1));
};

inline int shift_periodic_index(int idx, int lo, int hi) {
  if (idx > hi)
    idx -= hi - lo;

  if (idx < lo)
    idx += hi - lo;

  return idx;
};

template <class T> inline void a_cross_b(T (&a)[3], T (&b)[3], T (&c)[3]) {
  // c = a x b

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

  return;
}

template <class T> inline void zero_array(T* arr, int nSize) {
  for (int i = 0; i < nSize; i++)
    arr[i] = 0;
}

// rand1, rand2 are random numbers in [0,1].
inline void random_vector(amrex::Real rand1, amrex::Real rand2,
                          amrex::Real (&vec)[nDim3]) {
  amrex::Real lon = 2 * M_PI * rand1;
  amrex::Real costheta = 2 * rand2 - 1;
  amrex::Real sintheta = sqrt(1 - costheta * costheta);

  vec[ix_] = sintheta * cos(lon);
  vec[iy_] = sintheta * sin(lon);
  vec[iz_] = costheta;
}

inline amrex::Real l2_norm(amrex::Real* vec, int n) {
  amrex::Real sum = 0;
  for (int i = 0; i < n; i++)
    sum += vec[i] * vec[i];
  return sqrt(sum);
}

inline void linear_interpolation_coef(amrex::RealVect& dx,
                                      amrex::Real (&coef)[2][2][2]) {

  amrex::Real xi[2];
  amrex::Real eta[2];
  amrex::Real zeta[2];
  xi[0] = dx[0];
  xi[1] = 1 - xi[0];

  eta[0] = dx[1];
  eta[1] = 1 - eta[0];

  zeta[0] = nDim > 2 ? dx[2] : 0;
  zeta[1] = 1 - zeta[0];

  amrex::Real multi[2][2];
  multi[0][0] = xi[0] * eta[0];
  multi[0][1] = xi[0] * eta[1];
  multi[1][0] = xi[1] * eta[0];
  multi[1][1] = xi[1] * eta[1];

  coef[0][0][0] = multi[1][1] * zeta[1];
  coef[0][0][1] = multi[1][1] * zeta[0];
  coef[0][1][0] = multi[1][0] * zeta[1];
  coef[0][1][1] = multi[1][0] * zeta[0];
  coef[1][0][0] = multi[0][1] * zeta[1];
  coef[1][0][1] = multi[0][1] * zeta[0];
  coef[1][1][0] = multi[0][0] * zeta[1];
  coef[1][1][1] = multi[0][0] * zeta[0];
}

template <typename T>
inline T bound(const T& val, const T& xmin, const T& xmax) {
  if (val < xmin)
    return xmin;
  if (val > xmax)
    return xmax;
  return val;
}

template <typename T, int nRow, int nCol>
bool linear_solver_Gauss_Elimination(
    int m, int n, amrex::Array2D<T, 0, nRow - 1, 0, nCol - 1>& a,
    amrex::Vector<T>& x, amrex::Vector<T>& ref) {
  // a(m,n)
  if (m > nRow) {
    amrex::Abort("Error: m>nRow in linear_solver_Gauss_Elimination");
  }
  if (n > nCol) {
    amrex::Abort("Error: n>nCol in linear_solver_Gauss_Elimination");
  }

  for (int i = 0; i < m - 1; i++) {
    // Partial Pivoting
    for (int k = i + 1; k < m; k++) {
      // If diagonal element(absolute vallue) is smaller than
      // any of the terms below it
      if (fabs(a(i, i)) < fabs(a(k, i))) {
        // Swap the rows
        for (int j = 0; j < n; j++) {
          std::swap(a(i, j), a(k, j));
        }
      }
    }
    // Begin Gauss Elimination
    for (int k = i + 1; k < m; k++) {
      if (fabs(a(i, i)) < ref[i]) {
        return false;
      }
      double term = a(k, i) / a(i, i);
      for (int j = 0; j < n; j++) {
        a(k, j) = a(k, j) - term * a(i, j);
      }
    }
  }
  // Begin Back-substitution
  for (int i = m - 1; i >= 0; i--) {
    x[i] = a(i, n - 1);
    for (int j = i + 1; j < n - 1; j++) {
      x[i] = x[i] - a(i, j) * x[j];
    }

    if (fabs(a(i, i)) < ref[i]) {
      return false;
    }

    x[i] = x[i] / a(i, i);
  }
  return true;
};

// Fisher-Yates shuffle algorithm
template <class Vec, class Rand> void shuffle_fish_yates(Vec& arr, Rand& rd) {
  for (int i = arr.size() - 1; i > 0; --i) {
    int j = floor(rd() * i);
    if (i != j)
      std::swap(arr[i], arr[j]);
  }
};

// Randomly select m elements from an array with weights. These m elements can
// not repeat.
template <class Vec, class Rand>
std::vector<int> random_select_weighted_n(Vec arr, int m, Rand& rd) {

  const int n = arr.size();
  std::vector<int> idx;

  if (m >= n) {
    for (int i = 0; i < n; i++)
      idx.push_back(i);
    return idx;
  }

  amrex::Vector<double> w1;
  w1.resize(n);

  double wt = 0;
  for (int i = 0; i < arr.size(); i++) {
    wt += arr[i];
    w1[i] = wt;
  }

  const int iSelected = 1, iNotSelected = 0;
  amrex::Vector<int> selected(n, iNotSelected);

  while (idx.size() < size_t(m)) {
    double r = rd() * wt;
    int i = 0;
    while (w1[i] < r)
      i++;

    if (selected[i] == iNotSelected) {
      idx.push_back(i);
      selected[i] = iSelected;

      { // 'Remove' arr[i] so that it will not be selected again.
        arr[i] = 0;
        wt = 0;
        for (int j = 0; j < arr.size(); j++) {
          wt += arr[j];
          w1[j] = wt;
        }
      }
    }
  }

  return idx;
};

inline double read_mem_usage() {
  // This function returns the resident set size (RSS) of
  // this processor in unit MB.

  // From wiki:
  // RSS is the portion of memory occupied by a process that is
  // held in main memory (RAM).

  double rssMB = 0.0;

  std::ifstream stat_stream("/proc/self/stat", std::ios_base::in);

  if (!stat_stream.fail()) {
    // Dummy vars for leading entries in stat that we don't care about
    std::string pid, comm, state, ppid, pgrp, session, tty_nr;
    std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    std::string utime, stime, cutime, cstime, priority, nice;
    std::string O, itrealvalue, starttime;

    // Two values we want
    unsigned long vsize;
    unsigned long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >>
        tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >>
        stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >>
        starttime >> vsize >> rss; // Ignore the rest
    stat_stream.close();

    rssMB = rss * sysconf(_SC_PAGE_SIZE) / 1024.0 / 1024.0;
  }

  return rssMB;
}

#endif
