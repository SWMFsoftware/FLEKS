#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_iMultiFab.H>

#include "Constants.h"
#include "Timer.h"

// Only works for x>-8;
inline int fastfloor(amrex::Real x) { return (int)(x + 8) - 8; }

void curl_center_to_node(const amrex::MultiFab& centerMF,
                         amrex::MultiFab& nodeMF, const amrex::Real* invDx);

void curl_node_to_center(const amrex::MultiFab& nodeMF,
                         amrex::MultiFab& centerMF, const amrex::Real* invDx);

void lap_node_to_node(const amrex::MultiFab& srcMF, amrex::MultiFab& dstMF,
                      const amrex::DistributionMapping dm,
                      const amrex::Geometry& gm,
                      const amrex::iMultiFab& status);

void grad_node_to_center(const amrex::MultiFab& nodeMF,
                         amrex::MultiFab& centerMF, const amrex::Real* invDx,
                         const amrex::iMultiFab& status);

void grad_center_to_node(const amrex::MultiFab& centerMF,
                         amrex::MultiFab& nodeMF, const amrex::Real* invDx);

void div_center_to_node(const amrex::MultiFab& centerMF,
                        amrex::MultiFab& nodeMF, const amrex::Real* invDx);

void div_node_to_center(const amrex::MultiFab& nodeMF,
                        amrex::MultiFab& centerMF, const amrex::Real* invDx);

void div_center_to_center(const amrex::MultiFab& srcMF, amrex::MultiFab& dstMF,
                          const amrex::Real* invDx);

void average_center_to_node(const amrex::MultiFab& centerMF,
                            amrex::MultiFab& nodeMF);

void print_MultiFab(const amrex::iMultiFab& data, std::string tag,
                    int nshift = 0);

void print_MultiFab(const amrex::MultiFab& data, std::string tag,
                    const int iVarStart, const int iVarEnd, int nshift = 0);

void print_MultiFab(const amrex::MultiFab& data, std::string tag,
                    amrex::Geometry& gm, int nshift = 0);

template <class FAB>
void print_fab(const amrex::FabArray<FAB>& mf, std::string tag,
               const int iStart, const int nComp, int nshift = 0) {
  amrex::AllPrint() << "-----" << tag << " begin-----" << std::endl;
  amrex::Real sum = 0;
  amrex::Real sum2 = 0;

  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const FAB& fab = mf[mfi];
    const auto& box = mfi.validbox();
    const auto& data = fab.array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    amrex::AllPrint() << "------ box = " << box << std::endl;
    for (int i = lo.x - nshift; i <= hi.x + nshift; ++i)
      for (int j = lo.y - nshift; j <= hi.y + nshift; ++j)
        for (int k = lo.z; k <= hi.z; ++k)
          for (int iVar = iStart; iVar < iStart + nComp; iVar++) {
            amrex::AllPrint() << " i = " << i << " j = " << j << " k = " << k
                              << " iVar = " << iVar
                              << " data = " << data(i, j, k, iVar) << std::endl;
            sum += data(i, j, k, iVar);
            sum2 += pow(data(i, j, k, iVar), 2);
          }
  }
  amrex::AllPrint() << "sum = " << sum << " sum2 = " << sqrt(sum2)
                    << " on proc = " << amrex::ParallelDescriptor::MyProc()
                    << std::endl;
  amrex::AllPrint() << "-----" << tag << " end-----" << std::endl;
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

inline int get_local_node_or_cell_number(const amrex::MultiFab& MF) {

  int nTotal = 0;
  if (!MF.empty())
    for (amrex::MFIter mfi(MF); mfi.isValid(); ++mfi) {
      const amrex::Box& box = mfi.validbox();
      const auto lo = lbound(box);
      const auto hi = ubound(box);
      nTotal += (hi.x - lo.x + 1) * (hi.y - lo.y + 1) * (hi.z - lo.z + 1);
    }
  return nTotal;
}

template <class T> inline void zero_array(T* arr, int nSize) {
  for (int i = 0; i < nSize; i++)
    arr[i] = 0;
}

inline void linear_interpolation_coef(amrex::RealVect& dx,
                                      amrex::Real (&coef)[2][2][2]) {

  amrex::Real xi[2];
  amrex::Real eta[2];
  amrex::Real zeta[2];
  xi[0] = dx[0];
  eta[0] = dx[1];
  zeta[0] = dx[2];
  xi[1] = 1 - xi[0];
  eta[1] = 1 - eta[0];
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

inline amrex::Real get_value_at_loc(const amrex::MultiFab& mf,
                                    const amrex::MFIter& mfi,
                                    const amrex::Geometry& gm,
                                    const amrex::RealVect xyz, const int iVar) {
  const auto plo = gm.ProbLo();

  const auto invDx = gm.InvCellSize();

  int loIdx[nDim];
  amrex::Real dx[nDim];
  for (int i = 0; i < nDim; i++) {
    dx[i] = (xyz[i] - plo[i]) * invDx[i];
    loIdx[i] = fastfloor(dx[i]);
    dx[i] = dx[i] - loIdx[i];
  }

  amrex::Real coef[2][2][2];
  {
    amrex::Real xi[2];
    amrex::Real eta[2];
    amrex::Real zeta[2];
    xi[0] = dx[0];
    eta[0] = dx[1];
    zeta[0] = dx[2];
    xi[1] = 1 - xi[0];
    eta[1] = 1 - eta[0];
    zeta[1] = 1 - zeta[0];

    amrex::Real multi[2][2];
    multi[0][0] = xi[0] * eta[0];
    multi[0][1] = xi[0] * eta[1];
    multi[1][0] = xi[1] * eta[0];
    multi[1][1] = xi[1] * eta[1];

    // coef[k][j][i]: This is so wired. But is may be faster since it matches
    // the AMREX multifab data ordering.
    if (nDim == 2) {
      coef[0][0][0] = multi[1][1];
      coef[0][0][1] = multi[0][1];
      coef[0][1][0] = multi[1][0];
      coef[0][1][1] = multi[0][0];
    } else {
      coef[0][0][0] = multi[1][1] * zeta[1];
      coef[1][0][0] = multi[1][1] * zeta[0];
      coef[0][1][0] = multi[1][0] * zeta[1];
      coef[1][1][0] = multi[1][0] * zeta[0];
      coef[0][0][1] = multi[0][1] * zeta[1];
      coef[1][0][1] = multi[0][1] * zeta[0];
      coef[0][1][1] = multi[0][0] * zeta[1];
      coef[1][1][1] = multi[0][0] * zeta[0];
    }
  }

  const auto& arr = mf[mfi].array();
  amrex::Real val = 0;

  amrex::Box box = amrex::Box(amrex::IntVect(0), amrex::IntVect(1));

  amrex::ParallelFor(box,
                     [&] AMREX_GPU_DEVICE(int ii, int jj, int kk) noexcept {
                       int iIdx = loIdx[ix_] + ii;
                       int jIdx = loIdx[iy_] + jj;
                       int kIdx = nDim > 2 ? loIdx[iz_] + kk : 0;
                       val += arr(iIdx, jIdx, kIdx, iVar) * coef[kk][jj][ii];
                     });

  return val;
}

inline amrex::Real get_value_at_loc(const amrex::MultiFab& mf,
                                    const amrex::Geometry& gm,
                                    const amrex::RealVect xyz, const int iVar) {
  auto idx = gm.CellIndex(xyz.begin());

  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    // Cell box
    const amrex::Box& bx =
        amrex::convert(mfi.validbox(), { AMREX_D_DECL(0, 0, 0) });
    if (bx.contains(idx))
      return get_value_at_loc(mf, mfi, gm, xyz, iVar);
  }

  amrex::AllPrint() << "xyz = " << xyz << std::endl;

  amrex::Abort("Error: can not find this point!");
  return -1; // To suppress compiler warnings.
}

inline void add_to_mf(const amrex::Real& val, amrex::MultiFab& mf,
                      const amrex::MFIter& mfi, const amrex::Geometry& gm,
                      const amrex::Real x, const amrex::Real y,
                      const amrex::Real z, const int iVar) {
  const auto plo = gm.ProbLo();
  const amrex::RealVect loc = { AMREX_D_DECL(x, y, z) };

  const auto invDx = gm.InvCellSize();

  int loIdx[3];
  amrex::Real dx[3];
  for (int i = 0; i < 3; i++) {
    dx[i] = (loc[i] - plo[i]) * invDx[i];
    loIdx[i] = fastfloor(dx[i]);
    dx[i] = dx[i] - loIdx[i];
  }

  amrex::Real coef[2][2][2];
  {
    amrex::Real xi[2];
    amrex::Real eta[2];
    amrex::Real zeta[2];
    xi[0] = dx[0];
    eta[0] = dx[1];
    zeta[0] = dx[2];
    xi[1] = 1 - xi[0];
    eta[1] = 1 - eta[0];
    zeta[1] = 1 - zeta[0];

    amrex::Real multi[2][2];
    multi[0][0] = xi[0] * eta[0];
    multi[0][1] = xi[0] * eta[1];
    multi[1][0] = xi[1] * eta[0];
    multi[1][1] = xi[1] * eta[1];

    // coef[k][j][i]
    coef[0][0][0] = multi[1][1] * zeta[1];
    coef[1][0][0] = multi[1][1] * zeta[0];
    coef[0][1][0] = multi[1][0] * zeta[1];
    coef[1][1][0] = multi[1][0] * zeta[0];
    coef[0][0][1] = multi[0][1] * zeta[1];
    coef[1][0][1] = multi[0][1] * zeta[0];
    coef[0][1][1] = multi[0][0] * zeta[1];
    coef[1][1][1] = multi[0][0] * zeta[0];
  }

  const amrex::Array4<amrex::Real>& arr = mf[mfi].array();
  for (int kk = 0; kk < 2; kk++) {
    const int kIdx = loIdx[iz_] + kk;
    for (int jj = 0; jj < 2; jj++) {
      const int jIdx = loIdx[iy_] + jj;
      for (int ii = 0; ii < 2; ii++) {
        arr(loIdx[ix_] + ii, jIdx, kIdx, iVar) += val * coef[kk][jj][ii];
      }
    }
  }

  return;
}

template <typename T>
inline T bound(const T& val, const T& xmin, const T& xmax) {
  if (val < xmin)
    return xmin;
  if (val > xmax)
    return xmax;
  return val;
}

double read_mem_usage();

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

#endif
