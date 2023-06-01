#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_iMultiFab.H>

#include "Constants.h"
#include "Timer.h"

// Only works for x>-8;
inline int fastfloor(amrex::Real x) { return (int)(x + 8) - 8; }

inline bool test_bit(int i, int pos) { return i & (1 << pos); }

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
                    int nshift = 0);

void print_MultiFab(const amrex::MultiFab& data, std::string tag,
                    amrex::Geometry& gm, int nshift = 0);

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

inline amrex::Real get_value_at_node(const amrex::MultiFab& mf,
                                     const amrex::MFIter& mfi, const int i,
                                     const int j, const int k, const int iVar) {
  const auto& arr = mf[mfi].array();
  return arr(i, j, k, iVar);
}

inline amrex::Real get_value_at_loc(const amrex::MultiFab& mf,
                                    const amrex::MFIter& mfi,
                                    const amrex::Geometry& gm,
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

  const auto& arr = mf[mfi].array();
  amrex::Real val = 0;
  for (int kk = 0; kk < 2; kk++) {
    const int kIdx = loIdx[iz_] + kk;
    for (int jj = 0; jj < 2; jj++) {
      const int jIdx = loIdx[iy_] + jj;
      for (int ii = 0; ii < 2; ii++) {
        val += arr(loIdx[ix_] + ii, jIdx, kIdx, iVar) * coef[kk][jj][ii];
      }
    }
  }

  return val;
}

inline amrex::Real get_value_at_loc(const amrex::MultiFab& mf,
                                    const amrex::Geometry& gm,
                                    const amrex::Real x, const amrex::Real y,
                                    const amrex::Real z, const int iVar) {
  amrex::Real loc[3] = { x, y, z };
  auto idx = gm.CellIndex(loc);

  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    // Cell box
    const amrex::Box& bx =
        amrex::convert(mfi.validbox(), { AMREX_D_DECL(0, 0, 0) });
    if (bx.contains(idx))
      return get_value_at_loc(mf, mfi, gm, x, y, z, iVar);
  }

  amrex::AllPrint() << "loc = " << loc[ix_] << " " << loc[iy_] << " "
                    << loc[iz_] << " idx = " << idx << std::endl;

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
#endif
