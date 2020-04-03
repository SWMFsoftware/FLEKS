#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_iMultiFab.H>

#include "Constants.h"
#include "Timer.h"

inline int myfloor(amrex::Real x) { return (int)(x+8)-8; }

void curl_center_to_node(const amrex::MultiFab& centerMF,
                         amrex::MultiFab& nodeMF, const amrex::Real* invDx);

void curl_node_to_center(const amrex::MultiFab& nodeMF,
                         amrex::MultiFab& centerMF, const amrex::Real* invDx);

void lap_node_to_node(const amrex::MultiFab& srcMF, amrex::MultiFab& dstMF,
                      const amrex::DistributionMapping dm,
                      const amrex::Geometry& geom,
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
                    amrex::Geometry& geom, int nshift = 0);

inline int get_local_node_or_cell_number(const amrex::MultiFab& MF) {

  int nTotal = 0;
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

inline void linear_interpolation_coef(amrex::Real (&dx)[3],
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
  multi[0][0] = xi[0]*eta[0];
  multi[0][1] = xi[0]*eta[1];
  multi[1][0] = xi[1]*eta[0];
  multi[1][1] = xi[1]*eta[1];

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
                                    const amrex::Geometry& geom,
                                    const amrex::Real x, const amrex::Real y,
                                    const amrex::Real z, const int iVar) {
  const auto plo = geom.ProbLo();
  const amrex::Real loc[nDimMax] = { x, y, z };

  const auto invDx = geom.InvCellSize();

  int loIdx[3];
  amrex::Real dShift[3];
  for (int i = 0; i < 3; i++) {
    dShift[i] = (loc[i] - plo[i]) * invDx[i];
    loIdx[i] = myfloor(dShift[i]);
    dShift[i] = dShift[i] - loIdx[i];
  }

  amrex::Real coef[2][2][2];
  // Not a good name.
  linear_interpolation_coef(dShift, coef);

  amrex::Real val = 0;
  for (int kk = 0; kk < 2; kk++)
    for (int jj = 0; jj < 2; jj++)
      for (int ii = 0; ii < 2; ii++) {
        val += get_value_at_node(mf, mfi, loIdx[ix_] + ii, loIdx[iy_] + jj,
                                 loIdx[iz_] + kk, iVar) *
               coef[ii][jj][kk];
      }
  return val;
}

inline amrex::Real get_value_at_loc(const amrex::MultiFab& mf,
                                    const amrex::Geometry& geom,
                                    const amrex::Real x, const amrex::Real y,
                                    const amrex::Real z, const int iVar) {
  amrex::Real loc[3] = { x, y, z };
  auto idx = geom.CellIndex(loc);

  for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.validbox();
    if (bx.contains(idx))
      return get_value_at_loc(mf, mfi, geom, x, y, z, iVar);
  }

  amrex::AllPrint() << "loc = " << loc[ix_] << " " << loc[iy_] << " "
                    << loc[iz_] << " idx = " << idx << std::endl;

  amrex::Abort("Error: can not find this point!");
  return -1; // To suppress compiler warnings.
}

template <typename T>
inline T bound(const T& val, const T& xmin, const T& xmax) {
  if (val < xmin)
    return xmin;
  if (val > xmax)
    return xmax;
  return val;
}

#endif