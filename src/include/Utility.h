#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include "Constants.h"

void convert_1d_to_3d(const double* const p, amrex::MultiFab& MF,
                      amrex::Geometry& geom);

void convert_3d_to_1d(const amrex::MultiFab& MF, double* const p,
                      amrex::Geometry& geom);

void curl_center_to_node(const amrex::MultiFab& centerMF,
                         amrex::MultiFab& nodeMF, const amrex::Real* invDx);

void curl_node_to_center(const amrex::MultiFab& nodeMF,
                         amrex::MultiFab& centerMF, const amrex::Real* invDx);

void lap_node_to_node(const amrex::MultiFab& srcMF, amrex::MultiFab& dstMF,
                      const amrex::DistributionMapping dm,
                      const amrex::Geometry& geom);

void grad_node_to_center(const amrex::MultiFab& nodeMF,
                         amrex::MultiFab& centerMF, const amrex::Real* invDx);

void grad_center_to_node(const amrex::MultiFab& centerMF,
                         amrex::MultiFab& nodeMF, const amrex::Real* invDx);

void div_center_to_node(const amrex::MultiFab& centerMF,
                        amrex::MultiFab& nodeMF, const amrex::Real* invDx);

void div_node_to_center(const amrex::MultiFab& nodeMF,
                        amrex::MultiFab& centerMF, const amrex::Real* invDx);

void div_center_to_center(const amrex::MultiFab& srcMF,
                        amrex::MultiFab& dstMF, const amrex::Real* invDx);

void average_center_to_node(const amrex::MultiFab& centerMF,
                            amrex::MultiFab& nodeMF);

void print_MultiFab(amrex::MultiFab& data, std::string tag);

inline int get_fab_grid_points_number(const amrex::MultiFab& MF) {
  amrex::MFIter mfi(MF);
  const amrex::Box& box = mfi.validbox();
  const auto lo = lbound(box);
  const auto hi = ubound(box);

  return (hi.x - lo.x + 1) * (hi.y - lo.y + 1) * (hi.z - lo.z + 1);
}

template <class T> inline void zero_array(T* arr, int nSize) {
  for (int i = 0; i < nSize; i++)
    arr[i] = 0;
}

inline void part_grid_interpolation_coef(amrex::Real (&dx)[3],
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

  coef[0][0][0] = xi[1] * eta[1] * zeta[1];
  coef[0][0][1] = xi[1] * eta[1] * zeta[0];
  coef[0][1][0] = xi[1] * eta[0] * zeta[1];
  coef[0][1][1] = xi[1] * eta[0] * zeta[0];
  coef[1][0][0] = xi[0] * eta[1] * zeta[1];
  coef[1][0][1] = xi[0] * eta[1] * zeta[0];
  coef[1][1][0] = xi[0] * eta[0] * zeta[1];
  coef[1][1][1] = xi[0] * eta[0] * zeta[0];
}

#endif