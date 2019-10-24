#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <AMReX_REAL.H>

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