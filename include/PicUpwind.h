#ifndef _PIC_UPWIND_H_
#define _PIC_UPWIND_H_

#include <cmath>

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

#include "Constants.h"

// MUSCL slope limiter for the Pe advection and the #HYBRIDUPWIND / #UPWINDB B
// reconstruction.
// limType: 0 = upwind1, 1 = minmod, 2 = van leer, 3 = monotonized central (mc).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real pe_slope(
    const int limType, const amrex::Real dL, const amrex::Real dR) noexcept {
  if (limType == 0 || dL * dR <= 0.0)
    return 0.0;
  if (limType == 1) { // minmod
    return (dL > 0.0) ? amrex::min(dL, dR) : amrex::max(dL, dR);
  } else if (limType == 2) { // van leer
    return 2.0 * dL * dR / (dL + dR);
  } else if (limType == 3) { // mc
    const amrex::Real c = 0.5 * (dL + dR);
    const amrex::Real s = (dL > 0.0) ? 1.0 : -1.0;
    return s * amrex::min(2.0 * std::abs(dL),
                          amrex::min(2.0 * std::abs(dR), std::abs(c)));
  }
  return 0.0;
}

// Limited jump phi * (B^R - B^L) of component iVar of the MUSCL-reconstructed
// cell-centered B across the low face of cell (i,j,k) along DIR, i.e. between
// the cells at offsets -1 and 0, from the cached cell slopes (below). The
// upwind term at a node averages this over the 2 (2D) or 4 (3D) faces around
// the node in the transverse directions, the same cells curl_center_to_node
// averages over; computing each face once and averaging avoids reconstructing
// it again for every node that shares it.
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real upwind_face_jump(
    const amrex::Array4<amrex::Real const>& arrB,
    const amrex::Array4<amrex::Real const>& slope, const int i, const int j,
    const int k, const int iVar, const amrex::Real sensorFloor) noexcept {
  constexpr int di = (DIR == ix_) ? 1 : 0;
  constexpr int dj = (DIR == iy_) ? 1 : 0;
  constexpr int dk = (DIR == iz_) ? 1 : 0;
  const amrex::Real bM1 = arrB(i - di, j - dj, k - dk, iVar);
  const amrex::Real b0 = arrB(i, j, k, iVar);
  const amrex::Real bL = bM1 + 0.5 * slope(i - di, j - dj, k - dk, iVar);
  const amrex::Real bR = b0 - 0.5 * slope(i, j, k, iVar);
  amrex::Real phi = 1.0;
  if (sensorFloor < 1.0) {
    // Grid-scale sensor s = |third difference| / sum|first differences|:
    // ~ (k dx)^2 for smooth modes, O(1) for 2-3 cell noise and jumps.
    const amrex::Real bM2 = arrB(i - 2 * di, j - 2 * dj, k - 2 * dk, iVar);
    const amrex::Real bP1 = arrB(i + di, j + dj, k + dk, iVar);
    const amrex::Real dM = bM1 - bM2, d0 = b0 - bM1, dP = bP1 - b0;
    const amrex::Real den = std::abs(dM) + std::abs(d0) + std::abs(dP);
    const amrex::Real s =
        (den > 0.0) ? std::abs(dP - 2.0 * d0 + dM) / den : 0.0;
    phi = amrex::min(1.0, amrex::max(sensorFloor, s));
  }
  return phi * (bR - bL);
}

// Limited MUSCL slope of component iVar of the cell-centered B in cell (i,j,k)
// along DIR. Each cell's slope serves the faces on both of its sides, so it is
// computed once per direction instead of once per face.
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real upwind_cell_slope(
    const amrex::Array4<amrex::Real const>& arrB, const int i, const int j,
    const int k, const int iVar, const int limType) noexcept {
  constexpr int di = (DIR == ix_) ? 1 : 0;
  constexpr int dj = (DIR == iy_) ? 1 : 0;
  constexpr int dk = (DIR == iz_) ? 1 : 0;
  const amrex::Real bM1 = arrB(i - di, j - dj, k - dk, iVar);
  const amrex::Real b0 = arrB(i, j, k, iVar);
  const amrex::Real bP1 = arrB(i + di, j + dj, k + dk, iVar);
  return pe_slope(limType, b0 - bM1, bP1 - b0);
}

// Average of the face jumps around node (i,j,k) along DIR (see above).
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void upwind_node_jump(
    const amrex::Array4<amrex::Real const>& faceJump, const int i, const int j,
    const int k, amrex::Real jump[nDim3]) noexcept {
  constexpr int t1 = (DIR == ix_) ? iy_ : ix_;
  constexpr int t2 = (DIR == iz_) ? iy_ : iz_;
  const int o2Min = (nDim > 2 || t2 != iz_) ? -1 : 0;

  for (int iVar = 0; iVar < nDim3; ++iVar)
    jump[iVar] = 0.0;

  int nCell = 0;
  for (int o1 = -1; o1 <= 0; ++o1) {
    for (int o2 = o2Min; o2 <= 0; ++o2) {
      int c[nDim3] = { i, j, k };
      c[t1] += o1;
      c[t2] += o2;
      for (int iVar = 0; iVar < nDim3; ++iVar)
        jump[iVar] += faceJump(c[ix_], c[iy_], c[iz_], iVar);
      ++nCell;
    }
  }

  const amrex::Real invN = 1.0 / nCell;
  for (int iVar = 0; iVar < nDim3; ++iVar)
    jump[iVar] *= invN;
}

#endif // _PIC_UPWIND_H_
