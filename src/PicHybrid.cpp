#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_PhysBCFunct.H>

#include "GridUtility.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

namespace {

//==========================================================
// Helpers for the evolved electron pressure.

// Cell (i,j,k) shifted by `off` cells along direction DIR.
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
pe_at(const Array4<Real const>& arrPe, const int i, const int j, const int k,
      const int off) noexcept {
  if constexpr (DIR == ix_) {
    return arrPe(i + off, j, k);
  } else if constexpr (DIR == iy_) {
    return arrPe(i, j + off, k);
  } else {
    return arrPe(i, j, k + off);
  }
}

// Component DIR of the nodal velocity, shifted by `off` cells along DIR.
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
ue_at(const Array4<Real const>& arrUe, const int i, const int j, const int k,
      const int off) noexcept {
  if constexpr (DIR == ix_) {
    return arrUe(i + off, j, k, DIR);
  } else if constexpr (DIR == iy_) {
    return arrUe(i, j + off, k, DIR);
  } else {
    return arrUe(i, j, k + off, DIR);
  }
}

// MUSCL slope limiter for the Pe advection and the #HYBRIDUPWIND B
// reconstruction.
// limType: 0 = upwind1, 1 = minmod, 2 = van leer, 3 = monotonized central (mc).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real pe_slope(const int limType,
                                                       const Real dL,
                                                       const Real dR) noexcept {
  if (limType == 0 || dL * dR <= 0.0)
    return 0.0;
  if (limType == 1) { // minmod
    return (dL > 0.0) ? amrex::min(dL, dR) : amrex::max(dL, dR);
  } else if (limType == 2) { // van leer
    return 2.0 * dL * dR / (dL + dR);
  } else if (limType == 3) { // mc
    const Real c = 0.5 * (dL + dR);
    const Real s = (dL > 0.0) ? 1.0 : -1.0;
    return s * amrex::min(2.0 * std::abs(dL),
                          amrex::min(2.0 * std::abs(dR), std::abs(c)));
  }
  return 0.0;
}

// Add one direction's contribution to div(u_e Pe) and div(u_e). The caller
// skips a direction that carries no flux by testing invDx > 0 first.
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void pe_dir_flux(
    const Array4<Real const>& arrPe, const Array4<Real const>& arrUe,
    const int i, const int j, const int k, const Real invDx, const int limType,
    const Real peMin, Real& divFlux, Real& divUe) noexcept {
  const Real uLo = ue_at<DIR>(arrUe, i, j, k, 0);
  const Real uHi = ue_at<DIR>(arrUe, i, j, k, 1);

  Real peLo;
  if (limType == 0) {
    peLo = (uLo > 0.0) ? pe_at<DIR>(arrPe, i, j, k, -1)
                       : pe_at<DIR>(arrPe, i, j, k, 0);
  } else {
    const Real sLoL = pe_slope(
        limType,
        pe_at<DIR>(arrPe, i, j, k, -1) - pe_at<DIR>(arrPe, i, j, k, -2),
        pe_at<DIR>(arrPe, i, j, k, 0) - pe_at<DIR>(arrPe, i, j, k, -1));
    const Real sLoR = pe_slope(
        limType, pe_at<DIR>(arrPe, i, j, k, 0) - pe_at<DIR>(arrPe, i, j, k, -1),
        pe_at<DIR>(arrPe, i, j, k, 1) - pe_at<DIR>(arrPe, i, j, k, 0));
    const Real pL =
        amrex::max(pe_at<DIR>(arrPe, i, j, k, -1) + 0.5 * sLoL, peMin);
    const Real pR =
        amrex::max(pe_at<DIR>(arrPe, i, j, k, 0) - 0.5 * sLoR, peMin);
    peLo = (uLo > 0.0) ? pL : ((uLo < 0.0) ? pR : 0.5 * (pL + pR));
  }

  Real peHi;
  if (limType == 0) {
    peHi = (uHi > 0.0) ? pe_at<DIR>(arrPe, i, j, k, 0)
                       : pe_at<DIR>(arrPe, i, j, k, 1);
  } else {
    const Real sHiL = pe_slope(
        limType, pe_at<DIR>(arrPe, i, j, k, 0) - pe_at<DIR>(arrPe, i, j, k, -1),
        pe_at<DIR>(arrPe, i, j, k, 1) - pe_at<DIR>(arrPe, i, j, k, 0));
    const Real sHiR = pe_slope(
        limType, pe_at<DIR>(arrPe, i, j, k, 1) - pe_at<DIR>(arrPe, i, j, k, 0),
        pe_at<DIR>(arrPe, i, j, k, 2) - pe_at<DIR>(arrPe, i, j, k, 1));
    const Real pL =
        amrex::max(pe_at<DIR>(arrPe, i, j, k, 0) + 0.5 * sHiL, peMin);
    const Real pR =
        amrex::max(pe_at<DIR>(arrPe, i, j, k, 1) - 0.5 * sHiR, peMin);
    peHi = (uHi > 0.0) ? pL : ((uHi < 0.0) ? pR : 0.5 * (pL + pR));
  }

  divFlux += (uHi * peHi - uLo * peLo) * invDx;
  divUe += (uHi - uLo) * invDx;
}

// Limited jump phi * (B^R - B^L) of component iVar of the MUSCL-reconstructed
// cell-centered B across the low face of cell (i,j,k) along DIR, i.e. between
// the cells at offsets -1 and 0, from the cached cell slopes (below). The
// upwind term at a node averages this over the 2 (2D) or 4 (3D) faces around
// the node in the transverse directions, the same cells curl_center_to_node
// averages over; computing each face once and averaging avoids reconstructing
// it again for every node that shares it.
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
upwind_face_jump(const Array4<Real const>& arrB,
                 const Array4<Real const>& slope, const int i, const int j,
                 const int k, const int iVar, const Real sensorFloor) noexcept {
  constexpr int di = (DIR == ix_) ? 1 : 0;
  constexpr int dj = (DIR == iy_) ? 1 : 0;
  constexpr int dk = (DIR == iz_) ? 1 : 0;
  const Real bM1 = arrB(i - di, j - dj, k - dk, iVar);
  const Real b0 = arrB(i, j, k, iVar);
  const Real bL = bM1 + 0.5 * slope(i - di, j - dj, k - dk, iVar);
  const Real bR = b0 - 0.5 * slope(i, j, k, iVar);
  Real phi = 1.0;
  if (sensorFloor < 1.0) {
    // Grid-scale sensor s = |third difference| / sum|first differences|:
    // ~ (k dx)^2 for smooth modes, O(1) for 2-3 cell noise and jumps.
    const Real bM2 = arrB(i - 2 * di, j - 2 * dj, k - 2 * dk, iVar);
    const Real bP1 = arrB(i + di, j + dj, k + dk, iVar);
    const Real dM = bM1 - bM2, d0 = b0 - bM1, dP = bP1 - b0;
    const Real den = std::abs(dM) + std::abs(d0) + std::abs(dP);
    const Real s = (den > 0.0) ? std::abs(dP - 2.0 * d0 + dM) / den : 0.0;
    phi = amrex::min(1.0, amrex::max(sensorFloor, s));
  }
  return phi * (bR - bL);
}

// Limited MUSCL slope of component iVar of the cell-centered B in cell (i,j,k)
// along DIR. Each cell's slope serves the faces on both of its sides, so it is
// computed once per direction instead of once per face.
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
upwind_cell_slope(const Array4<Real const>& arrB, const int i, const int j,
                  const int k, const int iVar, const int limType) noexcept {
  constexpr int di = (DIR == ix_) ? 1 : 0;
  constexpr int dj = (DIR == iy_) ? 1 : 0;
  constexpr int dk = (DIR == iz_) ? 1 : 0;
  const Real bM1 = arrB(i - di, j - dj, k - dk, iVar);
  const Real b0 = arrB(i, j, k, iVar);
  const Real bP1 = arrB(i + di, j + dj, k + dk, iVar);
  return pe_slope(limType, b0 - bM1, bP1 - b0);
}

// Average of the face jumps around node (i,j,k) along DIR (see above).
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void upwind_node_jump(
    const Array4<Real const>& faceJump, const int i, const int j, const int k,
    Real jump[nDim3]) noexcept {
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

  const Real invN = 1.0 / nCell;
  for (int iVar = 0; iVar < nDim3; ++iVar)
    jump[iVar] *= invN;
}

// Algebraic polytropic electron closure the evolved equation starts from:
//   Pe = P0 * (rho/rho0)^gamma,  or Pe = Te*rho when gamma == 1.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
polytropic_closure(const Real rho, const Real p0, const Real invRho0,
                   const Real gamma, const Real Te) noexcept {
  if (gamma == 1.0) {
    return Te * rho;
  }
  return (rho > 0) ? p0 * std::pow(rho * invRho0, gamma) : 0.0;
}

} // namespace

//==========================================================
void Pic::assemble_ohm_E(const MultiFab& centerBin,
                         const MultiFab& centerBtimeAvg, MultiFab& Eout,
                         int iLev, Real hstep, bool includeAmbi) {
  BL_PROFILE("Pic::assemble_ohm_E");

  // Nodal total current J = curl(B)/(4*pi) from trial B (compact 1*dx stencil
  // from cell centers to nodes). Only needed for physical resistivity and Hall.
  const bool needJ =
      (etaResistivity > 0 || useHallTerm || hasRegionalResistivity_);
  if (needJ) {
    curl_center_to_node(centerBin, nodeJ[iLev], Geom(iLev).InvCellSize());
    nodeJ[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  // Magnetic field interpolated from cell centers to nodes for vector cross
  // products.
  average_center_to_node(centerBtimeAvg, nodeBstage[iLev]);
  nodeBstage[iLev].FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(node_status(iLev), nodeBstage[iLev], 0, 3, &Pic::get_node_B,
                   iLev, true);
  }
  if (is_body_conducting()) {
    project_body_B(nodeBstage[iLev], iLev);
    nodeBstage[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  // Moment time-interpolation weights:
  // For standard leapfrog: X(hstep) = (0.5-hstep)*X^{n-1/2} +
  // (0.5+hstep)*X^{n+1/2}. For predictor-corrector: X(hstep) = (1.0-hstep)*X^n
  // + hstep*X^{n+1}.
  const Real wPrev = usePredictorCorrector ? (1.0 - hstep) : (0.5 - hstep);
  const Real wCur = usePredictorCorrector ? hstep : (0.5 + hstep);
  const Real invFourPI = 1.0 / fourPI;

  // #HYBRIDUPWIND: directions with a single cell carry no jump.
  const bool doUpwind = useHybridUpwind;
  GpuArray<int, nDim3> upwindDir = { 0, 0, 0 };
  GpuArray<Real, nDim3> invDxUp = { 0.0, 0.0, 0.0 };
  for (int iDim = 0; iDim < nDim; ++iDim) {
    upwindDir[iDim] = (Geom(iLev).Domain().length(iDim) > 1) ? 1 : 0;
    invDxUp[iDim] = Geom(iLev).InvCellSize(iDim);
  }
  const int limTypeUp = upwindLimiterType;
  const Real sensorFloorUp = upwindSensorFloor;
  const Real coefWhistlerUp = useHallTerm ? upwindCoefWhistler : 0.0;
  const Real fixedSpeedUp = upwindFixedSpeed;
  const Real cs2Up =
      (electronTemperature > 0)
          ? ((electronGamma > 0) ? electronGamma : 1.0) * electronTemperature
          : 0.0;
  const bool useHall = useHallTerm;

  for (MFIter mfi(Eout); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrE = Eout[mfi].array();
    const Array4<Real const>& arrB = nodeBstage[iLev][mfi].array();
    const Array4<Real const>& moments = nodePlasma[nSpecies][iLev][mfi].array();
    const Array4<Real const>& momentsPrev =
        nodePlasmaPrev[nSpecies][iLev][mfi].array();
    const Array4<Real const> arrJ =
        needJ ? nodeJ[iLev][mfi].array() : Array4<Real const>();
    const Array4<Real const> arrEambi = (electronTemperature > 0)
                                            ? nodeEambi[iLev][mfi].array()
                                            : Array4<Real const>();
    const Array4<Real const> arrEtaReg =
        hasRegionalResistivity_ ? nodeEtaRegional[iLev][mfi].array()
                                : Array4<Real const>();
    // Upwind face jumps, one array per direction over the faces the nodes of
    // this box use: the node range along the direction, and one more cell on
    // the low side of each transverse direction (the cells c-1 and c around a
    // node). They read the cells i-2..i+1 along the direction, within the two
    // ghost layers of centerBin.
    FArrayBox faceJump[nDim3];
    Array4<Real const> arrFace[nDim3];
    if (doUpwind) {
      const Array4<Real const>& arrBin = centerBin[mfi].const_array();
      FArrayBox slopeFab;
      for (int iDir = 0; iDir < nDim; ++iDir) {
        if (!upwindDir[iDir])
          continue;
        Box faceBox(box.smallEnd(), box.bigEnd());
        for (int t = 0; t < nDim; ++t) {
          if (t != iDir)
            faceBox.growLo(t, 1);
        }
        // A face uses the slopes of the cells on both of its sides.
        Box slopeBox = faceBox;
        slopeBox.growLo(iDir, 1);
        slopeFab.resize(slopeBox, nDim3, The_Async_Arena());
        faceJump[iDir].resize(faceBox, nDim3, The_Async_Arena());
        const Array4<Real>& arrS = slopeFab.array();
        const Array4<Real const>& arrSc = slopeFab.const_array();
        const Array4<Real>& arrF = faceJump[iDir].array();
        if (iDir == ix_) {
          ParallelFor(slopeBox, nDim3,
                      [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                        arrS(i, j, k, iVar) = upwind_cell_slope<ix_>(
                            arrBin, i, j, k, iVar, limTypeUp);
                      });
          ParallelFor(faceBox, nDim3,
                      [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                        arrF(i, j, k, iVar) = upwind_face_jump<ix_>(
                            arrBin, arrSc, i, j, k, iVar, sensorFloorUp);
                      });
        } else if (iDir == iy_) {
          ParallelFor(slopeBox, nDim3,
                      [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                        arrS(i, j, k, iVar) = upwind_cell_slope<iy_>(
                            arrBin, i, j, k, iVar, limTypeUp);
                      });
          ParallelFor(faceBox, nDim3,
                      [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                        arrF(i, j, k, iVar) = upwind_face_jump<iy_>(
                            arrBin, arrSc, i, j, k, iVar, sensorFloorUp);
                      });
        } else {
          ParallelFor(slopeBox, nDim3,
                      [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                        arrS(i, j, k, iVar) = upwind_cell_slope<iz_>(
                            arrBin, i, j, k, iVar, limTypeUp);
                      });
          ParallelFor(faceBox, nDim3,
                      [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                        arrF(i, j, k, iVar) = upwind_face_jump<iz_>(
                            arrBin, arrSc, i, j, k, iVar, sensorFloorUp);
                      });
        }
        arrFace[iDir] = faceJump[iDir].const_array();
      }
    }
    const Array4<Real const> arrFaceX = arrFace[ix_];
    const Array4<Real const> arrFaceY = arrFace[iy_];
    const Array4<Real const> arrFaceZ = arrFace[iz_];

    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real rhoPrev = momentsPrev(i, j, k, iRho_);
      const Real rhoCur = moments(i, j, k, iRho_);
      const Real rho = wPrev * rhoPrev + wCur * rhoCur;
      const Real mx =
          wPrev * momentsPrev(i, j, k, iMx_) + wCur * moments(i, j, k, iMx_);
      const Real my =
          wPrev * momentsPrev(i, j, k, iMy_) + wCur * moments(i, j, k, iMy_);
      const Real mz =
          wPrev * momentsPrev(i, j, k, iMz_) + wCur * moments(i, j, k, iMz_);
      Real ui = 0, vi = 0, wi = 0;
      const Real invRhoEff =
          (rho > 0) ? (1.0 / amrex::max(rho, rhoMinOhm)) : 0.0;

      if (rho > 0) {
        ui = mx * invRhoEff;
        vi = my * invRhoEff;
        wi = mz * invRhoEff;
      }

      Real bx = arrB(i, j, k, ix_);
      Real by = arrB(i, j, k, iy_);
      Real bz = arrB(i, j, k, iz_);

      // Convection term: E = -U_i x B
      Real ex = -(vi * bz - wi * by);
      Real ey = -(wi * bx - ui * bz);
      Real ez = -(ui * by - vi * bx);

      // J = curl(B)/(4*pi)
      Real jx = 0.0, jy = 0.0, jz = 0.0;
      if (needJ) {
        jx = arrJ(i, j, k, ix_) * invFourPI;
        jy = arrJ(i, j, k, iy_) * invFourPI;
        jz = arrJ(i, j, k, iz_) * invFourPI;
      }

      // eta * J (global + regional resistivity)
      Real etaTotal = etaResistivity;
      if (hasRegionalResistivity_) {
        etaTotal += arrEtaReg(i, j, k);
      }
      if (etaTotal > 0.0) {
        ex += etaTotal * jx;
        ey += etaTotal * jy;
        ez += etaTotal * jz;
      }

      // Ambipolar electric field: E_ambi = -grad(p_e)/(e*n_e)
      // Analytically curl(E_ambi) == 0 for isothermal/polytropic electrons.
      if (includeAmbi && electronTemperature > 0) {
        ex += arrEambi(i, j, k, ix_);
        ey += arrEambi(i, j, k, iy_);
        ez += arrEambi(i, j, k, iz_);
      }

      // Hall term: (J x B) / rho_q
      if (rho > 0 && useHallTerm) {
        Real hall_x = (jy * bz - jz * by) * invRhoEff;
        Real hall_y = (jz * bx - jx * bz) * invRhoEff;
        Real hall_z = (jx * by - jy * bx) * invRhoEff;

        ex += hall_x;
        ey += hall_y;
        ez += hall_z;
      }

      // Upwind (Rusanov) dissipation: E += sum_d (a_d/2) e_d x (B^R - B^L).
      // The jumps use the evolved B only (no intrinsic B0), the speeds the
      // total B and the electron velocity u_e = u_i - J/rho.
      if (doUpwind && (rho > 0 || fixedSpeedUp > 0)) {
        Real cFast = 0.0, bAbs = 0.0;
        Real ue[nDim3] = { ui, vi, wi };
        if (fixedSpeedUp <= 0) {
          const Real b2 = bx * bx + by * by + bz * bz;
          bAbs = std::sqrt(b2);
          cFast = std::sqrt(b2 * invRhoEff * invFourPI + cs2Up);
          if (useHall) {
            ue[ix_] -= jx * invRhoEff;
            ue[iy_] -= jy * invRhoEff;
            ue[iz_] -= jz * invRhoEff;
          }
        }

        Real jump[nDim3];
        for (int iDir = 0; iDir < nDim; ++iDir) {
          if (!upwindDir[iDir])
            continue;
          if (iDir == ix_) {
            upwind_node_jump<ix_>(arrFaceX, i, j, k, jump);
          } else if (iDir == iy_) {
            upwind_node_jump<iy_>(arrFaceY, i, j, k, jump);
          } else {
            upwind_node_jump<iz_>(arrFaceZ, i, j, k, jump);
          }

          // Whistler phase speed at k = pi/dx: pi |B| / (4 pi rho dx).
          const Real a = (fixedSpeedUp > 0) ? fixedSpeedUp
                                            : std::abs(ue[iDir]) + cFast +
                                                  coefWhistlerUp * 0.25 * bAbs *
                                                      invRhoEff * invDxUp[iDir];
          const Real halfA = 0.5 * a;

          // e_d x jump
          if (iDir == ix_) {
            ey -= halfA * jump[iz_];
            ez += halfA * jump[iy_];
          } else if (iDir == iy_) {
            ex += halfA * jump[iz_];
            ez -= halfA * jump[ix_];
          } else {
            ex -= halfA * jump[iy_];
            ey += halfA * jump[ix_];
          }
        }
      }

      arrE(i, j, k, ix_) = ex;
      arrE(i, j, k, iy_) = ey;
      arrE(i, j, k, iz_) = ez;
    });
  }

  // Hyper-resistivity: E -= (eta_h / 4*pi) * curl(nabla^2 B).
  // centerLapB = Laplacian(centerBin);
  // nodeHyperE = curl_center_to_node(centerLapB).
  const bool doHyper = (etaHyperLev[iLev] > 0 || hasRegionalHyper_);
  if (doHyper) {
    lap_center_to_center(centerBin, centerLapB[iLev], Geom(iLev).InvCellSize());
    centerLapB[iLev].FillBoundary(Geom(iLev).periodicity());
    // The first ghost layer of B is Faraday-advanced, so the Laplacian
    // there (needed by the curl at the interface nodes) is computed
    // directly from it and must not be replaced by coarse interpolation.
    if (iLev > 0) {
      fill_fine_lev_outer_bny_from_coarse(
          centerLapB[iLev - 1], centerLapB[iLev], 0, centerLapB[iLev].nComp(),
          refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          *get_cell_interp(), faraday_ngrow(iLev));
    }
    // Dirichlet faces (fixed / inflow / coupled) must get the Laplacian of
    // the prescribed B, not B itself: copying B here leaves a jump of size B0
    // at the boundary cell, and its curl drives a spurious hyper-resistive E
    // at the nodes next to the face. The prescribed state is uniform, so zero.
    apply_field_bc(cell_status(iLev), centerLapB[iLev], 0,
                   centerLapB[iLev].nComp(), &Pic::get_zero, iLev, true);

    curl_center_to_node(centerLapB[iLev], nodeHyperE[iLev],
                        Geom(iLev).InvCellSize());
    nodeHyperE[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev > 0) {
      fill_fine_lev_bny_from_coarse(
          nodeHyperE[iLev - 1], nodeHyperE[iLev], 0, nodeHyperE[iLev].nComp(),
          refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), node_status(iLev),
          node_bilinear_interp);
    }
    // Likewise curl(nabla^2 B) of the prescribed state is zero, not the
    // prescribed (motional) E.
    apply_field_bc(node_status(iLev), nodeHyperE[iLev], 0,
                   nodeHyperE[iLev].nComp(), &Pic::get_zero, iLev, false);

    const Real fGlobal = etaHyperLev[iLev] / fourPI;

    if (!hasRegionalHyper_) {
      MultiFab::Saxpy(Eout, -fGlobal, nodeHyperE[iLev], 0, 0, nDim3, 0);
    } else {
      for (MFIter mfi(Eout); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real>& arrE = Eout[mfi].array();
        const Array4<Real const>& arrHyp = nodeHyperE[iLev][mfi].array();
        const Array4<Real const> arrHyperReg =
            nodeEtaHyperRegional[iLev][mfi].array();

        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          Real f = fGlobal + arrHyperReg(i, j, k) * invFourPI;
          if (f > 0.0) {
            arrE(i, j, k, ix_) -= f * arrHyp(i, j, k, ix_);
            arrE(i, j, k, iy_) -= f * arrHyp(i, j, k, iy_);
            arrE(i, j, k, iz_) -= f * arrHyp(i, j, k, iz_);
          }
        });
      }
    }
  }

  Eout.FillBoundary(Geom(iLev).periodicity());
  apply_field_bc(node_status(iLev), Eout, 0, nDim3, &Pic::get_node_E, iLev,
                 false);
  if (use_body()) {
    apply_body_E_bc(Eout, iLev);
  }
}

//==========================================================
// Ambipolar electric field: E_ambi = -grad(p_e) / (e * n_e).
// Precomputed once per PIC timestep outside the magnetic subcycling
// steps, avoiding redundant gradient and EOS operations during the
// high-frequency whistler integration.
void Pic::compute_ambipolar_E() {
  std::string nameFunc = "Pic::compute_ambipolar_E";
  timing_func(nameFunc);

  if (electronTemperature <= 0) {
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      nodeEambi[iLev].setVal(0.0);
    }
    return;
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    compute_ambipolar_E(iLev);
  }

  if (finestLevel() > 0) {
    for (int iLev = 1; iLev < n_lev(); ++iLev) {
      fill_fine_lev_bny_from_coarse(
          nodeEambi[iLev - 1], nodeEambi[iLev], 0, nDim3, refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    }
  }
}

//==========================================================
void Pic::compute_ambipolar_E(int iLev) {
  if (electronTemperature <= 0) {
    nodeEambi[iLev].setVal(0.0);
    return;
  }

  if (useElectronPressureEq) {
    // Pe is an evolved field: take it as it is instead of rebuilding it from
    // the density with the polytropic closure.
    MultiFab::Copy(centerPe[iLev], centerPeState[iLev], 0, 0, 1,
                   centerPe[iLev].nGrow());
    centerPe[iLev].FillBoundary(Geom(iLev).periodicity());
  } else {
    // Copy nodal ion density to nodeRhoTemp and average it to cell centers.
    compute_electron_density(iLev, nodeRhoTemp[iLev], centerPe[iLev], false);
    centerPe[iLev].FillBoundary(Geom(iLev).periodicity());

    // Evaluate electron pressure Pe at cell centers via EOS
    const Real p0 = electronDensity0 * electronTemperature;
    const Real invRho0 =
        (electronDensity0 > 0) ? (1.0 / electronDensity0) : 0.0;
    const Real gamma = electronGamma;
    const Real Te = electronTemperature;

    for (MFIter mfi(centerPe[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrPe = centerPe[iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        arrPe(i, j, k) =
            polytropic_closure(arrPe(i, j, k), p0, invRho0, gamma, Te);
      });
    }
    centerPe[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  if (iLev > 0) {
    fill_fine_lev_bny_from_coarse(
        centerPe[iLev - 1], centerPe[iLev], 0, 1, refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), cell_status(iLev), *get_cell_interp());
  }

  // Zero-gradient (Neumann / foextrap) BC across non-periodic domain boundaries
  if (!Geom(iLev).isAllPeriodic() && centerPe[iLev].nGrow() > 0) {
    apply_pe_zero_gradient_bc(iLev, centerPe[iLev]);
  }

  apply_fake2d_k_clamp(centerPe[iLev]);

  // Compute grad_center_to_node(Pe) directly into nodeEambi
  grad_center_to_node(centerPe[iLev], nodeEambi[iLev],
                      Geom(iLev).InvCellSize());

  // Scale in-place: E_ambi = -grad(Pe) / max(rho, rhoMinOhm)
  for (MFIter mfi(nodeEambi[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrEambi = nodeEambi[iLev][mfi].array();
    const Array4<Real const>& moments = nodePlasma[nSpecies][iLev][mfi].array();

    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real rho = moments(i, j, k, iRho_);
      if (rho > 0) {
        const Real invRhoEff = -1.0 / amrex::max(rho, rhoMinOhm);
        arrEambi(i, j, k, ix_) *= invRhoEff;
        arrEambi(i, j, k, iy_) *= invRhoEff;
        arrEambi(i, j, k, iz_) *= invRhoEff;
      } else {
        arrEambi(i, j, k, ix_) = 0.0;
        arrEambi(i, j, k, iy_) = 0.0;
        arrEambi(i, j, k, iz_) = 0.0;
      }
    });
  }

  // At physical conducting/reflecting walls, the normal ambipolar field
  // vanishes analytically (dPe/dn = 0)
  if (!Geom(iLev).isAllPeriodic()) {
    const BoundaryBounds bnd(Geom(iLev), nodeEambi[iLev].boxArray().ixType(),
                             &bcField);
    for (MFIter mfi(nodeEambi[iLev]); mfi.isValid(); ++mfi) {
      const Box& bxValid = mfi.validbox();
      Array4<Real> const& arr = nodeEambi[iLev][mfi].array();
      const Dim3 vLo = bxValid.smallEnd().dim3();
      const Dim3 vHi = bxValid.bigEnd().dim3();
      ParallelFor(bxValid, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const int ijk[3] = { i, j, k };
        const int vLoArr[3] = { vLo.x, vLo.y, vLo.z };
        const int vHiArr[3] = { vHi.x, vHi.y, vHi.z };
        for (int d = 0; d < nDim; ++d) {
          if (!bnd.isNode[d])
            continue;
          const bool onLoWall =
              (bnd.bcLo[d] == FieldBC::conducting) && (ijk[d] == bnd.loBnd[d]);
          const bool onHiWall =
              (bnd.bcHi[d] == FieldBC::conducting) && (ijk[d] == bnd.hiBnd[d]);
          if (onLoWall || onHiWall) {
            bool inValid = true;
            for (int od = 0; od < nDim; ++od) {
              if (od != d && (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
                inValid = false;
                break;
              }
            }
            if (inValid) {
              arr(i, j, k, d) = 0.0;
            }
          }
        }
      });
    }
  }

  if (is_fake_2d()) {
    for (amrex::MFIter mfi(nodeEambi[iLev]); mfi.isValid(); ++mfi) {
      const auto& vbox = mfi.validbox();
      const auto& fbox = mfi.fabbox();
      auto arr = nodeEambi[iLev][mfi].array();
      const int klo = vbox.smallEnd(2);
      const int khi = vbox.bigEnd(2);
      amrex::ParallelFor(fbox,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                           const int k_src = std::clamp(k, klo, khi);
                           if (k != k_src) {
                             for (int n = 0; n < nDim3; ++n) {
                               arr(i, j, k, n) = arr(i, j, k_src, n);
                             }
                           }
                         });
    }
  }

  nodeEambi[iLev].FillBoundary(Geom(iLev).periodicity());
  apply_field_bc(node_status(iLev), nodeEambi[iLev], 0, nDim3, &Pic::get_node_E,
                 iLev, false);
  if (use_body()) {
    apply_body_E_bc(nodeEambi[iLev], iLev);
  }
}

//==========================================================
void Pic::smooth_moments() {
  std::string nameFunc = "Pic::smooth_moments";
  timing_func(nameFunc);

  const bool shouldSmooth = useExplicitPIC
                                ? ((doSmoothJ && nSmoothJ > 0) ||
                                   (doSmoothMoments && nSmoothMoments > 0))
                                : (doSmoothMoments && nSmoothMoments > 0);
  if (!shouldSmooth)
    return;

  const int nSmooth = (useExplicitPIC && doSmoothJ) ? nSmoothJ : nSmoothMoments;
  const Real coefSmooth =
      (useExplicitPIC && doSmoothJ) ? coefSmoothJ : coefSmoothMoments;

  // Smooth the total plasma moments on the node grid.
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab& moments = nodePlasma[nSpecies][iLev];
    moments.FillBoundary(Geom(iLev).periodicity());
    if (iLev > 0) {
      fill_fine_lev_bny_from_coarse(
          nodePlasma[nSpecies][iLev - 1], nodePlasma[nSpecies][iLev], 0,
          nodePlasma[nSpecies][iLev].nComp(), refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    }
    for (int icount = 0; icount < nSmooth; ++icount) {
      smooth_multifab(moments, iLev, 1, coefSmooth);
      if (iLev > 0) {
        fill_fine_lev_bny_from_coarse(
            nodePlasma[nSpecies][iLev - 1], nodePlasma[nSpecies][iLev], 0,
            nodePlasma[nSpecies][iLev].nComp(), refRatio(iLev - 1),
            Geom(iLev - 1), Geom(iLev), node_status(iLev),
            node_bilinear_interp);
      }
    }
    moments.FillBoundary(Geom(iLev).periodicity());

    if (!Geom(iLev).isAllPeriodic()) {
      const BoundaryBounds bnd(Geom(iLev), moments.boxArray().ixType(),
                               &bcField);
      for (MFIter mfi(moments); mfi.isValid(); ++mfi) {
        const Box& bxValid = mfi.validbox();
        Array4<Real> const& arr = moments[mfi].array();
        const Dim3 vLo = bxValid.smallEnd().dim3();
        const Dim3 vHi = bxValid.bigEnd().dim3();
        ParallelFor(bxValid, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const int ijk[3] = { i, j, k };
          const int vLoArr[3] = { vLo.x, vLo.y, vLo.z };
          const int vHiArr[3] = { vHi.x, vHi.y, vHi.z };
          for (int d = 0; d < nDim; ++d) {
            if (!bnd.isNode[d])
              continue;
            const bool onLoWall = (bnd.bcLo[d] == FieldBC::conducting) &&
                                  (ijk[d] == bnd.loBnd[d]);
            const bool onHiWall = (bnd.bcHi[d] == FieldBC::conducting) &&
                                  (ijk[d] == bnd.hiBnd[d]);
            if (onLoWall || onHiWall) {
              bool inValid = true;
              for (int od = 0; od < nDim; ++od) {
                if (od != d && (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
                  inValid = false;
                  break;
                }
              }
              if (inValid) {
                if (d == 0) {
                  arr(i, j, k, iMx_) = 0.0;
                  arr(i, j, k, iPxy_) = 0.0;
                  arr(i, j, k, iPxz_) = 0.0;
                } else if (d == 1) {
                  arr(i, j, k, iMy_) = 0.0;
                  arr(i, j, k, iPxy_) = 0.0;
                  arr(i, j, k, iPyz_) = 0.0;
                } else if (d == 2) {
                  arr(i, j, k, iMz_) = 0.0;
                  arr(i, j, k, iPxz_) = 0.0;
                  arr(i, j, k, iPyz_) = 0.0;
                }
              }
            }
          }
        });
      }
    }

    if (use_body()) {
      mask_body(moments, node_status(iLev));
    }
  }

  if (finestLevel() > 0) {
    for (int iLev = 1; iLev < n_lev(); ++iLev) {
      fill_fine_lev_bny_from_coarse(
          nodePlasma[nSpecies][iLev - 1], nodePlasma[nSpecies][iLev], 0,
          nodePlasma[nSpecies][iLev].nComp(), refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    }
  }
}

//==========================================================
// Copy the current summed moment deposit into nodePlasmaPrev (J^{n-1/2})
// before a fresh deposit (J^{n+1/2}), so assemble_ohm_E can time-interpolate
// the two at the magnetic sub-step fraction hstep.
void Pic::save_current_moments_to_prev() {
  std::string nameFunc = "Pic::save_current_moments_to_prev";
  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::Copy(nodePlasmaPrev[nSpecies][iLev], nodePlasma[nSpecies][iLev],
                   0, 0, nHybridMomentsComps,
                   nodePlasma[nSpecies][iLev].nGrow());
    nodePlasmaPrev[nSpecies][iLev].FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
// Seed nodePlasmaPrev on the first hybrid step, where there is no previous
// deposit: initialise it from the current deposit so the time interpolation
// degrades to a plain average for that single step.
void Pic::seed_first_hybrid_step() {
  std::string nameFunc = "Pic::seed_first_hybrid_step";
  timing_func(nameFunc);

  save_current_moments_to_prev();
}

//==========================================================
// Zero-gradient (foextrap) ghost cells for the electron-pressure field and its
// scratch arrays. Same recipe as the one compute_ambipolar_E() uses for
// centerPe: Neumann on every non-periodic face, then FillBoundary.
void Pic::apply_pe_zero_gradient_bc(int iLev, MultiFab& mf) {
  if (Geom(iLev).isAllPeriodic() || mf.nGrow() == 0) {
    mf.FillBoundary(Geom(iLev).periodicity());
    return;
  }

  Vector<BCRec> bcr(mf.nComp());
  for (int c = 0; c < mf.nComp(); ++c) {
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      if (Geom(iLev).isPeriodic(d)) {
        bcr[c].setLo(d, BCType::int_dir);
        bcr[c].setHi(d, BCType::int_dir);
      } else {
        bcr[c].setLo(d, BCType::foextrap);
        bcr[c].setHi(d, BCType::foextrap);
      }
    }
  }
  GpuBndryFuncFab<FabFillNoOp> bfunc(FabFillNoOp{});
  PhysBCFunct<GpuBndryFuncFab<FabFillNoOp> > physbcf(Geom(iLev), bcr, bfunc);
  physbcf(mf, 0, mf.nComp(), mf.nGrowVect(), 0.0, 0);

  if (fi && fi->get_iPe() >= 0 && bcField.has(FieldBC::coupled)) {
    fill_ext_dir(cell_status(iLev), mf, 0, mf.nComp(), &Pic::get_center_Pe,
                 iLev, &bcField);
  }

  mf.FillBoundary(Geom(iLev).periodicity());
}

//==========================================================
// Fake-2D (one cell in z) ghost clamp: every z layer outside the valid range
// takes the value of the nearest valid layer. Single-component fields only; the
// multi-component variant lives inline in apply_centerB_BC().
void Pic::apply_fake2d_k_clamp(MultiFab& mf) {
  if (!is_fake_2d())
    return;

  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const auto& vbox = mfi.validbox();
    const auto& fbox = mfi.fabbox();
    auto arr = mf[mfi].array();
    const int klo = vbox.smallEnd(2);
    const int khi = vbox.bigEnd(2);
    ParallelFor(fbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int k_src = std::clamp(k, klo, khi);
      if (k != k_src)
        arr(i, j, k) = arr(i, j, k_src);
    });
  }
}

//==========================================================
void Pic::apply_centerPe_BC(int iLev) {
  apply_pe_zero_gradient_bc(iLev, centerPeState[iLev]);
  apply_fake2d_k_clamp(centerPeState[iLev]);
}

//==========================================================
void Pic::init_electron_pressure() {
  for (int iLev = 0; iLev < n_lev(); ++iLev)
    init_electron_pressure(iLev);
  peStateInitialized_ = true;
}

//==========================================================
void Pic::init_electron_pressure(int iLev) {
  // Start from the algebraic polytropic closure evaluated on the current ion
  // density, i.e. exactly the state the evolved equation replaces. Electrons
  // then depart from it through advection, compression and conduction.
  compute_electron_density(iLev, nodePeRho[iLev], centerPeRho[iLev], false);

  const Real p0 = electronDensity0 * electronTemperature;
  const Real invRho0 = (electronDensity0 > 0.0) ? 1.0 / electronDensity0 : 0.0;
  const Real gamma = electronGamma;
  const Real Te = electronTemperature;

  for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrPe = centerPeState[iLev][mfi].array();
    const Array4<Real const>& arrRho = centerPeRho[iLev][mfi].array();
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      arrPe(i, j, k) =
          polytropic_closure(arrRho(i, j, k), p0, invRho0, gamma, Te);
    });
  }
  apply_centerPe_BC(iLev);
}

//==========================================================
void Pic::fill_new_electron_pressure() {
  if (!useElectronPressureEq || finestLevel() <= 0)
    return;

  // Boxes created by a regrid have no history; interpolate them from the next
  // coarser level (boxes that already existed keep their copied values).
  auto& cellInterp = *get_cell_interp();
  for (int iLev = 1; iLev < n_lev(); ++iLev) {
    fill_fine_lev_new_from_coarse(centerPeState[iLev - 1], centerPeState[iLev],
                                  0, 1, refRatio(iLev - 1), Geom(iLev - 1),
                                  Geom(iLev), cell_status(iLev), cellInterp);
  }
}

//==========================================================
// Nodal ion density -> cell-centered n_e, into (nodalRho, cellRho).
// useGrownTile also covers the coarse-fine interface nodes, which
// compute_ambipolar_E() needs; the evolved-Pe path only reads valid nodes back.
void Pic::compute_electron_density(int iLev, MultiFab& nodalRho,
                                   MultiFab& cellRho, const bool useGrownTile) {
  nodalRho.setVal(0.0);
  for (int s = 0; s < nSpecies; ++s) {
    const Real qomFactor = (fi && s < fi->get_nS() && fi->get_MoMi(s) > 0.0)
                               ? (fi->get_QoQi(s) / fi->get_MoMi(s))
                               : 1.0;
    for (MFIter mfi(nodalRho); mfi.isValid(); ++mfi) {
      const Box& box = useGrownTile ? mfi.growntilebox() : mfi.validbox();
      const Array4<Real>& arrRho = nodalRho[mfi].array();
      const Array4<Real const>& moments = nodePlasma[s][iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        arrRho(i, j, k) += qomFactor * moments(i, j, k, iRho_);
      });
    }
  }
  nodalRho.FillBoundary(Geom(iLev).periodicity());
  average_node_to_center(nodalRho, cellRho);
}

//==========================================================
// Te = Pe/n_e at the cell centers, with fresh ghosts. Solver scratch for the
// conductivity, not a diagnostic source; see the Te block in
// write_amrex_field().
void Pic::compute_electron_temperature(int iLev) {
  const Real peMinLocal = peMin;
  const Real rhoFloor = rhoMinOhm;

  for (MFIter mfi(centerPeTe[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrTe = centerPeTe[iLev][mfi].array();
    const Array4<Real const>& arrPe = centerPeState[iLev][mfi].array();
    const Array4<Real const>& arrRho = centerPeRho[iLev][mfi].array();
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real ne = amrex::max(arrRho(i, j, k), rhoFloor);
      arrTe(i, j, k) = amrex::max(arrPe(i, j, k), peMinLocal) / ne;
    });
  }
  apply_pe_zero_gradient_bc(iLev, centerPeTe[iLev]);
}

//==========================================================
// Electron-ion collisional thermal equilibration (heat exchange) hook:
// dPe/dt = (Pi - Pe) / tau_eq following the BATSRUS point-implicit formulation.
void Pic::add_electron_ion_heating(int iLev, Real dt) {
  if (!useHeatExchange || collisionCoefEi <= 0.0)
    return;

  // 1) Compute scalar ion pressure Pi = Tr(P_i)/3 at nodes in nodePeRho
  for (MFIter mfi(nodePeRho[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrPiN = nodePeRho[iLev][mfi].array();
    const Array4<Real const>& moments = nodePlasma[nSpecies][iLev][mfi].array();
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      arrPiN(i, j, k) = (moments(i, j, k, iPxx_) + moments(i, j, k, iPyy_) +
                         moments(i, j, k, iPzz_)) /
                        3.0;
    });
  }
  nodePeRho[iLev].FillBoundary(Geom(iLev).periodicity());

  // 2) Average scalar ion pressure to cell centers in centerPe (scratch here)
  average_node_to_center(nodePeRho[iLev], centerPe[iLev]);
  apply_pe_zero_gradient_bc(iLev, centerPe[iLev]);

  // 3) Point-implicit energy exchange following BATSRUS:
  //    H = collisionCoefEi * ne / Te^1.5
  //    PePImpl = H / (1 + 2 * dt * H)
  //    DeltaPe = dt * PePImpl * (Pi - Pe)
  const Real cEi = collisionCoefEi;
  const Real peMinLocal = peMin;
  const Real rhoFloor = rhoMinOhm;

  for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrPe = centerPeState[iLev][mfi].array();
    const Array4<Real const>& arrPi = centerPe[iLev][mfi].array();
    const Array4<Real const>& arrRho = centerPeRho[iLev][mfi].array();

    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real ne = amrex::max(arrRho(i, j, k), rhoFloor);
      const Real pe = amrex::max(arrPe(i, j, k), peMinLocal);
      const Real pi = amrex::max(arrPi(i, j, k), 0.0);
      const Real te = pe / ne;

      if (te > 0.0) {
        const Real H = cEi * ne / std::pow(te, 1.5);
        const Real pePImpl = H / (1.0 + 2.0 * dt * H);
        const Real deltaPe = dt * pePImpl * (pi - pe);
        arrPe(i, j, k) = amrex::max(pe + deltaPe, peMinLocal);
      }
    });
  }
}

//==========================================================
void Pic::update_Pe_hybrid() {
  if (!useElectronPressureEq)
    return;

  std::string nameFunc = "Pic::update_Pe_hybrid";
  timing_func(nameFunc);

  const Real dt = tc->get_dt();

  if (!peStateInitialized_)
    init_electron_pressure();

  for (int iLev = 0; iLev < n_lev(); ++iLev)
    update_Pe_hybrid(iLev, dt);

  if (projectDownEmFields && finestLevel() > 0) {
    for (int iLev = finestLevel(); iLev > 0; iLev--)
      average_down(centerPeState[iLev], centerPeState[iLev - 1], 0, 1,
                   refRatio(iLev - 1));
  }
}

//==========================================================
void Pic::update_Pe_hybrid(int iLev, Real dt) {
  BL_PROFILE("Pic::update_Pe_hybrid");

  //--------------------------------------------------------------------
  // Operator-split advance of the scalar electron pressure on one level:
  //   dPe/dt + div(u_e Pe) + (gamma_e-1) Pe div(u_e)
  //       = (gamma_e-1) [ div(kappa_hat . grad(Te)) + H_ei ]
  // The order matters: each term acts on the density, velocity and Te that the
  // calls before it established.
  //--------------------------------------------------------------------

  electron_velocity_at_nodes(iLev);
  advect_electron_pressure(iLev, dt);

  compute_electron_density(iLev, nodePeRho[iLev], centerPeRho[iLev], false);
  compute_electron_temperature(iLev);

  // A no-op each when heatCondKappa0 is 0 / #ELECTRONCOLLISION is absent.
  apply_electron_heat_conduction(iLev, dt);
  add_electron_ion_heating(iLev, dt);

  // The three calls above write the valid box only, so the state ghost cells
  // are stale until now. Only the next step's advection reads them, which makes
  // this the one place that has to restore them.
  apply_centerPe_BC(iLev);
}

//==========================================================
// Electron velocity u_e = U_i - J/(e*n_e) on the NODES.
//
// A node at (i+1,j,k) is the x-face between cells i and i+1, so the nodal
// values are exactly the face-normal velocities the upwind fluxes of stage 2
// need, and their difference across a cell is div(u_e).
void Pic::electron_velocity_at_nodes(int iLev) {
  BL_PROFILE("Pic::electron_velocity_at_nodes");

  curl_center_to_node(centerB[iLev], nodeJ[iLev], Geom(iLev).InvCellSize());
  nodeJ[iLev].FillBoundary(Geom(iLev).periodicity());

  const Real invFourPILocal = 1.0 / fourPI;
  const Real rhoFloor = rhoMinOhm;

  for (MFIter mfi(nodePeVec[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrUe = nodePeVec[iLev][mfi].array();
    const Array4<Real const>& arrMom = nodePlasma[nSpecies][iLev][mfi].array();
    const Array4<Real const>& arrJ = nodeJ[iLev][mfi].array();
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real rho = arrMom(i, j, k, iRho_);
      if (rho > 0) {
        const Real invRho = 1.0 / amrex::max(rho, rhoFloor);
        for (int d = 0; d < nDim; ++d)
          arrUe(i, j, k, d) =
              (arrMom(i, j, k, iMx_ + d) - arrJ(i, j, k, d) * invFourPILocal) *
              invRho;
      } else {
        for (int d = 0; d < nDim; ++d)
          arrUe(i, j, k, d) = 0.0;
      }
    });
  }
  apply_pe_zero_gradient_bc(iLev, nodePeVec[iLev]);
}

//==========================================================
// Advection of Pe: dPe/dt = -div(u_e Pe) - (gamma_e-1) Pe div(u_e).
//
// TVD/MUSCL reconstructed face states (upwind1, minmod, vanleer, mc) with
// exponential or explicit compression. The result goes to the centerPe scratch
// and is copied back, because an in-place update would let a thread read a
// neighbour it has already overwritten.
void Pic::advect_electron_pressure(int iLev, Real dt) {
  BL_PROFILE("Pic::advect_electron_pressure");

  const Real* invDxGeom = Geom(iLev).InvCellSize();
  const IntVect domLen = Geom(iLev).Domain().length();
  // A dimension with a single cell carries no flux: zero its inverse spacing.
  const Real invDxX = (domLen[ix_] > 1) ? invDxGeom[ix_] : 0.0;
  const Real invDxY = (domLen[iy_] > 1) ? invDxGeom[iy_] : 0.0;
  const Real invDxZ = (nDim > 2 && domLen[iz_] > 1) ? invDxGeom[iz_] : 0.0;

  const Real gammaM1 = electronGamma - 1.0;
  const int limType = peLimiterType;
  const bool compExp = peCompressionExp;
  const Real peMinLocal = peMin;

  // The only stencil read of the state: the MUSCL reconstruction reaches two
  // cells away, so centerPeState needs both ghost layers.
  apply_centerPe_BC(iLev);

  for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real const>& arrPe = centerPeState[iLev][mfi].array();
    const Array4<Real>& arrPeNew = centerPe[iLev][mfi].array();
    const Array4<Real const>& arrUe = nodePeVec[iLev][mfi].array();

    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      Real divFlux = 0.0;
      Real divUe = 0.0;

      // Keep the x -> y -> z order; the summation order is part of the result.
      if (invDxX > 0.0) {
        pe_dir_flux<ix_>(arrPe, arrUe, i, j, k, invDxX, limType, peMinLocal,
                         divFlux, divUe);
      }
      if (invDxY > 0.0) {
        pe_dir_flux<iy_>(arrPe, arrUe, i, j, k, invDxY, limType, peMinLocal,
                         divFlux, divUe);
      }
      if constexpr (nDim > 2) {
        if (invDxZ > 0.0) {
          pe_dir_flux<iz_>(arrPe, arrUe, i, j, k, invDxZ, limType, peMinLocal,
                           divFlux, divUe);
        }
      }

      const Real pe = arrPe(i, j, k);
      if (compExp) {
        const Real pAdv = amrex::max(pe - dt * divFlux, peMinLocal);
        arrPeNew(i, j, k) =
            amrex::max(pAdv * std::exp(-gammaM1 * divUe * dt), peMinLocal);
      } else {
        arrPeNew(i, j, k) =
            amrex::max(pe - dt * (divFlux + gammaM1 * pe * divUe), peMinLocal);
      }
    });
  }
  MultiFab::Copy(centerPeState[iLev], centerPe[iLev], 0, 0, 1, 0);
}

//==========================================================
// Electron heat conduction:
//   dPe/dt = (gamma_e-1) div(kappa_hat . grad(Te)) = -(gamma_e-1) div(q),
// with q = -kappa_hat . grad(Te) and the Spitzer kappa = kappa0 Te^2.5.
//
// "point-implicit": runs nCondIter Jacobi sweeps of the point-implicit update
// "subcycle": takes explicit sub-steps bounded by a diffusion CFL estimate
void Pic::apply_electron_heat_conduction(int iLev, Real dt) {
  BL_PROFILE("Pic::apply_electron_heat_conduction");

  const Real* invDxGeom = Geom(iLev).InvCellSize();
  const IntVect domLen = Geom(iLev).Domain().length();
  const Real invDxX = (domLen[ix_] > 1) ? invDxGeom[ix_] : 0.0;
  const Real invDxY = (domLen[iy_] > 1) ? invDxGeom[iy_] : 0.0;
  const Real invDxZ = (nDim > 2 && domLen[iz_] > 1) ? invDxGeom[iz_] : 0.0;

  const Real gammaM1 = electronGamma - 1.0;
  // Fraction of the conductivity dyad that is field aligned: 1 = kappa*b*b,
  // 0 = isotropic.
  const Real fAlign = fieldAlignedConduction ? fieldAlignedFraction : 0.0;
  // m_p/m_e, needed by the free-streaming heat-flux limiter. |qomEl| is
  // |q/m| in units of e/m_p, so for the electron it is exactly m_p/m_e.
  const Real massRatioPe = std::abs(qomEl);

  if (heatCondKappa0 > 0.0) {
    const bool isExplicit = (heatCondMethod == "subcycle");
    const int nPasses = isExplicit ? 1 : nCondIter;
    Real dtRem = dt;
    int pass = 0;

    while (dtRem > 0.0 && pass < (isExplicit ? nCondSubcycleMax : nPasses)) {
      pass++;

      if (pass > 1) {
        for (MFIter mfi(centerPeTe[iLev]); mfi.isValid(); ++mfi) {
          const Box& box = mfi.validbox();
          const Array4<Real>& arrTe = centerPeTe[iLev][mfi].array();
          const Array4<Real const>& arrPe = centerPeState[iLev][mfi].array();
          const Array4<Real const>& arrRho = centerPeRho[iLev][mfi].array();
          ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const Real ne = amrex::max(arrRho(i, j, k), rhoMinOhm);
            arrTe(i, j, k) = amrex::max(arrPe(i, j, k), peMin) / ne;
          });
        }
        apply_pe_zero_gradient_bc(iLev, centerPeTe[iLev]);
      }

      average_center_to_node(centerPeTe[iLev], nodePeAux[iLev]);
      grad_center_to_node(centerPeTe[iLev], nodePeVec[iLev],
                          Geom(iLev).InvCellSize());
      average_center_to_node(centerB[iLev], nodeBstage[iLev]);
      nodeBstage[iLev].FillBoundary(Geom(iLev).periodicity());

      for (MFIter mfi(nodePeVec[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real>& arrQ = nodePeVec[iLev][mfi].array();
        const Array4<Real>& arrAux = nodePeAux[iLev][mfi].array();
        const Array4<Real const>& arrB = nodeBstage[iLev][mfi].array();
        const Array4<Real const>& arrRhoN = nodePeRho[iLev][mfi].array();

        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const Real teN = amrex::max(arrAux(i, j, k, 0), 0.0);
          Real kappa = (teN > 0.0) ? heatCondKappa0 * std::pow(teN, 2.5) : 0.0;

          // Unit vector along B, and a smooth weak-field weight: a field at
          // round-off level still has bmag > 0, but its direction is arbitrary,
          // so the conduction must not be locked to it. The weight blends the
          // dyad towards isotropic as |B| drops below fieldAlignedBMin.
          Real bhat[3] = { 0.0, 0.0, 0.0 };
          Real fEff = 0.0;
          if (fAlign > 0.0) {
            const Real bx = arrB(i, j, k, ix_);
            const Real by = arrB(i, j, k, iy_);
            const Real bz = arrB(i, j, k, iz_);
            const Real b2 = bx * bx + by * by + bz * bz;
            if (b2 > 0.0) {
              const Real bmag = std::sqrt(b2);
              bhat[ix_] = bx / bmag;
              bhat[iy_] = by / bmag;
              bhat[iz_] = bz / bmag;
              fEff = fAlign * b2 / (b2 + fieldAlignedBMin * fieldAlignedBMin);
            }
          }

          Real gT[3] = { 0.0, 0.0, 0.0 };
          for (int d = 0; d < nDim; ++d)
            gT[d] = arrQ(i, j, k, d);
          const Real bdotg =
              bhat[ix_] * gT[ix_] + bhat[iy_] * gT[iy_] + bhat[iz_] * gT[iz_];

          Real q[3] = { 0.0, 0.0, 0.0 };
          for (int d = 0; d < 3; ++d) {
            if (d >= nDim || domLen[d] <= 1) {
              q[d] = 0.0;
              continue;
            }
            q[d] = -kappa * (fEff * bhat[d] * bdotg + (1.0 - fEff) * gT[d]);
          }

          // Free-streaming saturation: |q| <= f * n_e * Te * v_th,e with
          // v_th,e = sqrt(Te * m_p/m_e) in code units.
          if (heatFluxLimiter > 0.0) {
            const Real qmag =
                std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
            if (qmag > 0.0) {
              const Real ne = amrex::max(arrRhoN(i, j, k), rhoMinOhm);
              const Real qsat =
                  heatFluxLimiter * ne * teN * std::sqrt(teN * massRatioPe);
              if (qmag > qsat) {
                const Real scale = qsat / qmag;
                for (int d = 0; d < 3; ++d)
                  q[d] *= scale;
                kappa *= scale;
              }
            }
          }

          for (int d = 0; d < 3; ++d) {
            arrQ(i, j, k, d) = q[d];
            arrAux(i, j, k, 1 + d) =
                kappa * (fEff * bhat[d] * bhat[d] + (1.0 - fEff));
          }
        });
      }
      apply_pe_zero_gradient_bc(iLev, nodePeVec[iLev]);
      apply_pe_zero_gradient_bc(iLev, nodePeAux[iLev]);

      div_node_to_center(nodePeVec[iLev], centerPe[iLev],
                         Geom(iLev).InvCellSize());

      Real subDt = dt / static_cast<Real>(nPasses);
      bool useExplicitThisStep = isExplicit;

      if (isExplicit) {
        Real maxDiffRate = 0.0;
        for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
          const Box& box = mfi.validbox();
          const Array4<Real const>& arrRho = centerPeRho[iLev][mfi].array();
          const Array4<Real const>& arrK = nodePeAux[iLev][mfi].array();
          Real localMaxRate = 0.0;
          amrex::Loop(box, [&](int i, int j, int k) {
            const Real ne = amrex::max(arrRho(i, j, k), rhoMinOhm);
            Real lam =
                (arrK(i + 1, j, k, 1) + arrK(i, j, k, 1)) * invDxX * invDxX;
            lam += (arrK(i, j + 1, k, 2) + arrK(i, j, k, 2)) * invDxY * invDxY;
            if (nDim > 2)
              lam +=
                  (arrK(i, j, k + 1, 3) + arrK(i, j, k, 3)) * invDxZ * invDxZ;
            localMaxRate = amrex::max(localMaxRate, gammaM1 * lam / ne);
          });
          maxDiffRate = amrex::max(maxDiffRate, localMaxRate);
        }
        amrex::ParallelDescriptor::ReduceRealMax(maxDiffRate);

        if (maxDiffRate > 0.0) {
          const Real dtCFL = 0.45 / maxDiffRate;
          subDt = amrex::min(dtRem, dtCFL);
        } else {
          subDt = dtRem;
        }

        if (pass == nCondSubcycleMax && subDt < dtRem) {
          subDt = dtRem;
          useExplicitThisStep = false;
        }
      }

      for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real>& arrPe = centerPeState[iLev][mfi].array();
        const Array4<Real const>& arrDivQ = centerPe[iLev][mfi].array();
        const Array4<Real const>& arrRho = centerPeRho[iLev][mfi].array();
        const Array4<Real const>& arrK = nodePeAux[iLev][mfi].array();

        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const Real ne = amrex::max(arrRho(i, j, k), rhoMinOhm);

          Real lambda =
              (arrK(i + 1, j, k, 1) + arrK(i, j, k, 1)) * invDxX * invDxX;
          lambda += (arrK(i, j + 1, k, 2) + arrK(i, j, k, 2)) * invDxY * invDxY;
          if (nDim > 2)
            lambda +=
                (arrK(i, j, k + 1, 3) + arrK(i, j, k, 3)) * invDxZ * invDxZ;

          const Real src = -arrDivQ(i, j, k);
          const Real denom =
              useExplicitThisStep ? 1.0 : (1.0 + subDt * gammaM1 * lambda / ne);

          arrPe(i, j, k) =
              amrex::max(arrPe(i, j, k) + subDt * gammaM1 * src / denom, peMin);
        });
      }

      dtRem -= subDt;
      apply_centerPe_BC(iLev);
    }
  }
}

//==========================================================
// Boundary conditions for cell-centered B, applied to the RK trial states and
// to the new state at the end of each B sub-step.
//
// For refined levels (iLev > 0), ghost cells at the coarse-fine boundary are
// interpolated from the coarse level, and physical boundary conditions are
// applied. When faraday_ngrow > 0, the first fine ghost layer is updated by
// the Faraday stage instead of being overwritten by coarse interpolation.
void Pic::apply_centerB_BC(int iLev) { apply_centerB_BC(iLev, centerB[iLev]); }

void Pic::apply_centerB_BC(int iLev, amrex::MultiFab& mfB) {
  mfB.FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(cell_status(iLev), mfB, 0, mfB.nComp(), &Pic::get_center_B,
                   iLev, true);
  } else {
    MultiFab& coarseB = (&mfB == &centerBstage[iLev])  ? centerBstage[iLev - 1]
                        : (&mfB == &centerBstar[iLev]) ? centerBstar[iLev - 1]
                                                       : centerB[iLev - 1];
    if (faraday_ngrow(iLev) > 0) {
      // The first ghost layer is advanced by the Faraday update with the
      // synchronized E; only the outer layers are interpolated.
      fill_fine_lev_outer_bny_from_coarse(
          coarseB, mfB, 0, mfB.nComp(), refRatio(iLev - 1), Geom(iLev - 1),
          Geom(iLev), cell_status(iLev), *get_cell_interp(),
          faraday_ngrow(iLev));
      if (is_bz_div_free() && mfB.nComp() > iz_) {
        // Without a z extent only Bx, By enter div(B). A free-running ghost Bz
        // is not needed and piles up where the outflow crosses the interface,
        // so it is taken from the coarse level in the whole ghost region.
        // (nDim alone is not enough: a fake-2D deck in a 3D build has nDim ==
        // 3 but the same single z cell.)
        fill_fine_lev_bny_from_coarse(coarseB, mfB, iz_, 1, refRatio(iLev - 1),
                                      Geom(iLev - 1), Geom(iLev),
                                      cell_status(iLev), *get_cell_interp());
      }
    } else {
      fill_fine_lev_bny_from_coarse(
          coarseB, mfB, 0, mfB.nComp(), refRatio(iLev - 1), Geom(iLev - 1),
          Geom(iLev), cell_status(iLev), *get_cell_interp());
    }
    apply_field_bc(cell_status(iLev), mfB, 0, mfB.nComp(), &Pic::get_center_B,
                   iLev, true);
  }

  if (is_body_conducting()) {
    project_body_B(mfB, iLev);
  }

  if (is_fake_2d()) {
    const int nComp = mfB.nComp();
    for (amrex::MFIter mfi(mfB); mfi.isValid(); ++mfi) {
      const auto& vbox = mfi.validbox();
      const auto& fbox = mfi.fabbox();
      auto arr = mfB[mfi].array();
      const int klo = vbox.smallEnd(2);
      const int khi = vbox.bigEnd(2);
      amrex::ParallelFor(fbox,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                           const int k_src = std::clamp(k, klo, khi);
                           if (k != k_src) {
                             for (int n = 0; n < nComp; ++n) {
                               arr(i, j, k, n) = arr(i, j, k_src, n);
                             }
                           }
                         });
    }
  }
}

//==========================================================
// Synchronize nodal E at the coarse-fine interfaces, so that both levels
// advance the cells that share an interface face with the same E, which keeps
// div(B) at round-off across the interface:
//   1) the coarse nodes on the interface take the E of the coincident fine
//      nodes (finest level first);
//   2) the fine interface nodes take the coarse E, linearly interpolated along
//      the coarse edges (fill_fine_lev_edge_from_coarse). On the coincident
//      nodes this returns the fine E of step 1; the midpoints get the average
//      of their two neighbors, so the circulation of E along every coarse
//      interface edge is the same on both levels.
// (The ghost nodes of the fine levels are filled afterwards in
// fill_ghost_emf_from_coarse.)
//
// Why the interface E comes from the fine level: with the coarse E there, the
// covered-cell relaxation closes a feedback loop. Each particle step it
// copies the fine B onto the covered coarse cells; the coarse grid propagates
// a whistler that is 3 coarse cells long with the wrong dispersion and, through
// the interface E, drives the fine level with the wrong phase. Parallel
// whistlers crossing an interface normal to B then grow at a rate set by the
// particle step, not the field sub-step (a linear 1D model of the scheme gives
// 38 Omega_ci for n = 0.2 and dx_coarse = 0.1 d_i, wavelength ~6 fine cells,
// as measured). With the fine E on the interface, the fine level only sees the
// coarse level through its ghost cells and the loop is neutral.
//
// The covered coarse nodes inside the interface keep the coarse E. Injecting
// the coincident fine E into them does not make the covered coarse B follow
// the fine one: curl_node_to_center is a 2*dx-wide stencil, so the coarse
// curl of any nodal E differs from the cell average of the fine curl at O(1)
// relative error (0.5-0.8 measured on the AMR tests). The covered coarse B is
// kept close to the fine average by relax_covered_B_to_fine instead. That
// injection was also harmful: it copies the grid-scale part of the fine E onto
// the coarse nodes, where it drives a 2*dx_coarse mode that the coarse level
// barely damped (the grid-mode hyper-resistivity was set by the finest dx). On
// a 2-level GEM deck (base dx = 0.1 d_i) that mode grew along the interface
// from the first cycle and blew the run up within an ion gyroperiod.
void Pic::sync_emf_at_interface(Vector<MultiFab>& nodeEmf) {
  BL_PROFILE("Pic::sync_emf_at_interface");
  for (int iLev = finestLevel(); iLev >= 1; --iLev) {
    // A node shared by two fine boxes can hold two values (each box computes
    // E from its own copy of the evolved ghost cells); take the owner's, so
    // that every coarse box receives the same one.
    nodeEmf[iLev].OverrideSync(Geom(iLev).periodicity());
    MultiFab& crse = nodeEmf[iLev - 1];
    // Restrict onto the coarsened fine layout, then copy with the coarse
    // periodicity: when the fine level touches a periodic face on one side
    // only, the coarse node on the opposite face is the same interface node
    // and must get the fine value too.
    MultiFab crseOfFine(
        amrex::coarsen(nodeEmf[iLev].boxArray(), refRatio(iLev - 1)),
        nodeEmf[iLev].DistributionMap(), nDim3, 0);
    average_down_nodal(nodeEmf[iLev], crseOfFine, refRatio(iLev - 1));
    MultiFab fromFine(crse.boxArray(), crse.DistributionMap(), nDim3, 0);
    fromFine.setVal(0.0);
    fromFine.ParallelCopy(crseOfFine, 0, 0, nDim3, IntVect(0), IntVect(0),
                          Geom(iLev - 1).periodicity());
    for (MFIter mfi(crse); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& e = crse[mfi].array();
      const auto& ef = fromFine[mfi].const_array();
      const auto& cell = cell_status(iLev - 1)[mfi].const_array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        // An interface node touches a covered cell and an uncovered cell
        // inside the domain.
        bool hasCovered = false, hasUncovered = false;
        const int kmin = nDim > 2 ? k - 1 : k;
        for (int kk = kmin; kk <= k; ++kk)
          for (int jj = j - 1; jj <= j; ++jj)
            for (int ii = i - 1; ii <= i; ++ii) {
              const int st = cell(ii, jj, kk);
              if (bit::is_refined(st))
                hasCovered = true;
              else if (!bit::is_domain_boundary(st))
                hasUncovered = true;
            }
        if (hasCovered && hasUncovered)
          for (int iVar = 0; iVar < nDim3; ++iVar)
            e(i, j, k, iVar) = ef(i, j, k, iVar);
      });
    }
    crse.OverrideSync(Geom(iLev - 1).periodicity());
    crse.FillBoundary(Geom(iLev - 1).periodicity());
  }

  for (int iLev = 1; iLev <= finestLevel(); ++iLev) {
    fill_fine_lev_edge_from_coarse(
        nodeEmf[iLev - 1], nodeEmf[iLev], 0, nDim3, refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    nodeEmf[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_field_bc(node_status(iLev), nodeEmf[iLev], 0, nDim3, &Pic::get_node_E,
                   iLev, false);
    if (use_body()) {
      apply_body_E_bc(nodeEmf[iLev], iLev);
    }
  }
}

//==========================================================
// Fill the ghost nodes of the fine-level E from the (synchronized) coarse E.
// Together with the interface nodes (already set from the coarse E), linear
// interpolation along the coarse edges makes the ghost face fluxes sum to the
// coarse face fluxes.
void Pic::fill_ghost_emf_from_coarse(Vector<MultiFab>& nodeEmf) {
  BL_PROFILE("Pic::fill_ghost_emf_from_coarse");
  for (int iLev = 1; iLev <= finestLevel(); ++iLev) {
    fill_fine_lev_bny_from_coarse(
        nodeEmf[iLev - 1], nodeEmf[iLev], 0, nDim3, refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    apply_field_bc(node_status(iLev), nodeEmf[iLev], 0, nDim3, &Pic::get_node_E,
                   iLev, false);
    if (use_body()) {
      apply_body_E_bc(nodeEmf[iLev], iLev);
    }
  }
}

//==========================================================
// kStage[iLev][iK] = curl(E_Ohm(Bin)) on all levels. E is assembled on
// every level first, then synchronized across the coarse-fine interfaces, and
// only then is the curl taken, so that every cell of a level (including the
// covered coarse cells and the first fine ghost layer) is
// advanced by one single nodal E field.
void Pic::faraday_stage_rhs(Vector<MultiFab>& Bin, int iK, Real hstep,
                            bool includeAmbi) {
  BL_PROFILE("Pic::faraday_stage_rhs");
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    assemble_ohm_E(Bin[iLev], total_center_B(Bin[iLev], iLev), nodeEstage[iLev],
                   iLev, hstep, includeAmbi);
  }

  if (finestLevel() > 0) {
    sync_emf_at_interface(nodeEstage);
    fill_ghost_emf_from_coarse(nodeEstage);
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    // curl_node_to_center also fills the first ghost layer, which is used by
    // the stage update when faraday_ngrow(iLev) > 0.
    curl_node_to_center(nodeEstage[iLev], kStage[iLev][iK],
                        Geom(iLev).InvCellSize());
    if (is_body_interior_frozen())
      mask_body_interior(kStage[iLev][iK], cell_status(iLev));
  }
}

//==========================================================
// Compute max |avg(B_fine) - B_coarse| over the covered cells of iLev,
// evaluated before divergence-free relaxation updates them.
void Pic::measure_covered_B_drift(int iLev) {
  BL_PROFILE("Pic::measure_covered_B_drift");
  MultiFab& target = relaxTarget[iLev];
  MultiFab::Copy(target, centerB[iLev], 0, 0, nDim3, 0);
  average_down(centerB[iLev + 1], target, 0, nDim3, refRatio(iLev));
  MultiFab::Subtract(target, centerB[iLev], 0, 0, nDim3, 0);

  Real drift = 0;
  for (MFIter mfi(target); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const auto& d = target[mfi].const_array();
    const auto& st = cell_status(iLev)[mfi].const_array();
    amrex::Loop(box, [&](int i, int j, int k) {
      if (!bit::is_refined(st(i, j, k)))
        return;
      for (int iVar = 0; iVar < nDim3; ++iVar)
        drift = amrex::max(drift, std::abs(d(i, j, k, iVar)));
    });
  }
  ParallelDescriptor::ReduceRealMax(drift);
  coveredBDrift[iLev] = drift;
}

//==========================================================
// Relax covered coarse B toward the cell-average of fine B using
// divergence-free corrections, B += curl(A), where A is a nodal potential
// supported on nodes whose surrounding cells are all covered. Uncovered cells
// and interface fluxes remain untouched, preserving coarse-level div(B) to
// machine precision. A is obtained by solving min |avg(B_fine) - B_coarse -
// curl(A)|^2 with CGLS.
//
// Why this is divergence-free:
//   div_node_to_center o average_center_to_node o curl_node_to_center == 0
// identically in 2D and 3D. Any curl correction strictly preserves div(B)
// regardless of CGLS convergence.
//
// Two variants, selected by is_bz_div_free():
//   * 2D / fake-2D: Bz does not enter div(B), so it is restricted directly from
//     fine cells, while (Bx, By) is relaxed via the scalar potential A_z.
//   * True 3D: all three components are relaxed with the nodal vector potential
//   A.
//     Gauge modes A -> A + grad(phi) lie in the null space of curl, so CGLS
//     does not amplify them and explicit gauge fixing is unnecessary.
void Pic::relax_covered_B_to_fine(int iLev) {
  BL_PROFILE("Pic::relax_covered_B_to_fine");

  MultiFab& B = centerB[iLev];
  MultiFab& target = relaxTarget[iLev];
  MultiFab::Copy(target, B, 0, 0, nDim3, 0);
  average_down(centerB[iLev + 1], target, 0, nDim3, refRatio(iLev));

  const bool bzFree = is_bz_div_free();
  // Number of cell components taking part in the relaxation.
  const int nRelax = bzFree ? 2 : nDim3;
  if (bzFree) {
    MultiFab::Copy(B, target, iz_, iz_, 1, 0);
  }

  MultiFab& r = relaxResidual[iLev];
  MultiFab& q = relaxCurl[iLev];
  MultiFab& s = relaxCurlT[iLev];
  MultiFab& p = relaxDir[iLev];

  // Face-averaged difference operators P_d of curl_node_to_center and their
  // transposes P_d^T, with c = 1/4 in 3D and 1/2 without a z extent. Writing
  // D+ for f(.. +1 ..) - f(..) and S- for f(..) + f(.. -1 ..), each acting on
  // the direction given by its subscript:
  //   P_x = c/dx * D+_x S+_y S+_z,   P_x^T = c/dx * D-_x S-_y S-_z
  // and likewise for y and z, with the sums over the *other* two directions.
  // curl(A) = (P_y A_z - P_z A_y, P_z A_x - P_x A_z, P_x A_y - P_y A_x); its
  // adjoint collects the coefficient of each A component:
  //   s_x = P_z^T r_y - P_y^T r_z, s_y = P_x^T r_z - P_z^T r_x,
  //   s_z = P_y^T r_x - P_x^T r_y.
  const auto invDx = Geom(iLev).InvCellSizeArray();
  const bool useZ = !bzFree;
  const Real cFac = useZ ? 0.25 : 0.5;
  const Real cInvDx = cFac * invDx[ix_];
  const Real cInvDy = cFac * invDx[iy_];
  const Real cInvDz = useZ ? cFac * invDx[iz_] : 0.0;
  // 1 in 3D, 0 without a z extent: the transverse sums then skip that index.
  const int kOff = useZ ? 1 : 0;
  constexpr int maxIter = 30;
  constexpr Real relTol = 1e-8;
  const auto period = Geom(iLev).periodicity();
  const auto owner = s.OwnerMask(period);

  // Precompute inner-node mask once before CGLS: nodes surrounded by covered
  // cells are inner (1), others are boundary/outer (0).
  iMultiFab innerMask(s.boxArray(), s.DistributionMap(), 1, 0);
  for (MFIter mfi(innerMask); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const auto& mask = innerMask[mfi].array();
    const auto& st = cell_status(iLev)[mfi].const_array();
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      bool inner = true;
      for (int kk = 0; kk <= kOff && inner; ++kk)
        for (int jj = 0; jj <= 1 && inner; ++jj)
          for (int ii = 0; ii <= 1 && inner; ++ii)
            if (!bit::is_refined(st(i - ii, j - jj, k - kk)))
              inner = false;
      mask(i, j, k) = inner ? 1 : 0;
    });
  }

  // s = curl^T r on the nodes surrounded by covered cells, zero elsewhere.
  // This is the exact adjoint of apply_curl below.
  auto apply_curlT = [&]() {
    r.FillBoundary(period);
    for (MFIter mfi(s); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& sa = s[mfi].array();
      const auto& d = r[mfi].const_array();
      const auto& mask = innerMask[mfi].const_array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (!mask(i, j, k)) {
          for (int iVar = 0; iVar < nDim3; ++iVar)
            sa(i, j, k, iVar) = 0.0;
          return;
        }

        // (P_z^T f)(i,j,k) = c/dz * sum_{ii,jj} [ f(i-ii,j-jj,k-1) -
        //                                         f(i-ii,j-jj,k)   ]
        Real pzT_y = 0.0, pzT_x = 0.0;
        if (useZ) {
          for (int jj = 0; jj <= 1; ++jj)
            for (int ii = 0; ii <= 1; ++ii) {
              pzT_y +=
                  d(i - ii, j - jj, k - 1, iy_) - d(i - ii, j - jj, k, iy_);
              pzT_x +=
                  d(i - ii, j - jj, k - 1, ix_) - d(i - ii, j - jj, k, ix_);
            }
          pzT_y *= cInvDz;
          pzT_x *= cInvDz;
        }

        // (P_y^T f)(i,j,k) = c/dy * sum_{ii,kk} [ f(i-ii,j-1,k-kk) -
        //                                         f(i-ii,j,  k-kk) ]
        Real pyT_x = 0.0, pyT_z = 0.0;
        for (int kk = 0; kk <= kOff; ++kk)
          for (int ii = 0; ii <= 1; ++ii) {
            pyT_x += d(i - ii, j - 1, k - kk, ix_) - d(i - ii, j, k - kk, ix_);
            pyT_z += d(i - ii, j - 1, k - kk, iz_) - d(i - ii, j, k - kk, iz_);
          }
        pyT_x *= cInvDy;
        pyT_z *= cInvDy;

        // (P_x^T f)(i,j,k) = c/dx * sum_{jj,kk} [ f(i-1,j-jj,k-kk) -
        //                                         f(i,  j-jj,k-kk) ]
        Real pxT_z = 0.0, pxT_y = 0.0;
        for (int kk = 0; kk <= kOff; ++kk)
          for (int jj = 0; jj <= 1; ++jj) {
            pxT_z += d(i - 1, j - jj, k - kk, iz_) - d(i, j - jj, k - kk, iz_);
            pxT_y += d(i - 1, j - jj, k - kk, iy_) - d(i, j - jj, k - kk, iy_);
          }
        pxT_z *= cInvDx;
        pxT_y *= cInvDx;

        sa(i, j, k, ix_) = pzT_y - pyT_z;
        sa(i, j, k, iy_) = pxT_z - pzT_x;
        sa(i, j, k, iz_) = pyT_x - pxT_y;
      });
    }
  };

  // q = curl(A) with A = p. Only non-zero on cells whose four (eight) corner
  // nodes are all inner, because p vanishes on every node touching an
  // uncovered cell. In 3D this calls the very operator the Faraday update
  // uses, so the correction is discretely a curl of the same stencil.
  auto apply_curl = [&]() {
    if (useZ) {
      curl_node_to_center(p, q, Geom(iLev).InvCellSize());
      return;
    }
    for (MFIter mfi(q); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& qa = q[mfi].array();
      const auto& pa = p[mfi].const_array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real pyAz =
            cInvDy * (pa(i, j + 1, k, iz_) + pa(i + 1, j + 1, k, iz_) -
                      pa(i, j, k, iz_) - pa(i + 1, j, k, iz_));
        const Real pxAz =
            cInvDx * (pa(i + 1, j, k, iz_) + pa(i + 1, j + 1, k, iz_) -
                      pa(i, j, k, iz_) - pa(i, j + 1, k, iz_));
        qa(i, j, k, ix_) = pyAz;
        qa(i, j, k, iy_) = -pxAz;
        qa(i, j, k, iz_) = 0.0;
      });
    }
  };

  // r = avg(Bf) - Bc on the covered cells, zero elsewhere. The components not
  // taking part in the relaxation stay zero, so the adjoint never drives their
  // potential components.
  r.setVal(0.0);
  for (MFIter mfi(B); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const auto& d = r[mfi].array();
    const auto& b = B[mfi].const_array();
    const auto& t = target[mfi].const_array();
    const auto& st = cell_status(iLev)[mfi].const_array();
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (bit::is_refined(st(i, j, k))) {
        for (int iVar = 0; iVar < nRelax; ++iVar)
          d(i, j, k, iVar) = t(i, j, k, iVar) - b(i, j, k, iVar);
      }
    });
  }

  // CGLS for min |r - curl(A)|^2: B directly accumulates the curl corrections
  // and r tracks the residual.
  apply_curlT();
  MultiFab::Copy(p, s, 0, 0, nDim3, 0);
  Real gamma = MultiFab::Dot(*owner, s, 0, s, 0, nDim3, 0);
  const Real gamma0 = gamma;

  for (int iter = 0; iter < maxIter && gamma > relTol * relTol * gamma0;
       ++iter) {
    apply_curl();
    const Real qq = MultiFab::Dot(q, 0, q, 0, nRelax, 0);
    if (qq <= 0)
      break;
    const Real a = gamma / qq;
    MultiFab::Saxpy(B, a, q, 0, 0, nRelax, 0);
    MultiFab::Saxpy(r, -a, q, 0, 0, nRelax, 0);
    apply_curlT();
    const Real gammaNew = MultiFab::Dot(*owner, s, 0, s, 0, nDim3, 0);
    MultiFab::Xpay(p, gammaNew / gamma, s, 0, 0, nDim3, 0);
    gamma = gammaNew;
  }
  B.FillBoundary(period);
}

//==========================================================
// Max |divB| per level, split into: cells next to the coarse-fine interface
// coarse cells, and the remaining interior cells. Also reports the drift
// between the covered coarse B and the average of the fine B, and the drift
// of the first fine ghost layer from the coarse interpolation. One line per
// cycle:
//   divB-AMR n=<cycle> L<lev> if/cv/in=<iface>/<covered>/<interior>/dom=<edge>
//            dCov=<max|Bc-avg(Bf)|> dGhost=<max|Bghost-interp(Bc)|>
//            |B|=<max|B|>
// (cv is omitted on the finest level; a non-finite divB shows as inf. 'dom'
// holds the cells on the physical domain boundary, which the field BC sets
// rather than the update, so they are not part of the interior bucket.)
void Pic::report_divB_amr() {
  BL_PROFILE("Pic::report_divB_amr");
  std::ostringstream line;
  line << std::scientific << std::setprecision(1);
  line << "divB-AMR n=" << tc->get_cycle();
  Real maxDrift = 0, maxGhostDrift = 0, maxBmag = 0;
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    if (divB[iLev].empty())
      continue;
    Real maxIface = 0, maxCovered = 0, maxInterior = 0, maxDomain = 0;
    const bool hasFiner = iLev < finestLevel();
    for (MFIter mfi(divB[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& arr = divB[iLev][mfi].array();
      const auto& status = cell_status(iLev)[mfi].array();
      amrex::Loop(box, [&](int i, int j, int k) {
        const int st = status(i, j, k);
        const Real raw = arr(i, j, k, 0);
        const Real v = std::isfinite(raw) ? std::abs(raw)
                                          : std::numeric_limits<Real>::max();
        // Cells on the physical domain boundary are set by the field BC, not
        // by the update, so they are reported separately: leaving them in the
        // interior bucket makes the solver's div(B) unreadable on decks with
        // outflow faces.
        if (bit::is_domain_edge(st)) {
          maxDomain = amrex::max(maxDomain, v);
        } else if (hasFiner && bit::is_refined(st)) {
          maxCovered = amrex::max(maxCovered, v);
        } else if ((iLev > 0 && bit::is_lev_edge(st)) ||
                   (hasFiner && bit::is_refined_neighbour(st))) {
          maxIface = amrex::max(maxIface, v);
        } else {
          maxInterior = amrex::max(maxInterior, v);
        }
      });
    }
    ParallelDescriptor::ReduceRealMax(maxIface);
    ParallelDescriptor::ReduceRealMax(maxCovered);
    ParallelDescriptor::ReduceRealMax(maxInterior);
    ParallelDescriptor::ReduceRealMax(maxDomain);

    // Sampled before the restriction in update_B_hybrid; measuring it here
    // instead would report zero whenever the plain average_down runs.
    const Real drift = (hasFiner && iLev < coveredBDrift.size())
                           ? coveredBDrift[iLev]
                           : Real(0.0);

    // Drift of the first fine ghost layer from the coarse interpolation.
    Real ghostDrift = 0;
    if (iLev > 0) {
      MultiFab tmp(centerB[iLev].boxArray(), centerB[iLev].DistributionMap(),
                   nDim3, centerB[iLev].nGrow());
      MultiFab::Copy(tmp, centerB[iLev], 0, 0, nDim3, centerB[iLev].nGrow());
      fill_fine_lev_bny_from_coarse(
          centerB[iLev - 1], tmp, 0, nDim3, refRatio(iLev - 1), Geom(iLev - 1),
          Geom(iLev), cell_status(iLev), *get_cell_interp());
      for (MFIter mfi(tmp); mfi.isValid(); ++mfi) {
        const Box box = amrex::grow(mfi.validbox(), 1);
        const auto& a = tmp[mfi].array();
        const auto& b = centerB[iLev][mfi].array();
        const auto& status = cell_status(iLev)[mfi].array();
        amrex::Loop(box, [&](int i, int j, int k) {
          const int st = status(i, j, k);
          if (!bit::is_lev_boundary(st) || bit::is_domain_boundary(st))
            return;
          for (int iVar = 0; iVar < nDim3; ++iVar)
            ghostDrift = amrex::max(
                ghostDrift, std::abs(a(i, j, k, iVar) - b(i, j, k, iVar)));
        });
      }
      ParallelDescriptor::ReduceRealMax(ghostDrift);
    }

    // Reference magnitude for the dGhost warning below.
    Real bMag = 0;
    for (int iVar = 0; iVar < nDim3; ++iVar)
      bMag = amrex::max(bMag, centerB[iLev].norm0(iVar, 0, false));

    auto fmt = [](Real v) {
      std::ostringstream s;
      s << std::scientific << std::setprecision(1);
      if (v >= std::numeric_limits<Real>::max())
        s << "inf";
      else
        s << v;
      return s.str();
    };
    line << " L" << iLev << (hasFiner ? " if/cv/in=" : " if/in=")
         << fmt(maxIface) << "/";
    if (hasFiner)
      line << fmt(maxCovered) << "/";
    line << fmt(maxInterior) << "/dom=" << fmt(maxDomain);
    maxDrift = amrex::max(maxDrift, drift);
    maxGhostDrift = amrex::max(maxGhostDrift, ghostDrift);
    maxBmag = amrex::max(maxBmag, bMag);
  }
  line << " dCov=" << maxDrift << " dGhost=" << maxGhostDrift
       << "|B|=" << maxBmag;

  // The first fine ghost layer is advanced by the Faraday update rather than
  // re-interpolated from the coarse level. Its difference from the coarse
  // interpolation typically saturates at a benign level without impacting
  // div(B). To avoid false positives on strong gradients, warn only when drift
  // exceeds a threshold ratio and exhibits continued growth past its previous
  // maximum.
  constexpr Real dGhostRunaway = 6.0; // threshold ratio above benign saturation
  constexpr Real dGrossGrowth = 2.0;  // factor that counts as sustained growth
  if (maxBmag > 0) {
    static Real maxRatioSeen = 0.0;
    const Real ratio = maxGhostDrift / maxBmag;
    if (ratio > dGhostRunaway && ratio > dGrossGrowth * maxRatioSeen) {
      amrex::Print()
          << printPrefix
          << "Warning: the evolved first fine ghost layer of B is "
             "running away from the coarse interpolation: dGhost/|B| = "
          << ratio << " (dGhost=" << maxGhostDrift << ", |B|=" << maxBmag
          << "), up from " << maxRatioSeen
          << " previously. Saturation near 5x|B| is benign; sustained "
             "growth is not.\n";
      maxRatioSeen = ratio;
    }
  }
  amrex::Print() << printPrefix << line.str() << std::endl;
}

//==========================================================

//==========================================================
void Pic::update_B_hybrid() {
  std::string nameFunc = "Pic::update_B_hybrid";
  timing_func(nameFunc);

  if (finestLevel() > 0 && nDim == 3 && !is_fake_2d()) {
    static bool isWarned = false;
    if (!isWarned) {
      amrex::Print()
          << printPrefix
          << "Note: AMR hybrid PIC on a 3D grid: the covered "
             "coarse B drifts from the average of the fine B, since the "
             "2dx-wide curl operator does not telescope between dx and "
             "dx/2 (divergence-free relaxation removes the drift each step). "
             "div(B) itself is preserved in 3D as well.\n";
      isWarned = true;
    }
  }

  const Real dt = tc->get_dt();
  const Real subDt = dt / nBSubcycle;

  // Grid-mode hyper-resistivity: every level is damped at its own grid scale,
  // eta_h = 4*pi*C_h*dx_min^4/dt with the dx of that level. The damping rate
  // of the grid-scale mode, C_h*(k dx)^4/dt, and the hyper-resistive CFL
  // number are then the same on every level. With the finest dx on all levels
  // the coarse Nyquist mode was damped 2^4 = 16 times less than the fine one,
  // and its noise reached the fine level through the coarse-fine interface E.
  if (etaHyperMode == "grid" && etaHyperCh > 0) {
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      const auto dxLev = Geom(iLev).CellSizeArray();
      Real dxMin = dxLev[0];
      for (int d = 1; d < nDim; ++d)
        dxMin = amrex::min(dxMin, dxLev[d]);
      etaHyperLev[iLev] = fourPI * etaHyperCh * std::pow(dxMin, 4) / dt;
    }
  }

  // Grid-mode localized regional hyper-resistivity
  if (hasRegionalHyper_) {
    update_regional_hyper_grid_mode(dt);
  }

  // CFL stability check for explicit resistive, hyper-resistive, and whistler
  // terms.
  const Real cflLimit = useRK4 ? 2.785 : 2.513;
  Real maxRegionalEta = 0.0;
  for (const auto& cfg : regionalResistivityConfigs) {
    maxRegionalEta = amrex::max(maxRegionalEta, cfg.etaCode);
  }
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    Real maxRegionalHyperLev = 0.0;
    for (const auto& cfg : regionalHyperResistivityConfigs) {
      if (iLev < static_cast<int>(cfg.etaLev.size())) {
        maxRegionalHyperLev = amrex::max(maxRegionalHyperLev, cfg.etaLev[iLev]);
      }
    }
    const Real maxEta = amrex::max(etaResistivity, maxRegionalEta);
    const Real maxHyper = amrex::max(etaHyperLev[iLev], maxRegionalHyperLev);
    if (maxEta <= 0 && maxHyper <= 0 && !useHallTerm && !useHybridUpwind)
      continue;

    const auto dx = Geom(iLev).CellSizeArray();
    const Box& domBox = Geom(iLev).Domain();
    Real sMax = 0.0, lMax = 0.0;
    for (int iDim = 0; iDim < nDim; ++iDim) {
      const int nCellDim = domBox.length(iDim);
      if (nCellDim < 2)
        continue;
      const Real invDx2 = 1.0 / (dx[iDim] * dx[iDim]);
      lMax += 4.0 * invDx2;
      if (nCellDim >= 3)
        sMax += invDx2;
    }
    const Real cflEta = (maxEta / fourPI) * sMax * subDt;
    const Real cflHyper = (maxHyper / fourPI) * sMax * lMax * subDt;
    if (cflEta > cflLimit)
      amrex::Print()
          << "  [CFL warning] resistivity: eta*kmax^2*dt_sub/(4pi) = " << cflEta
          << " (> " << cflLimit << ", explicit diffusion may be unstable)\n";
    if (cflHyper > cflLimit)
      amrex::Print()
          << "  [CFL warning] hyper-resistivity: eta_h*kmax^4*dt_sub/(4pi) = "
          << cflHyper << " (> " << cflLimit
          << ", explicit 4th-order diffusion may be unstable)\n";

    if (!useHallTerm && !useHybridUpwind)
      continue;

    Real bMaxLev = 0.0;
    const MultiFab& bForCfl = total_center_B(centerB[iLev], iLev);
    for (int d = 0; d < 3; ++d) {
      bMaxLev = amrex::max(bMaxLev, bForCfl.norm0(d, 0, false));
    }
    const Real rhoMinLev = nodePlasma[nSpecies][iLev].min(iRho_, 0, false);
    const Real rhoEff = amrex::max(rhoMinLev, rhoMinOhm);

    if (useHybridUpwind) {
      // Wave speeds only: the bulk-flow part of a_d is not known here.
      const Real cflUpwind =
          upwind_max_rate(iLev, 0.0, bMaxLev, rhoEff) * subDt;
      if (cflUpwind > cflLimit)
        amrex::Print() << "  [CFL warning] upwind dissipation: "
                          "sum_d 2*a_d/dx_d*dt_sub = "
                       << cflUpwind << " (> " << cflLimit
                       << ", explicit upwind term may be unstable)\n";
    }

    if (useHallTerm) {
      const Real cflWhistler = (sMax * bMaxLev / (fourPI * rhoEff)) * subDt;
      if (cflWhistler > cflLimit) {
        const int nSubReq =
            static_cast<int>(std::ceil(nBSubcycle * (cflWhistler / cflLimit)));
        amrex::Print()
            << "  [CFL warning] whistler wave: "
               "(kmax^2*Bmax)/(4pi*rho_eff)*dt_sub = "
            << cflWhistler << " (> " << cflLimit
            << ", whistler wave is linearly unstable in RK4! rho_min="
            << rhoMinLev << ", rho_eff=" << rhoEff << ", B_max=" << bMaxLev
            << "). Recommend nBSubcycle >= " << nSubReq
            << " or increasing rhoMinOhm.\n";
      }
    }
  }

  // Precalculate the ambipolar electric field E_ambi = -grad(p_e)/(e*n_e)
  // once per PIC timestep outside the magnetic subcycling steps.
  compute_ambipolar_E();

  // Maintain transverse force balance (J x B = grad(Pe)) in Faraday stages.
  // When ambiInStages was false, the uncompensated diamagnetic Hall field Ey
  // produced an artificial jump across coarse-fine boundaries, generating
  // unphysical horizontal dipole stripes of Bz in Harris reconnection.
  const bool ambiInStages = (electronTemperature > 0);

  // Ensure B ghost cells are fresh across all levels before subcycling.
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    apply_centerB_BC(iLev);
  }

  const Real invSubcycle = 1.0 / static_cast<Real>(nBSubcycle);

  // Each stage first evaluates curl(E) on ALL levels (faraday_stage_rhs, which
  // also synchronizes E across the coarse-fine interface), then updates every
  // level, then refreshes the ghost cells. nG = faraday_ngrow(iLev) is 1 on
  // refined levels: the first ghost layer is advanced with the same curl(E)
  // as the valid cells.
  for (int subStep = 0; subStep < nBSubcycle; ++subStep) {
    // Moment time-interpolation weights hstep for RK stages within the
    // sub-step.
    const Real g = static_cast<Real>(subStep) * invSubcycle;
    const Real hstepStart = g;
    const Real hstepHalf = g + 0.5 * invSubcycle;
    const Real hstepEnd = g + invSubcycle;

    if (useRK4) {
      // Classical RK4 on dB/dt = -curl(E_Ohm(B)) across all levels.
      const Real dtSixth = -subDt / 6.0;
      const Real dtThird = -subDt / 3.0;

      // Stages 1-3: trial state B_s = B_n - c_s * subDt * k_{s-1}.
      const Real cStage[3] = { 0.5, 0.5, 1.0 };
      const Real hStage[4] = { hstepStart, hstepHalf, hstepHalf, hstepEnd };
      for (int s = 0; s < 3; ++s) {
        faraday_stage_rhs(s == 0 ? centerB : centerBstage, s, hStage[s],
                          ambiInStages);
        for (int iLev = 0; iLev < n_lev(); ++iLev) {
          const int nG = faraday_ngrow(iLev);
          MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0,
                            -cStage[s] * subDt, kStage[iLev][s], 0, 0, nDim3,
                            nG);
        }
        for (int iLev = 0; iLev < n_lev(); ++iLev) {
          apply_centerB_BC(iLev, centerBstage[iLev]);
        }
      }

      // Stage 4 and the final RK4 combination.
      faraday_stage_rhs(centerBstage, 3, hStage[3], ambiInStages);
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        const int nG = faraday_ngrow(iLev);
        MultiFab::Saxpy(centerB[iLev], dtSixth, kStage[iLev][0], 0, 0, nDim3,
                        nG);
        MultiFab::Saxpy(centerB[iLev], dtThird, kStage[iLev][1], 0, 0, nDim3,
                        nG);
        MultiFab::Saxpy(centerB[iLev], dtThird, kStage[iLev][2], 0, 0, nDim3,
                        nG);
        MultiFab::Saxpy(centerB[iLev], dtSixth, kStage[iLev][3], 0, 0, nDim3,
                        nG);
      }
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        apply_centerB_BC(iLev);
      }
      continue;
    }

    if (!useRK4) {
      // Strong-stability-preserving RK3 (Shu-Osher) across all levels.
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        MultiFab::Copy(centerBstart[iLev], centerB[iLev], 0, 0, nDim3,
                       centerBstart[iLev].nGrow());
      }

      // Stage 1: B1 = B_n - subDt * curl(E(B_n))
      faraday_stage_rhs(centerB, 0, hstepStart, ambiInStages);
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        const int nG = faraday_ngrow(iLev);
        MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0, -subDt,
                          kStage[iLev][0], 0, 0, nDim3, nG);
      }
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        apply_centerB_BC(iLev, centerBstage[iLev]);
      }

      // Stage 2: B2 = (3/4)*B_n + (1/4)*(B1 - subDt * curl(E(B1)))
      faraday_stage_rhs(centerBstage, 1, hstepEnd, ambiInStages);
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        const int nG = faraday_ngrow(iLev);
        MultiFab::LinComb(centerBstage[iLev], 0.25, centerBstage[iLev], 0, 0.75,
                          centerBstart[iLev], 0, 0, nDim3, nG);
        MultiFab::Saxpy(centerBstage[iLev], -0.25 * subDt, kStage[iLev][1], 0,
                        0, nDim3, nG);
      }
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        apply_centerB_BC(iLev, centerBstage[iLev]);
      }

      // Stage 3: B^{n+1} = (1/3)*B_n + (2/3)*(B2 - subDt * curl(E(B2)))
      faraday_stage_rhs(centerBstage, 2, hstepHalf, ambiInStages);
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        const int nG = faraday_ngrow(iLev);
        MultiFab::LinComb(centerB[iLev], 2.0 / 3.0, centerBstage[iLev], 0,
                          1.0 / 3.0, centerBstart[iLev], 0, 0, nDim3, nG);
        MultiFab::Saxpy(centerB[iLev], (-2.0 / 3.0) * subDt, kStage[iLev][2], 0,
                        0, nDim3, nG);
      }
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        apply_centerB_BC(iLev);
      }
      continue;
    }
  }

  // Restriction. The covered coarse cells have already been advanced by the
  // coarse Faraday update, which keeps the coarse level divergence-consistent;
  // average_down of the cell-centered B would not (its face weights are
  // (F0+2F1+F2)/4 instead of the CT value (F0+F2)/2).
  // Left alone, however, the covered coarse B drifts away from the fine
  // average (the fine-scale part of E is invisible to the coarse curl), so it
  // is relaxed toward it with divergence-free corrections only.
  if (finestLevel() > 0) {
    // How far the covered coarse B has drifted from the fine average. This is
    // the quantity the relaxation has to remove, so it must be sampled before
    // the relaxation updates the covered cells.
    for (int iLev = finestLevel() - 1; iLev >= 0; iLev--) {
      measure_covered_B_drift(iLev);
    }
    for (int iLev = finestLevel() - 1; iLev >= 0; iLev--) {
      relax_covered_B_to_fine(iLev);
    }
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    apply_centerB_BC(iLev);
  }

  // div(B) control, which hybrid otherwise has none: correct_B() is only
  // reachable through update_B(), which is gated by solveEM.
  if (useHyperbolicCleaning) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      average_center_to_node(centerB[iLev], nodeB[iLev]);
      nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
      if (iLev == 0) {
        apply_field_bc(node_status(iLev), nodeB[iLev], 0, nDim3,
                       &Pic::get_node_B, iLev, true);
      }
      ensure_divB(iLev);
      div_node_to_center(nodeB[iLev], divB[iLev], Geom(iLev).InvCellSize());
      correct_B(iLev);
      centerB[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }

  // Update nodeB from centerB on all levels.
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    average_center_to_node(centerB[iLev], nodeB[iLev]);
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(node_status(iLev), nodeB[iLev], 0, nDim3, &Pic::get_node_B,
                     iLev, true);
    }
    if (is_body_conducting()) {
      project_body_B(nodeB[iLev], iLev);
      nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }

  // Fill coarse-fine interface ghost cells for nodeB.
  if (finestLevel() > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(
          nodeB[iLev - 1], nodeB[iLev], 0, nDim3, refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    }
  }

  // The 'divB' diagnostic, which hybrid has no other place to fill.
  if (need_divB()) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      compute_divB(iLev);
    }
    // The .out writer samples a cell-centered variable at the cell whose low
    // corner is the output node, so the nodes on the high side of a fine box
    // read a ghost cell. Fill those from the coarse level, so the plot shows
    // the coarse value there instead of the div of interpolated ghost B.
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(
          divB[iLev - 1], divB[iLev], 0, divB[iLev].nComp(), refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), cell_status(iLev), cell_bilinear_interp);
      fill_fine_lev_bny_from_coarse(
          centerDivB[iLev - 1], centerDivB[iLev], 0, centerDivB[iLev].nComp(),
          refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          cell_bilinear_interp);
    }
    if (need_divB() && n_lev() > 1 && doReport)
      report_divB_amr();
  }

  // Evaluate E^{n+1} into nodeE for the next push.
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    assemble_ohm_E(centerB[iLev], total_center_B(centerB[iLev], iLev),
                   nodeE[iLev], iLev, 1.0);
  }

  // Suppress grid-scale E component if enabled. This has to happen *before*
  // the coarse-fine synchronization: smoothing afterwards would leave the two
  // levels with different E at the interface nodes again.
  if (doSmoothE) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
      smooth_E(nodeE[iLev], iLev);
    }
  }

  if (finestLevel() > 0) {
    // Smoothing also changes the wall and body nodes; restore their
    // constraints on every level before sync_emf_at_interface copies the fine
    // E onto the coarse interface nodes.
    if (doSmoothE) {
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
        apply_field_bc(node_status(iLev), nodeE[iLev], 0, nDim3,
                       &Pic::get_node_E, iLev, false);
        if (use_body()) {
          apply_body_E_bc(nodeE[iLev], iLev);
        }
      }
    }
    sync_emf_at_interface(nodeE);
  }

  // Fill coarse-fine interface ghost cells for nodeE.
  if (finestLevel() > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(
          nodeE[iLev - 1], nodeE[iLev], 0, nDim3, refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    }
  }
}

//==========================================================
// The first-order upwind term damps the grid-scale mode of direction d at
// 2 a_d / dx_d (a real eigenvalue), so the explicit field advance needs
// dt_sub * sum_d 2 a_d / dx_d below the real-axis stability limit of the RK
// scheme. Limited reconstructions only damp less.
Real Pic::upwind_max_rate(int iLev, Real uMax, Real bMax, Real rhoEff) const {
  if (!useHybridUpwind)
    return 0.0;

  Real cFast = 0.0, cWhistlerDx = 0.0;
  if (rhoEff > 0.0) {
    Real cs2 = 0.0;
    if (electronTemperature > 0.0)
      cs2 = ((electronGamma > 0.0) ? electronGamma : 1.0) * electronTemperature;
    cFast = std::sqrt(bMax * bMax / (fourPI * rhoEff) + cs2);
    if (useHallTerm)
      cWhistlerDx = upwindCoefWhistler * 0.25 * bMax / rhoEff;
  }

  Real rate = 0.0;
  for (int iDim = 0; iDim < nDim; ++iDim) {
    if (Geom(iLev).Domain().length(iDim) < 2)
      continue;
    const Real invDx = Geom(iLev).InvCellSize(iDim);
    const Real a = (upwindFixedSpeed > 0.0)
                       ? upwindFixedSpeed
                       : uMax + cFast + cWhistlerDx * invDx;
    rate += 2.0 * a * invDx;
  }
  return rate;
}

//==========================================================
void Pic::calc_hybrid_dt_and_subcycle(Real& dtNext, int& nSubNext,
                                      bool doReport) {
  std::string nameFunc = "Pic::calc_hybrid_dt_and_subcycle";
  timing_func(nameFunc);

  const Real cflLimit = useRK4 ? 2.785 : 2.513;
  const Real fieldSafety = 0.85;

  Real dtMacroMin = std::numeric_limits<Real>::max();
  Real dtSubMin = std::numeric_limits<Real>::max();
  std::string localMacroReason = "ion kinetics";
  std::string localSubReason = "field advance";

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    const auto dx = Geom(iLev).CellSizeArray();
    Real dxMin = dx[0];
    for (int d = 1; d < nDim; ++d) {
      dxMin = amrex::min(dxMin, dx[d]);
    }

    // 1. Kinetic ion advection limit (including bulk flow + thermal spread)
    Real uMaxLev = 0.0;
    if (fixedUMax >= 0) {
      uMaxLev = fixedUMax;
    } else {
      Vector<Real> uMaxSpecies(nSpecies, 0.0);
      for (int i = 0; i < nSpecies; ++i) {
        amrex::MultiFab& momMF = nodePlasma[i][iLev];
        auto const& ma = momMF.const_arrays();
        Real uSpMax = ParReduce(
            TypeList<ReduceOpMax>{}, TypeList<Real>{}, momMF, IntVect(0),
            [=] AMREX_GPU_DEVICE(int nb, int i, int j, int k)
                noexcept -> GpuTuple<Real> {
                  Array4<Real const> const& arr = ma[nb];
                  Real rho = arr(i, j, k, iRho_);
                  if (rho <= 0.0)
                    return 0.0;
                  Real invRho = 1.0 / rho;
                  Real uBulk2 = (arr(i, j, k, iUx_) * arr(i, j, k, iUx_) +
                                 arr(i, j, k, iUy_) * arr(i, j, k, iUy_) +
                                 arr(i, j, k, iUz_) * arr(i, j, k, iUz_)) *
                                (invRho * invRho);
                  Real p = (arr(i, j, k, iPxx_) + arr(i, j, k, iPyy_) +
                            arr(i, j, k, iPzz_)) *
                           (1.0 / 3.0);
                  Real vth = (p > 0.0) ? sqrt(p * invRho) : 0.0;
                  return sqrt(uBulk2) + 3.0 * vth;
                });
        uMaxSpecies[i] = uSpMax;
      }
      ParallelDescriptor::ReduceRealMax(uMaxSpecies.data(), nSpecies);
      for (int i = 0; i < nSpecies; ++i) {
        if (uMaxSpecies[i] > uMaxLev) {
          uMaxLev = uMaxSpecies[i];
        }
      }
    }

    // 2. Maximum magnetic field magnitude on this level (including intrinsic
    // B0)
    const MultiFab& bForCfl = total_center_B(centerB[iLev], iLev);
    auto const& bArrays = bForCfl.const_arrays();
    Real bMaxLev = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, bForCfl,
                             IntVect(0),
                             [=] AMREX_GPU_DEVICE(int nb, int i, int j, int k)
                                 noexcept -> GpuTuple<Real> {
                                   Array4<Real const> const& arr = bArrays[nb];
                                   Real b2 = arr(i, j, k, 0) * arr(i, j, k, 0) +
                                             arr(i, j, k, 1) * arr(i, j, k, 1) +
                                             arr(i, j, k, 2) * arr(i, j, k, 2);
                                   return sqrt(b2);
                                 });
    ParallelDescriptor::ReduceRealMax(bMaxLev);

    // Kinetic advection limit
    Real dtPart =
        (uMaxLev > 0.0) ? (dxMin / uMaxLev) : std::numeric_limits<Real>::max();

    // Ion gyrofrequency limit (Omega_ci * dt <= thetaGyro)
    Real dtGyro = std::numeric_limits<Real>::max();
    if (bMaxLev > 0.0) {
      Real maxQom = 0.0;
      for (int i = 0; i < nSpecies; ++i) {
        Real qom = std::abs(parts[i]->get_charge() / parts[i]->get_mass());
        maxQom = amrex::max(maxQom, qom);
      }
      if (maxQom > 0.0) {
        const Real thetaGyro = 0.35; // radians (~18 steps per gyroperiod)
        dtGyro = thetaGyro / (maxQom * bMaxLev);
      }
    }

    Real dtMacroLev = amrex::min(dtPart, dtGyro);
    if (dtMacroLev < dtMacroMin) {
      dtMacroMin = dtMacroLev;
      if (dtGyro < dtPart) {
        localMacroReason = "ion gyro-frequency";
      } else {
        localMacroReason = "ion particle CFL";
      }
    }

    // 3. Field subcycling limits
    const Box& domBox = Geom(iLev).Domain();
    Real sMax = 0.0, lMax = 0.0;
    for (int iDim = 0; iDim < nDim; ++iDim) {
      const int nCellDim = domBox.length(iDim);
      if (nCellDim < 2)
        continue;
      const Real invDx2 = 1.0 / (dx[iDim] * dx[iDim]);
      lMax += 4.0 * invDx2;
      if (nCellDim >= 3)
        sMax += invDx2;
    }

    if (sMax > 0.0) {
      const Real rhoMinLev = nodePlasma[nSpecies][iLev].min(iRho_, 0, false);
      const Real rhoEff = amrex::max(
          rhoMinLev, amrex::max(rhoMinOhm, static_cast<Real>(1e-10)));

      if (useHallTerm && bMaxLev > 0.0) {
        const Real omegaWhistler = (sMax * bMaxLev) / (fourPI * rhoEff);
        if (omegaWhistler > 0.0) {
          Real dtWhistler = (cflLimit * fieldSafety) / omegaWhistler;
          if (dtWhistler < dtSubMin) {
            dtSubMin = dtWhistler;
            localSubReason = "whistler wave";
          }
        }
      }

      // Fast magnetosonic / Alfvén wave limit (governs when useHallTerm is
      // false or on coarse cells)
      if (bMaxLev > 0.0 && rhoEff > 0.0) {
        Real vA = bMaxLev / std::sqrt(fourPI * rhoEff);
        Real cs = 0.0;
        if (electronTemperature > 0.0) {
          Real gammaE = (electronGamma > 0.0) ? electronGamma : 1.0;
          cs = std::sqrt(gammaE * electronTemperature);
        }
        Real vFast = std::sqrt(vA * vA + cs * cs);
        if (vFast > 0.0) {
          Real dtWave = (cflLimit * fieldSafety) * dxMin / vFast;
          if (dtWave < dtSubMin) {
            dtSubMin = dtWave;
            localSubReason =
                useHallTerm ? "Alfvén / fast wave" : "ideal-MHD wave";
          }
        }
      }

      if (etaResistivity > 0.0) {
        const Real omegaEta = (etaResistivity / fourPI) * sMax;
        if (omegaEta > 0.0) {
          Real dtEta = (cflLimit * fieldSafety) / omegaEta;
          if (dtEta < dtSubMin) {
            dtSubMin = dtEta;
            localSubReason = "resistive diffusion";
          }
        }
      }

      if (etaHyperLev.size() > iLev && etaHyperLev[iLev] > 0.0 && lMax > 0.0) {
        const Real omegaHyper = (etaHyperLev[iLev] / fourPI) * sMax * lMax;
        if (omegaHyper > 0.0) {
          Real dtHyper = (cflLimit * fieldSafety) / omegaHyper;
          if (dtHyper < dtSubMin) {
            dtSubMin = dtHyper;
            localSubReason = "hyper-resistivity";
          }
        }
      }

      if (useHybridUpwind) {
        const Real omegaUpwind =
            upwind_max_rate(iLev, uMaxLev, bMaxLev, rhoEff);
        if (omegaUpwind > 0.0) {
          Real dtUpwind = (cflLimit * fieldSafety) / omegaUpwind;
          if (dtUpwind < dtSubMin) {
            dtSubMin = dtUpwind;
            localSubReason = "upwind dissipation";
          }
        }
      }
    }
  }

  const Real userCFL = (tc->get_cfl() > 0.0) ? tc->get_cfl() : 0.2;
  const Real dtMacroTarget = userCFL * dtMacroMin;
  dtLimitingReason = localMacroReason;

  if (tc->get_cfl() <= 0.0) {
    // Fixed macro dt mode (e.g. SWMF coupling or user-fixed dt)
    Real dtFixed = tc->get_dt();
    if (dtFixed <= 0.0) {
      dtFixed = tc->get_next_dt();
    }
    dtNext = dtFixed;

    if (isAutoSubcycle && dtSubMin < std::numeric_limits<Real>::max()) {
      int nSubReq = static_cast<int>(std::ceil(dtFixed / dtSubMin));
      nSubReq = std::clamp(nSubReq, nBSubcycleMin, nBSubcycleMax);

      if (nSubReq > nBSubcycle) {
        // Immediate increase for numerical stability
        nSubNext = nSubReq;
        subcycleHoldCount = 0;
      } else if (nSubReq < nBSubcycle) {
        // Hysteresis with 15% safety buffer before decreasing
        if (dtFixed / (nBSubcycle - 1) <= 0.85 * dtSubMin) {
          subcycleHoldCount++;
          if (subcycleHoldCount >= 3) {
            nSubNext = nBSubcycle - 1;
            subcycleHoldCount = 0;
          } else {
            nSubNext = nBSubcycle;
          }
        } else {
          subcycleHoldCount = 0;
          nSubNext = nBSubcycle;
        }
      } else {
        nSubNext = nBSubcycle;
        subcycleHoldCount = 0;
      }
    } else {
      nSubNext = nBSubcycle;
    }
    dtLimitingReason =
        isAutoSubcycle ? (localSubReason + " (subcycling)") : "fixed dt";
  } else {
    // Adaptive macro dt mode
    if (!isAutoSubcycle) {
      Real dtFieldLimit = nBSubcycle * dtSubMin;
      Real dtTarget = std::min(dtMacroTarget, dtFieldLimit);
      if (dtTarget == dtFieldLimit &&
          dtSubMin < std::numeric_limits<Real>::max()) {
        dtLimitingReason = localSubReason + " (fixed subcycles)";
      }

      Real dtOld = tc->get_dt();
      if (dtOld > 0.0) {
        // Allow stability limits to force an immediate drop; limit growth to
        // 1.10x.
        dtNext = std::min(dtTarget, static_cast<Real>(1.10) * dtOld);
      } else {
        dtNext = dtTarget;
      }
      nSubNext = nBSubcycle;
    } else {
      int nIdeal = (dtSubMin < std::numeric_limits<Real>::max())
                       ? static_cast<int>(std::ceil(dtMacroTarget / dtSubMin))
                       : 1;
      Real dtTarget = dtMacroTarget;

      if (nIdeal > nBSubcycleMax &&
          dtSubMin < std::numeric_limits<Real>::max()) {
        nSubNext = nBSubcycleMax;
        dtTarget = nBSubcycleMax * dtSubMin;
        dtLimitingReason = localSubReason + " (capped at nBSubcycleMax)";
        subcycleHoldCount = 0;
      } else {
        int nTarget = std::clamp(nIdeal, nBSubcycleMin, nBSubcycleMax);
        if (nTarget > nBSubcycle) {
          // Immediate increase for numerical stability
          nSubNext = nTarget;
          subcycleHoldCount = 0;
        } else if (nTarget < nBSubcycle) {
          if (dtTarget / (nBSubcycle - 1) <= 0.85 * dtSubMin) {
            subcycleHoldCount++;
            if (subcycleHoldCount >= 3) {
              nSubNext = nBSubcycle - 1;
              subcycleHoldCount = 0;
            } else {
              nSubNext = nBSubcycle;
            }
          } else {
            subcycleHoldCount = 0;
            nSubNext = nBSubcycle;
          }
        } else {
          nSubNext = nBSubcycle;
          subcycleHoldCount = 0;
        }
      }

      Real dtOld = tc->get_dt();
      if (dtOld > 0.0) {
        // Allow stability limits to force an immediate drop; limit growth to
        // 1.10x.
        dtNext = std::min(dtTarget, static_cast<Real>(1.10) * dtOld);
      } else {
        dtNext = dtTarget;
      }
    }
  }

  if (dtMacroMin > 0.0 && dtMacroMin < std::numeric_limits<Real>::max()) {
    maxCFL = dtNext / dtMacroMin;
  }

  if (doReport && !domainParameters.doCompact) {
    amrex::Print() << printPrefix << "[Hybrid Adaptive Step] dt = " << dtNext
                   << ", nSub = " << nSubNext
                   << " (limited by: " << dtLimitingReason << ")\n";
  }
}
