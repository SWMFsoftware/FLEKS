#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_PhysBCFunct.H>

#include "GridUtility.h"
#include "Pic.h"
#include "PicUpwind.h"
#include "Timer.h"

using namespace amrex;

namespace {

//==========================================================
// Helpers for the evolved electron pressure.

// Face-averaged component `comp` of a nodal vector field on the low (off=0) or
// high (off=1) face of cell (i,j,k) along direction DIR. Averages the 4 corner
// nodes in 3D or the 2 edge nodes in 2D / fake-2D, matching div_node_to_center.
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real ue_face_comp(
    const Array4<Real const>& arrUe, const int comp, const int i, const int j,
    const int k, const int off, const bool isFake2D) noexcept {
  const int kOff = (nDim > 2 && !isFake2D) ? 1 : 0;
  if constexpr (DIR == ix_) {
    if (kOff > 0) {
      return 0.25 *
             (arrUe(i + off, j, k, comp) + arrUe(i + off, j + 1, k, comp) +
              arrUe(i + off, j, k + 1, comp) +
              arrUe(i + off, j + 1, k + 1, comp));
    }
    return 0.5 * (arrUe(i + off, j, k, comp) + arrUe(i + off, j + 1, k, comp));
  } else if constexpr (DIR == iy_) {
    if (kOff > 0) {
      return 0.25 *
             (arrUe(i, j + off, k, comp) + arrUe(i + 1, j + off, k, comp) +
              arrUe(i, j + off, k + 1, comp) +
              arrUe(i + 1, j + off, k + 1, comp));
    }
    return 0.5 * (arrUe(i, j + off, k, comp) + arrUe(i + 1, j + off, k, comp));
  } else {
    return 0.25 *
           (arrUe(i, j, k + off, comp) + arrUe(i + 1, j, k + off, comp) +
            arrUe(i, j + 1, k + off, comp) +
            arrUe(i + 1, j + 1, k + off, comp));
  }
}

template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
ue_face(const Array4<Real const>& arrUe, const int i, const int j, const int k,
        const int off, const bool isFake2D) noexcept {
  return ue_face_comp<DIR>(arrUe, DIR, i, j, k, off, isFake2D);
}

// Smooth, positivity-preserving Le et al. (2009) trapped-passing equations of
// state for collisionless magnetized electrons:
//   Pe_par  = p0 * n~ * (1 + x^3)^(1/3)
//   Pe_perp = p0 * n~ * (1 + B~ * x^2) / (1 + x^2)
// where n~ = n/n0, B~ = sqrt(B^2 + B_min^2)/B0, and x = (pi/6) * (n~/B~)^2.
// Asymptotes to Boltzmann (Pe_par = Pe_perp = p0*n~) as x -> 0 and trapped
// CGL (Pe_par = p0*(pi/6)*n~^3/B~^2, Pe_perp = p0*n~*B~) as x -> inf.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void le2009_closure(
    const Real rho, const Real bMag, const Real p0, const Real invRho0,
    const Real invB0, const Real bMin, const Real peMin,
    const bool useFirehoseLimit, const Real firehoseCoef, Real& peParOut,
    Real& pePerpOut) noexcept {
  if (rho <= 0.0 || p0 <= 0.0 || invRho0 <= 0.0 || invB0 <= 0.0) {
    peParOut = peMin;
    pePerpOut = peMin;
    return;
  }
  const Real nTilde = rho * invRho0;
  const Real bEff = std::sqrt(bMag * bMag + bMin * bMin);
  const Real bTilde = bEff * invB0;
  const Real ratio = nTilde / bTilde;
  const Real x = (dPI / 6.0) * ratio * ratio;
  const Real x2 = x * x;
  const Real x3 = x2 * x;
  Real pPar = p0 * nTilde * std::cbrt(1.0 + x3);
  Real pPerp = p0 * nTilde * (1.0 + bTilde * x2) / (1.0 + x2);
  if (useFirehoseLimit && bMag > 0.0) {
    const Real pB = (bMag * bMag) / (8.0 * dPI);
    const Real pTot = (pPar + 2.0 * pPerp) / 3.0;
    const Real pParMax = pTot + (2.0 / 3.0) * firehoseCoef * pB;
    if (pPar > pParMax) {
      pPar = pParMax;
      pPerp = 0.5 * (3.0 * pTot - pPar);
    }
  }
  peParOut = amrex::max(pPar, peMin);
  pePerpOut = amrex::max(pPerp, peMin);
}

// Add one direction's contribution to div(u_e P) and div(u_e). The caller
// skips a direction that carries no flux by testing invDx > 0 first.
template <int DIR>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void pe_dir_flux(
    const Array4<Real const>& arrPe, const int comp,
    const Array4<Real const>& arrUe, const int i, const int j, const int k,
    const Real invDx, const int limType, const Real peMin, const bool isFake2D,
    Real& divFlux, Real& divUe) noexcept {
  const Real uLo = ue_face<DIR>(arrUe, i, j, k, 0, isFake2D);
  const Real uHi = ue_face<DIR>(arrUe, i, j, k, 1, isFake2D);

  auto p_at = [&](const int off) noexcept -> Real {
    if constexpr (DIR == ix_) {
      return arrPe(i + off, j, k, comp);
    } else if constexpr (DIR == iy_) {
      return arrPe(i, j + off, k, comp);
    } else {
      return arrPe(i, j, k + off, comp);
    }
  };

  Real peLo;
  if (limType == 0) {
    peLo = (uLo > 0.0) ? p_at(-1) : p_at(0);
  } else {
    const Real sLoL =
        pe_slope(limType, p_at(-1) - p_at(-2), p_at(0) - p_at(-1));
    const Real sLoR = pe_slope(limType, p_at(0) - p_at(-1), p_at(1) - p_at(0));
    const Real pL = amrex::max(p_at(-1) + 0.5 * sLoL, peMin);
    const Real pR = amrex::max(p_at(0) - 0.5 * sLoR, peMin);
    peLo = (uLo > 0.0) ? pL : ((uLo < 0.0) ? pR : 0.5 * (pL + pR));
  }

  Real peHi;
  if (limType == 0) {
    peHi = (uHi > 0.0) ? p_at(0) : p_at(1);
  } else {
    const Real sHiL = pe_slope(limType, p_at(0) - p_at(-1), p_at(1) - p_at(0));
    const Real sHiR = pe_slope(limType, p_at(1) - p_at(0), p_at(2) - p_at(1));
    const Real pL = amrex::max(p_at(0) + 0.5 * sHiL, peMin);
    const Real pR = amrex::max(p_at(1) - 0.5 * sHiR, peMin);
    peHi = (uHi > 0.0) ? pL : ((uHi < 0.0) ? pR : 0.5 * (pL + pR));
  }

  divFlux += (uHi * peHi - uLo * peLo) * invDx;
  divUe += (uHi - uLo) * invDx;
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

// True when node (i,j,k) is a corner of at least one active (in-level,
// uncovered, outside-body) cell on this level.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool touches_active_cell(
    const Array4<int const>& cStat, const int i, const int j, const int k,
    const bool hasBody, const bool hasFiner) noexcept {
  const int kmin = (nDim > 2) ? k - 1 : k;
  const int jmin = (nDim > 1) ? j - 1 : j;
  for (int kk = kmin; kk <= k; ++kk) {
    for (int jj = jmin; jj <= j; ++jj) {
      for (int ii = i - 1; ii <= i; ++ii) {
        const int st = cStat(ii, jj, kk);
        if (!bit::is_lev_boundary(st) && !(hasBody && bit::is_body(st)) &&
            !(hasFiner && bit::is_refined(st))) {
          return true;
        }
      }
    }
  }
  return false;
}

// True when node (i,j,k) lies outside #BODY and touches at least one active
// (in-level, uncovered, outside-body) cell on this level.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool is_active_cfl_node(
    const Array4<int const>& cStat, const Array4<int const>& nStat, const int i,
    const int j, const int k, const bool hasBody,
    const bool hasFiner) noexcept {
  if (hasBody && bit::is_body(nStat(i, j, k)))
    return false;
  return touches_active_cell(cStat, i, j, k, hasBody, hasFiner);
}

// Minimum positive stage density rho_q(hstep) encountered by assemble_ohm_E
// across all stages of a hybrid step. Returns numeric_limits<Real>::max() when
// every active stage has rho_q <= 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
min_positive_stage_rho(const Real rCur, const Real rPrev, const bool hasPrev,
                       const bool isPC, const Real minSubStageWeight) noexcept {
  Real rMin = (rCur > 0.0) ? rCur : std::numeric_limits<Real>::max();
  if (!hasPrev)
    return rMin;
  const Real rMid = 0.5 * (rPrev + rCur);
  if (rMid > 0.0)
    rMin = amrex::min(rMin, rMid);
  if (isPC) {
    if (rPrev > 0.0)
      rMin = amrex::min(rMin, rPrev);
    if ((rPrev > 0.0) != (rCur > 0.0)) {
      const Real rEdge = minSubStageWeight * amrex::max(rPrev, rCur);
      if (rEdge > 0.0)
        rMin = amrex::min(rMin, rEdge);
    }
  }
  return rMin;
}

} // namespace

//==========================================================
void Pic::assemble_ohm_E(const MultiFab& centerBin,
                         const MultiFab& centerBtimeAvg, MultiFab& Eout,
                         int iLev, Real hstep, bool includeAmbi,
                         Real timeOverride) {
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
                   iLev, true, 0.0, timeOverride);
  }
  if (is_body_conducting()) {
    project_body_B(nodeBstage[iLev], iLev);
    nodeBstage[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  // Moment time-interpolation weights:
  // - Charge density rho_q is deposited from particle positions x^{n+1} at
  //   integer time t^{n+1} (in nodePlasma) and x^n at t^n (in nodePlasmaPrev),
  //   so rho_q(hstep) = (1 - hstep)*rho_q^n + hstep*rho_q^{n+1} for both
  //   single-pass and predictor-corrector modes (strictly convex for hstep in
  //   [0, 1], preventing negative density extrapolation).
  // - In single-pass mode (!usePredictorCorrector), ion momentum m = rho_q u_i
  //   is deposited from velocities v^{n+1/2} at t^{n+1/2} (hstep = 0.5) in
  //   nodePlasma and v^{n-1/2} at t^n - 0.5*dtPrev in nodePlasmaPrev (or v^0
  //   at t^0 on the first step). The normalized step separation is dhM =
  //   dtVel / dt, giving m(hstep) = ((0.5 - hstep)/dhM)*m^{prev} +
  //   (1 - (0.5 - hstep)/dhM)*m^{n+1/2}. Faraday sub-steps pass hstep = 0.5
  //   (wPrevM = 0, wCurM = 1) to hold the midpoint current constant across
  //   [t^n, t^{n+1}], avoiding the secant-ramp whistler instability at
  //   omega_w * dt ~ pi, while the final nodeE evaluation at t^{n+1} passes
  //   hstep = 1.0.
  const Real wPrevRho = 1.0 - hstep;
  const Real wCurRho = hstep;
  const Real dtCur = tc ? tc->get_dt() : 0.0;
  const Real dhM = (isFirstHybridStep || dtPrevHybrid <= 0.0 || dtCur <= 0.0)
                       ? 0.5
                       : (0.5 * (dtPrevHybrid + dtCur) / dtCur);
  const Real wPrevM =
      usePredictorCorrector ? (1.0 - hstep) : ((0.5 - hstep) / dhM);
  const Real wCurM = usePredictorCorrector ? hstep : (1.0 - wPrevM);
  const bool usePrevMoments = (wPrevRho != 0.0 || wPrevM != 0.0);
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
  const Real teUp = (useElectronPressureEq && peStateInitialized_)
                        ? max_evolved_Te(iLev)
                        : electronTemperature;
  const Real cs2Up =
      (teUp > 0) ? ((electronGamma > 0) ? electronGamma : 1.0) * teUp : 0.0;
  const bool useHall = useHallTerm;

  {
    struct UpwindScratch {
      FArrayBox faceJump[nDim3];
      FArrayBox slopeFab;
    };
    std::optional<UpwindScratch> upwindScratch;
    if (doUpwind)
      upwindScratch.emplace();
    for (MFIter mfi(Eout, stream_ordered_mfi()); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrE = Eout[mfi].array();
      const Array4<Real const>& arrB = nodeBstage[iLev][mfi].array();
      const Array4<Real const>& moments =
          nodePlasma[nSpecies][iLev][mfi].array();
      const Array4<Real const> momentsPrev =
          usePrevMoments ? nodePlasmaPrev[nSpecies][iLev][mfi].array()
                         : Array4<Real const>();
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
      Array4<Real const> arrFace[nDim3];
      if (doUpwind) {
        FArrayBox(&faceJump)[nDim3] = upwindScratch->faceJump;
        FArrayBox& slopeFab = upwindScratch->slopeFab;
        const Array4<Real const>& arrBin = centerBin[mfi].const_array();
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
            parallel_for_omp(
                slopeBox, nDim3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                  arrS(i, j, k, iVar) =
                      upwind_cell_slope<ix_>(arrBin, i, j, k, iVar, limTypeUp);
                });
            parallel_for_omp(
                faceBox, nDim3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                  arrF(i, j, k, iVar) = upwind_face_jump<ix_>(
                      arrBin, arrSc, i, j, k, iVar, sensorFloorUp);
                });
          } else if (iDir == iy_) {
            parallel_for_omp(
                slopeBox, nDim3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                  arrS(i, j, k, iVar) =
                      upwind_cell_slope<iy_>(arrBin, i, j, k, iVar, limTypeUp);
                });
            parallel_for_omp(
                faceBox, nDim3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                  arrF(i, j, k, iVar) = upwind_face_jump<iy_>(
                      arrBin, arrSc, i, j, k, iVar, sensorFloorUp);
                });
          } else {
            parallel_for_omp(
                slopeBox, nDim3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                  arrS(i, j, k, iVar) =
                      upwind_cell_slope<iz_>(arrBin, i, j, k, iVar, limTypeUp);
                });
            parallel_for_omp(
                faceBox, nDim3,
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

      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const Real rhoCur = moments(i, j, k, iRho_);
        const Real rhoPrev = usePrevMoments ? momentsPrev(i, j, k, iRho_) : 0.0;
        Real rho = wCurRho * rhoCur;
        if (usePrevMoments) {
          rho += wPrevRho * rhoPrev;
        }
        Real ui = 0.0, vi = 0.0, wi = 0.0;
        const Real invRhoEff =
            (rho > 0.0) ? (1.0 / amrex::max(rho, rhoMinOhm)) : 0.0;

        if (usePredictorCorrector) {
          Real mx = wCurM * moments(i, j, k, iMx_);
          Real my = wCurM * moments(i, j, k, iMy_);
          Real mz = wCurM * moments(i, j, k, iMz_);
          if (usePrevMoments) {
            mx += wPrevM * momentsPrev(i, j, k, iMx_);
            my += wPrevM * momentsPrev(i, j, k, iMy_);
            mz += wPrevM * momentsPrev(i, j, k, iMz_);
          }
          if (rho > 0.0) {
            ui = mx * invRhoEff;
            vi = my * invRhoEff;
            wi = mz * invRhoEff;
          }
        } else {
          // Single-pass hybrid leapfrog:
          // moments(iMx_) carries momentum deposited at post-drift positions
          // x^{n+1} with half-step velocities v^{n+1/2}, i.e. rho_cur *
          // u^{n+1/2}. momentsPrev(iMx_) carries rho_prev * u^{n-1/2} at
          // positions x^n. Directly interpolating raw momentum and dividing by
          // time-interpolated density gives u_i = U * rho_next / ((rho_prev +
          // rho_next)/2) != U for nonuniform density advecting at uniform U,
          // introducing an O(dt) error. Instead, reconstruct ion fluid velocity
          // u_i at target time hstep:
          Real uiCur = 0.0, viCur = 0.0, wiCur = 0.0;
          if (rhoCur > 0.0) {
            const Real invRhoCur = 1.0 / amrex::max(rhoCur, rhoMinOhm);
            uiCur = moments(i, j, k, iMx_) * invRhoCur;
            viCur = moments(i, j, k, iMy_) * invRhoCur;
            wiCur = moments(i, j, k, iMz_) * invRhoCur;
          }
          ui = uiCur;
          vi = viCur;
          wi = wiCur;
          if (usePrevMoments && wPrevM != 0.0) {
            Real uiPrev = 0.0, viPrev = 0.0, wiPrev = 0.0;
            if (rhoPrev > 0.0) {
              const Real invRhoPrev = 1.0 / amrex::max(rhoPrev, rhoMinOhm);
              uiPrev = momentsPrev(i, j, k, iMx_) * invRhoPrev;
              viPrev = momentsPrev(i, j, k, iMy_) * invRhoPrev;
              wiPrev = momentsPrev(i, j, k, iMz_) * invRhoPrev;
            }
            ui = wCurM * uiCur + wPrevM * uiPrev;
            vi = wCurM * viCur + wPrevM * viPrev;
            wi = wCurM * wiCur + wPrevM * wiPrev;
          }
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
            const Real a = (fixedSpeedUp > 0)
                               ? fixedSpeedUp
                               : std::abs(ue[iDir]) + cFast +
                                     coefWhistlerUp * 0.25 * bAbs * invRhoEff *
                                         invDxUp[iDir];
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
    if (doUpwind)
      sync_before_free();
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
                   centerLapB[iLev].nComp(), &Pic::get_zero, iLev, true, 0.0);

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
                   nodeHyperE[iLev].nComp(), &Pic::get_zero, iLev, false, 0.0);

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

        parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
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
                 false, 0.0, timeOverride);
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

  if (useAnisoPe || useLe2009Pe) {
    if (useLe2009Pe) {
      compute_electron_density(iLev, nodeRhoTemp[iLev], centerPe[iLev], false);
    }

    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
    const MultiFab& bTot = total_center_B(centerB[iLev], iLev);
    const Real p0 = electronDensity0 * electronTemperature;
    const Real invRho0 =
        (electronDensity0 > 0.0) ? (1.0 / electronDensity0) : 0.0;
    const Real invB0 = (le2009B0 > 0.0) ? (1.0 / le2009B0) : 0.0;
    const Real bMin = fieldAlignedBMin;
    const Real bMin2 = bMin * bMin;
    const Real peMinLocal = peMin;
    const bool useFh = usePeInstabilityLimit;
    const Real fhCoef = peFirehoseCoef;
    const bool isLe = useLe2009Pe;

    for (MFIter mfi(centerPeTensor[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrT = centerPeTensor[iLev][mfi].array();
      const Array4<Real const>& arrB = bTot[mfi].const_array();
      const Array4<Real>& arrPeDiag = centerPe[iLev][mfi].array();
      const Array4<Real const> arrState =
          useAnisoPe ? centerPeState[iLev][mfi].const_array()
                     : Array4<Real const>();

      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const Real bx = arrB(i, j, k, ix_);
        const Real by = arrB(i, j, k, iy_);
        const Real bz = arrB(i, j, k, iz_);
        const Real b2 = bx * bx + by * by + bz * bz;
        const Real bMag = std::sqrt(b2);

        Real pe, pePar, pePerp;
        if (isLe) {
          const Real rho = arrPeDiag(i, j, k);
          le2009_closure(rho, bMag, p0, invRho0, invB0, bMin, peMinLocal, useFh,
                         fhCoef, pePar, pePerp);
          pe = (pePar + 2.0 * pePerp) / 3.0;
          arrPeDiag(i, j, k) = pe;
        } else {
          pe = amrex::max(arrState(i, j, k, 0), peMinLocal);
          pePar = amrex::max(arrState(i, j, k, 1), peMinLocal);
          pePerp = amrex::max(0.5 * (3.0 * pe - pePar), peMinLocal);
        }

        Real bhx = 0.0, bhy = 0.0, bhz = 0.0, wB = 0.0;
        if (b2 > 0.0) {
          const Real invB = 1.0 / bMag;
          bhx = bx * invB;
          bhy = by * invB;
          bhz = bz * invB;
          wB = b2 / (b2 + bMin2);
        }

        const Real dP = wB * (pePar - pePerp);
        const Real pPerpEff = pe - dP / 3.0;
        arrT(i, j, k, 0) = pPerpEff + dP * bhx * bhx;
        arrT(i, j, k, 1) = pPerpEff + dP * bhy * bhy;
        arrT(i, j, k, 2) = pPerpEff + dP * bhz * bhz;
        arrT(i, j, k, 3) = dP * bhx * bhy;
        arrT(i, j, k, 4) = dP * bhx * bhz;
        arrT(i, j, k, 5) = dP * bhy * bhz;
      });
    }

    centerPeTensor[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev > 0) {
      fill_fine_lev_bny_from_coarse(
          centerPeTensor[iLev - 1], centerPeTensor[iLev], 0, 6,
          refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          *get_cell_interp());
    }
    if (!Geom(iLev).isAllPeriodic() && centerPeTensor[iLev].nGrow() > 0) {
      apply_pe_zero_gradient_bc(iLev, centerPeTensor[iLev]);
    }
    apply_fake2d_k_clamp(centerPeTensor[iLev]);

    div_sym_tensor_center_to_node(centerPeTensor[iLev], nodeEambi[iLev],
                                  Geom(iLev).InvCellSize());
  } else {
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
        parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
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

    // Zero-gradient (Neumann / foextrap) BC across non-periodic domain
    // boundaries
    if (!Geom(iLev).isAllPeriodic() && centerPe[iLev].nGrow() > 0) {
      apply_pe_zero_gradient_bc(iLev, centerPe[iLev]);
    }

    apply_fake2d_k_clamp(centerPe[iLev]);

    // Compute grad_center_to_node(Pe) directly into nodeEambi
    grad_center_to_node(centerPe[iLev], nodeEambi[iLev],
                        Geom(iLev).InvCellSize());
  }

  // Scale in-place: E_ambi = -grad(Pe) / max(rho, rhoMinOhm)
  for (MFIter mfi(nodeEambi[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrEambi = nodeEambi[iLev][mfi].array();
    const Array4<Real const>& moments = nodePlasma[nSpecies][iLev][mfi].array();

    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
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
      parallel_for_omp(bxValid, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
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
      parallel_for_omp(fbox,
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
                 iLev, false, 0.0);
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
        parallel_for_omp(bxValid, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
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
  for (int iLev = 1; iLev < n_lev(); ++iLev) {
    fill_fine_lev_bny_from_coarse(
        nodePlasmaPrev[nSpecies][iLev - 1], nodePlasmaPrev[nSpecies][iLev], 0,
        nodePlasmaPrev[nSpecies][iLev].nComp(), refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
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
// Zero-gradient (foextrap) ghost cells for electron-pressure scratch arrays.
// Pure Neumann on every non-periodic face, then FillBoundary.
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

  mf.FillBoundary(Geom(iLev).periodicity());
}

//==========================================================
// Fake-2D (one cell in z) ghost clamp: every z layer outside the valid range
// takes the value of the nearest valid layer across all components.
void Pic::apply_fake2d_k_clamp(MultiFab& mf) {
  if (!is_fake_2d())
    return;

  const int nComp = mf.nComp();
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const auto& vbox = mfi.validbox();
    const auto& fbox = mfi.fabbox();
    auto arr = mf[mfi].array();
    const int klo = vbox.smallEnd(2);
    const int khi = vbox.bigEnd(2);
    parallel_for_omp(fbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int k_src = std::clamp(k, klo, khi);
      if (k != k_src) {
        for (int n = 0; n < nComp; ++n) {
          arr(i, j, k, n) = arr(i, j, k_src, n);
        }
      }
    });
  }
}

//==========================================================
void Pic::apply_centerPe_BC(int iLev) {
  MultiFab& mf = centerPeState[iLev];
  mf.FillBoundary(Geom(iLev).periodicity());

  if (iLev > 0) {
    fill_fine_lev_bny_from_coarse(
        centerPeState[iLev - 1], mf, 0, nPeComps_, refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), cell_status(iLev), *get_cell_interp());
  }

  apply_pe_zero_gradient_bc(iLev, mf);

  if (fi && fi->is_init_from_swmf() && fi->get_useMhdPe() &&
      fi->get_iPe() >= 0 && bcField.has(FieldBC::coupled)) {
    fill_ext_dir(cell_status(iLev), mf, 0, 1, &Pic::get_center_Pe, iLev,
                 &bcField);
    if (nPeComps_ > 1) {
      // Coupled MHD boundary supplies scalar Pe; seed PePar = Pe on those
      // Dirichlet ghosts.
      for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& fbox = mfi.fabbox();
        const Box& vbox = mfi.validbox();
        const Array4<Real>& arr = mf[mfi].array();
        const Array4<int const>& stat = cell_status(iLev)[mfi].const_array();
        parallel_for_omp(fbox,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                           if (!vbox.contains(IntVect(AMREX_D_DECL(i, j, k))) &&
                               bit::is_domain_boundary(stat(i, j, k))) {
                             arr(i, j, k, 1) = arr(i, j, k, 0);
                           }
                         });
      }
    }
  }

  apply_fake2d_k_clamp(mf);
  mf.FillBoundary(Geom(iLev).periodicity());
}

//==========================================================
void Pic::init_electron_pressure() {
  for (int iLev = 0; iLev < n_lev(); ++iLev)
    init_electron_pressure(iLev);
  peStateInitialized_ = true;
  update_max_evolved_Te();
}

//==========================================================
void Pic::compute_node_electron_pressure(int iLev, const MultiFab& nodeNe,
                                         MultiFab& nodePe) {
  if (useElectronPressureEq && peStateInitialized_) {
    // Evolved Pe: average the cell-centered state to the nodes. The ghost
    // cells are filled as in compute_ambipolar_E: same level, coarse-fine
    // boundary from the coarser level, zero gradient at the domain boundary.
    MultiFab centerTmp(centerPeState[iLev].boxArray(),
                       centerPeState[iLev].DistributionMap(), 1, 1);
    MultiFab::Copy(centerTmp, centerPeState[iLev], 0, 0, 1, 1);
    centerTmp.FillBoundary(Geom(iLev).periodicity());
    if (iLev > 0)
      fill_fine_lev_bny_from_coarse(
          centerPeState[iLev - 1], centerTmp, 0, 1, refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), cell_status(iLev), *get_cell_interp());
    if (!Geom(iLev).isAllPeriodic())
      apply_pe_zero_gradient_bc(iLev, centerTmp);
    average_center_to_node(centerTmp, nodePe);
    return;
  }

  if (electronTemperature <= 0) {
    nodePe.setVal(0.0);
    return;
  }

  const Real p0 = electronDensity0 * electronTemperature;
  const Real invRho0 = (electronDensity0 > 0) ? (1.0 / electronDensity0) : 0.0;
  const Real gamma = electronGamma;
  const Real Te = electronTemperature;
  for (MFIter mfi(nodePe); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrPe = nodePe[mfi].array();
    const Array4<Real const>& arrNe = nodeNe[mfi].const_array();
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      arrPe(i, j, k) =
          polytropic_closure(arrNe(i, j, k), p0, invRho0, gamma, Te);
    });
  }
}

//==========================================================
void Pic::init_electron_pressure(int iLev) {
  // Start from the coupled MHD Pe (when available) or the algebraic polytropic
  // closure evaluated on the current ion density.
  compute_electron_density(iLev, nodePeRho[iLev], centerPeRho[iLev], false);

  const Real p0 = electronDensity0 * electronTemperature;
  const Real invRho0 = (electronDensity0 > 0.0) ? 1.0 / electronDensity0 : 0.0;
  const Real gamma = electronGamma;
  const Real Te = electronTemperature;
  const Real peMinLocal = peMin;

  const bool fromShockTube = (fi && fi->get_is_shock_tube());
  const bool fromCoupledMhd = !fromShockTube && fi && fi->is_init_from_swmf() &&
                              fi->get_useMhdPe() && fi->get_iPe() >= 0 &&
                              fi->is_node_fluid_ready();
  if (fromCoupledMhd) {
    const MultiFab nodeMhdPe(fi->get_node_fluid(iLev), amrex::make_alias,
                             fi->get_iPe(), 1);
    average_node_to_center(nodeMhdPe, centerPe[iLev]);
  }

  const Real peL = fromShockTube ? fi->get_shock_pe_l() : 0.0;
  const Real peR = fromShockTube ? fi->get_shock_pe_r() : 0.0;
  const Real xShock = fromShockTube ? fi->get_shock_position() : 0.0;
  const Real slope = fromShockTube ? fi->get_shock_slope() : 0.0;
  const auto& geom = Geom(iLev);
  const auto& prob_lo = geom.ProbLo();
  const auto& dx = geom.CellSize();

  const bool isAniso = useAnisoPe;
  const Real anisoRatio0 = peAnisoRatio0;

  for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrPe = centerPeState[iLev][mfi].array();
    const Array4<Real const>& arrRho = centerPeRho[iLev][mfi].array();
    const Array4<Real const> arrCoupled =
        fromCoupledMhd ? centerPe[iLev][mfi].const_array()
                       : Array4<Real const>();
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      Real peInit;
      if (fromShockTube) {
        Real xPlane = prob_lo[0] + (i + 0.5) * dx[0];
#if AMREX_SPACEDIM > 1
        if (slope != 0.0) {
          xPlane += slope * (prob_lo[1] + (j + 0.5) * dx[1]);
        }
#endif
        peInit = (xPlane < xShock) ? peL : peR;
      } else if (fromCoupledMhd) {
        peInit = arrCoupled(i, j, k);
      } else {
        peInit = polytropic_closure(arrRho(i, j, k), p0, invRho0, gamma, Te);
      }
      const Real peVal = amrex::max(peInit, peMinLocal);
      arrPe(i, j, k, 0) = peVal;
      if (isAniso) {
        arrPe(i, j, k, 1) = amrex::max(
            3.0 * peVal / (1.0 + 2.0 * anisoRatio0), peMinLocal);
      }
    });
  }
  apply_centerPe_BC(iLev);
}

//==========================================================
void Pic::init_le2009_b0() {
  if (!useLe2009Pe || le2009B0 > 0.0 || n_lev() <= 0 || centerB[0].empty())
    return;

  centerB[0].FillBoundary(Geom(0).periodicity());
  const MultiFab& bTot = total_center_B(centerB[0], 0);
  auto const& bArr = bTot.const_arrays();
  auto const& statArr = cell_status(0).const_arrays();
  const bool hasBody = use_body();

  auto [bSum, count] = ParReduce(
      TypeList<ReduceOpSum, ReduceOpSum>{}, TypeList<Real, long>{}, bTot,
      IntVect(0),
      [=] AMREX_GPU_DEVICE(int nb, int i, int j,
                           int k) noexcept -> GpuTuple<Real, long> {
        if (hasBody && bit::is_body(statArr[nb](i, j, k)))
          return { 0.0, 0L };
        const Real bx = bArr[nb](i, j, k, ix_);
        const Real by = bArr[nb](i, j, k, iy_);
        const Real bz = bArr[nb](i, j, k, iz_);
        return { std::sqrt(bx * bx + by * by + bz * bz), 1L };
      });
  ParallelDescriptor::ReduceRealSum(bSum);
  ParallelDescriptor::ReduceLongSum(count);

  if (count <= 0 || bSum <= 0.0) {
    Abort("Pic::init_le2009_b0: mean |B| on level 0 is zero; provide a "
          "positive le2009B0SI in #ELECTRONANISO.");
  }
  le2009B0 = bSum / static_cast<Real>(count);
  Print() << "  electron pressure anisotropy (le2009): auto reference B0 = "
          << le2009B0 << " [code units]\n";
}

//==========================================================
void Pic::eval_le2009_cell_pressures(int iLev, const MFIter& mfi, IntVect ijk,
                                     Real& peOut, Real& peParOut,
                                     Real& pePerpOut, Real* neOut) const {
  peOut = 0.0;
  peParOut = 0.0;
  pePerpOut = 0.0;
  if (neOut) {
    *neOut = 0.0;
  }
  if (!useLe2009Pe || !usePIC || iLev >= static_cast<int>(centerB.size()) ||
      centerB[iLev].empty()) {
    return;
  }

  // Reconstruct quasi-neutral electron density at cell center from surrounding
  // 8 (4 in 2D) corner nodes of nodePlasma[s].
  Real ne = 0.0;
  const int i = ijk[ix_];
  const int j = ijk[iy_];
  const int k = (nDim > 2) ? ijk[iz_] : 0;
  const int kMax = (nDim > 2 && !is_fake_2d()) ? 1 : 0;
  const Real invN = (kMax > 0) ? 0.125 : 0.25;
  for (int s = 0; s < nSpecies; ++s) {
    if (s >= static_cast<int>(nodePlasma.size()) ||
        iLev >= static_cast<int>(nodePlasma[s].size()) ||
        nodePlasma[s][iLev].empty()) {
      continue;
    }
    const Real qomFactor = (fi && s < fi->get_nS() && fi->get_MoMi(s) > 0.0)
                               ? (fi->get_QoQi(s) / fi->get_MoMi(s))
                               : 1.0;
    const auto& mom = nodePlasma[s][iLev][mfi].const_array();
    Real rhoSum = 0.0;
    for (int dk = 0; dk <= kMax; ++dk)
      for (int dj = 0; dj <= 1; ++dj)
        for (int di = 0; di <= 1; ++di)
          rhoSum += mom(i + di, j + dj, k + dk, iRho_);
    ne += qomFactor * rhoSum * invN;
  }
  if (neOut) {
    *neOut = ne;
  }

  const auto& b1 = centerB[iLev][mfi].const_array();
  Real bx = b1(i, j, k, ix_);
  Real by = b1(i, j, k, iy_);
  Real bz = b1(i, j, k, iz_);
  if (use_intrinsic_B() && iLev < static_cast<int>(centerB0.size()) &&
      !centerB0[iLev].empty()) {
    const auto& b0 = centerB0[iLev][mfi].const_array();
    bx += b0(i, j, k, ix_);
    by += b0(i, j, k, iy_);
    bz += b0(i, j, k, iz_);
  }
  const Real bMag = std::sqrt(bx * bx + by * by + bz * bz);
  const Real p0 = electronDensity0 * electronTemperature;
  const Real invRho0 = (electronDensity0 > 0.0) ? (1.0 / electronDensity0) : 0.0;
  const Real invB0 = (le2009B0 > 0.0) ? (1.0 / le2009B0) : 0.0;
  le2009_closure(ne, bMag, p0, invRho0, invB0, fieldAlignedBMin, peMin,
                 usePeInstabilityLimit, peFirehoseCoef, peParOut, pePerpOut);
  peOut = (peParOut + 2.0 * pePerpOut) / 3.0;
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
                                  0, nPeComps_, refRatio(iLev - 1),
                                  Geom(iLev - 1), Geom(iLev), cell_status(iLev),
                                  cellInterp);
  }
  if (peStateInitialized_) {
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      compute_electron_density(iLev, nodePeRho[iLev], centerPeRho[iLev], false);
    }
    update_max_evolved_Te();
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
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
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
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real ne = amrex::max(arrRho(i, j, k), rhoFloor);
      arrTe(i, j, k) = amrex::max(arrPe(i, j, k, 0), peMinLocal) / ne;
    });
  }
  apply_pe_zero_gradient_bc(iLev, centerPeTe[iLev]);
  apply_fake2d_k_clamp(centerPeTe[iLev]);
}

//==========================================================
Real Pic::max_evolved_Te(int iLev) const {
  if (useElectronPressureEq && peStateInitialized_ &&
      iLev < static_cast<int>(maxEvolvedTe_.size()) &&
      maxEvolvedTe_[iLev] > 0.0) {
    return maxEvolvedTe_[iLev];
  }
  return electronTemperature;
}

//==========================================================
void Pic::update_max_evolved_Te() {
  if (!useElectronPressureEq || !peStateInitialized_)
    return;
  if (static_cast<int>(maxEvolvedTe_.size()) < n_lev()) {
    maxEvolvedTe_.resize(n_lev(), 0.0);
  }
  const Real peMinLocal = peMin;
  const Real rhoFloor = amrex::max(rhoMinOhm, static_cast<Real>(1e-10));
  const bool hasBody = use_body();

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    if (centerPeState[iLev].empty() || centerPeRho[iLev].empty()) {
      maxEvolvedTe_[iLev] = electronTemperature;
      continue;
    }
    const bool hasFiner = (iLev < finestLevel());
    auto const& peArrays = centerPeState[iLev].const_arrays();
    auto const& rhoArrays = centerPeRho[iLev].const_arrays();
    auto const& cStatArrays = cell_status(iLev).const_arrays();
    Real teMax = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
                           centerPeState[iLev], IntVect(0),
                           [=] AMREX_GPU_DEVICE(int nb, int i, int j, int k)
                               noexcept -> GpuTuple<Real> {
                                 const int st = cStatArrays[nb](i, j, k);
                                 if ((hasBody && bit::is_body(st)) ||
                                     (hasFiner && bit::is_refined(st)))
                                   return 0.0;
                                 const Real ne = rhoArrays[nb](i, j, k);
                                 if (ne <= 0.0)
                                   return 0.0;
                                 const Real pe = amrex::max(
                                     peArrays[nb](i, j, k, 0), peMinLocal);
                                 return pe / amrex::max(ne, rhoFloor);
                               });
    ParallelDescriptor::ReduceRealMax(teMax);
    maxEvolvedTe_[iLev] = (teMax > 0.0) ? teMax : electronTemperature;
  }
}

//==========================================================
// Resistive Ohmic heating: dPe/dt = (gamma_e - 1) * eta_tot * |J|^2.
// Assumes nodeJ[iLev] already contains curl(centerB[iLev]) from
// electron_velocity_at_nodes(iLev).
void Pic::add_ohmic_heating(int iLev, Real dt) {
  if (!useOhmicHeating || (etaResistivity <= 0.0 && !hasRegionalResistivity_)) {
    return;
  }

  const Real invFourPI = 1.0 / fourPI;
  const Real etaGlobal = etaResistivity;
  const bool hasReg = hasRegionalResistivity_;

  for (MFIter mfi(nodePeRho[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrQ = nodePeRho[iLev][mfi].array();
    const Array4<Real const>& arrJ = nodeJ[iLev][mfi].const_array();
    const Array4<Real const> arrEtaReg =
        hasReg ? nodeEtaRegional[iLev][mfi].const_array()
               : Array4<Real const>();
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real etaTot =
          etaGlobal + (hasReg ? arrEtaReg(i, j, k) : static_cast<Real>(0.0));
      if (etaTot > 0.0) {
        const Real jx = arrJ(i, j, k, ix_) * invFourPI;
        const Real jy = arrJ(i, j, k, iy_) * invFourPI;
        const Real jz = arrJ(i, j, k, iz_) * invFourPI;
        arrQ(i, j, k) = etaTot * (jx * jx + jy * jy + jz * jz);
      } else {
        arrQ(i, j, k) = 0.0;
      }
    });
  }
  nodePeRho[iLev].FillBoundary(Geom(iLev).periodicity());
  average_node_to_center(nodePeRho[iLev], centerPe[iLev]);

  const Real gammaM1 = electronGamma - 1.0;
  const Real peMinLocal = peMin;
  const int nComps = nPeComps_;
  for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrPe = centerPeState[iLev][mfi].array();
    const Array4<Real const>& arrQ = centerPe[iLev][mfi].const_array();
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real dPe = dt * gammaM1 * arrQ(i, j, k);
      if (dPe > 0.0) {
        arrPe(i, j, k, 0) = amrex::max(arrPe(i, j, k, 0) + dPe, peMinLocal);
        if (nComps > 1) {
          arrPe(i, j, k, 1) = amrex::max(arrPe(i, j, k, 1) + dPe, peMinLocal);
        }
      }
    });
  }
}

//==========================================================
// Electron-ion collisional thermal equilibration (heat exchange) hook:
// dPe/dt = (Pi - Pe) / tau_eq with kinetic ion Pi held fixed during the
// electron update: exact point-implicit denominator is 1 + dt * H.
void Pic::add_electron_ion_heating(int iLev, Real dt) {
  if (!useHeatExchange || collisionCoefEi <= 0.0)
    return;

  // 1) Compute scalar ion pressure Pi = Tr(P_i)/3 at nodes in nodePeRho
  for (MFIter mfi(nodePeRho[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrPiN = nodePeRho[iLev][mfi].array();
    const Array4<Real const>& moments = nodePlasma[nSpecies][iLev][mfi].array();
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      arrPiN(i, j, k) = (moments(i, j, k, iPxx_) + moments(i, j, k, iPyy_) +
                         moments(i, j, k, iPzz_)) /
                        3.0;
    });
  }
  nodePeRho[iLev].FillBoundary(Geom(iLev).periodicity());

  // 2) Average scalar ion pressure to cell centers in centerPe (scratch here)
  average_node_to_center(nodePeRho[iLev], centerPe[iLev]);

  // 3) Point-implicit energy exchange:
  //    H = collisionCoefEi * ne / Te^1.5
  //    PePImpl = H / (1 + dt * H)
  //    DeltaPe = dt * PePImpl * (Pi - Pe)
  const Real cEi = collisionCoefEi;
  const Real peMinLocal = peMin;
  const Real rhoFloor = rhoMinOhm;
  const int nComps = nPeComps_;

  for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrPe = centerPeState[iLev][mfi].array();
    const Array4<Real const>& arrPi = centerPe[iLev][mfi].const_array();
    const Array4<Real const>& arrRho = centerPeRho[iLev][mfi].const_array();

    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real ne = amrex::max(arrRho(i, j, k), rhoFloor);
      const Real pe = amrex::max(arrPe(i, j, k, 0), peMinLocal);
      const Real pi = amrex::max(arrPi(i, j, k), 0.0);
      const Real te = pe / ne;

      if (te > 0.0) {
        const Real H = cEi * ne / std::pow(te, 1.5);
        const Real pePImpl = H / (1.0 + dt * H);
        const Real deltaPe = dt * pePImpl * (pi - pe);
        arrPe(i, j, k, 0) = amrex::max(pe + deltaPe, peMinLocal);
        if (nComps > 1) {
          const Real pePar = amrex::max(arrPe(i, j, k, 1), peMinLocal);
          const Real deltaPePar = dt * pePImpl * (pi - pePar);
          arrPe(i, j, k, 1) = amrex::max(pePar + deltaPePar, peMinLocal);
        }
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
      average_down(centerPeState[iLev], centerPeState[iLev - 1], 0, nPeComps_,
                   refRatio(iLev - 1));
  }

  update_max_evolved_Te();
}

//==========================================================
void Pic::update_Pe_hybrid(int iLev, Real dt) {
  BL_PROFILE("Pic::update_Pe_hybrid");

  //--------------------------------------------------------------------
  // Operator-split advance of the electron pressure on one level:
  //   dPe/dt + div(u_e Pe) + (gamma_e-1) Pe div(u_e)
  //       = (gamma_e-1) [ eta |J|^2 + div(kappa_hat . grad(Te)) + H_ei ]
  //--------------------------------------------------------------------

  electron_velocity_at_nodes(iLev);
  advect_electron_pressure(iLev, dt);
  add_ohmic_heating(iLev, dt);

  compute_electron_density(iLev, nodePeRho[iLev], centerPeRho[iLev], false);
  compute_electron_temperature(iLev);

  // A no-op each when heatCondKappa0 is 0 / #ELECTRONCOLLISION is absent.
  apply_electron_heat_conduction(iLev, dt);
  add_electron_ion_heating(iLev, dt);
  if (useAnisoPe) {
    apply_electron_isotropization(iLev, dt);
  }

  apply_centerPe_BC(iLev);
}

//==========================================================
// Electron velocity u_e = (J_i - J_tot) / (e * n_e) on the NODES, using
// charge-weighted sums over all kinetic ion species so multi-species runs
// stay consistent with compute_electron_density().
void Pic::electron_velocity_at_nodes(int iLev) {
  BL_PROFILE("Pic::electron_velocity_at_nodes");

  curl_center_to_node(centerB[iLev], nodeJ[iLev], Geom(iLev).InvCellSize());
  nodeJ[iLev].FillBoundary(Geom(iLev).periodicity());

  nodePeRho[iLev].setVal(0.0);
  nodePeVec[iLev].setVal(0.0);
  for (int s : kineticSpecies_) {
    const Real qomFactor = (fi && s < fi->get_nS() && fi->get_MoMi(s) > 0.0)
                               ? (fi->get_QoQi(s) / fi->get_MoMi(s))
                               : 1.0;
    for (MFIter mfi(nodePeRho[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrNe = nodePeRho[iLev][mfi].array();
      const Array4<Real>& arrJi = nodePeVec[iLev][mfi].array();
      const Array4<Real const>& arrMom = nodePlasma[s][iLev][mfi].const_array();
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        arrNe(i, j, k) += qomFactor * arrMom(i, j, k, iRho_);
        for (int d = 0; d < nDim; ++d) {
          arrJi(i, j, k, d) += qomFactor * arrMom(i, j, k, iMx_ + d);
        }
      });
    }
  }

  const Real invFourPILocal = 1.0 / fourPI;
  const Real rhoFloor = rhoMinOhm;

  for (MFIter mfi(nodePeVec[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrUe = nodePeVec[iLev][mfi].array();
    const Array4<Real const>& arrNe = nodePeRho[iLev][mfi].const_array();
    const Array4<Real const>& arrJ = nodeJ[iLev][mfi].const_array();
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real rho = arrNe(i, j, k);
      if (rho > 0.0) {
        const Real invRho = 1.0 / amrex::max(rho, rhoFloor);
        for (int d = 0; d < nDim; ++d) {
          arrUe(i, j, k, d) =
              (arrUe(i, j, k, d) - arrJ(i, j, k, d) * invFourPILocal) * invRho;
        }
      } else {
        for (int d = 0; d < nDim; ++d) {
          arrUe(i, j, k, d) = 0.0;
        }
      }
    });
  }
  apply_pe_zero_gradient_bc(iLev, nodePeVec[iLev]);
  apply_fake2d_k_clamp(nodePeVec[iLev]);
}

//==========================================================
// Advection of Pe: dPe/dt = -div(u_e Pe) - (gamma_e-1) Pe div(u_e).
//
// Face-averaged u_e with TVD/MUSCL reconstruction, CFL-subcycled 2-stage
// SSP-RK2 (Heun) on the flux divergence, and exponential or explicit
// compression.
void Pic::advect_electron_pressure(int iLev, Real dt) {
  BL_PROFILE("Pic::advect_electron_pressure");

  const Real* invDxGeom = Geom(iLev).InvCellSize();
  const IntVect domLen = Geom(iLev).Domain().length();
  const Real invDxX = (domLen[ix_] > 1) ? invDxGeom[ix_] : 0.0;
  const Real invDxY = (domLen[iy_] > 1) ? invDxGeom[iy_] : 0.0;
  const Real invDxZ = (nDim > 2 && domLen[iz_] > 1) ? invDxGeom[iz_] : 0.0;

  const Real gammaM1 = electronGamma - 1.0;
  const int limType = peLimiterType;
  const bool compExp = peCompressionExp;
  const Real peMinLocal = peMin;
  const bool isFake2D = is_fake_2d();

  // Electron advection CFL subcycling bound.
  auto const& ueArrays = nodePeVec[iLev].const_arrays();
  Real maxAdvRate = ParReduce(
      TypeList<ReduceOpMax>{}, TypeList<Real>{}, centerPeState[iLev],
      IntVect(0),
      [=] AMREX_GPU_DEVICE(int nb, int i, int j,
                           int k) noexcept -> GpuTuple<Real> {
        Real rate = 0.0;
        if (invDxX > 0.0) {
          const Real u0 = ue_face<ix_>(ueArrays[nb], i, j, k, 0, isFake2D);
          const Real u1 = ue_face<ix_>(ueArrays[nb], i, j, k, 1, isFake2D);
          rate += amrex::max(std::abs(u0), std::abs(u1)) * invDxX;
        }
        if (invDxY > 0.0) {
          const Real u0 = ue_face<iy_>(ueArrays[nb], i, j, k, 0, isFake2D);
          const Real u1 = ue_face<iy_>(ueArrays[nb], i, j, k, 1, isFake2D);
          rate += amrex::max(std::abs(u0), std::abs(u1)) * invDxY;
        }
        if constexpr (nDim > 2) {
          if (invDxZ > 0.0) {
            const Real u0 = ue_face<iz_>(ueArrays[nb], i, j, k, 0, isFake2D);
            const Real u1 = ue_face<iz_>(ueArrays[nb], i, j, k, 1, isFake2D);
            rate += amrex::max(std::abs(u0), std::abs(u1)) * invDxZ;
          }
        }
        return rate;
      });
  ParallelDescriptor::ReduceRealMax(maxAdvRate);

  const int nAdvSub =
      (maxAdvRate > 0.0)
          ? std::clamp(static_cast<int>(std::ceil(dt * maxAdvRate / 0.45)), 1,
                       20)
          : 1;
  const Real dtSub = dt / static_cast<Real>(nAdvSub);

  if (useAnisoPe) {
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
    const MultiFab& bTot = total_center_B(centerB[iLev], iLev);
    const BoxArray& ba = centerPeState[iLev].boxArray();
    const DistributionMapping& dm = centerPeState[iLev].DistributionMap();
    const Real bMin2 = fieldAlignedBMin * fieldAlignedBMin;

    // 2-component work arrays: comp 0 = PePar, comp 1 = PePerp.
    MultiFab peAnisoWork(ba, dm, 2, get_n_ghost());
    MultiFab peAniso0(ba, dm, 2, 0);
    MultiFab peAniso1(ba, dm, 2, 0);

    auto fill_aniso_ghosts = [&](MultiFab& mf) {
      mf.FillBoundary(Geom(iLev).periodicity());
      if (iLev > 0) {
        MultiFab coarseParPerp(centerPeState[iLev - 1].boxArray(),
                               centerPeState[iLev - 1].DistributionMap(), 2,
                               centerPeState[iLev - 1].nGrow());
        for (MFIter mfi(coarseParPerp); mfi.isValid(); ++mfi) {
          const Box& fbox = mfi.fabbox();
          const Array4<Real>& arrC = coarseParPerp[mfi].array();
          const Array4<Real const>& arrS =
              centerPeState[iLev - 1][mfi].const_array();
          parallel_for_omp(fbox,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                             const Real pe = amrex::max(arrS(i, j, k, 0), peMinLocal);
                             const Real pPar =
                                 amrex::max(arrS(i, j, k, 1), peMinLocal);
                             arrC(i, j, k, 0) = pPar;
                             arrC(i, j, k, 1) = amrex::max(
                                 0.5 * (3.0 * pe - pPar), peMinLocal);
                           });
        }
        fill_fine_lev_bny_from_coarse(coarseParPerp, mf, 0, 2,
                                      refRatio(iLev - 1), Geom(iLev - 1),
                                      Geom(iLev), cell_status(iLev),
                                      *get_cell_interp());
      }
      apply_pe_zero_gradient_bc(iLev, mf);
      apply_fake2d_k_clamp(mf);
    };

    apply_centerPe_BC(iLev);
    for (MFIter mfi(peAnisoWork); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrW = peAnisoWork[mfi].array();
      const Array4<Real const>& arrS = centerPeState[iLev][mfi].const_array();
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real pe = amrex::max(arrS(i, j, k, 0), peMinLocal);
        const Real pPar = amrex::max(arrS(i, j, k, 1), peMinLocal);
        arrW(i, j, k, 0) = pPar;
        arrW(i, j, k, 1) = amrex::max(0.5 * (3.0 * pe - pPar), peMinLocal);
      });
    }

    for (int iSub = 0; iSub < nAdvSub; ++iSub) {
      fill_aniso_ghosts(peAnisoWork);
      MultiFab::Copy(peAniso0, peAnisoWork, 0, 0, 2, 0);

      // Stage 1: advect both PePar and PePerp.
      for (MFIter mfi(peAnisoWork); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real const>& arrW = peAnisoWork[mfi].const_array();
        const Array4<Real>& arr1 = peAniso1[mfi].array();
        const Array4<Real const>& arrUe = nodePeVec[iLev][mfi].const_array();
        parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          for (int c = 0; c < 2; ++c) {
            Real divFlux = 0.0;
            Real divUe = 0.0;
            if (invDxX > 0.0)
              pe_dir_flux<ix_>(arrW, c, arrUe, i, j, k, invDxX, limType,
                               peMinLocal, isFake2D, divFlux, divUe);
            if (invDxY > 0.0)
              pe_dir_flux<iy_>(arrW, c, arrUe, i, j, k, invDxY, limType,
                               peMinLocal, isFake2D, divFlux, divUe);
            if constexpr (nDim > 2) {
              if (invDxZ > 0.0)
                pe_dir_flux<iz_>(arrW, c, arrUe, i, j, k, invDxZ, limType,
                                 peMinLocal, isFake2D, divFlux, divUe);
            }
            arr1(i, j, k, c) =
                amrex::max(arrW(i, j, k, c) - dtSub * divFlux, peMinLocal);
          }
        });
      }
      MultiFab::Copy(peAnisoWork, peAniso1, 0, 0, 2, 0);
      fill_aniso_ghosts(peAnisoWork);

      // Stage 2: SSP-RK2 advection + CGL exponential compression.
      for (MFIter mfi(peAnisoWork); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real const>& arr0 = peAniso0[mfi].const_array();
        const Array4<Real const>& arr1 = peAnisoWork[mfi].const_array();
        const Array4<Real>& arrOut = peAniso1[mfi].array();
        const Array4<Real const>& arrUe = nodePeVec[iLev][mfi].const_array();
        const Array4<Real const>& arrB = bTot[mfi].const_array();

        parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          Real pAdv[2] = { 0.0, 0.0 };
          Real divUeUnused = 0.0;
          for (int c = 0; c < 2; ++c) {
            Real divFlux1 = 0.0;
            if (invDxX > 0.0)
              pe_dir_flux<ix_>(arr1, c, arrUe, i, j, k, invDxX, limType,
                               peMinLocal, isFake2D, divFlux1, divUeUnused);
            if (invDxY > 0.0)
              pe_dir_flux<iy_>(arr1, c, arrUe, i, j, k, invDxY, limType,
                               peMinLocal, isFake2D, divFlux1, divUeUnused);
            if constexpr (nDim > 2) {
              if (invDxZ > 0.0)
                pe_dir_flux<iz_>(arr1, c, arrUe, i, j, k, invDxZ, limType,
                                 peMinLocal, isFake2D, divFlux1, divUeUnused);
            }
            pAdv[c] = amrex::max(
                0.5 * arr0(i, j, k, c) +
                    0.5 * (arr1(i, j, k, c) - dtSub * divFlux1),
                peMinLocal);
          }

          // Cell-centered Jacobian du_{b}/dx_{a} matching div_node_to_center.
          Real du[3][3] = { { 0.0, 0.0, 0.0 },
                            { 0.0, 0.0, 0.0 },
                            { 0.0, 0.0, 0.0 } };
          for (int b = 0; b < 3; ++b) {
            if (invDxX > 0.0) {
              du[0][b] =
                  (ue_face_comp<ix_>(arrUe, b, i, j, k, 1, isFake2D) -
                   ue_face_comp<ix_>(arrUe, b, i, j, k, 0, isFake2D)) *
                  invDxX;
            }
            if (invDxY > 0.0) {
              du[1][b] =
                  (ue_face_comp<iy_>(arrUe, b, i, j, k, 1, isFake2D) -
                   ue_face_comp<iy_>(arrUe, b, i, j, k, 0, isFake2D)) *
                  invDxY;
            }
            if constexpr (nDim > 2) {
              if (invDxZ > 0.0) {
                du[2][b] =
                    (ue_face_comp<iz_>(arrUe, b, i, j, k, 1, isFake2D) -
                     ue_face_comp<iz_>(arrUe, b, i, j, k, 0, isFake2D)) *
                    invDxZ;
              }
            }
          }
          const Real divUe = du[0][0] + du[1][1] + du[2][2];

          const Real bx = arrB(i, j, k, ix_);
          const Real by = arrB(i, j, k, iy_);
          const Real bz = arrB(i, j, k, iz_);
          const Real b2 = bx * bx + by * by + bz * bz;
          Real bbGradEff = divUe / 3.0;
          if (b2 > 0.0) {
            const Real invB = 1.0 / std::sqrt(b2);
            const Real bh[3] = { bx * invB, by * invB, bz * invB };
            Real bbGrad = 0.0;
            for (int a = 0; a < 3; ++a)
              for (int b = 0; b < 3; ++b)
                bbGrad += bh[a] * bh[b] * du[a][b];
            const Real wB = b2 / (b2 + bMin2);
            bbGradEff = wB * bbGrad + (1.0 - wB) * (divUe / 3.0);
          }

          if (compExp) {
            arrOut(i, j, k, 0) = amrex::max(
                pAdv[0] * std::exp(-2.0 * bbGradEff * dtSub), peMinLocal);
            arrOut(i, j, k, 1) = amrex::max(
                pAdv[1] * std::exp(-(divUe - bbGradEff) * dtSub), peMinLocal);
          } else {
            arrOut(i, j, k, 0) = amrex::max(
                pAdv[0] - dtSub * 2.0 * pAdv[0] * bbGradEff, peMinLocal);
            arrOut(i, j, k, 1) = amrex::max(
                pAdv[1] - dtSub * pAdv[1] * (divUe - bbGradEff), peMinLocal);
          }
        });
      }
      MultiFab::Copy(peAnisoWork, peAniso1, 0, 0, 2, 0);
    }

    for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrS = centerPeState[iLev][mfi].array();
      const Array4<Real const>& arrW = peAnisoWork[mfi].const_array();
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real pPar = amrex::max(arrW(i, j, k, 0), peMinLocal);
        const Real pPerp = amrex::max(arrW(i, j, k, 1), peMinLocal);
        arrS(i, j, k, 0) = (pPar + 2.0 * pPerp) / 3.0;
        arrS(i, j, k, 1) = pPar;
      });
    }
    return;
  }

  MultiFab centerPe0(centerPeState[iLev].boxArray(),
                     centerPeState[iLev].DistributionMap(), 1, 0);

  for (int iSub = 0; iSub < nAdvSub; ++iSub) {
    apply_centerPe_BC(iLev);
    MultiFab::Copy(centerPe0, centerPeState[iLev], 0, 0, 1, 0);

    // Stage 1: P1 = max(P0 - dtSub * div(u_e P0), peMin)
    for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real const>& arrPe = centerPeState[iLev][mfi].const_array();
      const Array4<Real>& arrPe1 = centerPe[iLev][mfi].array();
      const Array4<Real const>& arrUe = nodePeVec[iLev][mfi].const_array();

      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        Real divFlux = 0.0;
        Real divUe = 0.0;
        if (invDxX > 0.0) {
          pe_dir_flux<ix_>(arrPe, 0, arrUe, i, j, k, invDxX, limType,
                           peMinLocal, isFake2D, divFlux, divUe);
        }
        if (invDxY > 0.0) {
          pe_dir_flux<iy_>(arrPe, 0, arrUe, i, j, k, invDxY, limType,
                           peMinLocal, isFake2D, divFlux, divUe);
        }
        if constexpr (nDim > 2) {
          if (invDxZ > 0.0) {
            pe_dir_flux<iz_>(arrPe, 0, arrUe, i, j, k, invDxZ, limType,
                             peMinLocal, isFake2D, divFlux, divUe);
          }
        }
        arrPe1(i, j, k) =
            amrex::max(arrPe(i, j, k, 0) - dtSub * divFlux, peMinLocal);
      });
    }
    MultiFab::Copy(centerPeState[iLev], centerPe[iLev], 0, 0, 1, 0);
    apply_centerPe_BC(iLev);

    // Stage 2: P_adv = 0.5*P0 + 0.5*(P1 - dtSub*div(u_e P1)), then compress
    for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real const>& arrP0 = centerPe0[mfi].const_array();
      const Array4<Real const>& arrP1 = centerPeState[iLev][mfi].const_array();
      const Array4<Real>& arrPeNew = centerPe[iLev][mfi].array();
      const Array4<Real const>& arrUe = nodePeVec[iLev][mfi].const_array();

      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        Real divFlux1 = 0.0;
        Real divUe = 0.0;
        if (invDxX > 0.0) {
          pe_dir_flux<ix_>(arrP1, 0, arrUe, i, j, k, invDxX, limType,
                           peMinLocal, isFake2D, divFlux1, divUe);
        }
        if (invDxY > 0.0) {
          pe_dir_flux<iy_>(arrP1, 0, arrUe, i, j, k, invDxY, limType,
                           peMinLocal, isFake2D, divFlux1, divUe);
        }
        if constexpr (nDim > 2) {
          if (invDxZ > 0.0) {
            pe_dir_flux<iz_>(arrP1, 0, arrUe, i, j, k, invDxZ, limType,
                             peMinLocal, isFake2D, divFlux1, divUe);
          }
        }
        const Real pAdv = amrex::max(
            0.5 * arrP0(i, j, k) + 0.5 * (arrP1(i, j, k, 0) - dtSub * divFlux1),
            peMinLocal);
        if (compExp) {
          arrPeNew(i, j, k) =
              amrex::max(pAdv * std::exp(-gammaM1 * divUe * dtSub), peMinLocal);
        } else {
          arrPeNew(i, j, k) =
              amrex::max(pAdv - dtSub * gammaM1 * pAdv * divUe, peMinLocal);
        }
      });
    }
    MultiFab::Copy(centerPeState[iLev], centerPe[iLev], 0, 0, 1, 0);
  }
}

//==========================================================
// Conservative face-centered electron heat conduction solver:
//   dP/dt = -energyScale * div(q),   q = -kappaScale * kappa_hat . grad(T),
// with T = P / n_e, Spitzer kappa = heatCondKappa0 * T^2.5, compact 2-point
// normal gradients (checkerboard-free), MC-limited transverse gradients, and
// RKL2 / subcycle / symmetric point-implicit time integration.
void Pic::solve_heat_conduction_scalar(int iLev, Real dt, MultiFab& state,
                                       int iStateComp, const MultiFab& rho,
                                       const MultiFab& bTot, Real kappaScale,
                                       Real energyScale,
                                       const MultiFab* coarseState,
                                       int iCoarseComp) {
  if (heatCondKappa0 <= 0.0 || kappaScale <= 0.0 || energyScale <= 0.0 ||
      dt <= 0.0) {
    return;
  }

  const BoxArray& ba = state.boxArray();
  const DistributionMapping& dm = state.DistributionMap();
  const IntVect domLen = Geom(iLev).Domain().length();
  const auto domLo = Geom(iLev).Domain().smallEnd();
  const auto domHi = Geom(iLev).Domain().bigEnd();
  const Real* invDxGeom = Geom(iLev).InvCellSize();
  GpuArray<int, nDim3> activeDir = { 0, 0, 0 };
  GpuArray<Real, nDim3> invDxArr = { 0.0, 0.0, 0.0 };
  for (int d = 0; d < nDim; ++d) {
    if (domLen[d] > 1) {
      activeDir[d] = 1;
      invDxArr[d] = invDxGeom[d];
    }
  }
  if (!activeDir[ix_] && !activeDir[iy_] && (nDim <= 2 || !activeDir[iz_])) {
    return;
  }

  const Real peMinLocal = peMin;
  const Real rhoFloor = rhoMinOhm;
  const Real kappa0Eff = heatCondKappa0 * kappaScale;
  const Real fAlign = fieldAlignedConduction ? fieldAlignedFraction : 0.0;
  const Real bMin2 = fieldAlignedBMin * fieldAlignedBMin;
  const Real fluxLim = heatFluxLimiter;
  const Real massRatioPe = cProtonMassSI / cElectronMassSI;
  const bool coupledDirichlet = fi && fi->is_init_from_swmf() &&
                                fi->get_useMhdPe() && fi->get_iPe() >= 0 &&
                                bcField.has(FieldBC::coupled);

  // Ghost-filled density and total magnetic field (1 ghost layer is enough for
  // face averages).
  MultiFab rhoGhost(ba, dm, 1, 1);
  MultiFab::Copy(rhoGhost, rho, 0, 0, 1, 0);
  rhoGhost.FillBoundary(Geom(iLev).periodicity());
  if (iLev > 0 && &rho == &centerPeRho[iLev] &&
      !centerPeRho[iLev - 1].empty()) {
    fill_fine_lev_bny_from_coarse(
        centerPeRho[iLev - 1], rhoGhost, 0, 1, refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), cell_status(iLev), *get_cell_interp());
  }
  apply_pe_zero_gradient_bc(iLev, rhoGhost);
  apply_fake2d_k_clamp(rhoGhost);

  MultiFab bGhost(ba, dm, nDim3, 1);
  if (bTot.nGrow() >= 1) {
    MultiFab::Copy(bGhost, bTot, 0, 0, nDim3, 1);
    bGhost.FillBoundary(Geom(iLev).periodicity());
  } else {
    MultiFab::Copy(bGhost, bTot, 0, 0, nDim3, 0);
    apply_pe_zero_gradient_bc(iLev, bGhost);
    apply_fake2d_k_clamp(bGhost);
  }

  // Helper to fill ghosts of a 1-comp pressure MultiFab.
  auto fill_pe_ghosts = [&](MultiFab& mf) {
    mf.FillBoundary(Geom(iLev).periodicity());
    if (iLev > 0) {
      const MultiFab* cPtr =
          coarseState ? coarseState
                      : ((&state == &centerPeState[iLev])
                             ? &centerPeState[iLev - 1]
                             : nullptr);
      const int cComp = coarseState ? iCoarseComp : iStateComp;
      if (cPtr && !cPtr->empty()) {
        MultiFab coarseAlias(const_cast<MultiFab&>(*cPtr), amrex::make_alias,
                             cComp, 1);
        fill_fine_lev_bny_from_coarse(coarseAlias, mf, 0, 1, refRatio(iLev - 1),
                                      Geom(iLev - 1), Geom(iLev),
                                      cell_status(iLev), *get_cell_interp());
      }
    }
    apply_pe_zero_gradient_bc(iLev, mf);
    if (coupledDirichlet && &state == &centerPeState[iLev] && iStateComp == 0) {
      fill_ext_dir(cell_status(iLev), mf, 0, 1, &Pic::get_center_Pe, iLev,
                   &bcField);
    }
    apply_fake2d_k_clamp(mf);
    mf.FillBoundary(Geom(iLev).periodicity());
  };

  // Helper to convert cell-centered P -> T with 2 fresh ghost layers for
  // face-transverse slopes.
  auto fill_T_from_P = [&](const MultiFab& peMF, MultiFab& teMF) {
    for (MFIter mfi(teMF); mfi.isValid(); ++mfi) {
      const Box& fbox = mfi.fabbox();
      const Box gbox = amrex::grow(mfi.validbox(), 1);
      const Box box = fbox & gbox;
      const Array4<Real>& arrTe = teMF[mfi].array();
      const Array4<Real const>& arrPe = peMF[mfi].const_array();
      const Array4<Real const>& arrRho = rhoGhost[mfi].const_array();
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real ne = amrex::max(arrRho(i, j, k), rhoFloor);
        arrTe(i, j, k) = amrex::max(arrPe(i, j, k), peMinLocal) / ne;
      });
    }
    teMF.FillBoundary(Geom(iLev).periodicity());
    if (!Geom(iLev).isAllPeriodic()) {
      // Only fill outer physical domain boundaries, preserving the layer-1
      // coarse-fine values already computed from peMF/rhoGhost.
      Vector<BCRec> bcr(1);
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (Geom(iLev).isPeriodic(d)) {
          bcr[0].setLo(d, BCType::int_dir);
          bcr[0].setHi(d, BCType::int_dir);
        } else {
          bcr[0].setLo(d, BCType::foextrap);
          bcr[0].setHi(d, BCType::foextrap);
        }
      }
      GpuBndryFuncFab<FabFillNoOp> bfunc(FabFillNoOp{});
      PhysBCFunct<GpuBndryFuncFab<FabFillNoOp> > physbcf(Geom(iLev), bcr,
                                                         bfunc);
      physbcf(teMF, 0, 1, teMF.nGrowVect(), 0.0, 0);
    }
    apply_fake2d_k_clamp(teMF);
    teMF.FillBoundary(Geom(iLev).periodicity());
  };

  MultiFab peWork(ba, dm, 1, get_n_ghost());
  MultiFab::Copy(peWork, state, iStateComp, 0, 1, 0);
  fill_pe_ghosts(peWork);

  // Initial Te^n used to freeze face conductivities and saturation factor.
  MultiFab teWork(ba, dm, 1, get_n_ghost());
  fill_T_from_P(peWork, teWork);

  // Face-centered effective conductivities: (0) kN, (1) kT1, (2) kT2.
  amrex::Vector<MultiFab> faceCond(nDim);
  for (int d = 0; d < nDim; ++d) {
    if (!activeDir[d])
      continue;
    faceCond[d].define(amrex::convert(ba, IntVect::TheDimensionVector(d)), dm,
                       3, 0);
    const int t1 = (d + 1) % 3;
    const int t2 = (d + 2) % 3;
    const bool actT1 = (t1 < nDim && activeDir[t1]);
    const bool actT2 = (t2 < nDim && activeDir[t2]);
    const Real invD = invDxArr[d];
    const Real invT1 = actT1 ? invDxArr[t1] : 0.0;
    const Real invT2 = actT2 ? invDxArr[t2] : 0.0;
    const int di = (d == ix_) ? 1 : 0;
    const int dj = (d == iy_) ? 1 : 0;
    const int dk = (d == iz_) ? 1 : 0;
    const int t1i = (t1 == ix_) ? 1 : 0;
    const int t1j = (t1 == iy_) ? 1 : 0;
    const int t1k = (t1 == iz_) ? 1 : 0;
    const int t2i = (t2 == ix_) ? 1 : 0;
    const int t2j = (t2 == iy_) ? 1 : 0;
    const int t2k = (t2 == iz_) ? 1 : 0;

    const bool periodicD = Geom(iLev).isPeriodic(d);
    const bool maskLoWall =
        !periodicD &&
        !(coupledDirichlet && bcField.face(d, 0) == FieldBC::coupled);
    const bool maskHiWall =
        !periodicD &&
        !(coupledDirichlet && bcField.face(d, 1) == FieldBC::coupled);
    const int fDomLo = domLo[d];
    const int fDomHi = domHi[d] + 1;

    for (MFIter mfi(faceCond[d]); mfi.isValid(); ++mfi) {
      const Box& fbox = mfi.validbox();
      const Array4<Real>& arrK = faceCond[d][mfi].array();
      const Array4<Real const>& arrT = teWork[mfi].const_array();
      const Array4<Real const>& arrRho = rhoGhost[mfi].const_array();
      const Array4<Real const>& arrB = bGhost[mfi].const_array();

      parallel_for_omp(fbox, [=] AMREX_GPU_DEVICE(int i, int j,
                                                  int k) noexcept {
        const int idxD = (d == ix_) ? i : ((d == iy_) ? j : k);
        if ((maskLoWall && idxD <= fDomLo) || (maskHiWall && idxD >= fDomHi)) {
          arrK(i, j, k, 0) = 0.0;
          arrK(i, j, k, 1) = 0.0;
          arrK(i, j, k, 2) = 0.0;
          return;
        }

        const int iL = i - di, jL = j - dj, kL = k - dk;
        const Real teL = amrex::max(arrT(iL, jL, kL), static_cast<Real>(0.0));
        const Real teR = amrex::max(arrT(i, j, k), static_cast<Real>(0.0));
        const Real teF = 0.5 * (teL + teR);
        Real kappa = (teF > 0.0) ? kappa0Eff * std::pow(teF, 2.5)
                                 : static_cast<Real>(0.0);

        Real bhat[3] = { 0.0, 0.0, 0.0 };
        Real fEff = 0.0;
        if (fAlign > 0.0) {
          const Real bx = 0.5 * (arrB(iL, jL, kL, ix_) + arrB(i, j, k, ix_));
          const Real by = 0.5 * (arrB(iL, jL, kL, iy_) + arrB(i, j, k, iy_));
          const Real bz = 0.5 * (arrB(iL, jL, kL, iz_) + arrB(i, j, k, iz_));
          const Real b2 = bx * bx + by * by + bz * bz;
          if (b2 > 0.0) {
            const Real invB = 1.0 / std::sqrt(b2);
            bhat[ix_] = bx * invB;
            bhat[iy_] = by * invB;
            bhat[iz_] = bz * invB;
            fEff = fAlign * b2 / (b2 + bMin2);
          }
        }

        Real kN = kappa * (fEff * bhat[d] * bhat[d] + (1.0 - fEff));
        Real kT1 = kappa * fEff * bhat[d] * bhat[t1];
        Real kT2 = kappa * fEff * bhat[d] * bhat[t2];

        if (fluxLim > 0.0 && teF > 0.0 && kappa > 0.0) {
          const Real gN = (teR - teL) * invD;
          Real gT1 = 0.0, gT2 = 0.0;
          if (actT1) {
            const Real gL = 0.5 *
                            (arrT(iL + t1i, jL + t1j, kL + t1k) -
                             arrT(iL - t1i, jL - t1j, kL - t1k)) *
                            invT1;
            const Real gR = 0.5 *
                            (arrT(i + t1i, j + t1j, k + t1k) -
                             arrT(i - t1i, j - t1j, k - t1k)) *
                            invT1;
            gT1 = pe_slope(3, gL, gR);
          }
          if (actT2) {
            const Real gL = 0.5 *
                            (arrT(iL + t2i, jL + t2j, kL + t2k) -
                             arrT(iL - t2i, jL - t2j, kL - t2k)) *
                            invT2;
            const Real gR = 0.5 *
                            (arrT(i + t2i, j + t2j, k + t2k) -
                             arrT(i - t2i, j - t2j, k - t2k)) *
                            invT2;
            gT2 = pe_slope(3, gL, gR);
          }
          const Real qAbs = std::abs(kN * gN + kT1 * gT1 + kT2 * gT2);
          if (qAbs > 0.0) {
            const Real neF = amrex::max(
                0.5 * (arrRho(iL, jL, kL) + arrRho(i, j, k)), rhoFloor);
            const Real qSat =
                fluxLim * neF * teF * std::sqrt(teF * massRatioPe);
            if (qAbs > qSat) {
              const Real scale = qSat / qAbs;
              kN *= scale;
              kT1 *= scale;
              kT2 *= scale;
            }
          }
        }

        arrK(i, j, k, 0) = kN;
        arrK(i, j, k, 1) = kT1;
        arrK(i, j, k, 2) = kT2;
      });
    }
  }

  // Maximum diffusion rate across cells on this level.
  Real maxDiffRate = 0.0;
  {
    auto const& rhoArrays = rhoGhost.const_arrays();
    auto const& kX =
        activeDir[ix_] ? faceCond[ix_].const_arrays() : rhoGhost.const_arrays();
    auto const& kY =
        activeDir[iy_] ? faceCond[iy_].const_arrays() : rhoGhost.const_arrays();
    auto const& kZ = (nDim > 2 && activeDir[iz_]) ? faceCond[iz_].const_arrays()
                                                  : rhoGhost.const_arrays();
    const bool actX = activeDir[ix_];
    const bool actY = activeDir[iy_];
    const bool actZ = (nDim > 2 && activeDir[iz_]);
    const Real invX2 = invDxArr[ix_] * invDxArr[ix_];
    const Real invY2 = invDxArr[iy_] * invDxArr[iy_];
    const Real invZ2 = invDxArr[iz_] * invDxArr[iz_];

    maxDiffRate = ParReduce(
        TypeList<ReduceOpMax>{}, TypeList<Real>{}, peWork, IntVect(0),
        [=] AMREX_GPU_DEVICE(int nb, int i, int j, int k)
            noexcept -> GpuTuple<Real> {
              const Real ne = amrex::max(rhoArrays[nb](i, j, k), rhoFloor);
              Real lam = 0.0;
              if (actX)
                lam += (kX[nb](i, j, k, 0) + kX[nb](i + 1, j, k, 0)) * invX2;
              if (actY)
                lam += (kY[nb](i, j, k, 0) + kY[nb](i, j + 1, k, 0)) * invY2;
              if (actZ)
                lam += (kZ[nb](i, j, k, 0) + kZ[nb](i, j, k + 1, 0)) * invZ2;
              return energyScale * lam / ne;
            });
    ParallelDescriptor::ReduceRealMax(maxDiffRate);
  }
  if (maxDiffRate <= 0.0)
    return;

  // Evaluate L(P) = -energyScale * div(q(T(P))) into rhsMF.
  MultiFab fluxTmp[nDim3];
  for (int d = 0; d < nDim; ++d) {
    if (activeDir[d]) {
      fluxTmp[d].define(amrex::convert(ba, IntVect::TheDimensionVector(d)), dm,
                        1, 0);
    }
  }

  auto eval_rhs = [&](const MultiFab& peMF, MultiFab& rhsMF) {
    fill_T_from_P(peMF, teWork);
    for (int d = 0; d < nDim; ++d) {
      if (!activeDir[d])
        continue;
      const int t1 = (d + 1) % 3;
      const int t2 = (d + 2) % 3;
      const bool actT1 = (t1 < nDim && activeDir[t1]);
      const bool actT2 = (t2 < nDim && activeDir[t2]);
      const Real invD = invDxArr[d];
      const Real invT1 = actT1 ? invDxArr[t1] : 0.0;
      const Real invT2 = actT2 ? invDxArr[t2] : 0.0;
      const int di = (d == ix_) ? 1 : 0;
      const int dj = (d == iy_) ? 1 : 0;
      const int dk = (d == iz_) ? 1 : 0;
      const int t1i = (t1 == ix_) ? 1 : 0;
      const int t1j = (t1 == iy_) ? 1 : 0;
      const int t1k = (t1 == iz_) ? 1 : 0;
      const int t2i = (t2 == ix_) ? 1 : 0;
      const int t2j = (t2 == iy_) ? 1 : 0;
      const int t2k = (t2 == iz_) ? 1 : 0;

      for (MFIter mfi(fluxTmp[d]); mfi.isValid(); ++mfi) {
        const Box& fbox = mfi.validbox();
        const Array4<Real>& arrQ = fluxTmp[d][mfi].array();
        const Array4<Real const>& arrK = faceCond[d][mfi].const_array();
        const Array4<Real const>& arrT = teWork[mfi].const_array();
        parallel_for_omp(
            fbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              const Real kN = arrK(i, j, k, 0);
              const Real kT1 = arrK(i, j, k, 1);
              const Real kT2 = arrK(i, j, k, 2);
              if (kN == 0.0 && kT1 == 0.0 && kT2 == 0.0) {
                arrQ(i, j, k) = 0.0;
                return;
              }
              const int iL = i - di, jL = j - dj, kL = k - dk;
              const Real gN = (arrT(i, j, k) - arrT(iL, jL, kL)) * invD;
              Real gT1 = 0.0, gT2 = 0.0;
              if (actT1 && kT1 != 0.0) {
                const Real gL = 0.5 *
                                (arrT(iL + t1i, jL + t1j, kL + t1k) -
                                 arrT(iL - t1i, jL - t1j, kL - t1k)) *
                                invT1;
                const Real gR = 0.5 *
                                (arrT(i + t1i, j + t1j, k + t1k) -
                                 arrT(i - t1i, j - t1j, k - t1k)) *
                                invT1;
                gT1 = pe_slope(3, gL, gR);
              }
              if (actT2 && kT2 != 0.0) {
                const Real gL = 0.5 *
                                (arrT(iL + t2i, jL + t2j, kL + t2k) -
                                 arrT(iL - t2i, jL - t2j, kL - t2k)) *
                                invT2;
                const Real gR = 0.5 *
                                (arrT(i + t2i, j + t2j, k + t2k) -
                                 arrT(i - t2i, j - t2j, k - t2k)) *
                                invT2;
                gT2 = pe_slope(3, gL, gR);
              }
              arrQ(i, j, k) = -(kN * gN + kT1 * gT1 + kT2 * gT2);
            });
      }
    }

    const bool actX = activeDir[ix_];
    const bool actY = activeDir[iy_];
    const bool actZ = (nDim > 2 && activeDir[iz_]);
    const Real invX = invDxArr[ix_];
    const Real invY = invDxArr[iy_];
    const Real invZ = invDxArr[iz_];

    for (MFIter mfi(rhsMF); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrRhs = rhsMF[mfi].array();
      const Array4<Real const> qX =
          actX ? fluxTmp[ix_][mfi].const_array() : Array4<Real const>();
      const Array4<Real const> qY =
          actY ? fluxTmp[iy_][mfi].const_array() : Array4<Real const>();
      const Array4<Real const> qZ =
          actZ ? fluxTmp[iz_][mfi].const_array() : Array4<Real const>();
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        Real divQ = 0.0;
        if (actX)
          divQ += (qX(i + 1, j, k) - qX(i, j, k)) * invX;
        if (actY)
          divQ += (qY(i, j + 1, k) - qY(i, j, k)) * invY;
        if (actZ)
          divQ += (qZ(i, j, k + 1) - qZ(i, j, k)) * invZ;
        arrRhs(i, j, k) = -energyScale * divQ;
      });
    }
  };

  // Symmetric conservative face-flux point-implicit pass over subDt.
  auto apply_point_implicit = [&](Real subDt) {
    fill_T_from_P(peWork, teWork);
    for (int d = 0; d < nDim; ++d) {
      if (!activeDir[d])
        continue;
      const int t1 = (d + 1) % 3;
      const int t2 = (d + 2) % 3;
      const bool actT1 = (t1 < nDim && activeDir[t1]);
      const bool actT2 = (t2 < nDim && activeDir[t2]);
      const Real invD = invDxArr[d];
      const Real invD2 = invD * invD;
      const Real invT1 = actT1 ? invDxArr[t1] : 0.0;
      const Real invT2 = actT2 ? invDxArr[t2] : 0.0;
      const int di = (d == ix_) ? 1 : 0;
      const int dj = (d == iy_) ? 1 : 0;
      const int dk = (d == iz_) ? 1 : 0;
      const int t1i = (t1 == ix_) ? 1 : 0;
      const int t1j = (t1 == iy_) ? 1 : 0;
      const int t1k = (t1 == iz_) ? 1 : 0;
      const int t2i = (t2 == ix_) ? 1 : 0;
      const int t2j = (t2 == iy_) ? 1 : 0;
      const int t2k = (t2 == iz_) ? 1 : 0;

      for (MFIter mfi(fluxTmp[d]); mfi.isValid(); ++mfi) {
        const Box& fbox = mfi.validbox();
        const Array4<Real>& arrQ = fluxTmp[d][mfi].array();
        const Array4<Real const>& arrK = faceCond[d][mfi].const_array();
        const Array4<Real const>& arrT = teWork[mfi].const_array();
        const Array4<Real const>& arrRho = rhoGhost[mfi].const_array();
        parallel_for_omp(
            fbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              const Real kN = arrK(i, j, k, 0);
              const Real kT1 = arrK(i, j, k, 1);
              const Real kT2 = arrK(i, j, k, 2);
              if (kN == 0.0 && kT1 == 0.0 && kT2 == 0.0) {
                arrQ(i, j, k) = 0.0;
                return;
              }
              const int iL = i - di, jL = j - dj, kL = k - dk;
              const Real gN = (arrT(i, j, k) - arrT(iL, jL, kL)) * invD;
              Real gT1 = 0.0, gT2 = 0.0;
              if (actT1 && kT1 != 0.0) {
                const Real gL = 0.5 *
                                (arrT(iL + t1i, jL + t1j, kL + t1k) -
                                 arrT(iL - t1i, jL - t1j, kL - t1k)) *
                                invT1;
                const Real gR = 0.5 *
                                (arrT(i + t1i, j + t1j, k + t1k) -
                                 arrT(i - t1i, j - t1j, k - t1k)) *
                                invT1;
                gT1 = pe_slope(3, gL, gR);
              }
              if (actT2 && kT2 != 0.0) {
                const Real gL = 0.5 *
                                (arrT(iL + t2i, jL + t2j, kL + t2k) -
                                 arrT(iL - t2i, jL - t2j, kL - t2k)) *
                                invT2;
                const Real gR = 0.5 *
                                (arrT(i + t2i, j + t2j, k + t2k) -
                                 arrT(i - t2i, j - t2j, k - t2k)) *
                                invT2;
                gT2 = pe_slope(3, gL, gR);
              }
              const Real neL = amrex::max(arrRho(iL, jL, kL), rhoFloor);
              const Real neR = amrex::max(arrRho(i, j, k), rhoFloor);
              const Real lamF = kN * invD2 * (1.0 / neL + 1.0 / neR);
              const Real denom = 1.0 + subDt * energyScale * lamF;
              arrQ(i, j, k) = -(kN * gN + kT1 * gT1 + kT2 * gT2) / denom;
            });
      }
    }

    const bool actX = activeDir[ix_];
    const bool actY = activeDir[iy_];
    const bool actZ = (nDim > 2 && activeDir[iz_]);
    const Real invX = invDxArr[ix_];
    const Real invY = invDxArr[iy_];
    const Real invZ = invDxArr[iz_];

    for (MFIter mfi(peWork); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrP = peWork[mfi].array();
      const Array4<Real const> qX =
          actX ? fluxTmp[ix_][mfi].const_array() : Array4<Real const>();
      const Array4<Real const> qY =
          actY ? fluxTmp[iy_][mfi].const_array() : Array4<Real const>();
      const Array4<Real const> qZ =
          actZ ? fluxTmp[iz_][mfi].const_array() : Array4<Real const>();
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        Real divQ = 0.0;
        if (actX)
          divQ += (qX(i + 1, j, k) - qX(i, j, k)) * invX;
        if (actY)
          divQ += (qY(i, j + 1, k) - qY(i, j, k)) * invY;
        if (actZ)
          divQ += (qZ(i, j, k + 1) - qZ(i, j, k)) * invZ;
        arrP(i, j, k) =
            amrex::max(arrP(i, j, k) - subDt * energyScale * divQ, peMinLocal);
      });
    }
    fill_pe_ghosts(peWork);
  };

  Real dtRem = dt;

  if (heatCondMethod == "rkl2") {
    const int sReq = std::max(
        2, static_cast<int>(std::ceil(
               0.5 * (-1.0 + std::sqrt(9.0 + 20.0 * dt * maxDiffRate)))));
    const int s = std::min(sReq, nCondSubcycleMax);
    const Real dtRkl =
        (s < sReq) ? amrex::min(dt, 0.2 * static_cast<Real>(s * s + s - 2) /
                                        maxDiffRate)
                   : dt;

    if (s >= 2 && dtRkl > 0.0) {
      MultiFab Y0(ba, dm, 1, 0);
      MultiFab Yjm2(ba, dm, 1, 0);
      MultiFab Yjm1(ba, dm, 1, get_n_ghost());
      MultiFab L0(ba, dm, 1, 0);
      MultiFab Ljm1(ba, dm, 1, 0);

      MultiFab::Copy(Y0, peWork, 0, 0, 1, 0);
      MultiFab::Copy(Yjm2, Y0, 0, 0, 1, 0);
      eval_rhs(peWork, L0);

      const Real sReal = static_cast<Real>(s);
      const Real omega1 = 4.0 / (sReal * sReal + sReal - 2.0);
      auto b_coef = [](int j) -> Real {
        if (j <= 2)
          return 1.0 / 3.0;
        const Real jr = static_cast<Real>(j);
        return (jr * jr + jr - 2.0) / (2.0 * jr * (jr + 1.0));
      };

      const Real mu1Tilde = b_coef(1) * omega1;
      MultiFab::LinComb(Yjm1, 1.0, Y0, 0, mu1Tilde * dtRkl, L0, 0, 0, 1, 0);
      fill_pe_ghosts(Yjm1);

      for (int j = 2; j <= s; ++j) {
        eval_rhs(Yjm1, Ljm1);
        const Real jr = static_cast<Real>(j);
        const Real bj = b_coef(j);
        const Real bjm1 = b_coef(j - 1);
        const Real bjm2 = b_coef(j - 2);
        const Real ajm1 = 1.0 - bjm1;
        const Real muj = ((2.0 * jr - 1.0) / jr) * (bj / bjm1);
        const Real nuj = -((jr - 1.0) / jr) * (bj / bjm2);
        const Real mujTilde = muj * omega1;
        const Real gammajTilde = -ajm1 * mujTilde;

        for (MFIter mfi(peWork); mfi.isValid(); ++mfi) {
          const Box& box = mfi.validbox();
          const Array4<Real>& arrY = peWork[mfi].array();
          const Array4<Real const>& arrY0 = Y0[mfi].const_array();
          const Array4<Real const>& arrYm1 = Yjm1[mfi].const_array();
          const Array4<Real const>& arrYm2 = Yjm2[mfi].const_array();
          const Array4<Real const>& arrL0 = L0[mfi].const_array();
          const Array4<Real const>& arrLm1 = Ljm1[mfi].const_array();
          parallel_for_omp(
              box, [=] AMREX_GPU_DEVICE(int i, int jCell, int k) noexcept {
                arrY(i, jCell, k) = muj * arrYm1(i, jCell, k) +
                                    nuj * arrYm2(i, jCell, k) +
                                    (1.0 - muj - nuj) * arrY0(i, jCell, k) +
                                    mujTilde * dtRkl * arrLm1(i, jCell, k) +
                                    gammajTilde * dtRkl * arrL0(i, jCell, k);
              });
        }
        if (j < s) {
          MultiFab::Copy(Yjm2, Yjm1, 0, 0, 1, 0);
          MultiFab::Copy(Yjm1, peWork, 0, 0, 1, 0);
          fill_pe_ghosts(Yjm1);
        }
      }

      for (MFIter mfi(peWork); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real>& arrP = peWork[mfi].array();
        parallel_for_omp(
            box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              arrP(i, j, k) = amrex::max(arrP(i, j, k), peMinLocal);
            });
      }
      fill_pe_ghosts(peWork);
      dtRem -= dtRkl;
    }
  } else if (heatCondMethod == "subcycle") {
    MultiFab rhs(ba, dm, 1, 0);
    const Real dtCFL = 0.45 / maxDiffRate;
    int pass = 0;
    while (dtRem > 0.0 && pass < nCondSubcycleMax) {
      ++pass;
      const Real subDt = amrex::min(dtRem, dtCFL);
      if (pass == nCondSubcycleMax && subDt < dtRem) {
        break;
      }
      eval_rhs(peWork, rhs);
      for (MFIter mfi(peWork); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real>& arrP = peWork[mfi].array();
        const Array4<Real const>& arrRhs = rhs[mfi].const_array();
        parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j,
                                                   int k) noexcept {
          arrP(i, j, k) =
              amrex::max(arrP(i, j, k) + subDt * arrRhs(i, j, k), peMinLocal);
        });
      }
      fill_pe_ghosts(peWork);
      dtRem -= subDt;
    }
  }

  if (dtRem > 0.0) {
    const int nPasses = (heatCondMethod == "point-implicit") ? nCondIter : 1;
    const Real subDt = dtRem / static_cast<Real>(nPasses);
    for (int pass = 0; pass < nPasses; ++pass) {
      apply_point_implicit(subDt);
    }
  }

  MultiFab::Copy(state, peWork, 0, iStateComp, 1, 0);
}

//==========================================================
// Electron heat conduction:
//   dPe/dt = (gamma_e-1) div(kappa_hat . grad(Te)) = -(gamma_e-1) div(q),
// with q = -kappa_hat . grad(Te) and Spitzer kappa = kappa0 Te^2.5.
void Pic::apply_electron_heat_conduction(int iLev, Real dt) {
  BL_PROFILE("Pic::apply_electron_heat_conduction");
  if (heatCondKappa0 <= 0.0)
    return;

  centerB[iLev].FillBoundary(Geom(iLev).periodicity());
  const MultiFab& bTot = total_center_B(centerB[iLev], iLev);

  if (!useAnisoPe) {
    solve_heat_conduction_scalar(iLev, dt, centerPeState[iLev], 0,
                                 centerPeRho[iLev], bTot, 1.0,
                                 electronGamma - 1.0);
    return;
  }

  // Gyrotropic CGL heat conduction: evolve parallel thermal energy (PePar/2)
  // and perpendicular thermal energy (PePerp) with kappaScale = 2/3,
  // energyScale = 2 (for PePar) and 1 (for PePerp), conserving total thermal
  // energy (PePar + 2*PePerp)/2 = 3/2 * Pe.
  const BoxArray& ba = centerPeState[iLev].boxArray();
  const DistributionMapping& dm = centerPeState[iLev].DistributionMap();
  const Real peMinLocal = peMin;

  MultiFab pePerpMF(ba, dm, 1, get_n_ghost());
  for (MFIter mfi(pePerpMF); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrPerp = pePerpMF[mfi].array();
    const Array4<Real const>& arrS = centerPeState[iLev][mfi].const_array();
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const Real pe = amrex::max(arrS(i, j, k, 0), peMinLocal);
      const Real pPar = amrex::max(arrS(i, j, k, 1), peMinLocal);
      arrPerp(i, j, k) = amrex::max(0.5 * (3.0 * pe - pPar), peMinLocal);
    });
  }

  MultiFab pePerpCoarse;
  if (iLev > 0) {
    pePerpCoarse.define(centerPeState[iLev - 1].boxArray(),
                        centerPeState[iLev - 1].DistributionMap(), 1,
                        centerPeState[iLev - 1].nGrow());
    for (MFIter mfi(pePerpCoarse); mfi.isValid(); ++mfi) {
      const Box& fbox = mfi.fabbox();
      const Array4<Real>& arrC = pePerpCoarse[mfi].array();
      const Array4<Real const>& arrS =
          centerPeState[iLev - 1][mfi].const_array();
      parallel_for_omp(fbox,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         const Real pe = amrex::max(arrS(i, j, k, 0), peMinLocal);
                         const Real pPar =
                             amrex::max(arrS(i, j, k, 1), peMinLocal);
                         arrC(i, j, k) =
                             amrex::max(0.5 * (3.0 * pe - pPar), peMinLocal);
                       });
    }
  }

  solve_heat_conduction_scalar(
      iLev, dt, centerPeState[iLev], 1, centerPeRho[iLev], bTot, 2.0 / 3.0, 2.0,
      iLev > 0 ? &centerPeState[iLev - 1] : nullptr, 1);
  solve_heat_conduction_scalar(iLev, dt, pePerpMF, 0, centerPeRho[iLev], bTot,
                               2.0 / 3.0, 1.0,
                               iLev > 0 ? &pePerpCoarse : nullptr, 0);

  // Conservative mirror-force heat-flux exchange between PePar and PePerp:
  //   dPe_par/dt += +2 * q_perp * div(bhat),  dPe_perp/dt += -q_perp * div(bhat)
  // which leaves Pe = (PePar + 2*PePerp)/3 strictly invariant.
  const Real fAlign = fieldAlignedConduction ? fieldAlignedFraction : 0.0;
  if (fAlign > 0.0) {
    const IntVect domLen = Geom(iLev).Domain().length();
    const Real* invDxGeom = Geom(iLev).InvCellSize();
    const Real invDxX = (domLen[ix_] > 1) ? invDxGeom[ix_] : 0.0;
    const Real invDxY = (domLen[iy_] > 1) ? invDxGeom[iy_] : 0.0;
    const Real invDxZ = (nDim > 2 && domLen[iz_] > 1) ? invDxGeom[iz_] : 0.0;
    const Real bMin2 = fieldAlignedBMin * fieldAlignedBMin;
    const Real rhoFloor = rhoMinOhm;
    const Real kappa0Perp = heatCondKappa0 * (2.0 / 3.0);

    MultiFab bHatMF(ba, dm, nDim3, 1);
    MultiFab tePerpMF(ba, dm, 1, 1);
    for (MFIter mfi(bHatMF); mfi.isValid(); ++mfi) {
      const Box& fbox = mfi.fabbox();
      const Box gbox = amrex::grow(mfi.validbox(), 1);
      const Box box = fbox & gbox;
      const Array4<Real>& arrBh = bHatMF[mfi].array();
      const Array4<Real const>& arrB = bTot[mfi].const_array();
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real bx = arrB(i, j, k, ix_);
        const Real by = arrB(i, j, k, iy_);
        const Real bz = arrB(i, j, k, iz_);
        const Real b2 = bx * bx + by * by + bz * bz;
        if (b2 > 0.0) {
          const Real invB = 1.0 / std::sqrt(b2);
          const Real wB = fAlign * b2 / (b2 + bMin2);
          arrBh(i, j, k, ix_) = wB * bx * invB;
          arrBh(i, j, k, iy_) = wB * by * invB;
          arrBh(i, j, k, iz_) = wB * bz * invB;
        } else {
          arrBh(i, j, k, ix_) = 0.0;
          arrBh(i, j, k, iy_) = 0.0;
          arrBh(i, j, k, iz_) = 0.0;
        }
      });
    }
    bHatMF.FillBoundary(Geom(iLev).periodicity());
    apply_pe_zero_gradient_bc(iLev, bHatMF);
    apply_fake2d_k_clamp(bHatMF);

    for (MFIter mfi(tePerpMF); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrT = tePerpMF[mfi].array();
      const Array4<Real const>& arrP = pePerpMF[mfi].const_array();
      const Array4<Real const>& arrRho = centerPeRho[iLev][mfi].const_array();
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real ne = amrex::max(arrRho(i, j, k), rhoFloor);
        arrT(i, j, k) = amrex::max(arrP(i, j, k), peMinLocal) / ne;
      });
    }
    tePerpMF.FillBoundary(Geom(iLev).periodicity());
    apply_pe_zero_gradient_bc(iLev, tePerpMF);
    apply_fake2d_k_clamp(tePerpMF);

    for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrS = centerPeState[iLev][mfi].array();
      const Array4<Real>& arrPerp = pePerpMF[mfi].array();
      const Array4<Real const>& arrBh = bHatMF[mfi].const_array();
      const Array4<Real const>& arrTp = tePerpMF[mfi].const_array();
      const Array4<Real const>& arrRho = centerPeRho[iLev][mfi].const_array();

      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real ne = amrex::max(arrRho(i, j, k), rhoFloor);
        const Real pPar = amrex::max(arrS(i, j, k, 1), peMinLocal);
        const Real pPerp = amrex::max(arrPerp(i, j, k), peMinLocal);
        const Real tPar = pPar / ne;
        const Real kap = (tPar > 0.0) ? kappa0Perp * std::pow(tPar, 2.5) : 0.0;
        const Real gX = (invDxX > 0.0)
                            ? 0.5 * (arrTp(i + 1, j, k) - arrTp(i - 1, j, k)) *
                                  invDxX
                            : 0.0;
        const Real gY = (invDxY > 0.0)
                            ? 0.5 * (arrTp(i, j + 1, k) - arrTp(i, j - 1, k)) *
                                  invDxY
                            : 0.0;
        const Real gZ =
            (invDxZ > 0.0)
                ? 0.5 * (arrTp(i, j, k + 1) - arrTp(i, j, k - 1)) * invDxZ
                : 0.0;
        const Real bDotG = arrBh(i, j, k, ix_) * gX + arrBh(i, j, k, iy_) * gY +
                           arrBh(i, j, k, iz_) * gZ;
        const Real qPerp = -kap * bDotG;
        Real divBh = 0.0;
        if (invDxX > 0.0)
          divBh +=
              0.5 * (arrBh(i + 1, j, k, ix_) - arrBh(i - 1, j, k, ix_)) * invDxX;
        if (invDxY > 0.0)
          divBh +=
              0.5 * (arrBh(i, j + 1, k, iy_) - arrBh(i, j - 1, k, iy_)) * invDxY;
        if (invDxZ > 0.0)
          divBh +=
              0.5 * (arrBh(i, j, k + 1, iz_) - arrBh(i, j, k - 1, iz_)) * invDxZ;
        Real dE = dt * qPerp * divBh;
        // Clamp transfer so neither parallel nor perpendicular pressure crosses
        // peMinLocal, keeping Pe = (PePar + 2*PePerp)/3 exact.
        const Real dEMax = pPerp - peMinLocal;
        const Real dEMin = -0.5 * (pPar - peMinLocal);
        dE = std::clamp(dE, dEMin, dEMax);
        arrS(i, j, k, 1) = pPar + 2.0 * dE;
        arrPerp(i, j, k) = pPerp - dE;
      });
    }
  }

  for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrS = centerPeState[iLev][mfi].array();
    const Array4<Real const>& arrPerp = pePerpMF[mfi].const_array();
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const Real pPar = amrex::max(arrS(i, j, k, 1), peMinLocal);
      const Real pPerp = amrex::max(arrPerp(i, j, k), peMinLocal);
      arrS(i, j, k, 0) = (pPar + 2.0 * pPerp) / 3.0;
      arrS(i, j, k, 1) = pPar;
    });
  }
}

//==========================================================
// Marginal-stability (firehose + whistler/mirror) & baseline isotropization for
// gyrotropic CGL electrons. Leaves Pe = (PePar + 2*PePerp)/3 strictly invariant.
void Pic::apply_electron_isotropization(int iLev, Real dt) {
  if (!useAnisoPe)
    return;

  centerB[iLev].FillBoundary(Geom(iLev).periodicity());
  const MultiFab& bTot = total_center_B(centerB[iLev], iLev);
  const Real bMin2 = fieldAlignedBMin * fieldAlignedBMin;
  const Real peMinLocal = peMin;
  const bool useLimit = usePeInstabilityLimit;
  const Real fhCoef = peFirehoseCoef;
  const Real whCoef = peWhistlerCoef;
  const Real whExp = peWhistlerExp;
  const Real tauInst = peAnisoRelaxTime;
  const Real baseRate = peAnisoBaseRate;

  for (MFIter mfi(centerPeState[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrS = centerPeState[iLev][mfi].array();
    const Array4<Real const>& arrB = bTot[mfi].const_array();

    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const Real pe = amrex::max(arrS(i, j, k, 0), peMinLocal);
      Real pePar = amrex::max(arrS(i, j, k, 1), peMinLocal);
      const Real peParCeiling = amrex::max(peMinLocal, 3.0 * pe - 2.0 * peMinLocal);
      pePar = amrex::min(pePar, peParCeiling);

      const Real bx = arrB(i, j, k, ix_);
      const Real by = arrB(i, j, k, iy_);
      const Real bz = arrB(i, j, k, iz_);
      const Real b2 = bx * bx + by * by + bz * bz;
      if (b2 <= 0.0) {
        arrS(i, j, k, 0) = pe;
        arrS(i, j, k, 1) = pe;
        return;
      }

      const Real pB = b2 / (8.0 * dPI);
      const Real wB = b2 / (b2 + bMin2);

      if (useLimit && pB > 0.0) {
        // Firehose upper bound: beta_par - beta_perp <= fhCoef
        // => PePar - PePerp <= fhCoef * P_B => PePar <= Pe + (2/3)*fhCoef*P_B.
        const Real pParMax = amrex::min(pe + (2.0 / 3.0) * fhCoef * pB,
                                        peParCeiling);

        // Whistler/mirror lower bound: PePerp/PePar - 1 <= whCoef / beta_par^whExp
        // i.e. g(p) = 3*(Pe - p)/(2*p) - whCoef * (P_B / p)^whExp <= 0.
        Real pParMin = peMinLocal;
        if (whCoef > 0.0 && pePar < pe) {
          const Real pFloor = amrex::max(peMinLocal, 1.0e-30 * pe);
          auto g_wh = [&](const Real p) noexcept -> Real {
            return 1.5 * (pe - p) / p - whCoef * std::pow(pB / p, whExp);
          };
          if (g_wh(pePar) > 0.0) {
            // Bisection + Newton bracketed on [pFloor, pe].
            Real lo = pFloor;
            Real hi = pe;
            Real p = amrex::max(pePar, lo);
            if (g_wh(lo) <= 0.0) {
              pParMin = lo;
            } else {
              for (int it = 0; it < 12; ++it) {
                const Real g = g_wh(p);
                if (g > 0.0) {
                  lo = p;
                } else {
                  hi = p;
                }
                const Real dg =
                    -1.5 * pe / (p * p) +
                    whCoef * whExp * std::pow(pB / p, whExp) / p;
                Real pNext = (dg < 0.0) ? (p - g / dg) : (0.5 * (lo + hi));
                if (!(pNext > lo && pNext < hi)) {
                  pNext = 0.5 * (lo + hi);
                }
                p = pNext;
              }
              pParMin = std::clamp(p, peMinLocal, pe);
            }
          }
        }

        const Real pClamp = std::clamp(pePar, pParMin, pParMax);
        if (pClamp != pePar) {
          if (tauInst <= 0.0) {
            pePar = pClamp;
          } else {
            const Real r = dt / tauInst;
            pePar = (pePar + r * pClamp) / (1.0 + r);
          }
        }
      }

      if (baseRate > 0.0) {
        const Real rBase = dt * baseRate;
        pePar = (pePar + rBase * pe) / (1.0 + rBase);
      }

      // Smoothly isotropize as |B| drops below fieldAlignedBMin.
      pePar = pe + wB * (pePar - pe);
      pePar = std::clamp(pePar, peMinLocal, peParCeiling);

      arrS(i, j, k, 0) = pe;
      arrS(i, j, k, 1) = pePar;
    });
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
void Pic::apply_centerB_BC(int iLev, Real timeOverride, Real dtOverride) {
  apply_centerB_BC(iLev, centerB[iLev], timeOverride, dtOverride);
}

void Pic::apply_centerB_BC(int iLev, amrex::MultiFab& mfB, Real timeOverride,
                           Real dtOverride) {
  mfB.FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(cell_status(iLev), mfB, 0, mfB.nComp(), &Pic::get_center_B,
                   iLev, true, dtOverride, timeOverride);
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
                   iLev, true, dtOverride, timeOverride);
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
      parallel_for_omp(fbox,
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
void Pic::sync_emf_at_interface(Vector<MultiFab>& nodeEmf, Real timeOverride) {
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
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
                   iLev, false, 0.0, timeOverride);
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
void Pic::fill_ghost_emf_from_coarse(Vector<MultiFab>& nodeEmf,
                                     Real timeOverride) {
  BL_PROFILE("Pic::fill_ghost_emf_from_coarse");
  for (int iLev = 1; iLev <= finestLevel(); ++iLev) {
    fill_fine_lev_bny_from_coarse(
        nodeEmf[iLev - 1], nodeEmf[iLev], 0, nDim3, refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    nodeEmf[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_field_bc(node_status(iLev), nodeEmf[iLev], 0, nDim3, &Pic::get_node_E,
                   iLev, false, 0.0, timeOverride);
    if (use_body()) {
      apply_body_E_bc(nodeEmf[iLev], iLev);
    }
  }
}

void Pic::fill_interface_and_ghost_emf_from_coarse(Vector<MultiFab>& nodeEmf,
                                                   Real timeOverride) {
  BL_PROFILE("Pic::fill_interface_and_ghost_emf_from_coarse");
  for (int iLev = 1; iLev <= finestLevel(); ++iLev) {
    fill_fine_lev_edge_from_coarse(
        nodeEmf[iLev - 1], nodeEmf[iLev], 0, nDim3, refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    fill_fine_lev_bny_from_coarse(
        nodeEmf[iLev - 1], nodeEmf[iLev], 0, nDim3, refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    nodeEmf[iLev].OverrideSync(Geom(iLev).periodicity());
    nodeEmf[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_field_bc(node_status(iLev), nodeEmf[iLev], 0, nDim3, &Pic::get_node_E,
                   iLev, false, 0.0, timeOverride);
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
                            bool includeAmbi, Real timeOverride) {
  BL_PROFILE("Pic::faraday_stage_rhs");
  if (has_time_dep_B0()) {
    const Real tEnd = tc ? tc->get_time() : 0.0;
    const Real dtFull = tc ? tc->get_dt() : 0.0;
    const Real tStage = (timeOverride >= 0.0)
                            ? timeOverride
                            : std::max(0.0, tEnd - (1.0 - hstep) * dtFull);
    refresh_intrinsic_fields(tStage);
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    assemble_ohm_E(Bin[iLev], total_center_B(Bin[iLev], iLev), nodeEstage[iLev],
                   iLev, hstep, includeAmbi, timeOverride);
  }

  if (finestLevel() > 0) {
    sync_emf_at_interface(nodeEstage, timeOverride);
    fill_ghost_emf_from_coarse(nodeEstage, timeOverride);
  }

  if (has_active_inductive_E0()) {
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      // Ohm's law produces total E; subtract inductive E0 = -dA0/dt so
      // Faraday's law advances only B1 = B_tot - B0.
      if (!nodeE0[iLev].empty()) {
        const int nG = std::min(nodeEstage[iLev].nGrow(), nodeE0[iLev].nGrow());
        MultiFab::Subtract(nodeEstage[iLev], nodeE0[iLev], 0, 0, nDim3, nG);
      }
    }
    // Re-synchronize coarse-fine interface hanging edges and fine ghost nodes
    // on nodeEstage (= E_tot - E0) so coarse and fine curls share identical
    // interface edge circulations and preserve div(B1) = 0.
    for (int iLev = 1; iLev <= finestLevel(); ++iLev) {
      fill_fine_lev_edge_from_coarse(
          nodeEstage[iLev - 1], nodeEstage[iLev], 0, nDim3, refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
      fill_fine_lev_bny_from_coarse(
          nodeEstage[iLev - 1], nodeEstage[iLev], 0, nDim3, refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
      nodeEstage[iLev].OverrideSync(Geom(iLev).periodicity());
      nodeEstage[iLev].FillBoundary(Geom(iLev).periodicity());
    }
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
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
  Real gamma = mf_dot(s, 0, s, 0, nDim3, 0, owner.get());
  const Real gamma0 = gamma;

  for (int iter = 0; iter < maxIter && gamma > relTol * relTol * gamma0;
       ++iter) {
    apply_curl();
    const Real qq = mf_dot(q, 0, q, 0, nRelax, 0);
    if (qq <= 0)
      break;
    const Real a = gamma / qq;
    MultiFab::Saxpy(B, a, q, 0, 0, nRelax, 0);
    MultiFab::Saxpy(r, -a, q, 0, 0, nRelax, 0);
    apply_curlT();
    const Real gammaNew = mf_dot(s, 0, s, 0, nDim3, 0, owner.get());
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
  if (dt <= 0.0)
    return;
  if (nBSubcycle < 1)
    nBSubcycle = 1;
  Real subDt = dt / nBSubcycle;

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

  // CFL stability check for explicit resistive, hyper-resistive, upwind, and
  // whistler terms. Real-axis dissipative terms (resistivity,
  // hyper-resistivity, upwind) use the real-axis stability limit (2.785 for
  // RK4, 2.513 for SSPRK3), while the oscillatory Hall whistler mode lies on
  // the imaginary axis (2.785 for RK4, sqrt(3) ~ 1.732 for SSPRK3).
  const Real cflLimitReal = useRK4 ? 2.785 : 2.513;
  const Real cflLimitWhistler = useRK4 ? 2.785 : 1.732;
  const Real fieldSafety = 0.85;
  auto safe_ceil_int = [](Real val) -> int {
    if (!std::isfinite(val) ||
        val >= static_cast<Real>(std::numeric_limits<int>::max() / 2)) {
      return std::numeric_limits<int>::max() / 2;
    }
    if (val <= 1.0)
      return 1;
    return static_cast<int>(std::ceil(val));
  };

  Real maxRegionalEta = 0.0;
  for (const auto& cfg : regionalResistivityConfigs) {
    maxRegionalEta = amrex::max(maxRegionalEta, cfg.etaCode);
  }
  const Real maxEta = etaResistivity + maxRegionalEta;

  const Real tStepEnd = tc ? tc->get_time() : 0.0;
  const Real tStepStart = std::max(0.0, tStepEnd - dt);

  refresh_intrinsic_fields(tStepStart);

  // Ensure B ghost cells are fresh across all levels at tStepStart before CFL
  // checks and Faraday subcycling.
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    apply_centerB_BC(iLev, tStepStart);
  }

  const bool hasBody = use_body();

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    const bool hasFiner = (iLev < finestLevel());
    Real maxRegionalHyperLev = 0.0;
    for (const auto& cfg : regionalHyperResistivityConfigs) {
      if (iLev < static_cast<int>(cfg.etaLev.size())) {
        maxRegionalHyperLev = amrex::max(maxRegionalHyperLev, cfg.etaLev[iLev]);
      }
    }
    const Real globalHyperLev =
        (iLev < static_cast<int>(etaHyperLev.size())) ? etaHyperLev[iLev] : 0.0;
    const Real maxHyper = globalHyperLev + maxRegionalHyperLev;
    if (maxEta <= 0 && maxHyper <= 0 && !useHallTerm && !useHybridUpwind)
      continue;

    const auto dx = Geom(iLev).CellSizeArray();
    const Box& domBox = Geom(iLev).Domain();
    // Staggered curl-curl has multilinear Fourier symbol
    // sum_d (4/dx_d^2) alpha_d prod_{d'!=d}(1 - alpha_{d'}) on alpha_d in
    // [0,1], attaining its global maximum sMax = max_d (4/dx_d^2) on
    // axis-aligned Nyquist modes (for any active dimension with nCellDim >= 2).
    // The cell-centered Laplacian has maximum symbol lMax = sum_d (4/dx_d^2).
    Real sMax = 0.0, lMax = 0.0;
    for (int iDim = 0; iDim < nDim; ++iDim) {
      const int nCellDim = domBox.length(iDim);
      if (nCellDim < 2)
        continue;
      const Real invDx2 = 1.0 / (dx[iDim] * dx[iDim]);
      lMax += 4.0 * invDx2;
      sMax = amrex::max(sMax, 4.0 * invDx2);
    }

    Real bMaxLev = 0.0;
    Real rhoMinLev = 0.0;
    Real rhoEff = 0.0;
    Real bOverRhoMax = 0.0;
    Real bOverSqrtRhoMax = 0.0;
    Real uMaxLev = 0.0;
    bool hasPlasma = false;
    if (useHallTerm || useHybridUpwind) {
      const MultiFab& bForCfl = total_center_B(centerB[iLev], iLev);
      auto const& bArrays = bForCfl.const_arrays();
      auto const& cStatArrays = cell_status(iLev).const_arrays();
      bMaxLev = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, bForCfl,
                          IntVect(0),
                          [=] AMREX_GPU_DEVICE(int nb, int i, int j, int k)
                              noexcept -> GpuTuple<Real> {
                                const int st = cStatArrays[nb](i, j, k);
                                if ((hasBody && bit::is_body(st)) ||
                                    (hasFiner && bit::is_refined(st)))
                                  return 0.0;
                                Array4<Real const> const& arr = bArrays[nb];
                                Real b2 = arr(i, j, k, 0) * arr(i, j, k, 0) +
                                          arr(i, j, k, 1) * arr(i, j, k, 1) +
                                          arr(i, j, k, 2) * arr(i, j, k, 2);
                                return sqrt(b2);
                              });
      ParallelDescriptor::ReduceRealMax(bMaxLev);

      // assemble_ohm_E disables Hall and wave-speed upwind terms at nodes with
      // rho <= 0 (vacuum, particle-free tests, or inside #BODY), interpolates
      // cell-centered B to nodeBstage using the full 8-cell stencil (plus BC
      // and conducting-body projection), and evaluates Hall/upwind speeds from
      // stage-interpolated densities between nodePlasmaPrev and nodePlasma.
      average_center_to_node(bForCfl, nodeBstage[iLev]);
      nodeBstage[iLev].FillBoundary(Geom(iLev).periodicity());
      if (iLev == 0) {
        apply_field_bc(node_status(iLev), nodeBstage[iLev], 0, 3,
                       &Pic::get_node_B, iLev, true, 0.0);
      }
      if (is_body_conducting()) {
        project_body_B(nodeBstage[iLev], iLev);
        nodeBstage[iLev].FillBoundary(Geom(iLev).periodicity());
      }

      const Real rhoFloor = amrex::max(rhoMinOhm, static_cast<Real>(1e-10));
      auto const& nbArrays = nodeBstage[iLev].const_arrays();
      auto const& rhoArrays = nodePlasma[nSpecies][iLev].const_arrays();
      const bool hasPrevMoments =
          !nodePlasmaPrev.empty() && !nodePlasmaPrev[nSpecies][iLev].empty();
      auto const& rhoPrevArrays =
          hasPrevMoments ? nodePlasmaPrev[nSpecies][iLev].const_arrays()
                         : rhoArrays;
      const bool isPC = usePredictorCorrector;
      const int nSubBound =
          amrex::max(nBSubcycleMax, amrex::max(nBSubcycle, 1));
      const Real minSubW = 0.5 / static_cast<Real>(nSubBound);
      auto const& nStatArrays = node_status(iLev).const_arrays();
      const int iRhoComp = iRho_;
      auto redTuple = ParReduce(
          TypeList<ReduceOpMin, ReduceOpMax, ReduceOpMax>{},
          TypeList<Real, Real, Real>{}, nodePlasma[nSpecies][iLev], IntVect(0),
          [=] AMREX_GPU_DEVICE(int nb, int i, int j,
                               int k) noexcept -> GpuTuple<Real, Real, Real> {
            if (!is_active_cfl_node(cStatArrays[nb], nStatArrays[nb], i, j, k,
                                    hasBody, hasFiner))
              return { std::numeric_limits<Real>::max(), 0.0, 0.0 };
            const Real rCur = rhoArrays[nb](i, j, k, iRhoComp);
            const Real rPrev =
                hasPrevMoments ? rhoPrevArrays[nb](i, j, k, iRhoComp) : rCur;
            const Real r = min_positive_stage_rho(rCur, rPrev, hasPrevMoments,
                                                  isPC, minSubW);
            if (r >= std::numeric_limits<Real>::max())
              return { std::numeric_limits<Real>::max(), 0.0, 0.0 };
            Array4<Real const> const& nbArr = nbArrays[nb];
            const Real bNode =
                std::sqrt(nbArr(i, j, k, ix_) * nbArr(i, j, k, ix_) +
                          nbArr(i, j, k, iy_) * nbArr(i, j, k, iy_) +
                          nbArr(i, j, k, iz_) * nbArr(i, j, k, iz_));
            const Real rEff = amrex::max(r, rhoFloor);
            return { r, bNode / rEff, bNode / std::sqrt(rEff) };
          });
      Real rhoPosMin = amrex::get<0>(redTuple);
      bOverRhoMax = amrex::get<1>(redTuple);
      bOverSqrtRhoMax = amrex::get<2>(redTuple);
      ParallelDescriptor::ReduceRealMin(rhoPosMin);
      Real maxRatios[2] = { bOverRhoMax, bOverSqrtRhoMax };
      ParallelDescriptor::ReduceRealMax(maxRatios, 2);
      bOverRhoMax = maxRatios[0];
      bOverSqrtRhoMax = maxRatios[1];

      hasPlasma = (rhoPosMin < std::numeric_limits<Real>::max());
      rhoMinLev = hasPlasma ? rhoPosMin : 0.0;
      rhoEff = hasPlasma ? amrex::max(rhoMinLev, rhoFloor) : 0.0;

      // Bound stage electron velocity (bulk velocity when Hall is disabled) for
      // the upwind dissipation rate estimate across active solver stages.
      if (useHybridUpwind && upwindFixedSpeed <= 0.0 && hasPlasma) {
        if (useHallTerm) {
          curl_center_to_node(centerB[iLev], nodeJ[iLev],
                              Geom(iLev).InvCellSize());
          nodeJ[iLev].FillBoundary(Geom(iLev).periodicity());
          auto const& jArrays = nodeJ[iLev].const_arrays();
          const Real invFourPI = 1.0 / fourPI;
          uMaxLev = ParReduce(
              TypeList<ReduceOpMax>{}, TypeList<Real>{},
              nodePlasma[nSpecies][iLev], IntVect(0),
              [=] AMREX_GPU_DEVICE(int nb, int i, int j,
                                   int k) noexcept -> GpuTuple<Real> {
                if (!is_active_cfl_node(cStatArrays[nb], nStatArrays[nb], i, j,
                                        k, hasBody, hasFiner))
                  return 0.0;
                auto eval_ue2 = [&](const Real rStage, const Real wPrevM)
                                    noexcept -> Real {
                  if (rStage <= 0.0)
                    return 0.0;
                  const Real invRho = 1.0 / amrex::max(rStage, rhoMinOhm);
                  const Real wCurM = 1.0 - wPrevM;
                  Real u2 = 0.0;
                  for (int d = 0; d < 3; ++d) {
                    const Real mStage =
                        wCurM * rhoArrays[nb](i, j, k, iMx_ + d) +
                        wPrevM * rhoPrevArrays[nb](i, j, k, iMx_ + d);
                    const Real ue_d =
                        (mStage - jArrays[nb](i, j, k, d) * invFourPI) * invRho;
                    u2 += ue_d * ue_d;
                  }
                  return u2;
                };
                const Real rCur = rhoArrays[nb](i, j, k, iRhoComp);
                Real u2Max = eval_ue2(rCur, 0.0);
                if (hasPrevMoments) {
                  const Real rPrev = rhoPrevArrays[nb](i, j, k, iRhoComp);
                  const Real rMid = 0.5 * (rPrev + rCur);
                  u2Max = amrex::max(u2Max, eval_ue2(rMid, isPC ? 0.5 : 0.0));
                  if (isPC) {
                    u2Max = amrex::max(u2Max, eval_ue2(rPrev, 1.0));
                  }
                }
                return std::sqrt(u2Max);
              });
        } else {
          uMaxLev = ParReduce(
              TypeList<ReduceOpMax>{}, TypeList<Real>{},
              nodePlasma[nSpecies][iLev], IntVect(0),
              [=] AMREX_GPU_DEVICE(int nb, int i, int j,
                                   int k) noexcept -> GpuTuple<Real> {
                if (!is_active_cfl_node(cStatArrays[nb], nStatArrays[nb], i, j,
                                        k, hasBody, hasFiner))
                  return 0.0;
                auto eval_ui2 = [&](const Real rStage, const Real wPrevM)
                                    noexcept -> Real {
                  if (rStage <= 0.0)
                    return 0.0;
                  const Real invRho = 1.0 / amrex::max(rStage, rhoMinOhm);
                  const Real wCurM = 1.0 - wPrevM;
                  Real u2 = 0.0;
                  for (int d = 0; d < 3; ++d) {
                    const Real mStage =
                        wCurM * rhoArrays[nb](i, j, k, iMx_ + d) +
                        wPrevM * rhoPrevArrays[nb](i, j, k, iMx_ + d);
                    const Real ui_d = mStage * invRho;
                    u2 += ui_d * ui_d;
                  }
                  return u2;
                };
                const Real rCur = rhoArrays[nb](i, j, k, iRhoComp);
                Real u2Max = eval_ui2(rCur, 0.0);
                if (hasPrevMoments) {
                  const Real rPrev = rhoPrevArrays[nb](i, j, k, iRhoComp);
                  const Real rMid = 0.5 * (rPrev + rCur);
                  u2Max = amrex::max(u2Max, eval_ui2(rMid, isPC ? 0.5 : 0.0));
                  if (isPC) {
                    u2Max = amrex::max(u2Max, eval_ui2(rPrev, 1.0));
                  }
                }
                return std::sqrt(u2Max);
              });
        }
        ParallelDescriptor::ReduceRealMax(uMaxLev);
      }
    }

    const Real rateEtaDt = (maxEta / fourPI) * sMax * dt;
    const Real rateHyperDt = (maxHyper / fourPI) * sMax * lMax * dt;
    const Real rateUpwindDt =
        useHybridUpwind ? upwind_max_rate(iLev, uMaxLev, bOverSqrtRhoMax,
                                          bOverRhoMax, hasPlasma) *
                              dt
                        : 0.0;
    const Real rateWhistlerDt =
        (useHallTerm && hasPlasma) ? (sMax * bOverRhoMax / fourPI) * dt : 0.0;
    // Bound the summed dissipative rate on the real axis (resistive +
    // hyper-resistive + upwind dissipation).
    const Real rateDissipativeDt = rateEtaDt + rateHyperDt + rateUpwindDt;

    // When auto-subcycling is enabled, immediately adapt nBSubcycle before the
    // Faraday sub-step loop if the freshly deposited moments or updated fields
    // would violate any CFL limit with the current nBSubcycle.
    if (isAutoSubcycle) {
      int nSubNeed = nBSubcycle;
      if (rateDissipativeDt / nBSubcycle > cflLimitReal) {
        nSubNeed =
            amrex::max(nSubNeed, safe_ceil_int(rateDissipativeDt /
                                               (cflLimitReal * fieldSafety)));
      }
      if (rateWhistlerDt / nBSubcycle > cflLimitWhistler) {
        nSubNeed = amrex::max(
            nSubNeed,
            safe_ceil_int(rateWhistlerDt / (cflLimitWhistler * fieldSafety)));
      }
      const int nSubClamped =
          std::clamp(nSubNeed, nBSubcycleMin, nBSubcycleMax);
      if (nSubClamped > nBSubcycle) {
        nBSubcycle = nSubClamped;
        nSubcycle = nSubClamped;
        subcycleHoldCount = 0;
        subDt = dt / nBSubcycle;
      }
    }

    const Real cflDissipative = rateDissipativeDt / nBSubcycle;
    const Real cflWhistler = rateWhistlerDt / nBSubcycle;
    if (cflDissipative > cflLimitReal) {
      const Real cflEta = rateEtaDt / nBSubcycle;
      const Real cflHyper = rateHyperDt / nBSubcycle;
      const Real cflUpwind = rateUpwindDt / nBSubcycle;
      amrex::Print() << "  [CFL warning] combined dissipation: (eta + hyper + "
                        "upwind)*dt_sub = "
                     << cflDissipative << " [eta=" << cflEta
                     << ", hyper=" << cflHyper << ", upwind=" << cflUpwind
                     << "] (> " << cflLimitReal
                     << ", explicit diffusion may be unstable)\n";
    }

    if (useHallTerm && hasPlasma) {
      if (cflWhistler > cflLimitWhistler) {
        const int nSubReq = safe_ceil_int(rateWhistlerDt / cflLimitWhistler);
        amrex::Print() << "  [CFL warning] whistler wave: "
                          "(kmax^2*Bmax)/(4pi*rho_eff)*dt_sub = "
                       << cflWhistler << " (> " << cflLimitWhistler
                       << ", whistler wave is linearly unstable in "
                       << fieldIntegrator << "! rho_min=" << rhoMinLev
                       << ", rho_eff=" << rhoEff << ", B_max=" << bMaxLev
                       << ", (B/rho_eff)_max=" << bOverRhoMax
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
  const bool ambiInStages = (electronTemperature > 0) &&
                            (!useElectronPressureEq || ambipolarInStages);

  const Real invSubcycle = 1.0 / static_cast<Real>(nBSubcycle);

  // Refresh cell-centered B ghosts at tStepStart before the first RK/SSPRK
  // stage so upwind and hyper-resistive stencils read stage-consistent ghosts.
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    apply_centerB_BC(iLev, tStepStart);
  }

  // Each stage first evaluates curl(E) on ALL levels (faraday_stage_rhs, which
  // also synchronizes E across the coarse-fine interface), then updates every
  // level, then refreshes the ghost cells. nG = faraday_ngrow(iLev) is 1 on
  // refined levels: the first ghost layer is advanced with the same curl(E)
  // as the valid cells.
  for (int subStep = 0; subStep < nBSubcycle; ++subStep) {
    // Moment time-interpolation weights hstep for RK stages within the
    // sub-step. In standard leapfrog (!usePredictorCorrector), the freshly
    // deposited moments J^{n+1/2} are already centered at the midpoint h = 0.5
    // of [t^n, t^{n+1}]; holding hstep = 0.5 across the Faraday sub-steps
    // preserves 2nd-order midpoint accuracy while eliminating the secant-ramp
    // whistler instability at omega_w * dt ~ pi.
    const Real g = static_cast<Real>(subStep) * invSubcycle;
    const Real hstepStart = usePredictorCorrector ? g : 0.5;
    const Real hstepHalf =
        usePredictorCorrector ? (g + 0.5 * invSubcycle) : 0.5;
    const Real hstepEnd = usePredictorCorrector ? (g + invSubcycle) : 0.5;
    const Real tSubStart = tStepStart + static_cast<Real>(subStep) * subDt;
    const Real tSubHalf = tSubStart + 0.5 * subDt;
    const Real tSubEnd = tSubStart + subDt;

    if (useRK4) {
      // Classical RK4 on dB/dt = -curl(E_Ohm(B)) across all levels.
      const Real dtSixth = -subDt / 6.0;
      const Real dtThird = -subDt / 3.0;

      // Stages 1-3: trial state B_s = B_n - c_s * subDt * k_{s-1}.
      const Real cStage[3] = { 0.5, 0.5, 1.0 };
      const Real hStage[4] = { hstepStart, hstepHalf, hstepHalf, hstepEnd };
      const Real tStageRhs[4] = { tSubStart, tSubHalf, tSubHalf, tSubEnd };
      const Real tStageOut[3] = { tSubHalf, tSubHalf, tSubEnd };
      for (int s = 0; s < 3; ++s) {
        faraday_stage_rhs(s == 0 ? centerB : centerBstage, s, hStage[s],
                          ambiInStages, tStageRhs[s]);
        for (int iLev = 0; iLev < n_lev(); ++iLev) {
          const int nG = faraday_ngrow(iLev);
          MultiFab::Copy(centerBstage[iLev], centerB[iLev], 0, 0, nDim3,
                         centerBstage[iLev].nGrow());
          MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0,
                            -cStage[s] * subDt, kStage[iLev][s], 0, 0, nDim3,
                            nG);
        }
        for (int iLev = 0; iLev < n_lev(); ++iLev) {
          apply_centerB_BC(iLev, centerBstage[iLev], tStageOut[s],
                           cStage[s] * subDt);
        }
      }

      // Stage 4 and the final RK4 combination.
      faraday_stage_rhs(centerBstage, 3, hStage[3], ambiInStages, tStageRhs[3]);
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
        apply_centerB_BC(iLev, tSubEnd, subDt);
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
      faraday_stage_rhs(centerB, 0, hstepStart, ambiInStages, tSubStart);
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        const int nG = faraday_ngrow(iLev);
        MultiFab::Copy(centerBstage[iLev], centerB[iLev], 0, 0, nDim3,
                       centerBstage[iLev].nGrow());
        MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0, -subDt,
                          kStage[iLev][0], 0, 0, nDim3, nG);
      }
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        apply_centerB_BC(iLev, centerBstage[iLev], tSubEnd, subDt);
      }

      // Stage 2: B2 = (3/4)*B_n + (1/4)*(B1 - subDt * curl(E(B1)))
      faraday_stage_rhs(centerBstage, 1, hstepEnd, ambiInStages, tSubEnd);
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        const int nG = faraday_ngrow(iLev);
        MultiFab::Saxpy(centerBstage[iLev], -subDt, kStage[iLev][1], 0, 0,
                        nDim3, nG);
      }
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        apply_centerB_BC(iLev, centerBstage[iLev], tSubEnd, subDt);
        MultiFab::LinComb(centerBstage[iLev], 0.25, centerBstage[iLev], 0, 0.75,
                          centerBstart[iLev], 0, 0, nDim3,
                          centerBstage[iLev].nGrow());
      }
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        apply_centerB_BC(iLev, centerBstage[iLev], tSubHalf, -2.0);
      }

      // Stage 3: B^{n+1} = (1/3)*B_n + (2/3)*(B2 - subDt * curl(E(B2)))
      faraday_stage_rhs(centerBstage, 2, hstepHalf, ambiInStages, tSubHalf);
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        const int nG = faraday_ngrow(iLev);
        MultiFab::Saxpy(centerBstage[iLev], -subDt, kStage[iLev][2], 0, 0,
                        nDim3, nG);
      }
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        apply_centerB_BC(iLev, centerBstage[iLev], tSubEnd, subDt);
        MultiFab::LinComb(centerB[iLev], 2.0 / 3.0, centerBstage[iLev], 0,
                          1.0 / 3.0, centerBstart[iLev], 0, 0, nDim3,
                          centerB[iLev].nGrow());
      }
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        apply_centerB_BC(iLev, tSubEnd, -2.0);
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

  refresh_intrinsic_fields(tStepEnd);

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
                       &Pic::get_node_B, iLev, true, 0.0);
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
                     iLev, true, 0.0);
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
                       &Pic::get_node_E, iLev, false, 0.0);
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
Real Pic::upwind_max_rate(int iLev, Real uMax, Real bOverSqrtRhoMax,
                          Real bOverRhoMax, bool hasPlasma) const {
  if (!useHybridUpwind || (!hasPlasma && upwindFixedSpeed <= 0.0))
    return 0.0;

  Real cFast = 0.0, cWhistlerDx = 0.0;
  if (hasPlasma && upwindFixedSpeed <= 0.0) {
    const Real teUp = (useElectronPressureEq && peStateInitialized_)
                          ? max_evolved_Te(iLev)
                          : electronTemperature;
    Real cs2 = 0.0;
    if (teUp > 0.0)
      cs2 = ((electronGamma > 0.0) ? electronGamma : 1.0) * teUp;
    cFast = std::sqrt(bOverSqrtRhoMax * bOverSqrtRhoMax / fourPI + cs2);
    if (useHallTerm)
      cWhistlerDx = upwindCoefWhistler * 0.25 * bOverRhoMax;
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

  const Real cflLimitReal = useRK4 ? 2.785 : 2.513;
  const Real cflLimitWhistler = useRK4 ? 2.785 : 1.732;
  const Real fieldSafety = 0.85;
  auto safe_ceil_int = [](Real val) -> int {
    if (!std::isfinite(val) ||
        val >= static_cast<Real>(std::numeric_limits<int>::max() / 2)) {
      return std::numeric_limits<int>::max() / 2;
    }
    if (val <= 1.0)
      return 1;
    return static_cast<int>(std::ceil(val));
  };

  Real maxRegionalEta = 0.0;
  for (const auto& cfg : regionalResistivityConfigs) {
    maxRegionalEta = amrex::max(maxRegionalEta, cfg.etaCode);
  }
  const Real maxEta = etaResistivity + maxRegionalEta;

  Real dtMacroMin = std::numeric_limits<Real>::max();
  Real dtSubMin = std::numeric_limits<Real>::max();
  int nSubHyperGridMin = 1;
  std::string localMacroReason = "ion kinetics";
  std::string localSubReason = "field advance";

  struct LevSubReq {
    Real omegaDiss = 0.0;
    Real cflHyperMacro = 0.0;
    Real dtWhistler = std::numeric_limits<Real>::max();
    Real dtWave = std::numeric_limits<Real>::max();
  };
  Vector<LevSubReq> levReqs(n_lev());
  const bool hasBody = use_body();
  if (doReport) {
    if (fixedUMax >= 0) {
      maxUthSpecies.clear();
    } else {
      maxUthSpecies.assign(n_lev() * nSpecies, 0.0);
    }
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    LevSubReq& req = levReqs[iLev];
    const bool hasFiner = (iLev < finestLevel());
    const auto dx = Geom(iLev).CellSizeArray();
    Real dxMin = dx[0];
    for (int d = 1; d < nDim; ++d) {
      dxMin = amrex::min(dxMin, dx[d]);
    }

    auto const& cStatArrays = cell_status(iLev).const_arrays();
    auto const& nStatArrays = node_status(iLev).const_arrays();

    // 1. Kinetic ion advection limit (including bulk flow + thermal spread)
    Real uMaxLev = 0.0;
    if (fixedUMax >= 0) {
      uMaxLev = fixedUMax;
    } else {
      Vector<Real> uMaxAndUth(doReport ? 2 * nSpecies : nSpecies, 0.0);
      for (int i = 0; i < nSpecies; ++i) {
        amrex::MultiFab& momMF = nodePlasma[i][iLev];
        auto const& ma = momMF.const_arrays();
        auto [uSpMax, uthSpMax] = ParReduce(
            TypeList<ReduceOpMax, ReduceOpMax>{}, TypeList<Real, Real>{}, momMF,
            IntVect(0),
            [=] AMREX_GPU_DEVICE(int nb, int i, int j, int k)
                noexcept -> GpuTuple<Real, Real> {
                  Array4<Real const> const& arr = ma[nb];
                  Real rho = arr(i, j, k, iRho_);
                  if (rho <= 0.0 ||
                      !is_active_cfl_node(cStatArrays[nb], nStatArrays[nb], i,
                                          j, k, hasBody, hasFiner))
                    return { 0.0, 0.0 };
                  Real invRho = 1.0 / rho;
                  Real uBulk2 = (arr(i, j, k, iUx_) * arr(i, j, k, iUx_) +
                                 arr(i, j, k, iUy_) * arr(i, j, k, iUy_) +
                                 arr(i, j, k, iUz_) * arr(i, j, k, iUz_)) *
                                (invRho * invRho);
                  Real p = (arr(i, j, k, iPxx_) + arr(i, j, k, iPyy_) +
                            arr(i, j, k, iPzz_)) *
                           (1.0 / 3.0);
                  Real vth = (p > 0.0) ? sqrt(p * invRho) : 0.0;
                  return { sqrt(uBulk2) + 3.0 * vth, vth };
                });
        uMaxAndUth[i] = uSpMax;
        if (doReport) {
          uMaxAndUth[nSpecies + i] = uthSpMax;
        }
      }
      ParallelDescriptor::ReduceRealMax(uMaxAndUth.data(), uMaxAndUth.size());
      for (int i = 0; i < nSpecies; ++i) {
        if (uMaxAndUth[i] > uMaxLev) {
          uMaxLev = uMaxAndUth[i];
        }
        if (doReport) {
          const Real uSp = uMaxAndUth[nSpecies + i];
          maxUthSpecies[iLev * nSpecies + i] = uSp;
          if (!doCompact) {
            amrex::Print() << printPrefix << std::setprecision(5) << "lev "
                           << iLev << " Species " << i << ": max(uth) = " << uSp
                           << std::endl;
          }
        }
      }
    }

    // 2. Maximum magnetic field magnitude on active plasma cells and their
    // corner nodes on this level (including intrinsic B0 and nodal stencils
    // that reach into #BODY or masked neighbor cells).
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
    const MultiFab& bForCfl = total_center_B(centerB[iLev], iLev);
    average_center_to_node(bForCfl, nodeBstage[iLev]);
    nodeBstage[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(node_status(iLev), nodeBstage[iLev], 0, 3,
                     &Pic::get_node_B, iLev, true, 0.0);
    }
    if (is_body_conducting()) {
      project_body_B(nodeBstage[iLev], iLev);
      nodeBstage[iLev].FillBoundary(Geom(iLev).periodicity());
    }

    auto const& bArrays = bForCfl.const_arrays();
    auto const& nbArrays = nodeBstage[iLev].const_arrays();
    const bool hasMoverNodeB = !nodeB[iLev].empty();
    const bool hasMoverNodeB0 = use_intrinsic_B() && !nodeB0[iLev].empty();
    auto const& moverBArrays =
        hasMoverNodeB ? nodeB[iLev].const_arrays() : nbArrays;
    auto const& moverB0Arrays =
        hasMoverNodeB0 ? nodeB0[iLev].const_arrays() : nbArrays;
    Real bMaxLev = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, bForCfl,
                             IntVect(0),
                             [=] AMREX_GPU_DEVICE(int nb, int i, int j, int k)
                                 noexcept -> GpuTuple<Real> {
                                   const int st = cStatArrays[nb](i, j, k);
                                   if ((hasBody && bit::is_body(st)) ||
                                       (hasFiner && bit::is_refined(st)))
                                     return 0.0;
                                   Array4<Real const> const& arr = bArrays[nb];
                                   Real b2 = arr(i, j, k, 0) * arr(i, j, k, 0) +
                                             arr(i, j, k, 1) * arr(i, j, k, 1) +
                                             arr(i, j, k, 2) * arr(i, j, k, 2);
                                   return sqrt(b2);
                                 });
    Real bGyroNodeMaxLev = ParReduce(
        TypeList<ReduceOpMax>{}, TypeList<Real>{}, nodeBstage[iLev], IntVect(0),
        [=] AMREX_GPU_DEVICE(int nb, int i, int j,
                             int k) noexcept -> GpuTuple<Real> {
          if (!touches_active_cell(cStatArrays[nb], i, j, k, hasBody, hasFiner))
            return 0.0;
          Array4<Real const> const& nbArr = nbArrays[nb];
          Real b2 = nbArr(i, j, k, ix_) * nbArr(i, j, k, ix_) +
                    nbArr(i, j, k, iy_) * nbArr(i, j, k, iy_) +
                    nbArr(i, j, k, iz_) * nbArr(i, j, k, iz_);
          if (hasMoverNodeB) {
            const Real bx =
                moverBArrays[nb](i, j, k, ix_) +
                (hasMoverNodeB0 ? moverB0Arrays[nb](i, j, k, ix_) : 0.0);
            const Real by =
                moverBArrays[nb](i, j, k, iy_) +
                (hasMoverNodeB0 ? moverB0Arrays[nb](i, j, k, iy_) : 0.0);
            const Real bz =
                moverBArrays[nb](i, j, k, iz_) +
                (hasMoverNodeB0 ? moverB0Arrays[nb](i, j, k, iz_) : 0.0);
            b2 = amrex::max(b2, bx * bx + by * by + bz * bz);
          }
          return std::sqrt(b2);
        });
    Real bMaxPair[2] = { bMaxLev, bGyroNodeMaxLev };
    ParallelDescriptor::ReduceRealMax(bMaxPair, 2);
    bMaxLev = bMaxPair[0];
    const Real bGyroMaxLev = amrex::max(bMaxLev, bMaxPair[1]);

    // Kinetic advection limit
    Real dtPart =
        (uMaxLev > 0.0) ? (dxMin / uMaxLev) : std::numeric_limits<Real>::max();

    // Ion gyrofrequency limit (Omega_ci * dt <= thetaGyro)
    Real dtGyro = std::numeric_limits<Real>::max();
    if (bGyroMaxLev > 0.0) {
      Real maxQom = 0.0;
      for (int i = 0; i < nSpecies; ++i) {
        Real qom = std::abs(parts[i]->get_charge() / parts[i]->get_mass());
        maxQom = amrex::max(maxQom, qom);
      }
      if (maxQom > 0.0) {
        const Real thetaGyro = 0.35; // radians (~18 steps per gyroperiod)
        dtGyro = thetaGyro / (maxQom * bGyroMaxLev);
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
      sMax = amrex::max(sMax, 4.0 * invDx2);
    }

    if (sMax > 0.0) {
      const Real rhoFloor = amrex::max(rhoMinOhm, static_cast<Real>(1e-10));
      auto const& rhoArrays = nodePlasma[nSpecies][iLev].const_arrays();
      const bool hasPrevMoments = !isFirstHybridStep &&
                                  !nodePlasmaPrev.empty() &&
                                  !nodePlasmaPrev[nSpecies][iLev].empty();
      auto const& rhoPrevArrays =
          hasPrevMoments ? nodePlasmaPrev[nSpecies][iLev].const_arrays()
                         : rhoArrays;
      const bool isPC = usePredictorCorrector;
      const int nSubBound =
          amrex::max(nBSubcycleMax, amrex::max(nBSubcycle, 1));
      const Real minSubW = 0.5 / static_cast<Real>(nSubBound);
      const int iRhoComp = iRho_;
      auto redTuple = ParReduce(
          TypeList<ReduceOpMin, ReduceOpMax, ReduceOpMax>{},
          TypeList<Real, Real, Real>{}, nodePlasma[nSpecies][iLev], IntVect(0),
          [=] AMREX_GPU_DEVICE(int nb, int i, int j,
                               int k) noexcept -> GpuTuple<Real, Real, Real> {
            if (!is_active_cfl_node(cStatArrays[nb], nStatArrays[nb], i, j, k,
                                    hasBody, hasFiner))
              return { std::numeric_limits<Real>::max(), 0.0, 0.0 };
            const Real rCur = rhoArrays[nb](i, j, k, iRhoComp);
            const Real rPrev =
                hasPrevMoments ? rhoPrevArrays[nb](i, j, k, iRhoComp) : rCur;
            const Real r = min_positive_stage_rho(rCur, rPrev, hasPrevMoments,
                                                  isPC, minSubW);
            if (r >= std::numeric_limits<Real>::max())
              return { std::numeric_limits<Real>::max(), 0.0, 0.0 };
            Array4<Real const> const& nbArr = nbArrays[nb];
            const Real bNode =
                std::sqrt(nbArr(i, j, k, ix_) * nbArr(i, j, k, ix_) +
                          nbArr(i, j, k, iy_) * nbArr(i, j, k, iy_) +
                          nbArr(i, j, k, iz_) * nbArr(i, j, k, iz_));
            const Real rEff = amrex::max(r, rhoFloor);
            return { r, bNode / rEff, bNode / std::sqrt(rEff) };
          });
      Real rhoPosMin = amrex::get<0>(redTuple);
      Real bOverRhoMax = amrex::get<1>(redTuple);
      Real bOverSqrtRhoMax = amrex::get<2>(redTuple);
      ParallelDescriptor::ReduceRealMin(rhoPosMin);
      Real maxRatios[2] = { bOverRhoMax, bOverSqrtRhoMax };
      ParallelDescriptor::ReduceRealMax(maxRatios, 2);
      bOverRhoMax = maxRatios[0];
      bOverSqrtRhoMax = maxRatios[1];

      const bool hasPlasma = (rhoPosMin < std::numeric_limits<Real>::max());

      if (useHallTerm && hasPlasma && bOverRhoMax > 0.0) {
        const Real omegaWhistler = (sMax * bOverRhoMax) / fourPI;
        if (omegaWhistler > 0.0) {
          req.dtWhistler = (cflLimitWhistler * fieldSafety) / omegaWhistler;
          if (req.dtWhistler < dtSubMin) {
            dtSubMin = req.dtWhistler;
            localSubReason = "whistler wave";
          }
        }
      }

      // Fast magnetosonic / Alfvén wave limit (governs when useHallTerm is
      // false or on coarse cells)
      const Real teWave = (useElectronPressureEq && peStateInitialized_)
                              ? max_evolved_Te(iLev)
                              : electronTemperature;
      if (hasPlasma && (bOverSqrtRhoMax > 0.0 || teWave > 0.0)) {
        Real vA = bOverSqrtRhoMax / std::sqrt(fourPI);
        Real cs = 0.0;
        if (teWave > 0.0) {
          Real gammaE = (electronGamma > 0.0) ? electronGamma : 1.0;
          cs = std::sqrt(gammaE * teWave);
        }
        Real vFast = std::sqrt(vA * vA + cs * cs);
        if (vFast > 0.0) {
          req.dtWave = (cflLimitReal * fieldSafety) * dxMin / vFast;
          if (req.dtWave < dtSubMin) {
            dtSubMin = req.dtWave;
            localSubReason =
                useHallTerm ? "Alfvén / fast wave" : "ideal-MHD wave";
          }
        }
      }

      // Combined real-axis dissipation (resistivity + hyper-resistivity +
      // upwinding)
      if (maxEta > 0.0) {
        req.omegaDiss += (maxEta / fourPI) * sMax;
      }

      if (lMax > 0.0) {
        // Physical ("si" mode) hyper-resistivity sets a physical time scale
        // dtHyper; "grid" mode (eta_h = 4*pi*C_h*dx_min^4/dt) has dimensionless
        // per-macro-step rate C_h*dx_min^4*sMax*lMax that directly constrains
        // nSubcycle. Regional and global hyper-resistivity are added together
        // in assemble_ohm_E, so bound their combined contribution additively.
        Real maxHyperSI = (etaHyperMode != "grid" &&
                           iLev < static_cast<int>(etaHyperLev.size()))
                              ? etaHyperLev[iLev]
                              : 0.0;
        Real maxRegionalHyperSI = 0.0;
        Real maxRegionalHyperCh = 0.0;
        for (const auto& cfg : regionalHyperResistivityConfigs) {
          if (cfg.mode == "grid") {
            maxRegionalHyperCh = amrex::max(maxRegionalHyperCh, cfg.ch);
          } else if (iLev < static_cast<int>(cfg.etaLev.size())) {
            maxRegionalHyperSI =
                amrex::max(maxRegionalHyperSI, cfg.etaLev[iLev]);
          }
        }
        maxHyperSI += maxRegionalHyperSI;
        const Real globalGridCh =
            (etaHyperMode == "grid" && etaHyperCh > 0.0) ? etaHyperCh : 0.0;
        const Real maxHyperChGrid = globalGridCh + maxRegionalHyperCh;

        if (maxHyperSI > 0.0) {
          req.omegaDiss += (maxHyperSI / fourPI) * sMax * lMax;
        }
        if (maxHyperChGrid > 0.0) {
          req.cflHyperMacro = maxHyperChGrid * std::pow(dxMin, 4) * sMax * lMax;
        }
      }

      if (useHybridUpwind) {
        req.omegaDiss += upwind_max_rate(iLev, uMaxLev, bOverSqrtRhoMax,
                                         bOverRhoMax, hasPlasma);
      }

      if (req.omegaDiss > 0.0) {
        const Real dtDiss = (cflLimitReal * fieldSafety) / req.omegaDiss;
        if (dtDiss < dtSubMin) {
          dtSubMin = dtDiss;
          localSubReason = "dissipative diffusion";
        }
      }
      if (req.cflHyperMacro > 0.0) {
        const int nReqGrid =
            safe_ceil_int(req.cflHyperMacro / (cflLimitReal * fieldSafety));
        if (nReqGrid > nSubHyperGridMin) {
          nSubHyperGridMin = nReqGrid;
          localSubReason = "hyper-resistivity";
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

    if (isAutoSubcycle) {
      int nSubReq = 1;
      std::string activeSubReason = "field advance";
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        const auto& req = levReqs[iLev];
        const Real rateDissDt = req.omegaDiss * dtFixed + req.cflHyperMacro;
        const int nReqDiss =
            safe_ceil_int(rateDissDt / (cflLimitReal * fieldSafety));
        const int nReqWhist =
            (req.dtWhistler < std::numeric_limits<Real>::max())
                ? safe_ceil_int(dtFixed / req.dtWhistler)
                : 1;
        const int nReqWave = (req.dtWave < std::numeric_limits<Real>::max())
                                 ? safe_ceil_int(dtFixed / req.dtWave)
                                 : 1;
        if (nReqDiss > nSubReq) {
          nSubReq = nReqDiss;
          activeSubReason = "dissipative diffusion";
        }
        if (nReqWhist > nSubReq) {
          nSubReq = nReqWhist;
          activeSubReason = "whistler wave";
        }
        if (nReqWave > nSubReq) {
          nSubReq = nReqWave;
          activeSubReason =
              useHallTerm ? "Alfvén / fast wave" : "ideal-MHD wave";
        }
      }
      nSubReq = std::clamp(nSubReq, nBSubcycleMin, nBSubcycleMax);

      if (nSubReq > nBSubcycle) {
        // Immediate increase for numerical stability
        nSubNext = nSubReq;
        subcycleHoldCount = 0;
        localSubReason = activeSubReason;
      } else if (nSubReq < nBSubcycle) {
        // Hysteresis with 15% safety buffer before decreasing
        const int testSub = nBSubcycle - 1;
        bool canDecreaseSub = (testSub >= nBSubcycleMin);
        for (int iLev = 0; iLev < n_lev() && canDecreaseSub; ++iLev) {
          const auto& req = levReqs[iLev];
          const Real rateDissSub =
              (req.omegaDiss * dtFixed + req.cflHyperMacro) / testSub;
          if (rateDissSub > 0.85 * cflLimitReal * fieldSafety) {
            canDecreaseSub = false;
          }
          if (req.dtWhistler < std::numeric_limits<Real>::max() &&
              dtFixed / testSub > 0.85 * req.dtWhistler) {
            canDecreaseSub = false;
          }
          if (req.dtWave < std::numeric_limits<Real>::max() &&
              dtFixed / testSub > 0.85 * req.dtWave) {
            canDecreaseSub = false;
          }
        }
        if (canDecreaseSub) {
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
      dtLimitingReason = localSubReason + " (subcycling)";
    } else {
      nSubNext = nBSubcycle;
      dtLimitingReason = "fixed dt";
    }
  } else {
    // Adaptive macro dt mode
    if (!isAutoSubcycle) {
      Real dtFieldLimit = std::numeric_limits<Real>::max();
      std::string fieldReason = localMacroReason;
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        const auto& req = levReqs[iLev];
        if (req.omegaDiss > 0.0) {
          const Real maxDissRate = nBSubcycle * cflLimitReal * fieldSafety;
          const Real dtLevDiss =
              (maxDissRate > req.cflHyperMacro)
                  ? ((maxDissRate - req.cflHyperMacro) / req.omegaDiss)
                  : 0.0;
          if (dtLevDiss < dtFieldLimit) {
            dtFieldLimit = dtLevDiss;
            fieldReason = "dissipative diffusion";
          }
        }
        if (req.dtWhistler < std::numeric_limits<Real>::max()) {
          const Real dtLevWhist = nBSubcycle * req.dtWhistler;
          if (dtLevWhist < dtFieldLimit) {
            dtFieldLimit = dtLevWhist;
            fieldReason = "whistler wave";
          }
        }
        if (req.dtWave < std::numeric_limits<Real>::max()) {
          const Real dtLevWave = nBSubcycle * req.dtWave;
          if (dtLevWave < dtFieldLimit) {
            dtFieldLimit = dtLevWave;
            fieldReason = useHallTerm ? "Alfvén / fast wave" : "ideal-MHD wave";
          }
        }
      }

      Real dtTarget = std::min(dtMacroTarget, dtFieldLimit);
      if (dtTarget == dtFieldLimit &&
          dtFieldLimit < std::numeric_limits<Real>::max()) {
        dtLimitingReason = fieldReason + " (fixed subcycles)";
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
      Real dtOld = tc->get_dt();
      Real dtTarget = (dtOld > 0.0) ? std::min(dtMacroTarget,
                                               static_cast<Real>(1.10) * dtOld)
                                    : dtMacroTarget;
      int nIdeal = 1;
      std::string activeSubReason = localMacroReason;
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        const auto& req = levReqs[iLev];
        const Real rateDissDt = req.omegaDiss * dtTarget + req.cflHyperMacro;
        const int nReqDiss =
            safe_ceil_int(rateDissDt / (cflLimitReal * fieldSafety));
        const int nReqWhist =
            (req.dtWhistler < std::numeric_limits<Real>::max())
                ? safe_ceil_int(dtTarget / req.dtWhistler)
                : 1;
        const int nReqWave = (req.dtWave < std::numeric_limits<Real>::max())
                                 ? safe_ceil_int(dtTarget / req.dtWave)
                                 : 1;
        if (nReqDiss > nIdeal) {
          nIdeal = nReqDiss;
          activeSubReason = "dissipative diffusion";
        }
        if (nReqWhist > nIdeal) {
          nIdeal = nReqWhist;
          activeSubReason = "whistler wave";
        }
        if (nReqWave > nIdeal) {
          nIdeal = nReqWave;
          activeSubReason =
              useHallTerm ? "Alfvén / fast wave" : "ideal-MHD wave";
        }
      }

      if (nIdeal > nBSubcycleMax) {
        nSubNext = nBSubcycleMax;
        for (int iLev = 0; iLev < n_lev(); ++iLev) {
          const auto& req = levReqs[iLev];
          if (req.omegaDiss > 0.0) {
            const Real maxDissRate = nBSubcycleMax * cflLimitReal * fieldSafety;
            if (maxDissRate > req.cflHyperMacro) {
              dtTarget = std::min(dtTarget, (maxDissRate - req.cflHyperMacro) /
                                                req.omegaDiss);
            }
          }
          if (req.dtWhistler < std::numeric_limits<Real>::max()) {
            dtTarget = std::min(dtTarget, nBSubcycleMax * req.dtWhistler);
          }
          if (req.dtWave < std::numeric_limits<Real>::max()) {
            dtTarget = std::min(dtTarget, nBSubcycleMax * req.dtWave);
          }
        }
        dtLimitingReason = activeSubReason + " (capped at nBSubcycleMax)";
        subcycleHoldCount = 0;
      } else {
        int nTarget = std::clamp(nIdeal, nBSubcycleMin, nBSubcycleMax);
        if (nTarget > nBSubcycle) {
          // Immediate increase for numerical stability
          nSubNext = nTarget;
          subcycleHoldCount = 0;
          dtLimitingReason = activeSubReason + " (subcycling)";
        } else if (nTarget < nBSubcycle) {
          const int testSub = nBSubcycle - 1;
          bool canDecreaseSub = (testSub >= nBSubcycleMin);
          for (int iLev = 0; iLev < n_lev() && canDecreaseSub; ++iLev) {
            const auto& req = levReqs[iLev];
            const Real rateDissSub =
                (req.omegaDiss * dtTarget + req.cflHyperMacro) / testSub;
            if (rateDissSub > 0.85 * cflLimitReal * fieldSafety) {
              canDecreaseSub = false;
            }
            if (req.dtWhistler < std::numeric_limits<Real>::max() &&
                dtTarget / testSub > 0.85 * req.dtWhistler) {
              canDecreaseSub = false;
            }
            if (req.dtWave < std::numeric_limits<Real>::max() &&
                dtTarget / testSub > 0.85 * req.dtWave) {
              canDecreaseSub = false;
            }
          }
          if (canDecreaseSub) {
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

      dtNext = dtTarget;
    }
  }

  if (dtMacroMin > 0.0 && dtMacroMin < std::numeric_limits<Real>::max()) {
    maxCFL = dtNext / dtMacroMin;
  }

  if (doReport && !doCompact) {
    amrex::Print() << printPrefix << "[Hybrid Adaptive Step] dt = " << dtNext
                   << ", nSub = " << nSubNext
                   << " (limited by: " << dtLimitingReason << ")\n";
  }
}
