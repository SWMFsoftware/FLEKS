#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <AMReX_ParReduce.H>

#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

namespace {

// Charge-conserving (Esirkepov) current of one particle of charge q moving
// from x0 to x1 (at most one cell per dimension) with the node-CIC shape S,
// deposited onto the nodes of jArr (components x, y, z):
//   face current  F_d = -q dx_d / (dt V) cumsum_d W_d   (faces between nodes)
//   node current  J_d = (F_d(i - 1/2) + F_d(i + 1/2)) / 2 along d.
// With the cell-centered charge rho_c = (node-CIC charge averaged over the
// 2^D nodes of each cell), rho_c^{new} - rho_c^{old} = -dt
// div_node_to_center(J) holds to round-off. In 2D the out-of-plane component is
// q vOut / V times the time-averaged shape. Returns false if the particle moved
// more than one cell.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool deposit_esirkepov(
    const Array4<Real>& jArr, const Real* x0, const Real* x1, Real q, Real vOut,
    const Real* plo, const Real* invDx, const Real* dx, Real invDt,
    Real invVol) {
  // A particle that stays in its cell touches 2 nodes per dimension, one that
  // crosses a cell face 3.
  constexpr int nS = 3;
  int base[3] = { 0, 0, 0 };
  int n[3] = { 1, 1, 1 };
  Real S0[3][nS] = {}, DS[3][nS] = {};
  for (int d = 0; d < 3; ++d) {
    if (d >= nDim) {
      S0[d][0] = 1.0; // a single node with no motion
      continue;
    }
    const Real s0 = (x0[d] - plo[d]) * invDx[d];
    const Real s1 = (x1[d] - plo[d]) * invDx[d];
    const int i0 = static_cast<int>(std::floor(s0));
    const int i1 = static_cast<int>(std::floor(s1));
    if (amrex::Math::abs(i1 - i0) > 1)
      return false;
    base[d] = amrex::min(i0, i1);
    n[d] = 2 + amrex::Math::abs(i1 - i0);
    const Real f0 = s0 - i0, f1 = s1 - i1;
    Real S1[nS] = {};
    S0[d][i0 - base[d]] += 1.0 - f0;
    S0[d][i0 - base[d] + 1] += f0;
    S1[i1 - base[d]] += 1.0 - f1;
    S1[i1 - base[d] + 1] += f1;
    for (int a = 0; a < n[d]; ++a)
      DS[d][a] = S1[a] - S0[d][a];
  }

  // Shape weight transverse to direction d: product of the time-averaged
  // shapes of the other two directions (Esirkepov's 1/2 and 1/3 terms).
  auto transverse = [&](int e, int f, int ie, int jf) {
    return S0[e][ie] * S0[f][jf] + 0.5 * DS[e][ie] * S0[f][jf] +
           0.5 * S0[e][ie] * DS[f][jf] + DS[e][ie] * DS[f][jf] / 3.0;
  };

  for (int d = 0; d < nDim; ++d) {
    const int e = (d + 1) % 3, f = (d + 2) % 3;
    const Real coef = -q * dx[d] * invDt * invVol;
    for (int ie = 0; ie < n[e]; ++ie)
      for (int jf = 0; jf < n[f]; ++jf) {
        const Real wt = transverse(e, f, ie, jf);
        Real fPrev = 0.0, fCum = 0.0;
        for (int a = 0; a < n[d]; ++a) {
          fCum += coef * DS[d][a] * wt;
          const Real jNode = 0.5 * (fCum + fPrev);
          fPrev = fCum;
          if (jNode == 0.0)
            continue;
          int idx[3];
          idx[d] = base[d] + a;
          idx[e] = base[e] + ie;
          idx[f] = base[f] + jf;
          deposit_add(jArr.ptr(AMREX_D_DECL(idx[0], idx[1], idx[2]), d), jNode);
        }
      }
  }

  if (nDim == 2) {
    const Real coef = q * vOut * invVol;
    for (int a = 0; a < n[0]; ++a)
      for (int b = 0; b < n[1]; ++b) {
        const Real w = transverse(0, 1, a, b);
        if (w != 0.0)
          deposit_add(jArr.ptr(AMREX_D_DECL(base[0] + a, base[1] + b, 0), iz_),
                      coef * w);
      }
  }
  return true;
}

// E at a particle with the shape of the charge-conserving current: along
// component d, the average of the two nodes of the particle's cell; across d
// (and for the out-of-plane component in 2D), linear. loIdx and dShift come
// from find_node_interpolation (dShift = weight of the upper node).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gather_E_matched_to_current(
    const Array4<Real const>& EArr, const IntVect& loIdx,
    const RealVect& dShift, Real* ep) {
  for (int c = 0; c < 3; ++c) {
    Real w[3][2];
    for (int d = 0; d < 3; ++d) {
      if (d >= nDim) {
        w[d][0] = 1.0;
        w[d][1] = 0.0;
      } else if (d == c) {
        w[d][0] = 0.5;
        w[d][1] = 0.5;
      } else {
        w[d][0] = 1.0 - dShift[d];
        w[d][1] = dShift[d];
      }
    }
    Real sum = 0.0;
    for (int k = 0; k < (nDim > 2 ? 2 : 1); ++k)
      for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i) {
          const IntVect idx(AMREX_D_DECL(loIdx[0] + i, loIdx[1] + j,
                                         nDim > 2 ? loIdx[nDim - 1] + k : 0));
          sum += w[0][i] * w[1][j] * w[2][k] * EArr(idx, c);
        }
    ep[c] = sum;
  }
}
} // namespace

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::update_position_to_half_stage(
    const MultiFab& nodeEMF, const MultiFab& nodeBMF, Real dt) {
  timing_func("Pts::update_position_to_half_stage");

  const Real dtLoc = 0.5 * dt;

  const int iLev = 0;
  const LevelGeomBox lb = make_level_geom_box(iLev);
  const ParticleBoundary pb = particle_boundary();
  int lastFab = -1;
  ActiveRegionBox ab;

  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    ParticleType* const pData = pti.GetArrayOfStructs()().data();

    if (pti.index() != lastFab) {
      lastFab = pti.index();
      const Box& bx = cell_status(iLev)[pti].box();
      ab = make_active_region_box(lb, pti.validbox(), bx);
    }

    const Array4<int const> status = cell_status(iLev)[pti].const_array();

    ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE(int ip) {
      ParticleType& p = pData[ip];
      if (p.id() < 0)
        return;

      for (int iDim = 0; iDim < nDim; iDim++) {
        p.pos(iDim) += p.rdata(iup_ + iDim) * dtLoc;
      }

      // Mark for deletion
      if (pb.reflect_or_delete_particle(p, status, ab)) {
        p.id() = -1;
      }
    });
  } // for pti

  redistribute_particles();
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::mover(
    const Vector<MultiFab>& nodeE, const Vector<MultiFab>& nodeB,
    const Vector<MultiFab>& nodeB0, const Vector<MultiFab>& eBg,
    const Vector<MultiFab>& uBg, Real dt, Real dtNext, int pushStage,
    Vector<MultiFab>* currentMF, bool trialOnly) {
  if (is_neutral()) {
    if (currentMF) {
      for (auto& mf : *currentMF)
        mf.setVal(0.0);
    }
    if (!trialOnly)
      neutral_mover(dt);
  } else {
    charged_particle_mover(nodeE, nodeB, nodeB0, eBg, uBg, dt, dtNext,
                           pushStage, currentMF, trialOnly);
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::charged_particle_mover(
    const Vector<MultiFab>& nodeE, const Vector<MultiFab>& nodeB,
    const Vector<MultiFab>& nodeB0, const Vector<MultiFab>& eBg,
    const Vector<MultiFab>& uBg, Real dt, Real dtNext, int pushStage,
    Vector<MultiFab>* currentMF, bool trialOnly) {
  timing_func("Pts::charged_particle_mover");

  // nodeB0 is empty when no intrinsic field is configured.
  const bool useB0 = nodeB0.size() == nodeB.size() && !nodeB0[0].empty();

  const Real qdto2mc = charge / mass * 0.5 * dt;
  Real dtLoc = 0.5 * (dt + dtNext);

  // Current of the step (Esirkepov), deposited as the particles move.
  const Real qFactor = charge * qomSign;
  const Real invDt = 1.0 / dt;
  const bool depositCurrent = (currentMF != nullptr);
  const bool savePitchAngle = (pMode == PartMode::PIC && imu_ < NStructReal);
  const ParticleBoundary pb = particle_boundary();

  // Set by a kernel if a particle moved more than one cell.
  Gpu::DeviceVector<int> movedTooFarFlag(1, 0);
  int* const movedTooFar = movedTooFarFlag.data();

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    // Copies, so that the kernels capture values, not host pointers.
    Real ploLoc[3] = { 0, 0, 0 }, invDxLoc[3] = { 0, 0, 0 };
    Real dxLoc[3] = { 0, 0, 0 };
    for (int d = 0; d < nDim; ++d) {
      ploLoc[d] = plo[iLev][d];
      invDxLoc[d] = invDx[iLev][d];
      dxLoc[d] = Geom(iLev).CellSize(d);
    }
    const Real invVolLev = invVol[iLev];
    const LevelGeomBox lb = make_level_geom_box(iLev);
    int lastFab = -1;
    ActiveRegionBox ab;

    if (currentMF)
      (*currentMF)[iLev].setVal(0.0);

    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      const Array4<Real> jArr =
          currentMF ? (*currentMF)[iLev][pti].array() : Array4<Real>();
      const Array4<Real const> EArr = nodeE[iLev][pti].const_array();
      const Array4<Real const> BArr = nodeB[iLev][pti].const_array();
      const Array4<Real const> B0Arr =
          useB0 ? nodeB0[iLev][pti].const_array() : Array4<Real const>();

      if (pti.index() != lastFab) {
        lastFab = pti.index();
        const Box& bx = cell_status(iLev)[pti].box();
        ab = make_active_region_box(lb, pti.validbox(), bx);
      }

      const Array4<int const> status = cell_status(iLev)[pti].const_array();

      ParticleType* const pData = pti.GetArrayOfStructs()().data();

      const Dim3 lo = init_dim3(0);
      const Dim3 hi = init_dim3(1);

      ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE(int ip) {
        ParticleType& p = pData[ip];
        if (p.id() < 0)
          return;

        const Real u_start = p.rdata(iup_);
        const Real v_start = p.rdata(ivp_);
        const Real w_start = p.rdata(iwp_);
        Real up = u_start;
        Real vp = v_start;
        Real wp = w_start;
        const Real xp = p.pos(ix_);
        const Real yp = p.pos(iy_);
        const Real zp = nDim > 2 ? p.pos(iz_) : 0;

        //-----calculate interpolate coef begin-------------
        // The stencil is centered on the node that contains the particle: plain
        // trilinear weights over the surrounding 2x2x2 nodes.
        IntVect loIdx;
        RealVect dShift;

        Real coef[2][2][2];
        RealVect posEval = p.pos();
        if (pushStage == 2) {
          posEval[ix_] += 0.5 * u_start * dt;
          posEval[iy_] += 0.5 * v_start * dt;
          if (nDim > 2)
            posEval[iz_] += 0.5 * w_start * dt;
        }
        find_node_interpolation(posEval, ploLoc, invDxLoc, loIdx, dShift, coef);
        //-----calculate interpolate coef end-------------

        Real bp[3] = { 0, 0, 0 };
        Real ep[3] = { 0, 0, 0 };
        Real u0p[3] = { 0, 0, 0 };
        Real b0p[3] = { 0, 0, 0 };
        if (useB0) {
          const Array4<Real const> fields[3] = { BArr, EArr, B0Arr };
          Real* values[3] = { bp, ep, b0p };
          interpolate_vector_fields(fields, loIdx, coef, lo, hi, values);
          bp[ix_] += b0p[ix_];
          bp[iy_] += b0p[iy_];
          bp[iz_] += b0p[iz_];
        } else {
          const Array4<Real const> fields[2] = { BArr, EArr };
          Real* values[2] = { bp, ep };
          interpolate_vector_fields(fields, loIdx, coef, lo, hi, values);
        }

        // The charge-conserving current of E_d has a zeroth-order shape along
        // d (the average of the two nodes of the cell the particle is in).
        // Gather E with the same shape, the transpose of the deposit, so the
        // work of the field on the particles matches the field energy loss.
        if (depositCurrent)
          gather_E_matched_to_current(EArr, loIdx, dShift, ep);

        up = up - u0p[ix_];
        vp = vp - u0p[iy_];
        wp = wp - u0p[iz_];

        const Real omx = qdto2mc * bp[ix_];
        const Real omy = qdto2mc * bp[iy_];
        const Real omz = qdto2mc * bp[iz_];

        // end interpolation
        const Real velocity[3] = { up, vp, wp };
        const Real ut = up + qdto2mc * ep[ix_];
        const Real vt = vp + qdto2mc * ep[iy_];
        const Real wt = wp + qdto2mc * ep[iz_];
        const Real electricVelocity[3] = { ut, vt, wt };
        const Real omega[3] = { omx, omy, omz };
        Real updatedVelocity[3];
        boris_push_nonrelativistic(velocity, electricVelocity, omega,
                                   updatedVelocity);

        const Real unp1 = updatedVelocity[ix_] + u0p[ix_];
        const Real vnp1 = updatedVelocity[iy_] + u0p[iy_];
        const Real wnp1 = updatedVelocity[iz_] + u0p[iz_];

        // Drift velocity: v^{n+1/2} (leapfrog) or the average of the start
        // and end velocities (predictor-corrector stages).
        const Real uDrift = pushStage > 0 ? 0.5 * (u_start + unp1) : unp1;
        const Real vDrift = pushStage > 0 ? 0.5 * (v_start + vnp1) : vnp1;
        const Real wDrift = pushStage > 0 ? 0.5 * (w_start + wnp1) : wnp1;
        const Real xOld[3] = { xp, yp, zp };
        const Real xNew[3] = { xp + uDrift * dtLoc, yp + vDrift * dtLoc,
                               zp + wDrift * dtLoc };

        // J = q dx/dt for every component; in 2D the out-of-plane dx is
        // wDrift * dtLoc (dtLoc != dt when the time step changes).
        if (depositCurrent &&
            !deposit_esirkepov(jArr, xOld, xNew, qFactor * p.rdata(iqp_),
                               wDrift * dtLoc * invDt, ploLoc, invDxLoc, dxLoc,
                               invDt, invVolLev))
          *movedTooFar = 1;

        // A trial push only deposits the current of the trial state.
        if (trialOnly)
          return;

        p.rdata(iup_) = unp1;
        p.rdata(ivp_) = vnp1;
        p.rdata(iwp_) = wnp1;

        if (savePitchAngle) {
          // Note: bp should be calculated at the new position. Now, bp at the
          // old position is used to save the calculation.
          p.rdata(imu_) = cosine(p, bp);
        }

        p.pos(ix_) = xNew[ix_];
        p.pos(iy_) = xNew[iy_];
        if (nDim > 2)
          p.pos(iz_) = xNew[iz_];

        // Apply boundary condition (absorb: delete; reflect: mirror).
        if (pb.reflect_or_delete_particle(p, status, ab)) {
          p.id() = -1;
        }
      });
    } // for pti

    int movedTooFarHost = 0;
    Gpu::copy(Gpu::deviceToHost, movedTooFarFlag.begin(), movedTooFarFlag.end(),
              &movedTooFarHost);
    if (movedTooFarHost)
      Abort("charged_particle_mover: a particle moved more than one cell in a "
            "step; the charge-conserving current needs a smaller time step.");

    if (currentMF) {
      fold_wall_ghost_nodes((*currentMF)[iLev], iLev, true);
      (*currentMF)[iLev].SumBoundary(Geom(iLev).periodicity(),
                                     /*deterministic=*/true);
    }
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::neutral_mover(Real dt) {
  timing_func("Pts::neutral_mover");

  const ParticleBoundary pb = particle_boundary();
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const LevelGeomBox lb = make_level_geom_box(iLev);
    int lastFab = -1;
    ActiveRegionBox ab;

    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      ParticleType* const pData = pti.GetArrayOfStructs()().data();

      if (pti.index() != lastFab) {
        lastFab = pti.index();
        const Box& bx = cell_status(iLev)[pti].box();
        ab = make_active_region_box(lb, pti.validbox(), bx);
      }

      const Array4<int const> status = cell_status(iLev)[pti].const_array();

      ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE(int ip) {
        ParticleType& p = pData[ip];
        if (p.id() < 0)
          return;

        const Real up = p.rdata(iup_);
        const Real vp = p.rdata(ivp_);
        const Real wp = p.rdata(iwp_);
        const Real xp = p.pos(ix_);
        const Real yp = p.pos(iy_);
        const Real zp = p.pos(iz_);

        p.pos(ix_) = xp + up * dt;
        p.pos(iy_) = yp + vp * dt;
        p.pos(iz_) = zp + wp * dt;

        // Apply boundary condition (absorb: delete; reflect: mirror).
        if (pb.reflect_or_delete_particle(p, status, ab)) {
          p.id() = -1;
        }
      });
    } // for pti
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::divE_correct_position(
    const amrex::Vector<MultiFab>& phiMF, int iLev) {
  timing_func("Pts:divE_correct_position");

  const Real sign = charge / fabs(charge);
  const Real epsLimit = 0.1;

  // Copies, so that the kernels capture values, not host pointers.
  Real ploLoc[3] = { 0, 0, 0 }, invDxLoc[3] = { 0, 0, 0 };
  Real dxLoc[3] = { 0, 0, 0 };
  for (int d = 0; d < nDim; ++d) {
    ploLoc[d] = plo[iLev][d];
    invDxLoc[d] = invDx[iLev][d];
    dxLoc[d] = dx[iLev][d];
  }
  const Real invVolLev = invVol[iLev];
  const bool useBody = (grid != nullptr && grid->use_body());
  const bool multiLevel = n_lev() > 1;
  const int nD = get_dim();
  const ParticleBoundary pb = particle_boundary();
  const LevelGeomBox lb = make_level_geom_box(iLev);
  int lastFab = -1;
  ActiveRegionBox ab;

  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const phiArr = phiMF[iLev][pti].const_array();
    const Array4<int const> status = cell_status(iLev)[pti].const_array();

    ParticleType* const pData = pti.GetArrayOfStructs()().data();

    if (pti.index() != lastFab) {
      lastFab = pti.index();
      const Box& bx = cell_status(iLev)[pti].box();
      ab = make_active_region_box(lb, pti.validbox(), bx);
    }

    ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE(int ip) {
      ParticleType& p = pData[ip];
      if (p.id() == -1 || pb.is_outside_active_region(p, status, ab)) {
        p.id() = -1;
        return;
      }

      if (skip_particle_for_dive_cleaning(p.pos(), ploLoc, invDxLoc, iLev,
                                          status) &&
          multiLevel) {
        return;
      }
      IntVect loIdx;
      RealVect dShift;
      find_cell_index(p.pos(), ploLoc, invDxLoc, loIdx, dShift);

      // Since the boundary condition for solving phi is not perfect,
      // correcting particles that are close to the boundaries may produce
      // artificial oscillations, which are seen in Earth's magnetotail
      // simulations. So, it is better to skip the boundary physical cells.
      bool isBoundaryPhysicalCell = false;
      for (int iz = 0; iz <= 1; iz++)
        for (int iy = 0; iy <= 1; iy++)
          for (int ix = 0; ix <= 1; ix++) {
            IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + ix, loIdx[iy_] + iy,
                                         loIdx[iz_] + iz) };
            if (bit::is_lev_boundary(status(ijk)) ||
                (useBody && bit::is_body(status(ijk))))
              isBoundaryPhysicalCell = true;
          }
      if (isBoundaryPhysicalCell && iLev == 0)
        return;

      {
        Real weights_IIID[2][2][2][nDim3];
        //----- Mass matrix calculation begin--------------
        const Real xi0 = dShift[ix_] * dxLoc[ix_];
        const Real eta0 = dShift[iy_] * dxLoc[iy_];
        const Real zeta0 = nDim > 2 ? dShift[iz_] * dxLoc[iz_] : 0;
        const Real xi1 = dxLoc[ix_] - xi0;
        const Real eta1 = dxLoc[iy_] - eta0;
        const Real zeta1 = nDim > 2 ? dxLoc[iz_] - zeta0 : 1;

        const Real zeta02Vol = zeta0 * invVolLev;
        const Real zeta12Vol = zeta1 * invVolLev;
        const Real eta02Vol = eta0 * invVolLev;
        const Real eta12Vol = eta1 * invVolLev;

        weights_IIID[1][1][1][ix_] = eta0 * zeta02Vol;
        weights_IIID[1][1][1][iy_] = xi0 * zeta02Vol;
        weights_IIID[1][1][1][iz_] = xi0 * eta02Vol;

        // xi0*eta0*zeta1*invVolLev;
        weights_IIID[1][1][0][ix_] = eta0 * zeta12Vol;
        weights_IIID[1][1][0][iy_] = xi0 * zeta12Vol;
        weights_IIID[1][1][0][iz_] = -xi0 * eta02Vol;

        // xi0*eta1*zeta0*invVolLev;
        weights_IIID[1][0][1][ix_] = eta1 * zeta02Vol;
        weights_IIID[1][0][1][iy_] = -xi0 * zeta02Vol;
        weights_IIID[1][0][1][iz_] = xi0 * eta12Vol;

        // xi0*eta1*zeta1*invVolLev;
        weights_IIID[1][0][0][ix_] = eta1 * zeta12Vol;
        weights_IIID[1][0][0][iy_] = -xi0 * zeta12Vol;
        weights_IIID[1][0][0][iz_] = -xi0 * eta12Vol;

        // xi1*eta0*zeta0*invVolLev;
        weights_IIID[0][1][1][ix_] = -eta0 * zeta02Vol;
        weights_IIID[0][1][1][iy_] = xi1 * zeta02Vol;
        weights_IIID[0][1][1][iz_] = xi1 * eta02Vol;

        // xi1*eta0*zeta1*invVolLev;
        weights_IIID[0][1][0][ix_] = -eta0 * zeta12Vol;
        weights_IIID[0][1][0][iy_] = xi1 * zeta12Vol;
        weights_IIID[0][1][0][iz_] = -xi1 * eta02Vol;

        // xi1*eta1*zeta0*invVolLev;
        weights_IIID[0][0][1][ix_] = -eta1 * zeta02Vol;
        weights_IIID[0][0][1][iy_] = -xi1 * zeta02Vol;
        weights_IIID[0][0][1][iz_] = xi1 * eta12Vol;

        // xi1*eta1*zeta1*invVolLev;
        weights_IIID[0][0][0][ix_] = -eta1 * zeta12Vol;
        weights_IIID[0][0][0][iy_] = -xi1 * zeta12Vol;
        weights_IIID[0][0][0][iz_] = -xi1 * eta12Vol;

        RealVect eps_D = { AMREX_D_DECL(0, 0, 0) };

        // Do not shift along z direction for both 2D and fake 2D cases (nD).

#if AMREX_SPACEDIM > 2
        constexpr int kEnd = 1;
#else
        constexpr int kEnd = 0;
#endif
        for (int k = 0; k <= kEnd; ++k) {
          for (int j = 0; j <= 1; ++j) {
            for (int i = 0; i <= 1; ++i) {
              IntVect ijk = { AMREX_D_DECL(i, j, k) };
              const Real coef = phiArr(loIdx + ijk);
              for (int iDim = 0; iDim < nD; iDim++) {
                eps_D[iDim] += coef * weights_IIID[i][j][k][iDim];
              }
            }
          }
        }

        for (int iDim = 0; iDim < nDim; iDim++)
          eps_D[iDim] *= sign * fourPI;

        Real eps_D_dot_invDx_Max = 0.0;
        for (int iDim = 0; iDim < nDim; iDim++) {
          if (fabs(eps_D[iDim] * invDxLoc[iDim]) > eps_D_dot_invDx_Max) {
            eps_D_dot_invDx_Max = fabs(eps_D[iDim] * invDxLoc[iDim]);
          }
        }

        if (eps_D_dot_invDx_Max > epsLimit) {
          // If eps_D is too large, the underlying assumption of the particle
          // correction method will be not valid. Comparing each exp_D
          // component instead of the length dl saves the computational time.
          const Real dl = eps_D.vectorLength();
          const Real ratio = epsLimit * dxLoc[ix_] / dl;
          for (int iDim = 0; iDim < nDim; iDim++)
            eps_D[iDim] *= ratio;
        }

        for (int iDim = 0; iDim < nDim; iDim++) {
          p.pos(iDim) += eps_D[iDim];
        }

        if (pb.is_outside_active_region(p, status, ab)) {
          // Do not allow moving particles from physical cells to ghost cells
          // during divE correction.
          for (int iDim = 0; iDim < nDim; iDim++) {
            p.pos(iDim) -= eps_D[iDim];
          }

          // p.id() = -1;
        }
      }
    });
  }
}

//==========================================================

#define INSTANTIATE_PARTICLES_MOVER(T)                                         \
  template void T::update_position_to_half_stage(const MultiFab&,              \
                                                 const MultiFab&, Real);       \
  template void T::mover(const Vector<MultiFab>&, const Vector<MultiFab>&,     \
                         const Vector<MultiFab>&, const Vector<MultiFab>&,     \
                         const Vector<MultiFab>&, Real, Real, int,             \
                         Vector<MultiFab>*, bool);                             \
  template void T::charged_particle_mover(                                     \
      const Vector<MultiFab>&, const Vector<MultiFab>&,                        \
      const Vector<MultiFab>&, const Vector<MultiFab>&,                        \
      const Vector<MultiFab>&, Real, Real, int, Vector<MultiFab>*, bool);      \
  template void T::neutral_mover(Real);                                        \
  template void T::divE_correct_position(const Vector<MultiFab>&, int);

INSTANTIATE_PARTICLES_MOVER(PicParticles)
INSTANTIATE_PARTICLES_MOVER(PTParticles)
#undef INSTANTIATE_PARTICLES_MOVER
