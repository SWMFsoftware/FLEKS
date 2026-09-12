#include <cstdlib>

#include <AMReX_ParReduce.H>

#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::update_position_to_half_stage(
    const MultiFab& nodeEMF, const MultiFab& nodeBMF, Real dt) {
  timing_func("Pts::update_position_to_half_stage");

  Real dtLoc = 0.5 * dt;

  const int iLev = 0;
  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    AoS& particles = pti.GetArrayOfStructs();

    const Box& bx = cell_status(iLev)[pti].box();
    const Array4<int const>& status = cell_status(iLev)[pti].array();

    const IntVect lowCorner = bx.smallEnd();
    const IntVect highCorner = bx.bigEnd();

    for (auto& p : particles) {
      if (p.id() < 0)
        continue;

      for (int iDim = 0; iDim < nDim; iDim++) {
        p.pos(iDim) += p.rdata(iup_ + iDim) * dtLoc;
      }

      // Mark for deletion
      if (reflect_or_delete_particle(p, status, lowCorner, highCorner, iLev)) {
        p.id() = -1;
      }
    } // for p
  } // for pti

  redistribute_particles();
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::mover(const Vector<MultiFab>& nodeE,
                                               const Vector<MultiFab>& nodeB,
                                               const Vector<MultiFab>& eBg,
                                               const Vector<MultiFab>& uBg,
                                               Real dt, Real dtNext) {
  if (is_neutral()) {
    neutral_mover(dt);
  } else {
    charged_particle_mover(nodeE, nodeB, eBg, uBg, dt, dtNext);
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::mover_cell_centered(
    const Vector<MultiFab>& centerE, const Vector<MultiFab>& centerB,
    const Vector<MultiFab>& eBg, const Vector<MultiFab>& uBg, Real dt,
    Real dtNext) {
  if (is_neutral()) {
    neutral_mover(dt);
  } else {
    charged_particle_mover_cell_centered(centerE, centerB, eBg, uBg, dt,
                                         dtNext);
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::charged_particle_mover(
    const Vector<MultiFab>& nodeE, const Vector<MultiFab>& nodeB,
    const Vector<MultiFab>& eBg, const Vector<MultiFab>& uBg, Real dt,
    Real dtNext) {
  timing_func("Pts::charged_particle_mover");

  const Real qdto2mc = charge / mass * 0.5 * dt;
  Real dtLoc = 0.5 * (dt + dtNext);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      const Array4<Real const>& nodeEArr = nodeE[iLev][pti].array();
      const Array4<Real const>& nodeBArr = nodeB[iLev][pti].array();

      const Box& bx = cell_status(iLev)[pti].box();
      const Array4<int const>& status = cell_status(iLev)[pti].array();

      const IntVect lowCorner = bx.smallEnd();
      const IntVect highCorner = bx.bigEnd();

      AoS& particles = pti.GetArrayOfStructs();

      const Dim3 lo = init_dim3(0);
      const Dim3 hi = init_dim3(1);

      for (auto& p : particles) {
        if (p.id() < 0)
          continue;

        Real up = p.rdata(iup_);
        Real vp = p.rdata(ivp_);
        Real wp = p.rdata(iwp_);
        const Real xp = p.pos(ix_);
        const Real yp = p.pos(iy_);
        const Real zp = nDim > 2 ? p.pos(iz_) : 0;

        //-----calculate interpolate coef begin-------------
        IntVect loIdx;
        RealVect dShift;

        find_node_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                        loIdx, dShift);

        Real coef[2][2][2];
        linear_interpolation_coef(dShift, coef);
        //-----calculate interpolate coef end-------------

        Real bp[3] = { 0, 0, 0 };
        Real ep[3] = { 0, 0, 0 };
        Real u0p[3] = { 0, 0, 0 };
        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + i, loIdx[iy_] + j,
                                           loIdx[iz_] + k) };

              const Real& c0 = coef[i][j][k];
              for (int iDim = 0; iDim < nDim3; iDim++) {
                bp[iDim] += nodeBArr(ijk, iDim) * c0;
                ep[iDim] += nodeEArr(ijk, iDim) * c0;
              }
            }

        up = up - u0p[ix_];
        vp = vp - u0p[iy_];
        wp = wp - u0p[iz_];

        const Real omx = qdto2mc * bp[ix_];
        const Real omy = qdto2mc * bp[iy_];
        const Real omz = qdto2mc * bp[iz_];

        // end interpolation
        const Real omsq = (omx * omx + omy * omy + omz * omz);
        const Real denom = 1.0 / (1.0 + omsq);
        // solve the position equation
        const Real ut = up + qdto2mc * ep[ix_];
        const Real vt = vp + qdto2mc * ep[iy_];
        const Real wt = wp + qdto2mc * ep[iz_];
        // const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
        const Real udotOm = ut * omx + vt * omy + wt * omz;
        // solve the velocity equation
        const Real uavg = (ut + (vt * omz - wt * omy + udotOm * omx)) * denom;
        const Real vavg = (vt + (wt * omx - ut * omz + udotOm * omy)) * denom;
        const Real wavg = (wt + (ut * omy - vt * omx + udotOm * omz)) * denom;

        Real unp1 = 2.0 * uavg - up + u0p[ix_];
        Real vnp1 = 2.0 * vavg - vp + u0p[iy_];
        Real wnp1 = 2.0 * wavg - wp + u0p[iz_];

        p.rdata(iup_) = unp1;
        p.rdata(ivp_) = vnp1;
        p.rdata(iwp_) = wnp1;

        if (pMode == PartMode::PIC && imu_ < NStructReal) {
          // Note: bp should be calculated at the new position. Now, bp at the
          // old position is used to save the calculation.
          p.rdata(imu_) = cosine(p, bp);
        }

        p.pos(ix_) = xp + unp1 * dtLoc;
        p.pos(iy_) = yp + vnp1 * dtLoc;
        if (nDim > 2)
          p.pos(iz_) = zp + wnp1 * dtLoc;

        // Apply boundary condition (absorb: delete; reflect: mirror).
        if (reflect_or_delete_particle(p, status, lowCorner, highCorner,
                                       iLev)) {
          p.id() = -1;
        }
      } // for p
    } // for pti
  }
}

//==========================================================
// Cell-centred Boris push. The E and B are gathered from cell fields. The
// gather is a plain cell-centred trilinear interpolation.

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::charged_particle_mover_cell_centered(
    const Vector<MultiFab>& centerE, const Vector<MultiFab>& centerB,
    const Vector<MultiFab>& eBg, const Vector<MultiFab>& uBg, Real dt,
    Real dtNext) {
  timing_func("Pts::charged_particle_mover_cell_centered");

  const Real qdto2mc = charge / mass * 0.5 * dt;
  Real dtLoc = 0.5 * (dt + dtNext);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      const Array4<Real const>& centerEArr = centerE[iLev][pti].array();
      const Array4<Real const>& centerBArr = centerB[iLev][pti].array();

      const Box& bx = cell_status(iLev)[pti].box();
      const Array4<int const>& status = cell_status(iLev)[pti].array();

      const IntVect lowCorner = bx.smallEnd();
      const IntVect highCorner = bx.bigEnd();

      AoS& particles = pti.GetArrayOfStructs();

      const Dim3 lo = init_dim3(0);
      const Dim3 hi = init_dim3(1);

      for (auto& p : particles) {
        if (p.id() < 0)
          continue;

        Real up = p.rdata(iup_);
        Real vp = p.rdata(ivp_);
        Real wp = p.rdata(iwp_);
        const Real xp = p.pos(ix_);
        const Real yp = p.pos(iy_);
        const Real zp = nDim > 2 ? p.pos(iz_) : 0;

        //-----calculate interpolate coef begin-------------
        IntVect loIdx;
        RealVect dShift;
        find_cell_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                        loIdx, dShift);

        // Plain cell-centred trilinear gather. The linear weights couple cells
        // loIdx and loIdx+1 (offsets 0 and 1); the 3x3x3 coef array is zero for
        // the unused offset-2 entry.
        Real coef[3][3][3];
        Real coefLin[2][2][2];
        linear_interpolation_coef(dShift, coefLin);
        for (int k = 0; k <= 2; ++k)
          for (int j = 0; j <= 2; ++j)
            for (int i = 0; i <= 2; ++i)
              coef[i][j][k] =
                  (i <= 1 && j <= 1 && k <= 1) ? coefLin[i][j][k] : 0.0;
        //-----calculate interpolate coef end-------------

        Real bp[3] = { 0, 0, 0 };
        Real ep[3] = { 0, 0, 0 };
        Real u0p[3] = { 0, 0, 0 };
        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + i, loIdx[iy_] + j,
                                           loIdx[iz_] + k) };

              const Real& c0 = coef[i - lo.x][j - lo.y][k - lo.z];
              for (int iDim = 0; iDim < nDim3; iDim++) {
                bp[iDim] += centerBArr(ijk, iDim) * c0;
                ep[iDim] += centerEArr(ijk, iDim) * c0;
              }
            }

        up = up - u0p[ix_];
        vp = vp - u0p[iy_];
        wp = wp - u0p[iz_];

        const Real omx = qdto2mc * bp[ix_];
        const Real omy = qdto2mc * bp[iy_];
        const Real omz = qdto2mc * bp[iz_];

        // end interpolation
        const Real omsq = (omx * omx + omy * omy + omz * omz);
        const Real denom = 1.0 / (1.0 + omsq);
        // solve the position equation
        const Real ut = up + qdto2mc * ep[ix_];
        const Real vt = vp + qdto2mc * ep[iy_];
        const Real wt = wp + qdto2mc * ep[iz_];
        // const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
        const Real udotOm = ut * omx + vt * omy + wt * omz;
        // solve the velocity equation
        const Real uavg = (ut + (vt * omz - wt * omy + udotOm * omx)) * denom;
        const Real vavg = (vt + (wt * omx - ut * omz + udotOm * omy)) * denom;
        const Real wavg = (wt + (ut * omy - vt * omx + udotOm * omz)) * denom;

        Real unp1 = 2.0 * uavg - up + u0p[ix_];
        Real vnp1 = 2.0 * vavg - vp + u0p[iy_];
        Real wnp1 = 2.0 * wavg - wp + u0p[iz_];

        p.rdata(iup_) = unp1;
        p.rdata(ivp_) = vnp1;
        p.rdata(iwp_) = wnp1;

        if (pMode == PartMode::PIC && imu_ < NStructReal) {
          // Note: bp should be calculated at the new position. Now, bp at the
          // old position is used to save the calculation.
          p.rdata(imu_) = cosine(p, bp);
        }

        p.pos(ix_) = xp + unp1 * dtLoc;
        p.pos(iy_) = yp + vnp1 * dtLoc;
        if (nDim > 2)
          p.pos(iz_) = zp + wnp1 * dtLoc;

        // Apply boundary condition (absorb: delete; reflect: mirror).
        if (reflect_or_delete_particle(p, status, lowCorner, highCorner,
                                       iLev)) {
          p.id() = -1;
        }
      } // for p
    } // for pti
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::neutral_mover(Real dt) {
  timing_func("Pts::neutral_mover");

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      AoS& particles = pti.GetArrayOfStructs();

      const Box& bx = cell_status(iLev)[pti].box();
      const Array4<int const>& status = cell_status(iLev)[pti].array();

      const IntVect lowCorner = bx.smallEnd();
      const IntVect highCorner = bx.bigEnd();
      for (auto& p : particles) {
        if (p.id() < 0)
          continue;

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
        if (reflect_or_delete_particle(p, status, lowCorner, highCorner,
                                       iLev)) {
          p.id() = -1;
        }
      } // for p
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
  Real epsMax = 0;

  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const& phiArr = phiMF[iLev][pti].array();
    const Array4<int const>& status = cell_status(iLev)[pti].array();

    AoS& particles = pti.GetArrayOfStructs();

    const Box& bx = cell_status(iLev)[pti].box();
    const IntVect lowCorner = bx.smallEnd();
    const IntVect highCorner = bx.bigEnd();

    for (auto& p : particles) {
      if (p.id() == -1 ||
          is_outside_active_region(p, status, lowCorner, highCorner, iLev)) {
        p.id() = -1;
        continue;
      }

      if (skip_particle_for_dive_cleaning(p.pos(), Geom(iLev), iLev, status) &&
          n_lev() > 1) {
        continue;
      }

      IntVect loIdx;
      RealVect dShift;
      find_cell_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                      loIdx, dShift);

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
            if (bit::is_lev_boundary(status(ijk)))
              isBoundaryPhysicalCell = true;
          }
      if (isBoundaryPhysicalCell && iLev == 0)
        continue;

      {
        Real weights_IIID[2][2][2][nDim3];
        //----- Mass matrix calculation begin--------------
        const Real xi0 = dShift[ix_] * dx[iLev][ix_];
        const Real eta0 = dShift[iy_] * dx[iLev][iy_];
        const Real zeta0 = nDim > 2 ? dShift[iz_] * dx[iLev][iz_] : 0;
        const Real xi1 = dx[iLev][ix_] - xi0;
        const Real eta1 = dx[iLev][iy_] - eta0;
        const Real zeta1 = nDim > 2 ? dx[iLev][iz_] - zeta0 : 1;

        const Real zeta02Vol = zeta0 * invVol[iLev];
        const Real zeta12Vol = zeta1 * invVol[iLev];
        const Real eta02Vol = eta0 * invVol[iLev];
        const Real eta12Vol = eta1 * invVol[iLev];

        weights_IIID[1][1][1][ix_] = eta0 * zeta02Vol;
        weights_IIID[1][1][1][iy_] = xi0 * zeta02Vol;
        weights_IIID[1][1][1][iz_] = xi0 * eta02Vol;

        // xi0*eta0*zeta1*invVol[iLev];
        weights_IIID[1][1][0][ix_] = eta0 * zeta12Vol;
        weights_IIID[1][1][0][iy_] = xi0 * zeta12Vol;
        weights_IIID[1][1][0][iz_] = -xi0 * eta02Vol;

        // xi0*eta1*zeta0*invVol[iLev];
        weights_IIID[1][0][1][ix_] = eta1 * zeta02Vol;
        weights_IIID[1][0][1][iy_] = -xi0 * zeta02Vol;
        weights_IIID[1][0][1][iz_] = xi0 * eta12Vol;

        // xi0*eta1*zeta1*invVol[iLev];
        weights_IIID[1][0][0][ix_] = eta1 * zeta12Vol;
        weights_IIID[1][0][0][iy_] = -xi0 * zeta12Vol;
        weights_IIID[1][0][0][iz_] = -xi0 * eta12Vol;

        // xi1*eta0*zeta0*invVol[iLev];
        weights_IIID[0][1][1][ix_] = -eta0 * zeta02Vol;
        weights_IIID[0][1][1][iy_] = xi1 * zeta02Vol;
        weights_IIID[0][1][1][iz_] = xi1 * eta02Vol;

        // xi1*eta0*zeta1*invVol[iLev];
        weights_IIID[0][1][0][ix_] = -eta0 * zeta12Vol;
        weights_IIID[0][1][0][iy_] = xi1 * zeta12Vol;
        weights_IIID[0][1][0][iz_] = -xi1 * eta02Vol;

        // xi1*eta1*zeta0*invVol[iLev];
        weights_IIID[0][0][1][ix_] = -eta1 * zeta02Vol;
        weights_IIID[0][0][1][iy_] = -xi1 * zeta02Vol;
        weights_IIID[0][0][1][iz_] = xi1 * eta12Vol;

        // xi1*eta1*zeta1*invVol[iLev];
        weights_IIID[0][0][0][ix_] = -eta1 * zeta12Vol;
        weights_IIID[0][0][0][iy_] = -xi1 * zeta12Vol;
        weights_IIID[0][0][0][iz_] = -xi1 * eta12Vol;

        RealVect eps_D = { AMREX_D_DECL(0, 0, 0) };

        // Do not shift along z direction for both 2D and fake 2D cases.
        int nD = isFake2D ? 2 : nDim;

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
          if (fabs(eps_D[iDim] * invDx[iLev][iDim]) > eps_D_dot_invDx_Max) {
            eps_D_dot_invDx_Max = fabs(eps_D[iDim] * invDx[iLev][iDim]);
          }
        }

        if (eps_D_dot_invDx_Max > epsLimit) {
          // If eps_D is too large, the underlying assumption of the particle
          // correction method will be not valid. Comparing each exp_D
          // component instead of the length dl saves the computational time.
          const Real dl = eps_D.vectorLength();
          const Real ratio = epsLimit * dx[iLev][ix_] / dl;
          for (int iDim = 0; iDim < nDim; iDim++)
            eps_D[iDim] *= ratio;
        }

        for (int iDim = 0; iDim < nDim; iDim++) {
          if (fabs(eps_D[iDim] * invDx[iLev][iDim]) > epsMax)
            epsMax = fabs(eps_D[iDim] * invDx[iLev][iDim]);

          p.pos(iDim) += eps_D[iDim];
        }

        if (is_outside_active_region(p, status, lowCorner, highCorner, iLev)) {
          // Do not allow moving particles from physical cells to ghost cells
          // during divE correction.
          for (int iDim = 0; iDim < nDim; iDim++) {
            p.pos(iDim) -= eps_D[iDim];
          }

          // p.id() = -1;
        }
      }
    } // for p
  }
}

//==========================================================

#define INSTANTIATE_PARTICLES_MOVER(T)                                         \
  template void T::update_position_to_half_stage(const MultiFab&,              \
                                                 const MultiFab&, Real);       \
  template void T::mover(const Vector<MultiFab>&, const Vector<MultiFab>&,     \
                         const Vector<MultiFab>&, const Vector<MultiFab>&,     \
                         Real, Real);                                          \
  template void T::mover_cell_centered(                                        \
      const Vector<MultiFab>&, const Vector<MultiFab>&,                        \
      const Vector<MultiFab>&, const Vector<MultiFab>&, Real, Real);           \
  template void T::charged_particle_mover(                                     \
      const Vector<MultiFab>&, const Vector<MultiFab>&,                        \
      const Vector<MultiFab>&, const Vector<MultiFab>&, Real, Real);           \
  template void T::charged_particle_mover_cell_centered(                       \
      const Vector<MultiFab>&, const Vector<MultiFab>&,                        \
      const Vector<MultiFab>&, const Vector<MultiFab>&, Real, Real);           \
  template void T::neutral_mover(Real);                                        \
  template void T::divE_correct_position(const Vector<MultiFab>&, int);

INSTANTIATE_PARTICLES_MOVER(PicParticles)
INSTANTIATE_PARTICLES_MOVER(PTParticles)
#undef INSTANTIATE_PARTICLES_MOVER
