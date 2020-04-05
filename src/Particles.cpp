#include <algorithm>

#include "Particles.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

//==========================================================
Particles::Particles(const Geometry& geom, const DistributionMapping& dm,
                     const BoxArray& ba, TimeCtr* const tcIn,
                     const int speciesIDIn, const Real chargeIn,
                     const Real massIn, const IntVect& nPartPerCellIn)
    : ParticleContainer<4, 0, 0, 0>(geom, dm, ba),
      tc(tcIn),
      speciesID(speciesIDIn),
      charge(chargeIn),
      mass(massIn),
      nPartPerCell(nPartPerCellIn) {
  do_tiling = true;

  qom = charge / mass;
  qomSign = qom > 0 ? 1 : -1;

  for (int i = 0; i < nDimMax; i++)
    tile_size[i] = 1;
}

//==========================================================
void Particles::add_particles_cell(const MFIter& mfi,
                                   const FluidInterface& fluidInterface, int i,
                                   int j, int k) {
  int ig, jg, kg, nxcg, nycg, nzcg, iCycle, npcel, nRandom = 7;
  // Why +1? for comparison with iPIC3D.-----

  ig = i + 2;
  jg = j + 2;
  kg = k;
  if (fluidInterface.getnDim() > 2)
    kg = kg + 2; // just for comparison with iPIC3D;
  //----------------------------------------

  nxcg = fluidInterface.getFluidNxc() + 2;
  nycg = fluidInterface.getFluidNyc() + 2;
  nzcg = fluidInterface.getFluidNzc();
  if (fluidInterface.getnDim() > 2)
    nzcg += 2;

  iCycle = tc->get_cycle();
  npcel = nPartPerCell[ix_] * nPartPerCell[iy_] * nPartPerCell[iz_];

  // What if the seed overflow?
  const long seed =
      (speciesID + 3) * nRandom * npcel *
      (nxcg * nycg * nzcg * iCycle + nycg * nzcg * ig + nzcg * jg + kg);
  randNum.set_seed(seed);

  Real x, y, z; // Particle location.

  auto dx = Geom(0).CellSize();
  auto plo = Geom(0).ProbLo();

  const Real vol = dx[ix_] * dx[iy_] * dx[iz_];
  const Real vol2Npcel = qomSign * vol / npcel;

  const int lev = 0;
  auto& particles =
      GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

  //----------------------------------------------------------
  // Calculate the coefficient to correct the thermal velocity so that the
  // variation of the velocity space distribution is unbiased.
  int nPartEffective;
  {
    const Box& gbx = convert(Geom(0).Domain(), { 0, 0, 0 });
    const bool is2D = gbx.bigEnd(iz_) == gbx.smallEnd(iz_);
    int nCellContribute = is2D ? 4 : 8;
    const int nx = nPartPerCell[ix_];
    const int ny = nPartPerCell[iy_];
    const int nz = nPartPerCell[iz_];
    const Real coefCorrection = 27. / 8 * (nx + 1) * (ny + 1) * (nz + 1) /
                                ((2 * nx + 1) * (2 * ny + 1) * (2 * nz + 1));
    nPartEffective = nCellContribute * npcel * coefCorrection;
  }
  const Real coefSD = sqrt(Real(nPartEffective) / (nPartEffective - 1));
  //-----------------------------------------------------------

  // loop over particles inside grid cell i, j, k
  for (int ii = 0; ii < nPartPerCell[ix_]; ii++)
    for (int jj = 0; jj < nPartPerCell[iy_]; jj++)
      for (int kk = 0; kk < nPartPerCell[iz_]; kk++) {

        x = (ii + randNum()) * (dx[ix_] / nPartPerCell[ix_]) + i * dx[ix_] +
            plo[ix_];
        y = (jj + randNum()) * (dx[iy_] / nPartPerCell[iy_]) + j * dx[iy_] +
            plo[iy_];
        z = (kk + randNum()) * (dx[iz_] / nPartPerCell[iz_]) + k * dx[iz_] +
            plo[iz_];

        double q = vol2Npcel *
                   fluidInterface.get_number_density(mfi, x, y, z, speciesID);
        if (q != 0) {
          double rand;
          Real u, v, w;
          double rand1 = randNum();
          double rand2 = randNum();
          double rand3 = randNum();
          double rand4 = randNum();

          if (fluidInterface.getUseAnisoP() &&
              (speciesID > 0 || fluidInterface.get_useElectronFluid())) {
            fluidInterface.set_particle_uth_aniso(mfi, x, y, z, &u, &v, &w,
                                                  rand1, rand2, rand3, rand4,
                                                  speciesID);
          } else {
            fluidInterface.set_particle_uth_iso(mfi, x, y, z, &u, &v, &w, rand1,
                                                rand2, rand3, rand4, speciesID);
          }

          // Increase the thermal velocity a little so that the variation of the
          // velocity space distribution is unbiased.
          u *= coefSD;
          v *= coefSD;
          w *= coefSD;

          u += fluidInterface.get_ux(mfi, x, y, z, speciesID);
          v += fluidInterface.get_uy(mfi, x, y, z, speciesID);
          w += fluidInterface.get_uz(mfi, x, y, z, speciesID);

          ParticleType p;
          if (ParticleType::the_next_id >= amrex::LastParticleID) {
            // id should not larger than LastParticleID. This is a bad solution,
            // since the ID becomes nonunique. --Yuxi
            p.id() = amrex::LastParticleID;
          } else {
            p.id() = ParticleType::NextID();
          }
          p.cpu() = ParallelDescriptor::MyProc();
          p.pos(ix_) = x; // + plo[ix_];
          p.pos(iy_) = y; // + plo[iy_];
          p.pos(iz_) = z; // + plo[iz_];
          p.rdata(iup_) = u;
          p.rdata(ivp_) = v;
          p.rdata(iwp_) = w;
          p.rdata(iqp_) = q;
          particles.push_back(p);
          // Print() << p << std::endl;
        }
      }
}

//==========================================================
void Particles::add_particles_domain(const FluidInterface& fluidInterface,
                                     const iMultiFab& cellStatus) {
  Timer funcTimer("Particles::add_particles");

  const int lev = 0;
  for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const auto& status = cellStatus[mfi].array();
    const Box& tile_box = mfi.validbox();
    const auto lo = amrex::lbound(tile_box);
    const auto hi = amrex::ubound(tile_box);

    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    int iMin = lo.x, jMin = lo.y, kMin = lo.z;

    for (int i = iMin; i <= iMax; ++i)
      for (int j = jMin; j <= jMax; ++j)
        for (int k = kMin; k <= kMax; ++k) {
          if (status(i, j, k) == iOnNew_) {
            add_particles_cell(mfi, fluidInterface, i, j, k);
          }
        }
  }
}

//==========================================================
void Particles::inject_particles_at_boundary(
    const FluidInterface& fluidInterface, const iMultiFab& cellStatus) {
  Timer funcTimer("Particles::inject_particles_at_boundary");

  // Only inject nGstInject layers.
  const int nGstInject = 1;

  const int lev = 0;

  for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const auto& status = cellStatus[mfi].array();
    const Box& bx = mfi.validbox();
    const IntVect lo = IntVect(bx.loVect());
    const IntVect hi = IntVect(bx.hiVect());
    // IntVect mid = (lo + hi) / 2;

    IntVect idxMin = lo, idxMax = hi;

    for (int iDim = 0; iDim < 3; iDim++) {
      if (!Geom(0).isPeriodic(iDim)) {
        idxMin[iDim] -= nGstInject;
        idxMax[iDim] += nGstInject;
      }
    }

    for (int i = idxMin[ix_]; i <= idxMax[ix_]; ++i)
      for (int j = idxMin[iy_]; j <= idxMax[iy_]; ++j)
        for (int k = idxMin[iz_]; k <= idxMax[iz_]; ++k) {
          if (do_inject_particles_for_this_cell(bx, status, i, j, k)) {
            add_particles_cell(mfi, fluidInterface, i, j, k);
          }
        }
  }
}

//==========================================================
void Particles::sum_to_center(amrex::MultiFab& netChargeMF,
                              amrex::UMultiFab<RealCMM>& centerMM,
                              bool doNetChargeOnly) {
  Timer funcTimer("Particles::sum_to_center");

  const auto& plo = Geom(0).ProbLo();

  const auto& dx = Geom(0).CellSize();
  const auto& invDx = Geom(0).InvCellSize();
  const Real invVol = invDx[ix_] * invDx[iy_] * invDx[iz_];

  const int lev = 0;
  for (ParticlesIter pti(*this, lev); pti.isValid(); ++pti) {
    Array4<Real> const& chargeArr = netChargeMF[pti].array();
    Array4<RealCMM> const& mmArr = centerMM[pti].array();

    const auto& particles = pti.GetArrayOfStructs();

    for (const auto& p : particles) {
      // Print() << "particle = " << p << std::endl;

      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      int loIdx[3];
      Real dShift[3];
      for (int i = 0; i < 3; i++) {
        // plo is the corner location => -0.5
        dShift[i] = (p.pos(i) - plo[i]) * invDx[i] - 0.5;
        loIdx[i] = fastfloor(dShift[i] + 10) - 10; // floor() is slow.
        dShift[i] = dShift[i] - loIdx[i];
      }
      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

      const Real cTmp = qp * invVol;
      for (int kk = 0; kk < 2; kk++)
        for (int jj = 0; jj < 2; jj++)
          for (int ii = 0; ii < 2; ii++) {
            chargeArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk) +=
                coef[ii][jj][kk] * cTmp;
          }

      if (!doNetChargeOnly) {
        Real weights_IIID[2][2][2][3];
        //----- Mass matrix calculation begin--------------
        const Real xi0 = dShift[ix_] * dx[ix_];
        const Real eta0 = dShift[iy_] * dx[iy_];
        const Real zeta0 = dShift[iz_] * dx[iz_];
        const Real xi1 = dx[ix_] - xi0;
        const Real eta1 = dx[iy_] - eta0;
        const Real zeta1 = dx[iz_] - zeta0;

        weights_IIID[1][1][1][ix_] = eta0 * zeta0 * invVol;
        weights_IIID[1][1][1][iy_] = xi0 * zeta0 * invVol;
        weights_IIID[1][1][1][iz_] = xi0 * eta0 * invVol;

        // xi0*eta0*zeta1*invVOL;
        weights_IIID[1][1][0][ix_] = eta0 * zeta1 * invVol;
        weights_IIID[1][1][0][iy_] = xi0 * zeta1 * invVol;
        weights_IIID[1][1][0][iz_] = -xi0 * eta0 * invVol;

        // xi0*eta1*zeta0*invVOL;
        weights_IIID[1][0][1][ix_] = eta1 * zeta0 * invVol;
        weights_IIID[1][0][1][iy_] = -xi0 * zeta0 * invVol;
        weights_IIID[1][0][1][iz_] = xi0 * eta1 * invVol;

        // xi0*eta1*zeta1*invVOL;
        weights_IIID[1][0][0][ix_] = eta1 * zeta1 * invVol;
        weights_IIID[1][0][0][iy_] = -xi0 * zeta1 * invVol;
        weights_IIID[1][0][0][iz_] = -xi0 * eta1 * invVol;

        // xi1*eta0*zeta0*invVOL;
        weights_IIID[0][1][1][ix_] = -eta0 * zeta0 * invVol;
        weights_IIID[0][1][1][iy_] = xi1 * zeta0 * invVol;
        weights_IIID[0][1][1][iz_] = xi1 * eta0 * invVol;

        // xi1*eta0*zeta1*invVOL;
        weights_IIID[0][1][0][ix_] = -eta0 * zeta1 * invVol;
        weights_IIID[0][1][0][iy_] = xi1 * zeta1 * invVol;
        weights_IIID[0][1][0][iz_] = -xi1 * eta0 * invVol;

        // xi1*eta1*zeta0*invVOL;
        weights_IIID[0][0][1][ix_] = -eta1 * zeta0 * invVol;
        weights_IIID[0][0][1][iy_] = -xi1 * zeta0 * invVol;
        weights_IIID[0][0][1][iz_] = xi1 * eta1 * invVol;

        // xi1*eta1*zeta1*invVOL;
        weights_IIID[0][0][0][ix_] = -eta1 * zeta1 * invVol;
        weights_IIID[0][0][0][iy_] = -xi1 * zeta1 * invVol;
        weights_IIID[0][0][0][iz_] = -xi1 * eta1 * invVol;

        const int iMin = loIdx[ix_];
        const int jMin = loIdx[iy_];
        const int kMin = loIdx[iz_];
        const int iMax = iMin + 1;
        const int jMax = jMin + 1;
        const int kMax = kMin + 1;

        const Real coef = fabs(qp) * invVol;
        Real wg_D[3];
        for (int k1 = kMin; k1 <= kMax; k1++)
          for (int j1 = jMin; j1 <= jMax; j1++)
            for (int i1 = iMin; i1 <= iMax; i1++) {

              for (int iDim = 0; iDim < 3; iDim++) {
                wg_D[iDim] =
                    coef * weights_IIID[i1 - iMin][j1 - jMin][k1 - kMin][iDim];
              }

              Real* const data = mmArr(i1, j1, k1).data;
              // Real weights[27] = { 0 };
              for (int i2 = iMin; i2 <= iMax; i2++) {
                int ip = i2 - i1 + 1;
                const int gp0 = ip * 9;
                for (int j2 = jMin; j2 <= jMax; j2++) {
                  int jp = j2 - j1 + 1;
                  const int gp1 = gp0 + jp * 3;
                  for (int k2 = kMin; k2 <= kMax; k2++) {
                    const Real(&wg1_D)[3] =
                        weights_IIID[i2 - iMin][j2 - jMin][k2 - kMin];

                    // const int kp = k2 - k1 + 1;
                    const int gp = gp1 + k2 - k1 + 1;

                    data[gp] += wg_D[ix_] * wg1_D[ix_] +
                                wg_D[iy_] * wg1_D[iy_] + wg_D[iz_] * wg1_D[iz_];
                    ;
                  }
                }
              }
            }
      } // if doChargeOnly

    } // for p
  }
}

//==========================================================
PartInfo Particles::sum_moments(MultiFab& momentsMF, UMultiFab<RealMM>& nodeMM,
                                MultiFab& nodeBMF, Real dt) {
  Timer funcTimer("Particles::sum_moments");
  const auto& plo = Geom(0).ProbLo();

  momentsMF.setVal(0.0);

  Real qdto2mc = charge / mass * 0.5 * dt;

  const auto& invDx = Geom(0).InvCellSize();
  const Real invVol = invDx[ix_] * invDx[iy_] * invDx[iz_];

  PartInfo pinfo;
  const int lev = 0;
  for (ParticlesIter pti(*this, lev); pti.isValid(); ++pti) {
    Array4<Real> const& momentsArr = momentsMF[pti].array();
    Array4<Real const> const& nodeBArr = nodeBMF[pti].array();
    Array4<RealMM> const& mmArr = nodeMM[pti].array();

    const auto& particles = pti.GetArrayOfStructs();

    for (const auto& p : particles) {

      // Print()<<"p = "<<p<<std::endl;
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      int loIdx[3];
      Real dShift[3];
      for (int i = 0; i < 3; i++) {
        dShift[i] = (p.pos(i) - plo[i]) * invDx[i];
        loIdx[i] = fastfloor(dShift[i]);
        dShift[i] = dShift[i] - loIdx[i];
      }

      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

      //----- Mass matrix calculation begin--------------
      Real Bxl = 0, Byl = 0, Bzl = 0; // should be bp[3];

      for (int kk = 0; kk < 2; kk++)
        for (int jj = 0; jj < 2; jj++)
          for (int ii = 0; ii < 2; ii++) {
            Bxl += nodeBArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk,
                            ix_) *
                   coef[ii][jj][kk];
            Byl += nodeBArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk,
                            iy_) *
                   coef[ii][jj][kk];
            Bzl += nodeBArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk,
                            iz_) *
                   coef[ii][jj][kk];
          }

      const Real Omx = qdto2mc * Bxl;
      const Real Omy = qdto2mc * Byl;
      const Real Omz = qdto2mc * Bzl;

      // end interpolation
      const Real omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
      const Real denom = 1.0 / (1.0 + omsq);
      const Real udotOm = up * Omx + vp * Omy + wp * Omz;

      const Real c0 = denom * invVol * qp * qdto2mc;
      Real alpha[9];
      alpha[0] = (1 + Omx * Omx) * c0;
      alpha[1] = (Omz + Omx * Omy) * c0;
      alpha[2] = (-Omy + Omx * Omz) * c0;
      alpha[3] = (-Omz + Omx * Omy) * c0;
      alpha[4] = (1 + Omy * Omy) * c0;
      alpha[5] = (Omx + Omy * Omz) * c0;
      alpha[6] = (Omy + Omx * Omz) * c0;
      alpha[7] = (-Omx + Omy * Omz) * c0;
      alpha[8] = (1 + Omz * Omz) * c0;

      // int count = 0;
      // int ip, jp, kp;         // Indexes for g'
      // double wg, wgp, weight; // W_pg, W_pg'

      const int iMin = loIdx[ix_];
      const int jMin = loIdx[iy_];
      const int kMin = loIdx[iz_];
      const int iMax = iMin + 2;
      const int jMax = jMin + 2;
      const int kMax = kMin + 2;

      for (int k1 = kMin; k1 < kMax; k1++)
        for (int j1 = jMin; j1 < jMax; j1++)
          for (int i1 = iMin; i1 < iMax; i1++) {
            const Real wg = coef[i1 - iMin][j1 - jMin][k1 - kMin];
            Real* const data0 = mmArr(i1, j1, k1).data;
            for (int k2 = kMin; k2 < kMax; k2++) {
              const int kp = k2 - k1 + 1;
              if (kp > 0) {
                for (int j2 = jMin; j2 < jMax; j2++) {
                  const int jp = j2 - j1 + 1;
                  for (int i2 = iMin; i2 < iMax; i2++) {
                    const Real weight =
                        wg * coef[i2 - iMin][j2 - jMin][k2 - kMin];
                    const int idx0 = kp * 81 + jp * 27 + (i2 - i1 + 1) * 9;

                    Real* const data = &(data0[idx0]);
                    for (int idx = 0; idx < 9; idx++) {
                      data[idx] += alpha[idx] * weight;
                    }
                  } // k2

                } // j2
              }   // if (ip > 0)
            }     // i2
          }       // k1

      //----- Mass matrix calculation end--------------

      //-------Moments begin---------
      Real pMoments[nMoments];

      pMoments[iNum_] = 1;
      pMoments[iRho_] = qp;

      {
        const Real mx = qp * up;
        const Real my = qp * vp;
        const Real mz = qp * wp;
        pMoments[iMx_] = mx;
        pMoments[iMy_] = my;
        pMoments[iMz_] = mz;

        pMoments[iPxx_] = mx * up;
        pMoments[iPyy_] = my * vp;
        pMoments[iPzz_] = mz * wp;

        pMoments[iPxy_] = mx * vp;
        pMoments[iPxz_] = mx * wp;
        pMoments[iPyz_] = my * wp;
      }

      {
        const Real coef1 = denom * qp;
        pMoments[iJhx_] = (up + (vp * Omz - wp * Omy + udotOm * Omx)) * coef1;
        pMoments[iJhy_] = (vp + (wp * Omx - up * Omz + udotOm * Omy)) * coef1;
        pMoments[iJhz_] = (wp + (up * Omy - vp * Omx + udotOm * Omz)) * coef1;
      }

      for (int iVar = 0; iVar < nMoments; iVar++)
        for (int kk = 0; kk < 2; kk++)
          for (int jj = 0; jj < 2; jj++)
            for (int ii = 0; ii < 2; ii++) {
              momentsArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk,
                         iVar) += coef[ii][jj][kk] * pMoments[iVar];
            }

      //-------Moments end---------

      pinfo.energy += qp * (up * up + vp * vp + wp * wp);
    } // for p
  }

  for (MFIter mfi(nodeMM); mfi.isValid(); ++mfi) {
    // Finalize the mass matrix calculation.
    const Box box = mfi.validbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Array4<RealMM> const& mmArr = nodeMM[mfi].array();

    // We only need the mass matrix on the physical nodes. But the first layer
    // of the ghost nodes may contributes to the physical nodes below (ghost
    // node constributes as a sender). So, we need the '-1' and '+1' staff.
    const int iMin = lo.x - 1, jMin = lo.y - 1, kMin = lo.z - 1;
    const int iMax = hi.x + 1, jMax = hi.y + 1, kMax = hi.z + 1;

    int gps, gpr; // gp_send, gp_receive
    for (int k1 = kMin; k1 <= kMax; k1++)
      for (int j1 = jMin; j1 <= jMax; j1++)
        for (int i1 = iMin; i1 <= iMax; i1++) {
          const int kp = 2, kpr = 0;
          const int kr = k1 + kp - 1;
          const Real* const datas0 = mmArr(i1, j1, k1).data;
          for (int jp = 0; jp < 3; jp++) {
            const int jr = j1 + jp - 1;
            const int jpr = 2 - jp;
            for (int ip = 0; ip < 3; ip++) {
              const int ir = i1 + ip - 1;
              const int ipr = 2 - ip;
              gpr = jpr * 3 + ipr;
              gps = 18 + jp * 3 + ip; // gps = kp*9+jp*3+kp

              Real* const datar = &(mmArr(ir, jr, kr).data[gpr * 9]);
              const Real* const datas = &(datas0[gps * 9]);
              for (int idx = 0; idx < 9; idx++) {
                datar[idx] = datas[idx];
              } // idx
            }   // kp
          }     // jp
        }       // k1
  }

  // Exclude the number density.
  momentsMF.mult(invVol, 0, nMoments - 1, momentsMF.nGrow());

  momentsMF.SumBoundary(Geom(0).periodicity());

  // FillBoundary seems unnecessary. --Yuxi
  momentsMF.FillBoundary(Geom(0).periodicity());

  // This function should be called before 'convert_to_fluid_moments'
  pinfo.uMax = calc_max_thermal_velocity(momentsMF);

  convert_to_fluid_moments(momentsMF);

  // Calculate the total particle energy--------------
  pinfo.energy *= 0.5 / get_qom();
  ParallelDescriptor::ReduceRealSum(pinfo.energy,
                                    ParallelDescriptor::IOProcessorNumber());
  if (!ParallelDescriptor::IOProcessor())
    pinfo.energy = 0;
  //------------------------------------------------------

  return pinfo;
}

//==========================================================
Real Particles::calc_max_thermal_velocity(MultiFab& momentsMF) {

  Real uthMax = 0;
  const Real c1over3 = 1. / 3;
  for (MFIter mfi(momentsMF); mfi.isValid(); ++mfi) {
    FArrayBox& fab = momentsMF[mfi];
    const Box& box = mfi.validbox();
    const Array4<Real>& arr = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    // Do not calculate for edges twice.
    for (int k = lo.z; k <= hi.z - 1; ++k)
      for (int j = lo.y; j <= hi.y - 1; ++j)
        for (int i = lo.x; i <= hi.x - 1; ++i) {
          Real rho = arr(i, j, k, iRho_);
          if (rho == 0)
            continue;

          Real p = (arr(i, j, k, iPxx_) + arr(i, j, k, iPyy_) +
                    arr(i, j, k, iPzz_)) *
                   c1over3;

          Real uth = sqrt(p / rho);
          if (uth > uthMax)
            uthMax = uth;
        }
  }

  ParallelDescriptor::ReduceRealMax(uthMax,
                                    ParallelDescriptor::IOProcessorNumber());

  if (!ParallelDescriptor::IOProcessor())
    uthMax = 0;

  return uthMax;
}

//==========================================================
void Particles::convert_to_fluid_moments(MultiFab& momentsMF) {
  MultiFab tmpMF(momentsMF, make_alias, iRho_, iPyz_ - iRho_ + 1);
  tmpMF.mult(1.0 / get_qom(), tmpMF.nGrow());

  for (MFIter mfi(momentsMF); mfi.isValid(); ++mfi) {
    FArrayBox& fab = momentsMF[mfi];
    const Box& box = mfi.fabbox();
    const Array4<Real>& arr = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          const Real rho = arr(i, j, k, iRho_);
          if (rho > 1e-99) {
            const Real ux = arr(i, j, k, iUx_) / rho;
            const Real uy = arr(i, j, k, iUy_) / rho;
            const Real uz = arr(i, j, k, iUz_) / rho;
            arr(i, j, k, iPxx_) = arr(i, j, k, iPxx_) - rho * ux * ux;
            arr(i, j, k, iPyy_) = arr(i, j, k, iPyy_) - rho * uy * uy;
            arr(i, j, k, iPzz_) = arr(i, j, k, iPzz_) - rho * uz * uz;

            arr(i, j, k, iPxy_) = arr(i, j, k, iPxy_) - rho * ux * uy;
            arr(i, j, k, iPxz_) = arr(i, j, k, iPxz_) - rho * ux * uz;
            arr(i, j, k, iPyz_) = arr(i, j, k, iPyz_) - rho * uy * uz;
          }
        }
  }
}

//==========================================================
void Particles::mover(const amrex::MultiFab& nodeEMF,
                      const amrex::MultiFab& nodeBMF, amrex::Real dt) {
  Timer funcTimer("Particles::mover");

  const auto& plo = Geom(0).ProbLo();

  const Real dtLoc = dt;
  const Real qdto2mc = charge / mass * 0.5 * dt;

  const auto& invDx = Geom(0).InvCellSize();

  const int lev = 0;
  for (ParticlesIter pti(*this, lev); pti.isValid(); ++pti) {
    const Array4<Real const>& nodeEArr = nodeEMF[pti].array();
    const Array4<Real const>& nodeBArr = nodeBMF[pti].array();

    auto& particles = pti.GetArrayOfStructs();
    for (auto& p : particles) {
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);
      const Real xp = p.pos(ix_);
      const Real yp = p.pos(iy_);
      const Real zp = p.pos(iz_);

      //-----calculate interpolate coef begin-------------
      int loIdx[3];
      Real dShift[3];
      for (int i = 0; i < 3; i++) {
        dShift[i] = (p.pos(i) - plo[i]) * invDx[i];
        loIdx[i] = fastfloor(dShift[i]);
        dShift[i] = dShift[i] - loIdx[i];
      }

      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

      Real Bxl = 0, Byl = 0, Bzl = 0; // should be bp[3];
      Real Exl = 0, Eyl = 0, Ezl = 0;
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++) {
            const int iNodeX = loIdx[ix_] + ii;
            const int iNodeY = loIdx[iy_] + jj;
            const int iNodeZ = loIdx[iz_] + kk;
            const Real& c0 = coef[ii][jj][kk];
            Bxl += nodeBArr(iNodeX, iNodeY, iNodeZ, ix_) * c0;
            Byl += nodeBArr(iNodeX, iNodeY, iNodeZ, iy_) * c0;
            Bzl += nodeBArr(iNodeX, iNodeY, iNodeZ, iz_) * c0;

            Exl += nodeEArr(iNodeX, iNodeY, iNodeZ, ix_) * c0;
            Eyl += nodeEArr(iNodeX, iNodeY, iNodeZ, iy_) * c0;
            Ezl += nodeEArr(iNodeX, iNodeY, iNodeZ, iz_) * c0;
          }

      const double Omx = qdto2mc * Bxl;
      const double Omy = qdto2mc * Byl;
      const double Omz = qdto2mc * Bzl;

      // end interpolation
      const Real omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
      const Real denom = 1.0 / (1.0 + omsq);
      // solve the position equation
      const Real ut = up + qdto2mc * Exl;
      const Real vt = vp + qdto2mc * Eyl;
      const Real wt = wp + qdto2mc * Ezl;
      // const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
      const Real udotOm = ut * Omx + vt * Omy + wt * Omz;
      // solve the velocity equation
      const Real uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
      const Real vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
      const Real wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;

      const double unp1 = 2.0 * uavg - up;
      const double vnp1 = 2.0 * vavg - vp;
      const double wnp1 = 2.0 * wavg - wp;

      p.rdata(ix_) = unp1;
      p.rdata(iy_) = vnp1;
      p.rdata(iz_) = wnp1;

      p.pos(ix_) = xp + unp1 * dtLoc;
      p.pos(iy_) = yp + vnp1 * dtLoc;
      p.pos(iz_) = zp + wnp1 * dtLoc;

      // Mark for deletion
      if (is_outside_ba(p)) {
        p.id() = -1;
      }
    } // for p
  }   // for pti

  // This function distributes particles to proper processors and apply
  // periodic boundary conditions if needed.
  Redistribute();
}

//==========================================================
void Particles::divE_correct_position(const amrex::MultiFab& phiMF) {
  Timer funcTimer("Particles::divE_correct_position");

  const auto& plo = Geom(0).ProbLo();

  const auto& dx = Geom(0).CellSize();
  const auto& invDx = Geom(0).InvCellSize();
  const Real invVol = invDx[ix_] * invDx[iy_] * invDx[iz_];

  const Real coef = charge / fabs(charge);
  const Real epsLimit = 0.1;
  Real epsMax = 0;

  const int lev = 0;
  for (ParticlesIter pti(*this, lev); pti.isValid(); ++pti) {
    Array4<Real const> const& phiArr = phiMF[pti].array();

    auto& particles = pti.GetArrayOfStructs();

    for (auto& p : particles) {

      const Real qp = p.rdata(iqp_);

      int loIdx[3];
      Real dShift[3];
      for (int i = 0; i < 3; i++) {
        // plo is the corner location => -0.5
        dShift[i] = (p.pos(i) - plo[i]) * invDx[i] - 0.5;
        loIdx[i] = fastfloor(dShift[i]);
        dShift[i] = dShift[i] - loIdx[i];
      }

      {
        Real weights_IIID[2][2][2][3];
        //----- Mass matrix calculation begin--------------
        const Real xi0 = dShift[ix_] * dx[ix_];
        const Real eta0 = dShift[iy_] * dx[iy_];
        const Real zeta0 = dShift[iz_] * dx[iz_];
        const Real xi1 = dx[ix_] - xi0;
        const Real eta1 = dx[iy_] - eta0;
        const Real zeta1 = dx[iz_] - zeta0;

        const Real zeta02Vol = zeta0 * invVol;
        const Real zeta12Vol = zeta1 * invVol;
        const Real eta02Vol = eta0 * invVol;
        const Real eta12Vol = eta1 * invVol;

        weights_IIID[1][1][1][ix_] = eta0 * zeta02Vol;
        weights_IIID[1][1][1][iy_] = xi0 * zeta02Vol;
        weights_IIID[1][1][1][iz_] = xi0 * eta02Vol;

        // xi0*eta0*zeta1*invVOL;
        weights_IIID[1][1][0][ix_] = eta0 * zeta12Vol;
        weights_IIID[1][1][0][iy_] = xi0 * zeta12Vol;
        weights_IIID[1][1][0][iz_] = -xi0 * eta02Vol;

        // xi0*eta1*zeta0*invVOL;
        weights_IIID[1][0][1][ix_] = eta1 * zeta02Vol;
        weights_IIID[1][0][1][iy_] = -xi0 * zeta02Vol;
        weights_IIID[1][0][1][iz_] = xi0 * eta12Vol;

        // xi0*eta1*zeta1*invVOL;
        weights_IIID[1][0][0][ix_] = eta1 * zeta12Vol;
        weights_IIID[1][0][0][iy_] = -xi0 * zeta12Vol;
        weights_IIID[1][0][0][iz_] = -xi0 * eta12Vol;

        // xi1*eta0*zeta0*invVOL;
        weights_IIID[0][1][1][ix_] = -eta0 * zeta02Vol;
        weights_IIID[0][1][1][iy_] = xi1 * zeta02Vol;
        weights_IIID[0][1][1][iz_] = xi1 * eta02Vol;

        // xi1*eta0*zeta1*invVOL;
        weights_IIID[0][1][0][ix_] = -eta0 * zeta12Vol;
        weights_IIID[0][1][0][iy_] = xi1 * zeta12Vol;
        weights_IIID[0][1][0][iz_] = -xi1 * eta02Vol;

        // xi1*eta1*zeta0*invVOL;
        weights_IIID[0][0][1][ix_] = -eta1 * zeta02Vol;
        weights_IIID[0][0][1][iy_] = -xi1 * zeta02Vol;
        weights_IIID[0][0][1][iz_] = xi1 * eta12Vol;

        // xi1*eta1*zeta1*invVOL;
        weights_IIID[0][0][0][ix_] = -eta1 * zeta12Vol;
        weights_IIID[0][0][0][iy_] = -xi1 * zeta12Vol;
        weights_IIID[0][0][0][iz_] = -xi1 * eta12Vol;

        const int iMin = loIdx[ix_];
        const int jMin = loIdx[iy_];
        const int kMin = loIdx[iz_];

        Real eps_D[3] = { 0, 0, 0 };

        for (int k = 0; k < 2; k++)
          for (int j = 0; j < 2; j++)
            for (int i = 0; i < 2; i++) {
              const Real coef = phiArr(iMin + i, jMin + j, kMin + k);
              for (int iDim = 0; iDim < 3; iDim++) {
                eps_D[iDim] += coef * weights_IIID[i][j][k][iDim];
              }
            }

        for (int iDim = 0; iDim < 3; iDim++)
          eps_D[iDim] *= coef * fourPI;

        if (fabs(eps_D[ix_] * invDx[ix_]) > epsLimit ||
            fabs(eps_D[iy_] * invDx[iy_]) > epsLimit ||
            fabs(eps_D[iz_] * invDx[iz_]) > epsLimit) {
          // If eps_D is too large, the underlying assumption of the particle
          // correction method will be not valid. Comparing each exp_D
          // component instead of the length dl saves the computational time.
          const Real dl = sqrt(pow(eps_D[ix_], 2) + pow(eps_D[iy_], 2) +
                               pow(eps_D[iz_], 2));
          const Real ratio = epsLimit * dx[ix_] / dl;
          for (int iDim = 0; iDim < 3; iDim++)
            eps_D[iDim] *= ratio;
        }

        for (int iDim = 0; iDim < 3; iDim++) {
          if (fabs(eps_D[iDim] * invDx[iDim]) > epsMax)
            epsMax = fabs(eps_D[iDim] * invDx[iDim]);

          p.pos(iDim) += eps_D[iDim];

          if (is_outside_ba(p)) {
            p.id() = -1;
          }
        }
      }

    } // for p
  }
}

//==========================================================
void Particles::split_particles(Real limit) {
  Timer funcTimer("Particles::split_particles");
  const int nPartGoal =
      nPartPerCell[ix_] * nPartPerCell[iy_] * nPartPerCell[iz_] * limit;

  IntVect iv = { 1, 1, 1 };
  if (!(do_tiling && tile_size == iv))
    return;

  const int lev = 0;
  Real dl = 0.5 * Geom(0).CellSize()[ix_] / nPartPerCell.max();

  for (ParticlesIter pti(*this, lev); pti.isValid(); ++pti) {

    auto& particles = pti.GetArrayOfStructs();

    const int nPartOrig = particles.size();

    const int nNew =
        nPartGoal - nPartOrig > nPartOrig ? nPartOrig : nPartGoal - nPartOrig;

    if (nNew <= 0)
      continue;

    // Find the 'heaviest' nNew particles by sorting the weight (charge).-----

    // Sort the particles by x first to make sure the results are the same for
    // different number of processors
    std::sort(particles.begin(), particles.end(),
              [ix_ = ix_](const auto& pl, const auto& pr) {
                return pl.rdata(ix_) > pr.rdata(ix_);
              });

    std::sort(particles.begin(), particles.end(),
              [ix_ = ix_](const auto& pl, const auto& pr) {
                return fabs(pl.rdata(iqp_)) > fabs(pr.rdata(iqp_));
              });
    //----------------------------------------------------------------

    const auto lo = lbound(pti.tilebox());
    const auto hi = ubound(pti.tilebox());

    const Real xMin =
                   Geom(0).LoEdge(lo.x, ix_) + Geom(0).CellSize()[ix_] * 1e-10,
               xMax =
                   Geom(0).HiEdge(hi.x, ix_) - Geom(0).CellSize()[ix_] * 1e-10;

    const Real yMin =
                   Geom(0).LoEdge(lo.y, iy_) + Geom(0).CellSize()[iy_] * 1e-10,
               yMax =
                   Geom(0).HiEdge(hi.y, iy_) - Geom(0).CellSize()[iy_] * 1e-10;

    const Real zMin =
                   Geom(0).LoEdge(lo.z, iz_) + Geom(0).CellSize()[iz_] * 1e-10,
               zMax =
                   Geom(0).HiEdge(hi.z, iz_) - Geom(0).CellSize()[iz_] * 1e-10;

    for (int ip = 0; ip < nNew; ip++) {
      auto& p = particles[ip];
      Real qp1 = p.rdata(iqp_);
      Real xp1 = p.pos(ix_);
      Real yp1 = p.pos(iy_);
      Real zp1 = p.pos(iz_);
      Real up1 = p.rdata(iup_);
      Real vp1 = p.rdata(ivp_);
      Real wp1 = p.rdata(iwp_);

      const Real coef = dl / sqrt(up1 * up1 + vp1 * vp1 + wp1 * wp1);
      const Real dpx = coef * up1;
      const Real dpy = coef * vp1;
      const Real dpz = coef * wp1;

      Real xp2 = xp1 + dpx;
      Real yp2 = yp1 + dpy;
      Real zp2 = zp1 + dpz;

      xp1 -= dpx;
      yp1 -= dpy;
      zp1 -= dpz;

      xp1 = bound(xp1, xMin, xMax);
      yp1 = bound(yp1, yMin, yMax);
      zp1 = bound(zp1, zMin, zMax);

      xp2 = bound(xp2, xMin, xMax);
      yp2 = bound(yp2, yMin, yMax);
      zp2 = bound(zp2, zMin, zMax);

      p.rdata(iqp_) = (qp1 / 2);
      p.pos(ix_) = xp1;
      p.pos(iy_) = yp1;
      p.pos(iz_) = zp1;

      ParticleType pnew;
      if (ParticleType::the_next_id >= amrex::LastParticleID) {
        // id should not larger than LastParticleID. This is a bad solution,
        // since the ID becomes nonunique. --Yuxi
        pnew.id() = amrex::LastParticleID;
      } else {
        pnew.id() = ParticleType::NextID();
      }

      pnew.cpu() = ParallelDescriptor::MyProc();
      pnew.pos(ix_) = xp2;
      pnew.pos(iy_) = yp2;
      pnew.pos(iz_) = zp2;
      pnew.rdata(iup_) = up1;
      pnew.rdata(ivp_) = vp1;
      pnew.rdata(iwp_) = wp1;
      pnew.rdata(iqp_) = qp1 / 2;
      particles.push_back(pnew);
    }
  }
}

//==========================================================
void Particles::combine_particles(Real limit) {
  Timer funcTimer("Particles::combine_particles");
  IntVect iv = { 1, 1, 1 };
  if (!(do_tiling && tile_size == iv))
    return;

  const int nPartGoal =
      nPartPerCell[ix_] * nPartPerCell[iy_] * nPartPerCell[iz_] * limit;

  const int lev = 0;

  for (ParticlesIter pti(*this, lev); pti.isValid(); ++pti) {

    auto& particles = pti.GetArrayOfStructs();

    const int nPartOrig = particles.size();

    if (nPartOrig <= nPartGoal)
      continue;
    const int nCombineGoal = nPartOrig - nPartGoal;

    // Phase space cell number in one direction.
    // The const 0.8 is choosen by experience.
    const int nDim = 3;
    const int nCell = 0.8 * pow(nPartGoal, 1. / nDim);
    if (nCell < 1)
      continue;

    //----------------------------------------------------------------
    // Estimate the bulk velocity and thermal velocity.
    Real uAv = 0, vAv = 0, wAv = 0;
    for (int pid = 0; pid < nPartOrig; pid++) {
      auto& pcl = particles[pid];
      uAv += pcl.rdata(ix_);
      vAv += pcl.rdata(iy_);
      wAv += pcl.rdata(iz_);
    }
    uAv /= nPartOrig;
    vAv /= nPartOrig;
    wAv /= nPartOrig;

    Real thVel = 0, thVel2 = 0;
    for (int pid = 0; pid < nPartOrig; pid++) {
      auto& pcl = particles[pid];
      thVel2 += pow(pcl.rdata(ix_) - uAv, 2) + pow(pcl.rdata(iy_) - vAv, 2) +
                pow(pcl.rdata(iz_) - wAv, 2);
    }

    thVel2 /= nPartOrig;
    thVel = sqrt(thVel2);

    //----------------------------------------------------------------

    // Storing the particle indices in the corresponding phase space cell.
    Vector<int> phasePartIdx_III[nCell][nCell][nCell];

    const int u_ = 0, v_ = 1, w_ = 2;
    Real velMin_D[nDim], velMax_D[nDim], inv_dVel_D[nDim];
    const Real r0 = 1.0;
    velMin_D[u_] = -r0 * thVel + uAv;
    velMax_D[u_] = r0 * thVel + uAv;
    velMin_D[v_] = -r0 * thVel + vAv;
    velMax_D[v_] = r0 * thVel + vAv;
    velMin_D[w_] = -r0 * thVel + wAv;
    velMax_D[w_] = r0 * thVel + wAv;
    for (int iDim = 0; iDim < nDim; iDim++) {
      inv_dVel_D[iDim] = nCell / (velMax_D[iDim] - velMin_D[iDim]);
    }

    Real vel_D[nDim];
    int iCell_D[nDim];
    for (int pid = 0; pid < nPartOrig; pid++) {
      auto& pcl = particles[pid];
      vel_D[u_] = pcl.rdata(ix_);
      vel_D[v_] = pcl.rdata(iy_);
      vel_D[w_] = pcl.rdata(iz_);

      bool isOutside = false;
      for (int iDim = 0; iDim < nDim; iDim++) {
        if (vel_D[iDim] < velMin_D[iDim] || vel_D[iDim] > velMax_D[iDim])
          isOutside = true;
      }
      if (isOutside)
        continue;

      for (int iDim = 0; iDim < nDim; iDim++) {
        iCell_D[iDim] =
            fastfloor((vel_D[iDim] - velMin_D[iDim]) * inv_dVel_D[iDim]);
      }

      phasePartIdx_III[iCell_D[u_]][iCell_D[v_]][iCell_D[w_]].push_back(pid);
    }

    // Sorting the particles indexes so that the results change with different
    // number of processors.
    for (int iCell = 0; iCell < nCell; iCell++)
      for (int jCell = 0; jCell < nCell; jCell++)
        for (int kCell = 0; kCell < nCell; kCell++) {
          std::sort(phasePartIdx_III[iCell][jCell][kCell].begin(),
                    phasePartIdx_III[iCell][jCell][kCell].end(),
                    [& particles = particles, ix_ = ix_](const int& idl,
                                                         const int& idr) {
                      return particles[idl].rdata(ix_) >
                             particles[idr].rdata(ix_);
                    });
        }

    const int nPartCombine = 6;
    int iCount = 0;
    int nAvailableCombine = 0;
    for (int iu = 0; iu < nCell; iu++)
      for (int iv = 0; iv < nCell; iv++)
        for (int iw = 0; iw < nCell; iw++) {
          iCount += phasePartIdx_III[iu][iv][iw].size();
          nAvailableCombine +=
              fastfloor(phasePartIdx_III[iu][iv][iw].size() / nPartCombine);
        }

    Real ratioCombine;
    if (nAvailableCombine < nCombineGoal) {
      ratioCombine = 1;
    } else {
      ratioCombine = Real(nCombineGoal) / nAvailableCombine;
    }

    const auto lo = lbound(pti.tilebox());
    const auto hi = ubound(pti.tilebox());

    const Real xMin = Geom(0).LoEdge(lo.x, ix_),
               xMax = Geom(0).HiEdge(hi.x, ix_);

    const Real yMin = Geom(0).LoEdge(lo.y, iy_),
               yMax = Geom(0).HiEdge(hi.y, iy_);

    const Real zMin = Geom(0).LoEdge(lo.z, iz_),
               zMax = Geom(0).HiEdge(hi.z, iz_);

    const Real dx = Geom(0).CellSize(ix_);
    const Real dy = Geom(0).CellSize(iy_);
    const Real dz = Geom(0).CellSize(iz_);

    const Real invdl2 = 1.0 / (dx * dx + dy * dy + dz * dz);
    const Real invthVel2 = 1.0 / thVel2;
    for (int iu = 0; iu < nCell; iu++)
      for (int iv = 0; iv < nCell; iv++)
        for (int iw = 0; iw < nCell; iw++) {
          int nCombineCell =
              fastfloor(ratioCombine * phasePartIdx_III[iu][iv][iw].size() /
                      nPartCombine);
          for (int iCombine = 0; iCombine < nCombineCell; iCombine++) {
            /*
                Delete 1 particle out of 6 particles:
                1) Choose two particles that are closest to each other.
                2) Delete the lighter one.
                3) Distribute its weights to another 5 particles to conserve
                   mass, momentum and energy.
                4) Adjust the location of the particles to conserve
                   mass center (optional. Turned off by default.)
             */

            int idx_I[nPartCombine];
            for (int ip = 0; ip < nPartCombine; ip++) {
              // Pop three particle indices.
              idx_I[ip] = phasePartIdx_III[iu][iv][iw].back();
              phasePartIdx_III[iu][iv][iw].pop_back();
            }

            // Find the pairs close to each other in phase space.
            int pair1 = 0, pair2 = 0;
            Real dis2Min = 2;
            for (int ip1 = 0; ip1 < nPartCombine - 1; ip1++)
              for (int ip2 = ip1 + 1; ip2 < nPartCombine; ip2++) {
                const Real dup = particles[idx_I[ip1]].rdata(ix_) -
                                 particles[idx_I[ip2]].rdata(ix_);
                const Real dvp = particles[idx_I[ip1]].rdata(iy_) -
                                 particles[idx_I[ip2]].rdata(iy_);
                const Real dwp = particles[idx_I[ip1]].rdata(iz_) -
                                 particles[idx_I[ip2]].rdata(iz_);
                const Real dxp = particles[idx_I[ip1]].pos(ix_) -
                                 particles[idx_I[ip2]].pos(ix_);
                const Real dyp = particles[idx_I[ip1]].pos(iy_) -
                                 particles[idx_I[ip2]].pos(iy_);
                const Real dzp = particles[idx_I[ip1]].pos(iz_) -
                                 particles[idx_I[ip2]].pos(iz_);

                const Real dis2 =
                    (dup * dup + dvp * dvp + dwp * dwp) * invthVel2 +
                    (dxp * dxp + dyp * dyp + dzp * dzp) * invdl2;

                if (dis2 < dis2Min) {
                  dis2Min = dis2;
                  pair1 = ip1;
                  pair2 = ip2;
                }
              }

            // Delete the lighter one.
            int iPartDel = pair1, iPartKeep = pair2;
            if (fabs(particles[idx_I[pair1]].rdata(iqp_)) >
                fabs(particles[idx_I[pair2]].rdata(iqp_))) {
              iPartDel = pair2;
              iPartKeep = pair1;
            }

            if (iPartDel != nPartCombine - 1) {
              const int idxTmp = idx_I[iPartDel];
              idx_I[iPartDel] = idx_I[nPartCombine - 1];
              idx_I[nPartCombine - 1] = idxTmp;
            }

            //---------------Solve the new particle weights---------------------
            const int nVar = 5;
            const int iq_ = 0, iu_ = 1, iv_ = 2, iw_ = 3, ie_ = 4;
            Real a[nVar][nVar + 1];
            Real x[nVar];
            for (int i = 0; i < nVar; i++) {
              x[i] = 0;
              for (int j = 0; j < nVar + 1; j++) {
                a[i][j] = 0;
              }
            }

            for (int ip = 0; ip < nPartCombine; ip++) {
              const Real qp = particles[idx_I[ip]].rdata(iqp_);
              const Real up = particles[idx_I[ip]].rdata(ix_);
              const Real vp = particles[idx_I[ip]].rdata(iy_);
              const Real wp = particles[idx_I[ip]].rdata(iz_);
              const Real v2 = (pow(up, 2) + pow(vp, 2) + pow(wp, 2));

              if (ip < nVar) {
                a[iq_][ip] = 1;
                a[iu_][ip] = up;
                a[iv_][ip] = vp;
                a[iw_][ip] = wp;
                a[ie_][ip] = v2;
              }

              a[iq_][nVar] += qp;
              a[iu_][nVar] += qp * up;
              a[iv_][nVar] += qp * vp;
              a[iw_][nVar] += qp * wp;
              a[ie_][nVar] += qp * v2;
            }

            auto linear_solver_Gauss_Elimination = [&a, &x, &nVar]() {
              for (int i = 0; i < nVar - 1; i++) {
                if (a[i][i] == 0)
                  return false;
                for (int k = i + 1; k < nVar; k++) {

                  Real t = a[k][i] / a[i][i];
                  for (int j = 0; j <= nVar; j++)
                    a[k][j] = a[k][j] - t * a[i][j];
                }
              }

              for (int i = nVar - 1; i >= 0; i--) {
                x[i] = a[i][nVar];
                for (int j = i + 1; j < nVar; j++) {
                  if (j != i)
                    x[i] = x[i] - a[i][j] * x[j];
                }
                if (a[i][i] == 0) {
                  return false;
                } else {
                  x[i] = x[i] / a[i][i];
                }
              }

              return true;
            };

            bool isSolved = linear_solver_Gauss_Elimination();

            if (isSolved) {
              // All the particle weights should have the same sign.
              Real qt = a[iq_][nVar];
              for (int ip = 0; ip < nPartCombine - 1; ip++) {
                if (qt * x[ip] < 0) {
                  isSolved = false;
                  break;
                }
              }
            }
            if (!isSolved)
              continue;
            //----------------------------------------------

            //----------------------------------------------

            // 1) Adjust the location to conserve the mass center.
            // 2) Adjust the weight to conserve mass and energy.
            const bool doConserveMassCenter = false;
            if (doConserveMassCenter) {
              Real centerxOld = 0, centeryOld = 0, centerzOld = 0;
              for (int ip = 0; ip < nPartCombine; ip++) {
                centerxOld += particles[idx_I[ip]].rdata(iqp_) *
                              particles[idx_I[ip]].pos(ix_);
                centeryOld += particles[idx_I[ip]].rdata(iqp_) *
                              particles[idx_I[ip]].pos(iy_);
                centerzOld += particles[idx_I[ip]].rdata(iqp_) *
                              particles[idx_I[ip]].pos(iz_);
              }

              Real centerxNew = 0, centeryNew = 0, centerzNew = 0;
              Real qtotal = 0;
              for (int ip = 0; ip < nPartCombine - 1; ip++) {
                centerxNew += x[ip] * particles[idx_I[ip]].pos(ix_);
                centeryNew += x[ip] * particles[idx_I[ip]].pos(iy_);
                centerzNew += x[ip] * particles[idx_I[ip]].pos(iz_);
                qtotal += x[ip];
              }

              const Real invQtotal = 1. / qtotal;
              const Real dPartX = (centerxOld - centerxNew) * invQtotal;
              const Real dPartY = (centeryOld - centeryNew) * invQtotal;
              const Real dPartZ = (centerzOld - centerzNew) * invQtotal;

              for (int ip = 0; ip < nPartCombine - 1; ip++) {
                Real xpNew = particles[idx_I[ip]].pos(ix_) + dPartX;
                Real ypNew = particles[idx_I[ip]].pos(iy_) + dPartY;
                Real zpNew = particles[idx_I[ip]].pos(iz_) + dPartZ;

                xpNew = bound(xpNew, xMin, xMax);
                ypNew = bound(ypNew, xMin, xMax);
                zpNew = bound(zpNew, xMin, xMax);

                particles[idx_I[ip]].pos(ix_) = xpNew;
                particles[idx_I[ip]].pos(iy_) = ypNew;
                particles[idx_I[ip]].pos(iz_) = zpNew;
              }
            } // if doConserveMassCenter

            // Adjust weight.
            for (int ip = 0; ip < nPartCombine - 1; ip++) {
              auto& p = particles[idx_I[ip]];
              p.rdata(iqp_) = x[ip];
            }
            // Mark for deletion
            particles[idx_I[nPartCombine - 1]].id() = -1;
          }
        }
  }
}

//==========================================================
bool Particles::do_inject_particles_for_this_cell(
    const amrex::Box& bx, const amrex::Array4<const int>& status, const int i,
    const int j, const int k) {

  // This cell should be a boundary cell at least.
  if (status(i, j, k) != iBoundary_)
    return false;

  for (int iloop = 1; iloop <= 3; iloop++) {
    // iloop==1: loop through faces;
    // iloop==2: loop through edges;
    // iloop==3: loop through corners;

    for (int di = -1; di <= 1; di++)
      for (int dj = -1; dj <= 1; dj++)
        for (int dk = -1; dk <= 1; dk++) {
          const int sum = abs(di) + abs(dj) + abs(dk);
          if (iloop != sum)
            continue;

          if (status(i + di, j + dj, k + dk) != iBoundary_) {
            // The first neighbor cell that is NOT a boundary cell.
            if (bx.contains({ i + di, j + dj, k + dk })) {
              return true;
            } else {
              return false;
            }
          }
        }
  }
  Abort("do_inject_particles_for_this_cell:something is wrong!");
  return false; // to suppress compilation warning.
}
