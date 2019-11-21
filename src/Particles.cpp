#include "Particles.h"
#include "Utility.h"

using namespace amrex;

Particles::Particles(const Geometry& geom, const DistributionMapping& dm,
                     const BoxArray& ba, const int speciesIDIn,
                     const Real chargeIn, const Real massIn,
                     const IntVect& nPartPerCellIn)
    : ParticleContainer<4, 0, 0, 0>(geom, dm, ba),
      speciesID(speciesIDIn),
      charge(chargeIn),
      mass(massIn),
      nPartPerCell(nPartPerCellIn) {}

void Particles::add_particles_cell(const MFIter& mfi,
                                   const FluidInterface& fluidInterface, int i,
                                   int j, int k) {
  int ig, jg, kg, nxcg, nycg, nzcg, iCycle, npcel, nRandom = 7;
  // Why +1? for comparison with iPIC3D.-----
  ig = i + 1;
  jg = j + 1;
  kg = k;
  if (fluidInterface.getnDim() > 2)
    kg++; // just for comparison with iPIC3D;
  //----------------------------------------

  nxcg = fluidInterface.getFluidNxc();
  nycg = fluidInterface.getFluidNyc();
  nzcg = fluidInterface.getFluidNzc();
  iCycle = fluidInterface.getCycle();
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

  const int lev = 0;
  auto& particles =
      GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

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

        double q = (charge / mass / fabs(charge / mass)) *
                   (fluidInterface.get_number_density(mfi, x, y, z, speciesID) /
                    npcel) *
                   vol;

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
          u += fluidInterface.get_ux(mfi, x, y, z, speciesID);
          v += fluidInterface.get_uy(mfi, x, y, z, speciesID);
          w += fluidInterface.get_uz(mfi, x, y, z, speciesID);

          ParticleType p;
          p.id() = ParticleType::NextID();
          p.cpu() = ParallelDescriptor::MyProc();
          p.pos(ix_) = x; // + plo[ix_];
          p.pos(iy_) = y; // + plo[iy_];
          p.pos(iz_) = z; // + plo[iz_];
          p.rdata(iup_) = u;
          p.rdata(ivp_) = v;
          p.rdata(iwp_) = w;
          p.rdata(iqp_) = q;
          particles.push_back(p);
          Print() << "i = " << i << " j = " << j << " k = " << k << p
                  << std::endl;
        }
      }
}

void Particles::add_particles_domain(const FluidInterface& fluidInterface) {
  BL_PROFILE("Particles::add_particles");

  const int lev = 0;
  for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const Box& tile_box = mfi.validbox();
    const auto lo = amrex::lbound(tile_box);
    const auto hi = amrex::ubound(tile_box);

    for (int i = lo.x; i <= hi.x; ++i)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int k = lo.z; k <= hi.z; ++k) {
          add_particles_cell(mfi, fluidInterface, i, j, k);
        }
  }
}

void Particles::inject_particles_at_boundary(
    const FluidInterface& fluidInterface) {

  if (nVirGst != 1)
    Abort("The function assumes nVirGst==1!!");

  const int lev = 0;

  const Box& gbx = Geom(0).Domain();

  for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    // Since nVirGst==1, fill in the cells in the 'validbox'.
    const Box& tile_box = mfi.validbox();
    const auto lo = amrex::lbound(tile_box);
    const auto hi = amrex::ubound(tile_box);

    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    int iMin = lo.x, jMin = lo.y, kMin = lo.z;

    bool doFillLeftX = false, doFillLeftY = false, doFillLeftZ = false;
    bool doFillRightX = false, doFillRightY = false, doFillRightZ = false;

    if (!(Geom(0).isPeriodic(ix_)) && gbx.bigEnd(ix_) == hi.x) {
      doFillRightX = true;
    }
    if ((!Geom(0).isPeriodic(iy_)) && gbx.bigEnd(iy_) == hi.y) {
      doFillRightY = true;
    }
    if ((!Geom(0).isPeriodic(iz_)) && gbx.bigEnd(iz_) == hi.z) {
      doFillRightZ = true;
    }

    if (!Geom(0).isPeriodic(ix_) && gbx.smallEnd(ix_) == lo.x) {
      doFillLeftX = true;
    }
    if (!Geom(0).isPeriodic(iy_) && gbx.smallEnd(iy_) == lo.y) {
      doFillLeftY = true;
    }
    if (!Geom(0).isPeriodic(iz_) && gbx.smallEnd(iz_) == lo.z) {
      doFillLeftZ = true;
    }

    if (doFillLeftX) {
      for (int i = iMin; i <= iMin - 1 + nVirGst; ++i)
        for (int j = jMin; j <= jMax; ++j)
          for (int k = kMin; k <= kMax; ++k) {

            add_particles_cell(mfi, fluidInterface, i, j, k);
          }
      iMin += nVirGst;
    }

    if (doFillRightX) {
      for (int i = iMax + 1 - nVirGst; i <= iMax; ++i)
        for (int j = jMin; j <= jMax; ++j)
          for (int k = kMin; k <= kMax; ++k) {

            add_particles_cell(mfi, fluidInterface, i, j, k);
          }
      iMax -= nVirGst;
    }

    if (doFillLeftY) {
      for (int i = iMin; i <= iMax; ++i)
        for (int j = jMin; j <= jMin - 1 + nVirGst; ++j)
          for (int k = kMin; k <= kMax; ++k) {

            add_particles_cell(mfi, fluidInterface, i, j, k);
          }
      jMin += nVirGst;
    }

    if (doFillRightY) {
      for (int i = iMin; i <= iMax; ++i)
        for (int j = jMax + 1 - nVirGst; j <= jMax; ++j)
          for (int k = kMin; k <= kMax; ++k) {

            add_particles_cell(mfi, fluidInterface, i, j, k);
          }
      jMax -= nVirGst;
    }

    if (doFillLeftZ) {
      for (int i = iMin; i <= iMax; ++i)
        for (int j = jMin; j <= jMax; ++j)
          for (int k = kMin; k <= kMin - 1 + nVirGst; ++k) {
            add_particles_cell(mfi, fluidInterface, i, j, k);
          }
      kMin += nVirGst;
    }

    if (doFillRightZ) {
      for (int i = iMin; i <= iMax; ++i)
        for (int j = jMin; j <= jMax; ++j)
          for (int k = kMax + 1 - nVirGst; k <= kMax; ++k) {

            add_particles_cell(mfi, fluidInterface, i, j, k);
          }
      kMax -= nVirGst;
    }
  }

  Print() << " i = " << speciesID << " number of particles after injection = "
          << TotalNumberOfParticles(false) << std::endl;
}

void Particles::sum_to_center(amrex::MultiFab& netChargeMF,
                              amrex::UMultiFab<RealCMM>& centerMM,
                              bool doNetChargeOnly) {
  BL_PROFILE("Particles::sum_to_center");

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
        loIdx[i] = floor(dShift[i]);
        dShift[i] = dShift[i] - loIdx[i];
      }
      Real coef[2][2][2];
      part_grid_interpolation_coef(dShift, coef);
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

        Real wg_D[3];

        for (int k1 = kMin; k1 <= kMax; k1++)
          for (int j1 = jMin; j1 <= jMax; j1++)
            for (int i1 = iMin; i1 <= iMax; i1++) {

              for (int iDim = 0; iDim < 3; iDim++) {
                wg_D[iDim] =
                    invVol *
                    weights_IIID[i1 - iMin][j1 - jMin][k1 - kMin][iDim];
              }

              for (int j2 = jMin; j2 <= jMax; j2++)
                for (int i2 = iMin; i2 <= iMax; i2++)
                  for (int k2 = kMin; k2 <= kMax; k2++) {
                    Real weight = 0;
                    for (int iDim = 0; iDim < 3; iDim++) {
                      weight +=
                          wg_D[iDim] *
                          weights_IIID[i2 - iMin][j2 - jMin][k2 - kMin][iDim];
                    }

                    const int ip = i2 - i1 + 1;
                    const int jp = j2 - j1 + 1;
                    const int kp = k2 - k1 + 1;
                    const int gp = ip * 9 + jp * 3 + kp;
                    mmArr(i1, j1, k1).data[gp] += weight * fabs(qp);
                  }
            }
      } // if doChargeOnly

    } // for p
  }
}

void Particles::sum_moments(MultiFab& momentsMF, UMultiFab<RealMM>& nodeMM,
                            MultiFab& nodeBMF, Real dt) {
  BL_PROFILE("Particles::sum_moments");
  const auto& plo = Geom(0).ProbLo();

  momentsMF.setVal(0.0);

  Real qdto2mc = charge / mass * 0.5 * dt;

  const auto& invDx = Geom(0).InvCellSize();
  const Real invVol = invDx[ix_] * invDx[iy_] * invDx[iz_];

  const int lev = 0;
  for (ParticlesIter pti(*this, lev); pti.isValid(); ++pti) {
    Array4<Real> const& momentsArr = momentsMF[pti].array();
    Array4<Real const> const& nodeBArr = nodeBMF[pti].array();
    Array4<RealMM> const& mmArr = nodeMM[pti].array();

    // Particle container is defined based on a center box.
    // So Convert it into a node box here.
    const Box& box = convert(pti.validbox(), { 1, 1, 1 });
    // Print() << " box = " << box << std::endl;

    const auto& particles = pti.GetArrayOfStructs();
    for (const auto& p : particles) {
      Print() << "particle = " << p << std::endl;

      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      int loIdx[3];
      Real dShift[3];
      for (int i = 0; i < 3; i++) {
        dShift[i] = (p.pos(i) - plo[i]) * invDx[i];
        loIdx[i] = floor(dShift[i]);
        dShift[i] = dShift[i] - loIdx[i];
      }

      Real coef[2][2][2];
      part_grid_interpolation_coef(dShift, coef);
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

      int count = 0;
      int ip, jp, kp;         // Indexes for g'
      double wg, wgp, weight; // W_pg, W_pg'

      const int iMin = loIdx[ix_];
      const int jMin = loIdx[iy_];
      const int kMin = loIdx[iz_];
      const int iMax = iMin + 2;
      const int jMax = jMin + 2;
      const int kMax = kMin + 2;

      int gp = -1;
      for (int k1 = kMin; k1 < kMax; k1++)
        for (int j1 = jMin; j1 < jMax; j1++)
          for (int i1 = iMin; i1 < iMax; i1++) {
            wg = coef[i1 - iMin][j1 - jMin][k1 - kMin];
            for (int k2 = kMin; k2 < kMax; k2++) {
              kp = k2 - k1 + 1;
              if (kp > 0) {
                for (int j2 = jMin; j2 < jMax; j2++) {
                  jp = j2 - j1 + 1;
                  for (int i2 = iMin; i2 < iMax; i2++) {
                    ip = i2 - i1 + 1;
                    weight = wg * coef[i2 - iMin][j2 - jMin][k2 - kMin];
                    gp = kp * 9 + jp * 3 + ip;
                    const int idx0 = gp * 9;
                    for (int idx = 0; idx < 9; idx++) {
                      mmArr(i1, j1, k1).data[idx0 + idx] += alpha[idx] * weight;
                    }
                  } // k2
                }   // j2
              }     // if (ip > 0)
            }       // i2
          }         // k1

      //----- Mass matrix calculation end--------------

      //-------Moments begin---------
      Real pMoments[nMoments];

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

      for (int kk = 0; kk < 2; kk++)
        for (int jj = 0; jj < 2; jj++)
          for (int ii = 0; ii < 2; ii++)
            for (int iVar = 0; iVar < nMoments; iVar++) {
              momentsArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk,
                         iVar) += coef[ii][jj][kk] * pMoments[iVar];
            }

      //-------Moments end---------

    } // for p

    for (MFIter mfi(nodeMM); mfi.isValid(); ++mfi) {
      // Finalize the mass matrix calculation.
      const Box& box = mfi.validbox();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      Array4<RealMM> const& mmArr = nodeMM[mfi].array();

      int gps, gpr; // gp_send, gp_receive
      for (int k1 = lo.z; k1 <= hi.z; k1++)
        for (int j1 = lo.y; j1 <= hi.y; j1++)
          for (int i1 = lo.x; i1 <= hi.x; i1++) {
            const int kp = 2, kpr = 0;
            const int kr = k1 + kp - 1;
            for (int jp = 0; jp < 3; jp++) {
              const int jr = j1 + jp - 1;
              const int jpr = 2 - jp;
              for (int ip = 0; ip < 3; ip++) {
                const int ir = i1 + ip - 1;
                const int ipr = 2 - ip;
                gpr = jpr * 3 + ipr;
                gps = 18 + jp * 3 + ip; // gps = kp*9+jp*3+kp
                for (int idx = 0; idx < 9; idx++) {
                  mmArr(ir, jr, kr).data[gpr * 9 + idx] =
                      mmArr(i1, j1, k1).data[gps * 9 + idx];
                } // idx
              }   // kp
            }     // jp
          }       // k1
    }
  }

  momentsMF.mult(invVol);

  momentsMF.SumBoundary(Geom(0).periodicity());
  // Print()<<momentsMF<<std::endl;

  // FillBoundary seems unnecessary. --Yuxi
  momentsMF.FillBoundary(Geom(0).periodicity());

  convert_to_fluid_moments(momentsMF);
}

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

void Particles::mover(const amrex::MultiFab& nodeEMF,
                      const amrex::MultiFab& nodeBMF, amrex::Real dt) {
  BL_PROFILE("Particles::mover");

  const auto& plo = Geom(0).ProbLo();

  const Real dtLoc = dt;
  const Real qdto2mc = charge / mass * 0.5 * dt;

  const auto& invDx = Geom(0).InvCellSize();
  // const Real invVol = invDx[ix_] * invDx[iy_] * invDx[iz_];
  // Print() << " invVol = " << invVol << std::endl;

  const int lev = 0;
  for (ParticlesIter pti(*this, lev); pti.isValid(); ++pti) {
    const Array4<Real const>& nodeEArr = nodeEMF[pti].array();
    const Array4<Real const>& nodeBArr = nodeBMF[pti].array();

    // Particle container is defined based on a center box.
    // So Convert it into a node box here.
    const Box& box = convert(pti.validbox(), { 1, 1, 1 });
    // Print() << " box = " << box << std::endl;

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
        loIdx[i] = floor(dShift[i]);
        dShift[i] = dShift[i] - loIdx[i];
      }

      Real coef[2][2][2];
      part_grid_interpolation_coef(dShift, coef);
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
      if (is_outside(p)) {
        p.id() = -1;
      }

    } // for p
  }   // for pti

  Print() << " i = " << speciesID << " number of particles before deletion = "
          << TotalNumberOfParticles(false) << std::endl;

  // This function distributes particles to proper processors and apply
  // periodic boundary conditions if needed.
  Redistribute();

  Print() << " i = " << speciesID << " number of particles after deletion = "
          << TotalNumberOfParticles(false) << std::endl;
}

void Particles::divE_correct_position(const amrex::MultiFab& phiMF) {

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
        loIdx[i] = floor(dShift[i]);
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

        Real eps_D[3] = { 0, 0, 0 };

        for (int k = 0; k < 2; k++)
          for (int j = 0; j < 2; j++)
            for (int i = 0; i < 2; i++)
              for (int iDim = 0; iDim < 3; iDim++) {
                eps_D[iDim] += fourPI * phiArr(iMin + i, jMin + j, kMin + k) *
                               weights_IIID[i][j][k][iDim];
              }

        for (int iDim = 0; iDim < 3; iDim++)
          eps_D[iDim] *= coef;

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
        }
      }

    } // for p
  }

  ParallelDescriptor::ReduceRealMax(epsMax,
                                    ParallelDescriptor::IOProcessorNumber());

  // A reduction is required here...
  Print() << "Particle position correction: epsMax = " << epsMax
          << " for species " << speciesID << std::endl;

  Redistribute();
}
