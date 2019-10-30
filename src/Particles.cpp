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
                                   const FluidPicInterface& fluidInterface,
                                   int iBlock, int i, int j, int k, int loi,
                                   int loj, int lok) {
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

  double x, y, z; // Particle location.

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

        x = (ii + randNum()) * (dx[ix_] / nPartPerCell[ix_]) + i * dx[ix_];
        y = (jj + randNum()) * (dx[iy_] / nPartPerCell[iy_]) + j * dx[iy_];
        z = (kk + randNum()) * (dx[iz_] / nPartPerCell[iz_]) + k * dx[iz_];

        double q =
            (charge / mass / fabs(charge / mass)) *
            (fluidInterface.getPICRhoNum(iBlock, x, y, z, speciesID) / npcel) *
            vol;

        if (q != 0) {
          double rand;
          double u, v, w;
          double rand1 = randNum();
          double rand2 = randNum();
          double rand3 = randNum();
          double rand4 = randNum();

          if (fluidInterface.getUseAnisoP() &&
              (speciesID > 0 || fluidInterface.get_useElectronFluid())) {
            fluidInterface.setPICAnisoUth(iBlock, x, y, z, &u, &v, &w, rand1,
                                          rand2, rand3, rand4, speciesID);
          } else {
            fluidInterface.setPICIsoUth(iBlock, x, y, z, &u, &v, &w, rand1,
                                        rand2, rand3, rand4, speciesID);
          }
          u += fluidInterface.getPICUx(iBlock, x, y, z, speciesID);
          v += fluidInterface.getPICUy(iBlock, x, y, z, speciesID);
          w += fluidInterface.getPICUz(iBlock, x, y, z, speciesID);

          ParticleType p;
          p.id() = ParticleType::NextID();
          p.cpu() = ParallelDescriptor::MyProc();
          p.pos(ix_) = x + plo[ix_];
          p.pos(iy_) = y + plo[iy_];
          p.pos(iz_) = z + plo[iz_];
          p.rdata(iup_) = u;
          p.rdata(ivp_) = v;
          p.rdata(iwp_) = w;
          p.rdata(iqp_) = q;
          particles.push_back(p);
        }
      }
}

void Particles::add_particles_domain(const FluidPicInterface& fluidInterface) {
  BL_PROFILE("Particles::add_particles");

  const int lev = 0;
  const auto dx = Geom(lev).CellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();

  int iBlock = 0;
  for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const Box& tile_box = mfi.validbox();
    const auto lo = amrex::lbound(tile_box);
    const auto hi = amrex::ubound(tile_box);

    for (int i = lo.x; i <= hi.x; ++i)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int k = lo.z; k <= hi.z; ++k) {
          add_particles_cell(mfi, fluidInterface, iBlock, i, j, k, lo.x, lo.y,
                             lo.z);
        }

    iBlock++;
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
      // Print() << "particle = " << p << std::endl;

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
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++) {
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
      for (int i1 = iMin; i1 < iMax; i1++)
        for (int j1 = jMin; j1 < jMax; j1++)
          for (int k1 = kMin; k1 < kMax; k1++) {
            wg = coef[i1 - iMin][j1 - jMin][k1 - kMin];
            for (int i2 = iMin; i2 < iMax; i2++) {
              ip = i2 - i1 + 1;
              if (ip > 0) {
                for (int j2 = jMin; j2 < jMax; j2++) {
                  jp = j2 - j1 + 1;
                  for (int k2 = kMin; k2 < kMax; k2++) {
                    kp = k2 - k1 + 1;
                    weight = wg * coef[i2 - iMin][j2 - jMin][k2 - kMin];
                    gp = ip * 9 + jp * 3 + kp;
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

      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
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

      int gps, gpr;                         // gp_send, gp_receive
      for (int i1 = lo.x; i1 <= hi.x; i1++) // Change the order here!!!!--Yuxi
        for (int j1 = lo.y; j1 <= hi.y; j1++)
          for (int k1 = lo.z; k1 <= hi.z; k1++) {
            const int ip = 2, ipr = 0;
            const int ir = i1 + ip - 1;
            for (int jp = 0; jp < 3; jp++) {
              const int jr = j1 + jp - 1;
              const int jpr = 2 - jp;
              for (int kp = 0; kp < 3; kp++) {
                const int kr = k1 + kp - 1;
                const int kpr = 2 - kp;
                gpr = jpr * 3 + kpr;    // gpr = ipr*9 + jpr*3 + kpr
                gps = 18 + jp * 3 + kp; // gps = ip*9+jp*3+kp

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

    } // for p
  }   // for pti

  this->Redistribute();
}
