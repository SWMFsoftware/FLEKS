#include <algorithm>
#include <cstdlib>

#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;
using namespace std;

//==========================================================
template <int NStructReal, int NStructInt>
Particles<NStructReal, NStructInt>::Particles(
    amrex::AmrCore* amrcore, FluidInterface* const fluidIn, TimeCtr* const tcIn,
    const int speciesIDIn, const Real chargeIn, const Real massIn,
    const IntVect& nPartPerCellIn, TestCase tcase)
    : AmrParticleContainer<NStructReal, NStructInt>(amrcore),
      fi(fluidIn),
      tc(tcIn),
      speciesID(speciesIDIn),
      charge(chargeIn),
      mass(massIn),
      nPartPerCell(nPartPerCellIn),
      testCase(tcase) {
  do_tiling = true;

  qom = charge / mass;
  qomSign = qom >= 0 ? 1 : -1;

  nLev = maxLevel() + 1;

  cellStatus.resize(nLev);

  plo.resize(nLev);
  phi.resize(nLev);
  dx.resize(nLev);
  invDx.resize(nLev);
  invVol.resize(nLev);

  for (int iLev = 0; iLev < nLev; iLev++) {
    invVol[iLev] = 1;
    for (int i = 0; i < nDim; i++) {
      tile_size[i] = 1;
      plo[iLev][i] = Geom(iLev).ProbLo(i);
      phi[iLev][i] = Geom(iLev).ProbHi(i);
      dx[iLev][i] = Geom(iLev).CellSize(i);
      invDx[iLev][i] = Geom(iLev).InvCellSize(i);
      invVol[iLev] *= invDx[iLev][i];
    }
  }

  // The following line is used to avoid an MPI bug (feature?) on Frontera. It
  // should be removed after the bug being fixed.
  SetUseUnlink(false);
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::outflow_bc(const amrex::MFIter& mfi,
                                                    const int ig, const int jg,
                                                    const int kg, const int ip,
                                                    const int jp,
                                                    const int kp) {
  const int iLev = 0;

  IntVect idxGst(ig, jg, kg);
  ParticleTileType& pGst = get_particle_tile(iLev, mfi, idxGst);

  IntVect idxPhy(ip, jp, kp);
  ParticleTileType& pPhy = get_particle_tile(iLev, mfi, idxGst);
  auto& phyParts = pPhy.GetArrayOfStructs();

  Real dxshift[3] = { 0, 0, 0 };
  for (int i = 0; i < nDim; i++) {
    dxshift[i] = Geom(iLev).CellSize(i) * (idxGst[i] - idxPhy[i]);
  }

  Vector<ParticleType> pList;
  for (const auto& p : phyParts) {
    IntVect iv = Index(p, iLev);
    // Q: Why do we need to check if the physical domain contains the particle?
    // A: Even if tiling with tile_size=1 is used, it seems the ghost cells
    // still share the the same tile with a physical cell. Therefore, we need to
    // make sure a particle in a "physical tile" is actually inside the physical
    // domain.
    if (mfi.validbox().contains(IntVect(iv))) {
      // TODO: Check NextID
      ParticleType pNew = p;

      pNew.id() = ParticleType::NextID();
      pNew.cpu() = ParallelDescriptor::MyProc();

      for (int i = 0; i < nDim; i++) {
        pNew.pos(i) = p.pos(i) + dxshift[i];
      }

      pList.push_back(pNew);
    }
  }

  // Q: Why do not push the new particles into pGst inside previous loop?
  // A: Sometimes, if not always, pPhy and pGst share the same tile. Previous
  // for-loop loops through al particles in pPhy. If we push the new particles
  // into pGst, which is the same as pPhy sometimes, the loop behavior is not
  // well defined.
  for (auto& p : pList) {
    pGst.push_back(p);
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_particles_cell(
    const int iLev, const MFIter& mfi, const int i, const int j, const int k,
    const FluidInterface& interface, IntVect ppc, const Vel tpVel, Real dt) {

  // If true, initialize the test particles with user defined velocities instead
  // of from fluid.
  bool userState = (tpVel.tag == speciesID);

  // If dt >0, it suggests the 'density' obtained from interface is actually the
  // density changing rate.
  if (dt <= 0)
    dt = 1;

  IntVect nPPC = nPartPerCell;
  if (!(ppc == 0)) {
    nPPC = ppc;
  }

  set_random_seed(iLev, i, j, k, nPPC);

  Real x, y, z; // Particle location

  const Real vol = dx[iLev][ix_] * dx[iLev][iy_] * dx[iLev][iz_];
  const int npcel = nPPC[ix_] * nPPC[iy_] * nPPC[iz_];
  const Real vol2Npcel = qomSign * vol / npcel;

  ParticleTileType& particles = get_particle_tile(iLev, mfi, i, j, k);

  int icount = 0;
  // Loop over particles inside grid cell i, j, k
  for (int ii = 0; ii < nPPC[ix_]; ii++)
    for (int jj = 0; jj < nPPC[iy_]; jj++)
      for (int kk = 0; kk < nPPC[iz_]; kk++) {

        x = (ii + randNum()) * (dx[iLev][ix_] / nPPC[ix_]) + i * dx[iLev][ix_] +
            plo[iLev][ix_];
        y = (jj + randNum()) * (dx[iLev][iy_] / nPPC[iy_]) + j * dx[iLev][iy_] +
            plo[iLev][iy_];
        z = (kk + randNum()) * (dx[iLev][iz_] / nPPC[iz_]) + k * dx[iLev][iz_] +
            plo[iLev][iz_];

        // If the particle weight is sampled in a random location, the sum of
        // particle mass is NOT the same as the integral of the grid density.
        // It is more convenient for debugging if mass is exactly conserved. For
        // a production run, it makes little difference.
        Real x0 = (ii + 0.5) * (dx[iLev][ix_] / nPPC[ix_]) + i * dx[iLev][ix_] +
                  plo[iLev][ix_];
        Real y0 = (jj + 0.5) * (dx[iLev][iy_] / nPPC[iy_]) + j * dx[iLev][iy_] +
                  plo[iLev][iy_];
        Real z0 = (kk + 0.5) * (dx[iLev][iz_] / nPPC[iz_]) + k * dx[iLev][iz_] +
                  plo[iLev][iz_];

        double q = vol2Npcel * interface.get_number_density(mfi, x0, y0, z0,
                                                            speciesID, iLev);
        if (q != 0) {
          Real u, v, w;
          double rand1 = randNum();
          double rand2 = randNum();
          double rand3 = randNum();
          double rand4 = randNum();

          double uth = (userState ? tpVel.vth : -1);
          if (!is_neutral() && interface.get_UseAnisoP() &&
              (speciesID > 0 || interface.get_useElectronFluid())) {
            interface.set_particle_uth_aniso(iLev, mfi, x, y, z, &u, &v, &w,
                                             rand1, rand2, rand3, rand4,
                                             speciesID, uth, uth);
          } else {
            interface.set_particle_uth_iso(iLev, mfi, x, y, z, &u, &v, &w,
                                           rand1, rand2, rand3, rand4,
                                           speciesID, uth);
          }

          Real uBulk = userState
                           ? tpVel.vx
                           : interface.get_ux(mfi, x, y, z, speciesID, iLev);
          Real vBulk = userState
                           ? tpVel.vy
                           : interface.get_uy(mfi, x, y, z, speciesID, iLev);
          Real wBulk = userState
                           ? tpVel.vz
                           : interface.get_uz(mfi, x, y, z, speciesID, iLev);

          if (testCase == TwoStream && qom < 0 && icount % 2 == 0) {
            // Electron only (qom<0)
            uBulk = -uBulk;
            vBulk = -vBulk;
            wBulk = -wBulk;
          }
          // printf("p u=%e, v=%e, w=%e, ubulk=%e, vbulk=%e, wbulk=%e \n", u, v,
          // w,
          //        uBulk, vBulk, wBulk);
          u += uBulk;
          v += vBulk;
          w += wBulk;

          ParticleType p;
          if (ParticleType::the_next_id >= amrex::LastParticleID) {
            // id should not be larger than LastParticleID. This is a bad
            // solution, since the ID becomes nonunique. --Yuxi
            p.id() = amrex::LastParticleID;
          } else {
            p.id() = ParticleType::NextID();
          }
          p.cpu() = ParallelDescriptor::MyProc();
          p.pos(ix_) = x;
          p.pos(iy_) = y;
          p.pos(iz_) = z;
          p.rdata(iup_) = u;
          p.rdata(ivp_) = v;
          p.rdata(iwp_) = w;
          // Convert 'density changing rate' to 'density' if necessary.
          p.rdata(iqp_) = q * dt;

          if (NStructInt > 0) {
            // For test particle only.
            p.idata(iRecordCount_) = 0;
          }

          particles.push_back(p);
          // AllPrint() << "p=" << p << std::endl;
          icount++;
        }
      }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_particles_source(
    const FluidInterface& interface, const FluidInterface* const stateOH,
    Real dt, IntVect ppc, const bool doSelectRegion) {
  timing_func("Particles::add_particles_source");

  // 1. Inject particles for physical cells.
  const int iLev = 0;
  for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
    const Box& tile_box = mfi.validbox();
    const auto lo = amrex::lbound(tile_box);
    const auto hi = amrex::ubound(tile_box);

    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    int iMin = lo.x, jMin = lo.y, kMin = lo.z;

    for (int i = iMin; i <= iMax; ++i)
      for (int j = jMin; j <= jMax; ++j)
        for (int k = kMin; k <= kMax; ++k) {
          bool doAdd = true;
#ifdef _PT_COMPONENT_
          if (stateOH && doSelectRegion) {
            const int iFluid = 0;
            const int iRegion =
                stateOH->get_neu_source_region(mfi, i, j, k, iFluid, iLev);
            doAdd = (iRegion == speciesID);
          }
#endif
          if (doAdd)
            add_particles_cell(iLev, mfi, i, j, k, interface, ppc, Vel(), dt);
        }
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_particles_domain(
    const amrex::Vector<amrex::iMultiFab>& iRefinement) {
  timing_func("Particles::add_particles_domain");

  for (int iLev = 0; iLev <= finestLevel(); iLev++) {
    for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
      const auto& status = cellStatus[iLev][mfi].array();
      const auto& iRefine = iRefinement[iLev][mfi].array();
      const Box& bx = mfi.validbox();
      const auto lo = amrex::lbound(bx);
      const auto hi = amrex::ubound(bx);

      int iMax = hi.x, jMax = hi.y, kMax = hi.z;
      int iMin = lo.x, jMin = lo.y, kMin = lo.z;

      //      Print() << "iLev = " << iLev << " bx = " << bx << std::endl;

      for (int i = iMin; i <= iMax; ++i)
        for (int j = jMin; j <= jMax; ++j)
          for (int k = kMin; k <= kMax; ++k) {
            if (test_bit(status(i, j, k), iDigitNew_) &&
                iRefine(i, j, k) == iNotRefined) {
              // printf("add particles iLev = %d, i = %d, j = %d, k = %d\n",
              // iLev,
              //        i, j, k);
              add_particles_cell(iLev, mfi, i, j, k, *fi);
            }
          }
    }
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::inject_particles_at_boundary(
    const FluidInterface* fiIn, Real dt, IntVect ppc) {
  timing_func("Particles::inject_particles_at_boundary");

  // Only inject nGstInject layers.
  const int nGstInject = 1;

  // By default, use fi for injecting particles.
  const FluidInterface* fiTmp = fi;
  if (fiIn)
    fiTmp = fiIn;

  // Assume the boundary cells aren't and cannot be refined. --Yuxi
  int iLev = 0;

  for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
    const auto& status = cellStatus[iLev][mfi].array();
    const Box& bx = mfi.validbox();
    const IntVect lo = IntVect(bx.loVect());
    const IntVect hi = IntVect(bx.hiVect());
    // IntVect mid = (lo + hi) / 2;

    IntVect idxMin = lo, idxMax = hi;

    for (int iDim = 0; iDim < fiTmp->get_fluid_dimension(); iDim++) {
      idxMin[iDim] -= nGstInject;
      idxMax[iDim] += nGstInject;
    }

    for (int i = idxMin[ix_]; i <= idxMax[ix_]; ++i)
      for (int j = idxMin[iy_]; j <= idxMax[iy_]; ++j)
        for (int k = idxMin[iz_]; k <= idxMax[iz_]; ++k) {
          int isrc, jsrc, ksrc;
          if (do_inject_particles_for_this_cell(bx, status, i, j, k, isrc, jsrc,
                                                ksrc)) {
            if (((bc.lo[ix_] == bc.outflow) && i < lo[ix_]) ||
                ((bc.hi[ix_] == bc.outflow) && i > hi[ix_]) ||
                ((bc.lo[iy_] == bc.outflow) && j < lo[iy_]) ||
                ((bc.hi[iy_] == bc.outflow) && j > hi[iy_]) ||
                ((bc.lo[iz_] == bc.outflow) && k < lo[iz_]) ||
                ((bc.hi[iz_] == bc.outflow) && k > hi[iz_])) {
              outflow_bc(mfi, i, j, k, isrc, jsrc, ksrc);
            } else {
              add_particles_cell(iLev, mfi, i, j, k, *fiTmp, ppc, Vel(), dt);
            }
          }
        }
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::sum_to_center(
    amrex::MultiFab& netChargeMF, amrex::UMultiFab<RealCMM>& centerMM,
    bool doNetChargeOnly) {
  timing_func("Particles::sum_to_center");

  const int iLev = 0;
  for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    Array4<Real> const& chargeArr = netChargeMF[pti].array();
    Array4<RealCMM> const& mmArr = centerMM[pti].array();

    const auto& particles = pti.GetArrayOfStructs();

    for (const auto& p : particles) {
      // Print() << "particle = " << p << std::endl;

      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      IntVect loIdx;
      RealVect dShift;
      for (int i = 0; i < nDim; i++) {
        // plo is the corner location => -0.5
        dShift[i] = (p.pos(i) - plo[iLev][i]) * invDx[iLev][i] - 0.5;
        loIdx[i] = fastfloor(dShift[i]); // floor() is slow.
        dShift[i] = dShift[i] - loIdx[i];
      }
      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

      const Real cTmp = qp * invVol[iLev];
      for (int kk = 0; kk < 2; kk++)
        for (int jj = 0; jj < 2; jj++)
          for (int ii = 0; ii < 2; ii++) {
            chargeArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk) +=
                coef[ii][jj][kk] * cTmp;
          }

      if (!doNetChargeOnly) {
        Real weights_IIID[2][2][2][nDim];
        //----- Mass matrix calculation begin--------------
        const Real xi0 = dShift[ix_] * dx[iLev][ix_];
        const Real eta0 = dShift[iy_] * dx[iLev][iy_];
        const Real zeta0 = dShift[iz_] * dx[iLev][iz_];
        const Real xi1 = dx[iLev][ix_] - xi0;
        const Real eta1 = dx[iLev][iy_] - eta0;
        const Real zeta1 = dx[iLev][iz_] - zeta0;

        weights_IIID[1][1][1][ix_] = eta0 * zeta0 * invVol[iLev];
        weights_IIID[1][1][1][iy_] = xi0 * zeta0 * invVol[iLev];
        weights_IIID[1][1][1][iz_] = xi0 * eta0 * invVol[iLev];

        // xi0*eta0*zeta1*invVol[iLev];
        weights_IIID[1][1][0][ix_] = eta0 * zeta1 * invVol[iLev];
        weights_IIID[1][1][0][iy_] = xi0 * zeta1 * invVol[iLev];
        weights_IIID[1][1][0][iz_] = -xi0 * eta0 * invVol[iLev];

        // xi0*eta1*zeta0*invVol[iLev];
        weights_IIID[1][0][1][ix_] = eta1 * zeta0 * invVol[iLev];
        weights_IIID[1][0][1][iy_] = -xi0 * zeta0 * invVol[iLev];
        weights_IIID[1][0][1][iz_] = xi0 * eta1 * invVol[iLev];

        // xi0*eta1*zeta1*invVol[iLev];
        weights_IIID[1][0][0][ix_] = eta1 * zeta1 * invVol[iLev];
        weights_IIID[1][0][0][iy_] = -xi0 * zeta1 * invVol[iLev];
        weights_IIID[1][0][0][iz_] = -xi0 * eta1 * invVol[iLev];

        // xi1*eta0*zeta0*invVol[iLev];
        weights_IIID[0][1][1][ix_] = -eta0 * zeta0 * invVol[iLev];
        weights_IIID[0][1][1][iy_] = xi1 * zeta0 * invVol[iLev];
        weights_IIID[0][1][1][iz_] = xi1 * eta0 * invVol[iLev];

        // xi1*eta0*zeta1*invVol[iLev];
        weights_IIID[0][1][0][ix_] = -eta0 * zeta1 * invVol[iLev];
        weights_IIID[0][1][0][iy_] = xi1 * zeta1 * invVol[iLev];
        weights_IIID[0][1][0][iz_] = -xi1 * eta0 * invVol[iLev];

        // xi1*eta1*zeta0*invVol[iLev];
        weights_IIID[0][0][1][ix_] = -eta1 * zeta0 * invVol[iLev];
        weights_IIID[0][0][1][iy_] = -xi1 * zeta0 * invVol[iLev];
        weights_IIID[0][0][1][iz_] = xi1 * eta1 * invVol[iLev];

        // xi1*eta1*zeta1*invVol[iLev];
        weights_IIID[0][0][0][ix_] = -eta1 * zeta1 * invVol[iLev];
        weights_IIID[0][0][0][iy_] = -xi1 * zeta1 * invVol[iLev];
        weights_IIID[0][0][0][iz_] = -xi1 * eta1 * invVol[iLev];

        const int iMin = loIdx[ix_];
        const int jMin = loIdx[iy_];
        const int kMin = loIdx[iz_];
        const int iMax = iMin + 1;
        const int jMax = jMin + 1;
        const int kMax = kMin + 1;

        const Real coef = fabs(qp) * invVol[iLev];
        RealVect wg_D;
        for (int k1 = kMin; k1 <= kMax; k1++)
          for (int j1 = jMin; j1 <= jMax; j1++)
            for (int i1 = iMin; i1 <= iMax; i1++) {

              for (int iDim = 0; iDim < nDim; iDim++) {
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
                  const int gp1 = gp0 + jp * nDim;
                  for (int k2 = kMin; k2 <= kMax; k2++) {
                    const Real(&wg1_D)[nDim] =
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

    }   // for p
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
std::array<Real, 5> Particles<NStructReal, NStructInt>::total_moments(
    bool localOnly) {
  timing_func("Particles::total_moments");

  std::array<Real, 5> sum = { 0, 0, 0, 0, 0 };

  for (int i = 0; i < 5; i++)
    sum[i] = 0;

  const int iLev = 0;
  for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    const auto& particles = pti.GetArrayOfStructs();
    for (const auto& p : particles) {
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      sum[0] += qp;
      sum[1] += qp * up;
      sum[2] += qp * vp;
      sum[3] += qp * wp;
      sum[4] += 0.5 * qp * (up * up + vp * vp + wp * wp);
    }
  }

  for (int i = 0; i < 5; i++)
    sum[i] *= qomSign * get_mass();

  if (!localOnly) {
    ParallelDescriptor::ReduceRealSum(sum.data(), sum.size(),
                                      ParallelDescriptor::IOProcessorNumber());
  }

  return sum;
}

//==========================================================
template <int NStructReal, int NStructInt>
Real Particles<NStructReal, NStructInt>::sum_moments(
    Vector<MultiFab>& momentsMF, Vector<MultiFab>& nodeBMF, Real dt) {
  timing_func("Particles::sum_moments");

  Real energy = 0;
  for (int iLev = 0; iLev <= finestLevel(); iLev++) {
    momentsMF[iLev].setVal(0.0);
    for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
         ++pti) {
      Array4<Real> const& momentsArr = momentsMF[iLev][pti].array();

      const auto& particles = pti.GetArrayOfStructs();

      // Print() << "iLev = " << iLev << std::endl;
      for (const auto& p : particles) {

        // Print() << "p = " << p << std::endl;
        const Real up = p.rdata(iup_);
        const Real vp = p.rdata(ivp_);
        const Real wp = p.rdata(iwp_);
        const Real qp = p.rdata(iqp_);

        //-----calculate interpolate coef begin-------------
        IntVect loIdx;
        RealVect dShift;
        for (int i = 0; i < nDim; i++) {
          dShift[i] = (p.pos(i) - plo[iLev][i]) * invDx[iLev][i];
          loIdx[i] = fastfloor(dShift[i]);
          dShift[i] = dShift[i] - loIdx[i];
        }

        Real coef[2][2][2];
        linear_interpolation_coef(dShift, coef);
        //-----calculate interpolate coef end-------------

        //-------nodePlasma begin---------
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

        for (int iVar = 0; iVar < nMoments; iVar++)
          for (int kk = 0; kk < 2; kk++)
            for (int jj = 0; jj < 2; jj++)
              for (int ii = 0; ii < 2; ii++) {
                const int i0 = loIdx[ix_] + ii;
                const int j0 = loIdx[iy_] + jj;
                const int k0 = loIdx[iz_] + kk;
                momentsArr(i0, j0, k0, iVar) +=
                    coef[ii][jj][kk] * pMoments[iVar];
              }

        //-------nodePlasma end---------

        energy += qp * (up * up + vp * vp + wp * wp);
      } // for p
    }

    // Exclude the number density.
    momentsMF[iLev].mult(invVol[iLev], 0, nMoments - 1,
                         momentsMF[iLev].nGrow());

    momentsMF[iLev].SumBoundary(Geom(iLev).periodicity());
  }

  for (int iLev = 0; iLev < finestLevel(); iLev++) {
    sum_two_lev_interface_nodes(momentsMF[iLev], momentsMF[iLev + 1], 0,
                                momentsMF[iLev].nComp(), get_ref_ratio(iLev));
  }

  energy *= 0.5 * qomSign * get_mass();

  return energy;
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_mass_matrix(
    UMultiFab<RealMM>& nodeMM, MultiFab& jHat, MultiFab& nodeBMF, Real dt) {
  timing_func("Particles::calc_mass_matrix");

  Real qdto2mc = charge / mass * 0.5 * dt;

  const int iLev = 0;
  for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    Array4<Real const> const& nodeBArr = nodeBMF[pti].array();
    Array4<Real> const& jArr = jHat[pti].array();
    Array4<RealMM> const& mmArr = nodeMM[pti].array();

    const auto& particles = pti.GetArrayOfStructs();

    for (const auto& p : particles) {

      // Print()<<"p = "<<p<<std::endl;
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      IntVect loIdx;
      RealVect dShift;
      for (int i = 0; i < nDim; i++) {
        dShift[i] = (p.pos(i) - plo[iLev][i]) * invDx[iLev][i];
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

      const Real c0 = denom * invVol[iLev] * qp * qdto2mc;
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

      {
        // jHat
        Real currents[3];

        {
          const Real coef1 = denom * qp;
          currents[ix_] = (up + (vp * Omz - wp * Omy + udotOm * Omx)) * coef1;
          currents[iy_] = (vp + (wp * Omx - up * Omz + udotOm * Omy)) * coef1;
          currents[iz_] = (wp + (up * Omy - vp * Omx + udotOm * Omz)) * coef1;
        }

        for (int iVar = 0; iVar < 3; iVar++)
          for (int kk = 0; kk < 2; kk++)
            for (int jj = 0; jj < 2; jj++)
              for (int ii = 0; ii < 2; ii++) {
                jArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk, iVar) +=
                    coef[ii][jj][kk] * currents[iVar];
              }
      }

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

                }   // j2
              }     // if (ip > 0)
            }       // i2
          }         // k1

      //----- Mass matrix calculation end--------------

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
          const int kp = 2;
          const int kr = k1 + kp - 1;
          if (kr > kMax || kr < kMin)
            continue;
          const Real* const datas0 = mmArr(i1, j1, k1).data;
          for (int jp = 0; jp < 3; jp++) {
            const int jr = j1 + jp - 1;
            if (jr > jMax || jr < jMin)
              continue;
            const int jpr = 2 - jp;
            for (int ip = 0; ip < 3; ip++) {
              const int ir = i1 + ip - 1;
              if (ir > iMax || ir < iMin)
                continue;
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
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_jhat(MultiFab& jHat,
                                                   MultiFab& nodeBMF, Real dt) {
  timing_func("Particles::calc_jhat");

  Real qdto2mc = charge / mass * 0.5 * dt;

  const int iLev = 0;
  for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    Array4<Real const> const& nodeBArr = nodeBMF[pti].array();
    Array4<Real> const& jArr = jHat[pti].array();

    const auto& particles = pti.GetArrayOfStructs();

    for (const auto& p : particles) {

      // Print()<<"p = "<<p<<std::endl;
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      IntVect loIdx;
      RealVect dShift;
      for (int i = 0; i < nDim; i++) {
        dShift[i] = (p.pos(i) - plo[iLev][i]) * invDx[iLev][i];
        loIdx[i] = fastfloor(dShift[i]);
        dShift[i] = dShift[i] - loIdx[i];
      }

      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

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

      {
        // jHat
        Real currents[3];

        {
          const Real coef1 = denom * qp;
          currents[ix_] = (up + (vp * Omz - wp * Omy + udotOm * Omx)) * coef1;
          currents[iy_] = (vp + (wp * Omx - up * Omz + udotOm * Omy)) * coef1;
          currents[iz_] = (wp + (up * Omy - vp * Omx + udotOm * Omz)) * coef1;
        }

        for (int iVar = 0; iVar < nDim; iVar++)
          for (int kk = 0; kk < 2; kk++)
            for (int jj = 0; jj < 2; jj++)
              for (int ii = 0; ii < 2; ii++) {
                jArr(loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk, iVar) +=
                    coef[ii][jj][kk] * currents[iVar];
              }
      }

    } // for p
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
Real Particles<NStructReal, NStructInt>::calc_max_thermal_velocity(
    MultiFab& momentsMF) {

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

  return uthMax;
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::convert_to_fluid_moments(
    Vector<MultiFab>& momentsMF) {

  for (int iLev = 0; iLev <= finestLevel(); iLev++) {
    MultiFab tmpMF(momentsMF[iLev], make_alias, iRho_, iPyz_ - iRho_ + 1);
    tmpMF.mult(qomSign * get_mass(), tmpMF.nGrow());

    for (MFIter mfi(momentsMF[iLev]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = momentsMF[iLev][mfi];
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
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::update_position_to_half_stage(
    const amrex::MultiFab& nodeEMF, const amrex::MultiFab& nodeBMF,
    amrex::Real dt) {
  timing_func("Particles::update_position_to_half_stage");

  Real dtLoc = 0.5 * dt;

  const int iLev = 0;
  for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    const Box& bx = cellStatus[iLev][pti].box();
    const IntVect lowCorner = bx.smallEnd();
    const IntVect highCorner = bx.bigEnd();

    const Array4<int const>& status = cellStatus[iLev][pti].array();

    auto& particles = pti.GetArrayOfStructs();
    for (auto& p : particles) {
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real xp = p.pos(ix_);
      const Real yp = p.pos(iy_);
      const Real zp = p.pos(iz_);

      p.pos(ix_) = xp + up * dtLoc;
      p.pos(iy_) = yp + vp * dtLoc;
      p.pos(iz_) = zp + wp * dtLoc;

      // Mark for deletion
      if (is_outside_ba(p, status, lowCorner, highCorner)) {
        p.id() = -1;
      }
    } // for p
  }   // for pti

  // This function distributes particles to proper processors and apply
  // periodic boundary conditions if needed.
  Redistribute();
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::mover(const Vector<MultiFab>& nodeE,
                                               const Vector<MultiFab>& nodeB,
                                               Real dt, Real dtNext) {
  if (is_neutral()) {
    neutral_mover(dt);
  } else {
    charged_particle_mover(nodeE, nodeB, dt, dtNext);
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::charged_particle_mover(
    const Vector<MultiFab>& nodeE, const Vector<MultiFab>& nodeB, Real dt,
    Real dtNext) {
  timing_func("Particles::charged_particle_mover");

  const Real qdto2mc = charge / mass * 0.5 * dt;
  Real dtLoc = 0.5 * (dt + dtNext);

  if (particlePosition == NonStaggered) {
    // Update location from x^{n+1/2} to x^{n+1} for nonstaggered case.
    dtLoc = 0.5 * dt;
  }

  for (int iLev = 0; iLev <= finestLevel(); iLev++) {
    for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
         ++pti) {
      const Array4<Real const>& nodeEArr = nodeE[iLev][pti].array();
      const Array4<Real const>& nodeBArr = nodeB[iLev][pti].array();

      const Array4<int const>& status = cellStatus[iLev][pti].array();
      // cellStatus[iLev][pti] is a FAB, and the box returned from the box()
      // method already contains the ghost cells.
      const Box& bx = cellStatus[iLev][pti].box();
      const IntVect lowCorner = bx.smallEnd();
      const IntVect highCorner = bx.bigEnd();

      auto& particles = pti.GetArrayOfStructs();
      for (auto& p : particles) {
        const Real up = p.rdata(iup_);
        const Real vp = p.rdata(ivp_);
        const Real wp = p.rdata(iwp_);
        const Real xp = p.pos(ix_);
        const Real yp = p.pos(iy_);
        const Real zp = p.pos(iz_);

        //-----calculate interpolate coef begin-------------
        IntVect loIdx;
        RealVect dShift;
        for (int i = 0; i < nDim; i++) {
          dShift[i] = (p.pos(i) - plo[iLev][i]) * invDx[iLev][i];
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
        if (is_outside_ba(p, status, lowCorner, highCorner)) {
          p.id() = -1;
        }
      } // for p
    }   // for pti
  }

  // This function distributes particles to proper processors and apply
  // periodic boundary conditions if needed.
  Redistribute();
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::neutral_mover(amrex::Real dt) {
  timing_func("Particles::neutral_mover");

  for (int iLev = 0; iLev <= finestLevel(); iLev++) {
    for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
         ++pti) {

      const Array4<int const>& status = cellStatus[iLev][pti].array();
      // cellStatus[iLev][pti] is a FAB, and the box returned from the box()
      // method already contains the ghost cells.
      const Box& bx = cellStatus[iLev][pti].box();
      const IntVect lowCorner = bx.smallEnd();
      const IntVect highCorner = bx.bigEnd();

      auto& particles = pti.GetArrayOfStructs();
      for (auto& p : particles) {
        const Real up = p.rdata(iup_);
        const Real vp = p.rdata(ivp_);
        const Real wp = p.rdata(iwp_);
        const Real xp = p.pos(ix_);
        const Real yp = p.pos(iy_);
        const Real zp = p.pos(iz_);

        p.pos(ix_) = xp + up * dt;
        p.pos(iy_) = yp + vp * dt;
        p.pos(iz_) = zp + wp * dt;

        // Mark for deletion
        if (is_outside_ba(p, status, lowCorner, highCorner)) {
          p.id() = -1;
        }
      } // for p
    }   // for pti
  }

  // This function distributes particles to proper processors and apply
  // periodic boundary conditions if needed.
  Redistribute();
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::divE_correct_position(
    const amrex::MultiFab& phiMF) {
  timing_func("Particles::divE_correct_position");

  const Real coef = charge / fabs(charge);
  const Real epsLimit = 0.1;
  Real epsMax = 0;

  const int iLev = 0;
  for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    Array4<Real const> const& phiArr = phiMF[pti].array();

    const Array4<int const>& status = cellStatus[iLev][pti].array();
    const Box& bx = cellStatus[iLev][pti].box();
    const IntVect lowCorner = bx.smallEnd();
    const IntVect highCorner = bx.bigEnd();

    auto& particles = pti.GetArrayOfStructs();

    for (auto& p : particles) {
      if (p.id() == -1 || is_outside_ba(p, status, lowCorner, highCorner)) {
        p.id() = -1;
        continue;
      }

      IntVect loIdx;
      RealVect dShift;
      for (int i = 0; i < nDim; i++) {
        // plo is the corner location => -0.5
        dShift[i] = (p.pos(i) - plo[iLev][i]) * invDx[iLev][i] - 0.5;
        loIdx[i] = fastfloor(dShift[i]);
        dShift[i] = dShift[i] - loIdx[i];
      }

      // Since the boundary condition for solving phi is not perfect,
      // correcting particles that are close to the boundaries may produce
      // artificial oscillations, which are seen in Earth's magnetotail
      // simulations. So, it is better to skip the boundary physical cells.
      bool isBoundaryPhysicalCell = false;
      for (int ix = 0; ix <= 1; ix++)
        for (int iy = 0; iy <= 1; iy++)
          for (int iz = 0; iz <= 1; iz++) {
            if (test_bit(
                    status(loIdx[ix_] + ix, loIdx[iy_] + iy, loIdx[iz_] + iz),
                    iDigitBny_))
              isBoundaryPhysicalCell = true;
          }
      if (isBoundaryPhysicalCell)
        continue;

      {
        Real weights_IIID[2][2][2][nDim];
        //----- Mass matrix calculation begin--------------
        const Real xi0 = dShift[ix_] * dx[iLev][ix_];
        const Real eta0 = dShift[iy_] * dx[iLev][iy_];
        const Real zeta0 = dShift[iz_] * dx[iLev][iz_];
        const Real xi1 = dx[iLev][ix_] - xi0;
        const Real eta1 = dx[iLev][iy_] - eta0;
        const Real zeta1 = dx[iLev][iz_] - zeta0;

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

        const int iMin = loIdx[ix_];
        const int jMin = loIdx[iy_];
        const int kMin = loIdx[iz_];

        RealVect eps_D = { AMREX_D_DECL(0, 0, 0) };

        for (int k = 0; k < 2; k++)
          for (int j = 0; j < 2; j++)
            for (int i = 0; i < 2; i++) {
              const Real coef = phiArr(iMin + i, jMin + j, kMin + k);
              for (int iDim = 0; iDim < nDim; iDim++) {
                eps_D[iDim] += coef * weights_IIID[i][j][k][iDim];
              }
            }

        for (int iDim = 0; iDim < nDim; iDim++)
          eps_D[iDim] *= coef * fourPI;

        if (fabs(eps_D[ix_] * invDx[iLev][ix_]) > epsLimit ||
            fabs(eps_D[iy_] * invDx[iLev][iy_]) > epsLimit ||
            fabs(eps_D[iz_] * invDx[iLev][iz_]) > epsLimit) {
          // If eps_D is too large, the underlying assumption of the particle
          // correction method will be not valid. Comparing each exp_D
          // component instead of the length dl saves the computational time.
          const Real dl = sqrt(pow(eps_D[ix_], 2) + pow(eps_D[iy_], 2) +
                               pow(eps_D[iz_], 2));
          const Real ratio = epsLimit * dx[iLev][ix_] / dl;
          for (int iDim = 0; iDim < nDim; iDim++)
            eps_D[iDim] *= ratio;
        }

        for (int iDim = 0; iDim < nDim; iDim++) {
          if (fabs(eps_D[iDim] * invDx[iLev][iDim]) > epsMax)
            epsMax = fabs(eps_D[iDim] * invDx[iLev][iDim]);

          p.pos(iDim) += eps_D[iDim];
        }

        if (is_outside_ba(p, status, lowCorner, highCorner)) {
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
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::split_particles(Real limit) {
  timing_func("Particles::split_particles");

  const int nInitial =
      nPartPerCell[ix_] * nPartPerCell[iy_] * nPartPerCell[iz_];
  const int nLowerLimit = nInitial * limit;

  const int nGoal = nLowerLimit > nInitial ? nLowerLimit : nInitial;

  IntVect iv = { AMREX_D_DECL(1, 1, 1) };
  if (!(do_tiling && tile_size == iv))
    return;

  for (int iLev = 0; iLev <= finestLevel(); iLev++) {

    const Real dl = 0.1 * Geom(iLev).CellSize()[ix_] / nPartPerCell.max();

    for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
         ++pti) {

      auto& particles = pti.GetArrayOfStructs();

      const int nPartOrig = particles.size();

      if (nPartOrig > nLowerLimit)
        continue;

      const int nNew =
          nGoal - nPartOrig > nPartOrig ? nPartOrig : nGoal - nPartOrig;

      // Find the 'heaviest' nNew particles by sorting the weight (charge).-----

      // Sort the particles by the location first to make sure the results
      // are the same for different number of processors
      std::sort(particles.begin(), particles.end(),
                [](const ParticleType& pl, const ParticleType& pr) {
                  return pl.pos(ix_) + pl.pos(iy_) + pl.pos(iz_) >
                         pr.pos(ix_) + pr.pos(iy_) + pr.pos(iz_);
                });

      // Sort the particles by the weight in decending order.
      std::sort(particles.begin(), particles.end(),
                [](const ParticleType& pl, const ParticleType& pr) {
                  const Real ql = fabs(pl.rdata(iqp_));
                  const Real qr = fabs(pr.rdata(iqp_));
                  // Q: Why use '1e-6*ql' instead of `0'?
                  // A: If it is sorted by 'ql > qr', and all the particles in
                  // this cell
                  //   have the same weight, the particles are essentially
                  //   sorted by the last digit, which is random. The threshold
                  //   '1e-6*ql' is introduced to avoid the randomness.
                  return ql - qr > 1e-6 * ql;
                });
      //----------------------------------------------------------------

      const auto lo = lbound(pti.tilebox());
      const auto hi = ubound(pti.tilebox());

      const Real xMin = Geom(iLev).LoEdge(lo.x, ix_) +
                        Geom(iLev).CellSize()[ix_] * 1e-10,
                 xMax = Geom(iLev).HiEdge(hi.x, ix_) -
                        Geom(iLev).CellSize()[ix_] * 1e-10;

      const Real yMin = Geom(iLev).LoEdge(lo.y, iy_) +
                        Geom(iLev).CellSize()[iy_] * 1e-10,
                 yMax = Geom(iLev).HiEdge(hi.y, iy_) -
                        Geom(iLev).CellSize()[iy_] * 1e-10;

      const Real zMin = Geom(iLev).LoEdge(lo.z, iz_) +
                        Geom(iLev).CellSize()[iz_] * 1e-10,
                 zMax = Geom(iLev).HiEdge(hi.z, iz_) -
                        Geom(iLev).CellSize()[iz_] * 1e-10;

      for (int ip = 0; ip < nNew; ip++) {
        auto& p = particles[ip];
        Real qp1 = p.rdata(iqp_);
        Real xp1 = p.pos(ix_);
        Real yp1 = p.pos(iy_);
        Real zp1 = p.pos(iz_);
        Real up1 = p.rdata(iup_);
        Real vp1 = p.rdata(ivp_);
        Real wp1 = p.rdata(iwp_);

        const Real u2 = up1 * up1 + vp1 * vp1 + wp1 * wp1;

        Real coef = (u2 < 1e-13) ? 0 : dl / sqrt(u2);
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
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::merge_particles(Real limit) {
  timing_func("Particles::merge_particles");
  IntVect iv = { AMREX_D_DECL(1, 1, 1) };
  if (!(do_tiling && tile_size == iv))
    return;

  const int nPartGoal =
      nPartPerCell[ix_] * nPartPerCell[iy_] * nPartPerCell[iz_] * limit;

  const int nPartCombine = 6, nPartNew = 5;

  const Real coefVel = 1, coefPos = 1;

  // The range of the velocity domain: [-r0,r0]*thermal_velocity+bulk_velocity
  const Real r0 = 1.0;

  int nAvailableCombines = 0, nEqs = 0, nSolved = 0;

  for (int iLev = 0; iLev <= finestLevel(); iLev++) {
    for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
         ++pti) {

      auto& particles = pti.GetArrayOfStructs();

      const int nPartOrig = particles.size();

      if (nPartOrig <= nPartGoal)
        continue;

      // Phase space cell number in one direction.
      // The const 0.8 is choosen by experience.
      int nCell = ceil(0.8 * r0 * pow(nPartOrig, 1. / nDim));

      if (nCell < 3)
        continue;

      if (fastMerge) {
        nCell = nCell > 5 ? 5 : nCell;
      }

      // Sort the particles by the location first to make sure the results
      // are the same for different number of processors
      std::sort(particles.begin(), particles.end(),
                [](const ParticleType& pl, const ParticleType& pr) {
                  return pl.pos(ix_) + pl.pos(iy_) + pl.pos(iz_) >
                         pr.pos(ix_) + pr.pos(iy_) + pr.pos(iz_);
                });

      // Print() << "nPartOrig = " << nPartOrig << " nPartGoal = " << nPartGoal
      //         << " nCell = " << nCell << std::endl;

      // One particle may belong to more than one velocity bins, but it can be
      // only merged at most once.
      std::vector<bool> merged;
      merged.resize(nPartOrig, false);

      //----------------------------------------------------------------
      // Estimate the bulk velocity and thermal velocity.
      Real uBulk[nDimVel] = { 0, 0, 0 };
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];
        for (int iDir = 0; iDir < 3; iDir++) {
          uBulk[iDir] += pcl.rdata(iDir);
        }
      }

      for (int iDir = 0; iDir < nDimVel; iDir++) {
        uBulk[iDir] /= nPartOrig;
      }

      Real thVel = 0, thVel2 = 0;
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];
        for (int iDir = 0; iDir < nDimVel; iDir++) {
          thVel2 += pow(pcl.rdata(iDir) - uBulk[iDir], 2);
        }
      }

      thVel2 /= nPartOrig;
      thVel = sqrt(thVel2);

      // The coef 0.5 if choosen by experience.
      const Real velNorm = (thVel < 1e-13) ? 0 : 1.0 / (0.5 * thVel);
      //----------------------------------------------------------------

      //----------------------------------------------------------------
      // Assign the particle IDs to the corresponding velocity space cells.
      Vector<int> phasePartIdx_III[nCell][nCell][nCell];

      // Velocity domain range.
      Real velMin_D[nDimVel], velMax_D[nDimVel];
      for (int iDir = 0; iDir < nDimVel; iDir++) {
        velMin_D[iDir] = -r0 * thVel + uBulk[iDir];
        velMax_D[iDir] = r0 * thVel + uBulk[iDir];
      }

      Real dv = (2.0 * r0 * thVel) / nCell;
      Real invDv = (dv < 1e-13) ? 0 : 1.0 / dv;

      int iCell_D[nDimVel];
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];

        bool isOutside = false;
        for (int iDim = 0; iDim < nDimVel; iDim++) {
          if (pcl.rdata(iDim) < velMin_D[iDim] ||
              pcl.rdata(iDim) > velMax_D[iDim])
            isOutside = true;
        }
        if (isOutside)
          continue;

        for (int iDim = 0; iDim < nDimVel; iDim++) {
          iCell_D[iDim] = fastfloor((pcl.rdata(iDim) - velMin_D[iDim]) * invDv);
        }

        // One particle may belong to multiple bins when each bin has a buffer
        // region.
        for (int xCell = iCell_D[ix_] - 1; xCell <= iCell_D[ix_] + 1; xCell++)
          for (int yCell = iCell_D[iy_] - 1; yCell <= iCell_D[iy_] + 1; yCell++)
            for (int zCell = iCell_D[iz_] - 1; zCell <= iCell_D[iz_] + 1;
                 zCell++) {

              if (xCell < 0 || xCell >= nCell || yCell < 0 || yCell >= nCell ||
                  zCell < 0 || zCell >= nCell)
                continue;

              IntVect cellIdx = { AMREX_D_DECL(xCell, yCell, zCell) };

              Real binMin_D[nDimVel], binMax_D[nDimVel];

              for (int iDim = 0; iDim < nDimVel; iDim++) {
                binMin_D[iDim] =
                    velMin_D[iDim] + (cellIdx[iDim] - velBinBufferSize) * dv;

                binMax_D[iDim] = velMin_D[iDim] +
                                 (cellIdx[iDim] + 1 + velBinBufferSize) * dv;
              }

              bool isInside = true;
              for (int iDim = 0; iDim < nDimVel; iDim++) {
                if (pcl.rdata(iDim) < binMin_D[iDim] ||
                    pcl.rdata(iDim) > binMax_D[iDim])
                  isInside = false;
              }

              if (isInside) {
                phasePartIdx_III[cellIdx[ix_]][cellIdx[iy_]][cellIdx[iz_]]
                    .push_back(pid);
              }
            }
      }
      //----------------------------------------------------------------

      for (int iu = 0; iu < nCell; iu++)
        for (int iv = 0; iv < nCell; iv++)
          for (int iw = 0; iw < nCell; iw++) {
            Vector<int> partIdx;
            for (int i = 0; i < phasePartIdx_III[iu][iv][iw].size(); i++) {
              auto& pIdx = phasePartIdx_III[iu][iv][iw];
              int pid = pIdx[i];
              if (!merged[pid]) {
                partIdx.push_back(pid);
              }
            }

            if (partIdx.size() < nPartCombine)
              continue;

            Vector<Real> distance;
            distance.resize(partIdx.size(), 0);

            // Find the center of the particles, and sort the particles based on
            // its distance to the 6-D center.
            //----------------------------------------------------------
            Real middle[6] = { 0, 0, 0, 0, 0, 0 };
            for (int pID : partIdx) {
              for (int iDir = ix_; iDir <= iz_; iDir++) {
                middle[iDir] += particles[pID].pos(iDir);
                middle[nDim + iDir] += particles[pID].rdata(iDir);
              }
            }

            for (int iDir = 0; iDir < 2 * nDim; iDir++) {
              middle[iDir] /= partIdx.size();
            }

            auto calc_distance2_to_center = [&, this](int pID) {
              Real dl2 = 0, dvel2 = 0;
              for (int iDir = ix_; iDir <= iz_; iDir++) {
                Real pos = particles[pID].pos(iDir);
                dl2 += pow((pos - middle[iDir]) * invDx[iLev][iDir], 2);

                Real v = particles[pID].rdata(iDir);
                dvel2 += pow((v - middle[nDim + iDir]) * velNorm, 2);
              }
              return coefPos * dl2 + coefVel * dvel2;
            };

            std::sort(partIdx.begin(), partIdx.end(),
                      [this, &particles, calc_distance2_to_center](
                          const int& idl, const int& idr) {
                        return calc_distance2_to_center(idl) <
                               calc_distance2_to_center(idr);
                      });
            //----------------------------------------------------------

            if (partIdx.size() >= nPartCombine) {
              /*
                  Delete 1 particle out of 6 particles:
                  1) Choose two particles that are closest to each other.
                  2) Delete the lighter one.
                  3) Distribute its weights to another 5 particles to conserve
                     mass, momentum and energy.
               */
              int idx_I[nPartCombine];
              for (int ip = 0; ip < nPartCombine; ip++) {
                idx_I[ip] = partIdx[ip];
              }

              // Calculate the center of the particles for combination
              for (int i = 0; i < 2 * nDim; i++) {
                middle[i] = 0;
              }
              for (int pID : idx_I) {
                for (int iDir = ix_; iDir <= iz_; iDir++) {
                  middle[iDir] += particles[pID].pos(iDir);
                  middle[nDim + iDir] += particles[pID].rdata(iDir);
                }
              }
              for (int iDir = 0; iDir < 2 * nDim; iDir++) {
                middle[iDir] /= nPartCombine;
              }

              bool doCombine = true;
              // Print()<<"--------------------"<<std::endl;
              for (int pID : idx_I) {
                Real distance = sqrt(calc_distance2_to_center(pID));
                // Print()<<"distance = "<<distance<<std::endl;
                if (distance > mergeThresholdDistance) {
                  doCombine = false;
                }
              }
              // Print()<<"--------------------"<<std::endl;
              nAvailableCombines++;
              if (!doCombine)
                continue;

              int nTry = 1;
              if (fastMerge) {
                nTry = nPartCombine;
              }

              int tryDeleteOrder[nPartCombine];
              for (int i = 0; i < nPartCombine; i++) {
                tryDeleteOrder[i] = idx_I[i];
              }

              if (nTry == 1) {
                // Find the pair that is closest to each other in phase space
                int pair1 = 0, pair2 = 0;
                Real dis2Min = 2;
                for (int ip1 = 0; ip1 < nPartCombine - 1; ip1++)
                  for (int ip2 = ip1 + 1; ip2 < nPartCombine; ip2++) {

                    // Distance between two particles in 6D space.
                    Real dl2 = 0, dv2 = 0;
                    for (int iDir = 0; iDir < nDim; iDir++) {
                      Real dv = velNorm * (particles[idx_I[ip1]].rdata(iDir) -
                                           particles[idx_I[ip2]].rdata(iDir));
                      dv2 += pow(dv, 2);

                      Real dx =
                          invDx[iLev][iDir] * (particles[idx_I[ip1]].pos(iDir) -
                                               particles[idx_I[ip2]].pos(iDir));
                      dv2 += pow(dx, 2);
                    }

                    const Real dis2 = dv2 * coefVel + dl2 * coefPos;

                    if (dis2 < dis2Min) {
                      dis2Min = dis2;
                      pair1 = ip1;
                      pair2 = ip2;
                    }
                  }
                //-------------------------------

                // Delete the lighter one.
                int iPartDel = pair1;
                // Q: Why is it 'l>(1+1e-6)*r' instead of 'l>r'?
                // A: The particle weights can be the same for some cases. 'l>r'
                // may return random results due to the truncation error.
                if (fabs(particles[idx_I[pair1]].rdata(iqp_)) >
                    (1 + 1e-6) * fabs(particles[idx_I[pair2]].rdata(iqp_))) {
                  iPartDel = pair2;
                }
                if (iPartDel != 0) {
                  tryDeleteOrder[0] = idx_I[iPartDel];
                  tryDeleteOrder[iPartDel] = idx_I[0];
                }
              } else {
                // Sort the particles by weights.
                std::sort(tryDeleteOrder, tryDeleteOrder + nPartCombine,
                          [&particles](const int& idLeft, const int& idRight) {
                            const Real ql = fabs(particles[idLeft].rdata(iqp_));
                            const Real qr =
                                fabs(particles[idRight].rdata(iqp_));
                            return ql - qr < 1e-6 * ql;
                          });
              }

              nEqs++;
              bool isSolved;
              const int nVar = nPartNew;
              Real x[nVar];

              for (int iPartDel = 0; iPartDel < nTry; iPartDel++) {
                for (int i = 0; i < nPartCombine; i++) {
                  idx_I[i] = tryDeleteOrder[i];
                }

                if (iPartDel != nPartCombine - 1) {
                  const int idxTmp = idx_I[iPartDel];
                  idx_I[iPartDel] = idx_I[nPartCombine - 1];
                  idx_I[nPartCombine - 1] = idxTmp;
                }

                //-----------Solve the new particle weights---------------
                // const int nVar = nPartNew;
                const int iq_ = 0, iu_ = 1, iv_ = 2, iw_ = 3, ie_ = 4;
                Real a[nVar][nVar + 1];
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

                //------------------------------------------
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
                //------------------------------------------

                isSolved = linear_solver_Gauss_Elimination();

                if (isSolved) {
                  // All the particle weights should have the same sign.
                  Real qt = a[iq_][nVar];
                  for (int ip = 0; ip < nPartNew; ip++) {
                    if (qt * x[ip] <= 0) {
                      isSolved = false;
                      break;
                    }
                  }
                }

                // Print() << " ipartdel = " << iPartDel << " w = "
                //         << particles[idx_I[nPartCombine - 1]].rdata(iqp_)
                //         << " isSolved = " << (isSolved ? 'T' : 'F') <<
                //         std::endl;
                if (isSolved) {
                  break;
                }
              }
              if (!isSolved)
                continue;
              //----------------------------------------------

              // Adjust weight.
              for (int ip = 0; ip < nPartNew; ip++) {
                auto& p = particles[idx_I[ip]];
                p.rdata(iqp_) = x[ip];
                merged[idx_I[ip]] = true;
              }
              // Mark for deletion
              for (int ip = nPartNew; ip < nPartCombine; ip++) {
                particles[idx_I[ip]].id() = -1;
                merged[idx_I[ip]] = true;
              }
              nSolved++;
            }
          }
    }
  }

  // if (nAvailableCombines > 0)
  //   AllPrint() << "Particle merging:: nAvailableCombine = "
  //              << nAvailableCombines << " nEquations = " << nEqs
  //              << " nSolved = " << nSolved << std::endl;
}

//==========================================================
template <int NStructReal, int NStructInt>
bool Particles<NStructReal, NStructInt>::do_inject_particles_for_this_cell(
    const amrex::Box& bx, const amrex::Array4<const int>& status, const int i,
    const int j, const int k, int& isrc, int& jsrc, int& ksrc) {

  // This cell should be a boundary cell at least.
  if (!test_bit(status(i, j, k), iDigitBny_))
    return false;

  for (int iloop = 1; iloop <= 3; iloop++) {
    // iloop==1: loop through faces;
    // iloop==2: loop through edges;
    // iloop==3: loop through corners;

    for (int di = -1; di <= 1; di++)
      for (int dj = -1; dj <= 1; dj++)
        for (int dk = -1; dk <= 1; dk++) {
          const int sum = std::abs(di) + std::abs(dj) + std::abs(dk);
          if (iloop != sum)
            continue;

          if (!test_bit(status(i + di, j + dj, k + dk), iDigitBny_)) {
            // The first neighbor cell that is NOT a boundary cell.
            if (bx.contains(IntVect{ AMREX_D_DECL(i + di, j + dj, k + dk) })) {
              isrc = i + di;
              jsrc = j + dj;
              ksrc = k + dk;
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

IOParticles::IOParticles(Particles& other, AmrCore* amrcore, Real no2outL,
                         Real no2outV, Real no2outM, RealBox IORange)
    : Particles(amrcore, nullptr, nullptr, other.get_speciesID(),
                other.get_charge(), other.get_mass(),
                amrex::IntVect(AMREX_D_DECL(-1, -1, -1))) {
  const int iLev = 0;
  no2outM *= qomSign * get_mass();

  const bool doLimit = IORange.ok();
  const auto& plevelOther = other.GetParticles(iLev);
  auto& plevel = GetParticles(iLev);
  for (MFIter mfi = other.MakeMFIter(iLev); mfi.isValid(); ++mfi) {
    auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());

    if (plevelOther.find(index) == plevelOther.end())
      continue;

    const auto& tileOther = plevelOther.at(index);

    if (tileOther.numParticles() == 0)
      continue;

    const auto& aosOther = tileOther.GetArrayOfStructs();

    const Box& bx = other.get_cell_status(iLev)[mfi].box();
    const IntVect lowCorner = bx.smallEnd();
    const IntVect highCorner = bx.bigEnd();
    const Array4<int const>& status = other.get_cell_status(iLev)[mfi].array();

    for (auto p : aosOther) {
      if (other.is_outside_ba(p, status, lowCorner, highCorner)) {
        // Redistribute() may fail if the ghost cell particles' IDs are not
        // -1 (marked for deletion);
        p.id() = -1;
      }

      for (int iDim = 0; iDim < nDim; iDim++) {
        p.pos(ix_ + iDim) = no2outL * p.pos(ix_ + iDim);
      }

      if (doLimit && !IORange.contains(RealVect(
                         AMREX_D_DECL(p.pos(ix_), p.pos(iy_), p.pos(iz_)))))
        continue;

      for (int iDim = 0; iDim < nDim; iDim++) {
        p.rdata(iup_ + iDim) = no2outV * p.rdata(iup_ + iDim);
      }
      p.rdata(iqp_) = no2outM * p.rdata(iqp_);

      plevel[index].push_back(p);
    }
  }
  Redistribute();
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::charge_exchange(
    Real dt, FluidInterface* stateOH, FluidInterface* sourcePT2OH,
    SourceInterface* source) {
  string nameFunc = "Particles::charge_exchange";

  timing_func(nameFunc);

  if (dt <= 0)
    return;

  const int iLev = 0;
  for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev); pti.isValid();
       ++pti) {
    auto& particles = pti.GetArrayOfStructs();
    for (auto& p : particles) {
      const Real xp = p.pos(ix_);
      const Real yp = p.pos(iy_);
      const Real zp = p.pos(iz_);

      double cs2Neu = 0, uNeu[3], rhoNeu;
      double cs2Ion, uIon[3], rhoIon;
      double ion2neu[5], neu2ion[5];
      const int iRho_ = 0, iUx_ = 1, iUy_ = 2, iUz_ = 3, iP_ = 4;
      const int iRhoUx_ = iUx_, iRhoUy_ = iUy_, iRhoUz_ = iUz_, iE_ = iP_;

      // amu/m^3
      rhoNeu = qomSign * p.rdata(iqp_) * get_mass() * invVol[iLev] *
               stateOH->get_No2SiRho() / cProtonMassSI;

      for (int i = 0; i < nDim; i++) {
        uNeu[i] = p.rdata(iup_ + i) * stateOH->get_No2SiV();
      }

      // MHD fluid index.
      const int fluidID = 0;
      // amu/m^3
      rhoIon = stateOH->get_fluid_mass_density(pti, xp, yp, zp, fluidID) *
               stateOH->get_No2SiRho() / cProtonMassSI;

      // cs = sqrt(P/n); m/s
      // Assume p = pi + pe = 2pi, so divide by sqrt(2.0).
      double cs = stateOH->get_fluid_uth(pti, xp, yp, zp, fluidID) *
                  stateOH->get_No2SiV() / sqrt(2.0);

      // cs2Ion = 2*P/n. The definition of thermal speed in get_uth_iso() is
      // different from the requirement in OH_get_charge_exchange_wrapper().
      // See page 92 of Adam Michael's thesis.
      cs2Ion = 2 * pow(cs, 2);

      uIon[ix_] = stateOH->get_fluid_ux(pti, xp, yp, zp, fluidID) *
                  stateOH->get_No2SiV();
      uIon[iy_] = stateOH->get_fluid_uy(pti, xp, yp, zp, fluidID) *
                  stateOH->get_No2SiV();
      uIon[iz_] = stateOH->get_fluid_uz(pti, xp, yp, zp, fluidID) *
                  stateOH->get_No2SiV();

      OH_get_charge_exchange_wrapper(&rhoIon, &cs2Ion, uIon, &rhoNeu, &cs2Neu,
                                     uNeu, ion2neu, neu2ion);

      // The function above returns number density changing rate.
      ion2neu[iRho_] *= cProtonMassSI;
      neu2ion[iRho_] *= cProtonMassSI;

      Real dtSI = dt * stateOH->get_No2SiT();
      // Print() << "rhoion = " << rhoIon << " cs2Ion = " << cs2Ion
      //         << " rhoNeu = " << rhoNeu << " cs2Neu = " << cs2Neu
      //         << " dtSI = " << dtSI << std::endl;
      for (int i = iRho_; i <= iP_; i++) {
        ion2neu[i] *= dtSI;
        neu2ion[i] *= dtSI;
        // Print() << " i = " << i << " ion2neu = " << ion2neu[i]
        //         << " neu2ion = " << neu2ion[i] << std::endl;
      }

      Real massExchange =
          neu2ion[iRho_] * stateOH->get_Si2NoRho() / invVol[iLev];

      if (massExchange == 0) {
        // It can happen for some special cases. For example, when neutral
        // density is zero.
        continue;
      }

      // Print() << "nden = " << p.rdata(iqp_)
      //         << " massExchange = " << massExchange << std::endl;
      if (p.rdata(iqp_) - massExchange <= 0) {
        // Mark for deletion
        p.id() = -1;

        // Reduce the sources accordingly to conserve total masses.
        const Real ratio = p.rdata(iqp_) / massExchange;
        for (int i = iRho_; i <= iP_; i++) {
          ion2neu[i] *= ratio;
          neu2ion[i] *= ratio;
        }
      } else {
        // Reduce particle mass due to charge exchange
        p.rdata(iqp_) = p.rdata(iqp_) - massExchange;
      }

      {
        // Q: Why is (neu2ion-ion2neu) divided by rhoIon?
        // A: What passed between PT and OH is 'source per ion density'
        // instead of source. The ion density will be multiplied back in OH
        // ModUser.f90
        sourcePT2OH->add_rho_to_loc((neu2ion[iRho_] - ion2neu[iRho_]) / rhoIon,
                                    pti, xp, yp, zp, fluidID);
        sourcePT2OH->add_mx_to_loc((neu2ion[iRhoUx_] - ion2neu[iRhoUx_]) /
                                       rhoIon,
                                   pti, xp, yp, zp, fluidID);
        sourcePT2OH->add_my_to_loc((neu2ion[iRhoUy_] - ion2neu[iRhoUy_]) /
                                       rhoIon,
                                   pti, xp, yp, zp, fluidID);
        sourcePT2OH->add_mz_to_loc((neu2ion[iRhoUz_] - ion2neu[iRhoUz_]) /
                                       rhoIon,
                                   pti, xp, yp, zp, fluidID);
        sourcePT2OH->add_p_to_loc((neu2ion[iP_] - ion2neu[iP_]) / rhoIon, pti,
                                  xp, yp, zp, fluidID);
      }

      if (ion2neu[iRho_] > 0) { // Add source to nodes.
        Real si2no_v[5];
        si2no_v[iRho_] = source->get_Si2NoRho();
        si2no_v[iRhoUx_] = source->get_Si2NoV() * si2no_v[iRho_];
        si2no_v[iRhoUy_] = si2no_v[iRhoUx_];
        si2no_v[iRhoUz_] = si2no_v[iRhoUx_];
        si2no_v[iP_] = source->get_Si2NoP();

        Real m2 = 0;
        for (int i = iRhoUx_; i <= iRhoUz_; i++) {
          m2 += pow(ion2neu[i], 2);
        }

        const Real gamma = 5. / 3;
        // P = (gamma-1)*(E - 0.5*rho*u2)
        ion2neu[iP_] = (gamma - 1) * (ion2neu[iE_] - 0.5 * m2 / ion2neu[iRho_]);

        if (ion2neu[iP_] < 0) {
          ion2neu[iP_] = 0;
        }

        // source saves changing rate (density/s...).
        source->add_rho_to_loc(ion2neu[iRho_] * si2no_v[iRho_] / dt, pti, xp,
                               yp, zp, speciesID);
        source->add_mx_to_loc(ion2neu[iRhoUx_] * si2no_v[iRhoUx_] / dt, pti, xp,
                              yp, zp, speciesID);
        source->add_my_to_loc(ion2neu[iRhoUy_] * si2no_v[iRhoUy_] / dt, pti, xp,
                              yp, zp, speciesID);
        source->add_mz_to_loc(ion2neu[iRhoUz_] * si2no_v[iRhoUz_] / dt, pti, xp,
                              yp, zp, speciesID);
        source->add_p_to_loc(ion2neu[iP_] * si2no_v[iP_] / dt, pti, xp, yp, zp,
                             speciesID);
      }
    } // for p
  }   // for pti
}

// Since Particles is a template, it is necessary to explicitly instantiate
// with template arguments.
template class Particles<nPicPartReal>;
template class Particles<nPTPartReal, nPTPartInt>;
