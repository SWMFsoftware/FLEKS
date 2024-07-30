#include <cstdlib>

#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

//==========================================================
template <int NStructReal, int NStructInt>
Particles<NStructReal, NStructInt>::Particles(
    Grid* gridIn, FluidInterface* const fluidIn, TimeCtr* const tcIn,
    const int speciesIDIn, const Real chargeIn, const Real massIn,
    const IntVect& nPartPerCellIn, const PartMode pModeIn, TestCase tcase)
    : AmrParticleContainer<NStructReal, NStructInt>(gridIn),
      grid(gridIn),
      fi(fluidIn),
      tc(tcIn),
      pMode(pModeIn),
      speciesID(speciesIDIn),
      charge(chargeIn),
      mass(massIn),
      nPartPerCell(nPartPerCellIn),
      testCase(tcase) {

  isParticleLocationRandom = gridIn->is_particle_location_random();
  do_tiling = true;

  qom = charge / mass;
  qomSign = qom >= 0 ? 1 : -1;

  plo.resize(n_lev_max());
  phi.resize(n_lev_max());
  dx.resize(n_lev_max());
  invDx.resize(n_lev_max());
  invVol.resize(n_lev_max());

  for (int iLev = 0; iLev < n_lev_max(); iLev++) {
    for (int i = 0; i < nDim; ++i) {
      tile_size[i] = 1;
      plo[iLev][i] = Geom(iLev).ProbLo(i);
      phi[iLev][i] = Geom(iLev).ProbHi(i);
      dx[iLev][i] = Geom(iLev).CellSize(i);
      invDx[iLev][i] = Geom(iLev).InvCellSize(i);
    }
    invVol[iLev] = invDx[iLev].product();
  }

  // The following line is used to avoid an MPI bug (feature?) on Frontera. It
  // should be removed after the bug being fixed.
  SetUseUnlink(false);
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::outflow_bc(const MFIter& mfi,
                                                    const IntVect ijkGst,
                                                    const IntVect ijkPhy) {
  const int iLev = 0;

  ParticleTileType& pGst = get_particle_tile(iLev, mfi, ijkGst);

  ParticleTileType& pPhy = get_particle_tile(iLev, mfi, ijkPhy);

  AoS& phyParts = pPhy.GetArrayOfStructs();

  RealVect dxshift;
  for (int i = 0; i < nDim; ++i) {
    dxshift[i] = Geom(iLev).CellSize(i) * (ijkGst[i] - ijkPhy[i]);
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
      ParticleType pNew = p;
      set_ids(pNew);
      for (int i = 0; i < nDim; ++i) {
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
    const int iLev, const MFIter& mfi, IntVect ijk,
    const FluidInterface* interface, bool doVacuumLimit, IntVect ppc,
    const Vel tpVel, Real dt) {

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

  // nPPC = nPPC / pow(2, iLev);

  set_random_seed(iLev, ijk, nPPC);

  const Real vol = dx[iLev].product();
  const int npcel = product(nPPC);

  const Real vol2Npcel = qomSign * vol / npcel;

  ParticleTileType& particles = get_particle_tile(iLev, mfi, ijk);

  int icount = 0;
  // Loop over particles inside grid cell i, j, k

  int kmax = 1;
  if (nDim > 2)
    kmax = nPPC[iz_];

  for (int ii = 0; ii < nPPC[ix_]; ++ii)
    for (int jj = 0; jj < nPPC[iy_]; ++jj)
      for (int kk = 0; kk < kmax; ++kk) {
        RealVect xyz, xyz0;

        IntVect ijk0 = { AMREX_D_DECL(ii, jj, kk) };

        for (int iDim = 0; iDim < nDim; iDim++) {
          xyz[iDim] = (ijk0[iDim] + randNum()) * (dx[iLev][iDim] / nPPC[iDim]) +
                      ijk[iDim] * dx[iLev][iDim] + plo[iLev][iDim];

          // If the particle weight is sampled in a random location, the sum of
          // particle mass is NOT the same as the integral of the grid density.
          // It is more convenient for debugging if mass is exactly conserved.
          // For a production run, it makes little difference.
          xyz0[iDim] = (ijk0[iDim] + 0.5) * (dx[iLev][iDim] / nPPC[iDim]) +
                       ijk[iDim] * dx[iLev][iDim] + plo[iLev][iDim];
        }

        if (nDim == 2) {
          // For comparison with the 3D case only.
          randNum();
        }

        if (!isParticleLocationRandom) {
          xyz = xyz0;
        }

        const double nDens =
            interface->get_number_density(mfi, xyz0, speciesID, iLev);

        if (doVacuumLimit && nDens * dt < vacuum)
          continue;

        double q = vol2Npcel * nDens;

        if (q != 0) {
          Real u, v, w;
          double rand1 = randNum();
          double rand2 = randNum();
          double rand3 = randNum();
          double rand4 = randNum();

          double uth = (userState ? tpVel.vth : -1);
          if (!is_neutral() && interface->get_UseAnisoP() &&
              (speciesID > 0 || interface->get_useElectronFluid())) {
            interface->set_particle_uth_aniso(iLev, mfi, xyz, &u, &v, &w, rand1,
                                              rand2, rand3, rand4, speciesID,
                                              uth, uth);
          } else {
            interface->set_particle_uth_iso(iLev, mfi, xyz, &u, &v, &w, rand1,
                                            rand2, rand3, rand4, speciesID,
                                            uth);
          }

          Real uBulk = userState ? tpVel.vx
                                 : interface->get_ux(mfi, xyz, speciesID, iLev);
          Real vBulk = userState ? tpVel.vy
                                 : interface->get_uy(mfi, xyz, speciesID, iLev);
          Real wBulk = userState ? tpVel.vz
                                 : interface->get_uz(mfi, xyz, speciesID, iLev);

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

          if (moveparticleswithconstantvelocity) {
            u = 0.0;
            v = 0.4;
            w = 0.0;
          }

          ParticleType p;
          set_ids(p);
          for (int iDim = 0; iDim < nDim; iDim++) {
            p.pos(iDim) = xyz[iDim];
          }
          p.rdata(iup_) = u;
          p.rdata(ivp_) = v;
          p.rdata(iwp_) = w;
          // Convert 'density changing rate' to 'density' if necessary.
          p.rdata(iqp_) = q * dt;

          if (NStructInt > iRecordCount_) {
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
    const FluidInterface* interface, const FluidInterface* const stateOH,
    Real dt, IntVect ppc, const bool doSelectRegion, const bool adaptivePPC) {
  timing_func("Pts::add_particles_source");

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
      const Box& tile_box = mfi.validbox();
      const auto lo = lbound(tile_box);
      const auto hi = ubound(tile_box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            const auto& status = cell_status(iLev)[mfi].array();
            if (bit::is_refined(status(i, j, k)))
              continue;

            IntVect ijk = { AMREX_D_DECL(i, j, k) };

            bool doAdd = true;
#ifdef _PT_COMPONENT_
            if (stateOH && doSelectRegion) {
              const int iFluid = 0;
              const int iRegion =
                  stateOH->get_neu_source_region(mfi, ijk, iFluid, iLev);
              doAdd = (iRegion == speciesID);
            }
#endif
            if (doAdd) {
              if (adaptivePPC) {
                // Adjust ppc so that the weight of the
                // source particles is not too small.
                const int initPPC = product(nPartPerCell);
                const int sourcePPC = product(ppc);

                Real rho = fi->get_number_density(mfi, ijk, speciesID, iLev);
                Real rhoSource =
                    interface->get_number_density(mfi, ijk, speciesID, iLev);
                if (dt > 0)
                  rhoSource *= dt;

                Real avgInitW = rho / initPPC;
                Real avgSourceW = rhoSource / sourcePPC;

                Real targetSourceW = avgInitW * 0.1;

                if (avgSourceW < targetSourceW) {
                  Real ratio = pow(avgSourceW / targetSourceW, 1.0 / nDim);
                  for (int iDim = 0; iDim < nDim; iDim++) {
                    ppc[iDim] = std::max(1, int(ppc[iDim] * ratio));
                  }
                }

                // Print() << "ppc = " << ppc << " rho = " << rho
                //         << " rhoSource = " << rhoSource
                //         << " rho/rhosource = " << (rho / rhoSource)
                //         << " dt = " << dt << std::endl;
              }

              add_particles_cell(iLev, mfi, ijk, interface, false, ppc, Vel(),
                                 dt);
            }
          }
    }
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_particles_domain() {
  timing_func("Pts::add_particles_domain");
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {

      const auto& status = cell_status(iLev)[mfi].array();
      ParallelFor(mfi.validbox(), [&](int i, int j, int k) noexcept {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };
        if (bit::is_new(status(ijk)) && !bit::is_refined(status(ijk))) {
          add_particles_cell(iLev, mfi, ijk, fi, true);
        }
      });
    }
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::inject_particles_at_boundary() {
  timing_func("Pts::inject_particles_at_boundary");

  // Only inject nGstInject layers.
  const int nGstInject = 1;

  // Only launch particles to the base grid boundary cells. The particle moments
  // of the domain edge nodes can be corrected by calling
  // interp_from_coarse_to_fine_for_domain_edge() in order from coarest level to
  // finest level.
  int iLev = 0;

  for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
    const auto& status = cell_status(iLev)[mfi].array();
    const Box& bx = mfi.validbox();
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    Box bxGst = bx;
    for (int iDim = 0; iDim < fi->get_fluid_dimension(); iDim++) {
      bxGst.grow(iDim, nGstInject);
    }

    ParallelFor(bxGst, [&](int i, int j, int k) noexcept {
      IntVect ijk = { AMREX_D_DECL(i, j, k) };
      IntVect ijksrc;
      if (do_inject_particles_for_this_cell(bx, status, ijk, ijksrc)) {
        if (((bc.lo[ix_] == bc.outflow) && i < lo.x) ||
            ((bc.hi[ix_] == bc.outflow) && i > hi.x) ||
            (nDim > 1 && (bc.lo[iy_] == bc.outflow) && j < lo.y) ||
            (nDim > 1 && (bc.hi[iy_] == bc.outflow) && j > hi.y) ||
            (nDim > 2 && (bc.lo[iz_] == bc.outflow) && k < lo.z) ||
            (nDim > 2 && (bc.hi[iz_] == bc.outflow) && k > hi.z)) {
          outflow_bc(mfi, ijk, ijksrc);
        } else if (((bc.lo[ix_] == bc.vacume) && i < lo.x) ||
                   ((bc.hi[ix_] == bc.vacume) && i > hi.x) ||
                   (nDim > 1 && (bc.lo[iy_] == bc.vacume) && j < lo.y) ||
                   (nDim > 1 && (bc.hi[iy_] == bc.vacume) && j > hi.y) ||
                   (nDim > 2 && (bc.lo[iz_] == bc.vacume) && k < lo.z) ||
                   (nDim > 2 && (bc.hi[iz_] == bc.vacume) && k > hi.z)) {
          // pass
        } else {
          add_particles_cell(iLev, mfi, ijk, fi, true, IntVect(), Vel(), -1);
        }
      }
    });
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::sum_to_center(
    MultiFab& netChargeMF, UMultiFab<RealCMM>& centerMM, bool doNetChargeOnly) {
  timing_func("Pts::sum_to_center");

  const int iLev = 0;
  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real> const& chargeArr = netChargeMF[pti].array();
    Array4<RealCMM> const& mmArr = centerMM[pti].array();

    const AoS& particles = pti.GetArrayOfStructs();

    const Dim3 lo = init_dim3(0);
    const Dim3 hi = init_dim3(1);

    for (const auto& p : particles) {
      /*
      Q: Why do not check p.id() < 0?
      A: IDs of ghost cell particles are set to -1 inside
      divE_correct_position(), but these particles should be take into account
      here.
      */

      // Print() << "particle = " << p << std::endl;

      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      IntVect loIdx;
      RealVect dShift;
      find_cell_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                      loIdx, dShift);
      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

      const Real cTmp = qp * invVol[iLev];
      for (int kk = lo.z; kk <= hi.z; ++kk)
        for (int jj = lo.y; jj <= hi.y; ++jj)
          for (int ii = lo.x; ii <= hi.x; ++ii) {
            const IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + ii, loIdx[iy_] + jj,
                                               loIdx[iz_] + kk) };
            chargeArr(ijk) += coef[ii][jj][kk] * cTmp;
          }

      if (!doNetChargeOnly) {
        Real weights_IIID[2][2][2][nDim3];
        //----- Mass matrix calculation begin--------------
        const Real xi0 = dShift[ix_] * dx[iLev][ix_];
        const Real eta0 = dShift[iy_] * dx[iLev][iy_];
        const Real zeta0 = nDim > 2 ? dShift[iz_] * dx[iLev][iz_] : 0;
        const Real xi1 = dx[iLev][ix_] - xi0;
        const Real eta1 = dx[iLev][iy_] - eta0;
        const Real zeta1 = nDim > 2 ? dx[iLev][iz_] - zeta0 : 1;

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
        const int kMin = nDim > 2 ? loIdx[iz_] : 0;
        const int iMax = iMin + 1;
        const int jMax = jMin + 1;
        const int kMax = nDim > 2 ? kMin + 1 : 0;

        const Real coef = fabs(qp) * invVol[iLev];
        RealVect wg_D;
        for (int k1 = kMin; k1 <= kMax; k1++)
          for (int j1 = jMin; j1 <= jMax; j1++)
            for (int i1 = iMin; i1 <= iMax; i1++) {

              for (int iDim = 0; iDim < nDim; iDim++) {
                wg_D[iDim] =
                    coef * weights_IIID[i1 - iMin][j1 - jMin][k1 - kMin][iDim];
              }

              auto& data = mmArr(i1, j1, k1);
              // Real weights[27] = { 0 };
              for (int i2 = iMin; i2 <= iMax; i2++) {
                int ip = i2 - i1 + 1;
                const int gp0 = ip * 9;
                for (int j2 = jMin; j2 <= jMax; j2++) {
                  int jp = j2 - j1 + 1;
                  const int gp1 = gp0 + jp * nDim3;
                  for (int k2 = kMin; k2 <= kMax; k2++) {
                    const Real(&wg1_D)[nDim3] =
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
template <int NStructReal, int NStructInt>
std::array<Real, 5> Particles<NStructReal, NStructInt>::total_moments(
    bool localOnly) {
  timing_func("Pts::total_moments");

  std::array<Real, 5> sum = { 0, 0, 0, 0, 0 };

  for (int i = 0; i < 5; ++i)
    sum[i] = 0;

  const int iLev = 0;
  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    const AoS& particles = pti.GetArrayOfStructs();
    for (const auto& p : particles) {
      if (p.id() < 0)
        continue;

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

  for (int i = 0; i < 5; ++i)
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
  timing_func("Pts::sum_moments");

  Real energy = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    timing_func("Pts::sum_moments_1");
    momentsMF[iLev].setVal(0.0);
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      Array4<Real> const& momentsArr = momentsMF[iLev][pti].array();

      const AoS& particles = pti.GetArrayOfStructs();

      const Dim3 lo = init_dim3(0);
      const Dim3 hi = init_dim3(1);

      // Print() << "iLev = " << iLev << std::endl;
      for (const auto& p : particles) {
        if (p.id() < 0)
          continue;

        // Print() << "p = " << p << std::endl;
        const Real up = p.rdata(iup_);
        const Real vp = p.rdata(ivp_);
        const Real wp = p.rdata(iwp_);
        const Real qp = p.rdata(iqp_);

        //-----calculate interpolate coef begin-------------
        IntVect loIdx;
        RealVect dShift;
        find_node_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                        loIdx, dShift);
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
          for (int kk = lo.z; kk <= hi.z; ++kk)
            for (int jj = lo.y; jj <= hi.y; ++jj)
              for (int ii = lo.x; ii <= hi.x; ++ii) {
                const IntVect ijk = { AMREX_D_DECL(
                    loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk) };
                momentsArr(ijk, iVar) += coef[ii][jj][kk] * pMoments[iVar];
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

  for (int iLev = n_lev() - 2; iLev >= 0; iLev--) {
    timing_func("Pts::sum_moments_2");
    sum_two_lev_interface_node(
        momentsMF[iLev], momentsMF[iLev + 1], 0, momentsMF[iLev].nComp(),
        get_ref_ratio(iLev), Geom(iLev), Geom(iLev + 1), node_status(iLev + 1));
  }

  // Correct domain edge nodes
  for (int iLev = 0; iLev < n_lev() - 1; iLev++) {
    timing_func("Pts::sum_moments_3");
    interp_from_coarse_to_fine_for_domain_edge(
        momentsMF[iLev], momentsMF[iLev + 1], 0, momentsMF[iLev].nComp(),
        get_ref_ratio(iLev), Geom(iLev), Geom(iLev + 1), node_status(iLev + 1));
  }

  energy *= 0.5 * qomSign * get_mass();

  return energy;
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_mass_matrix(
    UMultiFab<RealMM>& nodeMM, MultiFab& jHat, MultiFab& nodeBMF,
    MultiFab& u0MF, Real dt, int iLev, bool solveInCoMov) {
  timing_func("Pts::calc_mass_matrix");

  Real qdto2mc = charge / mass * 0.5 * dt;

  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const& nodeBArr = nodeBMF[pti].array();
    Array4<Real> const& jArr = jHat[pti].array();
    Array4<RealMM> const& mmArr = nodeMM[pti].array();

    Array4<Real const> const& u0Arr = u0MF[pti].array();

    const AoS& particles = pti.GetArrayOfStructs();

    const Dim3 lo = init_dim3(0);
    const Dim3 hi = init_dim3(1);

    for (const auto& p : particles) {
      if (p.id() < 0)
        continue;

      // Print()<<"p = "<<p<<std::endl;
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      IntVect loIdx;
      RealVect dShift;

      find_node_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                      loIdx, dShift);

      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

      //----- Mass matrix calculation begin--------------
      Real u0[3] = { 0, 0, 0 };
      Real bp[3] = { 0, 0, 0 };

      for (int kk = lo.z; kk <= hi.z; ++kk)
        for (int jj = lo.y; jj <= hi.y; ++jj)
          for (int ii = lo.x; ii <= hi.x; ++ii) {
            const IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + ii, loIdx[iy_] + jj,
                                               loIdx[iz_] + kk) };
            for (int iDim = 0; iDim < nDim3; iDim++) {
              bp[iDim] += nodeBArr(ijk, iDim) * coef[ii][jj][kk];

              if (solveInCoMov)
                u0[iDim] += u0Arr(ijk, iDim) * coef[ii][jj][kk];
            }
          }

      const Real Omx = qdto2mc * bp[ix_];
      const Real Omy = qdto2mc * bp[iy_];
      const Real Omz = qdto2mc * bp[iz_];

      // end interpolation
      const Real omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
      const Real denom = 1.0 / (1.0 + omsq);

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

        const Real up1 = up - u0[0];
        const Real vp1 = vp - u0[1];
        const Real wp1 = wp - u0[2];

        const Real udotOm1 = up1 * Omx + vp1 * Omy + wp1 * Omz;

        {
          const Real coef1 = denom * qp;
          currents[ix_] =
              (up1 + (vp1 * Omz - wp1 * Omy + udotOm1 * Omx)) * coef1;
          currents[iy_] =
              (vp1 + (wp1 * Omx - up1 * Omz + udotOm1 * Omy)) * coef1;
          currents[iz_] =
              (wp1 + (up1 * Omy - vp1 * Omx + udotOm1 * Omz)) * coef1;
        }

        for (int iVar = 0; iVar < 3; iVar++)
          for (int kk = lo.z; kk <= hi.z; ++kk)
            for (int jj = lo.y; jj <= hi.y; ++jj)
              for (int ii = lo.x; ii <= hi.x; ++ii) {
                IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + ii, loIdx[iy_] + jj,
                                             loIdx[iz_] + kk) };
                jArr(ijk, iVar) += coef[ii][jj][kk] * currents[iVar];
              }
      }

      const int iMin = loIdx[ix_];
      const int jMin = loIdx[iy_];
      const int kMin = nDim > 2 ? loIdx[iz_] : 0;
      const int iMax = iMin + 1;
      const int jMax = jMin + 1;
      const int kMax = nDim > 2 ? kMin + 1 : 0;

      for (int k1 = kMin; k1 <= kMax; k1++)
        for (int j1 = jMin; j1 <= jMax; j1++)
          for (int i1 = iMin; i1 <= iMax; i1++) {
            const Real wg = coef[i1 - iMin][j1 - jMin][k1 - kMin];
            auto& data0 = mmArr(i1, j1, k1);
            for (int k2 = kMin; k2 <= kMax; k2++) {
              const int kp = k2 - k1 + 1;
              if (kp > 0) {
                for (int j2 = jMin; j2 <= jMax; j2++) {
                  const int jp = j2 - j1 + 1;
                  for (int i2 = iMin; i2 <= iMax; i2++) {
                    const Real weight =
                        wg * coef[i2 - iMin][j2 - jMin][k2 - kMin];
                    const int idx0 = kp * 81 + jp * 27 + (i2 - i1 + 1) * 9;

                    Real* const data = &(data0[idx0]);
                    for (int idx = 0; idx < 9; idx++) {
                      data[idx] += alpha[idx] * weight;
                    }
                  } // k2

                } // j2
              } // if (ip > 0)
            } // i2
          } // k1

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
    const int iMin = lo.x - 1, jMin = lo.y - 1, kMin = nDim > 2 ? lo.z - 1 : 0;
    const int iMax = hi.x + 1, jMax = hi.y + 1, kMax = nDim > 2 ? hi.z + 1 : 0;

    int gps, gpr; // gp_send, gp_receive
    for (int k1 = kMin; k1 <= kMax; k1++)
      for (int j1 = jMin; j1 <= jMax; j1++)
        for (int i1 = iMin; i1 <= iMax; i1++) {
          const int kp = 2;
          const int kr = nDim > 2 ? k1 + kp - 1 : 0;
          if (kr > kMax || kr < kMin)
            continue;
          auto& datas0 = mmArr(i1, j1, k1);
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

              Real* const datar = &(mmArr(ir, jr, kr)[gpr * 9]);
              const Real* const datas = &(datas0[gps * 9]);
              for (int idx = 0; idx < 9; idx++) {
                datar[idx] = datas[idx];
              } // idx
            } // kp
          } // jp
        } // k1
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_mass_matrix_new(
    amrex::Vector<UMultiFab<RealMM> >& nodeMM, UMultiFab<RealMM>& nmmt,
    amrex::Vector<MultiFab>& jHat, MultiFab& jht,
    amrex::Vector<MultiFab>& nodeBMF, amrex::Vector<MultiFab>& u0MF, Real dt,
    int iLev, bool solveInCoMov, amrex::Vector<amrex::iMultiFab>& nodestatus) {
  timing_func("Pts::calc_mass_matrix");

  Real qdto2mc = charge / mass * 0.5 * dt;

  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const& nodeBArr = nodeBMF[iLev][pti].array();
    Array4<Real> const& jArr = jHat[iLev][pti].array();
    Array4<RealMM> const& mmArr = nodeMM[iLev][pti].array();
    Array4<Real const> const& u0Arr = u0MF[iLev][pti].array();

    const AoS& particles = pti.GetArrayOfStructs();

    const Dim3 lo = init_dim3(0);
    const Dim3 hi = init_dim3(1);

    /////////////////////////////////////////////////////////////////
    Array4<Real> const& jArrt = jht[pti].array();
    Array4<RealMM> const& mmArrt = nmmt[pti].array();
    /////////////////////////////////////////////////////////////////

    for (const auto& p : particles) {
      if (p.id() < 0)
        continue;

      // Print()<<"p = "<<p<<std::endl;
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      IntVect loIdx;
      IntVect cloIdx;
      IntVect floIdx;
      RealVect dShift;
      RealVect cdShift;
      RealVect fdShift;
      Real coef[2][2][2];
      Real coef_finer[2][2][2];
      Real coef_coarser[2][2][2];

      find_node_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                      loIdx, dShift);
      linear_interpolation_coef(dShift, coef);

      if (iLev > 0) {
        find_node_index(p.pos(), Geom(iLev - 1).ProbLo(),
                        Geom(iLev - 1).InvCellSize(), cloIdx, cdShift);
        linear_interpolation_coef(cdShift, coef_coarser);
      }
      if (iLev == 0) {
        find_node_index(p.pos(), Geom(iLev + 1).ProbLo(),
                        Geom(iLev + 1).InvCellSize(), floIdx, fdShift);
        linear_interpolation_coef(fdShift, coef_finer);
      }

      //-----calculate interpolate coef end-------------

      //----- Mass matrix calculation begin--------------
      Real u0[3] = { 0, 0, 0 };
      Real bp[3] = { 0, 0, 0 };

      for (int kk = lo.z; kk <= hi.z; ++kk)
        for (int jj = lo.y; jj <= hi.y; ++jj)
          for (int ii = lo.x; ii <= hi.x; ++ii) {
            const IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + ii, loIdx[iy_] + jj,
                                               loIdx[iz_] + kk) };
            for (int iDim = 0; iDim < nDim3; iDim++) {
              bp[iDim] += nodeBArr(ijk, iDim) * coef[ii][jj][kk];

              if (solveInCoMov)
                u0[iDim] += u0Arr(ijk, iDim) * coef[ii][jj][kk];
            }
          }

      const Real Omx = qdto2mc * bp[ix_];
      const Real Omy = qdto2mc * bp[iy_];
      const Real Omz = qdto2mc * bp[iz_];

      // end interpolation
      const Real omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
      const Real denom = 1.0 / (1.0 + omsq);

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

        const Real up1 = up - u0[0];
        const Real vp1 = vp - u0[1];
        const Real wp1 = wp - u0[2];

        const Real udotOm1 = up1 * Omx + vp1 * Omy + wp1 * Omz;

        {
          const Real coef1 = denom * qp;
          currents[ix_] =
              (up1 + (vp1 * Omz - wp1 * Omy + udotOm1 * Omx)) * coef1;
          currents[iy_] =
              (vp1 + (wp1 * Omx - up1 * Omz + udotOm1 * Omy)) * coef1;
          currents[iz_] =
              (wp1 + (up1 * Omy - vp1 * Omx + udotOm1 * Omz)) * coef1;
        }

        for (int iVar = 0; iVar < 3; iVar++)
          for (int kk = lo.z; kk <= hi.z; ++kk)
            for (int jj = lo.y; jj <= hi.y; ++jj)
              for (int ii = lo.x; ii <= hi.x; ++ii) {
                IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + ii, loIdx[iy_] + jj,
                                             loIdx[iz_] + kk) };
                jArr(ijk, iVar) += coef[ii][jj][kk] * currents[iVar];

                if (iLev > 0) {
                  ijk = { AMREX_D_DECL(cloIdx[ix_] + ii, cloIdx[iy_] + jj,
                                       cloIdx[iz_] + kk) };
                  jArrt(ijk, iVar) += coef_coarser[ii][jj][kk] * currents[iVar];
                }
                if (iLev == 0) {
                  ijk = { AMREX_D_DECL(floIdx[ix_] + ii, floIdx[iy_] + jj,
                                       floIdx[iz_] + kk) };
                  jArrt(ijk, iVar) += coef_finer[ii][jj][kk] * currents[iVar];
                }
              }
      }

      const int iMin = loIdx[ix_];
      const int jMin = loIdx[iy_];
      const int kMin = nDim > 2 ? loIdx[iz_] : 0;
      const int iMax = iMin + 1;
      const int jMax = jMin + 1;
      const int kMax = nDim > 2 ? kMin + 1 : 0;

      for (int k1 = kMin; k1 <= kMax; k1++)
        for (int j1 = jMin; j1 <= jMax; j1++)
          for (int i1 = iMin; i1 <= iMax; i1++) {
            const Real wg = coef[i1 - iMin][j1 - jMin][k1 - kMin];
            auto& data0 = mmArr(i1, j1, k1);
            for (int k2 = kMin; k2 <= kMax; k2++) {
              const int kp = k2 - k1 + 1;
              if (kp > 0) {
                for (int j2 = jMin; j2 <= jMax; j2++) {
                  const int jp = j2 - j1 + 1;
                  for (int i2 = iMin; i2 <= iMax; i2++) {
                    const Real weight =
                        wg * coef[i2 - iMin][j2 - jMin][k2 - kMin];
                    const int idx0 = kp * 81 + jp * 27 + (i2 - i1 + 1) * 9;

                    Real* const data = &(data0[idx0]);
                    for (int idx = 0; idx < 9; idx++) {
                      data[idx] += alpha[idx] * weight;
                    }
                  } // k2

                } // j2
              } // if (ip > 0)
            } // i2
          } // k1

      if (iLev > 0) {
        const int iMin = cloIdx[ix_];
        const int jMin = cloIdx[iy_];
        const int kMin = nDim > 2 ? cloIdx[iz_] : 0;
        const int iMax = iMin + 1;
        const int jMax = jMin + 1;
        const int kMax = nDim > 2 ? kMin + 1 : 0;

        for (int k1 = kMin; k1 <= kMax; k1++)
          for (int j1 = jMin; j1 <= jMax; j1++)
            for (int i1 = iMin; i1 <= iMax; i1++) {
              const Real wg = coef_coarser[i1 - iMin][j1 - jMin][k1 - kMin];
              auto& data0 = mmArrt(i1, j1, k1);
              for (int k2 = kMin; k2 <= kMax; k2++) {
                const int kp = k2 - k1 + 1;
                if (kp > 0) {
                  for (int j2 = jMin; j2 <= jMax; j2++) {
                    const int jp = j2 - j1 + 1;
                    for (int i2 = iMin; i2 <= iMax; i2++) {
                      const Real weight =
                          wg * coef_coarser[i2 - iMin][j2 - jMin][k2 - kMin];
                      const int idx0 = kp * 81 + jp * 27 + (i2 - i1 + 1) * 9;

                      Real* const data = &(data0[idx0]);
                      for (int idx = 0; idx < 9; idx++) {
                        data[idx] += alpha[idx] *
                                     (invVol[iLev - 1] / invVol[iLev]) * weight;
                      }
                    } // k2

                  } // j2
                } // if (ip > 0)
              } // i2
            } // k1
      }

      if (iLev == 0) {
        const int iMin = floIdx[ix_];
        const int jMin = floIdx[iy_];
        const int kMin = nDim > 2 ? floIdx[iz_] : 0;
        const int iMax = iMin + 1;
        const int jMax = jMin + 1;
        const int kMax = nDim > 2 ? kMin + 1 : 0;

        for (int k1 = kMin; k1 <= kMax; k1++)
          for (int j1 = jMin; j1 <= jMax; j1++)
            for (int i1 = iMin; i1 <= iMax; i1++) {
              const Real wg = coef_finer[i1 - iMin][j1 - jMin][k1 - kMin];
              auto& data0 = mmArrt(i1, j1, k1);
              for (int k2 = kMin; k2 <= kMax; k2++) {
                const int kp = k2 - k1 + 1;
                if (kp > 0) {
                  for (int j2 = jMin; j2 <= jMax; j2++) {
                    const int jp = j2 - j1 + 1;
                    for (int i2 = iMin; i2 <= iMax; i2++) {
                      const Real weight =
                          wg * coef_finer[i2 - iMin][j2 - jMin][k2 - kMin];
                      const int idx0 = kp * 81 + jp * 27 + (i2 - i1 + 1) * 9;

                      Real* const data = &(data0[idx0]);
                      for (int idx = 0; idx < 9; idx++) {
                        data[idx] += alpha[idx] *
                                     (invVol[iLev + 1] / invVol[iLev]) * weight;
                      }
                    } // k2

                  } // j2
                } // if (ip > 0)
              } // i2
            } // k1
      }

      //----- Mass matrix calculation end--------------

    } // for p
  }

  for (MFIter mfi(nodeMM[iLev]); mfi.isValid(); ++mfi) {
    // Finalize the mass matrix calculation.
    const Box box = mfi.validbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Array4<RealMM> const& mmArr = nodeMM[iLev][mfi].array();

    // We only need the mass matrix on the physical nodes. But the first layer
    // of the ghost nodes may contributes to the physical nodes below (ghost
    // node constributes as a sender). So, we need the '-1' and '+1' staff.
    const int iMin = lo.x - 1, jMin = lo.y - 1, kMin = nDim > 2 ? lo.z - 1 : 0;
    const int iMax = hi.x + 1, jMax = hi.y + 1, kMax = nDim > 2 ? hi.z + 1 : 0;

    int gps, gpr; // gp_send, gp_receive
    for (int k1 = kMin; k1 <= kMax; k1++)
      for (int j1 = jMin; j1 <= jMax; j1++)
        for (int i1 = iMin; i1 <= iMax; i1++) {
          const int kp = 2;
          const int kr = nDim > 2 ? k1 + kp - 1 : 0;
          if (kr > kMax || kr < kMin)
            continue;
          auto& datas0 = mmArr(i1, j1, k1);
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

              Real* const datar = &(mmArr(ir, jr, kr)[gpr * 9]);
              const Real* const datas = &(datas0[gps * 9]);
              for (int idx = 0; idx < 9; idx++) {
                datar[idx] = datas[idx];
              } // idx
            } // kp
          } // jp
        } // k1
  }
}
//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calc_jhat(MultiFab& jHat,
                                                   MultiFab& nodeBMF, Real dt) {
  timing_func("Pts::calc_jhat");

  Real qdto2mc = charge / mass * 0.5 * dt;

  const int iLev = 0;
  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const& nodeBArr = nodeBMF[pti].array();
    Array4<Real> const& jArr = jHat[pti].array();

    const AoS& particles = pti.GetArrayOfStructs();

    for (const auto& p : particles) {
      if (p.id() < 0)
        continue;

      // Print()<<"p = "<<p<<std::endl;
      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real qp = p.rdata(iqp_);

      //-----calculate interpolate coef begin-------------
      IntVect loIdx;
      RealVect dShift;

      find_node_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                      loIdx, dShift);

      Real coef[2][2][2];
      linear_interpolation_coef(dShift, coef);
      //-----calculate interpolate coef end-------------

      Real Bxl = 0, Byl = 0, Bzl = 0; // should be bp[3];

      for (int kk = 0; kk < 2; ++kk)
        for (int jj = 0; jj < 2; ++jj)
          for (int ii = 0; ii < 2; ++ii) {
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
          for (int kk = 0; kk < 2; ++kk)
            for (int jj = 0; jj < 2; ++jj)
              for (int ii = 0; ii < 2; ++ii) {
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

    ParallelFor(box, [&](int i, int j, int k) {
      Real rho = arr(i, j, k, iRho_);
      if (rho == 0)
        return;

      Real p =
          (arr(i, j, k, iPxx_) + arr(i, j, k, iPyy_) + arr(i, j, k, iPzz_)) *
          c1over3;

      Real uth = sqrt(p / rho);
      if (uth > uthMax)
        uthMax = uth;
    });
  }

  return uthMax;
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::convert_to_fluid_moments(
    Vector<MultiFab>& momentsMF) {

  for (int iLev = 0; iLev < n_lev(); iLev++) {
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
      if (is_outside_active_region(p, status, lowCorner, highCorner)) {
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
                                               Real dt, Real dtNext,
                                               bool solveInCoMov) {
  if (is_neutral()) {
    neutral_mover(dt);
  } else {
    charged_particle_mover(nodeE, nodeB, eBg, uBg, dt, dtNext, solveInCoMov);
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::charged_particle_mover(
    const Vector<MultiFab>& nodeE, const Vector<MultiFab>& nodeB,
    const Vector<MultiFab>& eBg, const Vector<MultiFab>& uBg, Real dt,
    Real dtNext, bool solveInCoMov) {
  timing_func("Pts::charged_particle_mover");

  const Real qdto2mc = charge / mass * 0.5 * dt;
  Real dtLoc = 0.5 * (dt + dtNext);

  if (particlePosition == NonStaggered) {
    // Update location from x^{n+1/2} to x^{n+1} for nonstaggered case.
    dtLoc = 0.5 * dt;
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      const Array4<Real const>& nodeEArr = nodeE[iLev][pti].array();
      const Array4<Real const>& nodeBArr = nodeB[iLev][pti].array();

      const Array4<Real const>& E0Arr = eBg[iLev][pti].array();
      const Array4<Real const>& U0Arr = uBg[iLev][pti].array();

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

                if (solveInCoMov) {
                  ep[iDim] += -E0Arr(ijk, iDim) * c0;
                  u0p[iDim] += U0Arr(ijk, iDim) * c0;
                }
              }
            }

        up = up - u0p[ix_];
        vp = vp - u0p[iy_];
        wp = wp - u0p[iz_];

        const double Omx = qdto2mc * bp[ix_];
        const double Omy = qdto2mc * bp[iy_];
        const double Omz = qdto2mc * bp[iz_];

        // end interpolation
        const Real omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
        const Real denom = 1.0 / (1.0 + omsq);
        // solve the position equation
        const Real ut = up + qdto2mc * ep[ix_];
        const Real vt = vp + qdto2mc * ep[iy_];
        const Real wt = wp + qdto2mc * ep[iz_];
        // const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
        const Real udotOm = ut * Omx + vt * Omy + wt * Omz;
        // solve the velocity equation
        const Real uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
        const Real vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
        const Real wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;

        double unp1 = 2.0 * uavg - up + u0p[ix_];
        double vnp1 = 2.0 * vavg - vp + u0p[iy_];
        double wnp1 = 2.0 * wavg - wp + u0p[iz_];

        if (moveparticleswithconstantvelocity) {
          unp1 = up;
          vnp1 = vp;
          wnp1 = wp;
        }
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

        // Mark for deletion
        if (is_outside_active_region(p, status, lowCorner, highCorner)) {
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

        // Mark for deletion
        if (is_outside_active_region(p, status, lowCorner, highCorner)) {
          p.id() = -1;
        }
      } // for p
    } // for pti
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::divE_correct_position(
    const MultiFab& phiMF) {
  timing_func("Pts:divE_correct_position");

  const Real coef = charge / fabs(charge);
  const Real epsLimit = 0.1;
  Real epsMax = 0;

  const int iLev = 0;
  for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
    Array4<Real const> const& phiArr = phiMF[pti].array();

    const Array4<int const>& status = cell_status(iLev)[pti].array();

    AoS& particles = pti.GetArrayOfStructs();

    const Box& bx = cell_status(iLev)[pti].box();
    const IntVect lowCorner = bx.smallEnd();
    const IntVect highCorner = bx.bigEnd();

    for (auto& p : particles) {
      if (p.id() == -1 ||
          is_outside_active_region(p, status, lowCorner, highCorner)) {
        p.id() = -1;
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
      if (isBoundaryPhysicalCell)
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

        const int iMin = loIdx[ix_];
        const int jMin = loIdx[iy_];
        const int kMin = nDim > 2 ? loIdx[iz_] : 0;

        RealVect eps_D = { AMREX_D_DECL(0, 0, 0) };

        for (int k = 0; k < 2; ++k)
          for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i) {
              IntVect ijk = { AMREX_D_DECL(iMin + i, jMin + j, kMin + k) };
              const Real coef = phiArr(ijk);
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

        if (is_outside_active_region(p, status, lowCorner, highCorner)) {
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
void Particles<NStructReal, NStructInt>::limit_weight(Real maxRatio,
                                                      bool seperateVelocity) {
  timing_func("Pts::limit_weight");

  if (maxRatio <= 1)
    return;

  IntVect iv(1);
  if (!(do_tiling && tile_size == iv))
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {

      Vector<ParticleType> newparticles;

      auto& pTile = get_particle_tile(iLev, pti);
      AoS& particles = pti.GetArrayOfStructs();

      // Sort the particles first to make sure the results
      // are the same for different number of processors
      std::sort(particles.begin(), particles.end(), compare_two_parts);

      Real totalMass = 0;
      Real totalMoment[nDim3] = { 0, 0, 0 };
      for (auto& p : particles) {
        totalMass += fabs(p.rdata(iqp_));
        for (int i = 0; i < nDim3; ++i)
          totalMoment[i] += fabs(p.rdata(iqp_) * p.rdata(iup_ + i));
      }
      Real avg = totalMass / particles.size();

      // Real maxWeight = avg + maxRatio * vars;
      Real maxWeight = avg * maxRatio;

      if (seperateVelocity) {
        Box bx = pti.tilebox();
        set_random_seed(iLev, bx.smallEnd(), IntVect(444));
        Vector<ParticleType*> pold;
        for (size_t ip = 0; ip < particles.size(); ip++) {
          Real qp1 = particles[ip].rdata(iqp_);
          if (fabs(qp1) < maxWeight)
            continue;
          pold.push_back(&(particles[ip]));
        }
        split_particles_by_velocity(pold, newparticles);
      } else {

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

        if (is_neutral()) {
          Box bx = pti.tilebox();
          set_random_seed(iLev, bx.smallEnd(), IntVect(999));
        }

        for (auto& p : particles) {
          Real qp1 = p.rdata(iqp_);
          if (fabs(qp1) < maxWeight)
            continue;

          Real up1 = p.rdata(iup_);
          Real vp1 = p.rdata(ivp_);
          Real wp1 = p.rdata(iwp_);

          Real xp1 = p.pos(ix_);
          Real yp1 = p.pos(iy_);
          Real zp1 = p.pos(iz_);

          int nNew = is_neutral() ? 7 : 1;

          p.rdata(iqp_) = qp1 / (nNew + 1);

          for (int iNew = 0; iNew < nNew; iNew++) {
            ParticleType pnew;
            set_ids(pnew);

            Real xp2 = xp1 + (xMax - xMin) * (randNum() - 0.5);
            Real yp2 = yp1 + (yMax - yMin) * (randNum() - 0.5);
            Real zp2 = zp1 + (zMax - zMin) * (randNum() - 0.5);

            xp2 = bound(xp2, xMin, xMax);
            yp2 = bound(yp2, yMin, yMax);
            zp2 = bound(zp2, zMin, zMax);

            pnew.pos(ix_) = xp2;
            pnew.pos(iy_) = yp2;
            pnew.pos(iz_) = zp2;
            pnew.rdata(iup_) = up1;
            pnew.rdata(ivp_) = vp1;
            pnew.rdata(iwp_) = wp1;

            pnew.rdata(iqp_) = qp1 / (nNew + 1);
            newparticles.push_back(pnew);
          }
        }
      }

      for (auto& p : newparticles) {
        pTile.push_back(p);
      }
    }
  }
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::split_particles_by_velocity(
    Vector<ParticleType*>& plist, Vector<ParticleType>& newparticles) {

  if (plist.size() < 2)
    return;

  const int nCell = pow(2, 3);
  // Assign the particle IDs to the corresponding velocity space cells.
  Vector<int> phasePartIdx_III[nCell][nCell][nCell];

  // Velocity domain range.
  Real velMin_D[nDim3] = { 1e9, 1e9, 1e9 },
       velMax_D[nDim3] = { -1e9, -1e9, -1e9 };

  for (int pid = 0; pid < plist.size(); pid++) {
    auto& pcl = *plist[pid];
    for (int iDir = 0; iDir < nDim3; iDir++) {
      if (pcl.rdata(iup_ + iDir) > velMax_D[iDir])
        velMax_D[iDir] = pcl.rdata(iup_ + iDir);
      if (pcl.rdata(iup_ + iDir) < velMin_D[iDir])
        velMin_D[iDir] = pcl.rdata(iup_ + iDir);
    }
  }

  Real dvMax = 0;
  for (int iDir = 0; iDir < nDim3; iDir++) {
    Real dv = velMax_D[iDir] - velMin_D[iDir];
    const Real vref = 0.5 * (fabs(velMax_D[iDir]) + fabs(velMin_D[iDir]));
    if (dv < 1e-6 * vref)
      dv = 1e-6 * vref;
    velMax_D[iDir] += 1e-3 * dv;
    velMin_D[iDir] -= 1e-3 * dv;
    dv = velMax_D[iDir] - velMin_D[iDir];

    if (dv > dvMax)
      dvMax = dv;
  }

  Real dvCell = dvMax == 0 ? 1e-9 : dvMax / nCell;
  Real invDv = 1 / dvCell;

  int iCell_D[nDim3];
  for (int pid = 0; pid < plist.size(); pid++) {
    auto& pcl = *plist[pid];
    for (int iDim = 0; iDim < nDim3; iDim++) {
      iCell_D[iDim] = fastfloor((pcl.rdata(iDim) - velMin_D[iDim]) * invDv);
    }

    phasePartIdx_III[iCell_D[ix_]][iCell_D[iy_]][iCell_D[iz_]].push_back(pid);
  }

  Vector<std::array<int, 3> > morton_idx(pow(nCell, 3));

  for (int iu = 0; iu < nCell; iu++)
    for (int iv = 0; iv < nCell; iv++)
      for (int iw = 0; iw < nCell; iw++) {
        morton_idx[encode_morton_3d(iu, iv, iw)] = { iu, iv, iw };
      }

  Vector<ParticleType*> p_morton;

  for (int i = 0; i < morton_idx.size(); ++i) {
    int iu = morton_idx[i][0];
    int iv = morton_idx[i][1];
    int iw = morton_idx[i][2];

    // printf("1 iu = %d iv = %d iw = %d\n", iu, iv, iw);
    for (int ip = 0; ip < phasePartIdx_III[iu][iv][iw].size(); ip++) {
      p_morton.push_back(plist[phasePartIdx_III[iu][iv][iw][ip]]);
    }
  }

  int nPair = floor(p_morton.size() / 2.0);
  for (int ip = 0; ip < nPair * 2; ip += 2) {
    ParticleType& p1 = *p_morton[ip];
    ParticleType& p2 = *p_morton[ip + 1];
    ParticleType p3, p4;

    Real du = p1.rdata(iup_) - p2.rdata(iup_);
    Real dv = p1.rdata(ivp_) - p2.rdata(ivp_);
    Real dw = p1.rdata(iwp_) - p2.rdata(iwp_);
    Real dspeed = sqrt(du * du + dv * dv + dw * dw);

    if (dspeed / dvCell > 2)
      continue;

    bool doSucceed = split_by_seperate_velocity(p1, p2, p3, p4);
    if (doSucceed) {
      newparticles.push_back(p3);
      newparticles.push_back(p4);
    }
  }
}

template <int NStructReal, int NStructInt>
bool Particles<NStructReal, NStructInt>::split_by_seperate_velocity(
    ParticleType& p1, ParticleType& p2, ParticleType& p3, ParticleType& p4) {
  // AllPrint() << "Old: p1 = " << p1 << std::endl;
  // AllPrint() << "Old: p2 = " << p2 << std::endl;

  Real mt = p1.rdata(iqp_) + p2.rdata(iqp_);
  Real wavg = mt / 4.0;

  // Calculate the average velocity and total energy.
  Real et = 0, uavg2 = 0;
  Real u[nDim3], du1[nDim3], du2[nDim3];
  for (int i = 0; i < nDim3; ++i) {
    u[i] = (p1.rdata(iqp_) * p1.rdata(iup_ + i) +
            p2.rdata(iqp_) * p2.rdata(iup_ + i)) /
           mt;

    et += 0.5 * p1.rdata(iqp_) * pow(p1.rdata(iup_ + i), 2);
    et += 0.5 * p2.rdata(iqp_) * pow(p2.rdata(iup_ + i), 2);

    uavg2 += pow(u[i], 2);

    // Get the direction of du1.
    du1[i] = p1.rdata(iup_ + i) - p2.rdata(iup_ + i);
  }

  Real du1Amp = l2_norm(du1, nDim3);
  if (du1Amp < 1e-16) {
    // If p1 and p2 have essentially the same velocity, do not split them. But
    // why the velocity difference can be so small?
    // A: with billions of particles, it can happen. I have done some
    // investigation, and it does not look like a bug.
    return false;
  }

  // The amplitude of du1 and du2.
  Real duAmp2 = et / (2 * wavg) - uavg2;
  if (duAmp2 < 0) {
    // Q: Why duAmp2 can be negative?
    // A: Rounding error.
    return false;
  }
  Real duAmp = sqrt(duAmp2);

  // Scale the amplitude of du1
  Real scale = duAmp / du1Amp;
  for (int i = 0; i < nDim3; ++i) {
    du1[i] *= scale;
  }

  {
    // Get the direction of du2
    Real utmp[nDim3];
    random_vector(randNum(), randNum(), utmp);

    // Correct the amplitide of du2
    for (int i = 0; i < nDim3; ++i) {
      du2[i] = utmp[i] * duAmp;
    }
  }

  // auto p_energy = [](const ParticleType& p) {
  //   Real energy = 0;
  //   for (int i = 0; i < nDim3; ++i) {
  //     energy += 0.5 * p.rdata(iqp_) * pow(p.rdata(iup_ + i), 2);
  //   }
  //   return energy;
  // };

  // Real eold = p_energy(p1) + p_energy(p2);
  // Real mold[3];
  // for (int i = 0; i < nDim3; ++i) {
  //   mold[i] = p1.rdata(iqp_) * p1.rdata(iup_ + i) +
  //             p2.rdata(iqp_) * p2.rdata(iup_ + i);
  // }

  set_ids(p3);
  set_ids(p4);

  p1.rdata(iqp_) = wavg;
  p2.rdata(iqp_) = wavg;
  p3.rdata(iqp_) = wavg;
  p4.rdata(iqp_) = wavg;

  for (int i = 0; i < nDim3; ++i) {
    p1.rdata(iup_ + i) = u[i] + du1[i];
    p2.rdata(iup_ + i) = u[i] - du1[i];

    p3.rdata(iup_ + i) = u[i] + du2[i];
    p4.rdata(iup_ + i) = u[i] - du2[i];
  }

  for (int i = 0; i < nDim; ++i) {
    p3.pos(i) = p1.pos(i);
    p4.pos(i) = p2.pos(i);
  }

  // Real enew = p_energy(p1) + p_energy(p2) + p_energy(p3) + p_energy(p4);
  // Real mnew[3];
  // for (int i = 0; i < nDim3; ++i) {
  //   mnew[i] = p1.rdata(iqp_) * p1.rdata(iup_ + i) +
  //             p2.rdata(iqp_) * p2.rdata(iup_ + i) +
  //             p3.rdata(iqp_) * p3.rdata(iup_ + i) +
  //             p4.rdata(iqp_) * p4.rdata(iup_ + i);
  // }

  // AllPrint() << "eold = " << eold << " enew = " << enew
  //            << " eold - enew = " << eold - enew << " mold - mnew "
  //            << mold[0] - mnew[0] << " " << mold[1] - mnew[1] << " "
  //            << mold[2] - mnew[2] << std::endl;

  // AllPrint() << "New: p1 = " << p1 << std::endl;
  // AllPrint() << "New: p2 = " << p2 << std::endl;
  // AllPrint() << "New: p3 = " << p3 << std::endl;
  // AllPrint() << "New: p4 = " << p4 << std::endl;

  return true;
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::split(Real limit,
                                               bool seperateVelocity) {
  timing_func("Pts::split");

  const int nInitial = product(nPartPerCell);

  IntVect iv = { AMREX_D_DECL(1, 1, 1) };
  if (!(do_tiling && tile_size == iv))
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {

    const Real dl = 0.1 * Geom(iLev).CellSize()[ix_] / nPartPerCell.max();

    const int nLowerLimit = nInitial * limit * pow(pLevRatio, iLev);

    const int nGoal = nLowerLimit > nInitial ? nLowerLimit : nInitial;

    const Real vol = dx[iLev].product();

    const Real vacuumMass = vacuum * vol;

    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {

      Vector<ParticleType> newparticles;

      auto& pTile = get_particle_tile(iLev, pti);
      AoS& particles = pTile.GetArrayOfStructs();

      const int nPartOrig = particles.size();

      if (nPartOrig > nLowerLimit)
        continue;

      const int nSplit =
          nGoal - nPartOrig > nPartOrig ? nPartOrig : nGoal - nPartOrig;

      Real totalMass = 0;
      for (auto& p : particles) {
        // So far, the vacuum limit is designed for OH-PT neutrals only. It is
        // not clear how should it be done for PIC, where electrons and ions
        // have different mass. --Yuxi
        totalMass += qomSign * p.rdata(iqp_);
      }
      if (totalMass < vacuumMass)
        continue;

      // Find the 'heaviest' nNew particles by sorting the weight
      // (charge).-----

      // Sort the particles by the location first to make sure the results
      // are the same for different number of processors
      std::sort(particles.begin(), particles.end(), compare_two_parts);

      const Real invLx = 1. / (phi[iLev][ix_] - plo[iLev][ix_]);
      const Real plox = plo[iLev][ix_];

      // Sort the particles by the weight in decending order.
      std::sort(
          particles.begin(), particles.end(),
          [&plox, &invLx](const ParticleType& pl, const ParticleType& pr) {
            const Real ql = fabs(pl.rdata(iqp_));
            const Real qr = fabs(pr.rdata(iqp_));
            if (fabs(ql - qr) > 1e-9 * (ql + qr)) {
              return ql > qr;
            }

            if (fabs(pl.pos(ix_) - pr.pos(ix_)) >
                1e-9 * (fabs(pl.pos(ix_)) + fabs(pr.pos(ix_)))) {
              return pl.pos(ix_) > pr.pos(ix_);
            }
            return false;
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

      if (is_neutral() || seperateVelocity) {
        Box bx = pti.tilebox();
        set_random_seed(iLev, bx.smallEnd(), IntVect(888));
      }

      if (seperateVelocity) {
        Vector<ParticleType*> pold;
        for (int ip = 0; ip < nSplit; ip++) {
          pold.push_back(&(particles[ip]));
        }
        split_particles_by_velocity(pold, newparticles);
      } else {
        for (int ip = 0; ip < nSplit; ip++) {
          auto& p = particles[ip];
          Real qp1 = p.rdata(iqp_);
          Real xp1 = p.pos(ix_);
          Real yp1 = p.pos(iy_);
          Real zp1 = nDim > 2 ? p.pos(iz_) : 0;
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

          int nNew = is_neutral() ? 7 : 1;

          p.rdata(iqp_) = qp1 / (nNew + 1.0);

          for (int iNew = 0; iNew < nNew; iNew++) {

            if (is_neutral()) {
              xp2 = xp1 + (xMax - xMin) * (randNum() - 0.5);
              yp2 = yp1 + (yMax - yMin) * (randNum() - 0.5);
              zp2 = zp1 + (zMax - zMin) * (randNum() - 0.5);
            } else {
              xp1 -= dpx;
              yp1 -= dpy;
              zp1 -= dpz;

              xp1 = bound(xp1, xMin, xMax);
              yp1 = bound(yp1, yMin, yMax);
              zp1 = bound(zp1, zMin, zMax);
              p.pos(ix_) = xp1;
              p.pos(iy_) = yp1;

              if (nDim > 2)
                p.pos(iz_) = zp1;
            }

            xp2 = bound(xp2, xMin, xMax);
            yp2 = bound(yp2, yMin, yMax);
            zp2 = bound(zp2, zMin, zMax);

            ParticleType pnew;
            set_ids(pnew);

            pnew.pos(ix_) = xp2;
            pnew.pos(iy_) = yp2;
            if (nDim > 2)
              pnew.pos(iz_) = zp2;
            pnew.rdata(iup_) = up1;
            pnew.rdata(ivp_) = vp1;
            pnew.rdata(iwp_) = wp1;
            pnew.rdata(iqp_) = qp1 / (nNew + 1.0);
            newparticles.push_back(pnew);
          }
        }
      }

      for (auto& p : newparticles) {
        pTile.push_back(p);
      }
    }
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
bool Particles<NStructReal, NStructInt>::merge_particles_accurate(
    int iLev, AoS& particles, Vector<int>& partIdx, Vector<int>& idx_I,
    int nPartCombine, int nPartNew, Vector<Real>& x, Real velNorm) {
  timing_func("Pts::merge_particles_accurate");

  constexpr int nVar = 5;
  constexpr int iq_ = 0, iu_ = 1, iv_ = 2, iw_ = 3, ie_ = 4;
  const Real coefVel = 1, coefPos = 1;

  Vector<Real> ref(nVar, 0);
  Array2D<Real, 0, nVar - 1, 0, nVar> a;
  for (int i = 0; i < nVar; ++i)
    for (int j = 0; j < nVar + 1; ++j) {
      a(i, j) = 0;
    }

  // Find the center of the particles, and sort the particles based
  // on its distance to the 6-D center.
  //----------------------------------------------------------
  Vector<Real> middle(nDim + nDim3, 0);
  for (int pID : partIdx) {
    for (int iDir = ix_; iDir <= iz_; iDir++) {
      if (iDir < nDim)
        middle[iDir] += particles[pID].pos(iDir);
      middle[nDim + iDir] += particles[pID].rdata(iDir);
    }
  }

  for (int i = 0; i < middle.size(); ++i) {
    middle[i] /= partIdx.size();
  }

  auto calc_distance2_to_center = [&, this](int pID) {
    Real dl2 = 0, dvel2 = 0;
    for (int iDir = ix_; iDir <= iz_; iDir++) {

      if (iDir < nDim) {
        Real pos = particles[pID].pos(iDir);
        dl2 += pow((pos - middle[iDir]) * invDx[iLev][iDir], 2);
      }

      Real v = particles[pID].rdata(iDir);
      dvel2 += pow((v - middle[nDim + iDir]) * velNorm, 2);
    }
    return coefPos * dl2 + coefVel * dvel2;
  };

  std::sort(partIdx.begin(), partIdx.end(),
            [this, &particles, calc_distance2_to_center](const int& idl,
                                                         const int& idr) {
              Real dll = calc_distance2_to_center(idl);
              Real dlr = calc_distance2_to_center(idr);

              if (fabs(dll - dlr) > 1e-9 * (dll + dlr)) {
                return dll < dlr;
              }
              return false;
            });

  /*
      Delete 1 particle out of 6 particles:
      1) Choose two particles that are closest to each other.
      2) Delete the lighter one.
      3) Distribute its weights to another 5 particles to conserve
         mass, momentum and energy.
   */

  idx_I.resize(nPartCombine, 0);
  for (int ip = 0; ip < nPartCombine; ip++) {
    idx_I[ip] = partIdx[ip];
  }

  // Calculate the center of the particles for combination
  for (int i = 0; i < middle.size(); ++i) {
    middle[i] = 0;
  }
  for (int pID : idx_I) {
    for (int iDir = ix_; iDir <= iz_; iDir++) {
      if (iDir < nDim)
        middle[iDir] += particles[pID].pos(iDir);

      middle[nDim + iDir] += particles[pID].rdata(iDir);
    }
  }
  for (int i = 0; i < middle.size(); ++i) {
    middle[i] /= nPartCombine;
  }

  bool doCombine = true;
  for (int pID : idx_I) {
    Real distance = sqrt(calc_distance2_to_center(pID));
    if (distance > mergeThresholdDistance) {
      // printf("Warning: distance=%e\n", distance);
      doCombine = false;
    }
  }

  if (!doCombine)
    return false;

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

        Real dx = invDx[iLev][iDir] * (particles[idx_I[ip1]].pos(iDir) -
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
  // Q: Why is it 'l>(1+1e-9)*r' instead of 'l>r'?
  // A: The particle weights can be the same for some cases. 'l>r'
  // may return random results due to the truncation error.
  if (fabs(particles[idx_I[pair1]].rdata(iqp_)) >
      (1 + 1e-9) * fabs(particles[idx_I[pair2]].rdata(iqp_))) {
    iPartDel = pair2;
  }

  std::swap(idx_I[iPartDel], idx_I[nPartCombine - 1]);

  //-----------Solve the new particle weights-------
  for (int ip = 0; ip < nPartCombine; ip++) {
    const Real qp = particles[idx_I[ip]].rdata(iqp_);
    const Real up = particles[idx_I[ip]].rdata(iup_);
    const Real vp = particles[idx_I[ip]].rdata(ivp_);
    const Real wp = particles[idx_I[ip]].rdata(iwp_);
    const Real v2 = (pow(up, 2) + pow(vp, 2) + pow(wp, 2));

    if (ip < nVar) {
      a(iq_, ip) = 1;
      a(iu_, ip) = up;
      a(iv_, ip) = vp;
      a(iw_, ip) = wp;
      a(ie_, ip) = v2;
    }

    a(iq_, nVar) += qp;
    a(iu_, nVar) += qp * up;
    a(iv_, nVar) += qp * vp;
    a(iw_, nVar) += qp * wp;
    a(ie_, nVar) += qp * v2;
  }

  const Real csmall = 1e-9;
  const Real tmp = csmall * fabs(1. / a(iq_, nVar));
  for (int i = iq_; i <= ie_; ++i) {
    ref[i] = fabs(a(i, nVar) * tmp);
  }

  x.resize(nVar, 0);
  bool isSolved = linear_solver_Gauss_Elimination<Real, nVar, nVar + 1>(
      nVar, nVar + 1, a, x, ref);

  if (isSolved) {
    // All the particle weights should have the same sign.
    Real qt = x[0];
    for (int ip = 0; ip < nPartNew; ip++) {
      if (qt * x[ip] <= 0) {
        isSolved = false;
        break;
      }
    }
  }

  return isSolved;
}

//==========================================================
template <int NStructReal, int NStructInt>
bool Particles<NStructReal, NStructInt>::merge_particles_fast(
    int iLev, AoS& particles, Vector<int>& partIdx, Vector<int>& idx_I,
    int nPartCombine, int nPartNew, Vector<Real>& x, long seed) {
  timing_func("Pts::merge_particles_fast");

  constexpr int iq_ = 0, iu_ = 1, iv_ = 2, iw_ = 3, ie_ = 4;
  constexpr int nPartNewMax = 16;
  constexpr int nVarMax = nPartNewMax + 5;

  int nVar = nPartNew + 5;

  Vector<Real> ref(nVar, 0);
  Array2D<Real, 0, nVarMax - 1, 0, nVarMax> a;
  for (int i = 0; i < nVar; ++i)
    for (int j = 0; j < nVar + 1; ++j) {
      a(i, j) = 0;
    }

  const Real invLx = 1. / (phi[iLev][ix_] - plo[iLev][ix_]);
  const Real plox = plo[iLev][ix_];
  // Q: Sort the particles by weights in ascending order.
  // But, why is it required here?
  // A: Eliminate randomness.
  std::sort(partIdx.begin(), partIdx.end(),
            [&particles, &invLx, &plox](int idLeft, int idRight) {
              const Real ql = fabs(particles[idLeft].rdata(iqp_));
              const Real qr = fabs(particles[idRight].rdata(iqp_));
              if (fabs(ql - qr) > 1e-9 * (ql + qr)) {
                return ql < qr;
              }

              Real xl = particles[idLeft].pos(ix_);
              Real xr = particles[idRight].pos(ix_);
              if (fabs(xl - xr) > 1e-9 * (fabs(xl) + fabs(xr))) {
                return xl < xr;
              }
              return false;
            });

  if (mergeLight) {
    idx_I.resize(nPartCombine, 0);
    for (int ip = 0; ip < nPartCombine; ip++) {
      idx_I[ip] = partIdx[ip];
    }

    Real plight = 1e99, pheavy = 0;
    for (int ip = 0; ip < nPartCombine; ip++) {
      auto& p = particles[idx_I[ip]];
      Real w = fabs(p.rdata(iqp_));
      if (w < plight)
        plight = w;
      if (w > pheavy)
        pheavy = w;
    }

    if (pheavy / plight > mergePartRatioMax)
      return false;

    randNum.set_seed(seed);
    shuffle_fish_yates(idx_I, randNum);

  } else {
    randNum.set_seed(seed);
    shuffle_fish_yates(partIdx, randNum);

    idx_I.resize(nPartCombine, 0);
    for (int ip = 0; ip < nPartCombine; ip++) {
      idx_I[ip] = partIdx[ip];
    }
  }

  // Sum the moments of all the old particles.
  for (int ip = 0; ip < nPartCombine; ip++) {
    const Real qp = particles[idx_I[ip]].rdata(iqp_);
    const Real up = particles[idx_I[ip]].rdata(iup_);
    const Real vp = particles[idx_I[ip]].rdata(ivp_);
    const Real wp = particles[idx_I[ip]].rdata(iwp_);
    const Real v2 = 0.5 * (pow(up, 2) + pow(vp, 2) + pow(wp, 2));
    a(nPartNew + iq_, nVar) += qp;
    a(nPartNew + iu_, nVar) += qp * up;
    a(nPartNew + iv_, nVar) += qp * vp;
    a(nPartNew + iw_, nVar) += qp * wp;
    a(nPartNew + ie_, nVar) += qp * v2;
  }

  const Real invAvg = 2 * nPartNew / a(nPartNew + iq_, nVar);
  for (int ip = 0; ip < nPartNew; ip++) {
    const Real qp = particles[idx_I[ip]].rdata(iqp_);
    const Real up = particles[idx_I[ip]].rdata(iup_);
    const Real vp = particles[idx_I[ip]].rdata(ivp_);
    const Real wp = particles[idx_I[ip]].rdata(iwp_);
    const Real v2 = 0.5 * (pow(up, 2) + pow(vp, 2) + pow(wp, 2));

    a(ip, nVar) = 2;

    a(ip, ip) = 2. / qp;
    a(ip, nPartNew + iq_) = 1;
    a(ip, nPartNew + iu_) = up;
    a(ip, nPartNew + iv_) = vp;
    a(ip, nPartNew + iw_) = wp;
    a(ip, nPartNew + ie_) = v2;

    a(nPartNew + iq_, ip) = 1;
    a(nPartNew + iu_, ip) = up;
    a(nPartNew + iv_, ip) = vp;
    a(nPartNew + iw_, ip) = wp;
    a(nPartNew + ie_, ip) = v2;
  }

  const Real csmall = 1e-9;
  const Real tmp = csmall * invAvg;
  for (int i = 0; i < nVar; ++i) {
    if (i < nPartNew) {
      ref[i] = tmp;
    } else {
      ref[i] = fabs(a(i, nVar) * tmp);
    }
  }

  x.resize(nVar, 0);
  bool isSolved = linear_solver_Gauss_Elimination<Real, nVarMax, nVarMax + 1>(
      nVar, nVar + 1, a, x, ref);

  if (isSolved) {
    // All the particle weights should have the same sign.
    Real qt = x[0];
    for (int ip = 0; ip < nPartNew; ip++) {
      if (qt * x[ip] <= 0) {
        isSolved = false;
        break;
      }
    }
  }

  if (isSolved) {
    for (int ip = 0; ip < nPartNew; ip++) {
      Real pold = particles[idx_I[ip]].rdata(iqp_);
      Real pnew = x[ip];

      Real c0 = Real(nPartCombine) / nPartNew * mergeRatioMax;
      if (pnew / pold > c0 || pold / pnew > c0) {
        isSolved = false;
        break;
      }
    }
  }

  return isSolved;
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::merge(Real limit) {
  timing_func("Pts::merge");
  IntVect iv = { AMREX_D_DECL(1, 1, 1) };
  if (!(do_tiling && tile_size == iv))
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {

    int nPartGoal = product(nPartPerCell) * limit * pow(pLevRatio, iLev);

    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {

      // It is assumed the tile size is 1x1x1.
      Box bx = pti.tilebox();
      long seed = set_random_seed(iLev, bx.smallEnd(), IntVect(777));

      AoS& particles = pti.GetArrayOfStructs();

      const int nPartOrig = particles.size();

      if (nPartOrig <= nPartGoal)
        continue;

      // The range of the velocity domain:
      // [-r0,r0]*thermal_velocity+bulk_velocity
      const Real r0 = fastMerge ? 2.0 : 1.0;

      // Phase space cell number in one direction.
      // The const 0.5/0.8 is choosen by experiments.
      int nCell = 0;
      if (fastMerge) {
        nCell = r0 * ceil(0.5 * pow(nPartOrig, 1. / nDim));
      } else {
        nCell = r0 * ceil(0.8 * pow(nPartOrig, 1. / nDim));
      }

      if (nCell < 3)
        continue;

      // Sort the particles by the location first to make sure the results
      // are the same for different number of processors
      std::sort(particles.begin(), particles.end(), compare_two_parts);

      // One particle may belong to more than one velocity bins, but it can be
      // only merged at most once.
      std::vector<bool> merged;
      merged.resize(nPartOrig, false);

      //----------------------------------------------------------------
      // Estimate the bulk velocity and thermal velocity.
      Real uBulk[nDim3] = { 0, 0, 0 };
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];
        for (int iDir = 0; iDir < 3; iDir++) {
          uBulk[iDir] += pcl.rdata(iDir);
        }
      }

      for (int iDir = 0; iDir < nDim3; iDir++) {
        uBulk[iDir] /= nPartOrig;
      }

      Real thVel = 0, thVel2 = 0;
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];
        for (int iDir = 0; iDir < nDim3; iDir++) {
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

      Real dv = (2.0 * r0 * thVel) / nCell;
      Real invDv = (dv < 1e-13) ? 0 : 1.0 / dv;

      // Velocity domain range.
      Real velMin_D[nDim3], velMax_D[nDim3];
      for (int iDir = 0; iDir < nDim3; iDir++) {
        Real dvshift = (randNum() - 0.5) * dv;
        velMin_D[iDir] = -r0 * thVel + uBulk[iDir] + dvshift;
        velMax_D[iDir] = r0 * thVel + uBulk[iDir] + dvshift;
      }

      int iCell_D[nDim3];
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];

        bool isOutside = false;
        for (int iDim = 0; iDim < nDim3; iDim++) {
          if (pcl.rdata(iDim) < velMin_D[iDim] ||
              pcl.rdata(iDim) > velMax_D[iDim])
            isOutside = true;
        }
        if (isOutside)
          continue;

        for (int iDim = 0; iDim < nDim3; iDim++) {
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

              Vector<int> cellIdx = { xCell, yCell, zCell };

              Real binMin_D[nDim3], binMax_D[nDim3];

              for (int iDim = 0; iDim < nDim3; iDim++) {
                binMin_D[iDim] =
                    velMin_D[iDim] + (cellIdx[iDim] - velBinBufferSize) * dv;

                binMax_D[iDim] = velMin_D[iDim] +
                                 (cellIdx[iDim] + 1 + velBinBufferSize) * dv;
              }

              bool isInside = true;
              for (int iDim = 0; iDim < nDim3; iDim++) {
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
            for (int i = 0; i < phasePartIdx_III[iu][iv][iw].size(); ++i) {
              auto& pIdx = phasePartIdx_III[iu][iv][iw];
              int pid = pIdx[i];
              if (!merged[pid]) {
                partIdx.push_back(pid);
              }
            }

            if (partIdx.size() < nPartNew + 1)
              continue;

            Vector<Real> x;
            Vector<int> idx_I;

            int nOld = nPartCombine;
            bool isSolved;
            if (fastMerge) {
              if (nOld > partIdx.size())
                nOld = partIdx.size();

              for (int iTry = 0; iTry < nMergeTry; iTry++) {
                long sd = seed + iu * 777 + iv * 77 + iw + iTry;
                isSolved = merge_particles_fast(iLev, particles, partIdx, idx_I,
                                                nOld, nPartNew, x, sd);
                if (isSolved)
                  break;
              }
            } else {
              isSolved = merge_particles_accurate(
                  iLev, particles, partIdx, idx_I, nOld, nPartNew, x, velNorm);
            }
            if (!isSolved)
              continue;

            //----------------------------------------------

            Real plight = 1e99, pheavy = 0;
            for (int ip = 0; ip < nOld; ip++) {
              auto& p = particles[idx_I[ip]];
              Real w = fabs(p.rdata(iqp_));
              if (w < plight)
                plight = w;
              if (w > pheavy)
                pheavy = w;
            }
            // AllPrint() << "pheavy/plight = " << pheavy / plight <<
            // std::endl;

            // Adjust weight.
            for (int ip = 0; ip < nPartNew; ip++) {
              auto& p = particles[idx_I[ip]];
              p.rdata(iqp_) = x[ip];
              merged[idx_I[ip]] = true;
            }
            // Mark for deletion
            for (int ip = nPartNew; ip < nOld; ip++) {
              particles[idx_I[ip]].id() = -1;
              particles[idx_I[ip]].rdata(iqp_) = 0;
              merged[idx_I[ip]] = true;
            }
          }
    }
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
bool Particles<NStructReal, NStructInt>::do_inject_particles_for_this_cell(
    const Box& bx, const Array4<const int>& status, const IntVect ijk,
    IntVect& ijksrc) {

  // This cell should be a boundary cell at least.
  if (!bit::is_lev_boundary(status(ijk)))
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

          IntVect ijk1 = ijk + IntVect{ AMREX_D_DECL(di, dj, dk) };
          if (!bit::is_lev_boundary(status(ijk1))) {
            // The first neighbor cell that is NOT a boundary cell.
            if (bx.contains(ijk1)) {
              ijksrc = ijk1;
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

IOParticles::IOParticles(Particles& other, Grid* gridIn, Real no2outL,
                         Real no2outV, Real no2outM, RealBox IORange)
    : Particles(gridIn, nullptr, nullptr, other.get_speciesID(),
                other.get_charge(), other.get_mass(),
                IntVect(AMREX_D_DECL(-1, -1, -1)), other.part_mode()) {

  no2outM *= qomSign * get_mass();

  const bool doLimit = IORange.ok();

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const auto& plevelOther = other.GetParticles(iLev);
    auto& plevel = GetParticles(iLev);
    for (MFIter mfi = other.MakeMFIter(iLev); mfi.isValid(); ++mfi) {
      auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());

      if (plevelOther.find(index) == plevelOther.end())
        continue;

      const auto& tileOther = plevelOther.at(index);

      if (tileOther.numParticles() == 0)
        continue;

      const AoS& aosOther = tileOther.GetArrayOfStructs();

      const Box& bx = other.cell_status(iLev)[mfi].box();
      const Array4<int const>& status = other.cell_status(iLev)[mfi].array();

      const IntVect lowCorner = bx.smallEnd();
      const IntVect highCorner = bx.bigEnd();

      for (auto p : aosOther) {
        if (other.is_outside_active_region(p, status, lowCorner, highCorner)) {
          // redistribute_particles() may fail if the ghost cell particles'
          // IDs are not -1 (marked for deletion);
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
  }
  redistribute_particles();
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::sample_charge_exchange(
    Real* vp, Real* vh, Real* up, Real vth, CrossSection cs) {

  timing_func("Pts::sample_charge_exchange");

  // M (holds normalization costants for both distributions), g(vp) is the
  // maxwellian distribution
  Real sup[3] = { up[0] - 3. * vth, up[1] - 3. * vth, up[2] - 3. * vth };
  Real M = charge_exchange_dis(sup, vh, up, vth, cs) /
           exp(-((sup[0] - up[0]) * (sup[0] - up[0]) +
                 (sup[1] - up[1]) * (sup[1] - up[1]) +
                 (sup[2] - up[2]) * (sup[2] - up[2])) /
               (vth * vth));

  bool accepted = false;
  while (!accepted) {

    {
      Real prob, theta, Uth;

      // u = X velocity
      prob = sqrt(-2.0 * log(1.0 - .999999999 * randNum()));
      theta = 2.0 * M_PI * randNum();
      Uth = vth / sqrt(2.0);
      vp[0] = Uth * prob * cos(theta) + up[0];
      // v = Y velocity
      vp[1] = Uth * prob * sin(theta) + up[1];

      // w = Z velocity
      prob = sqrt(-2.0 * log(1.0 - .999999999 * randNum()));
      theta = 2.0 * M_PI * randNum();
      vp[2] = Uth * prob * cos(theta) + up[2];
    }

    if (randNum() < charge_exchange_dis(vp, vh, up, vth, cs) /
                        (M * exp(-((vp[0] - up[0]) * (vp[0] - up[0]) +
                                   (vp[1] - up[1]) * (vp[1] - up[1]) +
                                   (vp[2] - up[2]) * (vp[2] - up[2])) /
                                 (vth * vth)))) {
      accepted = true;
    }
  }
}

template <int NStructReal, int NStructInt>
Real Particles<NStructReal, NStructInt>::charge_exchange_dis(Real* vp, Real* vh,
                                                             Real* up, Real vth,
                                                             CrossSection cs) {
  Real dv_D[3], dv2 = 0, dv = 0;
  for (int i = 0; i < 3; ++i) {
    dv_D[i] = vh[i] - vp[i];
    dv2 += dv_D[i] * dv_D[i];
  }

  if (dv2 == 0)
    return 0.0;

  dv = sqrt(dv2);

  Real erel = 0.5 * 1.674E-27 * dv2 * 6.2415E15; // in keV

  Real sigma = 0;
  if (cs == CrossSection::LS) {
    sigma = (4.15 - 0.531 * log(erel)) * (4.15 - 0.531 * log(erel)) *
            pow(1 - exp(-67.3 / erel), 4.5) * 1E-20; // cross section in m^2
  } else if (cs == CrossSection::MT) {
    Real dvcm = dv * 1E2;                             // velocity in cm/s
    sigma = pow(1.6 - 0.0695 * log(dvcm), 2) * 1e-18; // cross section in m^2
  }

  Real dvpup2 = 0;
  for (int i = 0; i < 3; ++i) {
    dvpup2 += pow(vp[i] - up[i], 2);
  }

  return dv * sigma * exp(-dvpup2 / (vth * vth));
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::get_analytic_ion_fluid(
    const RealVect xyz, Real& rhoIon, Real& cs2Ion,
    amrex::Real (&uIon)[nDim3]) {

  // The units of ionOH are assumed to be:
  // r: AU
  // rho: amu/cc
  // T: K
  // U: km/s

  Real r = max(xyz.vectorLength(), 1e-9 * Geom(0).CellSize(0));

  Real rAU = r * fi->get_No2SiL() / cAUSI;

  if (rAU > ionOH.rAnalytic) {
    Abort("Error: rAU > ionOH.rAnalytic");
  }

  if (ionOH.doGetFromOH) {
    Real xSI = xyz[ix_] * fi->get_No2SiL();
    Real ySI = xyz[iy_] * fi->get_No2SiL();
    Real zSI = xyz[iz_] * fi->get_No2SiL();

    Real temp, ur, b[nDim3];
    OH_get_solar_wind(&xSI, &ySI, &zSI, &rhoIon, &ur, &temp, b);

    // v_th = sqrt(2kT/m); m/s
    cs2Ion = 2 * cBoltzmannSI * temp / cProtonMassSI; // m^2/s^2

    for (int i = 0; i < nDim; ++i) {
      uIon[i] = ur * xyz[i] / r;
    }
  } else {
    Real r0 = 0;
    if (rAU < ionOH.rCutoff) {
      r0 = ionOH.rAnalytic / ionOH.rCutoff;
    } else {
      r0 = ionOH.rAnalytic / rAU;
    }

    rhoIon = ionOH.swRho * pow(r0, 2) * 1e6; // amu/cc -> amu/m^3

    // v_th = sqrt(2kT/m); m/s
    cs2Ion = 2 * cBoltzmannSI * ionOH.swT / cProtonMassSI; // m^2/s^2

    for (int i = 0; i < nDim; ++i) {
      uIon[i] = ionOH.swU * xyz[i] / r * 1e3; // km/s -> m/s
    }
  }

  // AllPrint() << "r = " << r << " xyz = " << xyz[0] << ", " << xyz[1] << ", "
  //            << xyz[2] << " rhoIon = " << rhoIon << ", cs2Ion = " << cs2Ion
  //            << ", uIon = " << uIon[0] << ", " << uIon[1] << ", " << uIon[2]
  //            << std::endl;
}

// Get the iFluid-th ion fluid properties at the location xyz
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::get_ion_fluid(
    FluidInterface* stateOH, PIter& pti, const int iLev, const int iFluid,
    const RealVect xyz, Real& rhoIon, Real& cs2Ion, Real (&uIon)[nDim3]) {

  Real rAU = xyz.vectorLength() * stateOH->get_No2SiL() / cAUSI;

  if (iFluid == 0 && rAU < ionOH.rAnalytic) {
    get_analytic_ion_fluid(xyz, rhoIon, cs2Ion, uIon);
    return;
  }

  // amu/m^3
  rhoIon = stateOH->get_fluid_mass_density(pti, xyz, iFluid, iLev) *
           stateOH->get_No2SiRho() / cProtonMassSI;

  // cs = sqrt(P/n); m/s
  // Assume p = pi + pe = 2pi, so divide by sqrt(2.0).
  double cs = stateOH->get_fluid_uth(pti, xyz, iFluid, iLev) *
              stateOH->get_No2SiV() / sqrt(2.0);

  // cs2Ion = 2*P/n. The definition of thermal speed in get_uth_iso() is
  // different from the requirement in OH_get_charge_exchange_wrapper().
  // See page 92 of Adam Michael's thesis.
  cs2Ion = 2 * pow(cs, 2);

  uIon[ix_] =
      stateOH->get_fluid_ux(pti, xyz, iFluid, iLev) * stateOH->get_No2SiV();
  uIon[iy_] =
      stateOH->get_fluid_uy(pti, xyz, iFluid, iLev) * stateOH->get_No2SiV();
  uIon[iz_] =
      stateOH->get_fluid_uz(pti, xyz, iFluid, iLev) * stateOH->get_No2SiV();
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::charge_exchange(
    Real dt, FluidInterface* stateOH, FluidInterface* sourcePT2OH,
    SourceInterface* source, bool kineticSource,
    Vector<std::unique_ptr<PicParticles> >& sourceParts, bool doSelectRegion,
    int nppc) {
  std::string nameFunc = "Pts::charge_exchange";

  timing_func(nameFunc);

  if (dt <= 0)
    return;

  // for (auto& ptr : sourceParts) {
  //   ptr->clearParticles();
  // }

  struct NeuPlasmaPair {
    Real q; // weight
    Real vp[3];
    Real vh[3];
    Real up[3];
    Real xyz[3];
    Real vth;
  };

  Real maxExchangeRatio = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      AoS& particles = pti.GetArrayOfStructs();

      if (kineticSource) {
        // Sort the particles by the location first to make sure the results
        // are the same for different number of processors
        std::sort(particles.begin(), particles.end(), compare_two_parts);
      }

      // It is assumed the tile size is 1x1x1.
      Box bx = pti.tilebox();
      auto cellIdx = bx.smallEnd();

      int iRegion = 0;
      if (kineticSource && doSelectRegion) {
        const int iFluid = 0;
        iRegion = stateOH->get_neu_source_region(pti, cellIdx, iFluid, iLev);
        if (iRegion < 0)
          continue;
      }
      // ParticleTileType
      auto& spTile =
          sourceParts[iRegion]->get_particle_tile(iLev, pti, cellIdx);

      Vector<NeuPlasmaPair> neuPlasmaPairs;

      if (kineticSource) {
        set_random_seed(iLev, cellIdx, IntVect(999));
      }

      for (auto& p : particles) {
        if (p.id() < 0)
          continue;

        RealVect xyz;
        for (int i = 0; i < nDim; ++i) {
          xyz[i] = p.pos(i);
        }

        double cs2Neu = 0, uNeu[3], rhoNeu;
        double cs2Ion, uIon[3], rhoIon;
        double ion2neu[5], neu2ion[5];
        const int iRho_ = 0, iUx_ = 1, iUy_ = 2, iUz_ = 3, iP_ = 4;
        const int iRhoUx_ = iUx_, iRhoUy_ = iUy_, iRhoUz_ = iUz_, iE_ = iP_;

        // amu/m^3
        rhoNeu = qomSign * p.rdata(iqp_) * get_mass() * invVol[iLev] *
                 stateOH->get_No2SiRho() / cProtonMassSI;

        for (int i = 0; i < nDim; ++i) {
          uNeu[i] = p.rdata(iup_ + i) * stateOH->get_No2SiV();
        }

        // A neutral particle interacts with all the ion fluids.
        for (int fluidID = 0; fluidID < stateOH->get_nFluid(); fluidID++) {
          get_ion_fluid(stateOH, pti, iLev, fluidID, xyz, rhoIon, cs2Ion, uIon);

          OH_get_charge_exchange_wrapper(&rhoIon, &cs2Ion, uIon, &rhoNeu,
                                         &cs2Neu, uNeu, ion2neu, neu2ion);

          // The function above returns number density changing rate.
          ion2neu[iRho_] *= cProtonMassSI;
          neu2ion[iRho_] *= cProtonMassSI;

          Real dtSI = dt * stateOH->get_No2SiT();
          // Print() << "rhoion = " << rhoIon << " cs2Ion = " << cs2Ion
          //         << " rhoNeu = " << rhoNeu << " cs2Neu = " << cs2Neu
          //         << " dtSI = " << dtSI << std::endl;
          for (int i = iRho_; i <= iP_; ++i) {
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
            for (int i = iRho_; i <= iP_; ++i) {
              ion2neu[i] *= ratio;
              neu2ion[i] *= ratio;
            }
          } else {
            // Reduce particle mass due to charge exchange

            Real ratio = massExchange / p.rdata(iqp_);
            if (ratio > maxExchangeRatio)
              maxExchangeRatio = ratio;

            p.rdata(iqp_) = p.rdata(iqp_) - massExchange;
          }

          {
            int iFluidAddTo;
            Real rhoIonAddTo, cs2IonAddTo, uIonAddTo[3];

            const int iSW = 2; // Pop3
            switch (stateOH->get_nFluid()) {
              case 1:
                iFluidAddTo = fluidID;
                rhoIonAddTo = rhoIon;
                break;
              case 2:
                if (iRegion == iSW) {
                  iFluidAddTo = 1; // Pu3
                } else {
                  iFluidAddTo = 0; // SW
                }
                get_ion_fluid(stateOH, pti, iLev, iFluidAddTo, xyz, rhoIonAddTo,
                              cs2IonAddTo, uIonAddTo);
                break;
              default:
                Abort("Error: nFluid > 2 is not supported yet.");
            }

            // Q: Why is (neu2ion-ion2neu) divided by rhoIon?
            // A: What passed between PT and OH is 'source per ion density'
            // instead of source. The ion density will be multiplied back in OH
            // ModUser.f90
            sourcePT2OH->add_rho_to_loc((neu2ion[iRho_] - ion2neu[iRho_]) /
                                            rhoIonAddTo,
                                        pti, xyz, iFluidAddTo, iLev);
            sourcePT2OH->add_mx_to_loc((neu2ion[iRhoUx_] - ion2neu[iRhoUx_]) /
                                           rhoIonAddTo,
                                       pti, xyz, iFluidAddTo, iLev);
            sourcePT2OH->add_my_to_loc((neu2ion[iRhoUy_] - ion2neu[iRhoUy_]) /
                                           rhoIonAddTo,
                                       pti, xyz, iFluidAddTo, iLev);
            sourcePT2OH->add_mz_to_loc((neu2ion[iRhoUz_] - ion2neu[iRhoUz_]) /
                                           rhoIonAddTo,
                                       pti, xyz, iFluidAddTo, iLev);
            sourcePT2OH->add_p_to_loc((neu2ion[iP_] - ion2neu[iP_]) /
                                          rhoIonAddTo,
                                      pti, xyz, iFluidAddTo, iLev);
          }

          if (ion2neu[iRho_] > 0) { // Add source to nodes.

            if (kineticSource) {

              NeuPlasmaPair pair;
              pair.q = massExchange;
              for (int i = 0; i < nDim3; ++i) {
                pair.vh[i] = uNeu[i];
                pair.up[i] = uIon[i];
                pair.vth = sqrt(cs2Ion);
              }

              for (int i = 0; i < nDim; ++i) {
                pair.xyz[i] = p.pos(ix_ + i);
              }

              neuPlasmaPairs.push_back(pair);

            } else {

              Real si2no_v[5];
              si2no_v[iRho_] = source->get_Si2NoRho();
              si2no_v[iRhoUx_] = source->get_Si2NoV() * si2no_v[iRho_];
              si2no_v[iRhoUy_] = si2no_v[iRhoUx_];
              si2no_v[iRhoUz_] = si2no_v[iRhoUx_];
              si2no_v[iP_] = source->get_Si2NoP();

              Real m2 = 0;
              for (int i = iRhoUx_; i <= iRhoUz_; ++i) {
                m2 += pow(ion2neu[i], 2);
              }

              const Real gamma = 5. / 3;
              // P = (gamma-1)*(E - 0.5*rho*u2)
              ion2neu[iP_] =
                  (gamma - 1) * (ion2neu[iE_] - 0.5 * m2 / ion2neu[iRho_]);

              if (ion2neu[iP_] < 0) {
                ion2neu[iP_] = 0;
              }

              // source saves changing rate (density/s...).
              source->add_rho_to_loc(ion2neu[iRho_] * si2no_v[iRho_] / dt, pti,
                                     xyz, speciesID, iLev);
              source->add_mx_to_loc(ion2neu[iRhoUx_] * si2no_v[iRhoUx_] / dt,
                                    pti, xyz, speciesID, iLev);
              source->add_my_to_loc(ion2neu[iRhoUy_] * si2no_v[iRhoUy_] / dt,
                                    pti, xyz, speciesID, iLev);
              source->add_mz_to_loc(ion2neu[iRhoUz_] * si2no_v[iRhoUz_] / dt,
                                    pti, xyz, speciesID, iLev);
              source->add_p_to_loc(ion2neu[iP_] * si2no_v[iP_] / dt, pti, xyz,
                                   speciesID, iLev);
            }
          }
          // p.id() = -1;
        }
      } // for p

      if (kineticSource) {
        // Sample the velocity distribution function.

        Vector<NeuPlasmaPair> newPairs;

        Real wt = 0;

        Vector<Real> weights;
        weights.resize(neuPlasmaPairs.size());
        for (int i = 0; i < neuPlasmaPairs.size(); ++i) {
          weights[i] = neuPlasmaPairs[i].q;
          wt += neuPlasmaPairs[i].q;
        }

        std::vector<int> idx = random_select_weighted_n(weights, nppc, randNum);

        Real wtnew = 0;
        for (int i : idx) {
          newPairs.push_back(neuPlasmaPairs[i]);
          wtnew += neuPlasmaPairs[i].q;
        }
        Real scale = wt / wtnew;

        for (auto& pair : newPairs) {
          sample_charge_exchange(pair.vp, pair.vh, pair.up, pair.vth,
                                 CrossSection::MT);
          for (int i = 0; i < nDim3; ++i) {
            pair.vp[i] *= stateOH->get_Si2NoV();
          }

          PicParticle newp;
          newp.rdata(iqp_) = pair.q * scale;

          for (int i = 0; i < nDim3; ++i) {
            newp.rdata(iup_ + i) = pair.vp[i];
          }

          for (int i = 0; i < nDim; ++i) {
            newp.pos(ix_ + i) = pair.xyz[i];
          }

          spTile.push_back(newp);
        }
      }
    } // for pti
  }

  ParallelDescriptor::ReduceRealMax(maxExchangeRatio);

  Print() << "maxExchangeRatio = " << maxExchangeRatio << std::endl;
  if (maxExchangeRatio > 0.2) {
    Print() << "Warning: maybe the charge exchange rate is too high within one "
               "time step! Reducing the time step will help to slow down the "
               "charge exchange."
            << std::endl;
  }
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_source_particles(
    std::unique_ptr<PicParticles>& sourcePart, IntVect ppc,
    const bool adaptivePPC) {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      // It is assumed the tile size is 1x1x1.
      Box bx = pti.tilebox();
      auto cellIdx = bx.smallEnd();

      ParticleTileType& pTile = get_particle_tile(iLev, pti, cellIdx);
      AoS& particles = pti.GetArrayOfStructs();

      // ParticleTileType
      auto& spTile = sourcePart->get_particle_tile(iLev, pti, cellIdx);

      // AoS type
      auto& sps = spTile.GetArrayOfStructs();

      if (sps.size() == 0)
        continue;

      Real rhoSource = 0;
      for (auto& p : sps) {
        rhoSource += p.rdata(iqp_);
      }

      if (adaptivePPC) {
        set_random_seed(iLev, cellIdx, IntVect(787));
        // Adjust ppc so that the weight of the
        // source particles is not too small.

        Real rho = 0;
        for (auto& p : particles) {
          rho += p.rdata(iqp_);
        }

        Real avgInitW = rho / product(nPartPerCell);
        Real avgSourceW = rhoSource / product(ppc);

        Real targetSourceW = avgInitW * 0.1;

        if (avgSourceW < targetSourceW) {
          Real ratio = pow(avgSourceW / targetSourceW, 1.0 / nDim);
          for (int iDim = 0; iDim < nDim; iDim++) {
            ppc[iDim] = std::max(1, int(ppc[iDim] * ratio));
          }
        }
      }

      Vector<Real> weights;
      weights.resize(sps.size());
      for (size_t i = 0; i < sps.size(); ++i) {
        weights[i] = sps[i].rdata(iqp_);
      }

      std::vector<int> idx =
          random_select_weighted_n(weights, product(ppc), randNum);

      Real wTmp = 0;
      for (int i : idx) {
        wTmp += sps[i].rdata(iqp_);
      }
      Real scale = rhoSource / wTmp;

      for (int i : idx) {
        ParticleType newp;
        set_ids(newp);

        newp.rdata(iqp_) = sps[i].rdata(iqp_) * scale;
        for (int iDim = 0; iDim < nDim; iDim++) {
          newp.rdata(iup_ + iDim) = sps[i].rdata(iup_ + iDim);
          newp.pos(ix_ + iDim) = sps[i].pos(ix_ + iDim);
        }
        pTile.push_back(newp);
      }
    }
  }
}

template <int NStructReal, int NStructInt>
Real Particles<NStructReal, NStructInt>::sum_moments_new(
    Vector<MultiFab>& momentsMF, Vector<MultiFab>& nodeBMF, Real dt,
    amrex::Vector<amrex::iMultiFab>& nodestatus) {

  //////////////////// This part will be removed when Nodestatus - bit :
  /// is_lev_edge is fixed for coarse level
  amrex::iMultiFab ctmp(nodestatus[0].boxArray(),
                        nodestatus[0].DistributionMap(), 1, 0);
  amrex::iMultiFab ftmp(nodestatus[1].boxArray(),
                        nodestatus[1].DistributionMap(), 1, 0);
  for (amrex::MFIter mfi(nodestatus[1]); mfi.isValid(); ++mfi) {
    const auto& status = nodestatus[1][mfi].array();
    const auto& ft = ftmp[mfi].array();

    ParallelFor(mfi.validbox(), [&](int i, int j, int k) {
      const amrex::IntVect ijk = { AMREX_D_DECL(i, j, k) };

      if (bit::is_lev_edge(status(ijk))) {
        ft(ijk) = 1;
      } else {
        ft(ijk) = 0;
      }
    });
  }
  ctmp.setVal(0);
  amrex::IntVect ratio = { AMREX_D_DECL(2, 2, 2) };
  average_down_nodal(ftmp, ctmp, ratio);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  Real energy = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    momentsMF[iLev].setVal(0.0);
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      Array4<Real> const& momentsArr = momentsMF[iLev][pti].array();
      // const auto& statusArr = nodestatus[iLev][pti].array();
      const auto& cstatus = ctmp[pti].array();

      const AoS& particles = pti.GetArrayOfStructs();

      const Dim3 lo = init_dim3(0);
      const Dim3 hi = init_dim3(1);

      // Print() << "iLev = " << iLev << std::endl;
      for (const auto& p : particles) {
        if (p.id() < 0)
          continue;

        // Print() << "p = " << p << std::endl;
        const Real up = p.rdata(iup_);
        const Real vp = p.rdata(ivp_);
        const Real wp = p.rdata(iwp_);
        const Real qp = p.rdata(iqp_);

        //-----calculate interpolate coef begin-------------
        IntVect loIdx;
        RealVect dShift;
        find_node_index(p.pos(), Geom(iLev).ProbLo(), Geom(iLev).InvCellSize(),
                        loIdx, dShift);
        Real coef[2][2][2];
        Real coef_finer[2][2][2];

        linear_interpolation_coef(dShift, coef);
        linear_interpolation_coef_finer(dShift, coef_finer);

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
          for (int kk = lo.z; kk <= hi.z; ++kk)
            for (int jj = lo.y; jj <= hi.y; ++jj)
              for (int ii = lo.x; ii <= hi.x; ++ii) {
                const IntVect ijk = { AMREX_D_DECL(
                    loIdx[ix_] + ii, loIdx[iy_] + jj, loIdx[iz_] + kk) };
                if (iLev == 0 && cstatus(ijk) == 1)

                {
                  momentsArr(ijk, iVar) +=
                      coef_finer[ii][jj][kk] * pMoments[iVar];
                }

                else {
                  momentsArr(ijk, iVar) += coef[ii][jj][kk] * pMoments[iVar];
                }
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

  for (int iLev = n_lev() - 2; iLev >= 0; iLev--) {
    sum_two_lev_interface_node(
        momentsMF[iLev], momentsMF[iLev + 1], 0, momentsMF[iLev].nComp(),
        get_ref_ratio(iLev), Geom(iLev), Geom(iLev + 1), node_status(iLev + 1));
  }

  // Correct domain edge nodes
  for (int iLev = 0; iLev < n_lev() - 1; iLev++) {
    interp_from_coarse_to_fine_for_domain_edge(
        momentsMF[iLev], momentsMF[iLev + 1], 0, momentsMF[iLev].nComp(),
        get_ref_ratio(iLev), Geom(iLev), Geom(iLev + 1), node_status(iLev + 1));
  }

  energy *= 0.5 * qomSign * get_mass();

  return energy;
}

// Since Particles is a template, it is necessary to explicitly instantiate
// with template arguments.
template class Particles<nPicPartReal, nPicPartInt>;
template class Particles<nPTPartReal, nPTPartInt>;
