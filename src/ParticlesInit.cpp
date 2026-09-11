#include <cstdlib>

#include <AMReX_Loop.H>
#include <AMReX_ParReduce.H>

#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

//==========================================================
/**
 * @brief Adds particles to a specific cell in the grid.
 *
 * This function adds particles to a specific cell in the grid at a given level
 * (`iLev`). It initializes the particles based on the provided fluid interface,
 * user-defined velocities, and other parameters.
 *
 * @tparam NStructReal Number of real components in the particle structure.
 * @tparam NStructInt Number of integer components in the particle structure.
 * @param iLev The level at which to add the particles.
 * @param mfi The MultiFab iterator for the current tile.
 * @param ijk The cell index where particles are to be added.
 * @param interface Pointer to the fluid interface used for initializing
 * particles.
 * @param doVacuumLimit Flag indicating whether to apply vacuum limit.
 * @param ppc Particles per cell.
 * @param tpVel User-defined velocity for initializing test particles.
 * @param dt Time step for density change rate.
 */
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_particles_cell(
    const int iLev, const MFIter& mfi, IntVect ijk,
    const FluidInterface* interface, bool doVacuumLimit, IntVect ppc,
    const Vel& tpVel, Real dt) {

  // If true, initialize the test particles with user defined velocities instead
  // of from fluid.
  bool userState = (tpVel.tag == speciesID);

  // If dt>0, it suggests the 'density' obtained from interface is actually the
  // density changing rate.
  if (dt <= 0)
    dt = 1;

  IntVect nPPC = nPartPerCell;
  if (!(ppc == 0)) {
    nPPC = ppc;
  }

  if (nPPC == 0)
    return;

  if (isTargetPPCDefined && !isFake2D) {
    const auto tppc = target_PPC(iLev)[mfi].array();
    for (int i = 0; i < nDim; ++i) {
      if (nDim > 2) {
        nPPC[i] = cbrt(tppc(ijk));
      } else {
        nPPC[i] = sqrt(tppc(ijk));
      }
    }
  }
  set_random_seed(iLev, ijk, nPPC);

  const Real vol = dx[iLev].product();
  const int npcel = product(nPPC);

  const Real vol2Npcel = qomSign * vol / npcel;

  ParticleTileType& particles = get_particle_tile(iLev, mfi, ijk);

  int icount = 0;
  // Loop over particles inside grid cell i, j, k

  const int kmax = nDim > 2 ? nPPC[iz_] : 1;

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

        const Real zp = (nDim > 2) ? xyz[2] : 0.0;

        const Real nDens =
            (userState && tpVel.nDens >= 0.0)
                ? tpVel.nDens
                : interface->get_number_density(mfi, xyz0, speciesID, iLev);

        if (doVacuumLimit && nDens * dt < vacuum)
          continue;

        Real q = vol2Npcel * nDens;

        // Per-particle weight modification (e.g. ion-acoustic-wave density
        // perturbation). Routed through the InitialCondition plugin; runs
        // before the q != 0 guard so an IC can suppress a particle.
        if (ic_ && ic_->modifies_weights()) {
          ParticleICState pics;
          pics.iLev = iLev;
          pics.iSpec = speciesID;
          pics.iCount = icount;
          pics.nPerCell = npcel;
          pics.x = xyz[0];
          pics.y = xyz[1];
          pics.z = zp;
          pics.q = q;
          ic_->modify_particle_weight(pics);
          q = pics.q;
        }

        if (q != 0) {
          Real u, v, w;
          Real rand1 = randNum();
          Real rand2 = randNum();
          Real rand3 = randNum();
          Real rand4 = randNum();

          Real uth = (userState ? tpVel.vth : -1);

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

          // Per-particle THERMAL-velocity override (e.g. a bi-Maxwellian with
          // distinct T_perp / T_par for the proton-cyclotron anisotropy
          // instability). Routed through the InitialCondition plugin; runs
          // AFTER the isotropic thermal velocity is sampled and BEFORE the bulk
          // velocity is added, so the IC may fully replace the sampled thermal
          // velocity with its own anisotropic draw.
          if (ic_ && ic_->modifies_thermal_velocity()) {
            ParticleICState pics;
            pics.iLev = iLev;
            pics.iSpec = speciesID;
            pics.iCount = icount;
            pics.nPerCell = npcel;
            pics.x = xyz[0];
            pics.y = xyz[1];
            pics.z = zp;
            pics.uThermal = u;
            pics.vThermal = v;
            pics.wThermal = w;
            pics.rand[0] = rand1;
            pics.rand[1] = rand2;
            pics.rand[2] = rand3;
            pics.rand[3] = rand4;
            ic_->modify_particle_thermal_velocity(pics);
            u = pics.uThermal;
            v = pics.vThermal;
            w = pics.wThermal;
          }

          Real uBulk = userState ? tpVel.vx
                                 : interface->get_ux(mfi, xyz, speciesID, iLev);
          Real vBulk = userState ? tpVel.vy
                                 : interface->get_uy(mfi, xyz, speciesID, iLev);
          Real wBulk = userState ? tpVel.vz
                                 : interface->get_uz(mfi, xyz, speciesID, iLev);

          // Per-particle velocity / weight modification (e.g. beam bulk
          // override, or the hybrid-wave Alfven velocity kick). Routed through
          // the InitialCondition plugin; runs after the bulk velocity is
          // computed and before it is added to the thermal velocity.
          if (ic_ && ic_->modifies_velocities()) {
            ParticleICState pics;
            pics.iLev = iLev;
            pics.iSpec = speciesID;
            pics.iCount = icount;
            pics.nPerCell = npcel;
            pics.x = xyz[0];
            pics.y = xyz[1];
            pics.z = zp;
            pics.uBulk = uBulk;
            pics.vBulk = vBulk;
            pics.wBulk = wBulk;
            pics.qScale = 1.0;
            pics.charge = charge;
            pics.temperature =
                interface->get_uth_iso(mfi, xyz, speciesID, iLev);
            pics.temperature *= pics.temperature * mass;
            ic_->modify_particle_velocity(pics);
            uBulk = pics.uBulk;
            vBulk = pics.vBulk;
            wBulk = pics.wBulk;
            q *= pics.qScale;
          }
          // Wave boundary: add the velocity kick to particles in a wave cell.
          if (waveVelocityKick) {
            const int loX = mfi.validbox().smallEnd(ix_);
            const int hiX = mfi.validbox().bigEnd(ix_);
            const int loY = mfi.validbox().smallEnd(iy_);
            const int hiY = mfi.validbox().bigEnd(iy_);
            const int loZ = (nDim > 2) ? mfi.validbox().smallEnd(iz_) : 0;
            const int hiZ = (nDim > 2) ? mfi.validbox().bigEnd(iz_) : 0;
            // Driven by the field-side wave faces (set from Pic::bcField),
            // not by a particle-side spelling.
            const bool onWaveX = (isWaveFace[0] && ijk[ix_] < loX) ||
                                 (isWaveFace[1] && ijk[ix_] > hiX);
            const bool onWaveY = (isWaveFace[2] && ijk[iy_] < loY) ||
                                 (isWaveFace[3] && ijk[iy_] > hiY);
            const bool onWaveZ =
                (nDim > 2) && ((isWaveFace[4] && ijk[iz_] < loZ) ||
                               (isWaveFace[5] && ijk[iz_] > hiZ));
            if (onWaveX || onWaveY || onWaveZ) {
              const Real ppos[3] = { xyz[0], xyz[1], zp };
              const Real tNow = tc ? tc->get_time() : 0.0;
              Real dvx = 0, dvy = 0, dvz = 0;
              waveVelocityKick(ppos, tNow, dvx, dvy, dvz);
              uBulk += dvx;
              vBulk += dvy;
              wBulk += dvz;
            }
          }
          u += uBulk;
          v += vBulk;
          w += wBulk;

          auto p = make_particle();
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
            const auto& status = host_cell_status(iLev)[mfi].array();
            if (bit::is_refined(status(i, j, k)))
              continue;

            IntVect ijk = { AMREX_D_DECL(i, j, k) };

            bool doAdd = true;
#ifdef _PT_COMPONENT_
            if (stateOH && doSelectRegion) {
              const int iRegion =
                  stateOH->get_neu_source_region(mfi, ijk, iLev);
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
  int iLevMax = 0;
  if (tc->get_cycle() == 0) {
    iLevMax = n_lev() - 1;
  }
  for (int iLev = 0; iLev <= iLevMax; iLev++) {
    for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {

      const auto& status = host_cell_status(iLev)[mfi].array();

      // Host-only kernel: CPU particle allocation into ParticleContainer
      amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) noexcept {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };
        if (bit::is_new(status(ijk)) && !bit::is_refined(status(ijk))) {
          add_particles_cell(iLev, mfi, ijk, fi, true);
        }
      });
    }
  }
}

//==========================================================
// Return true if ghost-cell particle injection should be skipped for this BC.
// - outflow: Ghost particles would be folded into edge cells in
//   sum_moments_cell_centered(), causing double counting.
// - inflow: Flux is injected at the physical face by
//   inject_flux_at_inflow_faces().

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
        auto newp = make_particle();
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

// Explicit template instantiations.
template class Particles<nPicPartReal, nPicPartInt>;
template class Particles<nPTPartReal, nPTPartInt>;
