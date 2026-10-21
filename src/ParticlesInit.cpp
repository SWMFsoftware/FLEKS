#include <cstdlib>

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
 * @param weightScale Factor applied to every particle weight (the weight
 * floor of the source particles, see add_particles_source()).
 * @param allowBody Also fill a cell inside the inner body (the ionospheric
 * ghost layer of inject_particles_at_body()).
 */
template <int NStructReal, int NStructInt>
int Particles<NStructReal, NStructInt>::add_particles_cell(
    const int iLev, const MFIter& mfi, IntVect ijk,
    const FluidInterface* interface, bool doVacuumLimit, IntVect ppc,
    const Vel& tpVel, Real dt, Real weightScale, bool allowBody) {

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
    return 0;

  // No particle is created inside the absorbing inner body (see #BODY). The
  // test uses the cell center, i.e., the same cell-based staircase as the
  // particle absorption and the field mask. This covers the initial fill, the
  // source particles and the boundary (ghost) injection, which all funnel
  // through this function.
  if (grid != nullptr && grid->use_body() && !allowBody) {
    Real cellCenter[nDim];
    for (int iDim = 0; iDim < nDim; iDim++)
      cellCenter[iDim] = (ijk[iDim] + 0.5) * dx[iLev][iDim] + plo[iLev][iDim];

    if (grid->is_inside_body(cellCenter))
      return 0;
  }

  // The level target (#CONSTANTPPV, pre-splitting) applies only when the
  // caller does not ask for an explicit count (the source particles do).
  if (ppc == 0 && isTargetPPCDefined && !isFake2D) {
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
          p.rdata(iqp_) = q * dt * weightScale;

          if (NStructInt > iRecordCount_) {
            // For test particle only.
            p.idata(iRecordCount_) = 0;
          }

          particles.push_back(p);

          icount++;
        }
      }
  return icount;
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_particles_source(
    const FluidInterface* interface, const FluidInterface* const stateOH,
    Real dt, IntVect ppc, const bool doSelectRegion, const bool adaptivePPC,
    const Real minWeightRatio) {
  timing_func("Pts::add_particles_source");

  // Weight floor (#SOURCEWEIGHTFLOOR): a source particle is not lighter than
  // minWeightRatio times the weight of a particle of the reference plasma
  // (all ions of fi, i.e. #UNIFORMSTATE stand-alone). A cell whose whole
  // source of the step is below the floor gets one particle of the floor
  // weight with the probability source/floor ("Russian roulette"): the mass
  // is conserved on average, and the weak sources that cover most of an
  // extended exosphere no longer create a particle in every cell and step.
  const bool useWeightFloor = minWeightRatio > 0.0;
  const Real initPPCFloor = product(nPartPerCell);

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
              const int iRegion =
                  stateOH->get_neu_source_region(mfi, ijk, iLev);
              doAdd = (iRegion == speciesID);
            }
#endif
            if (doAdd) {
              IntVect ppcCell = ppc;
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
                    ppcCell[iDim] = std::max(1, int(ppc[iDim] * ratio));
                  }
                }
              }

              Real weightScale = 1.0;
              if (useWeightFloor) {
                RealVect xyzC;
                for (int iDim = 0; iDim < nDim; iDim++)
                  xyzC[iDim] =
                      (ijk[iDim] + 0.5) * dx[iLev][iDim] + plo[iLev][iDim];

                // Number density of the step's source and of the reference
                // plasma, in code units.
                Real nSource =
                    interface->get_number_density(mfi, xyzC, speciesID, iLev);
                if (dt > 0)
                  nSource *= dt;
                Real nRef = 0;
                for (int iS = 0; iS < fi->get_nS(); iS++)
                  if (fi->get_species_charge(iS) > 0)
                    nRef += fi->get_number_density(mfi, xyzC, iS, iLev);
                // The reference particle of the cell: with #CONSTANTPPV a
                // refined cell has fewer particles per cell than the base.
                const Real ppcRef =
                    (isTargetPPCDefined && !isFake2D)
                        ? static_cast<Real>(target_PPC(iLev)[mfi].array()(ijk))
                        : initPPCFloor;
                const Real wMin = minWeightRatio * nRef / ppcRef;

                if (nSource <= 0.0)
                  continue;
                if (nSource / product(ppcCell) < wMin) {
                  ppcCell = IntVect(1);
                  if (nSource < wMin) {
                    set_random_seed(iLev, ijk, IntVect(919));
                    if (randNum() * wMin >= nSource)
                      continue;
                    weightScale = wMin / nSource;
                  }
                }
              }

              add_particles_cell(iLev, mfi, ijk, interface, false, ppcCell,
                                 Vel(), dt, weightScale);
            }
          }
    }
  }
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::inject_emission_shell(
    const Vector<MultiFab>& nAdd, Real vth, int nPPC, const Real* center) {
  timing_func("Pts::inject_emission_shell");

  if (nPPC <= 0)
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const Real vol = dx[iLev].product();
    for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
      const auto add = nAdd[iLev][mfi].const_array();
      const Box& box = mfi.validbox();
      const auto lo = lbound(box);
      const auto hi = ubound(box);
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            const Real n = add(i, j, k);
            if (n <= 0)
              continue;
            IntVect ijk = { AMREX_D_DECL(i, j, k) };
            set_random_seed(iLev, ijk, IntVect(271));
            ParticleTileType& tile = get_particle_tile(iLev, mfi, ijk);
            const Real q = qomSign * n * vol / nPPC;
            for (int ip = 0; ip < nPPC; ip++) {
              auto p = make_particle();
              set_ids(p);
              Real r[3] = { 0, 0, 0 };
              for (int iDim = 0; iDim < nDim; iDim++) {
                p.pos(iDim) = (ijk[iDim] + randNum()) * dx[iLev][iDim] +
                              plo[iLev][iDim];
                r[iDim] = p.pos(iDim) - center[iDim];
              }
              Real v[3];
              for (int iDir = 0; iDir < 3; iDir += 2) {
                const Real prob =
                    std::sqrt(-2.0 * std::log(1.0 - .999999999 * randNum()));
                const Real theta = 2.0 * M_PI * randNum();
                v[iDir] = vth * prob * std::cos(theta);
                if (iDir == 0)
                  v[1] = vth * prob * std::sin(theta);
              }
              // Outward half-Maxwellian: the emitted ions move away from the
              // planet (rhybrid IonosphereInjector).
              if (v[0] * r[0] + v[1] * r[1] + v[2] * r[2] < 0)
                for (auto& vi : v)
                  vi = -vi;
              p.rdata(iup_) = v[0];
              p.rdata(ivp_) = v[1];
              p.rdata(iwp_) = v[2];
              p.rdata(iqp_) = q;
              if (NStructInt > iRecordCount_)
                p.idata(iRecordCount_) = 0;
              tile.push_back(p);
            }
            bodyInjectCount += nPPC;
            bodyInjectMass += n * vol * mass;
          }
    }
  }
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::inject_particles_at_body(
    const std::function<Real(const Real*)>& nDensAt, Real vth) {
  timing_func("Pts::inject_particles_at_body");

  if (grid == nullptr || !grid->use_body())
    return;

  const Real* center = grid->get_body_center();

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const Real vol = dx[iLev].product();
    for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
      const auto& status = cell_status(iLev)[mfi].array();
      const Box& box = mfi.validbox();
      const auto lo = lbound(box);
      const auto hi = ubound(box);
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            const int st = status(i, j, k);
            // The surface layer: body cells that touch the plasma.
            if (!bit::is_body(st) || bit::is_body_interior(st) ||
                bit::is_refined(st))
              continue;

            IntVect ijk = { AMREX_D_DECL(i, j, k) };
            Real xyz[3] = { 0.0, 0.0, 0.0 };
            for (int iDim = 0; iDim < nDim; iDim++)
              xyz[iDim] = (ijk[iDim] + 0.5) * dx[iLev][iDim] + plo[iLev][iDim] -
                          center[iDim];

            const Real n = nDensAt(xyz);
            if (n <= 0)
              continue;

            Vel state;
            state.tag = speciesID;
            state.nDens = n;
            state.vth = vth;
            bodyInjectCount += add_particles_cell(
                iLev, mfi, ijk, fi, false, IntVect(), state, -1, 1.0, true);
            // Positive mass, as in body_absorb_tally().
            bodyInjectMass += n * vol * mass;
          }
    }
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_particles_domain(bool fillAll) {
  timing_func("Pts::add_particles_domain");
  // Full-domain seeding occurs on initial particle filling (fillAll=true, e.g.
  // cold start or first fill following a fluid-only restart).
  // Dynamic AMR regrids rely on AMReX Redistribute() to migrate existing
  // particles to newly refined levels; re-seeding iLev > 0 would duplicate
  // particles. After initial fill, only newly activated level-0 base domain
  // cells (if any) are seeded.
  const int nLev = fillAll ? n_lev() : 1;
  for (int iLev = 0; iLev < nLev; iLev++) {
    for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {

      const auto& status = cell_status(iLev)[mfi].array();
      ParallelFor(mfi.validbox(), [&](int i, int j, int k) noexcept {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };
        // A fluid-only restart can retain mesh cells without restoring any
        // particles. The first particle fill must seed every leaf cell.
        if ((fillAll || bit::is_new(status(ijk))) &&
            !bit::is_refined(status(ijk))) {
          add_particles_cell(iLev, mfi, ijk, fi, true);
        }
      });
    }
  }
}

//==========================================================
// Return true if ghost-cell particle injection should be skipped for this BC.
// - outflow: Ghost particles would be folded into edge cells in
//   sum_moments(), causing double counting.
// - inflow: Flux is injected at the physical face by
//   inject_flux_at_inflow_faces().

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_source_particles(
    std::unique_ptr<PicParticles>& sourcePart, IntVect ppc,
    const bool adaptivePPC) {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PicParticles::PIter pti(*sourcePart, iLev); pti.isValid(); ++pti) {
      // It is assumed the tile size is 1x1x1.
      Box bx = pti.tilebox();
      auto cellIdx = bx.smallEnd();

      ParticleTileType& pTile = get_particle_tile(iLev, pti, cellIdx);
      auto& sps = pti.GetArrayOfStructs();

      if (sps.size() == 0)
        continue;

      Real rhoSource = 0;
      for (auto& p : sps) {
        rhoSource += p.rdata(iqp_);
      }

      set_random_seed(iLev, cellIdx, IntVect(787));

      IntVect ppcCell = ppc;
      if (adaptivePPC) {
        // Adjust ppc so that the weight of the
        // source particles is not too small.

        Real rho = 0;
        for (auto& p : pTile.GetArrayOfStructs()) {
          rho += p.rdata(iqp_);
        }

        Real avgInitW = rho / product(nPartPerCell);
        Real avgSourceW = rhoSource / product(ppc);

        Real targetSourceW = avgInitW * 0.1;

        if (avgSourceW < targetSourceW) {
          Real ratio = pow(avgSourceW / targetSourceW, 1.0 / nDim);
          for (int iDim = 0; iDim < nDim; iDim++) {
            ppcCell[iDim] = std::max(1, int(ppc[iDim] * ratio));
          }
        }
      }

      Vector<Real> weights;
      weights.resize(sps.size());
      for (size_t i = 0; i < sps.size(); ++i) {
        weights[i] = sps[i].rdata(iqp_);
      }

      std::vector<int> idx =
          random_select_weighted_n(weights, product(ppcCell), randNum);

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

#define INSTANTIATE_PARTICLES_INIT(T)                                          \
  template void T::add_particles_domain(bool);                                 \
  template int T::add_particles_cell(const int, const MFIter&, const IntVect,  \
                                     const FluidInterface*, bool, IntVect,     \
                                     const Vel&, Real, Real, bool);            \
  template void T::inject_emission_shell(const Vector<MultiFab>&, Real, int,  \
                                         const Real*);                         \
  template void T::inject_particles_at_body(                                   \
      const std::function<Real(const Real*)>&, Real);                          \
  template void T::add_particles_source(                                       \
      const FluidInterface*, const FluidInterface* const, Real, IntVect,       \
      const bool, const bool, const Real);                                     \
  template void T::add_source_particles(std::unique_ptr<PicParticles>&,        \
                                        IntVect, const bool);

INSTANTIATE_PARTICLES_INIT(PicParticles)
INSTANTIATE_PARTICLES_INIT(PTParticles)
#undef INSTANTIATE_PARTICLES_INIT
