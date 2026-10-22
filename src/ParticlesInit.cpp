#include <cstdlib>

#include <AMReX_ParReduce.H>

#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "ParticlesSeed.h"
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
  // A single cell: see add_particles_cells().
  add_particles_cells(
      iLev, mfi, Box(ijk, ijk),
      [=] AMREX_GPU_HOST_DEVICE(const IntVect&, IntVect& ppcCell) {
        ppcCell = ppc;
        return true;
      },
      interface, doVacuumLimit, tpVel, dt);
}

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_particles_source(
    const FluidInterface* interface, const FluidInterface* const stateOH,
    Real dt, IntVect ppc, const bool doSelectRegion, const bool adaptivePPC) {
  timing_func("Pts::add_particles_source");

  const int initPPC = product(nPartPerCell);
  const int sourcePPC = product(ppc);
  const int iSpecies = speciesID;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
      const Box& tile_box = mfi.validbox();
      const Array4<int const> status = cell_status(iLev)[mfi].const_array();

      // Cells of the source region of this species (OH-PT only). The region
      // test runs on the host.
      Vector<int> inRegion;
      const bool useRegion = (stateOH && doSelectRegion);
#ifdef _PT_COMPONENT_
      if (useRegion) {
        if (fleksUseGpu)
          Abort("add_particles_source: the source region selection is not "
                "available in GPU runs yet.");
        inRegion.resize(tile_box.numPts());
        int l = 0;
        LoopOnCpu(tile_box, [&](int i, int j, int k) {
          const IntVect ijk = { AMREX_D_DECL(i, j, k) };
          inRegion[l++] =
              (stateOH->get_neu_source_region(mfi, ijk, iLev) == speciesID);
        });
      }
#endif
      Gpu::AsyncArray<int> regionStore(inRegion.data(), inRegion.size());
      const int* const region = useRegion ? regionStore.data() : nullptr;
      const Dim3 lo = lbound(tile_box);
      const Dim3 len = length(tile_box);

      const FluidView fvPic = fi->fluid_view(mfi, iLev);
      const FluidView fvSource = interface->fluid_view(mfi, iLev);

      add_particles_cells(
          iLev, mfi, tile_box,
          [=] AMREX_GPU_HOST_DEVICE(const IntVect& ijk, IntVect& ppcCell) {
            if (bit::is_refined(status(ijk)))
              return false;
            if (region) {
              const int l =
                  (ijk[0] - lo.x) +
                  len.x * ((AMREX_SPACEDIM > 1 ? ijk[1] - lo.y : 0) +
                           len.y * (AMREX_SPACEDIM > 2 ? ijk[2] - lo.z : 0));
              if (!region[l])
                return false;
            }

            ppcCell = ppc;
            if (adaptivePPC) {
              // Adjust ppc so that the weight of the source particles is not
              // too small.
              Real rho = fvPic.get_number_density(fvPic.d, ijk, iSpecies);
              Real rhoSource =
                  fvSource.get_number_density(fvSource.d, ijk, iSpecies);
              if (dt > 0)
                rhoSource *= dt;

              Real avgInitW = rho / initPPC;
              Real avgSourceW = rhoSource / sourcePPC;

              Real targetSourceW = avgInitW * 0.1;

              if (avgSourceW < targetSourceW) {
                Real ratio = std::pow(avgSourceW / targetSourceW, 1.0 / nDim);
                for (int iDim = 0; iDim < nDim; iDim++) {
                  ppcCell[iDim] = amrex::max(1, int(ppc[iDim] * ratio));
                }
              }
            }
            return true;
          },
          interface, false, Vel(), dt);
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

      const Array4<int const> status = cell_status(iLev)[mfi].const_array();
      add_particles_cells(
          iLev, mfi, mfi.validbox(),
          [=] AMREX_GPU_HOST_DEVICE(const IntVect& ijk, IntVect&) {
            // A fluid-only restart can retain mesh cells without restoring any
            // particles. The first particle fill must seed every leaf cell.
            return (fillAll || bit::is_new(status(ijk))) &&
                   !bit::is_refined(status(ijk));
          },
          fi, true);
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
  template void T::add_particles_cell(const int, const MFIter&, const IntVect, \
                                      const FluidInterface*, bool, IntVect,    \
                                      const Vel&, Real);                       \
  template void T::add_particles_source(const FluidInterface*,                 \
                                        const FluidInterface* const, Real,     \
                                        IntVect, const bool, const bool);      \
  template void T::add_source_particles(std::unique_ptr<PicParticles>&,        \
                                        IntVect, const bool);

INSTANTIATE_PARTICLES_INIT(PicParticles)
INSTANTIATE_PARTICLES_INIT(PTParticles)
#undef INSTANTIATE_PARTICLES_INIT
