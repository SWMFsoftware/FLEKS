#include "Particles.h"

using namespace amrex;

Particles::Particles(const Geometry& a_geom, const DistributionMapping& a_dmap,
                     const BoxArray& a_ba, const int a_species_id,
                     const Real a_charge, const Real a_mass)
    : ParticleContainer<4, 0, 0, 0>(a_geom, a_dmap, a_ba),
      m_species_id(a_species_id),
      m_charge(a_charge),
      m_mass(a_mass) {}

void Particles::add_particles(const IntVect& a_num_particles_per_cell,
                              const Real a_thermal_momentum_std,
                              const Real a_thermal_momentum_mean,
                              const Real a_density, const RealBox& a_bounds,
                              const int a_problem) {
  BL_PROFILE("Particles::add_particles");

  
}

void Particles::sum_moments(const MultiFab& Ex, const MultiFab& Ey,
                                       const MultiFab& Ez, const MultiFab& Bx,
                                       const MultiFab& By, const MultiFab& Bz,
                                       MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                       Real dt) {
  BL_PROFILE("Particles::sum_moments");

}

void Particles::mover(const MultiFab& Ex, const MultiFab& Ey,
                                    const MultiFab& Ez, const MultiFab& Bx,
                                    const MultiFab& By, const MultiFab& Bz,
                                    Real dt) {
  BL_PROFILE("Particles::mover");

}
