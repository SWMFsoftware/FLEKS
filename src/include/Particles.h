#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include <AMReX_Particles.H>

class ParticlesIter : public amrex::ParIter<4, 0, 0, 0> {
public:
  using amrex::ParIter<4, 0, 0, 0>::ParIter;
};

class Particles : public amrex::ParticleContainer<4, 0, 0, 0> {

public:
  Particles(const amrex::Geometry& a_geom,
            const amrex::DistributionMapping& a_dmap,
            const amrex::BoxArray& a_ba, const int a_species_id,
            const amrex::Real a_charge, const amrex::Real a_mass);

  void add_particles(const amrex::IntVect& a_num_particles_per_cell,
                     const amrex::Real a_thermal_momentum_std,
                     const amrex::Real a_thermal_momentum_mean,
                     const amrex::Real a_density,
                     const amrex::RealBox& a_bounds, const int a_problem);

  void sum_moments(const amrex::MultiFab& Ex, const amrex::MultiFab& Ey,
                   const amrex::MultiFab& Ez, const amrex::MultiFab& Bx,
                   const amrex::MultiFab& By, const amrex::MultiFab& Bz,
                   amrex::MultiFab& jx, amrex::MultiFab& jy,
                   amrex::MultiFab& jz, amrex::Real dt);

  void mover(const amrex::MultiFab& Ex, const amrex::MultiFab& Ey,
             const amrex::MultiFab& Ez, const amrex::MultiFab& Bx,
             const amrex::MultiFab& By, const amrex::MultiFab& Bz,
             amrex::Real dt);

protected:
  int m_species_id;

  amrex::Real m_charge;
  amrex::Real m_mass;
};

#endif
