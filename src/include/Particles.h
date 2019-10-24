#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include <AMReX_Particles.H>

#include "Constants.h"
#include "FluidPicInterface.h"
#include "RandNum.h"

class ParticlesIter : public amrex::ParIter<4, 0, 0, 0> {
public:
  using amrex::ParIter<4, 0, 0, 0>::ParIter;
};

class Particles : public amrex::ParticleContainer<4, 0, 0, 0> {

public:
  Particles(const amrex::Geometry& geom, const amrex::DistributionMapping& dm,
            const amrex::BoxArray& ba, const int speciesID,
            const amrex::Real charge, const amrex::Real mass,
            const amrex::IntVect& nPartPerCellIn);

  void add_particles_domain(const FluidPicInterface& fluidInterface);
  void add_particles_cell(const amrex::MFIter& mfi,
                          const FluidPicInterface& fluidInterface, int iBlock,
                          int i, int j, int k, int loi, int loj, int lok);

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
  static const int iup_=0;
  static const int ivp_=1;
  static const int iwp_=2;
  static const int iqp_=3; 

  int speciesID;
  RandNum randNum;

  amrex::Real charge;
  amrex::Real mass;

  amrex::IntVect nPartPerCell;
};

#endif
