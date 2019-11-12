#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include <AMReX_Particles.H>

#include "Array1D.h"
#include "Constants.h"
#include "FluidPicInterface.h"
#include "RandNum.h"
#include "UMultiFab.h"

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

  void sum_moments(amrex::MultiFab& momentsMF, amrex::UMultiFab<RealMM>& nodeMM,
                   amrex::MultiFab& nodeBMF, amrex::Real dt);

  void sum_to_center(amrex::MultiFab& netChargeMF,
                     amrex::UMultiFab<RealCMM>& centerMM, bool doNetChargeOnly);

  void mover(const amrex::MultiFab& nodeEMF, const amrex::MultiFab& nodeBMF,
             amrex::Real dt);

  void divE_correct_position(const amrex::MultiFab& phiMF);

protected:
  static const int iup_ = 0;
  static const int ivp_ = 1;
  static const int iwp_ = 2;
  static const int iqp_ = 3;

  int speciesID;
  RandNum randNum;

  amrex::Real charge;
  amrex::Real mass;

  amrex::IntVect nPartPerCell;
};

#endif