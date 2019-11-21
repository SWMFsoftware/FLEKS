#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include <AMReX_Particles.H>

#include "Array1D.h"
#include "Constants.h"
#include "FluidInterface.h"
#include "RandNum.h"
#include "UMultiFab.h"
#include "TimeCtr.h"

class ParticlesIter : public amrex::ParIter<4, 0, 0, 0> {
public:
  using amrex::ParIter<4, 0, 0, 0>::ParIter;
};

class Particles : public amrex::ParticleContainer<4, 0, 0, 0> {

public:
  Particles(const amrex::Geometry& geom, const amrex::DistributionMapping& dm,
            const amrex::BoxArray& ba, TimeCtr *const tcIn, const int speciesID,
            const amrex::Real charge, const amrex::Real mass,
            const amrex::IntVect& nPartPerCellIn);

  void add_particles_domain(const FluidInterface& fluidInterface);
  void add_particles_cell(const amrex::MFIter& mfi,
                          const FluidInterface& fluidInterface, int i, int j,
                          int k);
  void inject_particles_at_boundary(const FluidInterface& fluidInterface);

  void sum_moments(amrex::MultiFab& momentsMF, amrex::UMultiFab<RealMM>& nodeMM,
                   amrex::MultiFab& nodeBMF, amrex::Real dt);

  void sum_to_center(amrex::MultiFab& netChargeMF,
                     amrex::UMultiFab<RealCMM>& centerMM, bool doNetChargeOnly);

  void mover(const amrex::MultiFab& nodeEMF, const amrex::MultiFab& nodeBMF,
             amrex::Real dt);

  void convert_to_fluid_moments(amrex::MultiFab& momentsMF);

  inline bool is_outside(const ParticleType& p) {
    const auto& plo = Geom(0).ProbLo();
    const auto& phi = Geom(0).ProbHi();
    const auto& dx = Geom(0).CellSize();

    for (int iDim = 0; iDim < nDimMax; iDim++) {
      if (!Geom(0).isPeriodic(iDim)) {
        if (p.pos(iDim) > phi[iDim] - nVirGst * dx[iDim] ||
            p.pos(iDim) < plo[iDim] + nVirGst * dx[iDim]) {
          return true;
        }
      }
    }

    return false;
  }

  void divE_correct_position(const amrex::MultiFab& phiMF);

  amrex::Real get_qom(){
    return charge/mass;
  }
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

  TimeCtr* tc; 

};

#endif
