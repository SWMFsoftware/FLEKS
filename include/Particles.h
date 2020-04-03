#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include <memory>

#include <AMReX_CoordSys.H>
#include <AMReX_Particles.H>

#include "Array1D.h"
#include "Constants.h"
#include "FluidInterface.h"
#include "RandNum.h"
#include "TimeCtr.h"
#include "UMultiFab.h"

class PartInfo {
public:
  amrex::Real energy;
  amrex::Real uMax;
  PartInfo() : energy(0), uMax(0) {}
};

class ParticlesIter : public amrex::ParIter<4, 0, 0, 0> {
public:
  using amrex::ParIter<4, 0, 0, 0>::ParIter;
};

class Particles : public amrex::ParticleContainer<4, 0, 0, 0> {
private:
  amrex::BoxArray regionBA;

  amrex::Vector<amrex::RealBox> boxRange_I;

public:
  Particles(const amrex::Geometry& geom, const amrex::DistributionMapping& dm,
            const amrex::BoxArray& ba, TimeCtr* const tcIn, const int speciesID,
            const amrex::Real charge, const amrex::Real mass,
            const amrex::IntVect& nPartPerCellIn);

  void set_region_ba(const amrex::BoxArray& in) {
    regionBA = in;
    boxRange_I.clear();
    for (int iBox = 0; iBox < regionBA.size(); iBox++) {
      amrex::RealBox rb(regionBA[iBox], Geom(0).CellSize(), Geom(0).Offset());
      boxRange_I.push_back(rb);
    }

    // for(auto&rb:boxRange_I){
    //   amrex::Print()<<"rb = "<<rb<<std::endl;
    // }
  }

  void add_particles_domain(const FluidInterface& fluidInterface,
                            const amrex::iMultiFab& cellStatus);
  void add_particles_cell(const amrex::MFIter& mfi,
                          const FluidInterface& fluidInterface, int i, int j,
                          int k);
  void inject_particles_at_boundary(const FluidInterface& fluidInterface,
                                    const amrex::iMultiFab& cellStatus);

  // 1) Only inject particles ONCE for one ghost cells. This function decides
  // which block injects particles. 2) bx should be a valid box 3) The cell
  // (i,j,k) can NOT be the outmost ghost cell layer!!!!
  bool do_inject_particles_for_this_cell(const amrex::Box& bx,
                                         const amrex::Array4<const int>& status,
                                         const int i, const int j, const int k);

  PartInfo sum_moments(amrex::MultiFab& momentsMF,
                       amrex::UMultiFab<RealMM>& nodeMM,
                       amrex::MultiFab& nodeBMF, amrex::Real dt);

  // It is real 'thermal velocity'. It is sqrt(sum(q*v2)/sum(q)).
  amrex::Real calc_max_thermal_velocity(amrex::MultiFab& momentsMF);

  void sum_to_center(amrex::MultiFab& netChargeMF,
                     amrex::UMultiFab<RealCMM>& centerMM, bool doNetChargeOnly);

  void mover(const amrex::MultiFab& nodeEMF, const amrex::MultiFab& nodeBMF,
             amrex::Real dt);

  void convert_to_fluid_moments(amrex::MultiFab& momentsMF);

  inline bool is_outside_domain(const ParticleType& p) {
    const auto& plo = Geom(0).ProbLo();
    const auto& phi = Geom(0).ProbHi();
    const auto& dx = Geom(0).CellSize();

    for (int iDim = 0; iDim < nDimMax; iDim++) {
      if (!Geom(0).isPeriodic(iDim)) {
        if (p.pos(iDim) > phi[iDim] || p.pos(iDim) < plo[iDim]) {
          return true;
        }
      }
    }

    return false;
  }

  inline bool is_outside_ba(const ParticleType& p) {
    if (is_outside_domain(p))
      return true;

    const auto& plo = Geom(0).ProbLo();
    const auto& phi = Geom(0).ProbHi();
    amrex::Real loc[3] = { 0, 0, 0 };
    for (int iDim = 0; iDim < 3; iDim++) {
      loc[iDim] = p.pos(iDim);
      if (Geom(0).isPeriodic(iDim)) {
        // Fix index/loc for periodic BC.
        while (loc[iDim] > phi[iDim])
          loc[iDim] -= phi[iDim] - plo[iDim];
        while (loc[iDim] < plo[iDim])
          loc[iDim] += phi[iDim] - plo[iDim];
      }
    }

    for (const auto& rb : boxRange_I) {
      if (rb.contains(loc))
        return false;
    }

    return true;
  }

  void label_particles_outside_ba() {
    const int lev = 0;
    for (ParticlesIter pti(*this, lev); pti.isValid(); ++pti) {
      auto& particles = pti.GetArrayOfStructs();
      for (auto& p : particles) {
        if (is_outside_ba(p)) {
          p.id() = -1;
          // amrex::Print()<<"particle outside ba = "<<p<<std::endl;
        }
      }
    }
  }

  void split_particles(amrex::Real limit);
  void combine_particles(amrex::Real limit);

  void divE_correct_position(const amrex::MultiFab& phiMF);

  amrex::Real get_qom() { return charge / mass; }

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
