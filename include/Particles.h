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

struct PID {
  int cpu;
  int id;
  bool flag;

  // This function is used by c++ STL algorithms.
  bool operator<(const PID& t) const {
    bool lt = cpu < t.cpu;

    if (cpu == t.cpu)
      lt = id < t.id;

    return lt;
  }

  bool operator==(const PID& t) const { return cpu == t.cpu && id == t.id; }
};

template <int NStructReal = nPicPartReal, int NStructInt = 0>
class ParticlesIter : public amrex::ParIter<NStructReal, NStructInt> {
public:
  using amrex::ParIter<NStructReal, NStructInt>::ParIter;
};

template <int NStructReal = nPicPartReal, int NStructInt = 0>
class Particles : public amrex::ParticleContainer<NStructReal, NStructInt> {
public:
  static ParticleStaggering particlePosition;

public:
  // Since this is a template, the compiler will not search names in the base
  // class by default, and the following 'using ' statements are required.
  using ParticleType = amrex::Particle<NStructReal, NStructInt>;
  using amrex::ParticleContainer<NStructReal, NStructInt>::Geom;
  using amrex::ParticleContainer<NStructReal, NStructInt>::do_tiling;
  using amrex::ParticleContainer<NStructReal, NStructInt>::tile_size;
  using amrex::ParticleContainer<NStructReal, NStructInt>::SetUseUnlink;
  using amrex::ParticleContainer<NStructReal, NStructInt>::GetParticles;
  using amrex::ParticleContainer<NStructReal, NStructInt>::MakeMFIter;
  using amrex::ParticleContainer<NStructReal, NStructInt>::Redistribute;

protected:
  FluidInterface* fluidInterface;
  TimeCtr* tc;

  int speciesID;
  RandNum randNum;

  amrex::Real charge;
  amrex::Real mass;

  amrex::Real qom;
  int qomSign;

  amrex::IntVect nPartPerCell;

  amrex::BoxArray regionBA;

  amrex::Vector<amrex::RealBox> boxRange_I;

  amrex::Real plo[nDim], phi[nDim], dx[nDim], invDx[nDim];
  amrex::Real invVol;
  bool isPeriodic[nDim];

  amrex::Real mergeThresholdDistance = 0.6;
  amrex::Real velBinBufferSize = 0.125;

  bool isRelativistic = false;

public:
  static const int iup_ = 0;
  static const int ivp_ = 1;
  static const int iwp_ = 2;
  static const int iqp_ = 3;

  TestCase testCase;

  // Index of the integer data.
  static const int iRecordCount_ = 0;

  amrex::iMultiFab cellStatus;

  Particles(const amrex::BoxArray& regionBAIn, const amrex::Geometry& geom,
            const amrex::DistributionMapping& dm, const amrex::BoxArray& ba,
            FluidInterface* fluidIn, TimeCtr* tcIn, const int speciesIDIn,
            const amrex::Real chargeIn, const amrex::Real massIn,
            const amrex::IntVect& nPartPerCellIn,
            TestCase tcase = RegularSimulation);

  amrex::BoxArray get_region_ba() const { return regionBA; }

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

  void add_particles_domain(const amrex::iMultiFab& cellStatus);
  void add_particles_cell(const amrex::MFIter& mfi, const int i, const int j,
                          const int k, const FluidInterface& interface,
                          amrex::Real ratio = 1);
  void inject_particles_at_boundary(const amrex::iMultiFab& cellStatus);

  void add_particles_source(const amrex::MultiFab& momentsMF,
                            const FluidInterface& interface);

  // 1) Only inject particles ONCE for one ghost cells. This function decides
  // which block injects particles. 2) bx should be a valid box 3) The cell
  // (i,j,k) can NOT be the outmost ghost cell layer!!!!
  bool do_inject_particles_for_this_cell(const amrex::Box& bx,
                                         const amrex::Array4<const int>& status,
                                         const int i, const int j, const int k);

  amrex::Real sum_moments(amrex::MultiFab& momentsMF,
                          amrex::UMultiFab<RealMM>& nodeMM,
                          amrex::MultiFab& nodeBMF, amrex::Real dt);

  void calc_mass_matrix(amrex::UMultiFab<RealMM>& nodeMM, amrex::MultiFab& jHat,
                        amrex::MultiFab& nodeBMF, amrex::Real dt);

  // It is real 'thermal velocity'. It is sqrt(sum(q*v2)/sum(q)).
  amrex::Real calc_max_thermal_velocity(amrex::MultiFab& momentsMF);

  void sum_to_center(amrex::MultiFab& netChargeMF,
                     amrex::UMultiFab<RealCMM>& centerMM, bool doNetChargeOnly);

  void mover(const amrex::MultiFab& nodeEMF, const amrex::MultiFab& nodeBMF,
             amrex::Real dt, amrex::Real dtNext);

  void update_position_to_half_stage(const amrex::MultiFab& nodeEMF,
                                     const amrex::MultiFab& nodeBMF,
                                     amrex::Real dt);

  void convert_to_fluid_moments(amrex::MultiFab& momentsMF);

  void set_ppc(amrex::IntVect& in) { nPartPerCell = in; };

  inline bool is_outside_ba(const ParticleType& p) {
    amrex::Real loc[3] = { 0, 0, 0 };
    for (int iDim = 0; iDim < 3; iDim++) {
      loc[iDim] = p.pos(iDim);
      if (isPeriodic[iDim]) {
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

  inline bool is_outside_ba(const ParticleType& p,
                            amrex::Array4<int const> const& status,
                            const amrex::IntVect& low,
                            const amrex::IntVect& high) {

    // Contains ghost cells.
    bool isInsideBox = true;
    int cellIdx[3];
    amrex::Real dShift[3];
    for (int i = 0; i < 3; i++) {
      dShift[i] = (p.pos(i) - plo[i]) * invDx[i];
      cellIdx[i] = fastfloor(dShift[i]);
      if (cellIdx[i] > high[i] || cellIdx[i] < low[i]) {
        isInsideBox = false;
        break;
      }
    }

    if (isInsideBox) {
      return status(cellIdx[ix_], cellIdx[iy_], cellIdx[iz_]) == iBoundary_;
    } else {
      return is_outside_ba(p);
    }
  }

  void label_particles_outside_ba() {
    const int lev = 0;
    for (ParticlesIter<NStructReal, NStructInt> pti(*this, lev); pti.isValid();
         ++pti) {
      auto& particles = pti.GetArrayOfStructs();
      const amrex::Array4<int const>& status = cellStatus[pti].array();
      const amrex::Box& bx = cellStatus[pti].box();
      const amrex::IntVect lowCorner = bx.smallEnd();
      const amrex::IntVect highCorner = bx.bigEnd();
      for (auto& p : particles) {
        if (is_outside_ba(p, status, lowCorner, highCorner)) {
          p.id() = -1;
          // amrex::Print()<<"particle outside ba = "<<p<<std::endl;
        }
      }
    }
  }

  void label_particles_outside_ba_general() {
    const int lev = 0;
    for (ParticlesIter<NStructReal, NStructInt> pti(*this, lev); pti.isValid();
         ++pti) {
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

  int get_speciesID() const { return speciesID; }
  amrex::Real get_charge() const { return charge; }
  amrex::Real get_mass() const { return mass; }

  void set_merge_threshold(amrex::Real in) { mergeThresholdDistance = in; }
  void set_merge_velocity_bin_buffer(amrex::Real in) { velBinBufferSize = in; }

  void set_relativistic(const bool& in) { isRelativistic = in; }
};

template <int NStructReal, int NStructInt>
ParticleStaggering Particles<NStructReal, NStructInt>::particlePosition =
    Staggered;

class IOParticles : public Particles<nPicPartReal> {
public:
  IOParticles() = delete;

  IOParticles(Particles<>& other, amrex::Geometry geomIO,
              amrex::Real no2outL = 1, amrex::Real no2outV = 1,
              amrex::Real no2OutM = 1,
              amrex::RealBox IORange = amrex::RealBox());
  ~IOParticles() = default;
};
#endif
