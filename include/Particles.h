#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include <memory>

#include <AMReX_AmrCore.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_CoordSys.H>

#include "Array1D.h"
#include "BC.h"
#include "Bit.h"
#include "Constants.h"
#include "FluidInterface.h"
#include "GridUtility.h"
#include "RandNum.h"
#include "SourceInterface.h"
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

struct Vel {
  amrex::Real vth;
  amrex::Real vx;
  amrex::Real vy;
  amrex::Real vz;
  // tag is usually the species ID
  int tag;

  Vel() {
    vth = 0;
    vx = 0;
    vy = 0;
    vz = 0;
    tag = -1;
  }
};

template <int NStructReal = nPicPartReal, int NStructInt = 0>
class ParticlesIter : public amrex::ParIter<NStructReal, NStructInt> {
public:
  using amrex::ParIter<NStructReal, NStructInt>::ParIter;
};

template <int NStructReal = nPicPartReal, int NStructInt = 0>
class Particles : public amrex::AmrParticleContainer<NStructReal, NStructInt> {
public:
  static ParticleStaggering particlePosition;

public:
  // Since this is a template, the compiler will not search names in the base
  // class by default, and the following 'using ' statements are required.
  using ParticleType = amrex::Particle<NStructReal, NStructInt>;
  using ParticleTileType = amrex::ParticleTile<NStructReal, NStructInt, 0, 0>;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::Geom;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::do_tiling;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::tile_size;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::SetUseUnlink;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::GetParticles;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::MakeMFIter;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::Redistribute;
  using amrex::AmrParticleContainer<NStructReal,
                                    NStructInt>::NumberOfParticlesAtLevel;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::finestLevel;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::Checkpoint;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::Index;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::ParticlesAt;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::maxLevel;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::GetParGDB;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::ParticleBoxArray;
  using amrex::AmrParticleContainer<NStructReal,
                                    NStructInt>::ParticleDistributionMap;

protected:
  FluidInterface* fi;
  TimeCtr* tc;

  int speciesID;
  RandNum randNum;

  int nLev;

  amrex::Real charge;
  amrex::Real mass;

  amrex::Real qom;
  int qomSign;

  amrex::IntVect nPartPerCell;

  amrex::Vector<amrex::RealBox> activeRegions;

  amrex::Vector<amrex::RealVect> plo, phi, dx, invDx;
  amrex::Vector<amrex::Real> invVol;

  amrex::Real mergeThresholdDistance = 0.6;
  amrex::Real velBinBufferSize = 0.125;

  // If fastMerge == false: find the particle pair that is closest to each other
  // in the phase space and try to delete the lighter one.
  // If fastMerge == true: sort the particles by weights from light to heavy,
  // and try to delete one of them one by one.
  bool fastMerge = false;

  bool isRelativistic = false;

  BC bc; // boundary condition

  amrex::Vector<amrex::iMultiFab> cellStatus;

public:
  static const int iup_ = 0;
  static const int ivp_ = 1;
  static const int iwp_ = 2;
  static const int iqp_ = 3;

  TestCase testCase;

  // Index of the integer data.
  static const int iRecordCount_ = 0;

  Particles(amrex::AmrCore* amrcore, FluidInterface* fluidIn, TimeCtr* tcIn,
            const int speciesIDIn, const amrex::Real chargeIn,
            const amrex::Real massIn, const amrex::IntVect& nPartPerCellIn,
            TestCase tcase = RegularSimulation);

  void set_region_range(const amrex::BoxArray& ba) {
    activeRegions.clear();
    for (int iBox = 0; iBox < ba.size(); iBox++) {
      amrex::RealBox rb(ba[iBox], Geom(0).CellSize(), Geom(0).Offset());
      activeRegions.push_back(rb);
    }
  }

  void add_particles_domain(const amrex::Vector<amrex::iMultiFab>& iRefinement);
  void add_particles_cell(const int iLev, const amrex::MFIter& mfi, const int i,
                          const int j, const int k,
                          const FluidInterface& interface,
                          amrex::IntVect ppc = amrex::IntVect(),
                          const Vel tpVel = Vel(), amrex::Real dt = -1);
  void inject_particles_at_boundary(const FluidInterface* fiIn = nullptr,
                                    amrex::Real dt = -1,
                                    amrex::IntVect ppc = amrex::IntVect());

  void add_particles_source(const FluidInterface& interface,
                            const FluidInterface* const stateOH = nullptr,
                            amrex::Real dt = -1,
                            amrex::IntVect ppc = amrex::IntVect(),
                            const bool doSelectRegion = false);

  // Copy particles from (ip,jp,kp) to (ig, jg, kg) and shift boundary
  // particle's coordinates accordingly.
  void outflow_bc(const amrex::MFIter& mfi, const int ig, const int jg,
                  const int kg, const int ip, const int jp, const int kp);

  // 1) Only inject particles ONCE for one ghost cells. This function decides
  // which block injects particles. 2) bx should be a valid box 3) The cell
  // (i,j,k) can NOT be the outmost ghost cell layer!!!!
  bool do_inject_particles_for_this_cell(const amrex::Box& bx,
                                         const amrex::Array4<const int>& status,
                                         const int i, const int j, const int k,
                                         int& isrc, int& jsrc, int& ksrc);

  amrex::Real sum_moments(amrex::Vector<amrex::MultiFab>& momentsMF,
                          amrex::Vector<amrex::MultiFab>& nodeBMF,
                          amrex::Real dt);

  std::array<amrex::Real, 5> total_moments(bool localOnly = false);

  void calc_mass_matrix(amrex::UMultiFab<RealMM>& nodeMM, amrex::MultiFab& jHat,
                        amrex::MultiFab& nodeBMF, amrex::Real dt);

  void calc_jhat(amrex::MultiFab& jHat, amrex::MultiFab& nodeBMF,
                 amrex::Real dt);

  // It is real 'thermal velocity'. It is sqrt(sum(q*v2)/sum(q)).
  amrex::Real calc_max_thermal_velocity(amrex::MultiFab& momentsMF);

  void sum_to_center(amrex::MultiFab& netChargeMF,
                     amrex::UMultiFab<RealCMM>& centerMM, bool doNetChargeOnly);

  void charge_exchange(amrex::Real dt, FluidInterface* stateOH,
                       FluidInterface* sourcePT2OH, SourceInterface* source);

  void mover(const amrex::Vector<amrex::MultiFab>& nodeE,
             const amrex::Vector<amrex::MultiFab>& nodeB, amrex::Real dt,
             amrex::Real dtNext);

  void charged_particle_mover(const amrex::Vector<amrex::MultiFab>& nodeE,
                              const amrex::Vector<amrex::MultiFab>& nodeB,
                              amrex::Real dt, amrex::Real dtNext);

  void neutral_mover(amrex::Real dt);

  void update_position_to_half_stage(const amrex::MultiFab& nodeEMF,
                                     const amrex::MultiFab& nodeBMF,
                                     amrex::Real dt);

  void convert_to_fluid_moments(amrex::Vector<amrex::MultiFab>& momentsMF);

  amrex::IntVect get_ref_ratio(const int iLev) const {
    const amrex::ParGDBBase* gdb = GetParGDB();
    return gdb->refRatio(iLev);
  }

  void set_random_seed(const int iLev, const int i, const int j, const int k,
                       const amrex::IntVect nPPC) {
    amrex::IntVect nCell = Geom(iLev).Domain().size();

    int nRandom = 7;

    int nxcg = nCell[ix_];
    int nycg = nCell[iy_];
    int nzcg = nCell[iz_];
    int iCycle = tc->get_cycle();
    int npcel = nPPC[ix_] * nPPC[iy_] * nPPC[iz_];

    // What if the seed overflow?
    const long seed =
        (speciesID + 3) * nRandom * npcel *
        (nxcg * nycg * nzcg * iCycle + nycg * nzcg * i + nzcg * j + k);

    randNum.set_seed(seed);
  }

  const amrex::Vector<amrex::iMultiFab>& get_cell_status() const {
    return cellStatus;
  }

  const amrex::iMultiFab& get_cell_status(int iLev) const {
    return cellStatus[iLev];
  }

  void update_cell_status(const amrex::Vector<amrex::iMultiFab>& in) {
    for (int iLev = 0; iLev < nLev; iLev++) {
      const int nGst = in[iLev].nGrow();
      distribute_FabArray(cellStatus[iLev], ParticleBoxArray(iLev),
                          ParticleDistributionMap(iLev), 1, nGst, false);

      if (!in[iLev].empty()) {
        amrex::iMultiFab::Copy(cellStatus[iLev], in[iLev], 0, 0,
                               in[iLev].nComp(), nGst);
      }
    }
  }

  ParticleTileType& get_particle_tile(int iLev, const amrex::MFIter& mfi, int i,
                                      int j, int k) {
    return get_particle_tile(iLev, mfi, amrex::IntVect(i, j, k));
  }

  ParticleTileType& get_particle_tile(int iLev, const amrex::MFIter& mfi,
                                      const amrex::IntVect& iv) {
    amrex::Box tileBox;
    const int tileIdx =
        getTileIndex(iv, mfi.validbox(), do_tiling, tile_size, tileBox);
    return GetParticles(iLev)[std::make_pair(mfi.index(), tileIdx)];
  }

  void set_ppc(amrex::IntVect& in) { nPartPerCell = in; };

  void set_bc(BC& bcIn) { bc = bcIn; }

  inline bool is_outside_ba(const ParticleType& p) {
    int iLev = 0;
    amrex::Real loc[3] = { 0, 0, 0 };
    for (int iDim = 0; iDim < 3; iDim++) {
      loc[iDim] = p.pos(iDim);
      if (Geom(iLev).isPeriodic(iDim)) {
        // Fix index/loc for periodic BC.
        while (loc[iDim] > phi[iLev][iDim])
          loc[iDim] -= phi[iLev][iDim] - plo[iLev][iDim];
        while (loc[iDim] < plo[iLev][iDim])
          loc[iDim] += phi[iLev][iDim] - plo[iLev][iDim];
      }
    }

    for (const auto& rb : activeRegions) {
      if (rb.contains(loc))
        return false;
    }

    return true;
  }

  inline bool is_outside_ba(const ParticleType& p,
                            amrex::Array4<int const> const& status,
                            const amrex::IntVect& low,
                            const amrex::IntVect& high) {
    int iLev = 0;
    // Contains ghost cells.
    bool isInsideBox = true;
    int cellIdx[3];
    amrex::Real dShift[3];
    for (int i = 0; i < 3; i++) {
      dShift[i] = (p.pos(i) - plo[iLev][i]) * invDx[iLev][i];
      cellIdx[i] = fastfloor(dShift[i]);
      if (cellIdx[i] > high[i] || cellIdx[i] < low[i]) {
        isInsideBox = false;
        break;
      }
    }

    if (isInsideBox) {
      return bit::is_boundary(status(cellIdx[ix_], cellIdx[iy_], cellIdx[iz_]));
    } else {
      return is_outside_ba(p);
    }
  }

  void label_particles_outside_ba() {
    const int iLev = 0;
    if (NumberOfParticlesAtLevel(iLev, true, true) > 0) {
      for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev);
           pti.isValid(); ++pti) {
        auto& particles = pti.GetArrayOfStructs();
        const amrex::Array4<int const>& status = cellStatus[iLev][pti].array();
        const amrex::Box& bx = cellStatus[iLev][pti].box();
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
  }

  void label_particles_outside_ba_general() {
    const int iLev = 0;
    if (NumberOfParticlesAtLevel(iLev, true, true) > 0) {
      for (ParticlesIter<NStructReal, NStructInt> pti(*this, iLev);
           pti.isValid(); ++pti) {
        auto& particles = pti.GetArrayOfStructs();
        for (auto& p : particles) {
          if (is_outside_ba(p)) {
            p.id() = -1;
            // amrex::Print()<<"particle outside ba = "<<p<<std::endl;
          }
        }
      }
    }
  }

  void split_particles(amrex::Real limit);
  void merge_particles(amrex::Real limit);

  void divE_correct_position(const amrex::MultiFab& phiMF);

  bool is_neutral() const { return charge == 0; };

  int get_speciesID() const { return speciesID; }
  amrex::Real get_charge() const { return charge; }
  amrex::Real get_mass() const { return mass; }

  void set_merge_threshold(amrex::Real in) { mergeThresholdDistance = in; }
  void set_merge_velocity_bin_buffer(amrex::Real in) { velBinBufferSize = in; }
  void fast_merge(bool in) { fastMerge = in; }

  void set_relativistic(const bool& in) { isRelativistic = in; }

  void Write_Paraview(std::string folder = "Particles",
                      std::string particletype = "1") {
    Redistribute();
    std::string command = "python "
                          "../util/AMREX/Tools/Py_util/amrex_particles_to_vtp/"
                          "amrex_binary_particles_to_vtp.py";
    Checkpoint(folder, particletype);
    command = command + " " + folder + " " + particletype;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::system(command.c_str());
    }
    command = "mv";
    command = command + " " + folder + ".vtp" + " " + folder + "_" +
              particletype + ".vtp";
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::system(command.c_str());
    }
    command = "rm -rf";
    command = command + " " + folder;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::system(command.c_str());
    }
  }

  void Write_Binary(std::string folder = "Particles",
                    std::string particletype = "1") {
    Checkpoint(folder, particletype);
  }
};

template <int NStructReal, int NStructInt>
ParticleStaggering Particles<NStructReal, NStructInt>::particlePosition =
    Staggered;

class IOParticles : public Particles<nPicPartReal> {
public:
  IOParticles() = delete;

  IOParticles(Particles<>& other, amrex::AmrCore* amrcore,
              amrex::Real no2outL = 1, amrex::Real no2outV = 1,
              amrex::Real no2OutM = 1,
              amrex::RealBox IORange = amrex::RealBox());
  ~IOParticles() = default;
};
#endif
