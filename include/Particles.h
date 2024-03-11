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

enum class CrossSection { LS = 0, MT };

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

struct IDs {
  int id;
  int supID;
};

struct ParticlesInfo {
  bool fastMerge = false;
  bool mergeLight = false;

  int nPartCombine = 6;
  int nPartNew = 5;
  int nMergeTry = 3;

  amrex::Real mergeThresholdDistance = 0.6;
  amrex::Real velBinBufferSize = 0.125;

  amrex::Real mergeRatioMax = 1.5;
  amrex::Real pLevRatio = 1.2;

  amrex::Real mergePartRatioMax = 10;

  // [amu/cc]
  amrex::Real vacuumIO = 0;
};

template <int NStructReal, int NStructInt>
class ParticlesIter : public amrex::ParIter<NStructReal, NStructInt> {
public:
  using amrex::ParIter<NStructReal, NStructInt>::ParIter;
};

/*
Q: How it the Grid* gridIn variable used inside the particle container?
A: The following codes show how gridIn is passed into the particle container.

  1.
    AmrParticleContainer (AmrCore* amr_core)
        : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt,
Allocator>(amr_core->GetParGDB()){ }

  2.
      ParticleContainer (ParGDBBase* gdb):
        ParticleContainerBase(gdb)

  3.
      ParticleContainerBase (ParGDBBase* gdb)
        :
        m_verbose(0),
        m_gdb(gdb)
    {}

But, what is amr_core->GetParGDB()? It returns a AmrParGDB pointer that is
pointing to AmrCore::m_gdb. AmrParGDB class (object AmrCore::m_gdb) contains a
member variable  AmrCore* m_amrcore, which is a copy of amr_core pointer. So,
the Geometries, DistributionMaps and BoxArrays of Grid can be accessed by the
particle container through m_gdb pointer.

In short, once the grids or distributions maps of Pic or ParticleTracker change,
the particle contains aware of the changes through the m_gdb pointer.
*/

using PicParticle = amrex::Particle<nPicPartReal, nPicPartInt>;
using nPTParticle = amrex::Particle<nPTPartReal, nPTPartInt>;

// Forward declaration.
template <int NStructReal, int NStructInt> class Particles;

using PicParticles = Particles<nPicPartReal, nPicPartInt>;
using PTParticles = Particles<nPTPartReal, nPTPartInt>;

template <int NStructReal, int NStructInt>
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
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::Checkpoint;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::Index;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::ParticlesAt;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::maxLevel;
  using amrex::AmrParticleContainer<NStructReal, NStructInt>::GetParGDB;

  using AoS = amrex::ArrayOfStructs<NStructReal, NStructInt>;

  using PIter = ParticlesIter<NStructReal, NStructInt>;

protected:
  Grid* grid = nullptr;

  FluidInterface* fi = nullptr;
  TimeCtr* tc = nullptr;

  int speciesID;
  RandNum randNum;

  amrex::Real charge;
  amrex::Real mass;

  amrex::Real qom;
  int qomSign;

  amrex::IntVect nPartPerCell;

  amrex::Vector<amrex::RealVect> plo, phi, dx, invDx;
  amrex::Vector<amrex::Real> invVol;

  // ------- Particle resampling begin -------
  amrex::Real mergeThresholdDistance = 0.6;
  amrex::Real velBinBufferSize = 0.125;

  // If fastMerge == false: find the particle pair that is closest to each other
  // in the phase space and try to delete the lighter one.
  // If fastMerge == true: merge nPartCombineMax particles into nPartNew with
  // Lagrange multiplier method.
  bool fastMerge = false;
  int nPartCombine = 6;
  int nPartNew = 5;
  int nMergeTry = 3;
  amrex::Real mergeRatioMax = 1.5;

  bool mergeLight = false;
  amrex::Real mergePartRatioMax = 10;
  // ------- Particle resampling end -------

  amrex::Real pLevRatio = 1.2;

  amrex::Real vacuum = 0;

  bool isRelativistic = false;

  bool isParticleLocationRandom = true;

  BC bc; // boundary condition

  // AMREX uses 40 bits(it is 40! Not a typo. See AMReX_Particle.H) to store
  // p.id(), but it is converted to a 32-bit integer when saving to disk. To
  // avoid the mismatch, FLEKS set the maximum value of p.id() to 2^31-1, and
  // introduce a new integer 'supID' to avoid the overflow of p.id(). See
  // set_ids() below. In short, a FLEKS particle is identified by p.cpu(),
  // p.id() and p.idata(iSupID_).
  int supID = 1;

public:
  static const int iup_ = 0;
  static const int ivp_ = 1;
  static const int iwp_ = 2;
  static const int iqp_ = 3;

  // mu = cos(theta), theta is the pitch angle.
  static const int imu_ = 4;

  TestCase testCase;

  // Index of the integer data.
  static const int iRecordCount_ = 1;

  Particles(Grid* gridIn, FluidInterface* fluidIn, TimeCtr* tcIn,
            const int speciesIDIn, const amrex::Real chargeIn,
            const amrex::Real massIn, const amrex::IntVect& nPartPerCellIn,
            TestCase tcase = RegularSimulation);

  int n_lev() const { return GetParGDB()->finestLevel() + 1; }

  int n_lev_max() const { return maxLevel() + 1; }

  void add_particles_domain();
  void add_particles_cell(const int iLev, const amrex::MFIter& mfi,
                          const amrex::IntVect ijk,
                          const FluidInterface* interface, bool doVacuumLimit,
                          amrex::IntVect ppc = amrex::IntVect(),
                          const Vel tpVel = Vel(), amrex::Real dt = -1);
  void inject_particles_at_boundary();

  void add_particles_source(const FluidInterface* interface,
                            const FluidInterface* const stateOH = nullptr,
                            amrex::Real dt = -1,
                            amrex::IntVect ppc = amrex::IntVect(),
                            const bool doSelectRegion = false,
                            const bool adaptivePPC = false);

  // Copy particles from (ip,jp,kp) to (ig, jg, kg) and shift boundary
  // particle's coordinates accordingly.
  void outflow_bc(const amrex::MFIter& mfi, const amrex::IntVect ijkGst,
                  const amrex::IntVect ijkPhy);

  // 1) Only inject particles ONCE for one ghost cells. This function decides
  // which block injects particles. 2) bx should be a valid box 3) The cell
  // (i,j,k) can NOT be the outmost ghost cell layer!!!!
  bool do_inject_particles_for_this_cell(const amrex::Box& bx,
                                         const amrex::Array4<const int>& status,
                                         const amrex::IntVect ijk,
                                         amrex::IntVect& ijksrc);

  amrex::Real sum_moments(amrex::Vector<amrex::MultiFab>& momentsMF,
                          amrex::Vector<amrex::MultiFab>& nodeBMF,
                          amrex::Real dt);

  std::array<amrex::Real, 5> total_moments(bool localOnly = false);

  void calc_mass_matrix(amrex::UMultiFab<RealMM>& nodeMM, amrex::MultiFab& jHat,
                        amrex::MultiFab& nodeBMF, amrex::MultiFab& u0MF,
                        amrex::Real dt, int iLev, bool solveInCoMov);

  void calc_jhat(amrex::MultiFab& jHat, amrex::MultiFab& nodeBMF,
                 amrex::Real dt);

  // It is real 'thermal velocity'. It is sqrt(sum(q*v2)/sum(q)).
  amrex::Real calc_max_thermal_velocity(amrex::MultiFab& momentsMF);

  void sum_to_center(amrex::MultiFab& netChargeMF,
                     amrex::UMultiFab<RealCMM>& centerMM, bool doNetChargeOnly);

  void charge_exchange(
      amrex::Real dt, FluidInterface* stateOH, FluidInterface* sourcePT2OH,
      SourceInterface* source, bool kineticSource,
      amrex::Vector<std::unique_ptr<PicParticles> >& sourceParts,
      bool doSelectRegion, int nppc);

  void get_ion_fluid(FluidInterface* stateOH, PIter& pti, const int iLev,
                     const int iFluid, const amrex::RealVect xyz,
                     amrex::Real& rhoIon, amrex::Real& cs2Ion,
                     amrex::Real (&uIon)[nDim3]);

  void add_source_particles(std::unique_ptr<PicParticles>& sourcePart,
                            amrex::IntVect ppc, const bool adaptivePPC);

  void mover(const amrex::Vector<amrex::MultiFab>& nodeE,
             const amrex::Vector<amrex::MultiFab>& nodeB,
             const amrex::Vector<amrex::MultiFab>& eBg,
             const amrex::Vector<amrex::MultiFab>& uBg, amrex::Real dt,
             amrex::Real dtNext, bool solveInCoMov);

  void charged_particle_mover(const amrex::Vector<amrex::MultiFab>& nodeE,
                              const amrex::Vector<amrex::MultiFab>& nodeB,
                              const amrex::Vector<amrex::MultiFab>& eBg,
                              const amrex::Vector<amrex::MultiFab>& uBg,
                              amrex::Real dt, amrex::Real dtNext,
                              bool solveInCoMov);

  // Both the input are in the SI unit: m/s
  amrex::Real charge_exchange_dis(amrex::Real* vp, amrex::Real* vh,
                                  amrex::Real* up, amrex::Real vth,
                                  CrossSection cs);

  void sample_charge_exchange(amrex::Real* vp, amrex::Real* vh, amrex::Real* up,
                              amrex::Real vth, CrossSection cs);

  void neutral_mover(amrex::Real dt);

  void update_position_to_half_stage(const amrex::MultiFab& nodeEMF,
                                     const amrex::MultiFab& nodeBMF,
                                     amrex::Real dt);

  void convert_to_fluid_moments(amrex::Vector<amrex::MultiFab>& momentsMF);

  void set_info(ParticlesInfo& pi) {
    fastMerge = pi.fastMerge;
    mergeLight = pi.mergeLight;
    nPartCombine = pi.nPartCombine;
    nPartNew = pi.nPartNew;
    nMergeTry = pi.nMergeTry;
    mergeThresholdDistance = pi.mergeThresholdDistance;
    velBinBufferSize = pi.velBinBufferSize;
    mergeRatioMax = pi.mergeRatioMax;
    pLevRatio = pi.pLevRatio;
    mergePartRatioMax = pi.mergePartRatioMax;
    vacuum = pi.vacuumIO * cProtonMassSI * 1e6 * fi->get_Si2NoRho();
  }

  IDs get_next_ids() {
    constexpr long idMax = 2147483647L;
    long id = ParticleType::NextID();
    if (id > idMax) {
      id = 1;
      ParticleType::NextID(id);
      supID++;
    }

    IDs ids = { static_cast<int>(id), supID };
    return ids;
  }

  // set p.id(), p.cpu() and p.idata(iSupID_)
  void set_ids(ParticleType& p) {
    auto ids = get_next_ids();
    p.id() = ids.id;
    p.idata(iSupID_) = ids.supID;
    p.cpu() = amrex::ParallelDescriptor::MyProc();
  }

  int sup_id() const { return supID; }
  void set_sup_id(int in) { supID = in; }

  amrex::IntVect get_ref_ratio(const int iLev) const {
    const amrex::ParGDBBase* gdb = GetParGDB();
    return gdb->refRatio(iLev);
  }

  // This function distributes particles to proper processors and apply
  // periodic boundary conditions if needed.
  void redistribute_particles() {
    const amrex::ParGDBBase* gdb = GetParGDB();

    if (!gdb->boxArray(0).empty()) {
      // It will crash if there is no active cells.
      Redistribute();
    }
  }

  long calc_random_seed(const int iLev, const amrex::IntVect ijk,
                        const amrex::IntVect nPPC) {
    amrex::IntVect nCell = Geom(iLev).Domain().size();

    int nRandom = 7;

    int nxcg = nCell[ix_];
    int nycg = nCell[iy_];
    int nzcg = 1;
    if (nDim > 2)
      nzcg = nCell[iz_];

    int iCycle = tc->get_cycle();

    int i = ijk[0];
    int j = ijk[1];
    int k = nDim > 2 ? ijk[2] : 0;

    // What if the seed overflow?
    const long seed =
        (speciesID + 3) * nRandom * product(nPPC) *
        (nxcg * nycg * nzcg * iCycle + nycg * nzcg * i + nzcg * j + k);
    return seed;
  }

  long set_random_seed(const int iLev, const amrex::IntVect ijk,
                       const amrex::IntVect nPPC) {
    long seed = calc_random_seed(iLev, ijk, nPPC);
    randNum.set_seed(seed);
    return seed;
  }

  const amrex::iMultiFab& cell_status(int iLev) const {
    return grid->cell_status(iLev);
  }

  const amrex::iMultiFab& node_status(int iLev) const {
    return grid->node_status(iLev);
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

  inline bool is_outside_active_region(const ParticleType& p) {
    int iLev = 0;
    amrex::RealVect loc;
    for (int iDim = 0; iDim < nDim; iDim++) {
      loc[iDim] = p.pos(iDim);
      if (Geom(iLev).isPeriodic(iDim)) {
        // Fix index/loc for periodic BC.
        while (loc[iDim] > phi[iLev][iDim])
          loc[iDim] -= phi[iLev][iDim] - plo[iLev][iDim];
        while (loc[iDim] < plo[iLev][iDim])
          loc[iDim] += phi[iLev][iDim] - plo[iLev][iDim];
      }
    }

    return !grid->is_inside_domain(loc.begin());
  }

  // validBox should NOT include ghost cells.
  inline bool is_outside_active_region(const ParticleType& p, const int iLev,
                                       const amrex::Box& validBox) {
    amrex::IntVect cellIdx = Geom(iLev).CellIndex(p.pos().begin());

    if (validBox.contains(cellIdx))
      return false;

    return is_outside_active_region(p);
  }

  inline void label_particles_outside_active_region() {
    for (int iLev = 0; iLev < n_lev(); iLev++)
      if (NumberOfParticlesAtLevel(iLev, true, true) > 0) {
        for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
          AoS& particles = pti.GetArrayOfStructs();
          if (cell_status(iLev).empty()) {
            for (auto& p : particles) {
              p.id() = -1;
            }
          } else {
            const amrex::Box& validBox = pti.validbox();
            for (auto& p : particles) {
              if (is_outside_active_region(p, iLev, validBox)) {
                p.id() = -1;
              }
            }
          }
        }
      }
  }

  inline void label_particles_outside_active_region_general() {
    for (int iLev = 0; iLev < n_lev(); iLev++)
      if (NumberOfParticlesAtLevel(iLev, true, true) > 0) {
        for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
          AoS& particles = pti.GetArrayOfStructs();
          for (auto& p : particles) {
            if (is_outside_active_region(p)) {
              p.id() = -1;
            }
          }
        }
      }
  }

  void limit_weight(amrex::Real maxRatio, bool seperateVelocity = false);
  void split(amrex::Real limit, bool seperateVelocity = false);
  void split_particles_by_velocity(amrex::Vector<ParticleType*>& plist,
                                   amrex::Vector<ParticleType>& newparticles);
  bool split_by_seperate_velocity(ParticleType& p1, ParticleType& p2,
                                  ParticleType& p3, ParticleType& p4);
  void merge(amrex::Real limit);
  bool merge_particles_fast(int iLev, AoS& particles,
                            amrex::Vector<int>& partIdx,
                            amrex::Vector<int>& idx_I, int nPartCombine,
                            int nPartNew, amrex::Vector<amrex::Real>& x,
                            long seed);

  bool merge_particles_accurate(int iLev, AoS& particles,
                                amrex::Vector<int>& partIdx,
                                amrex::Vector<int>& idx_I, int nPartCombine,
                                int nPartNew, amrex::Vector<amrex::Real>& x,
                                amrex::Real velNorm);

  void divE_correct_position(const amrex::MultiFab& phiMF);

  bool is_neutral() const { return charge == 0; };

  int get_speciesID() const { return speciesID; }
  amrex::Real get_charge() const { return charge; }
  amrex::Real get_mass() const { return mass; }

  void set_relativistic(const bool& in) { isRelativistic = in; }

  void Write_Paraview(std::string folder = "Particles",
                      std::string particletype = "1") {
    redistribute_particles();
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

class IOParticles : public PicParticles {
public:
  IOParticles() = delete;

  IOParticles(PicParticles& other, Grid* gridIn, amrex::Real no2outL = 1,
              amrex::Real no2outV = 1, amrex::Real no2OutM = 1,
              amrex::RealBox IORange = amrex::RealBox());
  ~IOParticles() = default;
};
#endif
