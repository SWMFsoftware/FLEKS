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

enum class PartMode { PIC = 0, Neutral, SEP };

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

struct OHIon {
  amrex::Real rAnalytic = 0;
  amrex::Real rCutoff = 0;
  amrex::Real swRho = 0;
  amrex::Real swT = 0;
  amrex::Real swU = 0;
  bool doGetFromOH = false;
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
  using amrex::AmrParticleContainer<NStructReal,
                                    NStructInt>::CreateGhostParticles;
  using amrex::AmrParticleContainer<NStructReal,
                                    NStructInt>::CreateVirtualParticles;
  using amrex::AmrParticleContainer<NStructReal,
                                    NStructInt>::AddParticlesAtLevel;

  using AoS = amrex::ArrayOfStructs<NStructReal, NStructInt>;

  using PIter = ParticlesIter<NStructReal, NStructInt>;

protected:
  Grid* grid = nullptr;

  FluidInterface* fi = nullptr;
  TimeCtr* tc = nullptr;

  PartMode pMode = PartMode::PIC;

  int speciesID;
  RandNum randNum;
  bool moveParticlesWithConstantVelocity = false;
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

  bool isPPVconstant = false;

  BC bc; // boundary condition

  // AMREX uses 40 bits(it is 40! Not a typo. See AMReX_Particle.H) to store
  // p.id(), but it is converted to a 32-bit integer when saving to disk. To
  // avoid the mismatch, FLEKS set the maximum value of p.id() to 2^31-1, and
  // introduce a new integer 'supID' to avoid the overflow of p.id(). See
  // set_ids() below. In short, a FLEKS particle is identified by p.cpu(),
  // p.id() and p.idata(iSupID_).
  int supID = 1;

  OHIon ionOH;

  bool isFake2D;

public:
  static constexpr int iup_ = 0;
  static constexpr int ivp_ = 1;
  static constexpr int iwp_ = 2;
  static constexpr int iqp_ = 3;

  // mu = cos(theta), theta is the pitch angle.
  static constexpr int imu_ = 4;

  TestCase testCase;

  // Index of the integer data.
  static constexpr int iRecordCount_ = 1;

  Particles(Grid* gridIn, FluidInterface* fluidIn, TimeCtr* tcIn,
            const int speciesIDIn, const amrex::Real chargeIn,
            const amrex::Real massIn, const amrex::IntVect& nPartPerCellIn,
            const PartMode pModeIn, TestCase tcase = RegularSimulation);

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

  amrex::Real sum_moments_new(amrex::Vector<amrex::MultiFab>& momentsMF,
                              amrex::Vector<amrex::MultiFab>& nodeBMF,
                              amrex::Real dt,
                              amrex::Vector<amrex::iMultiFab>& nodestatus);

  std::array<amrex::Real, 5> total_moments(bool localOnly = false);

  void calc_mass_matrix(amrex::UMultiFab<RealMM>& nodeMM, amrex::MultiFab& jHat,
                        amrex::MultiFab& nodeBMF, amrex::MultiFab& u0MF,
                        amrex::Real dt, int iLev, bool solveInCoMov);

  void calc_mass_matrix_amr(
      amrex::UMultiFab<RealMM>& nodeMM,
      amrex::Vector<amrex::Vector<amrex::UMultiFab<RealMM> > >& nmmc,
      amrex::Vector<amrex::UMultiFab<RealMM> >& nmmf, amrex::MultiFab& jHat,
      amrex::Vector<amrex::Vector<amrex::MultiFab> >& jhc,
      amrex::Vector<amrex::MultiFab>& jhf, amrex::MultiFab& nodeBMF,
      amrex::MultiFab& u0MF, amrex::Real dt, int iLev, bool solveInCoMov,
      amrex::Vector<amrex::iMultiFab>& cellstatus);

  void calc_mass_matrix_new(amrex::Vector<amrex::UMultiFab<RealMM> >& nodeMM,
                            amrex::UMultiFab<RealMM>& nmmc,
                            amrex::UMultiFab<RealMM>& nmmf,
                            amrex::Vector<amrex::MultiFab>& jHat,
                            amrex::MultiFab& jhc, amrex::MultiFab& jhf,
                            amrex::Vector<amrex::MultiFab>& nodeBMF,
                            amrex::Vector<amrex::MultiFab>& u0MF,
                            amrex::Real dt, int iLev, bool solveInCoMov,
                            amrex::Vector<amrex::iMultiFab>& nodestatus,
                            amrex::Vector<amrex::iMultiFab>& cellstatus);

  void calc_mass_matrix_new_optimized(
      amrex::Vector<amrex::UMultiFab<RealMM> >& nodeMM,
      amrex::UMultiFab<RealMM>& nmmc, amrex::UMultiFab<RealMM>& nmmf,
      amrex::Vector<amrex::MultiFab>& jHat, amrex::MultiFab& jhc,
      amrex::MultiFab& jhf, amrex::Vector<amrex::MultiFab>& nodeBMF,
      amrex::Vector<amrex::MultiFab>& u0MF, amrex::Real dt, int iLev,
      bool solveInCoMov, amrex::Vector<amrex::iMultiFab>& nodestatus,
      amrex::Vector<amrex::iMultiFab>& cellstatus);

  void calc_jhat(amrex::MultiFab& jHat, amrex::MultiFab& nodeBMF,
                 amrex::Real dt);

  // It is real 'thermal velocity'. It is sqrt(sum(q*v2)/sum(q)).
  amrex::Real calc_max_thermal_velocity(amrex::MultiFab& momentsMF);

  void sum_to_center(amrex::MultiFab& netChargeMF,
                     amrex::UMultiFab<RealCMM>& centerMM, bool doNetChargeOnly,
                     int iLev);

  void sum_to_center_new(amrex::MultiFab& netChargeMF, amrex::MultiFab& jc,
                         amrex::MultiFab& jf,
                         amrex::UMultiFab<RealCMM>& centerMM,
                         bool doNetChargeOnly, int iLev);

  void charge_exchange(
      amrex::Real dt, FluidInterface* stateOH, FluidInterface* sourcePT2OH,
      SourceInterface* source, bool kineticSource,
      amrex::Vector<std::unique_ptr<PicParticles> >& sourceParts,
      bool doSelectRegion, int nppc);

  void get_ion_fluid(FluidInterface* stateOH, PIter& pti, const int iLev,
                     const int iFluid, const amrex::RealVect xyz,
                     amrex::Real& rhoIon, amrex::Real& cs2Ion,
                     amrex::Real (&uIon)[nDim3]);

  // Input:
  // xyz in NO units.
  // Output:
  // rhoIon: amu/m^3
  // cs2Ion: (m/s)^2
  // uIon: m/s
  void get_analytic_ion_fluid(const amrex::RealVect xyz, amrex::Real& rhoIon,
                              amrex::Real& cs2Ion, amrex::Real (&uIon)[nDim3]);

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

  // select particles based on input supid and id
  void select_particle(amrex::Vector<std::array<int, 3> >& selectParticleIn);

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

  PartMode part_mode() const { return pMode; }

  void set_ion_fluid(const OHIon& in) { ionOH = in; }

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

  static inline bool compare_two_parts(const ParticleType& pl,
                                       const ParticleType& pr) {
    // It is non-trivial to compare floating point numbers. If there is
    // significant difference between the two floating point numbers, the
    // comparison is based on the floating point numbers. Otherwise, the
    // comparison is based on the integer numbers (particle ids). However,
    // different number of processors may have different results for id
    // comparison.
    if (fabs(pl.pos(ix_) - pr.pos(ix_)) >
        1e-9 * (fabs(pl.pos(ix_)) + fabs(pr.pos(ix_)))) {
      return pl.pos(ix_) > pr.pos(ix_);
    }

    if (fabs(pl.rdata(iup_) - pr.rdata(iup_)) >
        1e-9 * (fabs(pl.rdata(iup_)) + fabs(pr.rdata(iup_)))) {
      return pl.rdata(iup_) > pr.rdata(iup_);
    }

    return false;
  }

  amrex::Real cosine(ParticleType& p, amrex::Real (&bIn)[nDim3]) {
    amrex::Real u[nDim3];
    amrex::Real b[nDim3];
    for (int i = 0; i < nDim3; ++i) {
      u[i] = p.rdata(iup_ + i);
      b[i] = bIn[i];
    }
    amrex::Real mu = 0;
    for (int i = 0; i < nDim; ++i)
      mu += u[i] * b[i];

    const amrex::Real bNorm = l2_norm(b, nDim3);
    const amrex::Real uNorm = l2_norm(u, nDim3);

    const amrex::Real invB = bNorm > 1e-99 ? 1.0 / bNorm : 0;
    const amrex::Real invU = uNorm > 1e-99 ? 1.0 / uNorm : 0;

    return mu * invB * invU;
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

  /**
   * @brief Sets the IDs for a particle.
   *
   * This function sets the unique ID, supplementary ID, and CPU ID for the
   * given particle.
   *
   * @param p The particle for which the IDs are to be set.
   */
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

  ParticleTileType& get_particle_tile(int iLev, const amrex::MFIter& mfi) {
    return GetParticles(
        iLev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
  }

  void set_ppc(amrex::IntVect& in) { nPartPerCell = in; };

  void set_bc(BC& bcIn) { bc = bcIn; }

  inline bool is_outside_active_region(const ParticleType& p, int iLev) {
    if (iLev > 0) {
      return false;
    }
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

  inline bool is_outside_level(const ParticleType& p, int iLev,
                               amrex::Array4<int const> const& status,
                               const amrex::IntVect cellIdx) {
    bool isOutsideLevel = false;
    if (bit::is_refined(status(cellIdx)) ||
        bit::is_lev_boundary(status(cellIdx))) {
      isOutsideLevel = true;
    }

    return isOutsideLevel;
  }

  /**
   * @brief Checks if a particle is outside the active region at a given level.
   *
   * This function determines whether a particle is outside the active region
   * at a specified level (`iLev`). It takes into account periodic boundary
   * conditions and adjusts the particle's position accordingly.
   *
   * @param p The particle to check.
   * @param iLev The level at which to check the particle's position.
   * @return True if the particle is outside the active region, false otherwise.
   */
  inline bool is_outside_active_region(const ParticleType& p,
                                       amrex::Array4<int const> const& status,
                                       const amrex::IntVect& low,
                                       const amrex::IntVect& high, int iLev) {

    // TODO: It does not work with AMR.
    // Contains ghost cells.
    if (iLev > 0) {
      return false;
    }
    bool isInsideBox = true;
    amrex::IntVect cellIdx;
    amrex::RealVect dShift;
    for (int i = 0; i < nDim; ++i) {
      dShift[i] = (p.pos(i) - plo[iLev][i]) * invDx[iLev][i];
      cellIdx[i] = fastfloor(dShift[i]);
      if (cellIdx[i] > high[i] || cellIdx[i] < low[i]) {
        isInsideBox = false;
        break;
      }
    }

    if (isInsideBox) {
      return bit::is_domain_boundary(status(cellIdx));
    } else {
      return is_outside_active_region(p, iLev);
    }
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
            const amrex::Box& bx = cell_status(iLev)[pti].box();
            const amrex::Array4<int const>& status =
                cell_status(iLev)[pti].array();

            const amrex::IntVect lowCorner = bx.smallEnd();
            const amrex::IntVect highCorner = bx.bigEnd();

            for (auto& p : particles) {
              if (is_outside_active_region(p, status, lowCorner, highCorner,
                                           iLev)) {
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
            if (is_outside_active_region(p, iLev)) {
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

  void divE_correct_position(const amrex::Vector<amrex::MultiFab>& phiMF,
                             int iLev);

  bool is_neutral() const { return charge == 0; };

  int get_speciesID() const { return speciesID; }
  amrex::Real get_charge() const { return charge; }
  amrex::Real get_mass() const { return mass; }

  void set_relativistic(const bool& in) { isRelativistic = in; }

  void Write_Paraview(std::string folder = "Particles",
                      std::string particletype = "1") {
    // redistribute_particles();
    std::string command = "python "
                          "../util/AMREX/Tools/Py_util/amrex_particles_to_vtp/"
                          "amrex_binary_particles_to_vtp.py";
    Checkpoint(folder, particletype);
    command = command + " " + folder + " " + particletype;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      int result = std::system(command.c_str());
      if (result != 0) {
        std::cerr << "Error executing command: " << command << std::endl;
      }
    }
    command = "mv";
    command = command + " " + folder + ".vtp" + " " + folder + "_" +
              particletype + ".vtp";
    if (amrex::ParallelDescriptor::IOProcessor()) {
      int result = std::system(command.c_str());
      if (result != 0) {
        std::cerr << "Error executing command: " << command << std::endl;
      }
    }
    command = "rm -rf";
    command = command + " " + folder;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      int result = std::system(command.c_str());
      if (result != 0) {
        std::cerr << "Error executing command: " << command << std::endl;
      }
    }
  }

  void Write_Binary(std::string folder = "Particles",
                    std::string particletype = "1") {
    Checkpoint(folder, particletype);
  }

  void Generate_GhostParticles(int iLev, int nGhost) {
    ParticleTileType ptile;
    CreateGhostParticles(iLev - 1, nGhost, ptile);
    AddParticlesAtLevel(ptile, iLev, nGhost);
  }

  void Generate_VirtualParticles(int iLev) {
    ParticleTileType ptile;
    CreateVirtualParticles(iLev + 1, ptile);
    AddParticlesAtLevel(ptile, iLev);
  }

  void Exchange_VirtualParticles(int iLev) {
    ParticleTileType ptile;
    ParticleTileType ptile2;
    CreateVirtualParticles(iLev + 1, ptile);
    CreateGhostParticles(iLev, 1, ptile2);
    AddParticlesAtLevel(ptile, iLev);
    AddParticlesAtLevel(ptile2, iLev + 1, 1);
  }

  void delete_particles_from_refined_region(int iLev) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      AoS& particles = pti.GetArrayOfStructs();
      const auto& status = cell_status(iLev)[pti].array();
      for (auto& p : particles) {
        amrex::IntVect loIdx;
        amrex::RealVect dShift;
        find_cell_index_exp(p.pos(), Geom(iLev).ProbLo(),
                            Geom(iLev).InvCellSize(), loIdx, dShift);
        if (bit::is_refined(status(loIdx))) {
          p.id() = -1;
        }
      }
    }
  }

  void delete_particles_from_ghost_cells(int iLev) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      AoS& particles = pti.GetArrayOfStructs();
      const auto& status = cell_status(iLev)[pti].array();
      for (auto& p : particles) {
        amrex::IntVect loIdx;
        amrex::RealVect dShift;
        find_cell_index_exp(p.pos(), Geom(iLev).ProbLo(),
                            Geom(iLev).InvCellSize(), loIdx, dShift);
        if (bit::is_lev_boundary(status(loIdx))) {
          p.id() = -1;
        }
      }
    }
  }

  void shape_fix_DisplaceEqually4() {

    for (int iLev = 0; iLev < n_lev() - 1; iLev++) {
      amrex::Real dx = Geom(iLev).CellSize(iLev);
      amrex::Real ratio = 0.0;
      amrex::Real disp = dx * 0.25 * sqrt(2.0);
      amrex::Real theta = 45.0;
      amrex::Real PI = 3.14159265358979323846;
      theta = theta * PI / 180.0;
      for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
        // auto& pTile1 = get_particle_tile(iLev, pti);
        auto& pTile2 = get_particle_tile(iLev + 1, pti);
        AoS& particles = pti.GetArrayOfStructs();
        for (auto& p : particles) {
          if (p.id() < 0 || abs(p.rdata(iqp_)) < ((1 - ratio)) * (3.5e-4) / 1.0)
            continue;
          const amrex::Real xp = p.pos(ix_);
          const amrex::Real yp = p.pos(iy_);
          // const amrex::Real zp = nDim > 2 ? p.pos(iz_) : 0;
          if (abs(xp) <= 20.0 && abs(yp) <= 20.0) {
            amrex::Real up = p.rdata(iup_);
            amrex::Real vp = p.rdata(ivp_);
            amrex::Real wp = p.rdata(iwp_);
            amrex::Real qp = p.rdata(iqp_);
            p.id() = -1;
            amrex::Vector<ParticleType> newparticles;
            ParticleType pnew1;
            ParticleType pnew2;
            ParticleType pnew3;
            ParticleType pnew4;
            ParticleType pnew5;
            set_ids(pnew1);
            set_ids(pnew2);
            set_ids(pnew3);
            set_ids(pnew4);
            set_ids(pnew5);
            pnew1.pos(ix_) = xp + (disp * cos(theta));
            pnew1.pos(iy_) = yp + (disp * sin(theta));
            pnew1.pos(iz_) = 0.0;
            pnew1.rdata(iup_) = up;
            pnew1.rdata(ivp_) = vp;
            pnew1.rdata(iwp_) = wp;
            pnew1.rdata(iqp_) = (1.0 - ratio) * qp / 4.0;
            pnew2.pos(ix_) = xp - (disp * cos(theta));
            pnew2.pos(iy_) = yp - (disp * sin(theta));
            pnew2.pos(iz_) = 0.0;
            pnew2.rdata(iup_) = up;
            pnew2.rdata(ivp_) = vp;
            pnew2.rdata(iwp_) = wp;
            pnew2.rdata(iqp_) = (1.0 - ratio) * qp / 4.0;
            pnew3.pos(ix_) = xp + (disp * cos(theta + PI / 2));
            pnew3.pos(iy_) = yp + (disp * sin(theta + PI / 2));
            pnew3.pos(iz_) = 0.0;
            pnew3.rdata(iup_) = up;
            pnew3.rdata(ivp_) = vp;
            pnew3.rdata(iwp_) = wp;
            pnew3.rdata(iqp_) = (1.0 - ratio) * qp / 4.0;
            pnew4.pos(ix_) = xp - (disp * cos(theta + PI / 2));
            pnew4.pos(iy_) = yp - (disp * sin(theta + PI / 2));
            pnew4.pos(iz_) = 0.0;
            pnew4.rdata(iup_) = up;
            pnew4.rdata(ivp_) = vp;
            pnew4.rdata(iwp_) = wp;
            pnew4.rdata(iqp_) = (1.0 - ratio) * qp / 4.0;
            pnew5.pos(ix_) = xp;
            pnew5.pos(iy_) = yp;
            pnew5.pos(iz_) = 0.0;
            pnew5.rdata(iup_) = up;
            pnew5.rdata(ivp_) = vp;
            pnew5.rdata(iwp_) = wp;
            pnew5.rdata(iqp_) = ratio * qp;
            newparticles.push_back(pnew1);
            newparticles.push_back(pnew2);
            newparticles.push_back(pnew3);
            newparticles.push_back(pnew4);
            // newparticles.push_back(pnew5);

            for (auto& p : newparticles) {
              pTile2.push_back(p);
            }

            // pTile.push_back(pnew1);
            // pTile.push_back(pnew2);
            // pTile.push_back(pnew3);
            // pTile.push_back(pnew4);
          }
        }
      }
    }
  }
};

class IOParticles : public PicParticles {
public:
  IOParticles() = delete;

  IOParticles(PicParticles& other, Grid* gridIn, amrex::Real no2outL,
              amrex::Real no2outV, amrex::Real no2OutM, amrex::RealBox IORange);
  ~IOParticles() = default;
};
#endif
