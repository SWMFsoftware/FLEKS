#ifndef _PIC_H_
#define _PIC_H_

#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "Array1D.h"
#include "Bit.h"
#include "Constants.h"
#include "DomainParameters.h"
#include "FleksDistributionMap.h"
#include "FluidInterface.h"
#include "GridAccess.h"
#include "InitialCondition.h"
#include "IntrinsicBField.h"
#include "LinearSolver.h"
#include "MeshChange.h"
#include "OHInterface.h"
#include "Particles.h"
#include "ReadParam.h"
#include "Regions.h"
#include "Shape.h"
#include "SourceInterface.h"
#include "TimeCtr.h"
#include "WaveBC.h"

class ParticleTracker;
class Pic;

class FieldSolver {
public:
  amrex::Real theta;
  amrex::Real coefDiff;
  bool useLaggedLimiter;
  FieldSolverMode mode;
  FieldSolver() {
    theta = 0.51;
    coefDiff = 0.1;
    useLaggedLimiter = false;
    mode = FieldSolverMode::GMRES;
  }
};

typedef amrex::Real (Pic::*GETVALUE)(amrex::MFIter &mfi, amrex::IntVect ijk,
                                     int iVar, const int iLev);

typedef void (Pic::*PicWriteAmrex)(const std::string &filename,
                                   const std::string varName);

struct NodeMMCommData {
  struct LocTagEntry {
    int srcIndex;
    int dstIndex;
    amrex::Box sbox;
    amrex::Box dbox;
    std::size_t bufOffset;
  };

  struct RemoteTagEntry {
    int boxIndex;
    amrex::Box box;
    std::size_t bufOffset;
  };

  struct PeerComm {
    int rank;
    std::vector<RemoteTagEntry> tags;
    std::vector<RealMM> buf;
    std::size_t totalPts = 0;
  };

  bool is_initialized = false;
  amrex::FabArrayBase::BDKey bdkey;

  std::vector<LocTagEntry> loc_tags;
  std::vector<RealMM> local_buf;
  std::vector<std::vector<int> > box_to_loc_tags;

  std::vector<PeerComm> sends;
  std::vector<PeerComm> recvs;

#ifdef BL_USE_MPI
  std::vector<MPI_Request> recv_reqs;
  std::vector<MPI_Request> send_reqs;
#endif
};

// The grid is defined in DomainGrid. This class contains the data on the grid.
class Pic : public GridAccess {
  friend PlotWriter;
  friend ParticleTracker;
  // private variables
private:
  bool usePIC = true;
  bool solveEM = true;
  bool initEM = true;

  // ---- Hybrid PIC (kinetic ions + fluid electrons) solver ----
  bool useHybridPIC = false;
  // Resistive term eta * J. SI input [m^2/s], converted to code units.
  amrex::Real etaResistivitySI = 0.0;
  amrex::Real etaResistivity = 0.0;
  // Electron pressure gradient. Input [eV], converted to code units.
  amrex::Real electronTemperatureEV = 0.0;
  amrex::Real electronTemperature = 0.0;
  // Polytropic index for the adiabatic electron pressure closure.
  amrex::Real electronGamma = 1.0;
  // Reference charge density. Input [amu/cc], converted to code units.
  amrex::Real electronDensity0In = 1.0;
  amrex::Real electronDensity0 = 0.0;
  // Number of sub-steps for the B-field update within one dt.
  int nBSubcycle = 1;
  // Hall term in the generalized Ohm's law.
  bool useHallTerm = true;

  // Hyper-resistivity (fourth-order) term in the Ohm's law:
  //   E -= eta_h * nabla^2 J = -(eta_h/4*pi) * nabla x (nabla^2 B).
  // etaHyperMode selects how etaHyperSI is interpreted:
  //   "si"   -> direct physical value [m^4/s], converted to code units;
  //   "grid" -> CFL-scaled eta_h = C_h * dx^4 / dt_sub (dimless C_h).
  // etaHyperLev[iLev] (code units, 0 disables) is the value actually applied.
  amrex::Real etaHyperSI = 0.0;
  std::string etaHyperMode = "si";
  amrex::Real etaHyperCh = 0.01;
  amrex::Vector<amrex::Real> etaHyperLev;

  // Regional resistivity and hyper-resistivity
  struct RegionalResistivityConfig {
    std::string regionStr;
    amrex::Real etaSI = 0.0;
    amrex::Real etaCode = 0.0;
  };

  struct RegionalHyperResistivityConfig {
    std::string regionStr;
    amrex::Real etaSI = 0.0;
    std::string mode = "grid";
    amrex::Real ch = 0.0;
    amrex::Vector<amrex::Real> etaLev;
  };

  amrex::Vector<RegionalResistivityConfig> regionalResistivityConfigs;
  amrex::Vector<RegionalHyperResistivityConfig> regionalHyperResistivityConfigs;
  bool hasRegionalResistivity_ = false;
  bool hasRegionalHyper_ = false;
  amrex::Vector<std::shared_ptr<Shape> > regionShapes;
  amrex::Vector<amrex::MultiFab> nodeEtaRegional;
  amrex::Vector<amrex::MultiFab> nodeEtaHyperRegional;

  // Minimum charge density in the Hall and electron pressure gradient term.
  // <= 0 means auto: 1e-6 * electronDensity0.
  amrex::Real rhoMinOhm = 0.0;

  // ---- Evolved electron pressure (see #ELECTRONPRESSURE) ----
  bool useElectronPressureEq = false;
  // Spitzer kappa0 [W/(m K^(7/2))]; 9.2e-12 is for coulombLog = 20.
  amrex::Real heatCondKappa0SI = 9.2e-12;
  amrex::Real heatCondKappa0 = 0.0; // code units (convert_electron_heat_cond)
  amrex::Real coulombLog = 20.0;
  // kappa_hat = kappa [f b b + (1-f) I] with f = fieldAlignedFraction.
  bool fieldAlignedConduction = true;
  amrex::Real fieldAlignedFraction = 1.0;
  // |B| [T] under which the direction is treated as noise: use isotropic.
  amrex::Real fieldAlignedBMinSI = 1.0e-15;
  amrex::Real fieldAlignedBMin = 0.0;
  amrex::Real heatFluxLimiter = 0.0; // free-streaming limiter; 0 disables it
  amrex::Real peMin = 0.0;           // floor on the evolved Pe
  // Numerical scheme for Pe advection and compression (see #ELECTRONADVECTION)
  std::string peAdvectionLimiter = "vanleer";
  std::string peCompressionScheme = "exponential";
  int peLimiterType = 2; // 0: upwind1, 1: minmod, 2: vanleer, 3: mc
  bool peCompressionExp = true;
  // Solving strategy for Pe heat conduction (see #ELECTRONCONDUCTION)
  std::string heatCondMethod = "point-implicit";
  int nCondIter = 1;
  int nCondSubcycleMax = 100;
  // Electron-ion collisional heat exchange (#ELECTRONCOLLISION)
  bool useHeatExchange = false;
  amrex::Real collisionFactor = 1.0;
  amrex::Real collisionCoefEi = 0.0; // code units
  // Add the ambipolar E to the RK stages; needed once Pe is not polytropic.
  bool ambipolarInStages = true;
  bool peStateInitialized_ = false; // centerPeState holds a valid field

  bool useExplicitPIC = false;
  bool projectDownEmFields = true;
  bool skipMassMatrix = false;
  bool reportParticleQuality = false;

  PartMode pMode = PartMode::PIC;

  FluidInterface *fi = nullptr;
  FluidInterface *stateOH = nullptr;
  FluidInterface *sourcePT2OH = nullptr;
  SourceInterface *source = nullptr;

  // Wave-injection boundary-condition manager (#WAVEBC).
  WaveBoundaryManager waveBC;
  TimeCtr *tc = nullptr;

  const DomainParameters &domainParameters;

  amrex::Vector<amrex::MultiFab> nodeE;
  amrex::Vector<amrex::MultiFab> nodeEth;
  amrex::Vector<amrex::MultiFab> nodeB;
  // div(B) diagnostics, both cell-centred: divB is div of the nodal field,
  // the quantity the cleaning acts on, centerDivB is div of the cell-centred
  // field the B update advances. Allocated in distribute_arrays (full PIC)
  // or on demand by ensure_divB() / ensure_centerDivB() (hybrid PIC).
  amrex::Vector<amrex::MultiFab> divB;
  amrex::Vector<amrex::MultiFab> centerDivB;
  amrex::Vector<amrex::MultiFab> centerB;
  // Hybrid hyper-resistivity scratch fields.
  amrex::Vector<amrex::MultiFab> centerLapB; // nabla^2 B  (stage A)
  amrex::Vector<amrex::MultiFab> nodeHyperE; // nabla x (nabla^2 B) (stage B)
  // Hybrid RK method shared intermediate solver scratch.
  amrex::Vector<amrex::MultiFab> centerBstage;
  // kStage[iLev][0..3]: the stage curls curl(E_stage) for the level-iLev.
  amrex::Vector<amrex::Vector<amrex::MultiFab> > kStage;

  // RK persistent scratch.
  amrex::Vector<amrex::MultiFab> centerBstart;
  amrex::Vector<amrex::MultiFab> centerBstar; // time-centered state used by E

  amrex::Vector<amrex::MultiFab> dBdt;
  amrex::Vector<amrex::MultiFab> particleQuality;

  // Hyperbolic cleaning
  bool useHyperbolicCleaning = false;
  amrex::Vector<amrex::MultiFab> hypPhi;
  amrex::Real hypDecay = 0.1;
  // Fill divB even when no cleaning is requested, so 'divB' can be plotted.
  bool alwaysComputeDivB = false;
  // With the cleaning on, div(B) comes for free as its by-product.
  bool need_divB() const { return alwaysComputeDivB || useHyperbolicCleaning; }

  // Background velocity and electric field.
  amrex::Vector<amrex::MultiFab> uBg;
  amrex::Vector<amrex::MultiFab> eBg;

  // Mach number: u/v_th
  amrex::Vector<amrex::MultiFab> mMach;

  amrex::Vector<NodeMMFab> nodeMM;
  amrex::Vector<NodeMMCommData> nodeMM_comm_data;

  // ------divE correction--------------
  // Old @ t=t_{n-1/2}; N @ t=t_n; New @ t=t_{n+1/2}
  amrex::Vector<amrex::MultiFab> centerNetChargeOld, centerNetChargeN,
      centerNetChargeNew;
  amrex::Vector<amrex::MultiFab> centerDivE, centerPhi;
  amrex::Vector<amrex::MultiFab> divEInMF, divEOutMF;
  amrex::Vector<CenterMMFab> centerMM;
  const amrex::Real rhoTheta = 0.51;
  //--------------------------------------

  LinearSolver eSolver;
  LinearSolver divESolver;

  // Persistent scratch MultiFabs for implicit field solver (full-PIC)
  amrex::Vector<amrex::MultiFab> solverVecMF;
  amrex::Vector<amrex::MultiFab> solverMatvecMF;
  amrex::Vector<amrex::MultiFab> solverTempNode3;
  amrex::Vector<amrex::MultiFab> solverCenterLapMF;
  amrex::Vector<amrex::MultiFab> solverTempCenter3;
  amrex::Vector<amrex::MultiFab> solverTempCenter1;
  amrex::Vector<amrex::MultiFab> solverRhsNode1;
  amrex::Vector<amrex::MultiFab> solverRhsNode2;
  amrex::Vector<amrex::MultiFab> centerDB;
  amrex::Vector<amrex::MultiFab> smoothScratchMF;
  amrex::Vector<amrex::MultiFab> projectScratchMF;

  int nSpecies;
  int iTot;

  // Hybrid species IDs
  // iElectron_ (-1 if none); kineticSpecies_ = non-electron species.
  int iElectron_ = -1;
  std::vector<int> kineticSpecies_;
  amrex::Vector<amrex::Vector<amrex::MultiFab> > nodePlasma;
  // Ion moments at J^{n-1/2}; interpolated with current nodePlasma
  // by hstep inside assemble_ohm_E.
  amrex::Vector<amrex::Vector<amrex::MultiFab> > nodePlasmaPrev;
  // ---- Staggered hybrid solver fields ----
  amrex::Vector<amrex::MultiFab> nodeEstage; // E at a stage B (nodal)
  amrex::Vector<amrex::MultiFab> nodeJ;      // total current J = curl(B)/(4*pi)
                                             // (nodal)
  amrex::Vector<amrex::MultiFab> nodeBstage; // B interpolated to nodes at RK
                                             // stages
  // Per-step scratch reused across stages: the advected Pe, then div(q), then
  // the scalar ion pressure. Also the Pe the Ohm's law differentiates when
  // #ELECTRONPRESSURE is off, so never read it as the electron pressure --
  // diagnostics must read centerPeState.
  amrex::Vector<amrex::MultiFab> centerPe;
  amrex::Vector<amrex::MultiFab> nodeEambi;   // ambipolar electric field
                                              // -grad(Pe)/(e*ne) at nodes
  amrex::Vector<amrex::MultiFab> nodeRhoTemp; // scratch for time-interpolated
                                              // density
  // ---- Evolved electron pressure (#ELECTRONPRESSURE) ----
  // centerPeState is the evolved field (written to / read from the restart
  // files and is copied on a regrid); everything else is per-step scratch.
  // All are allocated with the full nGst (= 2) ghost layers, which the
  // advection MUSCL stencil needs: it reads the state two cells away in every
  // direction.
  amrex::Vector<amrex::MultiFab> centerPeState; // Pe at cell centers (state)
  amrex::Vector<amrex::MultiFab> centerPeRho;   // n_e; read pointwise only, so
                                                // it carries no ghost BC
  amrex::Vector<amrex::MultiFab> centerPeTe;    // Te; conductivity scratch, NOT
                                                // the reported Te
  amrex::Vector<amrex::MultiFab> nodePeVec;     // u_e, then grad(Te), then q
  amrex::Vector<amrex::MultiFab> nodePeRho; // nodal ion charge density, then
                                            // the collision stage's Pi
  amrex::Vector<amrex::MultiFab> nodePeAux; // (0) Te, (1..3) kappa_hat_dd
  amrex::Vector<amrex::Real> plasmaEnergy;

  bool isMomentsUpdated = false;

  amrex::Vector<amrex::MultiFab> jHat;

  amrex::Vector<std::unique_ptr<PicParticles> > parts;
  amrex::Vector<std::unique_ptr<PicParticles> > sourceParts;

  amrex::Real qomEl = -100;

  // Particle Per Cell (PPC) of source particles.
  amrex::IntVect nSourcePPC = { AMREX_D_DECL(0, 0, 0) };
  bool adaptiveSourcePPC = false;
  bool kineticSource = false;
  amrex::Real maxExchangeRatio = 0;
  amrex::Real maxExchangeRatioLimit = 1;

  FieldSolver fsolver;

  bool doCorrectDivE = true;
  int nDivECorrection = 3;

  bool doReSampling = true;
  amrex::Real reSamplingLowLimit = 0.8;
  amrex::Real reSamplingHighLimit = 1.5;
  amrex::Real maxWeightRatio = 1.0;

  bool solveFieldInCoMov = false;
  int nSmoothBackGroundU = 0;

  bool useUpwindE = false;
  amrex::Real limiterThetaE = 0;
  amrex::Real cMaxE = -1;
  bool useUpwindB = false;
  amrex::Real limiterThetaB = 0;
  // Override upwind velocity in correct_B(). 0 = use plasma background
  // velocity.
  amrex::Real fixedUpwindVel = 0.0;

  // Override uMax for CFL estimate. < 0 = estimate from particle thermal
  // velocity.
  amrex::Real fixedUMax = -1.0;

  bool doSmoothJ = false;
  int nSmoothJ = 0;
  amrex::Real coefSmoothJ = 0.5;

  // Smoothing of ion moments before the generalized Ohm's law.
  bool doSmoothMoments = false;
  int nSmoothMoments = 0;
  amrex::Real coefSmoothMoments = 0.5;

  std::string fieldIntegrator = "rk4"; // B integrator
  bool useRK4 = false;

  // Guard: true on the first hybrid step before nodePlasmaPrev is seeded.
  bool isFirstHybridStep = true;

  bool doSmoothE = false;
  int nSmoothE = 0;

  // Plug-in initial condition via #TESTCASE registry.
  std::unique_ptr<InitialCondition> ic_;

  ParticlesInfo pInfo;

  BoxBC<FieldBC::Type> bcField;
  bool hasConductingBC_ = false;
  bool hasAbsorbBC_ = false;
  bool hasInflowBC_ = false;

  // Field condition on the surface of the inner body (#BODYBOUNDARY); the
  // particle condition lives in Grid.
  BodyFieldBC::Type bodyFieldBC = BodyFieldBC::linetied;
  bool bodyBoundarySet_ = false;

  bool is_body_linetied() const {
    return useBody && bodyFieldBC == BodyFieldBC::linetied;
  }
  bool is_body_conducting() const {
    return useBody && bodyFieldBC == BodyFieldBC::conducting;
  }
  bool is_body_insulating() const {
    return useBody && bodyFieldBC == BodyFieldBC::insulating;
  }
  // Static intrinsic field of the planet (#DIPOLE / #CRUSTALFIELD), on the
  // same grids as the evolved field. Filled at init and after every regrid,
  // never advanced by Faraday's law.
  amrex::Vector<amrex::MultiFab> nodeB0;
  amrex::Vector<amrex::MultiFab> centerB0;
  // Scratch for the cell-centered total field B1 + B0 used by the convective
  // and Hall terms of the hybrid Ohm's law. The current is computed from B1
  // alone: the intrinsic field is current-free, and the discrete curl of
  // B1 + B0 would feed its truncation error into J.
  amrex::Vector<amrex::MultiFab> centerBtotal;

  // The interior of the body (the body minus its one-cell-thick surface layer)
  // is a cavity: E = 0 and B frozen at its initial value. In full-PIC, an
  // 'insulating' body lets EM waves propagate inside via the wave equation;
  // in hybrid PIC, vacuum cavities have no wave propagation and stay frozen.
  bool is_body_interior_frozen() const {
    return useBody &&
           (!useHybridPIC ? (bodyFieldBC != BodyFieldBC::insulating) : true);
  }

  void update_bc_flags() {
    hasConductingBC_ = bcField.has(FieldBC::conducting);
    hasAbsorbBC_ = bcField.has(FieldBC::absorb);
    hasInflowBC_ = bcField.has(FieldBC::inflow) || bcField.has(FieldBC::fixed);
  }

  // De-duplicated boundary-condition warnings
  std::set<std::string> bcWarnings_;

  // Record a boundary-condition warning; empty messages are ignored.
  void add_bc_warning(const std::string &msg) {
    if (!msg.empty())
      bcWarnings_.insert(msg);
  }

  bool fieldBCSet_ = false;

  // Characteristic speed for the absorbing BC; 0 = auto (light speed).
  amrex::Real absorbCharSpeed = 0.0;

  // Inflow boundary upstream state set by #INFLOW (in code units).
  bool inflowDefined_ = false;
  amrex::Real inflowRho_ = 0.0;
  amrex::Real inflowUx_ = 0.0, inflowUy_ = 0.0, inflowUz_ = 0.0;
  amrex::Real inflowT_ = 0.0;

  // Static intrinsic magnetic field B0 of the planet (#DIPOLE /
  // #CRUSTALFIELD). It never evolves: it is added to the evolved field B1
  // wherever a *total* magnetic field is needed.
  std::unique_ptr<IntrinsicBField> intrinsicB_;
  // True once either model is enabled and its coefficients are in code units.
  bool use_intrinsic_B() const {
    return intrinsicB_ != nullptr && intrinsicB_->is_active();
  }

  // select particle params
  bool doSelectParticle = false;
  std::string selectParticleInputFile;

  bool doReport = false;
  amrex::Real maxCFL = 0.0;
  int dnMemory = -1;

  std::string logFile;
  std::ofstream picLogStream;

protected:
  // The shared Grid reference and the mesh queries live in GridAccess.
  std::string tag = "pic";
  std::string gridName;
  std::string printPrefix;
  int gridID;
  int nGst;
  const bool &isFake2D;
  const bool &isGridEmpty;
  const int &finest_level;
  const amrex::Vector<amrex::BoxArray> &cGrids;
  const amrex::Vector<amrex::BoxArray> &nGrids;
  const amrex::Vector<amrex::IntVect> &ref_ratio;
  amrex::Vector<amrex::iMultiFab> &cellStatus;
  amrex::Vector<amrex::iMultiFab> &nodeStatus;
  amrex::Vector<amrex::MultiFab> &cellCost;
  const amrex::BoxArray &activeRegion;
  const bool &useBody;
  const amrex::Real &bodyRadius;
  const amrex::Real *bodyCenter;

  amrex::Vector<amrex::iMultiFab> targetPPC;
  bool isTargetPPCDefined = false;
  bool doNeedFillNewCell = false;

  // public methods
public:
  // Shared mesh queries (n_lev, Geom, DistributionMap, cell_status,
  // get_base_grid, lev_string, get_finest_lev, ...) are inherited from
  // GridAccess.  The queries below are PIC-specific.
  bool is_inside_domain(const amrex::Real *loc) const {
    return grid.is_inside_domain(loc);
  }
  const amrex::Vector<amrex::RealBox> &domain_range() const {
    return grid.domain_range();
  }
  bool use_body() const { return grid.use_body(); }
  amrex::Real get_body_radius() const { return grid.get_body_radius(); }
  const amrex::Real *get_body_center() const { return grid.get_body_center(); }
  int get_dim() const { return grid.get_dim(); }
  const amrex::Vector<amrex::MultiFab> &get_cost() const {
    return grid.get_cost();
  }
  const amrex::iMultiFab &target_PPC(int iLev) const { return targetPPC[iLev]; }
  amrex::Real get_cell_volume(int iLev) const {
    return grid.get_cell_volume(iLev);
  }
  bool is_inside_body(const amrex::Real *loc) const {
    return grid.is_inside_body(loc);
  }
  void set_body(const amrex::Real *center, const amrex::Real radius) {
    grid.set_body(center, radius);
  }
  ParticleBC::Type &body_particle_bc() { return grid.bodyParticleBC; }
  ParticleBC::Type body_particle_bc() const { return grid.bodyParticleBC; }

  Pic(Grid &gridIn, FluidInterface *fluidIn, TimeCtr *tcIn, int id,
      const DomainParameters &parameters)
      : GridAccess(gridIn),
        fi(fluidIn),
        tc(tcIn),
        domainParameters(parameters),
        gridID(id),
        nGst(gridIn.get_n_ghost()),
        isFake2D(gridIn.is_fake_2d_ref()),
        isGridEmpty(gridIn.is_grid_empty_ref()),
        finest_level(gridIn.get_finest_level_ref()),
        cGrids(gridIn.box_arrays()),
        nGrids(gridIn.node_box_arrays()),
        ref_ratio(gridIn.ref_ratios()),
        cellStatus(gridIn.cell_status()),
        nodeStatus(gridIn.node_status()),
        cellCost(gridIn.cell_cost()),
        activeRegion(gridIn.active_region_ref()),
        useBody(gridIn.use_body_ref()),
        bodyRadius(gridIn.get_body_radius_ref()),
        bodyCenter(gridIn.get_body_center()) {
    gridName = std::string("FLEKS") + std::to_string(gridID);
    printPrefix = gridName + " pic: ";
    eSolver.set_tol(1e-6);
    eSolver.set_nIter(200);

    divESolver.set_tol(0.01);
    divESolver.set_nIter(20);

    //-----------------------------------------------------
    centerB.resize(n_lev_max());
    nodeB.resize(n_lev_max());
    nodeB0.resize(n_lev_max());
    centerB0.resize(n_lev_max());
    centerBtotal.resize(n_lev_max());
    dBdt.resize(n_lev_max());
    nodeE.resize(n_lev_max());
    nodeEth.resize(n_lev_max());
    divB.resize(n_lev_max());
    hypPhi.resize(n_lev_max());
    centerLapB.resize(n_lev_max());
    nodeHyperE.resize(n_lev_max());
    centerBstage.resize(n_lev_max());
    nodeEstage.resize(n_lev_max());
    nodeJ.resize(n_lev_max());
    nodeBstage.resize(n_lev_max());
    centerPe.resize(n_lev_max());
    nodeEambi.resize(n_lev_max());
    nodeRhoTemp.resize(n_lev_max());
    centerBstart.resize(n_lev_max());
    centerBstar.resize(n_lev_max());
    kStage.resize(n_lev_max());
    for (int iL = 0; iL < n_lev_max(); ++iL)
      kStage[iL].resize(4);
    etaHyperLev.resize(n_lev_max(), 0.0);
    targetPPC.resize(n_lev_max());
    if (reportParticleQuality) {
      particleQuality.resize(n_lev_max());
    }
    eBg.resize(n_lev_max());
    uBg.resize(n_lev_max());

    mMach.resize(n_lev_max());

    centerNetChargeOld.resize(n_lev_max());
    centerNetChargeN.resize(n_lev_max());
    centerNetChargeNew.resize(n_lev_max());

    centerDivE.resize(n_lev_max());
    centerPhi.resize(n_lev_max());
    divEInMF.resize(n_lev_max());
    divEOutMF.resize(n_lev_max());

    nodeMM.resize(n_lev_max());
    centerMM.resize(n_lev_max());

    jHat.resize(n_lev_max());

    solverVecMF.resize(n_lev_max());
    solverMatvecMF.resize(n_lev_max());
    solverTempNode3.resize(n_lev_max());
    solverCenterLapMF.resize(n_lev_max());
    solverTempCenter3.resize(n_lev_max());
    solverTempCenter1.resize(n_lev_max());
    solverRhsNode1.resize(n_lev_max());
    solverRhsNode2.resize(n_lev_max());
    centerDB.resize(n_lev_max());
    smoothScratchMF.resize(n_lev_max());
    projectScratchMF.resize(n_lev_max());

#ifdef _PT_COMPONENT_
    kineticSource = true;
    initEM = false;
    solveEM = false;

    doCorrectDivE = false;

    pMode = PartMode::Neutral;
#endif
  };
  ~Pic() {
    if (picLogStream.is_open()) {
      picLogStream.close();
    }
  };

  void free_memory();

  void update(bool doReportIn = false);

  PicParticles *get_particle_pointer(int i) { return parts[i].get(); }
  // TODO: no longer needed if we fix the output variable ordering.
  bool get_useHybridPIC() const { return useHybridPIC; }

  // Returns the cell-centered coarse-to-fine interpolater to use.  For the
  // hybrid solver (useHybridPIC), uses CellConservativeLinear (lincc_interp,
  // 2nd-order conservative with slope limiting) for higher accuracy at
  // coarse-fine interfaces.  For the full-PIC solver, keeps CellBilinear
  // (cell_bilinear_interp).
  amrex::Interpolater *get_cell_interp() const {
    return useHybridPIC
               ? static_cast<amrex::Interpolater *>(&amrex::lincc_interp)
               : static_cast<amrex::Interpolater *>(
                     &amrex::cell_bilinear_interp);
  }

  void set_stateOH(OHInterface *in) { stateOH = in; }
  void set_sourceOH(OHInterface *in) { sourcePT2OH = in; }
  void set_fluid_source(SourceInterface *in) { source = in; }

  //--------------Initialization begin-------------------------------
  void pre_regrid();
  void post_regrid(MeshChangeReason reason = MeshChangeReason::Topology);

  void distribute_arrays(const amrex::Vector<amrex::BoxArray> &cGridsOld =
                             amrex::Vector<amrex::BoxArray>());

  void fill_new_cells();
  void fill_E_B_fields(bool fillAll = false);

  void fill_new_node_E(bool fillAll = false);

  void fill_new_node_B(bool fillAll = false);
  void fill_new_center_B(bool fillAll = false);

  void fill_particles();

  // Narrow facade used by InitialCondition plug-ins (see InitialCondition.h).
  // These expose only the field arrays an IC needs; the facade is NOT a friend
  // of Pic.
  amrex::MultiFab &get_node_E(int iLev) { return nodeE[iLev]; }
  amrex::MultiFab &get_node_B(int iLev) { return nodeB[iLev]; }
  amrex::MultiFab &get_center_B(int iLev) { return centerB[iLev]; }
  PicICFields ic_fields() { return PicICFields(*this); }

  void init_source(const FluidInterface &interfaceIn) {
    // To be implemented

    //   sourceInterface = interfaceIn;
  }

  //----------------Initialization end-------------------------------

  void charge_exchange();

  void sum_moments(bool updateDt = false);

  void calc_mach_number();
  // Convert SI input parameters to code units after normalization is finalized.
  void finalize_units_conversion();
  void convert_resistivity();
  void convert_electron_density0();
  void convert_electron_heat_conduction();
  void convert_electron_collision();
  void convert_inflow_state();
  // Turn the SI input of #DIPOLE / #CRUSTALFIELD into code units and read the
  // spherical harmonic coefficients.
  void convert_intrinsic_B();
  // Fill the frozen B0 arrays (nodeB0 / centerB0) from the analytic field.
  void fill_intrinsic_B();
  // Regional resistivity / hyper-resistivity setup and update methods.
  void set_region_shapes(const amrex::Vector<std::shared_ptr<Shape> > &shapes) {
    if (useHybridPIC && (hasRegionalResistivity_ || hasRegionalHyper_)) {
      regionShapes = shapes;
    }
  }
  void init_regional_fields();
  void init_regional_fields(int iLev);
  void fill_regional_resistivity_field(int iLev);
  void fill_regional_hyper_field(int iLev);
  void update_regional_hyper_grid_mode(amrex::Real dt);
  // dst <- dst + B0, for a MultiFab living on the same centering as dst.
  void add_intrinsic_B(amrex::MultiFab &dst, int iLev);
  // Cell-centered total field: returns `src + B0` built in a scratch array, or
  // `src` itself when no intrinsic field is configured.
  amrex::MultiFab &total_center_B(amrex::MultiFab &src, int iLev);
  const IntrinsicBField *get_intrinsic_B() const { return intrinsicB_.get(); }

  void calc_mass_matrix();
  void calc_mass_matrix_amr();
  void init_boundary_node_mm_comm(int iLev);
  void sum_boundary_node_mm(int iLev);

  void update_part_loc_to_half_stage();

  void particle_mover();

  void re_sampling();

  void fill_source_particles();

  void inject_particles_for_new_cells() {
    if (!usePIC)
      return;

    for (auto &pts : parts) {
      pts->add_particles_domain();
    }
  }

  void inject_particles_for_boundary_cells() {
    if (!usePIC)
      return;

    for (auto &pts : parts) {
      pts->inject_particles_at_boundary();
    }
  }

  //------------Coupler related begin--------------
  void update_cells_for_pt();
  void get_fluid_state_for_points(const int nDim, const int nPoint,
                                  const double *const xyz_I,
                                  double *const data_I, const int nVar);
  void read_param(const std::string &command, ReadParam &param);
  void post_process_param();

  void report_bc_warnings(const std::string &context);
  void apply_periodicity_autofill(const amrex::Geometry &gm);
  void validate_bc_pairing(const amrex::Geometry &gm);
  //------------Coupler related end--------------

  //-------------Electric field solver begin-------------
  void update_E();
  void update_E_impl();
  void update_E_expl();
  void solve_E_gmres(int iLev);
  void solve_E_newton_krylov(int iLev);
  void update_E_rhs(double *rhos, int iLev);
  void update_E_matvec(const double *vecIn, double *vecOut, int iLev,
                       const bool useZeroBC = true);
  void update_E_M_dot_E(const amrex::MultiFab &inMF, amrex::MultiFab &outMF,
                        int iLev);
  void convert_1d_to_3d(const double *const p, amrex::MultiFab &MF, int iLev);
  void convert_3d_to_1d(const amrex::MultiFab &MF, double *const p, int iLev);

  void smooth_E(amrex::MultiFab &mfE, int iLev);
  void project_down_E();

  void smooth_multifab(amrex::MultiFab &mf, int iLev, int di,
                       amrex::Real coef = 0.5);

  void update_U0_E0();

  //-------------Hybrid PIC solver (kinetic ions + fluid electrons)-------------
  void smooth_moments();
  void update_B_hybrid();
  // Apply periodic and physical boundary conditions (and coarse-fine interface
  // ghosts on refined levels) to the cell-centered B, e.g. for intermediate RK
  // trial states that need fresh ghosts for the Ohm's law stencils.
  void apply_centerB_BC(int iLev);
  void apply_centerB_BC(int iLev, amrex::MultiFab &mfB);
  // Evaluate the Ohm's law E = -U_i x B + eta J + (J x B)/rho_q -
  // grad(Pe)/rho_q at an off-member B state (J from `centerBin`,
  // Hall/convection B from `centerBtimeAvg`), writing E into `Eout`. Ion
  // moments are time-interpolated between nodePlasmaPrev (J^{n-1/2}) and
  // nodePlasma (J^{n+1/2}) at the sub-step fraction `hstep`: X =
  // (0.5-hstep)X^{n-1/2} + (0.5+hstep)X^{n+1/2}.
  void assemble_ohm_E(const amrex::MultiFab &centerBin,
                      const amrex::MultiFab &centerBtimeAvg,
                      amrex::MultiFab &Eout, int iLev, amrex::Real hstep,
                      bool includeAmbi = true);
  void compute_ambipolar_E();
  void compute_ambipolar_E(int iLev);
  void save_current_moments_to_prev();
  void seed_first_hybrid_step();

  //-------------Evolved electron pressure (#ELECTRONPRESSURE)-----------------
  // Advance the scalar electron pressure by one PIC step, operator split:
  //   dPe/dt + div(u_e Pe) + (gamma_e-1) Pe div(u_e)
  //       = (gamma_e-1) [ div(kappa_hat . grad(Te)) + H_ei ]
  // update_Pe_hybrid(iLev, dt) is the orchestrator; the terms of the split are
  // the methods below, in the order they are applied.
  void update_Pe_hybrid();
  void update_Pe_hybrid(int iLev, amrex::Real dt);
  // u_e = U_i - J/(e*n_e) evaluated on the nodes.
  void electron_velocity_at_nodes(int iLev);
  // TVD/MUSCL advection of Pe plus the compression (pdV) term.
  void advect_electron_pressure(int iLev, amrex::Real dt);
  // Spitzer electron heat conduction; a no-op when heatCondKappa0 is 0.
  void apply_electron_heat_conduction(int iLev, amrex::Real dt);
  // Seed Pe from the algebraic polytropic closure using the current density.
  void init_electron_pressure();
  void init_electron_pressure(int iLev);
  // Interpolate the state onto boxes created by a regrid.
  void fill_new_electron_pressure();
  // Zero-gradient ghosts, used by the pressure field and scratch.
  void apply_pe_zero_gradient_bc(int iLev, amrex::MultiFab &mf);
  void apply_centerPe_BC(int iLev);
  // Fake-2D (single z cell) ghost clamp for a single-component cell-centered
  // field.
  void apply_fake2d_k_clamp(amrex::MultiFab &mf);
  // Nodal ion density -> cell-centered n_e. useGrownTile also covers the
  // coarse-fine interface nodes, which compute_ambipolar_E needs.
  void compute_electron_density(int iLev, amrex::MultiFab &nodalRho,
                                amrex::MultiFab &cellRho, bool useGrownTile);
  // Te = Pe/n_e at the cell centers, as conductivity scratch for the solver.
  void compute_electron_temperature(int iLev);
  // Electron-ion collisional thermal equilibration (heat exchange) hook:
  // dPe/dt = (Pi - Pe) / tau_eq, point-implicit formulation from BATSRUS.
  void add_electron_ion_heating(int iLev, amrex::Real dt);

  //-------------Electric field solver end-------------

  void update_B();

  void correct_B(int iLev);

  void solve_hyp_phi(int iLev);

  //-------------div(B) diagnostic begin----------------
  void ensure_divB(int iLev);
  void ensure_centerDivB(int iLev);
  void ensure_hypPhi(int iLev);
  // Refresh both div(B) diagnostics from the current B.
  void compute_divB(int iLev);
  //-------------div(B) diagnostic end------------------

  //-------------div(E) correction begin----------------
  void divE_correction();
  void amr_divE_correction();
  void divE_accurate_matvec(const double *vecIn, double *vecOut, int iLev);
  void divE_correct_particle_position();
  void sum_to_center(bool isBeforeCorrection);
  void sum_to_center_amr(bool isBeforeCorrection, int iLev);
  void calculate_phi(LinearSolver &solver, int iLev, bool reportSolver = true);
  //-------------div(E) correction end----------------

  void report_load_balance(bool doReportSummary = true,
                           bool doReportDetail = false);

  void calc_cost_per_cell();

  //--------------- IO begin--------------------------------
  void find_output_list(const PlotWriter &writerIn, long int &nPointAllProc,
                        VectorPointList &pointList_II, amrex::RealVect &xMin_D,
                        amrex::RealVect &xMax_D);

  void get_field_var(const VectorPointList &pointList_II,
                     const std::vector<std::string> &sVar_I,
                     MDArray<double> &var_II);
  double get_var(std::string_view var, const int iLev, const amrex::IntVect ijk,
                 const amrex::MFIter &mfi, bool isValidMFI = true);
  void save_restart_header(std::ofstream &headerFile);
  void save_restart_data();
  amrex::Vector<std::array<int, 3> > read_select_particle_input();
  void read_restart();
  void write_log(bool doForce = false, bool doCreateFile = false);
  void write_plots(bool doForce = false);
  void write_amrex(const PlotWriter &pw, double const timeNow,
                   int const iCycle);
  void write_amrex_field(const PlotWriter &pw, double const timeNow,
                         int const iCycle,
                         const std::string plotVars = "X E B plasma",
                         const std::string filenameIn = std::string(),
                         const amrex::BoxArray baOut = amrex::BoxArray());
  void write_amrex_particle(const PlotWriter &pw, double const timeNow,
                            int const iCycle);

  void set_IO_geom(amrex::Vector<amrex::Geometry> &geomIO,
                   const PlotWriter &pw);
  //--------------- IO end--------------------------------

  //--------------- Boundary begin ------------------------
  // Dispatch boundary conditions for electromagnetic fields (`isB` selects B vs
  // E).
  void apply_field_bc(const amrex::iMultiFab &status, amrex::MultiFab &mf,
                      const int iStart, const int nComp, GETVALUE func,
                      const int iLev, const bool isB);

  // Generic fill: zero-gradient copy on open faces, or values from `func`.
  void apply_BC(const amrex::iMultiFab &status, amrex::MultiFab &mf,
                const int iStart, const int nComp, GETVALUE func,
                const int iLev, const BoxBC<FieldBC::Type> *bc = nullptr);

  // Conducting wall: zero normal B / tangential E, mirror remaining components.
  void apply_conducting_wall(const amrex::iMultiFab &status,
                             amrex::MultiFab &mf, const int iStart,
                             const int nComp, const int iLev,
                             const BoxBC<FieldBC::Type> &bc, bool isB);

  // Absorbing wall: matched-impedance ghost-cell blending.
  void apply_absorbing_wall(const amrex::iMultiFab &status, amrex::MultiFab &mf,
                            const int iStart, const int nComp, const int iLev,
                            const BoxBC<FieldBC::Type> &bc, bool isB);

  // Inflow wall: pin physical boundary face nodes and fill ghost cells/nodes.
  void apply_inflow_wall(const amrex::iMultiFab &status, amrex::MultiFab &mf,
                         const int iStart, const int nComp, const int iLev,
                         const BoxBC<FieldBC::Type> &bc, bool isB,
                         GETVALUE func = nullptr);

  //--- Inner body field boundary (see #BODY / #BODYBOUNDARY) ---
  // These act on the nodes/cells flagged with bit::iBody_ and use the radial
  // direction from the body center as the surface normal.
  //
  // The body is split into a one-cell-thick surface layer, where the field
  // boundary condition acts, and the interior, which is a cavity with E = 0 and
  // B frozen at its initial value.
  // linetied   : E = 0 on every body node.
  // conducting : E <- (E.n) n on the surface nodes (E_t = 0, E_r kept) and
  //              E = 0 in the interior;
  //              B <- B - (B.n) n on the surface cells/nodes (B_r = 0).
  // insulating : nothing (the fields pass through the body).
  void zero_body_E(amrex::MultiFab &mf, const int iLev);
  void zero_body_interior_E(amrex::MultiFab &mf, const int iLev);
  void project_body_E(amrex::MultiFab &mf, const int iLev);
  void project_body_B(amrex::MultiFab &mf, const int iLev);
  void fill_body_E_insulating(amrex::MultiFab &mf, const int iLev);

  // Dispatch the electric-field condition selected by #BODYBOUNDARY.
  void apply_body_E_bc(amrex::MultiFab &mf, const int iLev);

  // Inject wave source into boundary ghost cells (iField: 0 = B, 1 = E).
  void apply_wave_field(const amrex::iMultiFab &status, amrex::MultiFab &mf,
                        const int iStart, const int nComp, const int iLev,
                        const BoxBC<FieldBC::Type> &bc, int iField,
                        amrex::Real t, GETVALUE func = nullptr);

  // Compute wave velocity perturbation for particle injection.
  void wave_velocity_kick(const amrex::Real *pos, amrex::Real t,
                          amrex::Real &dvx, amrex::Real &dvy, amrex::Real &dvz);

  // Fill external Dirichlet (coupled / fixed) ghost cells from func.
  void fill_ext_dir(const amrex::iMultiFab &status, amrex::MultiFab &mf,
                    const int iStart, const int nComp, GETVALUE func,
                    const int iLev, const BoxBC<FieldBC::Type> *bc = nullptr);

  amrex::Real get_zero(amrex::MFIter &mfi, amrex::IntVect ijk, int iVar,
                       int iLev) {
    return 0.0;
  }

  inline amrex::Real get_node_fluid_u(amrex::MFIter &mfi, amrex::IntVect ijk,
                                      int iVar, const int iLev, int iFluid) {
    amrex::Real u;
    if (iVar == ix_)
      u = fi->get_fluid_ux(mfi, ijk, iFluid, iLev);
    if (iVar == iy_)
      u = fi->get_fluid_uy(mfi, ijk, iFluid, iLev);
    if (iVar == iz_)
      u = fi->get_fluid_uz(mfi, ijk, iFluid, iLev);

    return u;
  }

  inline amrex::Real get_node_E(amrex::MFIter &mfi, amrex::IntVect ijk,
                                int iVar, const int iLev) {
    amrex::Real e;
    if (iVar == ix_)
      e = fi->get_ex(mfi, ijk, iLev);
    if (iVar == iy_)
      e = fi->get_ey(mfi, ijk, iLev);
    if (iVar == iz_)
      e = fi->get_ez(mfi, ijk, iLev);

    return e;
  }

  inline amrex::Real get_node_B(amrex::MFIter &mfi, amrex::IntVect ijk,
                                int iVar, const int iLev) {
    amrex::Real b;
    if (iVar == ix_)
      b = fi->get_bx(mfi, ijk, iLev);
    if (iVar == iy_)
      b = fi->get_by(mfi, ijk, iLev);
    if (iVar == iz_)
      b = fi->get_bz(mfi, ijk, iLev);

    return b;
  }

  inline amrex::Real get_center_B(amrex::MFIter &mfi, amrex::IntVect ijk,
                                  int iVar, const int iLev) {
    return fi->get_center_b(mfi, ijk, iVar, iLev);
  }

  inline amrex::Real get_center_E(amrex::MFIter &mfi, amrex::IntVect ijk,
                                  int iVar, const int iLev) {
    amrex::Real e;
    if (iVar == ix_)
      e = fi->get_ex(mfi, ijk, iLev);
    if (iVar == iy_)
      e = fi->get_ey(mfi, ijk, iLev);
    if (iVar == iz_)
      e = fi->get_ez(mfi, ijk, iLev);

    return e;
  }

  //--------------- Boundary end ------------------------

  void WriteDivEErrorToParaView() {
    amrex::Vector<amrex::MultiFab> errorDivE;
    errorDivE.resize(n_lev());
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      errorDivE[iLev].define(cGrids[iLev], DistributionMap(iLev), 1, nGst);
      errorDivE[iLev].setVal(0.0);

      for (amrex::MFIter mfi(errorDivE[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box &box = mfi.validbox();
        const amrex::Array4<amrex::Real> &error = errorDivE[iLev][mfi].array();
        const amrex::Array4<amrex::Real const> divEcc =
            centerDivE[iLev][mfi].array();
        const amrex::Array4<amrex::Real const> qcc =
            centerNetChargeN[iLev][mfi].array();
        const auto &status = cell_status(iLev)[mfi].array();

        amrex::ParallelFor(box, [&](int i, int j, int k) {
          error(i, j, k) =
              sqrt(pow((4.0 * dPI * qcc(i, j, k) - 1.0 * divEcc(i, j, k)), 2));
          if (bit::is_refined(status(i, j, k))) {
            error(i, j, k) = 0;
          }
        });
      }
    }
    grid.WriteMF(errorDivE, finest_level, "errorDivE");
  }

  void SetTargetPPC(int npresplitcells) {
    if (pInfo.isPPVconstant && pInfo.doPreSplitting) {
      amrex::Abort(
          "ConstantPPV and PreSplitting cannot be true at the same time");
    }
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      for (amrex::MFIter mfi(targetPPC[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box &box = mfi.fabbox();
        const auto &ppcArr = targetPPC[iLev][mfi].array();
        amrex::ParallelFor(box, [&](int i, int j, int k) noexcept {
          amrex::IntVect ijk = { AMREX_D_DECL(i, j, k) };
          ppcArr(ijk, 0) = product(pInfo.nPartPerCell);
        });
      }
    }
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      for (amrex::MFIter mfi(targetPPC[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box &box = mfi.validbox();
        const auto &ppcArr = targetPPC[iLev][mfi].array();
        const auto &status = cell_status(iLev)[mfi].array();
        amrex::ParallelFor(box, [&](int i, int j, int k) noexcept {
          amrex::IntVect ijk = { AMREX_D_DECL(i, j, k) };
          if (pInfo.isPPVconstant) {
            int tmp = 1;
            for (int i = 0; i < nDim; i++) {
              tmp *=
                  (pInfo.nPartPerCell[i] / pow((ref_ratio[iLev].max()), iLev));
            }
            ppcArr(ijk, 0) = tmp;
          } else {
            ppcArr(ijk, 0) = product(pInfo.nPartPerCell);
          }
          if (pInfo.doPreSplitting) {
            for (int ii = -npresplitcells; ii <= npresplitcells; ii++) {
              for (int jj = -npresplitcells; jj <= npresplitcells; jj++) {
                for (int kk = -npresplitcells; kk <= npresplitcells; kk++) {
                  amrex::IntVect ijk2 =
                      ijk + amrex::IntVect{ AMREX_D_DECL(ii, jj, kk) };
                  if (bit::is_refined(status(ijk2)) &&
                      !bit::is_refined(status(ijk))) {
                    ppcArr(ijk, 0) = product(pInfo.nPartPerCell) *
                                     pow(ref_ratio[iLev].max(), nDim);
                  }
                }
              }
            }
          }
        });
      }
    }
  }

  void WriteParticleQualityToParaView() {
    parts[0]->calculate_particle_quality(particleQuality);
    grid.WriteMF(particleQuality, finest_level, "particleQuality0");
    parts[1]->calculate_particle_quality(particleQuality);
    grid.WriteMF(particleQuality, finest_level, "particleQuality1");
  }
  // private methods
private:
  amrex::Real calc_E_field_energy();
  amrex::Real calc_B_field_energy();
};

void find_output_list_caller(const PlotWriter &writerIn,
                             long int &nPointAllProc,
                             VectorPointList &pointList_II,
                             amrex::RealVect &xMin_D, amrex::RealVect &xMax_D);

void get_field_var_caller(const VectorPointList &pointList_II,
                          const std::vector<std::string> &sVar_I,
                          MDArray<double> &var_II);

#endif
