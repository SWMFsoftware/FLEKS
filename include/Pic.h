#ifndef _PIC_H_
#define _PIC_H_

#include <iostream>

#include "Array1D.h"
#include "Bit.h"
#include "Constants.h"
#include "DomainParameters.h"
#include "FleksDistributionMap.h"
#include "FluidInterface.h"
#include "Grid.h"
#include "InitialCondition.h"
#include "LinearSolver.h"
#include "OHInterface.h"
#include "Particles.h"
#include "ReadParam.h"
#include "SourceInterface.h"
#include "TimeCtr.h"
#include "UMultiFab.h"

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

// The grid is defined in DomainGrid. This class contains the data on the grid.
class Pic : public Grid {
  friend PlotWriter;
  friend ParticleTracker;
  // private variables
private:
  bool usePIC = true;
  bool solveEM = true;
  bool initEM = true;

  // ---- Hybrid PIC (kinetic ions + fluid electrons) solver ----
  // Generalized Ohm's law field solver instead of the implicit GM-PC (solveEM).
  bool useHybridPIC = false;
  // Resistive term eta * J. SI input [m^2/s], converted to code units.
  amrex::Real etaResistivitySI = 0.0;
  amrex::Real etaResistivity = 0.0;
  // Electron pressure gradient. Input [eV], converted to code units.
  amrex::Real electronTemperatureEV = 0.0;
  amrex::Real electronTemperature = 0.0;
  // Polytropic index for the adiabatic electron pressure closure.
  amrex::Real electronGamma = 1.0;
  // Reference (upstream) charge density, input in amu/cc.
  amrex::Real electronDensity0EV = 1.0;
  // Reference charge density in code units; converted lazily from
  // electronDensity0EV at the first hybrid field advance.
  amrex::Real electronDensity0 = 0.0;
  // True once electronDensity0 (code units) has been converted.
  bool electronDensity0Converted_ = false;
  // Number of sub-steps for the B-field update within one coarse dt.
  int nBSubcycle = 1;
  // Hall term (J x B) / rho in the generalized Ohm's law. Default on.
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

  // Minimum charge density floor for the 1/rho factors in the Hall term and
  // the electron-pressure-gradient term. <= 0 means auto: 1e-6 *
  // electronDensity0.
  amrex::Real rhoMinOhm = 0.0;

  bool useExplicitPIC = false;
  bool projectDownEmFields = true;
  bool skipMassMatrix = false;
  bool reportParticleQuality = false;

  PartMode pMode = PartMode::PIC;

  FluidInterface *fi = nullptr;
  FluidInterface *stateOH = nullptr;
  FluidInterface *sourcePT2OH = nullptr;
  SourceInterface *source = nullptr;
  TimeCtr *tc = nullptr;

  const DomainParameters &domainParameters;

  amrex::Vector<amrex::MultiFab> nodeE;
  amrex::Vector<amrex::MultiFab> nodeEth;
  amrex::Vector<amrex::MultiFab> nodeB;
  amrex::Vector<amrex::MultiFab> divB;
  amrex::Vector<amrex::MultiFab> centerB;
  // Hyper-resistivity scratch fields (allocated when useHybridPIC).
  amrex::Vector<amrex::MultiFab> centerLapB; // nabla^2 B  (hyper stage A)
  amrex::Vector<amrex::MultiFab> nodeHyperE; // nabla x (nabla^2 B) (hyper stage
                                             // B)
  // RK4 (fieldIntegrator="rk4") scratch fields (allocated when useHybridPIC).
  // centerB_RK4 / nodeB_RK4 hold trial B states at the sub-stages; nodeE_RK4 is
  // the electric field evaluated at a trial B; kRK4[0..3] are the four stage
  // curls curl(E_stage). All reused per sub-step; only level 0 is advanced
  // (fine levels follow from projection, exactly as in the ssprk3 path).
  amrex::Vector<amrex::MultiFab> centerB_RK4;
  amrex::Vector<amrex::MultiFab> nodeB_RK4;
  amrex::Vector<amrex::MultiFab> nodeE_RK4;
  // kRK4[iLev][0..3]: the four stage curls curl(E_stage) for the level-iLev RK4
  // sub-step (4 stages per level).
  amrex::Vector<amrex::Vector<amrex::MultiFab> > kRK4;

  // rk3/rk4 persistent scratch (allocated when useHybridPIC). centerBstart
  // holds the sub-step start B_n; centerBstar holds the time-centred (trial +
  // B_n)/2 state used by the rk3/rk4 time-centred-E stages.
  amrex::Vector<amrex::MultiFab> centerBstart_heun;
  amrex::Vector<amrex::MultiFab> centerBstar_heun;

  amrex::Vector<amrex::MultiFab> dBdt;
  amrex::Vector<amrex::MultiFab> particleQuality;

  // Running time-averaged magnetic field (EMA). Allocated when useHybridPIC;
  // only used when useAvgFieldB is set. B_avg is NOT divergence-clean and is
  // never fed into the Faraday update -- it is used only inside the generalized
  // Ohm's law and in the particle Boris push.
  amrex::Vector<amrex::MultiFab> centerBavg; // cell-centred <B>
  amrex::Vector<amrex::MultiFab> nodeBavg;   // node-centred <B>
  bool isBavgInit = false;                   // first-step copy flag for the EMA

  // Hyperbolic cleaning
  bool useHyperbolicCleaning = false;
  amrex::Vector<amrex::MultiFab> hypPhi;
  amrex::Real hypDecay = 0.1;

  // Background velocity and electric field.
  amrex::Vector<amrex::MultiFab> uBg;
  amrex::Vector<amrex::MultiFab> eBg;

  // Mach number: u/v_th
  amrex::Vector<amrex::MultiFab> mMach;

  amrex::Vector<amrex::UMultiFab<RealMM> > nodeMM;

  // ------divE correction--------------
  // Old @ t=t_{n-1/2}; N @ t=t_n; New @ t=t_{n+1/2}
  amrex::Vector<amrex::MultiFab> centerNetChargeOld, centerNetChargeN,
      centerNetChargeNew;
  amrex::Vector<amrex::MultiFab> centerDivE, centerPhi;
  amrex::Vector<amrex::UMultiFab<RealCMM> > centerMM;
  const amrex::Real rhoTheta = 0.51;
  //--------------------------------------

  LinearSolver eSolver;
  LinearSolver divESolver;

  int nSpecies;
  int iTot;

  // Electron / kinetic-ion identification for the hybrid solver.
  // iElectron_ is the index of the (explicit) electron species, or -1 if none.
  // kineticSpecies_ lists the indices of all non-electron (kinetic ion)
  // species; in standard PIC runs (no electron particle species) this equals
  // all species, so existing behaviour is unchanged.
  int iElectron_ = -1;
  std::vector<int> kineticSpecies_;
  amrex::Vector<amrex::Vector<amrex::MultiFab> > nodePlasma;
  // Second-order hybrid: previous-step ion moments (J^{n-1/2}, rho^{n-1/2})
  // deposited BEFORE the particle push. Combined with the current deposit
  // nodePlasma (J^{n+1/2}) by a linear time interpolation inside
  // assemble_ohm_E, indexed by the magnetic sub-step fraction hstep
  // (hstep=0 -> J^n = 1/2(J^{n-1/2}+J^{n+1/2}); hstep=1 -> J^{n+1} =
  // 3/2 J^{n+1/2} - 1/2 J^{n-1/2}). Same layout as nodePlasma.
  amrex::Vector<amrex::Vector<amrex::MultiFab> > nodePlasmaPrev;
  // ---- Cell-centred hybrid solver fields (Hybrid-VPIC-style layout) ----
  // When useHybridPIC is true the hybrid time step reads/writes only these
  // cell-centred fields (B, E, J, moments). The node-centred nodeE/nodeB/
  // nodePlasma are kept as write-only OUTPUT MIRRORS refreshed once per step
  // (node-sync bridges) so the plot / restart / tracker path sees correct
  // data. All six arrays are allocated only inside the useHybridPIC block of
  // Pic::distribute_arrays.
  amrex::Vector<amrex::MultiFab> centerEhybrid; // cell-centred E
  amrex::Vector<amrex::MultiFab> centerJ;       // cell-centred J
  amrex::Vector<amrex::MultiFab> centerEprev;   // cell-centred E^n
                                                // (time-centring)
  amrex::Vector<amrex::MultiFab> centerBprev;   // cell-centred B^n
                                                // (time-centring)
  amrex::Vector<amrex::MultiFab> centerE_RK4;  // cell-centred E trial (replaces
                                               // nodeE_RK4)
  amrex::Vector<amrex::MultiFab> centerHyperE; // cell-centred hyper-resistivity
                                               // E
  amrex::Vector<amrex::Vector<amrex::MultiFab> > centerPlasma; // per-species
                                                               // cell-centred
                                                               // moments
  amrex::Vector<amrex::Vector<amrex::MultiFab> >
      centerPlasmaSum; // summed cell-centred moments [nSpecies][iLev]
  amrex::Vector<amrex::Vector<amrex::MultiFab> > centerPlasmaPrev; // per-species
                                                                   // previous-step
                                                                   // moments
  amrex::Vector<amrex::Real> plasmaEnergy;

  bool isMomentsUpdated = false;
  // When true, nodePlasma (and mMach) are stale: the per-step nodePlasma output
  // bridge and calc_mach_number are skipped in sum_moments, and must be
  // materialized on demand by sync_node_plasma_output() before any node-plasma
  // / mach output or load-balancing is requested. Only used in the hybrid path.
  bool nodePlasmaStale = false;

  // True when nodeE (an output mirror of the live centerEhybrid) is stale.
  // Marked stale each update_B_hybrid and materialized by sync_node_E_output()
  // at plot time, so the per-step average_center_to_node cost is deferred out
  // of the loop.
  bool nodeEStale = false;

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
  // Override upwind velocity in correct_B(). 0 = use plasma background velocity.
  amrex::Real fixedUpwindVel = 0.0;

  // Override uMax for CFL estimate. < 0 = estimate from particle thermal velocity.
  amrex::Real fixedUMax = -1.0;

  bool doSmoothJ = false;
  int nSmoothJ = 0;
  amrex::Real coefSmoothJ = 0.5;

  // Digital-filter smoothing of ion moments before the generalized Ohm's law.
  bool doSmoothMoments = false;
  int nSmoothMoments = 0;
  amrex::Real coefSmoothMoments = 0.5;

  // B integrator: "rk4" (default) or "ssprk3".
  std::string fieldIntegrator = "rk4";
  bool useRK4 = false;

  // Guard: true on the first hybrid step before nodePlasmaPrev is seeded.
  bool isFirstHybridStep = true;

  // EMA-averaged B fed to Ohm's law and Boris push. Dampens shot noise.
  bool useAvgFieldB = false;
  int nAvgFieldB = 10;

  bool doSmoothE = false;
  int nSmoothE = 0;

  // Plug-in initial condition via #TESTCASE registry.
  std::unique_ptr<InitialCondition> ic_;

  ParticlesInfo pInfo;

  // Boundary conditions for fields.
  BC bcBField;

  // select particle params
  bool doSelectParticle = false;
  std::string selectParticleInputFile;

  bool doReport = false;
  int dnMemory = -1;

  std::string logFile;

  // public methods
public:
  Pic(amrex::Geometry const &gm, amrex::AmrInfo const &amrInfo, int nGst,
      FluidInterface *fluidIn, TimeCtr *tcIn, int id,
      const DomainParameters &parameters)
      : Grid(gm, amrInfo, nGst, id, "pic"),
        fi(fluidIn),
        tc(tcIn),
        domainParameters(parameters) {
    eSolver.set_tol(1e-6);
    eSolver.set_nIter(200);

    divESolver.set_tol(0.01);
    divESolver.set_nIter(20);

    //-----------------------------------------------------
    centerB.resize(n_lev_max());
    nodeB.resize(n_lev_max());
    dBdt.resize(n_lev_max());
    nodeE.resize(n_lev_max());
    nodeEth.resize(n_lev_max());
    divB.resize(n_lev_max());
    hypPhi.resize(n_lev_max());
    centerLapB.resize(n_lev_max());
    nodeHyperE.resize(n_lev_max());
    centerB_RK4.resize(n_lev_max());
    nodeB_RK4.resize(n_lev_max());
    nodeE_RK4.resize(n_lev_max());
    centerEhybrid.resize(n_lev_max());
    centerJ.resize(n_lev_max());
    centerEprev.resize(n_lev_max());
    centerBprev.resize(n_lev_max());
    centerE_RK4.resize(n_lev_max());
    centerHyperE.resize(n_lev_max());
    centerBavg.resize(n_lev_max());
    nodeBavg.resize(n_lev_max());
    centerBstart_heun.resize(n_lev_max());
    centerBstar_heun.resize(n_lev_max());
    kRK4.resize(n_lev_max());
    for (int iL = 0; iL < n_lev_max(); ++iL)
      kRK4[iL].resize(4);
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

    nodeMM.resize(n_lev_max());
    centerMM.resize(n_lev_max());

    jHat.resize(n_lev_max());

#ifdef _PT_COMPONENT_
    kineticSource = true;
    initEM = false;
    solveEM = false;

    doCorrectDivE = false;

    pMode = PartMode::Neutral;
#endif
  };
  ~Pic() {};

  void free_memory();

  void update(bool doReportIn = false);

  PicParticles *get_particle_pointer(int i) { return parts[i].get(); }
  // TODO: no longer needed if we fix the output variable ordering.
  bool get_useHybridPIC() const { return useHybridPIC; }

  void set_stateOH(OHInterface *in) { stateOH = in; }
  void set_sourceOH(OHInterface *in) { sourcePT2OH = in; }
  void set_fluid_source(SourceInterface *in) { source = in; }

  //--------------Initialization begin-------------------------------
  virtual void pre_regrid() override;
  virtual void post_regrid() override;

  void distribute_arrays(const amrex::Vector<amrex::BoxArray> &cGridsOld);

  void fill_new_cells();
  void fill_E_B_fields();

  void fill_new_node_E();

  void fill_new_node_B();
  void fill_new_center_B();

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
  // Rebuild the stale nodePlasma/nodeE output mirrors; calc_mach_number only if
  // needMach. Used by output / load balancing, not the hybrid solver.
  void sync_node_plasma_output(bool needMach = false);
  void sync_node_E_output();

  // Convert electronDensity0 (amu/cc) to code units and set the auto density
  // floor. Idempotent; run at the first hybrid field advance after
  // fi->post_process_param() finalizes Si2NoRho.
  void convert_electron_density0();

  void calc_mass_matrix();
  void calc_mass_matrix_amr();

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

  void smooth_E(amrex::MultiFab &mfE, int iLev);
  void project_down_E();

  void smooth_multifab(amrex::MultiFab &mf, int iLev, int di,
                       amrex::Real coef = 0.5);

  void update_U0_E0();

  //-------------Hybrid PIC solver (kinetic ions + fluid electrons)-------------
  void smooth_moments();
  void update_B_hybrid();
  void project_centerB_to_nodeB(int iLev);
  void apply_centerB_BC(int iLev);
  void project_centerB_to_nodeB_scratch(amrex::MultiFab &centerIn,
                                        amrex::MultiFab &nodeOut, int iLev);
  // Evaluate the Ohm's law E = -U_i x B + eta J + (J x B)/rho_q -
  // grad(Pe)/rho_q at an off-member B state (J from `centerBin`,
  // Hall/convection B from `centerBtimeAvg`), writing E into `Eout`. Ion
  // moments are time-interpolated between centerPlasmaPrev (J^{n-1/2}) and
  // centerPlasmaSum (J^{n+1/2}) at the sub-step fraction `hstep`: X =
  // (0.5-hstep)X^{n-1/2} + (0.5+hstep)X^{n+1/2}.
  void assemble_ohm_E(const amrex::MultiFab &centerBin,
                      const amrex::MultiFab &centerBtimeAvg,
                      amrex::MultiFab &Eout, int iLev, amrex::Real hstep);
  void save_current_moments_to_prev();
  void seed_first_hybrid_step();

  //-------------Electric field solver end-------------

  void update_B();

  void correct_B(int iLev);

  void solve_hyp_phi(int iLev);

  //-------------div(E) correction begin----------------
  void divE_correction();
  void amr_divE_correction();
  void divE_accurate_matvec(const double *vecIn, double *vecOut, int iLev);
  void divE_correct_particle_position();
  void sum_to_center(bool isBeforeCorrection);
  void sum_to_center_amr(bool isBeforeCorrection, int iLev);
  void calculate_phi(LinearSolver &solver, int iLev);
  //-------------div(E) correction end----------------

  void report_load_balance(bool doReportSummary = true,
                           bool doReportDetail = false);

  void calc_cost_per_cell();

  void convert_1d_to_3d(const double *const p, amrex::MultiFab &MF, int iLev);

  void convert_3d_to_1d(const amrex::MultiFab &MF, double *const p, int iLev);

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
  void apply_BC(const amrex::iMultiFab &status, amrex::MultiFab &mf,
                const int iStart, const int nComp, GETVALUE func,
                const int iLev, const BC *bc = nullptr);

  bool use_float(const int i, const int j, const int k, int &ip, int &jp,
                 int &kp, const BC &bc, const amrex::Box &bxValid) {
    bool useFloat = false;
    ip = i;
    jp = j;
    kp = k;
    if (i < bxValid.smallEnd(ix_) && bc.lo[ix_] == BC::outflow) {
      useFloat = true;
      ip = bxValid.smallEnd(ix_);
    }
    if (i > bxValid.bigEnd(ix_) && bc.hi[ix_] == BC::outflow) {
      useFloat = true;
      ip = bxValid.bigEnd(ix_);
    }

    if (j < bxValid.smallEnd(iy_) && bc.lo[iy_] == BC::outflow) {
      useFloat = true;
      jp = bxValid.smallEnd(iy_);
    }

    if (j > bxValid.bigEnd(iy_) && bc.hi[iy_] == BC::outflow) {
      useFloat = true;
      jp = bxValid.bigEnd(iy_);
    }

    if (nDim > 2) {
      if (k < bxValid.smallEnd(iz_) && bc.lo[iz_] == BC::outflow) {
        useFloat = true;
        kp = bxValid.smallEnd(iz_);
      }

      if (k > bxValid.bigEnd(iz_) && bc.hi[iz_] == BC::outflow) {
        useFloat = true;
        kp = bxValid.bigEnd(iz_);
      }
    }
    return useFloat;
  }

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

  // Make a new level from scratch using provided BoxArray and
  // DistributionMapping. Only used during initialization. overrides the pure
  // virtual function in AmrCore
  virtual void MakeNewLevelFromScratch(
      int iLev, amrex::Real time, const amrex::BoxArray &ba,
      const amrex::DistributionMapping &dm) override {
    std::string nameFunc = "Pic::MakeNewLevelFromScratch";
    amrex::Print() << printPrefix << nameFunc << " iLev = " << iLev
                   << std::endl;
  };

  // Make a new level using provided BoxArray and DistributionMapping and
  // fill with interpolated coarse level data.
  // overrides the pure virtual function in AmrCore
  virtual void MakeNewLevelFromCoarse(
      int iLev, amrex::Real time, const amrex::BoxArray &ba,
      const amrex::DistributionMapping &dm) override {
    std::string nameFunc = "Pic::MakeNewLevelFromCoarse";
    amrex::Print() << printPrefix << nameFunc << " iLev = " << iLev
                   << std::endl;
  };

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
    WriteMF(errorDivE, finest_level, "errorDivE");
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
    WriteMF(particleQuality, finest_level, "particleQuality0");
    parts[1]->calculate_particle_quality(particleQuality);
    WriteMF(particleQuality, finest_level, "particleQuality1");
  }
  // private methods
private:
  amrex::Real calc_E_field_energy();
  amrex::Real calc_B_field_energy();
  AMREX_EXPORT amrex::UNode_FourthOrder<amrex::Real> node_fourth_order_interp;
};

void find_output_list_caller(const PlotWriter &writerIn,
                             long int &nPointAllProc,
                             VectorPointList &pointList_II,
                             amrex::RealVect &xMin_D, amrex::RealVect &xMax_D);

void get_field_var_caller(const VectorPointList &pointList_II,
                          const std::vector<std::string> &sVar_I,
                          MDArray<double> &var_II);

#endif
