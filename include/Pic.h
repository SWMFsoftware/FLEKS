#ifndef _PIC_H_
#define _PIC_H_

#include <iostream>

#include "Array1D.h"
#include "Constants.h"
#include "FluidInterface.h"
#include "Grid.h"
#include "GridUtility.h"
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
  FieldSolver() {
    theta = 0.51;
    coefDiff = 0.1;
  }
};

typedef amrex::Real (Pic::*GETVALUE)(amrex::MFIter &mfi, int i, int j, int k,
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

  bool useExplicitPIC = false;

  // If there is neutral species (OH-PT coupling), do not solve for
  // EM fields.
  bool solveEM = true;

  FluidInterface *fi = nullptr;
  FluidInterface *stateOH = nullptr;
  FluidInterface *sourcePT2OH = nullptr;
  SourceInterface *source = nullptr;
  TimeCtr *tc = nullptr;

  amrex::Vector<amrex::MultiFab> nodeE;
  amrex::Vector<amrex::MultiFab> nodeEth;
  amrex::Vector<amrex::MultiFab> nodeB;
  amrex::Vector<amrex::MultiFab> centerB;

  amrex::Vector<amrex::MultiFab> U0;
  amrex::Vector<amrex::MultiFab> E0;

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
  amrex::Vector<amrex::Vector<amrex::MultiFab> > nodePlasma;
  amrex::Vector<amrex::Real> plasmaEnergy;

  bool isMomentsUpdated = false;

  amrex::Vector<amrex::MultiFab> jHat;

  amrex::Vector<std::unique_ptr<Particles<> > > parts;

  amrex::IntVect nPartPerCell = { AMREX_D_DECL(6, 6, 6) };
  amrex::Real qomEl = -100;

  // Particle Per Cell (PPC) of source particles.
  amrex::IntVect nSourcePPC = { AMREX_D_DECL(0, 0, 0) };
  bool adaptiveSourcePPC = false;

  FieldSolver fsolver;

  bool doCorrectDivE = true;
  int nDivECorrection = 3;

  bool doReSampling = true;
  amrex::Real reSamplingLowLimit = 0.8;
  amrex::Real reSamplingHighLimit = 1.5;
  amrex::Real maxWeightRatio = 0;

  bool doSplitJ = false;
  bool doSmoothE = false;
  int nSmoothE = 1;
  amrex::Real coefStrongSmooth = 0.5;
  amrex::Real coefWeakSmooth = 0;
  amrex::Real strongSmoothMach = 0.8;
  amrex::Real weakSmoothMach = 0.7;

  amrex::MultiFab nodeSmoothCoef;

  TestCase testCase = RegularSimulation;

  amrex::Real particleMergeThreshold = -1, particleMergeBinBuffer = -1;
  bool fastMerge = false;
  int nMergeOld = 0, nMergeNew = 0;
  int nMergeTry = 1;
  amrex::Real mergeRatioMax = 1.5;
  amrex::Real pLevRatio = 1.2;
  bool mergeLight = false;
  amrex::Real mergePartRatioMax = 0;

  bool middlePointSource = false;

  // [amu/cc]
  amrex::Real vacuum = 0;

  bool doReport = false;

  std::string logFile;

  // Boundary conditions for particles.
  amrex::Vector<BC> pBCs;

  // public methods
public:
  Pic(amrex::Geometry const &gm, amrex::AmrInfo const &amrInfo, int nGst,
      FluidInterface *fluidIn, TimeCtr *tcIn, int id = 0)
      : Grid(gm, amrInfo, nGst, id, "pic"), fi(fluidIn), tc(tcIn) {
    eSolver.set_tol(1e-6);
    eSolver.set_nIter(200);

    divESolver.set_tol(0.01);
    divESolver.set_nIter(20);

    Particles<>::particlePosition = Staggered;

    //-----------------------------------------------------
    centerB.resize(n_lev_max());
    nodeB.resize(n_lev_max());
    nodeE.resize(n_lev_max());
    nodeEth.resize(n_lev_max());

    E0.resize(n_lev_max());
    U0.resize(n_lev_max());

    centerNetChargeOld.resize(n_lev_max());
    centerNetChargeN.resize(n_lev_max());
    centerNetChargeNew.resize(n_lev_max());

    centerDivE.resize(n_lev_max());
    centerPhi.resize(n_lev_max());

    nodeMM.resize(n_lev_max());
    centerMM.resize(n_lev_max());

    jHat.resize(n_lev_max());

    // At most 10 species.
    pBCs.resize(10);
  };
  ~Pic(){};

  void update(bool doReportIn = false);

  Particles<> *get_particle_pointer(int i) { return parts[i].get(); }

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

  void init_source(const FluidInterface &interfaceIn) {
    // To be implemented

    //   sourceInterface = interfaceIn;
  }

  //----------------Initialization end-------------------------------

  void charge_exchange();

  void sum_moments(bool updateDt = false);

  void calc_mass_matrix();

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
  void update_E_rhs(double *rhos, int iLev);
  void update_E_matvec(const double *vecIn, double *vecOut, int iLev,
                       const bool useZeroBC = true);
  void update_E_M_dot_E(const amrex::MultiFab &inMF, amrex::MultiFab &outMF,
                        int iLev);

  void smooth_E(amrex::MultiFab &mfE);

  void smooth_multifab(amrex::MultiFab &mf, bool useFixedCoef = false,
                       double coefIn = 1);

  void calc_smooth_coef();

  void update_U0_E0();
  //-------------Electric field solver end-------------

  void update_B();

  //-------------div(E) correction begin----------------
  void divE_correction();
  void divE_accurate_matvec(const double *vecIn, double *vecOut);
  void divE_correct_particle_position();
  void sum_to_center(bool isBeforeCorrection);
  void calculate_phi(LinearSolver &solver);
  //-------------div(E) correction end----------------

  void report_load_balance(bool doReportSummary = true,
                           bool doReportDetail = false);

  void calc_cost_per_cell(BalanceStrategy balanceStrategy);

  void convert_1d_to_3d(const double *const p, amrex::MultiFab &MF, int iLev);

  void convert_3d_to_1d(const amrex::MultiFab &MF, double *const p, int iLev);

  //--------------- IO begin--------------------------------
  void find_output_list(const PlotWriter &writerIn, long int &nPointAllProc,
                        VectorPointList &pointList_II, amrex::RealVect &xMin_D,
                        amrex::RealVect &xMax_D);

  void get_field_var(const VectorPointList &pointList_II,
                     const std::vector<std::string> &sVar_I,
                     MDArray<double> &var_II);
  double get_var(std::string var, const int iLev, const int ix, const int iy,
                 const int iz, const amrex::MFIter &mfi,
                 bool isValidMFI = true);
  void save_restart_header(std::ofstream &headerFile);
  void save_restart_data();
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
                const int iLev);

  amrex::Real get_zero(amrex::MFIter &mfi, int i, int j, int k, int iVar,
                       int iLev) {
    return 0.0;
  }

  inline amrex::Real get_node_E(amrex::MFIter &mfi, int i, int j, int k,
                                int iVar, const int iLev) {
    amrex::Real e;
    if (iVar == ix_)
      e = fi->get_ex(mfi, i, j, k, iLev);
    if (iVar == iy_)
      e = fi->get_ey(mfi, i, j, k, iLev);
    if (iVar == iz_)
      e = fi->get_ez(mfi, i, j, k, iLev);

    return e;
  }

  inline amrex::Real get_node_B(amrex::MFIter &mfi, int i, int j, int k,
                                int iVar, const int iLev) {
    amrex::Real b;
    if (iVar == ix_)
      b = fi->get_bx(mfi, i, j, k, iLev);
    if (iVar == iy_)
      b = fi->get_by(mfi, i, j, k, iLev);
    if (iVar == iz_)
      b = fi->get_bz(mfi, i, j, k, iLev);

    return b;
  }

  inline amrex::Real get_center_B(amrex::MFIter &mfi, int i, int j, int k,
                                  int iVar, const int iLev) {
    return fi->get_center_b(mfi, i, j, k, iVar, iLev);
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
