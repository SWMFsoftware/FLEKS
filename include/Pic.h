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

  bool useExplicitPIC = false;
  bool projectDownEmFields = true;
  bool skipMassMatrix = false;

  PartMode pMode = PartMode::PIC;

  FluidInterface *fi = nullptr;
  FluidInterface *stateOH = nullptr;
  FluidInterface *sourcePT2OH = nullptr;
  SourceInterface *source = nullptr;
  TimeCtr *tc = nullptr;

  amrex::Vector<amrex::MultiFab> nodeE;
  amrex::Vector<amrex::MultiFab> nodeEth;
  amrex::Vector<amrex::MultiFab> nodeB;
  amrex::Vector<amrex::MultiFab> divB;
  amrex::Vector<amrex::MultiFab> centerB;
  amrex::Vector<amrex::MultiFab> dBdt;

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
  amrex::Vector<amrex::Vector<amrex::MultiFab> > nodePlasma;
  amrex::Vector<amrex::Real> plasmaEnergy;

  bool isMomentsUpdated = false;

  amrex::Vector<amrex::MultiFab> jHat;

  amrex::Vector<std::unique_ptr<PicParticles> > parts;
  amrex::Vector<std::unique_ptr<PicParticles> > sourceParts;

  amrex::IntVect nPartPerCell = { AMREX_D_DECL(6, 6, 6) };
  amrex::Real qomEl = -100;

  // Particle Per Cell (PPC) of source particles.
  amrex::IntVect nSourcePPC = { AMREX_D_DECL(0, 0, 0) };
  bool adaptiveSourcePPC = false;
  bool kineticSource = false;

  FieldSolver fsolver;

  bool doCorrectDivE = true;
  int nDivECorrection = 3;

  bool doReSampling = true;
  amrex::Real reSamplingLowLimit = 0.8;
  amrex::Real reSamplingHighLimit = 1.5;
  amrex::Real maxWeightRatio = 1.0;

  bool solveFieldInCoMov = false;
  bool solvePartInCoMov = false;
  int nSmoothBackGroundU = 0;

  bool useUpwindE = false;
  amrex::Real limiterThetaE = 0;
  amrex::Real cMaxE = 0.0;
  bool useUpwindB = false;

  bool doSmoothJ = false;
  int nSmoothJ = 0;
  amrex::Real smoothJThreshold = 1;

  bool doSmoothE = false;
  int nSmoothE = 0;
  bool doSmoothB = false;
  amrex::Real limiterTheta = 0;
  amrex::Real smoothBIso = 0.1;

  TestCase testCase = RegularSimulation;

  ParticlesInfo pInfo;

  OHIon ionOH;

  // Boundary conditions for particles.
  amrex::Vector<BC> pBCs;

  // Boundary conditions for fields.
  BC bcBField;

  amrex::Vector<int> supIDs;

  // select particle params
  bool doSelectParticle = false;
  std::string selectParticleInputFile;

  bool doReport = false;

  std::string logFile;

  // public methods
public:
  Pic(amrex::Geometry const &gm, amrex::AmrInfo const &amrInfo, int nGst,
      FluidInterface *fluidIn, TimeCtr *tcIn, int id = 0)
      : Grid(gm, amrInfo, nGst, id, "pic"), fi(fluidIn), tc(tcIn) {
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

    // At most 10 species.
    pBCs.resize(10);

#ifdef _PT_COMPONENT_
    kineticSource = true;
    initEM = false;
    solveEM = false;

    doCorrectDivE = false;

    pMode = PartMode::Neutral;
#endif
  };
  ~Pic() {};

  void update(bool doReportIn = false);

  PicParticles *get_particle_pointer(int i) { return parts[i].get(); }

  void set_stateOH(OHInterface *in) { stateOH = in; }
  void set_sourceOH(OHInterface *in) { sourcePT2OH = in; }
  void set_fluid_source(SourceInterface *in) { source = in; }

  //--------------Initialization begin-------------------------------
  virtual void pre_regrid() override;
  virtual void post_regrid() override;

  void distribute_arrays(const amrex::Vector<amrex::BoxArray> &cGridsOld);

  void fill_new_cells();
  void fill_E_B_fields();
  void fill_lightwaves(amrex::Real wavelength, int EorB = -1,
                       amrex::Real time = 0, int lev = -1);

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

  void calc_mach_number();

  void calc_mass_matrix();
  void calc_mass_matrix_amr();
  void calc_mass_matrix_new();

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

  void smooth_E(amrex::MultiFab &mfE, int iLev);
  void project_down_E();

  void smooth_multifab(amrex::MultiFab &mf, int iLev, int di,
                       amrex::Real coef = 0.5);

  void update_U0_E0();
  //-------------Electric field solver end-------------

  void update_B();

  void smooth_B(int iLev);

  void solve_hyp_phi(int iLev);

  //-------------div(E) correction begin----------------
  void divE_correction();
  void amr_divE_correction();
  void divE_accurate_matvec(const double *vecIn, double *vecOut, int iLev);
  void divE_correct_particle_position();
  void sum_to_center(bool isBeforeCorrection);
  void sum_to_center_new(bool isBeforeCorrection, int iLev);
  void calculate_phi(LinearSolver &solver, int iLev);
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
  double get_var(std::string var, const int iLev, const amrex::IntVect ijk,
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
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      for (amrex::MFIter mfi(targetPPC[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box &box = mfi.validbox();
        const auto &ppcArr = targetPPC[iLev][mfi].array();
        const auto &status = cell_status(iLev)[mfi].array();
        amrex::ParallelFor(box, [&](int i, int j, int k) noexcept {
          amrex::IntVect ijk = { AMREX_D_DECL(i, j, k) };
          ppcArr(ijk, 0) = product(nPartPerCell);
          ppcArr(ijk, 1) = nPartPerCell[ix_];
          ppcArr(ijk, 2) = nPartPerCell[iy_];
          ppcArr(ijk, 3) = 0;
          if (nDim == 3) {
            ppcArr(ijk, 3) = nPartPerCell[iz_];
          }
          if (doPreSplitting) {
            for (int ii = -npresplitcells; ii <= npresplitcells; ii++) {
              for (int jj = -npresplitcells; jj <= npresplitcells; jj++) {
                for (int kk = -npresplitcells; kk <= npresplitcells; kk++) {
                  amrex::IntVect ijk2 =
                      ijk + amrex::IntVect{ AMREX_D_DECL(ii, jj, kk) };
                  if (bit::is_refined(status(ijk2)) &&
                      !bit::is_refined(status(ijk))) {
                    ppcArr(ijk, 0) = product(nPartPerCell) *
                                     pow(ref_ratio[iLev].max(), nDim);
                    ppcArr(ijk, 1) = nPartPerCell[ix_] * ref_ratio[iLev].max();
                    ppcArr(ijk, 2) = nPartPerCell[iy_] * ref_ratio[iLev].max();
                    if (nDim == 3) {
                      ppcArr(ijk, 3) =
                          nPartPerCell[iz_] * ref_ratio[iLev].max();
                    }
                  }
                }
              }
            }
          }
        });
      }
      if (isPPVconstant) {
        for (int iLev = 0; iLev < n_lev(); iLev++) {
          targetPPC[iLev].mult(
              1.0 / pow((pow(ref_ratio[iLev].max(), nDim)), iLev), 0, 1, 0);
          targetPPC[iLev].mult(1.0 / pow(ref_ratio[iLev].max(), iLev), 1, 3, 0);
        }
      }
    }
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
