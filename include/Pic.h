#ifndef _PIC_H_
#define _PIC_H_

#include <iostream>

#include "Array1D.h"
#include "Constants.h"
#include "FluidInterface.h"
#include "Grid.h"
#include "LinearSolver.h"
#include "Particles.h"
#include "ReadParam.h"
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
                                     int iVar);

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

  std::shared_ptr<FluidInterface> fluidInterface;
  std::shared_ptr<TimeCtr> tc;

  amrex::MultiFab nodeE;
  amrex::MultiFab nodeEth;
  amrex::MultiFab nodeB;
  amrex::MultiFab centerB;

  amrex::UMultiFab<RealMM> nodeMM;

  amrex::MultiFab costMF;

  // ------divE correction--------------
  // Old @ t=t_{n-1/2}; N @ t=t_n; New @ t=t_{n+1/2}
  amrex::MultiFab centerNetChargeOld, centerNetChargeN, centerNetChargeNew;
  amrex::MultiFab centerDivE, centerPhi;
  amrex::UMultiFab<RealCMM> centerMM;
  const amrex::Real rhoTheta = 0.51;
  //--------------------------------------

  LinearSolver eSolver;
  LinearSolver divESolver;

  //------Temporary variables for field---
  amrex::MultiFab tempNode3;
  amrex::MultiFab tempCenter3;
  amrex::MultiFab tempCenter1;
  amrex::MultiFab tempCenter1_1;
  //--------------------------------------

  bool useSource = false;
  // FluidInterface sourceInterface;

  int nSpecies;
  int iTot;
  amrex::Vector<amrex::MultiFab> nodePlasma;
  amrex::Vector<amrex::Real> plasmaEnergy;

  amrex::MultiFab jHat;

  amrex::Vector<std::unique_ptr<Particles<> > > parts;

  amrex::IntVect nPartPerCell = { 6, 6, 6 };
  amrex::Real qomEl = -100;

  FieldSolver fsolver;

  bool doCorrectDivE = true;
  int nDivECorrection = 3;

  bool doReSampling = true;
  amrex::Real reSamplingLowLimit = 0.8;
  amrex::Real reSamplingHighLimit = 1.5;

  bool doSmoothE = false;
  int nSmoothE = 1;
  amrex::Real coefStrongSmooth = 0.5;
  amrex::Real coefWeakSmooth = 0;
  amrex::Real strongSmoothMach = 0.8;
  amrex::Real weakSmoothMach = 0.7;

  amrex::MultiFab nodeSmoothCoef;

  TestCase testCase = RegularSimulation;

  amrex::Real particleMergeThreshold = -1, particleMergeBinBuffer = -1;

  bool doReport = false;

  // public methods
public:
  Pic(amrex::Geometry const &gm, amrex::AmrInfo const &amrInfo,
      std::shared_ptr<FluidInterface> &fluidIn, std::shared_ptr<TimeCtr> &tcIn,
      int id = 0)
      : Grid(gm, amrInfo, id), tc(tcIn), fluidInterface(fluidIn) {
    eSolver.set_tol(1e-6);
    eSolver.set_nIter(200);

    divESolver.set_tol(0.01);
    divESolver.set_nIter(20);

    Particles<>::particlePosition = Staggered;
  };
  ~Pic(){};

  void update(bool doReportIn = false);

  Particles<> *get_particle_pointer(int i) { return parts[i].get(); }

  //--------------Initialization begin-------------------------------
  void set_geom(int nGstIn, const amrex::Geometry &geomIn);

  void regrid(const amrex::BoxArray &activeRegionBAIn,
              const amrex::BoxArray &centerBAIn,
              const amrex::DistributionMapping &dmIn);

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
      pts->add_particles_domain(cellStatus);
    }
  }

  void inject_particles_for_boundary_cells() {
    if (!usePIC)
      return;

    for (auto &pts : parts) {
      pts->inject_particles_at_boundary(cellStatus);
    }
  }

  //------------Coupler related begin--------------
  void set_state_var(double *data, int *index);
  int get_grid_nodes_number();
  void get_grid(double *pos_DI);
  void find_mpi_rank_for_points(const int nPoint, const double *const xyz_I,
                                int *const rank_I);
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
  void update_E_rhs(double *rhos);
  void update_E_matvec(const double *vecIn, double *vecOut,
                       const bool useZeroBC = true);
  void update_E_M_dot_E(const amrex::MultiFab &inMF, amrex::MultiFab &outMF);

  void smooth_E(amrex::MultiFab &mfE);

  void smooth_multifab(amrex::MultiFab &mf, bool useFixedCoef = false,
                       double coefIn = 1);

  void calc_smooth_coef();
  //-------------Electric field solver end-------------

  void update_B();

  //-------------div(E) correction begin----------------
  void divE_correction();
  void divE_accurate_matvec(const double *vecIn, double *vecOut);
  void divE_correct_particle_position();
  void sum_to_center(bool isBeforeCorrection);
  void calculate_phi(LinearSolver &solver);
  //-------------div(E) correction end----------------

  //---------------load balance begin-------------------
  void compute_cost();
  void load_balance();
  //---------------load balance end---------------------

  void report_load_balance();

  void set_nodeShare();

  void convert_1d_to_3d(const double *const p, amrex::MultiFab &MF);

  void convert_3d_to_1d(const amrex::MultiFab &MF, double *const p);

  //--------------- IO begin--------------------------------
  void find_output_list(const PlotWriter &writerIn, long int &nPointAllProc,
                        PlotWriter::VectorPointList &pointList_II,
                        std::array<double, nDim> &xMin_D,
                        std::array<double, nDim> &xMax_D);

  void get_field_var(const VectorPointList &pointList_II,
                     const std::vector<std::string> &sVar_I,
                     MDArray<double> &var_II);
  double get_var(std::string var, const int ix, const int iy, const int iz,
                 const amrex::MFIter &mfi, bool isValidMFI = true);
  void save_restart_header(std::ofstream &headerFile);
  void save_restart_data();
  void read_restart();
  std::string logFile;
  void write_log(bool doForce = false, bool doCreateFile = false);
  void write_plots(bool doForce = false);
  void write_amrex(const PlotWriter &pw, double const timeNow,
                   int const iCycle);
  void write_amrex_field(const PlotWriter &pw, double const timeNow,
                         int const iCycle,
                         const std::string plotVars = "X E B plasma",
                         const std::string filenameIn = std::string());
  void write_amrex_particle(const PlotWriter &pw, double const timeNow,
                            int const iCycle);

  void set_IO_geom(amrex::Geometry &geomIO, const PlotWriter &pw);
  //--------------- IO end--------------------------------

  //--------------- Boundary begin ------------------------
  void apply_BC(const amrex::iMultiFab &status, amrex::MultiFab &mf,
                const int iStart, const int nComp, GETVALUE func = nullptr);

  amrex::Real get_zero(amrex::MFIter &mfi, int i, int j, int k, int iVar) {
    return 0.0;
  }

  inline amrex::Real get_node_E(amrex::MFIter &mfi, int i, int j, int k,
                                int iVar) {
    amrex::Real e;
    if (iVar == ix_)
      e = fluidInterface->get_ex(mfi, i, j, k);
    if (iVar == iy_)
      e = fluidInterface->get_ey(mfi, i, j, k);
    if (iVar == iz_)
      e = fluidInterface->get_ez(mfi, i, j, k);

    return e;
  }

  inline amrex::Real get_node_B(amrex::MFIter &mfi, int i, int j, int k,
                                int iVar) {
    amrex::Real b;
    if (iVar == ix_)
      b = fluidInterface->get_bx(mfi, i, j, k);
    if (iVar == iy_)
      b = fluidInterface->get_by(mfi, i, j, k);
    if (iVar == iz_)
      b = fluidInterface->get_bz(mfi, i, j, k);

    return b;
  }

  inline amrex::Real get_center_B(amrex::MFIter &mfi, int i, int j, int k,
                                  int iVar) {
    return fluidInterface->get_center_b(mfi, i, j, k, iVar);
  }

  //--------------- Boundary end ------------------------

  // Make a new level from scratch using provided BoxArray and
  // DistributionMapping. Only used during initialization. overrides the pure
  // virtual function in AmrCore
  virtual void MakeNewLevelFromScratch(
      int lev, amrex::Real time, const amrex::BoxArray &ba,
      const amrex::DistributionMapping &dm) override {
    std::string nameFunc = "Pic::MakeNewLevelFromScratch";
    amrex::Print() << printPrefix << nameFunc << " lev = " << lev
                   << " ba = " << ba << std::endl;
  };

  // Make a new level using provided BoxArray and DistributionMapping and
  // fill with interpolated coarse level data.
  // overrides the pure virtual function in AmrCore
  virtual void MakeNewLevelFromCoarse(
      int lev, amrex::Real time, const amrex::BoxArray &ba,
      const amrex::DistributionMapping &dm) override {
    std::string nameFunc = "Pic::MakeNewLevelFromCoarse";
    amrex::Print() << printPrefix << nameFunc << " lev = " << lev
                   << " ba = " << ba << std::endl;
  };

  virtual void PostProcessBaseGrids(amrex::BoxArray &ba) const override {
    std::string nameFunc = "Pic::PostProcessBaseGrids";
    amrex::Print() << printPrefix << nameFunc << " is called." << std::endl;
    ba = cGrid;
  };

  // private methods
private:
  amrex::Real calc_E_field_energy();
  amrex::Real calc_B_field_energy();
};

void find_output_list_caller(const PlotWriter &writerIn,
                             long int &nPointAllProc,
                             PlotWriter::VectorPointList &pointList_II,
                             std::array<double, nDim> &xMin_D,
                             std::array<double, nDim> &xMax_D);

void get_field_var_caller(const VectorPointList &pointList_II,
                          const std::vector<std::string> &sVar_I,
                          MDArray<double> &var_II);

#endif
