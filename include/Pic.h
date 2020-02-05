#ifndef _PIC_H_
#define _PIC_H_

#include <iostream>

#include "Array1D.h"
#include "Constants.h"
#include "PicGrid.h"
#include "FluidInterface.h"
#include "LinearSolver.h"
#include "Particles.h"
#include "TimeCtr.h"
#include "UMultiFab.h"

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

// The grid is defined in DomaiGrid. This class contains the data on the grid.
class Pic : public PicGrid {
  friend PlotWriter;
  // public variables
public:
  // private variables

private:
  std::shared_ptr<FluidInterface> fluidInterface;
  std::shared_ptr<TimeCtr> tc;

  amrex::MultiFab nodeE;
  amrex::MultiFab nodeEth;
  amrex::MultiFab nodeB;
  amrex::MultiFab centerB;

  amrex::UMultiFab<RealMM> nodeMM;

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

  int nSpecies;
  int iTot;
  amrex::Vector<amrex::MultiFab> nodePlasma;
  amrex::Vector<amrex::Real> plasmaEnergy;

  amrex::Vector<std::unique_ptr<Particles> > parts;

  amrex::IntVect nPartPerCell;
  amrex::Real qomEl;

  FieldSolver fsolver;

  bool doRestart;
  bool doCorrectDivE;

  bool doReSampling;
  amrex::Real reSamplingLowLimit;
  amrex::Real reSamplingHighLimit;
  // public methods
public:
  Pic() {
    qomEl = -100;

    doRestart = false;
    doCorrectDivE = true;

    doReSampling = false;

    for (int iDim = 0; iDim < nDim; iDim++)
      nPartPerCell[iDim] = 6;

    eSolver.set_tol(1e-6);
    eSolver.set_nIter(200);

    divESolver.set_tol(0.01);
    divESolver.set_nIter(20);
  };
  ~Pic() {};

  void update();

  //--------------Initialization begin-------------------------------
  void init(amrex::Real timeIn, const std::string &paramString, int *paramInt,
            double *gridDim, double *paramReal,
            std::shared_ptr<FluidInterface> &fluidIn,
            std::shared_ptr<TimeCtr> &tcIn);

  void make_grid(int nGstIn, const amrex::BoxArray &centerBAIn,
                 const amrex::Geometry &geomIn);

  void make_data();
  void set_ic();
  void set_ic_field();
  void set_ic_particles();
  //----------------Initialization end-------------------------------

  void sum_moments();

  void particle_mover();

  //------------Coupler related begin--------------
  void set_state_var(double *data, int *index);
  int get_grid_nodes_number();
  void get_grid(double *pos_DI);
  void find_mpi_rank_for_points(const int nPoint, const double *const xyz_I,
                                int *const rank_I);
  void get_fluid_state_for_points(const int nDim, const int nPoint,
                                  const double *const xyz_I,
                                  double *const data_I, const int nVar);
  void read_param(const std::string& command, ReadParam& readParam);
  //------------Coupler related end--------------

  //-------------Electric field solver begin-------------
  void update_E();
  void update_E_rhs(double *rhos);
  void update_E_matvec(const double *vecIn, double *vecOut,
                       const bool useZeroBC = true);
  void update_E_M_dot_E(const amrex::MultiFab &inMF, amrex::MultiFab &outMF);
  //-------------Electric field solver end-------------

  void update_B();

  //-------------div(E) correction begin----------------
  void divE_correction();
  void divE_accurate_matvec(double *vecIn, double *vecOut);
  void divE_correct_particle_position();
  void sum_to_center(bool isBeforeCorrection);
  void calculate_phi(LinearSolver &solver);
  //-------------div(E) correction end----------------

  //---------------load balance begin-------------------
  void compute_cost();
  void load_balance();
  //---------------load balance end---------------------

  //--------------- IO begin--------------------------------
  void find_output_list(const PlotWriter &writerIn, long int &nPointAllProc,
                        PlotWriter::VectorPointList &pointList_II,
                        std::array<double, nDimMax> &xMin_D,
                        std::array<double, nDimMax> &xMax_D);

  void get_field_var(const VectorPointList &pointList_II,
                     const std::vector<std::string> &sVar_I,
                     MDArray<double> &var_II);
  double get_var(std::string var, const int ix, const int iy, const int iz,
                 const amrex::MFIter &mfi);
  
  void save_restart_header(std::ofstream &headerFile);
  void save_restart_data();
  void read_restart();
  std::string logFile;
  void write_log(bool doForce = false, bool doCreateFile = false);
  //--------------- IO end--------------------------------

  //--------------- Boundary begin ------------------------
  void apply_external_BC(amrex::MultiFab &mf, const int iStart, const int nComp,
                         GETVALUE func);
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

  // private methods
private:
  amrex::Real calc_E_field_energy();
  amrex::Real calc_B_field_energy();
};

void find_output_list_caller(const PlotWriter &writerIn,
                             long int &nPointAllProc,
                             PlotWriter::VectorPointList &pointList_II,
                             std::array<double, nDimMax> &xMin_D,
                             std::array<double, nDimMax> &xMax_D);

void get_field_var_caller(const VectorPointList &pointList_II,
                          const std::vector<std::string> &sVar_I,
                          MDArray<double> &var_II);

#endif
