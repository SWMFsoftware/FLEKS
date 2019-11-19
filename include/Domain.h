#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <iostream>

#include "Array1D.h"
#include "Constants.h"
#include "DomainGrid.h"
#include "FluidInterface.h"
#include "LinearSolver.h"
#include "Particles.h"
#include "TimeCtr.h"
#include "UMultiFab.h"

class Domain;

class FieldSolver {
public:
  amrex::Real theta;
  amrex::Real coefDiff;
  FieldSolver() {
    theta = 0.5;
    coefDiff = 0;
  }
};

typedef amrex::Real (Domain::*GETVALUE)(amrex::MFIter &mfi, int i, int j, int k,
                                        int iVar);

// The grid is defined in DomaiGrid. This class contains the data on the grid.
class Domain : public DomainGrid {
  friend PlotWriter;
  // public variables
public:
  // private variables
  TimeCtr tc;

private:
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

  amrex::Vector<std::unique_ptr<Particles> > parts;
  int *npcelx, *npcely, *npcelz;
  double *qom;

  FieldSolver fsolver;

  FluidInterface fluidInterface;

  bool doRestart;

  // public methods
public:
  Domain() {
    qom = nullptr;
    npcelx = nullptr;
    npcely = nullptr;
    npcelz = nullptr;

    doRestart = false;
  };
  ~Domain() {
    delete[] qom;
    delete[] npcelx;
    delete[] npcely;
    delete[] npcelz;
  };

  void update();

  //--------------Initialization begin-------------------------------
  void init(amrex::Real timeIn, const std::string &paramString, int *paramInt,
            double *gridDim, double *paramReal, int iDomain = 1);
  void make_grid();

  void make_data();
  void set_ic();
  void set_ic_field();
  void set_ic_particles();
  void init_time_ctr();
  //----------------Initialization end-------------------------------

  void sum_moments();

  void particle_mover();

  //------------Coupler related begin--------------
  void set_state_var(double *data, int *index);
  int get_grid_nodes_number();
  void get_grid(double *pos_DI);
  void read_param();
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

  void save_restart();
  void save_restart_header();
  void save_restart_data();
  void read_restart();
  //--------------- IO end--------------------------------

  //--------------- Boundary begin ------------------------
  void apply_external_BC(amrex::MultiFab &mf, const int iStart, const int nComp,
                         GETVALUE func);
  amrex::Real get_zero(amrex::MFIter &mfi, int i, int j, int k, int iVar) {
    return 0.0;
  }

  amrex::Real get_node_E(amrex::MFIter &mfi, int i, int j, int k, int iVar) {
    amrex::Real e;
    if (iVar == ix_)
      e = fluidInterface.get_ex(mfi, i, j, k);
    if (iVar == iy_)
      e = fluidInterface.get_ey(mfi, i, j, k);
    if (iVar == iz_)
      e = fluidInterface.get_ez(mfi, i, j, k);

    return e;
  }

  //--------------- Boundary end ------------------------

  // private methods
private:
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
