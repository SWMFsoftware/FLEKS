#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <iostream>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include "Array1D.h"
#include "Constants.h"
#include "FluidPicInterface.h"
#include "LinearSolver.h"
#include "Particles.h"
#include "TimeCtr.h"
#include "UMultiFab.h"

class FieldSolver {
public:
  amrex::Real theta;
  amrex::Real coefDiff;
  FieldSolver() {
    theta = 0.5;
    coefDiff = 0;
  }
};

class Domain {
  // public variables
public:
  // private variables
  TimeCtr tc;

private:
  int iProc;

  int nGst;

  amrex::IntVect nCell;
  amrex::IntVect nNode;
  int nCellBlockMax;
  int periodicity[nDim];
  amrex::IntVect centerBoxLo;
  amrex::IntVect centerBoxHi;
  amrex::Box centerBox;
  amrex::RealBox boxRange;

  int coord;
  amrex::Geometry geom;

  amrex::BoxArray centerBA;
  amrex::BoxArray nodeBA;

  amrex::DistributionMapping dm;

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

  int nSolveNode;
  int nSolveCenter;

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

  FluidPicInterface fluidInterface;

  // public methods
public:
  Domain() {
    qom = nullptr;
    npcelx = nullptr;
    npcely = nullptr;
    npcelz = nullptr;
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
  void define_domain();
  void set_ic();
  void init_field();
  void init_particles();
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
  void update_E_matvec(const double *vecIn, double *vecOut);
  void update_E_M_dot_E(const amrex::MultiFab &inMF, amrex::MultiFab &outMF);
  //-------------Electric field solver end-------------

  void update_B();

  //-------------div(E) correction begin----------------
  void divE_correction();
  void divE_accurate_matvec(double *vecIn, double *vecOut);
  void divE_correct_particle_position();
  void sum_to_center(bool isBeforeCorrection);
  void calculate_phi(MATVEC fMatvec, amrex::Real tol = 1e-2, int nIter = 20);
  //-------------div(E) correction end----------------

  void find_output_list(const PlotWriter &writerIn, long int &nPointAllProc,
                        PlotWriter::VectorPointList &pointList_II,
                        std::array<double, nDimMax> &xMin_D,
                        std::array<double, nDimMax> &xMax_D);

  void get_field_var(const VectorPointList &pointList_II,
                     const std::vector<std::string> &sVar_I,
                     MDArray<double> &var_II);

  // private methods
private:
};

#endif
