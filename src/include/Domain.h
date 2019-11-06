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
#include "Particles.h"
#include "UMultiFab.h"
#include "LinearSolver.h"

using namespace amrex;

struct Tag {
  int index;
  bool b1;
};

class FieldSolver {
public:
  Real theta;
  Real coefDiff;
  FieldSolver() {
    theta = 0.5;
    coefDiff = 0;
  }
};

class Domain {
  // public variables
public:
  // private variables
private:
  int iProc;

  int nGst;

  IntVect nCell;
  IntVect nNode;
  int nCellBlockMax;
  int periodicity[nDim];
  IntVect centerBoxLo;
  IntVect centerBoxHi;
  Box centerBox;
  RealBox boxRange;

  int coord;
  Geometry geom;

  BoxArray centerBA;
  BoxArray nodeBA;

  DistributionMapping dm;

  MultiFab nodeE;
  MultiFab nodeEth;
  MultiFab nodeB;
  MultiFab centerB;

  UMultiFab<RealMM> nodeMM;

  // ------divE correction--------------
  // Old @ t=t_{n-1/2}; N @ t=t_n; New @ t=t_{n+1/2}
  MultiFab centerNetChargeOld, centerNetChargeN, centerNetChargeNew;
  MultiFab centerDivE, centerPhi;
  UMultiFab<RealCMM> centerMM;
  const Real rhoTheta = 0.51; 
  //--------------------------------------

  int nSolveNode;
  int nSolveCenter; 

  //------Temporary variables for field---
  MultiFab tempNode3;
  MultiFab tempCenter3;
  MultiFab tempCenter1;
  MultiFab tempCenter1_1;  
  //--------------------------------------

  int nSpecies;
  int iTot;
  Vector<MultiFab> nodePlasma;

  Vector<std::unique_ptr<Particles> > parts;
  int *npcelx, *npcely, *npcelz;
  double *qom;

  FieldSolver fsolver;

  Real dtSI, dt;
  Real timeNow, timeNowSI;

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
  void init(Real timeIn, const std::string &paramString, int *paramInt,
            double *gridDim, double *paramReal, int iDomain = 1);
  void define_domain();
  void set_ic();
  void init_field();
  void init_particles();
  void sum_moments();

  void divE_correction();
  void divE_accurate_matvec(double *vecIn, double *vecOut);
  void divE_correct_particle_position();
  void sum_to_center(bool isBeforeCorrection);
  void calculate_phi(MATVEC fMatvec, Real tol=1e-2, int nIter=20);

  void update();

  void particle_mover();

  void set_dtSI(double dtIn) {
    dtSI = dtIn;
    dt = dtSI;
  }
  Real get_dtSI() const { return dtSI; }

  Real get_timeSI() const {
    Print() << "timeNowSI = " << timeNowSI << std::endl;
    return timeNowSI;
  }

  //------------Coupler related begin--------------
  void set_state_var(double *data, int *index);
  int get_grid_nodes_number();
  void get_grid(double *pos_DI);

  void read_param();
  //------------Coupler related end--------------

  void update_E();
  void update_E_rhs(double *rhos);
  void update_E_matvec(const double *vecIn, double *vecOut);
  void update_E_M_dot_E(const MultiFab &inMF, MultiFab &outMF);

  void update_B();

  // private methods
private:
};

#endif
