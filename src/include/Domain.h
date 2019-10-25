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

#include "Constants.h"
#include "Particles.h"
#include "FluidPicInterface.h"

using namespace amrex;

struct Tag{
  int index; 
  bool b1; 
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
  MultiFab nodeB;
  MultiFab centerB;
  MultiFab nodeMMatrix;

  int nSpecies;
  int iTot;
  Vector<MultiFab> nodePlasma;  

  Vector<std::unique_ptr<Particles>> parts;
  int *npcelx, *npcely, *npcelz;
  double *qom;  


  Real dtSI, dt;
  Real timeNow, timeNowSI;

  FluidPicInterface fluidInterface;

  // public methods
public:
  Domain(){};
  ~Domain(){};
  void init(Real timeIn, const std::string &paramString, int *paramInt,
            double *gridDim, double *paramReal, int iDomain = 1);
  void define_domain();
  void set_ic();
  void init_field();
  void init_particles();
  void sum_moments();
  void update();

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

  // private methods
private:
};

#endif
