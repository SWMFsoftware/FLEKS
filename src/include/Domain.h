#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <iostream>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>
#include <AMReX_IndexType.H>

#include "Constants.h"
#include "Particles.h"

using namespace amrex;

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
  Vector<MultiFab> nodePlasma;

  //Vector<MParticles> partVect;
  // public methods
public:
  Domain(){};
  ~Domain(){};
  void init();
  void define_domain();
  void init_field();
  void init_particles();
  
  //private methods
private:
};

#endif
