#ifndef _DOMAIN_H_
#define _DOMAIN_H_

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

#include "Constants.h"

using namespace amrex;

class Domain {
  // public variables
public:
  // private variables
private:
  int nGst;

  IntVect nCell;
  int nCellBlockMax;
  int periodicity[nDim];
  IntVect domainLo;
  IntVect domainHi;
  RealBox domainRange;
  Box domainBox;

  int coord;
  Geometry geom;

  BoxArray ba;

  DistributionMapping dm;

  MultiFab E;
  MultiFab B;

  int nSpecies;
  Vector<MultiFab> plasma;

  // public methods
public:
  Domain(){};
  ~Domain(){};
  void init();
  // private methods
private:
};

#endif
