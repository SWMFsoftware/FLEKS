#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include "Constants.h"

using namespace amrex;

struct DomainParams{  
  IntVect nCell; 
  int nCellBlockMax;
  int periodicity[nDim];
  RealBox realLo; 
  RealBox realHi; 
};

class Domain {
  // public variables
public:
  // private variables
private:
  DomainParams domainPar;
  // public methods
public:
  Domain(){};
  ~Domain(){};
  void init();
  // private methods
private:  
};

#endif
