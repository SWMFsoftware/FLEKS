#ifndef _DOMAINGRID_H_
#define _DOMAINGRID_H_

#include <AMReX_BCRec.H>
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
#include "GridInfo.h"

class DomainGrid {

protected:
  int nGst;

  amrex::IntVect nCell;

  amrex::IntVect maxBlockSize;
  int periodicity[nDim];
  amrex::IntVect centerBoxLo;
  amrex::IntVect centerBoxHi;
  amrex::Box centerBox;
  amrex::RealBox domainRange;

  const int coord = 0; // Cartesian grid
  amrex::Geometry gm;

  GridInfo gridInfo;

  int iGrid = 1;
  int iDecomp = 1;

  int domainID;
  std::string printPrefix;
  std::string domainName;

  bool isGridInitialized = false;

public:
  DomainGrid() {
    for (int i = 0; i < nDim; i++) {
      periodicity[i] = 0;
      maxBlockSize[i] = 8;
    }
  }
  ~DomainGrid() = default;

  int get_iGrid() const { return iGrid; }
  int get_iDecomp() const { return iDecomp; }
  void set_periodicity(const int iDir, const bool isPeriodic) {
    periodicity[iDir] = (isPeriodic ? 1 : 0);
  }
};
#endif
