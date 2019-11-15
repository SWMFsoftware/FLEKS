#ifndef _DOMAINGRID_H_
#define _DOMAINGRID_H_

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

class DomainGrid {
  // This class define the grid information, but NOT the data on the grid.

protected:
  int nGst;

  amrex::IntVect nCell;
  amrex::IntVect nNode;
  int maxBlockSize;
  int periodicity[nDim];
  amrex::IntVect centerBoxLo;
  amrex::IntVect centerBoxHi;
  amrex::Box centerBox;
  amrex::RealBox boxRange;

  const int coord = 0; // Cartesian grid
  amrex::Geometry geom;

  amrex::BoxArray centerBA;
  amrex::BoxArray nodeBA;

  amrex::DistributionMapping dm;

public:
  DomainGrid() {
    maxBlockSize = 8;
    for (int i = 0; i < nDim; i++) {
      periodicity[i] = 0;
    }
  }

  void init() {
    // nCell and boxRange should have been set before calling this function.

    for (int i = 0; i < nDim; i++) {
      centerBoxLo[i] = 0;
      centerBoxHi[i] = nCell[i] - 1;
    }

    centerBox.setSmall(centerBoxLo);
    centerBox.setBig(centerBoxHi);

    geom.define(centerBox, &boxRange, coord, periodicity);

    centerBA.define(centerBox);
    centerBA.maxSize(maxBlockSize);

    dm.define(centerBA);

    nodeBA = convert(centerBA, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });

    amrex::Print() << "Domain:: Domain range = " << boxRange << std::endl;
    amrex::Print() << "Domain:: Total block #  = " << nodeBA.size() << std::endl;
  }

  void set_nGst(const int nGstIn) { nGst = nGstIn; }
  void set_maxBlockSize(const int in) { maxBlockSize = in; }
  void set_nCell(const amrex::IntVect& in) { nCell = in; }
  void set_boxRange(const amrex::RealBox& in) { boxRange = in; }
  void set_periodicity(const int iDir, const bool isPeriodic) {
    periodicity[iDir] = (isPeriodic ? 1 : 0);
  }
};

#endif