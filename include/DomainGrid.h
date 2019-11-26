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

// This class define the grid information, but NOT the data on the grid.
class DomainGrid {

protected:
  int nGst;

  amrex::IntVect nCell;
  amrex::IntVect nNode;
  amrex::IntVect maxBlockSize;
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
    for (int i = 0; i < nDim; i++) {
      periodicity[i] = 0;
      maxBlockSize[i] = 8;
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

    amrex::Print() << "DomainGrid:: Domain range = " << boxRange << std::endl;
    amrex::Print() << "DomainGrid:: Total block #  = " << nodeBA.size()
                   << std::endl;
    //amrex::Print() << "DomainGrid:: centerBA = " << centerBA << std::endl;
  }

  void set_nGst(const int nGstIn) { nGst = nGstIn; }
  void set_maxBlockSize(int iDir, const int in) { maxBlockSize[iDir] = in; }
  void set_nCell(const amrex::IntVect& in) { nCell = in; }
  void set_boxRange(const amrex::RealBox& in) { boxRange = in; }
  void set_periodicity(const int iDir, const bool isPeriodic) {
    periodicity[iDir] = (isPeriodic ? 1 : 0);
  }

  inline int find_mpi_rank_from_coord(amrex::Real const x, amrex::Real const y,
                                      amrex::Real const z) const {
    amrex::Real loc[3] = { x, y, z };
    auto idx = geom.CellIndex(loc);
    return find_mpi_rank_from_cell_index(idx[ix_], idx[iy_], idx[iz_]);
  }

  inline int find_mpi_rank_from_cell_index(int const i, int const j,
                                           int const k) const {
    amrex::IntVect idx = { i, j, k };
    for (int ii = 0, n = centerBA.size(); ii < n; ii++) {
      const amrex::Box& bx = centerBA[ii];
      if (bx.contains(idx))
        return dm[ii];
    }

    amrex::Abort("Error: can not find this cell!");
    return -1; // To suppress compiler warnings.
  }
};

#endif