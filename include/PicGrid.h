#ifndef _PICGRID_H_
#define _PICGRID_H_

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
#include <AMReX_iMultiFab.H>

#include "Constants.h"

// This class define the grid information, but NOT the data on the grid.
class PicGrid {

protected:
  int nGst;

  // const int coord = 0; // Cartesian grid
  amrex::Geometry geom;

  // A collection of boxes to describe the PIC domain. The boxes have been
  // combined if possible. It covers the same region as centerBA, but usually
  // contains less boxes.
  amrex::BoxArray picRegionBA;

  amrex::BoxArray centerBAOld;
  amrex::BoxArray centerBA;

  amrex::BoxArray nodeBA;
  amrex::BoxArray nodeBAOld;

  amrex::DistributionMapping dm;

  amrex::MultiFab costMF;

  // The status of a cell could be: iBoundary_, iOnNew_, or iOnOld_.
  amrex::iMultiFab cellStatus;

  // The status of a Node could be: iBoundary_, iOnNew_, or iOnOld_.
  amrex::iMultiFab nodeStatus;

  // A node may be shared by a few blocks/boxes. Sometimes (such as the E field
  // solver) only one of the boexes needs to take care such nodes. The following
  // multifab is usually used for this purpose: if one node of one box is
  // 'iAssign_', this box needs to take care of this node (this node must be
  // 'iNotAssign_' on other boxes).
  amrex::iMultiFab nodeAssignment;
  constexpr static int iAssign_ = 1, iNotAssign_ = 0;

  bool doNeedFillNewCell;

  bool isGridInitialized;

  bool isGridEmpty;

public:
  PicGrid() {
    doNeedFillNewCell = true;
    isGridInitialized = false;
    isGridEmpty = false;
  };
  ~PicGrid() = default;

  void set_nGst(const int nGstIn) { nGst = nGstIn; }

  bool is_grid_empty() const { return isGridEmpty; }

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

    amrex::AllPrint() << "idx = " << idx << std::endl;
    amrex::Abort("Error: can not find this cell!");
    return -1; // To suppress compiler warnings.
  }
};

#endif