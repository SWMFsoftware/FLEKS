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
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include "Constants.h"


// This class define the grid information, but NOT the data on the grid.
class PicGrid {

protected:
  int nGst;

  amrex::Box centerBox;

  // const int coord = 0; // Cartesian grid
  amrex::Geometry geom;

  amrex::BoxArray centerBAOld;
  amrex::BoxArray centerBA;

  amrex::BoxArray nodeBA;

  amrex::DistributionMapping dm;

  amrex::MultiFab costMF;

  amrex::iMultiFab cellStatus;

public:
  PicGrid() = default;
  ~PicGrid() = default;

  void set_nGst(const int nGstIn) { nGst = nGstIn; }

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