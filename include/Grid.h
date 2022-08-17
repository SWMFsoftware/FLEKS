#ifndef _Grid_H_
#define _Grid_H_

#include <AMReX_AmrCore.H>
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
class Grid : public amrex::AmrCore {

protected:
  int nGst;

  // const int coord = 0; // Cartesian grid
  amrex::Geometry gm;

  // A collection of boxes to describe the simulation domain. The boxes have
  // been combined if possible. It covers the same region as centerBA, but
  // usually contains less boxes.
  amrex::BoxArray activeRegionBA;

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
  // 'iAssign_', this box needs to take care of this node; in fake 2D, a node
  // can be 'iAbandon_' to ignore one layer node; if a node is neither
  // 'iAssign_' nor 'iAbandon_', nodeShare stors the neighbor that handle
  // this node.
  amrex::iMultiFab nodeShare;
  constexpr static int iAssign_ = (1 << 4), iAbandon_ = -1;

  bool doNeedFillNewCell = true;

  bool isGridInitialized = false;

  bool isGridEmpty = false;

  bool isFake2D = false;

public:
  Grid(amrex::Geometry const& gm, amrex::AmrInfo const& amrInfo)
      : AmrCore(gm, amrInfo){};
  ~Grid() = default;

  void set_nGst(const int nGstIn) { nGst = nGstIn; }

  bool is_grid_empty() const { return isGridEmpty; }

  inline int find_mpi_rank_from_coord(amrex::Real const x, amrex::Real const y,
                                      amrex::Real const z) const {
    amrex::Real loc[3] = { x, y, z };
    auto idx = gm.CellIndex(loc);
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

  // Make a new level using provided BoxArray and DistributionMapping and
  // fill with interpolated coarse level data.
  // overrides the pure virtual function in AmrCore
  virtual void MakeNewLevelFromCoarse(
      int lev, amrex::Real time, const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm) override{
    // To be implemented
  };

  // Remake an existing level using provided BoxArray and DistributionMapping
  // and fill with existing fine and coarse data. overrides the pure virtual
  // function in AmrCore
  virtual void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray& ba,
                           const amrex::DistributionMapping& dm) override{
    // To be implemented
  };

  // Delete level data
  // overrides the pure virtual function in AmrCore
  virtual void ClearLevel(int lev) override{
    // To be implemented
  };

  // Make a new level from scratch using provided BoxArray and
  // DistributionMapping. Only used during initialization. overrides the pure
  // virtual function in AmrCore
  virtual void MakeNewLevelFromScratch(
      int lev, amrex::Real time, const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm) override{
    // To be implemented
  };

  // tag all cells for refinement
  // overrides the pure virtual function in AmrCore
  virtual void ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time,
                        int ngrow) override{
    // To be implemented
  };
};

#endif
