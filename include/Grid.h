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

  std::string printPrefix;
  std::string gridName;
  int gridID;

  bool useAMRGrid = false;

  // const int coord = 0; // Cartesian grid

  // A collection of boxes to describe the simulation domain. The boxes have
  // been combined if possible. It covers the same region as cGrid, but
  // usually contains less boxes.
  amrex::BoxArray activeRegionBA;

  // Cell center
  amrex::BoxArray cGrid;

  // Nodal
  amrex::BoxArray nGrid;

  // The status of a cell could be: iBoundary_, iOnNew_, or iOnOld_.
  amrex::iMultiFab cellStatus;

  // The status of a Node could be: iBoundary_, iOnNew_, or iOnOld_.
  amrex::iMultiFab nodeStatus;

  // A node may be shared by a few blocks/boxes. Sometimes (such as the E field
  // solver) only one of the boexes needs to take care such nodes. The following
  // multifab is usually used for the following purposes: (1) if one node of one
  // box is 'iAssign_', this box needs to take care of this node; (2) in fake
  // 2D, some nodes are set to 'iIgnore_' so that only one layer of nodes is
  // solved; (3) if a node is neither 'iAssign_' nor 'iIgnore_', nodeShare
  // stores the neighbor that handle this node.
  amrex::iMultiFab nodeShare;
  constexpr static int iAssign_ = (1 << 4), iIgnore_ = -1;

  bool doNeedFillNewCell = true;

  bool isGridInitialized = false;

  bool isGridEmpty = false;

  bool isFake2D = false;

  std::string tag;

public:
  Grid(amrex::Geometry const& gm, amrex::AmrInfo const& amrInfo, int nGstIn,
       int id, std::string tagIn = std::string())
      : AmrCore(gm, amrInfo), nGst(nGstIn), gridID(id) {
    tag = tagIn;
    gridName = std::string("FLEKS") + std::to_string(gridID);
    if (tag.empty()) {
      printPrefix = gridName + ": ";
    } else {
      printPrefix = gridName + " " + tag + ": ";
    }

    isFake2D = Geom(0).Domain().bigEnd(iz_) == Geom(0).Domain().smallEnd(iz_);
  };
  ~Grid() = default;

  void set_base_grid(const amrex::BoxArray& ba) { cGrid = ba; }

  bool is_grid_empty() const { return isGridEmpty; }

  inline int find_mpi_rank_from_coord(amrex::Real const x, amrex::Real const y,
                                      amrex::Real const z) const {
    amrex::Real loc[3] = { x, y, z };
    auto idx = Geom(0).CellIndex(loc);
    return find_mpi_rank_from_cell_index(idx[ix_], idx[iy_], idx[iz_]);
  }

  inline int find_mpi_rank_from_cell_index(int const i, int const j,
                                           int const k) const {
    amrex::IntVect idx = { i, j, k };
    for (int ii = 0, n = cGrid.size(); ii < n; ii++) {
      const amrex::Box& bx = cGrid[ii];
      if (bx.contains(idx))
        return DistributionMap(0)[ii];
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
      const amrex::DistributionMapping& dm) override {
    std::string nameFunc = "Grid::MakeNewLevelFromCoarse";
    amrex::Print() << printPrefix << nameFunc << " lev = " << lev << std::endl;
  };

  // Remake an existing level using provided BoxArray and DistributionMapping
  // and fill with existing fine and coarse data. overrides the pure virtual
  // function in AmrCore
  virtual void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray& ba,
                           const amrex::DistributionMapping& dm) override {
    std::string nameFunc = "Grid::RemakeLevel";
    amrex::Print() << printPrefix << nameFunc << " lev = " << lev << std::endl;
  };

  // Delete level data
  // overrides the pure virtual function in AmrCore
  virtual void ClearLevel(int lev) override {
    std::string nameFunc = "Grid::ClearLevel";
    amrex::Print() << printPrefix << nameFunc << " lev = " << lev << std::endl;
  };

  // Make a new level from scratch using provided BoxArray and
  // DistributionMapping. Only used during initialization. overrides the pure
  // virtual function in AmrCore
  virtual void MakeNewLevelFromScratch(
      int lev, amrex::Real time, const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm) override {
    std::string nameFunc = "Grid::MakeNewLevelFromScratch";
    amrex::Print() << printPrefix << nameFunc << " lev = " << lev << std::endl;
  };

  // tag all cells for refinement
  // overrides the pure virtual function in AmrCore
  virtual void ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time,
                        int ngrow) override {
    std::string nameFunc = "Grid::ErrorEst";
    amrex::Print() << printPrefix << nameFunc << " lev = " << lev << std::endl;
  };

  virtual void PostProcessBaseGrids(amrex::BoxArray& ba) const override {
    std::string nameFunc = "Grid::PostProcessBaseGrids";
    amrex::Print() << printPrefix << nameFunc << " is called." << std::endl;
    ba = cGrid;
  };
};
#endif
