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

  // const int coord = 0; // Cartesian grid

  // A collection of boxes to describe the simulation domain. The boxes have
  // been combined if possible. It covers the same region as cGrid, but
  // usually contains less boxes.
  amrex::BoxArray activeRegionBA;

  // Cell center
  amrex::BoxArray& cGrid = grids[0];
  amrex::Vector<amrex::BoxArray>& cGrids = grids;

  // Nodal
  amrex::Vector<amrex::BoxArray> nGrids;
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

private:
  // Here is the inheritance chain: AmrInfo -> AmrMesh -> AmrCore -> Grid. We
  // need to copy Grid object sometime, but the copy constructor of AmrCore is
  // deleted. I also do not know how to copy the information in AmrInfo. In
  // order to simplify the copy operation, the following variable is used to
  // install the initial AmrInfo.
  amrex::AmrInfo gridAmrInfo;

public:
  Grid(amrex::Geometry const& gm, amrex::AmrInfo const& amrInfo,
       const int nGstIn, int id, std::string tagIn = std::string())
      : AmrCore(gm, amrInfo), nGst(nGstIn), gridID(id) {
    gridAmrInfo = amrInfo;
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

  int get_n_ghost() const { return nGst; }

  const amrex::AmrInfo& get_amr_info() const { return gridAmrInfo; }

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

  void calc_node_grids() {
    nGrids.clear();
    nGrid.clear();

    if (cGrid.empty())
      return;

    nGrids.resize(finest_level + 1);
    for (int lev = 0; lev <= finest_level; lev++) {
      nGrids[lev] =
          amrex::convert(cGrids[lev], amrex::IntVect::TheNodeVector());
    }
    nGrid = nGrids[0];
  }

  void print_grid_info(bool printBoxes = false) {
    amrex::Print() << printPrefix << " =======Grid Info========" << std::endl;
    amrex::Print() << printPrefix << " max_level = " << max_level << std::endl;
    amrex::Print() << printPrefix << " finest_level = " << finest_level
                   << std::endl;

    for (int lev = 0; lev <= finest_level; lev++) {
      amrex::Print() << printPrefix << " lev = " << lev
                     << "\t # of boxes = " << std::setw(9) << cGrids[lev].size()
                     << "\t # of cells = " << std::setw(11) << CountCells(lev)
                     << "\t max_grid_size = " << max_grid_size[lev]
                     << std::endl;
    }

    if (printBoxes) {
      for (int lev = 0; lev <= finest_level; lev++) {
        amrex::Print() << printPrefix << " Boxes of lev = " << lev << std::endl;
        for (int ii = 0, n = cGrids[lev].size(); ii < n; ii++) {
          amrex::Print() << printPrefix << " box " << ii << " = "
                         << cGrids[lev][ii] << std::endl;
        }
      }
    }

    amrex::Print() << printPrefix << " =========================\n"
                   << std::endl;
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

    const int tagval = amrex::TagBox::SET;
    // return;
    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      const auto tagfab = tags.array(mfi);

      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (i >= 16 && i < 24 && j >= 8 && j < 16) {
              tagfab(i, j, k) = 1;
            }
          }
    }
  };

  virtual void PostProcessBaseGrids(amrex::BoxArray& ba) const override {
    std::string nameFunc = "Grid::PostProcessBaseGrids";
    amrex::Print() << printPrefix << nameFunc << " is called." << std::endl;
    ba = cGrid;
  };
};
#endif
