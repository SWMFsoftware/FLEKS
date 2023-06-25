#ifndef _Grid_H_
#define _Grid_H_

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>
#include <AMReX_iMultiFab.H>

#include "Bit.h"
#include "Constants.h"
#include "GridUtility.h"
#include "Utility.h"

// This class define the grid information, but NOT the data on the grid.
class Grid : public amrex::AmrCore {

protected:
  int nGst;

  std::string printPrefix;
  std::string gridName;
  int gridID;

  int nLev = 1;

  // const int coord = 0; // Cartesian grid

  // A collection of boxes to describe the simulation domain. The boxes have
  // been combined if possible. It covers the same region as cGrids[0], but
  // usually contains less boxes.
  amrex::BoxArray activeRegion;

  // Cell center
  amrex::Vector<amrex::BoxArray>& cGrids = grids;

  // Nodal
  amrex::Vector<amrex::BoxArray> nGrids;

  // Each bit of the integer represents a status of a cell/node. See Bit.h for
  // potential status.
  amrex::Vector<amrex::iMultiFab> cellStatus;
  amrex::Vector<amrex::iMultiFab> nodeStatus;

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

    nLev = max_level + 1;

    cellStatus.resize(nLev);
    nodeStatus.resize(nLev);
  };

  ~Grid() = default;

  int get_n_ghost() const { return nGst; }

  const amrex::AmrInfo& get_amr_info() const { return gridAmrInfo; }

  void set_base_grid(const amrex::BoxArray& ba) { activeRegion = ba; }

  bool is_grid_empty() const { return isGridEmpty; }

  // 1. Allocate memory for Fab declared in this class.
  // 2. Set cellStatus and nodeStatus. If cGridsOld is not empty, it will also
  // decide if a cell/node is new or not.
  void distribute_grid_arrays(const amrex::Vector<amrex::BoxArray>& cGridsOld =
                                  amrex::Vector<amrex::BoxArray>()) {

    for (int iLev = 0; iLev <= finest_level; iLev++) {
      distribute_FabArray(cellStatus[iLev], cGrids[iLev], DistributionMap(iLev),
                          1, nGst, false);

      distribute_FabArray(nodeStatus[iLev], nGrids[iLev], DistributionMap(iLev),
                          1, nGst, false);
    }

    update_cell_status(cGridsOld);

    update_node_status(cGridsOld);
  }

  const amrex::iMultiFab& cell_status(int iLev) const {
    return cellStatus[iLev];
  }

  const amrex::iMultiFab& node_status(int iLev) const {
    return nodeStatus[iLev];
  }

  // If cGridsOld is provided, it will also decide if a cell is new or not.
  void update_cell_status(const amrex::Vector<amrex::BoxArray>& cGridsOld =
                              amrex::Vector<amrex::BoxArray>()) {

    for (int iLev = 0; iLev < nLev; iLev++) {
      if (cellStatus[iLev].empty())
        continue;

      // Set default status for all cells.
      for (amrex::MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.fabbox();
        const auto& cellArr = cellStatus[iLev][mfi].array();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);

        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              bit::set_boundary(cellArr(i, j, k));
            }
      }

      // Set 'boundary', 'new' status.
      for (amrex::MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const amrex::Array4<int>& cellArr = cellStatus[iLev][mfi].array();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);

        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              // Not boundary cell
              bit::set_not_boundary(cellArr(i, j, k));

              // New active cell
              bit::set_new(cellArr(i, j, k));

              if (!cGridsOld.empty()) {
                if (cGridsOld[iLev].contains(
                        amrex::IntVect{ AMREX_D_DECL(i, j, k) })) {
                  bit::set_not_new(cellArr(i, j, k));
                }
              }
            }
      }

      // Set the 'refined' status
      if (iLev < max_level) {
        const int iRefined = 1, iNotRefined = 2;
        auto iRefine =
            amrex::makeFineMask(grids[iLev], dmap[iLev], grids[iLev + 1],
                                ref_ratio[iLev], iNotRefined, iRefined);

        for (amrex::MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
          const amrex::Box& box = mfi.validbox();
          const amrex::Array4<int>& cellArr = cellStatus[iLev][mfi].array();
          const auto& iRef = iRefine[mfi].array();
          const auto lo = amrex::lbound(box);
          const auto hi = amrex::ubound(box);

          for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
              for (int i = lo.x; i <= hi.x; ++i) {
                if (iRef(i, j, k) == iRefined) {
                  bit::set_refined(cellArr(i, j, k));
                }
              }
        }
      }

      cellStatus[iLev].FillBoundary(Geom(iLev).periodicity());

      // Set the edge cells.
      // Q: But what is the edge cell?
      // A: It is a physical cell that has one or more neighbor cells are
      // boundary cell.
      for (amrex::MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const amrex::Array4<int>& cellArr = cellStatus[iLev][mfi].array();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);

        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {

              for (int kk = k - 1; kk <= k + 1; kk++)
                for (int jj = j - 1; jj <= j + 1; jj++)
                  for (int ii = i - 1; ii <= i + 1; ii++) {
                    if (bit::is_boundary(cellArr(ii, jj, kk)))
                      bit::set_edge(cellArr(i, j, k));
                  }
            }
      }

      if (isFake2D) {
        // For the fake 2D cases, in the z-direction, only the first layer
        // ghost cells are filled in correctly by the method FillBoundary.
        if (!cellStatus[iLev].empty())
          for (amrex::MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.fabbox();
            const amrex::Array4<int>& cellArr = cellStatus[iLev][mfi].array();
            const auto lo = amrex::lbound(box);
            const auto hi = amrex::ubound(box);

            for (int k = lo.z; k <= hi.z; ++k)
              if (k < -1 || k > 1)
                for (int j = lo.y; j <= hi.y; ++j)
                  for (int i = lo.x; i <= hi.x; ++i) {
                    cellArr(i, j, k) = cellArr(i, j, 0);
                  }
          }
      }
    }
  }

  // If cGridsOld is provided, it will also decide if a node is new or not.
  void update_node_status(const amrex::Vector<amrex::BoxArray>& cGridsOld =
                              amrex::Vector<amrex::BoxArray>()) {
    for (int iLev = 0; iLev < nLev; iLev++) {
      if (nodeStatus[iLev].empty())
        continue;

      // Set default status for all nodes.
      for (amrex::MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.fabbox();
        const auto& nodeArr = nodeStatus[iLev][mfi].array();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);

        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              bit::set_boundary(nodeArr(i, j, k));
              bit::set_skip(nodeArr(i, j, k));
            }
      }

      amrex::BoxArray nodeBAOld;

      if (!cGridsOld.empty()) {
        nodeBAOld =
            convert(cGridsOld[iLev], amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });
      }

      // Set 'boundary', 'new' status.
      for (amrex::MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const auto& nodeArr = nodeStatus[iLev][mfi].array();

        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);

        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              // Not boundary cell
              bit::set_not_boundary(nodeArr(i, j, k));

              // New active cell
              bit::set_new(nodeArr(i, j, k));

              if (!nodeBAOld.empty()) {
                if (nodeBAOld.contains(
                        amrex::IntVect{ AMREX_D_DECL(i, j, k) })) {
                  bit::set_not_new(nodeArr(i, j, k));
                }
              }
            }
      }

      nodeStatus[iLev].FillBoundary(Geom(iLev).periodicity());

      // Set the 'edge' status
      // Q: But what is the edge node?
      // A: It is a node at the boundary of a level.
      for (amrex::MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const amrex::Array4<int>& nodeArr = nodeStatus[iLev][mfi].array();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);

        { // Set 'owner' and 'skip' status
          const auto& cellBox = convert(box, { AMREX_D_DECL(0, 0, 0) });
          const auto& cell = cellStatus[iLev][mfi].array();
          int diMax = 0, diMin = -1;
          int djMax = 0, djMin = -1;
          int dkMax = 0, dkMin = -1;
          if (isFake2D) {
            dkMin = 0;
          }
          // If this box is the owner of this node?
          auto is_the_box_owner = [&](int i, int j, int k) {
            for (int dk = dkMax; dk >= dkMin; dk--)
              for (int dj = djMax; dj >= djMin; dj--)
                for (int di = diMax; di >= diMin; di--) {
                  if (!bit::is_boundary(cell(i + di, j + dj, k + dk))) {
                    // Find the first CELL that shares this node.
                    if (cellBox.contains(amrex::IntVect{
                            AMREX_D_DECL(i + di, j + dj, k + dk) })) {
                      return true;
                    } else {
                      return false;
                    }
                  }
                }
            amrex::Abort("Error: something is wrong here!");
            return false;
          };

          for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
              for (int i = lo.x; i <= hi.x; ++i) {
                if (!isFake2D || k == lo.z) {
                  // for fake 2D , only use the layer of k=0
                  bit::set_not_skip(nodeArr(i, j, k));

                  if (i == lo.x || i == hi.x || j == lo.y || j == hi.y ||
                      (!isFake2D && (k == lo.z || k == hi.z))) {
                    // Block boundary nodes.
                    if (is_the_box_owner(i, j, k)) {
                      bit::set_owner(nodeArr(i, j, k));
                    } else {
                      bit::set_not_owner(nodeArr(i, j, k));
                    }

                  } else {
                    // Nodes indside the box.
                    bit::set_owner(nodeArr(i, j, k));
                  }
                }
              }
        }

        // Set 'edge' status
        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {

              for (int kk = k - 1; kk <= k + 1; kk++)
                for (int jj = j - 1; jj <= j + 1; jj++)
                  for (int ii = i - 1; ii <= i + 1; ii++) {
                    if (bit::is_boundary(nodeArr(ii, jj, kk)))
                      bit::set_edge(nodeArr(i, j, k));
                  }
            }
      }
    }
  }

  std::string lev_string(int iLev) {
    std::string sLev = "_lev_" + std::to_string(iLev);
    if (nLev == 1) {
      // Keep backward compatibility.
      sLev = "";
    }
    return sLev;
  }

  // Find the finest level that contains the cell (x, y, z).
  inline int get_finest_lev(amrex::Real const x, amrex::Real const y,
                            amrex::Real const z) const {
    amrex::Real loc[3] = { x, y, z };
    for (int iLev = finest_level; iLev >= 0; iLev--) {
      auto idx = Geom(iLev).CellIndex(loc);
      if (cGrids[iLev].contains(idx)) {
        return iLev;
      }
    }

    printf("Error: x = %f, y = %f, z = %f is out of the domain\n", loc[0],
           loc[1], loc[2]);
    amrex::Abort();
    return -1; // To suppress compiler warnings.
  }

  inline int find_mpi_rank_from_coord(amrex::Real const x, amrex::Real const y,
                                      amrex::Real const z) const {
    amrex::Real loc[3] = { x, y, z };

    int iLev = get_finest_lev(x, y, z);

    auto idx = Geom(iLev).CellIndex(loc);

    int rank =
        find_mpi_rank_from_cell_index(iLev, idx[ix_], idx[iy_], idx[iz_]);

    return rank;
  }

  inline int find_mpi_rank_from_cell_index(int const iLev, int const i,
                                           int const j, int const k) const {
    amrex::IntVect idx = { AMREX_D_DECL(i, j, k) };
    for (int ii = 0, n = cGrids[iLev].size(); ii < n; ii++) {
      const amrex::Box& bx = cGrids[iLev][ii];
      if (bx.contains(idx))
        return DistributionMap(iLev)[ii];
    }

    amrex::AllPrint() << "iLev = " << iLev << " idx = " << idx << std::endl;
    amrex::Abort("Error: can not find this cell!");
    return -1; // To suppress compiler warnings.
  }

  void calc_node_grids() {
    nGrids.clear();

    if (activeRegion.empty()) {
      // If there is no active cell, still push an empty box array as the base
      // grid
      nGrids.push_back(amrex::BoxArray());
    } else {
      nGrids.resize(finest_level + 1);
      for (int iLev = 0; iLev <= finest_level; iLev++) {
        nGrids[iLev] =
            amrex::convert(cGrids[iLev], amrex::IntVect::TheNodeVector());
      }
    }
  }

  void print_grid_info(bool printBoxes = false) {
    amrex::Print() << printPrefix << " =======Grid Info========" << std::endl;
    amrex::Print() << printPrefix << " nLev = " << nLev << std::endl;
    amrex::Print() << printPrefix << " finest_level = " << finest_level
                   << std::endl;

    for (int iLev = 0; iLev <= finest_level; iLev++) {
      amrex::Print() << printPrefix << " iLev = " << iLev
                     << "\t # of boxes = " << std::setw(9)
                     << cGrids[iLev].size()
                     << "\t # of cells = " << std::setw(11) << CountCells(iLev)
                     << "\t max_grid_size = " << max_grid_size[iLev]
                     << std::endl;
    }

    if (printBoxes) {
      for (int iLev = 0; iLev <= finest_level; iLev++) {
        amrex::Print() << printPrefix << " Boxes of iLev = " << iLev
                       << std::endl;
        for (int ii = 0, n = cGrids[iLev].size(); ii < n; ii++) {
          amrex::Print() << printPrefix << " box " << ii << " = "
                         << cGrids[iLev][ii] << std::endl;
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
      int iLev, amrex::Real time, const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm) override {
    std::string nameFunc = "Grid::MakeNewLevelFromCoarse";
    amrex::Print() << printPrefix << nameFunc << " iLev = " << iLev
                   << std::endl;
  };

  // Remake an existing level using provided BoxArray and DistributionMapping
  // and fill with existing fine and coarse data. overrides the pure virtual
  // function in AmrCore
  virtual void RemakeLevel(int iLev, amrex::Real time,
                           const amrex::BoxArray& ba,
                           const amrex::DistributionMapping& dm) override {
    std::string nameFunc = "Grid::RemakeLevel";
    amrex::Print() << printPrefix << nameFunc << " iLev = " << iLev
                   << std::endl;
  };

  // Delete level data
  // overrides the pure virtual function in AmrCore
  virtual void ClearLevel(int iLev) override {
    std::string nameFunc = "Grid::ClearLevel";
    amrex::Print() << printPrefix << nameFunc << " iLev = " << iLev
                   << std::endl;
  };

  // Make a new level from scratch using provided BoxArray and
  // DistributionMapping. Only used during initialization. overrides the pure
  // virtual function in AmrCore
  virtual void MakeNewLevelFromScratch(
      int iLev, amrex::Real time, const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm) override {
    std::string nameFunc = "Grid::MakeNewLevelFromScratch";
    amrex::Print() << printPrefix << nameFunc << " iLev = " << iLev
                   << std::endl;
  };

  // Tag cells for refinement.
  virtual void ErrorEst(int iLev, amrex::TagBoxArray& tags, amrex::Real time,
                        int ngrow) override {
    std::string nameFunc = "Grid::ErrorEst";
    amrex::Print() << printPrefix << nameFunc << " iLev = " << iLev
                   << std::endl;
#ifdef _AMR_DEV_
    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      const auto tagArr = tags.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
#ifdef _PT_COMPONENT_
            if (i >= 4 && j >= 4 && k >= 4) {
#else
            if (i >= 4 && i < 6 && j >= 4 && j < 6) {
#endif
              tagArr(i, j, k) = amrex::TagBox::SET;
            }
          }
    }
#endif
  };

  virtual void PostProcessBaseGrids(amrex::BoxArray& ba) const override {
    std::string nameFunc = "Grid::PostProcessBaseGrids";
    amrex::Print() << printPrefix << nameFunc << " is called." << std::endl;

    ba = activeRegion;
    const int iLev = 0;
    ba.maxSize(max_grid_size[iLev]);
  };

  void WriteMF(amrex::iMultiFab& MF, int nlev = -1, std::string st = "WriteMF",
               amrex::Vector<std::string> var = {}) {

    amrex::Vector<amrex::iMultiFab> tmf;
    tmf.resize(1);
    tmf[0] = std::move(MF);

    WriteMF(tmf, nlev, st, var);
  }

  void WriteMF(amrex::MultiFab& MF, int nlev = -1, std::string st = "WriteMF",
               amrex::Vector<std::string> var = {}) {

    amrex::Vector<amrex::MultiFab> tmf;
    tmf.resize(1);
    tmf[0] = std::move(MF);

    WriteMF(tmf, nlev, st, var);
  }

  void WriteMF(amrex::Vector<amrex::iMultiFab>& MF, int nlev = -1,
               std::string st = "WriteMF",
               amrex::Vector<std::string> var = {}) {

    amrex::Vector<amrex::MultiFab> tmf;
    tmf.resize(MF.size());
    for (int iLev = 0; iLev < MF.size(); iLev++) {

      tmf[iLev].define(MF[iLev].boxArray(), MF[iLev].DistributionMap(),
                       MF[iLev].nComp(), MF[iLev].nGrow());

      for (amrex::MFIter mfi(MF[iLev]); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.fabbox();
        const amrex::Array4<int>& fab = MF[iLev][mfi].array();
        const amrex::Array4<amrex::Real>& fab2 = tmf[iLev][mfi].array();
        const auto lo = lbound(box);
        const auto hi = ubound(box);

        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              fab2(i, j, k) = fab(i, j, k);
            }
      }
    }
    WriteMF(tmf, nlev, st, var);
  }

  void WriteMF(amrex::Vector<amrex::MultiFab>& MF, int nlev = -1,
               std::string st = "WriteMF",
               amrex::Vector<std::string> var = {}) {
    if (nlev == -1) {
      nlev = MF.size();
    } else {
      nlev = nlev + 1;
    }
    amrex::Vector<const amrex::MultiFab*> tMF;
    for (int i = 0; i < nlev; ++i) {
      tMF.push_back(&MF[i]);
    }
    amrex::Vector<int> tmpVint;
    if (var.empty()) {
      for (int i = 0; i < MF[0].nComp(); i++) {
        var.push_back(std::to_string(i + 1));
      }
    }
    for (int i = 0; i <= nlev; i++) {
      tmpVint.push_back(0);
    }
    amrex::WriteMultiLevelPlotfile(st, nlev, tMF, var, geom, 0.0, tmpVint,
                                   ref_ratio);
  };
};
#endif
