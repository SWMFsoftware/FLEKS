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
#include "Regions.h"
#include "Utility.h"

// This class define the grid information, but NOT the data on the grid.
class Grid : public amrex::AmrCore {

protected:
  int nGst;

  std::string printPrefix;
  std::string gridName;
  int gridID;

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

  amrex::Vector<Regions> refineRegions;

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

    cellStatus.resize(n_lev_max());
    nodeStatus.resize(n_lev_max());
  };

  ~Grid() = default;

  int n_lev() const { return finestLevel() + 1; }

  // n_lev_max() is usually only used for initialization. n_lev() shoudl be used
  // for most purposes.
  int n_lev_max() const { return maxLevel() + 1; }

  int get_n_ghost() const { return nGst; }

  const amrex::AmrInfo& get_amr_info() const { return gridAmrInfo; }

  void set_base_grid(const amrex::BoxArray& ba) { activeRegion = ba; }

  bool is_grid_empty() const { return isGridEmpty; }

  void update_refine_region(const amrex::Vector<Regions>& in) {
    refineRegions = in;
  }

  void print_grid_info(bool printBoxes = false);

  // 1. Allocate memory for Fab declared in this class.
  // 2. Set cellStatus and nodeStatus. If cGridsOld is not empty, it will also
  // decide if a cell/node is new or not.
  void distribute_grid_arrays(const amrex::Vector<amrex::BoxArray>& cGridsOld =
                                  amrex::Vector<amrex::BoxArray>());

  // If cGridsOld is provided, it will also decide if a cell is new or not.
  void update_cell_status(const amrex::Vector<amrex::BoxArray>& cGridsOld =
                              amrex::Vector<amrex::BoxArray>());

  // If cGridsOld is provided, it will also decide if a node is new or not.
  void update_node_status(const amrex::Vector<amrex::BoxArray>& cGridsOld =
                              amrex::Vector<amrex::BoxArray>());

  const amrex::iMultiFab& cell_status(int iLev) const {
    return cellStatus[iLev];
  }

  const amrex::iMultiFab& node_status(int iLev) const {
    return nodeStatus[iLev];
  }

  std::string lev_string(int iLev) {
    std::string sLev = "_lev_" + std::to_string(iLev);
    if (n_lev_max() == 1) {
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

  //===========================================================================
  inline int find_mpi_rank_from_coord(amrex::Real const x, amrex::Real const y,
                                      amrex::Real const z) const {
    amrex::Real loc[3] = { x, y, z };

    int iLev = get_finest_lev(x, y, z);

    auto idx = Geom(iLev).CellIndex(loc);

    int rank =
        find_mpi_rank_from_cell_index(iLev, idx[ix_], idx[iy_], idx[iz_]);

    return rank;
  }

  //===========================================================================
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

  //===========================================================================
  void calc_node_grids() {
    nGrids.clear();

    if (activeRegion.empty()) {
      // If there is no active cell, still push an empty box array as the base
      // grid
      nGrids.push_back(amrex::BoxArray());
    } else {
      nGrids.resize(finest_level + 1);
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        nGrids[iLev] =
            amrex::convert(cGrids[iLev], amrex::IntVect::TheNodeVector());
      }
    }
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
    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.validbox();
      const auto tagArr = tags.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            amrex::Real xyz[nDim];
            Geom(iLev).CellCenter({ AMREX_D_DECL(i, j, k) }, xyz);
            if (refineRegions[iLev].is_inside(xyz)) {
              tagArr(i, j, k) = amrex::TagBox::SET;
            }
          }
    }
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
