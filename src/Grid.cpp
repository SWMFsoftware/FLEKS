#include "Grid.h"

using namespace amrex;

void Grid::regrid(const BoxArray& region, const Grid* const grid) {
  std::string nameFunc = "Grid::regrid_base";

  if (grid) {
    refineRegions = grid->get_refine_regions();
    SetGridEff(grid->gridEff());
  }

  // Why need 'isGridInitialized'? See the explanation in Domain::regrid().
  if (region == activeRegion && isGridInitialized)
    return;

  pre_regrid();

  cGridsOld = cGrids;

  doNeedFillNewCell = true;

  activeRegion = region;

  domainRange.clear();
  for (int iBox = 0; iBox < activeRegion.size(); iBox++) {
    amrex::RealBox rb(activeRegion[iBox], Geom(0).CellSize(), Geom(0).Offset());
    domainRange.push_back(rb);
  }

  isGridEmpty = activeRegion.empty();

  if (isGridEmpty) {
    cGrids.clear();
    cGrids.push_back(BoxArray());
  } else {
    if (grid) {
      SetFinestLevel(grid->finestLevel());
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        // Q: Why is it required to set distribution map here?
        // A: fi and pic should have the same grids and distribution maps.
        // However, it seems AMReX is too smart that it will try to load balance
        // the box arrays so that the distribution maps can be different even
        // the grid is the same. So we need to set the distribution map here.

        SetBoxArray(iLev, grid->boxArray(iLev));
        SetDistributionMap(iLev, grid->DistributionMap(iLev));
      }
    } else {
      // This method will call MakeNewLevelFromScratch() and
      // PostProcessBaseGrids()
      InitFromScratch(0.0);
    }
  }

  // Print() << "dm = " << DistributionMap(0) << std::endl;

  calc_node_grids();

  print_grid_info();

  // If regrid() is called from from read_restart(), activeRegion is not
  // simplifed. Simplify it here.
  activeRegion = activeRegion.simplified();

  isGridInitialized = true;

  post_regrid();
}

//============================================================================//
void Grid::print_grid_info(bool printBoxes) {
  Print() << printPrefix << " =======Grid Info========" << std::endl;
  Print() << printPrefix << " n_lev_max = " << n_lev_max() << std::endl;
  Print() << printPrefix << " n_lev     = " << n_lev() << std::endl;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    Print() << printPrefix << " iLev = " << iLev
            << "\t # of boxes = " << std::setw(9) << cGrids[iLev].size()
            << "\t # of cells = " << std::setw(11) << CountCells(iLev)
            << "\t max_grid_size = " << max_grid_size[iLev] << std::endl;
  }

  if (printBoxes) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      Print() << printPrefix << " Boxes of iLev = " << iLev << std::endl;
      for (int ii = 0, n = cGrids[iLev].size(); ii < n; ii++) {
        Print() << printPrefix << " box " << ii << " = " << cGrids[iLev][ii]
                << std::endl;
      }
    }
  }

  Print() << printPrefix << " =========================\n" << std::endl;
}

//============================================================================//
void Grid::distribute_grid_arrays(const Vector<BoxArray>& cGridsOld) {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    distribute_FabArray(cellStatus[iLev], cGrids[iLev], DistributionMap(iLev),
                        1, nGst, false);

    distribute_FabArray(nodeStatus[iLev], nGrids[iLev], DistributionMap(iLev),
                        1, nGst, false);
  }

  update_cell_status(cGridsOld);

  update_node_status(cGridsOld);
}

//============================================================================//
void Grid::update_cell_status(const Vector<BoxArray>& cGridsOld) {

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (cellStatus[iLev].empty())
      continue;

    // Set default status for all cells.
    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const auto& cellArr = cellStatus[iLev][mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            bit::set_lev_boundary(cellArr(i, j, k));
          }
    }

    // Set 'boundary', 'new' status.
    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            // Not boundary cell
            bit::set_not_lev_boundary(cellArr(i, j, k));

            // New active cell
            bit::set_new(cellArr(i, j, k));

            if (!cGridsOld.empty()) {
              if (cGridsOld[iLev].contains(IntVect{ AMREX_D_DECL(i, j, k) })) {
                bit::set_not_new(cellArr(i, j, k));
              }
            }
          }
    }

    // Set the 'refined' status
    if (iLev < n_lev() - 1) {
      const int iRefined = 1, iNotRefined = 2;
      auto iRefine = makeFineMask(grids[iLev], dmap[iLev], grids[iLev + 1],
                                  ref_ratio[iLev], iNotRefined, iRefined);

      for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
        const auto& iRef = iRefine[mfi].array();
        const auto lo = lbound(box);
        const auto hi = ubound(box);

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
    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {

            for (int kk = k - 1; kk <= k + 1; kk++)
              for (int jj = j - 1; jj <= j + 1; jj++)
                for (int ii = i - 1; ii <= i + 1; ii++) {
                  if (bit::is_lev_boundary(cellArr(ii, jj, kk)))
                    bit::set_lev_edge(cellArr(i, j, k));
                }
          }
    }

    if (isFake2D) {
      // For the fake 2D cases, in the z-direction, only the first layer
      // ghost cells are filled in correctly by the method FillBoundary.
      if (!cellStatus[iLev].empty())
        for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
          const Box& box = mfi.fabbox();
          const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
          const auto lo = lbound(box);
          const auto hi = ubound(box);

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

//============================================================================//
void Grid::update_node_status(const Vector<BoxArray>& cGridsOld) {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (nodeStatus[iLev].empty())
      continue;

    // Set default status for all nodes.
    for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const auto& nodeArr = nodeStatus[iLev][mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            bit::set_lev_boundary(nodeArr(i, j, k));
          }
    }

    BoxArray nodeBAOld;

    if (!cGridsOld.empty()) {
      nodeBAOld = convert(cGridsOld[iLev], IntVect{ AMREX_D_DECL(1, 1, 1) });
    }

    // Set 'boundary', 'new' status.
    for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& nodeArr = nodeStatus[iLev][mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            // Not boundary cell
            bit::set_not_lev_boundary(nodeArr(i, j, k));

            // New active cell
            bit::set_new(nodeArr(i, j, k));

            if (!nodeBAOld.empty()) {
              if (nodeBAOld.contains(IntVect{ AMREX_D_DECL(i, j, k) })) {
                bit::set_not_new(nodeArr(i, j, k));
              }
            }
          }
    }

    nodeStatus[iLev].FillBoundary(Geom(iLev).periodicity());

    // Set the 'edge' status
    // Q: But what is the edge node?
    // A: It is a node at the boundary of a level.
    for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<int>& nodeArr = nodeStatus[iLev][mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      { // Set 'owner' status
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
                if (!bit::is_lev_boundary(cell(i + di, j + dj, k + dk))) {
                  // Find the first CELL that shares this node.
                  if (cellBox.contains(
                          IntVect{ AMREX_D_DECL(i + di, j + dj, k + dk) })) {
                    return true;
                  } else {
                    return false;
                  }
                }
              }
          Abort("Error: something is wrong here!");
          return false;
        };

        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              if (!isFake2D || k == lo.z) {

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
                  if (bit::is_lev_boundary(nodeArr(ii, jj, kk)))
                    bit::set_lev_edge(nodeArr(i, j, k));
                }
          }
    }
  }
}
