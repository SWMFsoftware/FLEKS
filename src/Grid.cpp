#include <cmath>
#include <sstream>

#include <AMReX_PlotFileUtil.H>

#include "Bit.h"
#include "FleksDistributionMap.h"
#include "Grid.h"
#include "GridUtility.h"

using namespace amrex;

bool Grid::set_refine_region(int iLev, const std::string& selector) {
  if (iLev < 0 || iLev >= static_cast<int>(refineRegionsStr.size()) - 1)
    amrex::Abort("Invalid refinement level " + std::to_string(iLev) +
                 " in #REFINEREGION: max allowed level is " +
                 std::to_string(static_cast<int>(refineRegionsStr.size()) - 2));

  std::stringstream input(selector);
  std::string token, normalized;
  bool hasNone = false;
  int tokenCount = 0;
  while (input >> token) {
    ++tokenCount;
    hasNone |= token == "none";
    if (!normalized.empty())
      normalized += ' ';
    normalized += token;
  }
  if (hasNone && tokenCount > 1)
    amrex::Abort("Cannot combine 'none' with other regions in #REFINEREGION");
  if (normalized == "none")
    normalized.clear();
  if (refineRegionsStr[iLev] == normalized)
    return false;
  refineRegionsStr[iLev] = normalized;
  return true;
}

void Grid::replace_refine_region(int iLev, const Regions& region) {
  if (iLev < 0 || iLev >= static_cast<int>(refineRegions.size()) - 1)
    amrex::Abort("Invalid refinement level in dynamic refinement update");
  refineRegions[iLev] = region;
  mark_refinement_modified();
}

void Grid::rebuild_refine_regions(const Vector<std::shared_ptr<Shape> >& shapes,
                                  const std::set<std::string>& shapeNames) {
  for (int i = 0; i < static_cast<int>(refineRegionsStr.size()); ++i) {
    std::stringstream input(refineRegionsStr[i]);
    std::string token;
    while (input >> token) {
      if (token.size() < 2 || (token[0] != '+' && token[0] != '-'))
        amrex::Abort("Invalid shape prefix in #REFINEREGION: '" + token +
                     "' (must start with '+' or '-')");
      if (shapeNames.count(token.substr(1)) == 0)
        amrex::Abort("Unknown shape in #REFINEREGION: '" + token.substr(1) +
                     "'");
    }
    refineRegions[i].define(shapes, refineRegionsStr[i]);
  }
}

Vector<DistributionMapping> Grid::calc_balanced_maps(bool doSplitLevs) {
  BL_PROFILE("calc_balanced_maps");

  Vector<DistributionMapping> dmap(n_lev_max());

  Vector<int> rankStart(n_lev(), 0);
  Vector<int> nProcEachLev(n_lev(), ParallelDescriptor::NProcs());

  if (doSplitLevs) {
    Real totalCost = 0;
    Vector<Real> levCost(n_lev());

    for (int iLev = 0; iLev < n_lev(); iLev++) {
      levCost[iLev] = cellCost[iLev].sum();
      totalCost += levCost[iLev];
    }

    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (iLev == 0) {
        rankStart[iLev] = 0;
      } else {
        rankStart[iLev] = nProcEachLev[iLev - 1] + rankStart[iLev - 1];
      }

      if (iLev < n_lev() - 1) {
        nProcEachLev[iLev] =
            floor(ParallelDescriptor::NProcs() * levCost[iLev] / totalCost);
      } else {
        nProcEachLev[iLev] = ParallelDescriptor::NProcs() - rankStart[iLev];
      }
      Print() << printPrefix << " Ranks from " << rankStart[iLev] << " to "
              << rankStart[iLev] + nProcEachLev[iLev] - 1
              << " are assigned to iLev = " << iLev << std::endl;
    }
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    Vector<int> ord(ParallelDescriptor::NProcs());
    for (int i = 0; i < nProcEachLev[iLev]; ++i) {
      ord[i] = i + rankStart[iLev];
    }

    Real eff;
    dmap[iLev] = FleksDistributionMap::make_balanced_map(
        BalanceMethod::SFC, cellCost[iLev], nProcEachLev[iLev], ord, eff);
  }

  return dmap;
}

//==========================================================
void Grid::load_balance(const Grid* other, bool doSplitLevs) {

  if (other) {
    regrid(other, true);
  } else {
    Vector<DistributionMapping> dmap = calc_balanced_maps(doSplitLevs);

    cGridsOld = cGrids;
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      SetDistributionMap(iLev, dmap[iLev]);
    }

    calc_node_grids();
    print_grid_info();
    isNewGrid = false;
    post_regrid();
  }
}

//==========================================================
void Grid::regrid(const BoxArray& region) {
  std::string nameFunc = "Grid::regrid_base";
  BL_PROFILE(nameFunc);

  // The base boxes can be identical while new refinement tags or solver regions
  // require a different hierarchy. isNewGrid also forces first initialization.
  if (region == activeRegion && !isNewGrid && !refineRegions.is_modified() &&
      !solverRegionsModified)
    return;

  cGridsOld = cGrids;
  activeRegion = region;

  isGridEmpty = activeRegion.empty();

  if (isGridEmpty) {
    cGrids.clear();
    cGrids.push_back(BoxArray());
  } else {
    // This method will call MakeNewLevelFromScratch() and
    // PostProcessBaseGrids()
    InitFromScratch(0.0);
  }

  finish_regrid();
}

void Grid::regrid(const Grid* const grid, bool doLoadBalance) {
  if (!grid)
    return;

  std::string nameFunc = "Grid::regrid_base";
  BL_PROFILE(nameFunc);

  if (!doLoadBalance) {
    refineRegionsStr = grid->refineRegionsStr;
    refineRegions = grid->refineRegions;
    defaultSolverMode = grid->defaultSolverMode;
    regionSolvers = grid->regionSolvers;
    SetGridEff(grid->gridEff());

    const BoxArray& region = grid->get_base_grid();
    if (region == activeRegion && !isNewGrid && !refineRegions.is_modified() &&
        !solverRegionsModified)
      return;

    cGridsOld = cGrids;
    activeRegion = region;
  } else {
    // Load balancing changes ownership, not the cells in the hierarchy.
    // Keep the old boxes so retained nodes are not marked as new.
    cGridsOld = cGrids;
  }

  isGridEmpty = activeRegion.empty();

  if (isGridEmpty) {
    cGrids.clear();
    cGrids.push_back(BoxArray());
  } else {
    // Q: Why is it required to set distribution map (inside init_grid())?
    // A: fi and pic should have the same grids and distribution maps.
    // However, it seems AMReX is too smart that it will try to load balance
    // the box arrays so that the distribution maps can be different even
    // the grid is the same. So we need to set the distribution map here.
    set_ba_and_dm(grid);
  }

  finish_regrid();
}

void Grid::restore_grid(const Vector<BoxArray>& boxes) {
  if (boxes.empty() || static_cast<int>(boxes.size()) > n_lev_max())
    Abort("Invalid restart grid level count.");

  cGridsOld = cGrids;
  cGridsOld.resize(n_lev_max());
  activeRegion = boxes[0];
  isGridEmpty = activeRegion.empty();
  SetFinestLevel(static_cast<int>(boxes.size()) - 1);
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    SetBoxArray(iLev, boxes[iLev]);
    SetDistributionMap(iLev, DistributionMapping(boxes[iLev]));
  }
  for (int iLev = n_lev(); iLev < n_lev_max(); ++iLev) {
    ClearBoxArray(iLev);
    ClearDistributionMap(iLev);
  }
  finish_regrid();
  clear_refinement_modified();
}

void Grid::finish_regrid() {

  // InitFromScratch lowers finest_level without clearing removed levels.
  // Do not retain their boxes: a later refinement would otherwise classify
  // recreated fine nodes as old and skip their fluid initialization.
  for (int iLev = n_lev();
       iLev < n_lev_max() && iLev < static_cast<int>(cGrids.size()); ++iLev) {
    ClearBoxArray(iLev);
    ClearDistributionMap(iLev);
  }

  // Split refined levels (iLev >= 1) along solver region boundaries so every
  // fine Box also belongs 100% to a single PicSolverMode.
  if (!regionSolvers.empty() && !isGridEmpty) {
    for (int iLev = 1; iLev < n_lev(); ++iLev) {
      if (cGrids[iLev].empty())
        continue;
      BoxArray splitBA = split_boxes_by_solver_regions(cGrids[iLev], iLev);
      if (splitBA != cGrids[iLev]) {
        SetBoxArray(iLev, splitBA);
        SetDistributionMap(iLev, DistributionMapping(splitBA));
      }
    }
  }

  calc_node_grids();

  print_grid_info();

  // If regrid() is called from from read_restart(), activeRegion is not
  // simplifed. Simplify it here.
  activeRegion = activeRegion.simplified();

  activeDomainRanges.clear();
  for (int iBox = 0; iBox < activeRegion.size(); iBox++) {
    RealBox rb(activeRegion[iBox], Geom(0).CellSize(), Geom(0).Offset());
    activeDomainRanges.push_back(rb);
  }

  isNewGrid = false;

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
      for (int ii = 0, n = cGrids[iLev].size(); ii < n; ++ii) {
        Print() << printPrefix << " box " << ii << " = " << cGrids[iLev][ii]
                << std::endl;
      }
    }
  }

  Print() << printPrefix << " =========================\n" << std::endl;
}

//============================================================================//
void Grid::distribute_grid_arrays() { distribute_grid_arrays(cGridsOld); }

//============================================================================//
void Grid::distribute_grid_arrays(const Vector<BoxArray>& cGridsOld) {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    distribute_FabArray(cellStatus[iLev], cGrids[iLev], DistributionMap(iLev),
                        1, nGst, false);

    distribute_FabArray(nodeStatus[iLev], nGrids[iLev], DistributionMap(iLev),
                        1, nGst, false);

    distribute_FabArray(cellCost[iLev], cGrids[iLev], DistributionMap(iLev), 1,
                        0, false);

    if (use_body() && bodyParticleBC == ParticleBC::ambipolar) {
      distribute_FabArray(bodyNetCharge[iLev], cGrids[iLev],
                          DistributionMap(iLev), 1, nGst, true, 0.0);
      distribute_FabArray(bodyNetChargeInc[iLev], cGrids[iLev],
                          DistributionMap(iLev), 1, nGst, false);
      bodyNetChargeInc[iLev].setVal(0.0);
    }
  }

  update_grid_status(cGridsOld);
  update_solver_masks();
}

//============================================================================//
void Grid::update_cell_status(const Vector<BoxArray>& cGridsOld) {
  std::string nameFunc = "Grid::update_cell_status";
  BL_PROFILE(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (cellStatus[iLev].empty())
      continue;

    // Set default status for all cells.
    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const auto& cellArr = cellStatus[iLev][mfi].array();
      ParallelFor(box, [&](int i, int j, int k) noexcept {
        bit::set_lev_boundary(cellArr(i, j, k));
        bit::set_not_domain_boundary(cellArr(i, j, k));
      });
    }
    // Set 'boundary', 'new' status.
    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
      ParallelFor(box, [&](int i, int j, int k) noexcept {
        // Not boundary cell
        bit::set_not_lev_boundary(cellArr(i, j, k));

        // New active cell
        bit::set_new(cellArr(i, j, k));
      });

      if (!cGridsOld.empty()) {
        for (int b = 0, nb = cGridsOld[iLev].size(); b < nb; ++b) {
          const Box isect = box & cGridsOld[iLev][b];
          if (isect.ok()) {
            ParallelFor(isect, [&](int i, int j, int k) noexcept {
              bit::set_not_new(cellArr(i, j, k));
            });
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
        ParallelFor(box, [&](int i, int j, int k) noexcept {
          if (iRef(i, j, k) == iRefined) {
            bit::set_refined(cellArr(i, j, k));
          }
        });
      }
    }

    cellStatus[iLev].FillBoundary(Geom(iLev).periodicity());

    // Find domain boundary cells: level boundary cells whose centers are
    // outside the active region. A cell center is never on a level-0 cell
    // face, so it is inside a box of the active region if and only if its
    // level-0 cell is. Testing that in index space, against only the active
    // boxes that overlap this fab, is equivalent to is_inside_domain(), which
    // loops over all boxes for every cell and is very slow for many boxes.
    IntVect ratioToLev0(1);
    for (int lev = 0; lev < iLev; ++lev)
      ratioToLev0 *= ref_ratio[lev];

    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
      // Find every active-region box that overlaps the coarsened version of
      // this fab.
      const auto isects =
          activeRegion.intersections(amrex::coarsen(box, ratioToLev0));
      ParallelFor(box, [&](int i, int j, int k) noexcept {
        if (bit::is_lev_boundary(cellArr(i, j, k))) {
          const IntVect iv0 =
              amrex::coarsen(IntVect(AMREX_D_DECL(i, j, k)), ratioToLev0);
          bool isInside = false;
          for (const auto& isect : isects) {
            if (isect.second.contains(iv0)) {
              isInside = true;
              break;
            }
          }
          if (!isInside) {
            bit::set_domain_boundary(cellArr(i, j, k));
          }
        }
      });
    }

    // Mark the cells that belong to the absorbing inner body (see #BODY).
    // The body is approximated by the union of the cells whose centers are
    // inside the sphere, i.e., a staircase boundary with an error of dx/2.
    if (useBody) {
      for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.fabbox();
        const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
        ParallelFor(box, [&](int i, int j, int k) noexcept {
          Real xyz[nDim];
          Geom(iLev).CellCenter({ AMREX_D_DECL(i, j, k) }, xyz);
          if (is_inside_body(xyz)) {
            bit::set_body(cellArr(i, j, k));
          } else {
            bit::set_not_body(cellArr(i, j, k));
          }
        });
      }

      // Split the body into a one-cell-thick surface layer, where the field
      // boundary condition acts, and the interior, where the fields are frozen
      // (E = 0 and B at its initial value).
      for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
        ParallelFor(box, [&](int i, int j, int k) noexcept {
          if (!bit::is_body(cellArr(i, j, k))) {
            bit::set_not_body_interior(cellArr(i, j, k));
            return;
          }

          bool isInterior = bit::is_body(cellArr(i - 1, j, k)) &&
                            bit::is_body(cellArr(i + 1, j, k));
          if (nDim > 1)
            isInterior = isInterior && bit::is_body(cellArr(i, j - 1, k)) &&
                         bit::is_body(cellArr(i, j + 1, k));
          if (nDim > 2)
            isInterior = isInterior && bit::is_body(cellArr(i, j, k - 1)) &&
                         bit::is_body(cellArr(i, j, k + 1));

          if (isInterior) {
            bit::set_body_interior(cellArr(i, j, k));
          } else {
            bit::set_not_body_interior(cellArr(i, j, k));
          }
        });
      }
      cellStatus[iLev].FillBoundary(Geom(iLev).periodicity());
    }

    // Set edge cells and find cells with 'is_refined' neighbors
    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
      ParallelFor(box, [&](int i, int j, int k) noexcept {
        const int kmin = nDim > 2 ? k - 1 : k;
        const int kmax = nDim > 2 ? k + 1 : k;
        const bool notRefined = !bit::is_refined(cellArr(i, j, k));

        for (int kk = kmin; kk <= kmax; ++kk) {
          for (int jj = j - 1; jj <= j + 1; ++jj) {
            for (int ii = i - 1; ii <= i + 1; ++ii) {
              const int neighbor = cellArr(ii, jj, kk);
              if (bit::is_lev_boundary(neighbor)) {
                bit::set_lev_edge(cellArr(i, j, k));

                if (bit::is_domain_boundary(neighbor)) {
                  bit::set_domain_edge(cellArr(i, j, k));
                }
              }

              if (notRefined && bit::is_refined(neighbor)) {
                bit::set_refined_neighbour(cellArr(i, j, k));
              }
            }
          }
        }
      });
    }

    if (isFake2D) {
      // For the fake 2D cases, in the z-direction, only the first layer
      // ghost cells are filled in correctly by the method FillBoundary.
      if (!cellStatus[iLev].empty())
        for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
          const Box& box = mfi.fabbox();
          const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
          ParallelFor(box, [&](int i, int j, int k) noexcept {
            if (k < -1 || k > 1)
              cellArr(i, j, k) = cellArr(i, j, 0);
          });
        }
    }
  }
}

//============================================================================//
void Grid::update_node_status(const Vector<BoxArray>& cGridsOld) {
  std::string nameFunc = "Grid::update_node_status()";
  BL_PROFILE(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (nodeStatus[iLev].empty())
      continue;

    // Set default status for all nodes.
    for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const auto& nodeArr = nodeStatus[iLev][mfi].array();
      ParallelFor(box, [&](int i, int j, int k) noexcept {
        bit::set_lev_boundary(nodeArr(i, j, k));
        bit::set_not_domain_boundary(nodeArr(i, j, k));
        bit::set_not_refined(nodeArr(i, j, k));
      });
    }

    BoxArray nodeBAOld;

    if (!cGridsOld.empty()) {
      nodeBAOld = convert(cGridsOld[iLev], IntVect(1));
    }

    // Set 'boundary', 'new' status.
    for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& nodeArr = nodeStatus[iLev][mfi].array();
      ParallelFor(box, [&](int i, int j, int k) noexcept {
        // Not boundary cell
        bit::set_not_lev_boundary(nodeArr(i, j, k));

        // New active cell
        bit::set_new(nodeArr(i, j, k));
      });

      if (!nodeBAOld.empty()) {
        for (int b = 0, nb = nodeBAOld.size(); b < nb; ++b) {
          const Box isect = box & nodeBAOld[b];
          if (isect.ok()) {
            ParallelFor(isect, [&](int i, int j, int k) noexcept {
              bit::set_not_new(nodeArr(i, j, k));
            });
          }
        }
      }
    }

    nodeStatus[iLev].FillBoundary(Geom(iLev).periodicity());

    // Find domain boundary cells
    for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const Array4<int>& nodeArr = nodeStatus[iLev][mfi].array();

      ParallelFor(box, [&](int i, int j, int k) noexcept {
        if (bit::is_lev_boundary(nodeArr(i, j, k))) {
          Real xyz[nDim];
          Geom(iLev).LoNode({ AMREX_D_DECL(i, j, k) }, xyz);
          if (!is_inside_domain(xyz)) {
            bit::set_domain_boundary(nodeArr(i, j, k));
          }
        }
      });
    }

    // Mark the nodes inside the absorbing inner body (see #BODY). The
    // electric field is pinned to zero on these nodes.
    if (useBody) {
      for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.fabbox();
        const Array4<int>& nodeArr = nodeStatus[iLev][mfi].array();
        ParallelFor(box, [&](int i, int j, int k) noexcept {
          Real xyz[nDim];
          Geom(iLev).LoNode({ AMREX_D_DECL(i, j, k) }, xyz);
          if (is_inside_body(xyz)) {
            bit::set_body(nodeArr(i, j, k));
          } else {
            bit::set_not_body(nodeArr(i, j, k));
          }
        });
      }

      // A body node belongs to the frozen interior when all the cells around
      // it belong to the body interior, so that it only sees frozen cells.
      for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<int>& nodeArr = nodeStatus[iLev][mfi].array();
        const Array4<int const>& cellArr = cellStatus[iLev][mfi].array();
        ParallelFor(box, [&](int i, int j, int k) noexcept {
          if (!bit::is_body(nodeArr(i, j, k))) {
            bit::set_not_body_interior(nodeArr(i, j, k));
            return;
          }

          bool isInterior = true;
          for (int kk = (nDim > 2 ? -1 : 0); kk <= 0; ++kk) {
            for (int jj = (nDim > 1 ? -1 : 0); jj <= 0; ++jj) {
              for (int ii = -1; ii <= 0; ++ii) {
                if (!bit::is_body_interior(cellArr(i + ii, j + jj, k + kk))) {
                  isInterior = false;
                }
              }
            }
          }

          if (isInterior) {
            bit::set_body_interior(nodeArr(i, j, k));
          } else {
            bit::set_not_body_interior(nodeArr(i, j, k));
          }
        });
      }
      nodeStatus[iLev].FillBoundary(Geom(iLev).periodicity());
    }

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
        const int activeDim = get_dim();
        if (activeDim == 2) {
          dkMin = 0;
        }
        // Is the box the owner of this node?
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

        ParallelFor(box, [&](int i, int j, int k) noexcept {
          if (!isFake2D || k == lo.z) {
            if (i == lo.x || i == hi.x || j == lo.y || j == hi.y ||
                (activeDim == 3 && (k == lo.z || k == hi.z))) {
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
        });
      }

      // Set the 'edge' status
      // Q: But what is the edge node?
      // A: It is a node at the boundary of a level.
      ParallelFor(box, [&](int i, int j, int k) noexcept {
        const int kmin = nDim > 2 ? k - 1 : k;
        const int kmax = nDim > 2 ? k + 1 : k;

        for (int kk = kmin; kk <= kmax; ++kk) {
          for (int jj = j - 1; jj <= j + 1; ++jj) {
            for (int ii = i - 1; ii <= i + 1; ++ii) {
              if (bit::is_lev_boundary(nodeArr(ii, jj, kk))) {
                bit::set_lev_edge(nodeArr(i, j, k));

                if (bit::is_domain_boundary(nodeArr(ii, jj, kk))) {
                  bit::set_domain_edge(nodeArr(i, j, k));
                }
              }
            }
          }
        }
      });

      // Set the 'refined' status for nodes
      const auto& cell = cellStatus[iLev][mfi].array();
      ParallelFor(box, [&](int i, int j, int k) noexcept {
        const int kmin = nDim > 2 ? k - 1 : k;
        for (int kk = kmin; kk <= k; ++kk) {
          for (int jj = j - 1; jj <= j; ++jj) {
            for (int ii = i - 1; ii <= i; ++ii) {
              if (bit::is_refined(cell(ii, jj, kk))) {
                bit::set_refined(nodeArr(i, j, k));
              }
            }
          }
        }
      });
    }
  }
}

void Grid::write_mf(Vector<MultiFab>& mf, int nlev, std::string st,
                    Vector<std::string> var) {
  if (nlev == -1) {
    nlev = finest_level + 1;
  } else {
    nlev = nlev + 1;
  }
  Vector<const MultiFab*> tMF;
  tMF.reserve(nlev);
  for (int i = 0; i < nlev; ++i) {
    tMF.push_back(&mf[i]);
  }
  if (var.empty()) {
    var.reserve(mf[0].nComp());
    for (int i = 0; i < mf[0].nComp(); ++i) {
      var.push_back(std::to_string(i + 1));
    }
  }
  Vector<int> tmpVint(nlev + 1, 0);
  WriteMultiLevelPlotfile(st, nlev, tMF, var, geom, 0.0, tmpVint, ref_ratio);
}

//============================================================================//
void Grid::rebuild_solver_regions(const Vector<std::shared_ptr<Shape> >& shapes,
                                  const std::set<std::string>& shapeNames) {
  for (auto& spec : regionSolvers) {
    std::stringstream input(spec.regionStr);
    std::string token;
    bool hasInclude = false;
    while (input >> token) {
      if (token.size() < 2 || (token[0] != '+' && token[0] != '-'))
        Abort("Invalid shape prefix in #REGIONSOLVER: '" + token +
              "' (must start with '+' or '-')");
      if (token[0] == '+')
        hasInclude = true;
      if (shapeNames.count(token.substr(1)) == 0)
        Abort("Unknown shape in #REGIONSOLVER: '" + token.substr(1) + "'");
    }
    if (!hasInclude)
      Abort("Invalid #REGIONSOLVER region '" + spec.regionStr +
            "': must specify at least one '+shape'.");
    spec.region.define(shapes, spec.regionStr);
  }
}

PicSolverMode Grid::eval_region_solver(const Real* xyz) const {
  PicSolverMode mode = defaultSolverMode;
  for (const auto& spec : regionSolvers) {
    if (spec.region.is_inside(xyz)) {
      mode = spec.solver;
    }
  }
  return mode;
}

bool Grid::has_configured_solver(PicSolverMode mode) const {
  if (defaultSolverMode == mode)
    return true;
  for (const auto& spec : regionSolvers) {
    if (spec.solver == mode)
      return true;
  }
  return false;
}

bool Grid::has_active_solver(PicSolverMode mode) const {
  if (cGrids.empty() || cGrids[0].empty() || boxSolver.empty() ||
      boxSolver[0].empty()) {
    return has_configured_solver(mode);
  }
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    if (iLev >= static_cast<int>(boxSolver.size()))
      break;
    for (PicSolverMode m : boxSolver[iLev]) {
      if (m == mode)
        return true;
    }
  }
  return false;
}

bool Grid::is_multi_solver_grid() const {
  int nActive = 0;
  for (PicSolverMode m : { PicSolverMode::Implicit, PicSolverMode::Explicit,
                           PicSolverMode::Hybrid, PicSolverMode::None }) {
    if (has_active_solver(m))
      ++nActive;
  }
  return nActive >= 2;
}

BoxArray Grid::split_boxes_by_solver_regions(const BoxArray& inputBA,
                                             int iLev) const {
  const int maxGridLev =
      std::min(iLev, static_cast<int>(max_grid_size.size()) - 1);
  if (regionSolvers.empty() || inputBA.empty()) {
    BoxArray ba = inputBA;
    if (!ba.empty())
      ba.maxSize(max_grid_size[maxGridLev]);
    return ba;
  }

  const int bfLev =
      std::min(iLev, static_cast<int>(blocking_factor.size()) - 1);
  IntVect bf = blocking_factor[bfLev];
  if (iLev > 0 && (iLev - 1) < static_cast<int>(ref_ratio.size())) {
    bf = amrex::max(bf, ref_ratio[iLev - 1]);
  }
  if (isFake2D && nDim > 2) {
    bf[iz_] = 1;
  }
  for (int d = 0; d < nDim; ++d) {
    IntVect test(1);
    test[d] = bf[d];
    while (test[d] > 1 && !inputBA.coarsenable(test)) {
      test[d] /= 2;
    }
    bf[d] = test[d];
  }

  BoxArray coarseBA = inputBA;
  coarseBA.coarsen(bf);

  const Real* dx = Geom(iLev).CellSize();
  const Real* probLo = Geom(iLev).ProbLo();
  Real dxCoarse[3] = { 1.0, 1.0, 1.0 };
  for (int d = 0; d < nDim; ++d) {
    dxCoarse[d] = dx[d] * bf[d];
  }

  BoxList bl;
  for (int iBox = 0; iBox < coarseBA.size(); ++iBox) {
    const Box& cbx = coarseBA[iBox];
    const IntVect lo = cbx.smallEnd();
    const IntVect hi = cbx.bigEnd();
    const int nx = hi[0] - lo[0] + 1;
    const int ny = (nDim > 1) ? (hi[1] - lo[1] + 1) : 1;
    const int nz = (nDim > 2) ? (hi[2] - lo[2] + 1) : 1;
    const int jLo = (nDim > 1) ? lo[1] : 0;
    const int kLo = (nDim > 2) ? lo[2] : 0;

    std::vector<int> modes(static_cast<std::size_t>(nx) * ny * nz);
    auto idx3 = [nx, ny](int ii, int jj, int kk) -> std::size_t {
      return (static_cast<std::size_t>(kk) * ny + jj) * nx + ii;
    };

    bool uniform = true;
    int firstMode = -1;
    for (int kk = 0; kk < nz; ++kk) {
      for (int jj = 0; jj < ny; ++jj) {
        for (int ii = 0; ii < nx; ++ii) {
          Real xyz[3] = { 0.0, 0.0, 0.0 };
          xyz[0] = probLo[0] + (lo[0] + ii + Real(0.5)) * dxCoarse[0];
          if (nDim > 1)
            xyz[1] = probLo[1] + (jLo + jj + Real(0.5)) * dxCoarse[1];
          if (nDim > 2)
            xyz[2] = probLo[2] + (kLo + kk + Real(0.5)) * dxCoarse[2];
          const int m = static_cast<int>(eval_region_solver(xyz));
          modes[idx3(ii, jj, kk)] = m;
          if (firstMode < 0)
            firstMode = m;
          else if (m != firstMode)
            uniform = false;
        }
      }
    }

    if (uniform) {
      Box fineBx = cbx;
      fineBx.refine(bf);
      bl.push_back(fineBx);
      continue;
    }

    // Partition non-uniform coarse box into maximal contiguous rectangular
    // super-boxes of homogeneous PicSolverMode, then refine by bf.
    std::vector<unsigned char> visited(modes.size(), 0);
    for (int kk = 0; kk < nz; ++kk) {
      for (int jj = 0; jj < ny; ++jj) {
        for (int ii = 0; ii < nx; ++ii) {
          if (visited[idx3(ii, jj, kk)])
            continue;
          const int target = modes[idx3(ii, jj, kk)];

          int iEnd = ii;
          while (iEnd + 1 < nx && !visited[idx3(iEnd + 1, jj, kk)] &&
                 modes[idx3(iEnd + 1, jj, kk)] == target) {
            ++iEnd;
          }

          int jEnd = jj;
          while (jEnd + 1 < ny) {
            bool rowOk = true;
            for (int x = ii; x <= iEnd; ++x) {
              if (visited[idx3(x, jEnd + 1, kk)] ||
                  modes[idx3(x, jEnd + 1, kk)] != target) {
                rowOk = false;
                break;
              }
            }
            if (!rowOk)
              break;
            ++jEnd;
          }

          int kEnd = kk;
          while (kEnd + 1 < nz) {
            bool planeOk = true;
            for (int y = jj; y <= jEnd && planeOk; ++y) {
              for (int x = ii; x <= iEnd; ++x) {
                if (visited[idx3(x, y, kEnd + 1)] ||
                    modes[idx3(x, y, kEnd + 1)] != target) {
                  planeOk = false;
                  break;
                }
              }
            }
            if (!planeOk)
              break;
            ++kEnd;
          }

          for (int z = kk; z <= kEnd; ++z)
            for (int y = jj; y <= jEnd; ++y)
              for (int x = ii; x <= iEnd; ++x)
                visited[idx3(x, y, z)] = 1;

          IntVect sLo(AMREX_D_DECL(lo[0] + ii, jLo + jj, kLo + kk));
          IntVect sHi(AMREX_D_DECL(lo[0] + iEnd, jLo + jEnd, kLo + kEnd));
          Box superBx(sLo, sHi);
          superBx.refine(bf);
          bl.push_back(superBx);
        }
      }
    }
  }

  BoxArray splitBA(std::move(bl));
  splitBA.maxSize(max_grid_size[maxGridLev]);
  return splitBA;
}

void Grid::update_solver_masks() {
  std::string nameFunc = "Grid::update_solver_masks";
  BL_PROFILE(nameFunc);

  const int maxLev = static_cast<int>(boxSolver.size());
  for (int iLev = n_lev(); iLev < maxLev; ++iLev) {
    boxSolver[iLev].clear();
    cellSolver[iLev].clear();
    hybridCollarMask[iLev].clear();
    nodeSolverStatus[iLev].clear();
    nodeBlendWeight[iLev].clear();
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    boxSolver[iLev].clear();
    if (cGrids[iLev].empty()) {
      cellSolver[iLev].clear();
      hybridCollarMask[iLev].clear();
      nodeSolverStatus[iLev].clear();
      nodeBlendWeight[iLev].clear();
      continue;
    }

    const int nBoxes = cGrids[iLev].size();
    boxSolver[iLev].resize(nBoxes, defaultSolverMode);
    if (!regionSolvers.empty()) {
      const Real* dx = Geom(iLev).CellSize();
      const Real* probLo = Geom(iLev).ProbLo();
      for (int iBox = 0; iBox < nBoxes; ++iBox) {
        const Box& bx = cGrids[iLev][iBox];
        Real center[3] = { 0.0, 0.0, 0.0 };
        for (int d = 0; d < nDim; ++d) {
          center[d] = probLo[d] +
                      Real(0.5) * (bx.smallEnd(d) + bx.bigEnd(d) + 1) * dx[d];
        }
        boxSolver[iLev][iBox] = eval_region_solver(center);
      }
    }

    distribute_FabArray(cellSolver[iLev], cGrids[iLev], DistributionMap(iLev),
                        1, nGst, false);
    distribute_FabArray(hybridCollarMask[iLev], cGrids[iLev],
                        DistributionMap(iLev), 1, nGst, false);
    distribute_FabArray(nodeSolverStatus[iLev], nGrids[iLev],
                        DistributionMap(iLev), 1, nGst, false);

    // Fill cellSolver valid cells from owning boxSolver, and initialize ghosts.
    for (MFIter mfi(cellSolver[iLev]); mfi.isValid(); ++mfi) {
      const Box& fbox = mfi.fabbox();
      const Box& vbox = mfi.validbox();
      const auto csArr = cellSolver[iLev][mfi].array();
      const auto cStat = cellStatus[iLev][mfi].array();
      const int boxMode = static_cast<int>(boxSolver[iLev][mfi.index()]);
      if (regionSolvers.empty()) {
        ParallelFor(fbox, [&](int i, int j, int k)
                              noexcept { csArr(i, j, k) = boxMode; });
      } else {
        ParallelFor(fbox, [&](int i, int j, int k) noexcept {
          if (vbox.contains(IntVect(AMREX_D_DECL(i, j, k)))) {
            csArr(i, j, k) = boxMode;
          } else if (bit::is_domain_boundary(cStat(i, j, k))) {
            csArr(i, j, k) = boxMode;
          } else {
            Real xyz[3] = { 0.0, 0.0, 0.0 };
            Geom(iLev).CellCenter({ AMREX_D_DECL(i, j, k) }, xyz);
            csArr(i, j, k) = static_cast<int>(eval_region_solver(xyz));
          }
        });
      }
    }
    cellSolver[iLev].FillBoundary(Geom(iLev).periodicity());
    if (isFake2D && !cellSolver[iLev].empty()) {
      for (MFIter mfi(cellSolver[iLev]); mfi.isValid(); ++mfi) {
        const Box& fbox = mfi.fabbox();
        const auto csArr = cellSolver[iLev][mfi].array();
        ParallelFor(fbox, [&](int i, int j, int k) noexcept {
          if (k != 0)
            csArr(i, j, k) = csArr(i, j, 0);
        });
      }
    }

    // Build nodeSolverStatus by inspecting active adjacent cells around each
    // node (valid up to nGst - 1 ghosts via cellSolver ghosts).
    const int activeDim = get_dim();
    for (MFIter mfi(nodeSolverStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& fbox = mfi.fabbox();
      const Box& vbox = mfi.validbox();
      const auto nsArr = nodeSolverStatus[iLev][mfi].array();
      const auto csArr = cellSolver[iLev][mfi].array();
      const auto cStat = cellStatus[iLev][mfi].array();
      const auto csBox = cellSolver[iLev][mfi].box();
      const Box safeBox = Box(vbox).grow(std::max(0, nGst - 1));
      const int kMinOff = (activeDim > 2) ? -1 : 0;
      const int jMinOff = (activeDim > 1) ? -1 : 0;

      ParallelFor(fbox, [&](int i, int j, int k) noexcept {
        if (!safeBox.contains(IntVect(AMREX_D_DECL(i, j, k)))) {
          nsArr(i, j, k) = 0;
          return;
        }
        int bits = 0;
        for (int kk = kMinOff; kk <= 0; ++kk) {
          for (int jj = jMinOff; jj <= 0; ++jj) {
            for (int ii = -1; ii <= 0; ++ii) {
              const IntVect cIdx(AMREX_D_DECL(i + ii, j + jj, k + kk));
              if (!csBox.contains(cIdx))
                continue;
              if (bit::is_domain_boundary(cStat(cIdx)))
                continue;
              const int m = csArr(cIdx);
              if (m == static_cast<int>(PicSolverMode::Implicit))
                bits |= solver_bit::kTouchImplicit;
              else if (m == static_cast<int>(PicSolverMode::Explicit))
                bits |= solver_bit::kTouchExplicit;
              else if (m == static_cast<int>(PicSolverMode::Hybrid))
                bits |= solver_bit::kTouchHybrid;
              else
                bits |= solver_bit::kTouchNone;
            }
          }
        }
        if ((bits & solver_bit::kTouchImplicit) &&
            (bits & (solver_bit::kTouchExplicit | solver_bit::kTouchHybrid))) {
          bits |= solver_bit::kImplicitSeamDirichlet;
        }
        nsArr(i, j, k) = bits;
      });
    }
    nodeSolverStatus[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  // Second pass over levels now that boxSolver is populated across all levels:
  // build hybridCollarMask and optional nodeBlendWeight.
  const bool multiSolver = is_multi_solver_grid();
  const int blendWidth = max_blend_cells();
  const int collarRadius = std::max(nGst, blendWidth + 2);
  const int activeDim = get_dim();

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    if (cGrids[iLev].empty())
      continue;

    if (!multiSolver || !has_active_solver(PicSolverMode::Hybrid)) {
      const int val = has_active_solver(PicSolverMode::Hybrid) ? 2 : 0;
      hybridCollarMask[iLev].setVal(val);
    } else {
      iMultiFab isFullPic(cGrids[iLev], DistributionMap(iLev), 1, collarRadius);
      isFullPic.setVal(0);
      for (MFIter mfi(isFullPic); mfi.isValid(); ++mfi) {
        const Box& vbox = mfi.validbox();
        const auto fpArr = isFullPic[mfi].array();
        const auto csArr = cellSolver[iLev][mfi].array();
        ParallelFor(vbox, [&](int i, int j, int k) noexcept {
          const int m = csArr(i, j, k);
          fpArr(i, j, k) = (m == static_cast<int>(PicSolverMode::Implicit) ||
                            m == static_cast<int>(PicSolverMode::Explicit))
                               ? 1
                               : 0;
        });
      }
      isFullPic.FillBoundary(Geom(iLev).periodicity());

      for (MFIter mfi(hybridCollarMask[iLev]); mfi.isValid(); ++mfi) {
        const Box& vbox = mfi.validbox();
        const auto hcArr = hybridCollarMask[iLev][mfi].array();
        const auto csArr = cellSolver[iLev][mfi].array();
        const auto fpArr = isFullPic[mfi].array();
        const Box fpBox = isFullPic[mfi].box();
        const int rK = (activeDim > 2) ? collarRadius : 0;
        const int rJ = (activeDim > 1) ? collarRadius : 0;

        ParallelFor(vbox, [&](int i, int j, int k) noexcept {
          if (csArr(i, j, k) != static_cast<int>(PicSolverMode::Hybrid)) {
            hcArr(i, j, k) = 0;
            return;
          }
          bool nearFullPic = false;
          for (int kk = k - rK; kk <= k + rK && !nearFullPic; ++kk) {
            for (int jj = j - rJ; jj <= j + rJ && !nearFullPic; ++jj) {
              for (int ii = i - collarRadius; ii <= i + collarRadius; ++ii) {
                const IntVect nIdx(AMREX_D_DECL(ii, jj, kk));
                if (fpBox.contains(nIdx) && fpArr(nIdx) == 1) {
                  nearFullPic = true;
                  break;
                }
              }
            }
          }
          hcArr(i, j, k) = nearFullPic ? 1 : 2;
        });
      }
      hybridCollarMask[iLev].FillBoundary(Geom(iLev).periodicity());
      if (isFake2D) {
        for (MFIter mfi(hybridCollarMask[iLev]); mfi.isValid(); ++mfi) {
          const Box& fbox = mfi.fabbox();
          const auto hcArr = hybridCollarMask[iLev][mfi].array();
          ParallelFor(fbox, [&](int i, int j, int k) noexcept {
            if (k != 0)
              hcArr(i, j, k) = hcArr(i, j, 0);
          });
        }
      }
    }

    if (multiSolver && blendWidth > 0) {
      distribute_FabArray(nodeBlendWeight[iLev], nGrids[iLev],
                          DistributionMap(iLev), 1, 1, false);
      const int secMode = static_cast<int>(
          has_active_solver(PicSolverMode::Hybrid) ? PicSolverMode::Hybrid
                                                   : PicSolverMode::Explicit);
      const int rad = blendWidth + 1;
      iMultiFab isSecondary(cGrids[iLev], DistributionMap(iLev), 1, rad);
      isSecondary.setVal(0);
      for (MFIter mfi(isSecondary); mfi.isValid(); ++mfi) {
        const Box& vbox = mfi.validbox();
        const auto secArr = isSecondary[mfi].array();
        const auto csArr = cellSolver[iLev][mfi].array();
        ParallelFor(vbox, [&](int i, int j, int k) noexcept {
          secArr(i, j, k) = (csArr(i, j, k) == secMode) ? 1 : 0;
        });
      }
      isSecondary.FillBoundary(Geom(iLev).periodicity());

      const Real pi = Real(4.0) * std::atan(Real(1.0));
      const int rK = (activeDim > 2) ? rad : 0;
      const int rJ = (activeDim > 1) ? rad : 0;
      for (MFIter mfi(nodeBlendWeight[iLev]); mfi.isValid(); ++mfi) {
        const Box& vbox = mfi.validbox();
        const auto wArr = nodeBlendWeight[iLev][mfi].array();
        const auto secArr = isSecondary[mfi].array();
        const Box secBox = isSecondary[mfi].box();
        nodeBlendWeight[iLev][mfi].setVal<RunOn::Device>(Real(1.0));

        ParallelFor(vbox, [&](int i, int j, int k) noexcept {
          int minDist = blendWidth;
          for (int kk = k - rK; kk <= k + rK - (activeDim > 2 ? 1 : 0); ++kk) {
            const int dk =
                (activeDim > 2) ? ((kk >= k) ? (kk - k) : (k - 1 - kk)) : 0;
            if (dk >= minDist)
              continue;
            for (int jj = j - rJ; jj <= j + rJ - (activeDim > 1 ? 1 : 0);
                 ++jj) {
              const int dj =
                  (activeDim > 1) ? ((jj >= j) ? (jj - j) : (j - 1 - jj)) : 0;
              const int djk = std::max(dk, dj);
              if (djk >= minDist)
                continue;
              for (int ii = i - rad; ii <= i + rad - 1; ++ii) {
                const int di = (ii >= i) ? (ii - i) : (i - 1 - ii);
                const int dist = std::max(djk, di);
                if (dist >= minDist)
                  continue;
                const IntVect cIdx(AMREX_D_DECL(ii, jj, kk));
                if (secBox.contains(cIdx) && secArr(cIdx) == 1) {
                  minDist = dist;
                }
              }
            }
          }
          if (minDist <= 0) {
            wArr(i, j, k) = Real(0.0);
          } else if (minDist >= blendWidth) {
            wArr(i, j, k) = Real(1.0);
          } else {
            wArr(i, j, k) =
                Real(0.5) *
                (Real(1.0) - std::cos(pi * Real(minDist) / Real(blendWidth)));
          }
        });
      }
      nodeBlendWeight[iLev].FillBoundary(Geom(iLev).periodicity());
    } else {
      nodeBlendWeight[iLev].clear();
    }
  }
}
