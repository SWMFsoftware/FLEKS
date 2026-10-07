#include <AMReX_PlotFileUtil.H>

#include "Bit.h"
#include "FleksDistributionMap.h"
#include "Grid.h"
#include "GridUtility.h"

using namespace amrex;

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
    regrid(other->get_base_grid(), other, true);
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
void Grid::regrid(const BoxArray& region, const Grid* const grid,
                  bool doLoadBalance) {
  std::string nameFunc = "Grid::regrid_base";
  BL_PROFILE(nameFunc);

  if (!doLoadBalance) {

    if (grid) {
      refineRegions = grid->get_refine_regions();
      SetGridEff(grid->gridEff());
    }

    // The base boxes can be identical while new refinement tags require a
    // different fine hierarchy. isNewGrid also forces first initialization.
    if (region == activeRegion && !isNewGrid &&
        !(refineRegions && refineRegions->is_modified()))
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
    if (grid) {
      // Q: Why is it required to set distribution map (inside init_grid())?
      // A: fi and pic should have the same grids and distribution maps.
      // However, it seems AMReX is too smart that it will try to load balance
      // the box arrays so that the distribution maps can be different even
      // the grid is the same. So we need to set the distribution map here.
      set_ba_and_dm(grid);
    } else {
      // This method will call MakeNewLevelFromScratch() and
      // PostProcessBaseGrids()
      InitFromScratch(0.0);
    }
  }

  finish_regrid();
}

void Grid::restore_grid(const Vector<BoxArray>& boxes) {
  if (boxes.empty() || boxes.size() > static_cast<std::size_t>(n_lev_max()))
    Abort("Invalid restart grid level count.");

  cGridsOld = cGrids;
  cGridsOld.resize(n_lev_max());
  activeRegion = boxes[0];
  isGridEmpty = activeRegion.empty();
  SetFinestLevel(static_cast<int>(boxes.size()) - 1);
  for (int iLev = 0; iLev < n_lev_max(); ++iLev) {
    const BoxArray ba = iLev < n_lev() ? boxes[iLev] : BoxArray();
    SetBoxArray(iLev, ba);
    SetDistributionMap(iLev, DistributionMapping(ba));
  }
  finish_regrid();
}

void Grid::finish_regrid() {

  calc_node_grids();

  print_grid_info();

  // If regrid() is called from from read_restart(), activeRegion is not
  // simplifed. Simplify it here.
  activeRegion = activeRegion.simplified();

  domainRange.clear();
  for (int iBox = 0; iBox < activeRegion.size(); iBox++) {
    RealBox rb(activeRegion[iBox], Geom(0).CellSize(), Geom(0).Offset());
    domainRange.push_back(rb);
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
  }

  update_grid_status(cGridsOld);
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

    // Find domain boundary cells
    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
      ParallelFor(box, [&](int i, int j, int k) noexcept {
        if (bit::is_lev_boundary(cellArr(i, j, k))) {
          Real xyz[nDim];
          Geom(iLev).CellCenter({ AMREX_D_DECL(i, j, k) }, xyz);
          if (!is_inside_domain(xyz)) {
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
  for (int i = 0; i < nlev; ++i) {
    tMF.push_back(&mf[i]);
  }
  Vector<int> tmpVint;
  if (var.empty()) {
    for (int i = 0; i < mf[0].nComp(); ++i) {
      var.push_back(std::to_string(i + 1));
    }
  }
  for (int i = 0; i <= nlev; ++i) {
    tmpVint.push_back(0);
  }
  WriteMultiLevelPlotfile(st, nlev, tMF, var, geom, 0.0, tmpVint, ref_ratio);
}
