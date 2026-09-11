
#include "Bit.h"
#include "FleksDistributionMap.h"
#include "Grid.h"
#include "GridUtility.h"

using namespace amrex;

Vector<DistributionMapping> Grid::calc_balanced_maps(bool doSplitLevs) {
  BL_PROFILE("calc_balanced_maps");

  Vector<DistributionMapping> dmap(n_lev_max());

  Vector<MultiFab> cost(n_lev_max());

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    distribute_FabArray(cost[iLev], cGrids[iLev], DistributionMap(iLev), 1, 0,
                        false);
    MultiFab::Copy(cost[iLev], cellCost[iLev], 0, 0, 1, 0);
  }

  Vector<int> rankStart(n_lev(), 0);
  Vector<int> nProcEachLev(n_lev(), ParallelDescriptor::NProcs());

  if (doSplitLevs) {
    Real totalCost = 0;
    Vector<Real> levCost(n_lev());

    for (int iLev = 0; iLev < n_lev(); iLev++) {
      levCost[iLev] = cost[iLev].sum();
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

  // Real localProcCost = 0;
  // Vector<Real> pcost(ParallelDescriptor::NProcs(), 0);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    Vector<int> ord(ParallelDescriptor::NProcs());
    for (int i = 0; i < nProcEachLev[iLev]; ++i) {
      ord[i] = i + rankStart[iLev];
    }

    Real eff;
    dmap[iLev] = FleksDistributionMap::make_balanced_map(
        BalanceMethod::SFC, cost[iLev], nProcEachLev[iLev], ord, eff);
    // Print() << printPrefix << " iLev = " << iLev
    //         << " load balance efficiency = " << std::setw(10) << eff
    //         << std::endl;

    distribute_FabArray(cost[iLev], cGrids[iLev], dmap[iLev], 1, 0, true);

    //   for (MFIter mfi(cost[iLev]); mfi.isValid(); ++mfi) {
    //     localProcCost += cost[iLev][mfi].sum<RunOn::Device>(mfi.validbox(),
    //     0);
    //   }

    //   ParallelDescriptor::Gather(&localProcCost, 1, pcost.data(), 1,
    //                              ParallelDescriptor::IOProcessorNumber());

    //   ParallelDescriptor::Bcast(pcost.data(), pcost.size(),
    //                             ParallelDescriptor::IOProcessorNumber());

    //   using LIpair = std::pair<Long, int>;

    //   Vector<LIpair> pair;
    //   pair.reserve(ParallelDescriptor::NProcs());

    //   for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
    //     pair.push_back(LIpair(pcost[i], i));
    //   }

    //   std::sort(pair.begin(), pair.end(),
    //             [](const LIpair& lhs, const LIpair& rhs) {
    //               return lhs.first > rhs.first;
    //             });

    //   for (int i = 0; i < pcost.size(); ++i) {
    //     ord[i] = pair[i].second;
    //   }
  }

  return dmap;
}

//==========================================================
void Grid::load_balance(const Grid* other, bool doSplitLevs) {

  if (other) {
    regrid(other->get_base_grid(), other, true);
  } else {
    Vector<DistributionMapping> dmap = calc_balanced_maps(doSplitLevs);

    Grid grid(Geom(0), get_amr_info(), nGst, gridID);
    grid.SetFinestLevel(n_lev() - 1);
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      grid.SetBoxArray(iLev, boxArray(iLev));
      grid.SetDistributionMap(iLev, dmap[iLev]);
    }

    regrid(grid.get_base_grid(), &grid, true);
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

    if (refineRegions.empty()) {
      refineRegions.resize(n_lev_max());
    }

    // Why need 'isNewGrid'? See the explanation in Domain::regrid().
    if (region == activeRegion && !isNewGrid)
      return;

    pre_regrid();

    cGridsOld = cGrids;

    doNeedFillNewCell = true;

    activeRegion = region;
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

  // Print() << "dm = " << DistributionMap(0) << std::endl;

  calc_node_grids();

  print_grid_info();

  // If regrid() is called from from read_restart(), activeRegion is not
  // simplifed. Simplify it here.
  activeRegion = activeRegion.simplified();

  domainRange.clear();
  d_domainRange.clear();
  for (int iBox = 0; iBox < activeRegion.size(); iBox++) {
    RealBox rb(activeRegion[iBox], Geom(0).CellSize(), Geom(0).Offset());
    domainRange.push_back(rb);
    d_domainRange.push_back(rb);
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
void Grid::distribute_grid_arrays(const Vector<BoxArray>& cGridsOld) {
  if (h_cellStatus.size() < static_cast<size_t>(n_lev())) {
    h_cellStatus.resize(n_lev());
  }
  if (h_nodeStatus.size() < static_cast<size_t>(n_lev())) {
    h_nodeStatus.resize(n_lev());
  }
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    distribute_FabArray(cellStatus[iLev], cGrids[iLev], DistributionMap(iLev),
                        1, nGst, false);

    h_cellStatus[iLev].define(cGrids[iLev], DistributionMap(iLev), 1, nGst,
                              amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));

    distribute_FabArray(nodeStatus[iLev], nGrids[iLev], DistributionMap(iLev),
                        1, nGst, false);

    h_nodeStatus[iLev].define(nGrids[iLev], DistributionMap(iLev), 1, nGst,
                              amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));

    distribute_FabArray(nodeOffsetMap[iLev], nGrids[iLev],
                        DistributionMap(iLev), 1, 0, false);

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
      const auto cellArr = cellStatus[iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        bit::set_lev_boundary(cellArr(i, j, k));
        bit::set_not_domain_boundary(cellArr(i, j, k));
      });
    }
    // Set 'boundary', 'new' status.
    const bool hasOldGrids = !cGridsOld.empty();
    amrex::Gpu::DeviceVector<Box> d_oldBoxes;
    const Box* pOldBoxes = nullptr;
    int nOldBoxes = 0;
    if (hasOldGrids && iLev < static_cast<int>(cGridsOld.size()) && !cGridsOld[iLev].empty()) {
      const auto& ba = cGridsOld[iLev];
      nOldBoxes = static_cast<int>(ba.size());
      std::vector<Box> h_boxes(nOldBoxes);
      for (int b = 0; b < nOldBoxes; ++b) {
        h_boxes[b] = ba[b];
      }
      d_oldBoxes.resize(nOldBoxes);
      amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_boxes.data(), h_boxes.data() + nOldBoxes, d_oldBoxes.data());
      pOldBoxes = d_oldBoxes.data();
    }

    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<int> cellArr = cellStatus[iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        // Not boundary cell
        bit::set_not_lev_boundary(cellArr(i, j, k));

        // New active cell
        bit::set_new(cellArr(i, j, k));

        if (pOldBoxes) {
          const IntVect iv{ AMREX_D_DECL(i, j, k) };
          for (int b = 0; b < nOldBoxes; ++b) {
            if (pOldBoxes[b].contains(iv)) {
              bit::set_not_new(cellArr(i, j, k));
              break;
            }
          }
        }
      });
    }

    // Set the 'refined' status
    if (iLev < n_lev() - 1) {
      const int iRefined = 1, iNotRefined = 2;
      auto iRefine = makeFineMask(grids[iLev], dmap[iLev], grids[iLev + 1],
                                  ref_ratio[iLev], iNotRefined, iRefined);

      for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<int> cellArr = cellStatus[iLev][mfi].array();
        const auto iRef = iRefine[mfi].array();
        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (iRef(i, j, k) == iRefined) {
            bit::set_refined(cellArr(i, j, k));
          }
        });
      }
    }

    cellStatus[iLev].FillBoundary(Geom(iLev).periodicity());

    const auto geomdata = Geom(iLev).data();
    const RealBox* d_ranges = device_domain_range();
    const int nRanges = domain_range_size();

    // Find domain boundary cells
    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const Array4<int> cellArr = cellStatus[iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (bit::is_lev_boundary(cellArr(i, j, k))) {
          Real xyz[3] = {
            geomdata.ProbLo(0) + (i + 0.5) * geomdata.CellSize(0),
            (AMREX_SPACEDIM > 1) ? (geomdata.ProbLo(1) + (j + 0.5) * geomdata.CellSize(1)) : 0.0,
            (AMREX_SPACEDIM > 2) ? (geomdata.ProbLo(2) + (k + 0.5) * geomdata.CellSize(2)) : 0.0
          };
          bool inside = false;
          for (int r = 0; r < nRanges; ++r) {
            if (d_ranges[r].contains(xyz)) {
              inside = true;
              break;
            }
          }
          if (!inside) {
            bit::set_domain_boundary(cellArr(i, j, k));
          }
        }
      });
    }

    // Set the edge cells.
    // Q: But what is the edge cell?
    // A: It is a physical cell that has one or more neighbor cells are
    // boundary cell.
    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<int>& cellArr = cellStatus[iLev][mfi].array();
      // Flatten inner subBox loop: check all 26 neighbors.
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        for (int kk = k - 1; kk <= k + 1; ++kk)
          for (int jj = j - 1; jj <= j + 1; ++jj)
            for (int ii = i - 1; ii <= i + 1; ++ii) {
              if (bit::is_lev_boundary(cellArr(ii, jj, kk))) {
                bit::set_lev_edge(cellArr(i, j, k));

                if (bit::is_domain_boundary(cellArr(ii, jj, kk))) {
                  bit::set_domain_edge(cellArr(i, j, k));
                }
              }
            }
      });
    }

    // Find cells with 'is_refined' neighbors
    for (MFIter mfi(cellStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto status = cellStatus[iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int kmin = nDim > 2 ? k - 1 : k;
        int kmax = nDim > 2 ? k + 1 : k;
        for (int ii = i - 1; ii <= i + 1; ii++) {
          for (int jj = j - 1; jj <= j + 1; jj++) {
            for (int kk = kmin; kk <= kmax; kk++) {
              if (bit::is_refined(status(ii, jj, kk)) &&
                  !bit::is_refined(status(i, j, k))) {
                bit::set_refined_neighbour(status(i, j, k));
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
          const Array4<int> cellArr = cellStatus[iLev][mfi].array();
          ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (k < -1 || k > 1)
              cellArr(i, j, k) = cellArr(i, j, 0);
          });
        }
    }

    h_cellStatus[iLev].ParallelCopy(cellStatus[iLev]);
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
      const auto nodeArr = nodeStatus[iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        bit::set_lev_boundary(nodeArr(i, j, k));
        bit::set_not_domain_boundary(nodeArr(i, j, k));
        bit::set_not_refined(nodeArr(i, j, k));
      });
    }

    const bool hasOldGrids = !cGridsOld.empty();
    amrex::Gpu::DeviceVector<Box> d_oldNodeBoxes;
    const Box* pOldNodeBoxes = nullptr;
    int nOldNodeBoxes = 0;
    if (hasOldGrids && iLev < static_cast<int>(cGridsOld.size()) && !cGridsOld[iLev].empty()) {
      BoxArray nodeBAOld = convert(cGridsOld[iLev], IntVect(1));
      nOldNodeBoxes = static_cast<int>(nodeBAOld.size());
      std::vector<Box> h_nodeBoxes(nOldNodeBoxes);
      for (int b = 0; b < nOldNodeBoxes; ++b) {
        h_nodeBoxes[b] = nodeBAOld[b];
      }
      d_oldNodeBoxes.resize(nOldNodeBoxes);
      amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_nodeBoxes.data(), h_nodeBoxes.data() + nOldNodeBoxes, d_oldNodeBoxes.data());
      pOldNodeBoxes = d_oldNodeBoxes.data();
    }

    // Set 'boundary', 'new' status.
    for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto nodeArr = nodeStatus[iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        // Not boundary cell
        bit::set_not_lev_boundary(nodeArr(i, j, k));

        // New active cell
        bit::set_new(nodeArr(i, j, k));

        if (pOldNodeBoxes) {
          const IntVect iv{ AMREX_D_DECL(i, j, k) };
          for (int b = 0; b < nOldNodeBoxes; ++b) {
            if (pOldNodeBoxes[b].contains(iv)) {
              bit::set_not_new(nodeArr(i, j, k));
              break;
            }
          }
        }
      });
    }

    nodeStatus[iLev].FillBoundary(Geom(iLev).periodicity());

    const auto geomdataNode = Geom(iLev).data();
    const RealBox* d_ranges = device_domain_range();
    const int nRanges = domain_range_size();

    // Find domain boundary cells
    for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const Array4<int> nodeArr = nodeStatus[iLev][mfi].array();

      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (bit::is_lev_boundary(nodeArr(i, j, k))) {
          Real xyz[3] = {
            geomdataNode.ProbLo(0) + i * geomdataNode.CellSize(0),
            (AMREX_SPACEDIM > 1) ? (geomdataNode.ProbLo(1) + j * geomdataNode.CellSize(1)) : 0.0,
            (AMREX_SPACEDIM > 2) ? (geomdataNode.ProbLo(2) + k * geomdataNode.CellSize(2)) : 0.0
          };
          bool inside = false;
          for (int r = 0; r < nRanges; ++r) {
            if (d_ranges[r].contains(xyz)) {
              inside = true;
              break;
            }
          }
          if (!inside) {
            bit::set_domain_boundary(nodeArr(i, j, k));
          }
        }
      });
    }

    for (MFIter mfi(nodeStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<int> nodeArr = nodeStatus[iLev][mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      { // Set 'owner' status
        const Box cellBox = convert(box, { AMREX_D_DECL(0, 0, 0) });
        const Array4<int const> cell = cellStatus[iLev][mfi].const_array();
        const bool isFake2D_val = isFake2D;
        const int nDim_val = nDim;
        int diMax = 0, diMin = -1;
        int djMax = 0, djMin = -1;
        int dkMax = 0, dkMin = -1;
        if (isFake2D_val || nDim_val == 2) {
          dkMin = 0;
        }

        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (!isFake2D_val || k == lo.z) {
            if (i == lo.x || i == hi.x || j == lo.y || j == hi.y ||
                (nDim_val == 3 && !isFake2D_val && (k == lo.z || k == hi.z))) {
              // Block boundary nodes: check if this box is the owner of this node.
              bool isOwner = false;
              for (int dk = dkMax; dk >= dkMin; dk--) {
                for (int dj = djMax; dj >= djMin; dj--) {
                  for (int di = diMax; di >= diMin; di--) {
                    if (!bit::is_lev_boundary(cell(i + di, j + dj, k + dk))) {
                      isOwner = cellBox.contains(
                          IntVect{ AMREX_D_DECL(i + di, j + dj, k + dk) });
                      goto done_owner;
                    }
                  }
                }
              }
            done_owner:
              if (isOwner) {
                bit::set_owner(nodeArr(i, j, k));
              } else {
                bit::set_not_owner(nodeArr(i, j, k));
              }
            } else {
              // Nodes inside the box.
              bit::set_owner(nodeArr(i, j, k));
            }
          }
        });
      }

      // Set the 'edge' status
      // Q: But what is the edge node?
      // A: It is a node at the boundary of a level.

      // Flatten inner subBox loop: check all 26 neighbors.
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        for (int kk = k - 1; kk <= k + 1; ++kk)
          for (int jj = j - 1; jj <= j + 1; ++jj)
            for (int ii = i - 1; ii <= i + 1; ++ii) {
              if (bit::is_lev_boundary(nodeArr(ii, jj, kk))) {
                bit::set_lev_edge(nodeArr(i, j, k));

                if (bit::is_domain_boundary(nodeArr(ii, jj, kk))) {
                  bit::set_domain_edge(nodeArr(i, j, k));
                }
              }
            }
      });

      // Flatten inner subBox(ijk-1, ijk) loop: check 2x2x2 cell stencil.
      const auto cell = cellStatus[iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        for (int kk = k - 1; kk <= k; ++kk)
          for (int jj = j - 1; jj <= j; ++jj)
            for (int ii = i - 1; ii <= i; ++ii) {
              if (bit::is_refined(cell(ii, jj, kk))) {
                bit::set_refined(nodeArr(i, j, k));
              }
            }
      });
    }

    nodeOffsetMap[iLev].setVal(-1);
    nOwnedNodes[iLev].assign(nodeStatus[iLev].local_size(), 0);

    h_nodeStatus[iLev].ParallelCopy(nodeStatus[iLev]);
    amrex::Gpu::streamSynchronize();

    amrex::iMultiFab h_offsetMap(nodeOffsetMap[iLev].boxArray(), nodeOffsetMap[iLev].DistributionMap(),
                                 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    h_offsetMap.setVal(-1);

    for (MFIter mfi(h_nodeStatus[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& nodeArr = h_nodeStatus[iLev][mfi].array();
      const auto& offsetArr = h_offsetMap[mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      int m = 0;
      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_owner(nodeArr(i, j, k))) {
              offsetArr(i, j, k) = m++;
            }
          }
        }
      }
      nOwnedNodes[iLev][mfi.LocalIndex()] = m;
    }

    nodeOffsetMap[iLev].ParallelCopy(h_offsetMap);
  }
}
