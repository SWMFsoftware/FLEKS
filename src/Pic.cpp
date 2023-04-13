#include <math.h>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include "GridUtility.h"
#include "LinearSolver.h"
#include "Pic.h"
#include "SimDomains.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

//==========================================================
void Pic::read_param(const std::string& command, ReadParam& param) {

  if (command == "#PIC") {
    param.read_var("usePIC", usePIC);
  } else if (command == "#DIVE") {
    param.read_var("doCorrectDivE", doCorrectDivE);
    if (doCorrectDivE) {
      param.read_var("nDivECorrection", nDivECorrection);
    }
  } else if (command == "#EXPLICITPIC") {
    param.read_var("useExplicitPIC", useExplicitPIC);
  } else if (command == "#EFIELDSOLVER") {
    Real tol;
    int nIter;
    param.read_var("tol", tol);
    param.read_var("nIter", nIter);
    eSolver.set_tol(tol);
    eSolver.set_nIter(nIter);
  } else if (command == "#PARTICLES") {
    param.read_var("npcelx", nPartPerCell[ix_]);
    param.read_var("npcely", nPartPerCell[iy_]);
    param.read_var("npcelz", nPartPerCell[iz_]);
  } else if (command == "#SOURCEPARTICLES") {
    param.read_var("npcelx", nSourcePPC[ix_]);
    param.read_var("npcely", nSourcePPC[iy_]);
    param.read_var("npcelz", nSourcePPC[iz_]);
  } else if (command == "#ELECTRON") {
    param.read_var("qom", qomEl);
  } else if (command == "#DISCRETIZE" || command == "#DISCRETIZATION") {
    param.read_var("theta", fsolver.theta);
    param.read_var("coefDiff", fsolver.coefDiff);
  } else if (command == "#SMOOTHE") {
    param.read_var("doSmoothE", doSmoothE);
    if (doSmoothE) {
      param.read_var("nSmoothE", nSmoothE);
      param.read_var("coefStrongSmooth", coefStrongSmooth);
      param.read_var("coefWeakSmooth", coefWeakSmooth);
      param.read_var("strongSmoothMach", strongSmoothMach);
      param.read_var("weakSmoothMach", weakSmoothMach);
    }
  } else if (command == "#RESAMPLING") {
    param.read_var("doReSampling", doReSampling);
    if (doReSampling) {
      param.read_var("reSamplingLowLimit", reSamplingLowLimit);
      param.read_var("reSamplingHighLimit", reSamplingHighLimit);
    }
  } else if (command == "#MERGEEFFICIENCY") {
    param.read_var("mergeThresholdDistance", particleMergeThreshold);
    param.read_var("fastMerge", fastMerge);
  } else if (command == "#TESTCASE") {
    std::string testcase;
    param.read_var("testCase", testcase);
    if (testcase == "TwoStream") {
      testCase = TwoStream;
    }
  }
}

//==========================================================
void Pic::post_process_param() {
  fi->set_plasma_charge_and_mass(qomEl);
  nSpecies = fi->get_nS();
}

//==========================================================
void Pic::fill_new_cells() {
  std::string nameFunc = "Pic::fill_new_cells";

  if (isGridEmpty)
    return;

  if (usePIC && !doNeedFillNewCell)
    return;

  timing_func(nameFunc);

  if (!usePIC) {
    // If this method is called when PIC component is off, it suggests the test
    // particle component is activated. The test particle component copies EM
    // field from PIC, so PIC EM field should be updated here.
    if (!nodeStatus.empty()) {
      nodeStatus.setVal(iBoundary_);
      nodeStatus.setVal(iOnNew_, 0);
    }

    if (!cellStatus.empty()) {
      cellStatus.setVal(iBoundary_);
      cellStatus.setVal(iOnNew_, 0);
    }
  }

  fill_E_B_fields();
  
  if (usePIC) {
    fill_particles();
    sum_moments(true);
    sum_to_center(false);
  }

  doNeedFillNewCell = false;
}

//==========================================================
void Pic::init_Pic() {
  const int nLev = max_level + 1;
  {
    centerB.resize(nLev);
    nodeB.resize(nLev);
    nodeE.resize(nLev);
    nodeEth.resize(nLev);
  }
}
//==========================================================
void Pic::distribute_arrays() {
  if (cGrids[0].empty()) {
    return;
  }

  for (int lev = 0; lev <= finest_level; lev++) {
    distribute_FabArray(centerB[lev], cGrids[lev], DistributionMap(lev), 3,
                        nGst);
    distribute_FabArray(nodeB[lev], nGrids[lev], DistributionMap(lev), 3, nGst);
    distribute_FabArray(nodeE[lev], nGrids[lev], DistributionMap(lev), 3, nGst);
    distribute_FabArray(nodeEth[lev], nGrids[lev], DistributionMap(lev), 3,
                        nGst);
  }
}

//==========================================================
void Pic::regrid(const BoxArray& region, const Grid* const grid) {
  std::string nameFunc = "Pic::regrid";

  timing_func(nameFunc);

  // Why need 'isGridInitialized'? See the explaination in Domain::regrid().
  if (region == activeRegion && isGridInitialized)
    return;

  activeRegion = region;

  isGridEmpty = activeRegion.empty();

  doNeedFillNewCell = true;

  BoxArray cGridOld = cGrids[0];

  if (isGridEmpty) {
    cGrids.clear();
    cGrids.push_back(amrex::BoxArray());
  } else {
    // This method will call MakeNewLevelFromScratch() and
    // PostProcessBaseGrids()
    InitFromScratch(tc->get_time());
    for (int iLev = 0; iLev <= max_level; iLev++) {
      // Q: Why is it required to set distribution map here?
      // A: fi and pic should have the same grids and distribution maps.
      // However, it seems AMReX is too smart that it will try to load balance
      // the box arrays so that the distribution maps can be different even the
      // grid is the same. So we need to set the distribution map here.
      SetDistributionMap(iLev, grid->DistributionMap(iLev));
    }
  }

  calc_node_grids();

  print_grid_info(false);

  distribute_arrays();

  //===========Move field data around begin====================
  // distribute_FabArray(nodeE, nGrids[0], DistributionMap(0), 3, nGst);
  // distribute_FabArray(nodeEth, nGrids[0], DistributionMap(0), 3, nGst);
  // distribute_FabArray(centerB[0], cGrids[0], DistributionMap(0), 3, nGst);
  // distribute_FabArray(nodeB[0], cGrids[0], DistributionMap(0), 3, nGst);

  distribute_FabArray(centerNetChargeOld, cGrids[0], DistributionMap(0), 1,
                      nGst);
  distribute_FabArray(centerNetChargeN, cGrids[0], DistributionMap(0), 1,
                      nGst); // false??
  distribute_FabArray(centerNetChargeNew, cGrids[0], DistributionMap(0), 1,
                      nGst); // false??

  distribute_FabArray(centerDivE, cGrids[0], DistributionMap(0), 1, nGst);
  distribute_FabArray(centerPhi, cGrids[0], DistributionMap(0), 1, nGst);

  {
    bool doMoveData = false;

    iTot = nSpecies;
    if (plasmaEnergy.empty()) {
      plasmaEnergy.resize(nSpecies + 1);
    }

    if (nodePlasma.empty()) {
      // The last one is the sum of all species.
      nodePlasma.resize(nSpecies + 1);
      for (auto& pl : nodePlasma) {
        pl.define(nGrids[0], DistributionMap(0), nMoments, nGst);
        pl.setVal(0.0);
      }
    } else {
      for (auto& pl : nodePlasma) {
        distribute_FabArray(pl, nGrids[0], DistributionMap(0), nMoments, nGst,
                            doMoveData);
      }
    }

    distribute_FabArray(jHat, nGrids[0], DistributionMap(0), 3, nGst,
                        doMoveData);

    if (!useExplicitPIC) {
      distribute_FabArray(nodeMM, nGrids[0], DistributionMap(0), 1, 1,
                          doMoveData);
    }
    distribute_FabArray(centerMM, cGrids[0], DistributionMap(0), 1, nGst,
                        doMoveData);

    distribute_FabArray(nodeSmoothCoef, nGrids[0], DistributionMap(0), 1, nGst,
                        doMoveData);
  }
  //===========Move field data around end====================

  { //===========Label cellStatus/nodeStatus ==========

    // The algorithm decides inject particles or not needs at least 2 ghost cell
    // layers.
    distribute_FabArray(cellStatus, cGrids[0], DistributionMap(0), 1,
                        nGst >= 2 ? nGst : 2, false);
    if (!cellStatus.empty()) {
      cellStatus.setVal(iBoundary_);
      cellStatus.setVal(iOnNew_, 0);

      if (usePIC)
        for (MFIter mfi(cellStatus); mfi.isValid(); ++mfi) {
          const Box& box = mfi.validbox();
          const Array4<int>& cellArr = cellStatus[mfi].array();
          const auto lo = lbound(box);
          const auto hi = ubound(box);

          for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
              for (int i = lo.x; i <= hi.x; ++i) {
                if (cellArr(i, j, k) == iOnNew_ &&
                    cGridOld.contains(IntVect{ i, j, k })) {
                  cellArr(i, j, k) = iOnOld_;
                }
              }
        }
    }

    cellStatus.FillBoundary(Geom(0).periodicity());

    if (isFake2D) {
      // For the fake 2D cases, in the z-direction, only the first layer ghost
      // cells are filled in correctly by the method FillBoundary.
      if (!cellStatus.empty())
        for (MFIter mfi(cellStatus); mfi.isValid(); ++mfi) {
          const Box& box = mfi.fabbox();
          const Array4<int>& cellArr = cellStatus[mfi].array();
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

    distribute_FabArray(nodeStatus, nGrids[0], DistributionMap(0), 1, nGst,
                        false);
    if (!nodeStatus.empty()) {
      nodeStatus.setVal(iBoundary_);
      nodeStatus.setVal(iOnNew_, 0);

      if (usePIC) {
        auto nodeBAOld =
            convert(cGridOld, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });
        for (MFIter mfi(nodeStatus); mfi.isValid(); ++mfi) {
          const Box& box = mfi.validbox();
          const auto& nodeArr = nodeStatus[mfi].array();

          const auto lo = lbound(box);
          const auto hi = ubound(box);

          for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
              for (int i = lo.x; i <= hi.x; ++i) {
                if (nodeArr(i, j, k) == iOnNew_ &&
                    nodeBAOld.contains(IntVect{ i, j, k })) {
                  nodeArr(i, j, k) = iOnOld_;
                }
              }
        }
      }
    }

    nodeStatus.FillBoundary(Geom(0).periodicity());

    distribute_FabArray(nodeShare, nGrids[0], DistributionMap(0), 1, 0, false);
    set_nodeShare();
  }

  //--------------particles-----------------------------------
  if (parts.empty()) {
    solveEM = true;
    for (int i = 0; i < nSpecies; i++) {
      auto ptr = std::unique_ptr<Particles<> >(new Particles<>(
          this, fi.get(), tc.get(), i, fi->get_species_charge(i),
          fi->get_species_mass(i), nPartPerCell, testCase));

      // If contains neutrals, assume this is OH-PT coupling, and do not solve
      // for EM fields.
      if (ptr->is_neutral())
        solveEM = false;

      ptr->set_region_range(activeRegion);

      //----- Set parameters------------
      if (particleMergeThreshold >= 0) {
        ptr->set_merge_threshold(particleMergeThreshold);
      }

      if (particleMergeBinBuffer >= 0) {
        ptr->set_merge_velocity_bin_buffer(particleMergeBinBuffer);
      }

      ptr->fast_merge(fastMerge);
      //----------------------------------

      parts.push_back(std::move(ptr));
    }
  } else {
    for (int i = 0; i < nSpecies; i++) {
      // Label the particles outside the OLD PIC region.
      parts[i]->label_particles_outside_ba();

      parts[i]->set_region_range(activeRegion);

      // Label the particles outside the NEW PIC region.
      parts[i]->label_particles_outside_ba_general();

      if (cGrids[0].size() > 0) {
        parts[i]->SetParticleBoxArray(0, cGrids[0]);
        parts[i]->SetParticleDistributionMap(0, DistributionMap(0));
      }
      parts[i]->Redistribute();
    }
  }

  { // Copy cellStatus to Particles objects.
    for (int i = 0; i < nSpecies; i++) {
      distribute_FabArray(parts[i]->cellStatus, cGrids[0], DistributionMap(0),
                          1, nGst >= 2 ? nGst : 2, false);

      if (!cellStatus.empty()) {
        iMultiFab::Copy(parts[i]->cellStatus, cellStatus, 0, 0,
                        cellStatus.nComp(), cellStatus.nGrow());
      }
    }
  }
  //--------------particles-----------------------------------

  {
    for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
      int n = get_local_node_or_cell_number(nodeE[iLevTest]);
      eSolver.init(n, nDim, nDim, matvec_E_solver);
    }
  }

  {
    int n = get_local_node_or_cell_number(centerDivE);
    divESolver.init(n, 1, nDim, matvec_divE_accurate);
  }

  isGridInitialized = true;
}

//==========================================================
void Pic::set_nodeShare() {
  if (!nodeShare.empty())
    nodeShare.setVal(iIgnore_);

  const Box& gbx = convert(Geom(0).Domain(), { 0, 0, 0 });

  if (!nodeShare.empty())
    for (MFIter mfi(nodeShare); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();

      const Box& cellBox = convert(box, { 0, 0, 0 });

      const auto& typeArr = nodeShare[mfi].array();
      const auto& statusArr = cellStatus[mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      int diMax = 0, diMin = -1;
      int djMax = 0, djMin = -1;
      int dkMax = 0, dkMin = -1;
      if (isFake2D) {
        dkMin = 0;
      }

      auto fHandle = [&](int i, int j, int k) {
        for (int dk = dkMax; dk >= dkMin; dk--)
          for (int dj = djMax; dj >= djMin; dj--)
            for (int di = diMax; di >= diMin; di--) {
              if (statusArr(i + di, j + dj, k + dk) != iBoundary_) {
                // Find the first CELL that shares this node.
                if (cellBox.contains(IntVect{ i + di, j + dj, k + dk })) {
                  return iAssign_;
                } else {
                  int ii = 0;
                  if (di == 0)
                    ii += 1;
                  if (dj == 0)
                    ii += 2;
                  if (!isFake2D && dk == 0)
                    ii += 4;

                  return ii;
                }
              }
            }
        Abort("Error: something is wrong here!");
        return 0;
      };

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (!isFake2D || k == lo.z) {
              // for 2D (1 cell in the z-direction), only handle the layer of
              // k=0
              if (i == lo.x || i == hi.x || j == lo.y || j == hi.y ||
                  (!isFake2D && (k == lo.z || k == hi.z))) {
                // Block boundary nodes.

                // Check if this block boundary node needs to be handled by this
                // block.

                typeArr(i, j, k) = fHandle(i, j, k);
              } else {
                typeArr(i, j, k) = iAssign_;
              }
            }
          }
    }
}

//==========================================================
void Pic::fill_new_node_E() {
  // for (int iLev = 0; iLev < nodeE.size(); iLev++) //-Talha-Eventually needs to be implemented - FluidInterface already capable - checked
  int iLev = 0;
  for (MFIter mfi(nodeE[iLev]); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodeE[iLev][mfi];
    const Box& box = mfi.validbox();
    const Array4<Real>& arrE = fab.array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const auto& status = nodeStatus[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          if (status(i, j, k) == iOnNew_) {
            arrE(i, j, k, ix_) = fi->get_ex(mfi, i, j, k, iLev);
            arrE(i, j, k, iy_) = fi->get_ey(mfi, i, j, k, iLev);
            arrE(i, j, k, iz_) = fi->get_ez(mfi, i, j, k, iLev);
          }
        }
  }

  if (max_level > 0) {
    InterpFromCoarseAllLevels(nodeE, finest_level);
  }
}

//==========================================================
void Pic::fill_new_node_B() {
  // for (int iLev = 0; iLev < nodeE.size(); iLev++)
  int iLev = 0;
  for (MFIter mfi(nodeB[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrB = nodeB[iLev][mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const auto& status = nodeStatus[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          if (status(i, j, k) == iOnNew_) {
            arrB(i, j, k, ix_) = fi->get_bx(mfi, i, j, k, iLev);
            arrB(i, j, k, iy_) = fi->get_by(mfi, i, j, k, iLev);
            arrB(i, j, k, iz_) = fi->get_bz(mfi, i, j, k, iLev);
          }
        }
  }
if (max_level > 0) {
    InterpFromCoarseAllLevels(nodeB, finest_level);
  }

}

//==========================================================
void Pic::fill_new_center_B() {
  // for (int iLev = 0; iLev < nodeE.size(); iLev++)
  int iLev = 0;
  for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& centerArr = centerB[iLev][mfi].array();
    const auto& nodeArr = nodeB[iLev][mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const auto& status = cellStatus[mfi].array();

    for (int iVar = 0; iVar < centerB[iLev].nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (status(i, j, k) != iOnOld_) {
              centerArr(i, j, k, iVar) = 0;
              for (int di = 0; di <= 1; di++)
                for (int dj = 0; dj <= 1; dj++)
                  for (int dk = 0; dk <= 1; dk++) {
                    centerArr(i, j, k, iVar) +=
                        0.125 * nodeArr(i + di, j + dj, k + dk, iVar);
                  }
            }
          }
  }

  if (max_level > 0) {
    InterpFromCoarseAllLevels(centerB, finest_level);
  }
}

//==========================================================
void Pic::fill_E_B_fields() {
  if (!solveEM)
    return;

  fill_new_node_E();
  fill_new_node_B();

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    nodeE[iLevTest].FillBoundary(Geom(iLevTest).periodicity());
    nodeB[iLevTest].FillBoundary(Geom(iLevTest).periodicity());
    apply_BC(nodeStatus, nodeB[iLevTest], 0, nDim, &Pic::get_node_B, iLevTest);
    apply_BC(nodeStatus, nodeE[iLevTest], 0, nDim, &Pic::get_node_E, iLevTest);
  }

  fill_new_center_B();

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    centerB[iLevTest].FillBoundary(Geom(iLevTest).periodicity());

    apply_BC(cellStatus, centerB[iLevTest], 0, centerB[iLevTest].nComp(),
             &Pic::get_center_B, iLevTest);
  }
}

//==========================================================
void Pic::fill_particles() {
  inject_particles_for_new_cells();
  inject_particles_for_boundary_cells();
}

void Pic::fill_source_particles() {
  for (int i = 0; i < nSpecies; i++) {
    parts[i]->add_particles_source(cellStatus, *source, tc->get_dt(),
                                   nSourcePPC);
  }
}

//==========================================================
void Pic::update_part_loc_to_half_stage() {
  std::string nameFunc = "Pic::update_part_loc_to_half_stage";

  timing_func(nameFunc);

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    for (int i = 0; i < nSpecies; i++) {
      parts[i]->update_position_to_half_stage(nodeEth[iLevTest],
                                              nodeB[iLevTest], tc->get_dt());
    }
  }

  inject_particles_for_boundary_cells();
}

//==========================================================
void Pic::re_sampling() {
  std::string nameFunc = "Pic::re_sampling";

  timing_func(nameFunc);

  if (doReSampling) {
    for (int i = 0; i < nSpecies; i++) {
      parts[i]->split_particles(reSamplingLowLimit);
      parts[i]->combine_particles(reSamplingHighLimit);
    }
  }
}

//==========================================================
void Pic::particle_mover() {
  std::string nameFunc = "Pic::mover";

  timing_func(nameFunc);

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    if (useExplicitPIC) {
      MultiFab tmpE(nGrids[iLevTest], DistributionMap(iLevTest), 3, nGst);
      // nodeE/nodeEth is at t_n/t_{n+1}, tmpE is at t_{n+0.5}
      MultiFab::LinComb(tmpE, 0.5, nodeEth[iLevTest], 0, 0.5, nodeE[iLevTest],
                        0, 0, nodeE[iLevTest].nComp(), nodeE[iLevTest].nGrow());
      for (int i = 0; i < nSpecies; i++) {
        parts[i]->mover(tmpE, nodeB[iLevTest], tc->get_dt(), tc->get_next_dt());
      }

    } else {

      for (int i = 0; i < nSpecies; i++) {
        parts[i]->mover(nodeEth[iLevTest], nodeB[iLevTest], tc->get_dt(),
                        tc->get_next_dt());
      }
    }
  }
}

//==========================================================
void Pic::calc_mass_matrix() {
  std::string nameFunc = "Pic::calc_mass_matrix";

  if (isGridEmpty)
    return;

  if (!solveEM)
    return;

  timing_func(nameFunc);

  jHat.setVal(0.0);

  if (!useExplicitPIC) {
    const RealMM mm0(0.0);
    nodeMM.setVal(mm0);
  }

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    for (int i = 0; i < nSpecies; i++) {
      if (useExplicitPIC) {
        parts[i]->calc_jhat(jHat, nodeB[iLevTest], tc->get_dt());
      } else {
        parts[i]->calc_mass_matrix(nodeMM, jHat, nodeB[iLevTest], tc->get_dt());
      }
    }
  }

  Real invVol = 1;
  for (int i = 0; i < nDim; i++) {
    invVol *= Geom(0).InvCellSize(i);
  }

  jHat.mult(invVol, 0, jHat.nComp(), jHat.nGrow());

  jHat.SumBoundary(Geom(0).periodicity());

  if (!useExplicitPIC) {
    nodeMM.SumBoundary(Geom(0).periodicity());
    nodeMM.FillBoundary(Geom(0).periodicity());
  }
}

//==========================================================
void Pic::sum_moments(bool updateDt) {
  std::string nameFunc = "Pic::sum_moments";

  if (isGridEmpty)
    return;

  timing_func(nameFunc);

  const auto& dx = Geom(0).CellSize();
  Real minDx = 1e99;
  for (int iDim = 0; iDim < nDim; iDim++) {
    if (minDx > dx[iDim])
      minDx = dx[iDim];
  }

  nodePlasma[nSpecies].setVal(0.0);
  plasmaEnergy[iTot] = 0;
  for (int i = 0; i < nSpecies; i++) {

    for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
      Real energy =
          parts[i]->sum_moments(nodePlasma[i], nodeB[iLevTest], tc->get_dt());
      plasmaEnergy[i] = energy;
      plasmaEnergy[iTot] += energy;
    }
  }

  if (updateDt) {
    Real uMax = 0;
    if (tc->get_cfl() > 0 || doReport) {
      for (int i = 0; i < nSpecies; i++) {
        Real uMaxSpecies = parts[i]->calc_max_thermal_velocity(nodePlasma[i]);
        ParallelDescriptor::ReduceRealMax(uMaxSpecies);

        if (doReport)
          Print() << printPrefix << std::setprecision(5) << "Species " << i
                  << ": max(uth) = " << uMaxSpecies << std::endl;

        if (uMaxSpecies > uMax)
          uMax = uMaxSpecies;
      }
    }

    if (tc->get_cfl() > 0) {
      Real dt = tc->get_cfl() * minDx / uMax;
      tc->set_next_dt(dt);

      if (tc->get_dt() < 0) {
        tc->set_dt(dt);
      }

      if (Particles<>::particlePosition == NonStaggered) {
        tc->set_dt(dt);
      }
    }

    if (doReport) {
      if (Particles<>::particlePosition == Staggered) {
        Print() << printPrefix << std::setprecision(5)
                << "dt = " << tc->get_dt_si()
                << " dtNext = " << tc->get_next_dt_si()
                << " CFL(dtNext) = " << tc->get_next_dt() * uMax / minDx
                << std::endl;
      } else {
        Print() << printPrefix << std::setprecision(5)
                << "dt = " << tc->get_dt_si()
                << " CFL = " << tc->get_dt() * uMax / minDx << std::endl;
      }
    }
  }

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->convert_to_fluid_moments(nodePlasma[i]);
    MultiFab::Add(nodePlasma[nSpecies], nodePlasma[i], 0, 0, nMoments, nGst);
  }
}

//==========================================================
void Pic::divE_correction() {
  if (!solveEM)
    return;

  std::string nameFunc = "Pic::divE_correction";

  timing_func(nameFunc);

  for (int iIter = 0; iIter < nDivECorrection; iIter++) {
    sum_to_center(true);

    if (doReport)
      Print() << "\n-----" << printPrefix << " div(E) correction at iter "
              << iIter << "----------" << std::endl;
    calculate_phi(divESolver);

    divE_correct_particle_position();
  }

  for (int i = 0; i < nSpecies; i++) {
    // The particles outside the simulation domain is marked for deletion
    // inside divE_correct_particle_position(). Redistribute() deletes these
    // particles. In order to get correct moments, re-inject particles in the
    // ghost cells.
    parts[i]->Redistribute();
  }
  inject_particles_for_boundary_cells();

  sum_to_center(false);
}

//==========================================================
void Pic::divE_correct_particle_position() {
  std::string nameFunc = "Pic::correct_position";

  timing_func(nameFunc);

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->divE_correct_position(centerPhi);
  }
}

//==========================================================
void Pic::calculate_phi(LinearSolver& solver) {
  std::string nameFunc = "Pic::calculate_phi";

  timing_func(nameFunc);

  const int iLev = 0;
  MultiFab residual(cGrids[iLev], DistributionMap(iLev), 1, nGst);

  solver.reset(get_local_node_or_cell_number(centerDivE));
  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    div_node_to_center(nodeE[iLevTest], residual, Geom(iLevTest).InvCellSize());
  }

  Real coef = 1;
  if (Particles<>::particlePosition == Staggered) {
    coef = 1.0 / rhoTheta;
  }

  MultiFab::LinComb(residual, coef, residual, 0, -fourPI * coef,
                    centerNetChargeN, 0, 0, residual.nComp(), residual.nGrow());

  convert_3d_to_1d(residual, solver.rhs);

  BL_PROFILE_VAR("Pic::phi_iterate", solve);
  solver.solve(doReport);
  BL_PROFILE_VAR_STOP(solve);

  convert_1d_to_3d(solver.xLeft, centerPhi);
  centerPhi.FillBoundary(Geom(0).periodicity());
}

//==========================================================
void Pic::divE_accurate_matvec(const double* vecIn, double* vecOut) {
  std::string nameFunc = "Pic::divE_matvec";
  timing_func(nameFunc);

  const int iLev = 0;

  zero_array(vecOut, divESolver.get_nSolve());

  MultiFab inMF(cGrids[iLev], DistributionMap(iLev), 1, nGst);

  convert_1d_to_3d(vecIn, inMF);
  inMF.FillBoundary(0, 1, IntVect(1), Geom(0).periodicity());

  MultiFab outMF(cGrids[iLev], DistributionMap(iLev), 1, nGst);
  outMF.setVal(0.0);

  for (amrex::MFIter mfi(inMF); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real>& lArr = outMF[mfi].array();
    const amrex::Array4<amrex::Real const>& rArr = inMF[mfi].array();
    const amrex::Array4<RealCMM>& mmArr = centerMM[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i)

          for (int i2 = i - 1; i2 <= i + 1; i2++)
            for (int j2 = j - 1; j2 <= j + 1; j2++)
              for (int k2 = k - 1; k2 <= k + 1; k2++) {
                const int gp = (i2 - i + 1) * 9 + (j2 - j + 1) * 3 + k2 - k + 1;
                lArr(i, j, k) += rArr(i2, j2, k2) * mmArr(i, j, k).data[gp];
              }
  }
  outMF.mult(fourPI * fourPI);
  convert_3d_to_1d(outMF, vecOut);
}

//==========================================================
void Pic::sum_to_center(bool isBeforeCorrection) {
  std::string nameFunc = "Pic::sum_to_center";

  timing_func(nameFunc);

  centerNetChargeNew.setVal(0.0);

  const RealCMM mm0(0.0);
  centerMM.setVal(mm0);

  bool doNetChargeOnly = !isBeforeCorrection;

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->sum_to_center(centerNetChargeNew, centerMM, doNetChargeOnly);
  }

  if (!doNetChargeOnly) {
    centerMM.SumBoundary(Geom(0).periodicity());
  }

  centerNetChargeNew.SumBoundary(Geom(0).periodicity());

  const int iLevTest = 0;
  apply_BC(cellStatus, centerNetChargeNew, 0, centerNetChargeNew.nComp(),
           &Pic::get_zero, iLevTest);

  if (Particles<>::particlePosition == NonStaggered) {
    MultiFab::Copy(centerNetChargeN, centerNetChargeNew, 0, 0,
                   centerNetChargeN.nComp(), centerNetChargeN.nGrow());
  } else {
    MultiFab::LinComb(centerNetChargeN, 1 - rhoTheta, centerNetChargeOld, 0,
                      rhoTheta, centerNetChargeNew, 0, 0,
                      centerNetChargeN.nComp(), centerNetChargeN.nGrow());

    if (!isBeforeCorrection) {
      MultiFab::Copy(centerNetChargeOld, centerNetChargeNew, 0, 0,
                     centerNetChargeOld.nComp(), centerNetChargeOld.nGrow());
    }
  }
}

//==========================================================
void Pic::update(bool doReportIn) {
  std::string nameFunc = "Pic::update";

  if (isGridEmpty || !usePIC)
    return;

  timing_func(nameFunc);

  doReport = doReportIn;

  Real tStart = second();

  if (Particles<>::particlePosition == NonStaggered) {
    update_part_loc_to_half_stage();
  }

  calc_mass_matrix();

  update_E();

  particle_mover();

  // Calling re_sampling after particle mover so that all the particles outside
  // the domain have been deleted.
  re_sampling();

  charge_exchange();

  if (source) {
    fill_source_particles();
  }

  inject_particles_for_boundary_cells();

  update_B();

  if (doCorrectDivE) {
    divE_correction();
  }

  tc->set_dt(tc->get_next_dt());

  sum_moments(true);

  if (doReport) {
    Real tEnd = second();
    Real nPoint = activeRegion.d_numPts();
    int nProc = amrex::ParallelDescriptor::NProcs();
    // The unit of the speed is (cell per processor per second)
    Real speed = nPoint / nProc / (tEnd - tStart);

    // speedNorm is a value obtained from tests.
    Real speedNorm = 1000;
    Print() << printPrefix
            << "Normalized PIC simulation speed = " << speed / speedNorm
            << " (performance is good if the value >> 1 and bad if <<1 )"
            << std::endl;

    report_load_balance();
  }
}

//==========================================================
void Pic::update_E() {
  if (!solveEM)
    return;

  if (useExplicitPIC) {
    update_E_expl();
  } else {
    update_E_impl();
  }
}

//==========================================================
void Pic::update_E_expl() {
  std::string nameFunc = "Pic::update_E_expl";

  timing_func(nameFunc);

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    MultiFab::Copy(nodeEth[iLevTest], nodeE[iLevTest], 0, 0,
                   nodeE[iLevTest].nComp(), nodeE[iLevTest].nGrow());
    apply_BC(cellStatus, centerB[iLevTest], 0, centerB[iLevTest].nComp(),
             &Pic::get_center_B, iLevTest);
  }
  const Real dt = tc->get_dt();
  Real dt2dx[nDimMax];
  for (int i = 0; i < nDimMax; i++) {
    dt2dx[i] = dt * Geom(0).InvCellSize(i);
  }
  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    curl_center_to_node(centerB[iLevTest], nodeE[iLevTest], dt2dx);
    MultiFab::Saxpy(nodeE[iLevTest], -fourPI * dt, jHat, 0, 0,
                    nodeE[iLevTest].nComp(), nodeE[iLevTest].nGrow());

    MultiFab::Add(nodeE[iLevTest], nodeEth[iLevTest], 0, 0,
                  nodeE[iLevTest].nComp(), nodeE[iLevTest].nGrow());

    nodeE[iLevTest].FillBoundary(Geom(iLevTest).periodicity());
    apply_BC(nodeStatus, nodeE[iLevTest], 0, nDim, &Pic::get_node_E, iLevTest);
  }
}

//==========================================================
void Pic::update_E_impl() {
  std::string nameFunc = "Pic::update_E_impl";

  timing_func(nameFunc);
  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    eSolver.reset(get_local_node_or_cell_number(nodeE[iLevTest]));

    update_E_rhs(eSolver.rhs);

    convert_3d_to_1d(nodeE[iLevTest], eSolver.xLeft);
  }

  update_E_matvec(eSolver.xLeft, eSolver.matvec, false);

  for (int i = 0; i < eSolver.get_nSolve(); i++) {
    eSolver.rhs[i] -= eSolver.matvec[i];
    eSolver.xLeft[i] = 0;
  }

  if (doReport)
    Print() << "\n-------" << printPrefix << " E solver ------------------"
            << std::endl;

  BL_PROFILE_VAR("Pic::E_iterate", eSolver);
  eSolver.solve(doReport);
  BL_PROFILE_VAR_STOP(eSolver);

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    nodeEth[iLevTest].setVal(0.0);
    convert_1d_to_3d(eSolver.xLeft, nodeEth[iLevTest]);
    nodeEth[iLevTest].SumBoundary(Geom(iLevTest).periodicity());
    nodeEth[iLevTest].FillBoundary(Geom(iLevTest).periodicity());
    MultiFab::Add(nodeEth[iLevTest], nodeE[iLevTest], 0, 0,
                  nodeEth[iLevTest].nComp(), nGst);

    MultiFab::LinComb(nodeE[iLevTest], -(1.0 - fsolver.theta) / fsolver.theta,
                      nodeE[iLevTest], 0, 1. / fsolver.theta, nodeEth[iLevTest],
                      0, 0, nodeE[iLevTest].nComp(), nGst);

    apply_BC(nodeStatus, nodeE[iLevTest], 0, nDim, &Pic::get_node_E, iLevTest);
    apply_BC(nodeStatus, nodeEth[iLevTest], 0, nDim, &Pic::get_node_E,
             iLevTest);
    if (doSmoothE) {
      calc_smooth_coef();
      smooth_E(nodeEth[iLevTest]);
      smooth_E(nodeE[iLevTest]);
    }
    div_node_to_center(nodeE[iLevTest], centerDivE,
                       Geom(iLevTest).InvCellSize());
  }
}

//==========================================================
void Pic::update_E_matvec(const double* vecIn, double* vecOut,
                          const bool useZeroBC) {
  std::string nameFunc = "Pic::E_matvec";
  timing_func(nameFunc);

  const int iLev = 0;

  zero_array(vecOut, eSolver.get_nSolve());

  MultiFab vecMF(nGrids[0], DistributionMap(0), 3, nGst);
  vecMF.setVal(0.0);

  MultiFab matvecMF(nGrids[0], DistributionMap(0), 3, 1);
  matvecMF.setVal(0.0);

  MultiFab tempCenter3(cGrids[0], DistributionMap(0), 3, nGst);

  MultiFab tempNode3(nGrids[iLev], DistributionMap(iLev), 3, nGst);
  tempNode3.setVal(0.0);

  MultiFab tempCenter1(cGrids[iLev], DistributionMap(iLev), 1, nGst);

  convert_1d_to_3d(vecIn, vecMF);

  // The right side edges should be filled in.
  vecMF.SumBoundary(Geom(0).periodicity());

  // M*E needs ghost cell information.
  vecMF.FillBoundary(Geom(0).periodicity());

  if (useZeroBC) {
    // The boundary nodes would not be filled in by convert_1d_3d. So, there
    // is not need to apply zero boundary conditions again here.
  } else {
    // Even after apply_BC(), the outmost layer node E is still
    // unknow. See FluidInterface::calc_current for detailed explaniation.
    apply_BC(nodeStatus, vecMF, 0, nDim, &Pic::get_node_E, iLev);
  }

  lap_node_to_node(vecMF, matvecMF, DistributionMap(0), Geom(0), cellStatus);

  Real delt2 = pow(fsolver.theta * tc->get_dt(), 2);
  matvecMF.mult(-delt2);

  { // grad(divE)
    div_node_to_center(vecMF, centerDivE, Geom(0).InvCellSize());

    if (fsolver.coefDiff > 0) {
      // Calculate cell center E for center-to-center divE.
      // The outmost boundary layer of tempCenter3 is not accurate.
      average_node_to_cellcenter(tempCenter3, 0, vecMF, 0, 3,
                                 tempCenter3.nGrow());

      //----The following comments are left here for reference------
      // Q: Why apply float BC for all boundary ghost nodes, instead of just
      // the outmost layer? A: For the example described in
      // FluidInterface::calc_current, cell (c+4, c-1) of tempCenter3-block1
      // is not accurate, so the values at (c+4, c-2) will be wrong if we only
      // apply float BC for the outmost layer.
      // apply_float_boundary(cellStatus, tempCenter3, Geom(0), 0,
      //                           tempCenter3.nComp());
      //------------------------------------------------------------

      div_center_to_center(tempCenter3, tempCenter1, Geom(0).InvCellSize());

      tempCenter1.FillBoundary(0, 1, IntVect(1), Geom(0).periodicity());

      // 1) The outmost boundary layer of tempCenter3 is not accurate.
      // 2) The 2 outmost boundary layers (all ghosts if there are 2 ghost
      // cells) of tempCenter1 are not accurate
      apply_BC(cellStatus, tempCenter1, 0, tempCenter1.nComp(), &Pic::get_zero,
               iLev);

      MultiFab::LinComb(centerDivE, 1 - fsolver.coefDiff, centerDivE, 0,
                        fsolver.coefDiff, tempCenter1, 0, 0, 1, 1);
    }

    grad_center_to_node(centerDivE, tempNode3, Geom(0).InvCellSize());

    tempNode3.mult(delt2);
    MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(),
                  matvecMF.nGrow());
  }

  tempNode3.setVal(0);
  update_E_M_dot_E(vecMF, tempNode3);

  MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(), 0);

  MultiFab::Add(matvecMF, vecMF, 0, 0, matvecMF.nComp(), 0);

  convert_3d_to_1d(matvecMF, vecOut);
}

//==========================================================
void Pic::update_E_M_dot_E(const MultiFab& inMF, MultiFab& outMF) {
  std::string nameFunc = "Pic::update_E_M_dot_E";
  timing_func(nameFunc);

  outMF.setVal(0.0);
  Real c0 = fourPI * fsolver.theta * tc->get_dt();
  for (amrex::MFIter mfi(outMF); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const amrex::Array4<amrex::Real const>& inArr = inMF[mfi].array();
    const amrex::Array4<amrex::Real>& ourArr = outMF[mfi].array();
    const amrex::Array4<RealMM>& mmArr = nodeMM[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          Real* const data0 = mmArr(i, j, k).data;
          for (int k2 = k - 1; k2 <= k + 1; k2++)
            for (int j2 = j - 1; j2 <= j + 1; j2++)
              for (int i2 = i - 1; i2 <= i + 1; i2++) {
                const int gp = (k2 - k + 1) * 9 + (j2 - j + 1) * 3 + i2 - i + 1;
                const int idx0 = gp * 9;

                Real* const M_I = &(data0[idx0]);

                const double& vctX =
                    inArr(i2, j2, k2, ix_); // vectX[i2][j2][k2];
                const double& vctY = inArr(i2, j2, k2, iy_);
                const double& vctZ = inArr(i2, j2, k2, iz_);
                ourArr(i, j, k, ix_) +=
                    (vctX * M_I[0] + vctY * M_I[1] + vctZ * M_I[2]) * c0;
                ourArr(i, j, k, iy_) +=
                    (vctX * M_I[3] + vctY * M_I[4] + vctZ * M_I[5]) * c0;
                ourArr(i, j, k, iz_) +=
                    (vctX * M_I[6] + vctY * M_I[7] + vctZ * M_I[8]) * c0;
              }
        }
  }
}

//==========================================================
void Pic::update_E_rhs(double* rhs) {
  std::string nameFunc = "Pic::update_E_rhs";
  timing_func(nameFunc);

  MultiFab tempNode(nGrids[0], DistributionMap(0), 3, nGst);
  tempNode.setVal(0.0);
  MultiFab temp2Node(nGrids[0], DistributionMap(0), 3, nGst);
  temp2Node.setVal(0.0);

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    apply_BC(cellStatus, centerB[iLevTest], 0, centerB[iLevTest].nComp(),
             &Pic::get_center_B, iLevTest);
    apply_BC(nodeStatus, nodeB[iLevTest], 0, nodeB[iLevTest].nComp(),
             &Pic::get_node_B, iLevTest);
  }

  const Real* invDx = Geom(0).InvCellSize();

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    curl_center_to_node(centerB[iLevTest], tempNode, invDx);
  }

  MultiFab::Saxpy(temp2Node, -fourPI, jHat, 0, 0, temp2Node.nComp(),
                  temp2Node.nGrow());

  MultiFab::Add(temp2Node, tempNode, 0, 0, tempNode.nComp(), temp2Node.nGrow());

  temp2Node.mult(fsolver.theta * tc->get_dt());
  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    MultiFab::Add(temp2Node, nodeE[iLevTest], 0, 0, nodeE[iLevTest].nComp(),
                  temp2Node.nGrow());
  }

  convert_3d_to_1d(temp2Node, rhs);
}

//==========================================================
void Pic::update_B() {
  if (!solveEM)
    return;

  std::string nameFunc = "Pic::update_B";
  timing_func(nameFunc);
  MultiFab dB(cGrids[0], DistributionMap(0), 3, nGst);

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    curl_node_to_center(nodeEth[iLevTest], dB, Geom(0).InvCellSize());
    MultiFab::Saxpy(centerB[iLevTest], -tc->get_dt(), dB, 0, 0,
                    centerB[iLevTest].nComp(), centerB[iLevTest].nGrow());
    centerB[iLevTest].FillBoundary(Geom(0).periodicity());

    apply_BC(cellStatus, centerB[iLevTest], 0, centerB[iLevTest].nComp(),
             &Pic::get_center_B, iLevTest);

    average_center_to_node(centerB[iLevTest], nodeB[iLevTest]);
    nodeB[iLevTest].FillBoundary(Geom(iLevTest).periodicity());
    apply_BC(nodeStatus, nodeB[iLevTest], 0, nodeB[iLevTest].nComp(),
             &Pic::get_node_B, iLevTest);
  }
}

//==========================================================
void Pic::calc_smooth_coef() {
  std::string nameFunc = "Pic::calc_smooth_coef";
  timing_func(nameFunc);

  nodeSmoothCoef.setVal(coefWeakSmooth);

  if (fabs(coefStrongSmooth - coefWeakSmooth) < 1e-3)
    return;

  Real gamma = 5. / 3;
  for (MFIter mfi(nodePlasma[nSpecies]); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodePlasma[nSpecies][mfi];
    const Box& box = mfi.fabbox();
    const Array4<Real>& arr = fab.array();

    const Array4<Real>& arrCoef = nodeSmoothCoef[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          const Real rho = arr(i, j, k, iRho_);

          if (rho > 1e-99) {
            const Real uBulk =
                sqrt(pow(arr(i, j, k, iUx_), 2) + pow(arr(i, j, k, iUy_), 2) +
                     pow(arr(i, j, k, iUz_), 2)) /
                rho;

            const Real uth = sqrt(gamma / 3.0 *
                                  (arr(i, j, k, iPxx_) + arr(i, j, k, iPyy_) +
                                   arr(i, j, k, iPzz_)) /
                                  rho);

            const Real mach = uBulk / uth;

            if (mach > strongSmoothMach) {
              arrCoef(i, j, k) = coefStrongSmooth;
            } else if (mach > weakSmoothMach) {
              Real r0 =
                  (mach - weakSmoothMach) / (strongSmoothMach - weakSmoothMach);
              arrCoef(i, j, k) =
                  coefWeakSmooth + (coefStrongSmooth - coefWeakSmooth) * r0;
            }
          }
        }
  }

  smooth_multifab(nodeSmoothCoef, true, 0.5);
}

//==========================================================
void Pic::smooth_multifab(amrex::MultiFab& mf, bool useFixedCoef,
                          double coefIn) {
  std::string nameFunc = "Pic::smooth_multifab";
  timing_func(nameFunc);

  MultiFab mfOld(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrow());

  auto smooth_dir = [&](int iDir) {
    int dIdx[3] = { 0, 0, 0 };
    dIdx[iDir] = 1;

    MultiFab::Copy(mfOld, mf, 0, 0, mf.nComp(), mf.nGrow());

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.validbox();

      amrex::Array4<amrex::Real> const& arrE = mf[mfi].array();
      amrex::Array4<amrex::Real> const& arrTmp = mfOld[mfi].array();
      amrex::Array4<amrex::Real> const& arrCoef = nodeSmoothCoef[mfi].array();

      const auto lo = IntVect(bx.loVect());
      const auto hi = IntVect(bx.hiVect());

      for (int iVar = 0; iVar < mf.nComp(); iVar++)
        for (int k = lo[iz_]; k <= hi[iz_]; k++)
          for (int j = lo[iy_]; j <= hi[iy_]; j++)
            for (int i = lo[ix_]; i <= hi[ix_]; i++) {
              Real coef = coefIn;
              if (!useFixedCoef) {
                coef = arrCoef(i, j, k);
              }

              const Real weightSelf = 1 - coef;
              const Real WeightNei = coef / 2.0;

              const Real neiSum =
                  arrTmp(i - dIdx[ix_], j - dIdx[iy_], k - dIdx[iz_], iVar) +
                  arrTmp(i + dIdx[ix_], j + dIdx[iy_], k + dIdx[iz_], iVar);
              arrE(i, j, k, iVar) =
                  weightSelf * arrE(i, j, k, iVar) + WeightNei * neiSum;
            }
    }

    mf.FillBoundary(Geom(0).periodicity());
  };

  smooth_dir(ix_);
  smooth_dir(iy_);
  smooth_dir(iz_);
}

//==========================================================
void Pic::smooth_E(MultiFab& mfE) {
  if (!doSmoothE)
    return;

  std::string nameFunc = "Pic::smooth_E";
  timing_func(nameFunc);

  for (int icount = 0; icount < nSmoothE; icount++) {
    smooth_multifab(mfE);
  }
}

//==========================================================
void Pic::apply_BC(const iMultiFab& status, MultiFab& mf, const int iStart,
                   const int nComp, GETVALUE func, const int iLev) {
  std::string nameFunc = "Pic::apply_BC";
  timing_func(nameFunc);

  if (Geom(iLev).isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  bool useFloatBC = (func == nullptr);

  // BoxArray ba = mf.boxArray();
  BoxArray ba = convert(activeRegion, mf.boxArray().ixType());

  const IntVect& ngrow = mf.nGrowVect();
  if (Geom(iLev).Domain().bigEnd(iz_) == Geom(iLev).Domain().smallEnd(iz_)) {
    ba.grow(iz_, ngrow[iz_]);
  }

  if (useFloatBC) {
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.fabbox();

      //! if there are cells not in the valid + periodic grown box
      //! we need to fill them here
      if (!ba.contains(bx)) {
        amrex::Array4<amrex::Real> const& arr = mf[mfi].array();

        const Array4<const int>& statusArr = status[mfi].array();

        const auto lo = IntVect(bx.loVect());
        const auto hi = IntVect(bx.hiVect());

        for (int k = lo[iz_] + 1; k <= hi[iz_] - 1; k++)
          for (int j = lo[iy_] + 1; j <= hi[iy_] - 1; j++)
            for (int i = lo[ix_] + 1; i <= hi[ix_] - 1; i++)
              if (statusArr(i, j, k, 0) == iBoundary_) {
                bool isNeiFound = false;

                // Find the neighboring physical cell
                for (int kk = -1; kk <= 1; kk++)
                  for (int jj = -1; jj <= 1; jj++)
                    for (int ii = -1; ii <= 1; ii++) {
                      if (!isNeiFound &&
                          statusArr(i + ii, j + jj, k + kk, 0) != iBoundary_) {
                        isNeiFound = true;
                        for (int iVar = iStart; iVar < iStart + nComp; iVar++) {
                          arr(i, j, k, iVar) =
                              arr(i + ii, j + jj, k + kk, iVar);
                        }
                      }
                    }

                continue;
              }
      }
    }
  } else {
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.fabbox();

      //! if there are cells not in the valid + periodic grown box
      //! we need to fill them here
      if (!ba.contains(bx)) {
        amrex::Array4<amrex::Real> const& arr = mf[mfi].array();

        const Array4<const int>& statusArr = status[mfi].array();

        const auto lo = IntVect(bx.loVect());
        const auto hi = IntVect(bx.hiVect());

        for (int iVar = iStart; iVar < iStart + nComp; iVar++)
          for (int k = lo[iz_] + 1; k <= hi[iz_] - 1; k++)
            for (int j = lo[iy_]; j <= hi[iy_]; j++)
              for (int i = lo[ix_]; i <= hi[ix_]; i++)
                if (statusArr(i, j, k, 0) == iBoundary_) {
                  arr(i, j, k, iVar) =
                      (this->*func)(mfi, i, j, k, iVar - iStart, iLev);
                }

        continue;
      }
    }
  }
}

//==========================================================
Real Pic::calc_E_field_energy() {
  Real sum = 0;
  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    for (MFIter mfi(nodeE[iLevTest]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = nodeE[iLevTest][mfi];
      const Box& box = mfi.validbox();
      const Array4<Real>& arr = fab.array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      // Do not count the right edges.
      Real sumLoc = 0;
      for (int k = lo.z; k <= hi.z - 1; ++k)
        for (int j = lo.y; j <= hi.y - 1; ++j)
          for (int i = lo.x; i <= hi.x - 1; ++i) {
            sumLoc += pow(arr(i, j, k, ix_), 2) + pow(arr(i, j, k, iy_), 2) +
                      pow(arr(i, j, k, iz_), 2);
          }

      const auto& dx = Geom(0).CellSize();
      const Real coef = 0.5 * dx[ix_] * dx[iy_] * dx[iz_] / fourPI;
      sum += sumLoc * coef;
    }
    ParallelDescriptor::ReduceRealSum(sum,
                                      ParallelDescriptor::IOProcessorNumber());

    if (!ParallelDescriptor::IOProcessor())
      sum = 0;
  }
  return sum;
}

//==========================================================
Real Pic::calc_B_field_energy() {
  Real sum = 0;

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    for (MFIter mfi(centerB[iLevTest]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = centerB[iLevTest][mfi];
      const Box& box = mfi.validbox();
      const Array4<Real>& arr = fab.array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      Real sumLoc = 0;
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            sumLoc += pow(arr(i, j, k, ix_), 2) + pow(arr(i, j, k, iy_), 2) +
                      pow(arr(i, j, k, iz_), 2);
          }

      const auto& dx = Geom(0).CellSize();
      const Real coef = 0.5 * dx[ix_] * dx[iy_] * dx[iz_] / fourPI;
      sum += sumLoc * coef;
    }
  }
  ParallelDescriptor::ReduceRealSum(sum,
                                    ParallelDescriptor::IOProcessorNumber());

  if (!ParallelDescriptor::IOProcessor())
    sum = 0;

  return sum;
}

//==========================================================
void Pic::convert_1d_to_3d(const double* const p, amrex::MultiFab& MF) {
  std::string nameFunc = "Pic::convert_1d_to_3d";
  timing_func(nameFunc);

  bool isCenter = MF.ixType().cellCentered();

  MF.setVal(0.0);

  int iCount = 0;
  for (amrex::MFIter mfi(MF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.tilebox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const amrex::Array4<amrex::Real>& arr = MF[mfi].array();

    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    int iMin = lo.x, jMin = lo.y, kMin = lo.z;

    const auto& nodeArr = nodeShare[mfi].array();
    for (int iVar = 0; iVar < MF.nComp(); iVar++)
      for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
          for (int i = iMin; i <= iMax; ++i)
            if (isCenter || nodeArr(i, j, k) == iAssign_) {
              arr(i, j, k, iVar) = p[iCount++];
            }
  }
}

//==========================================================
void Pic::convert_3d_to_1d(const amrex::MultiFab& MF, double* const p) {
  std::string nameFunc = "Pic::convert_3d_to_1d";
  timing_func(nameFunc);

  bool isCenter = MF.ixType().cellCentered();

  int iCount = 0;
  for (amrex::MFIter mfi(MF, doTiling); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.tilebox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    const amrex::Array4<amrex::Real const>& arr = MF[mfi].array();

    // Avoid double counting the share edges.
    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    int iMin = lo.x, jMin = lo.y, kMin = lo.z;

    const auto& nodeArr = nodeShare[mfi].array();
    for (int iVar = 0; iVar < MF.nComp(); iVar++)
      for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
          for (int i = iMin; i <= iMax; ++i)
            if (isCenter || nodeArr(i, j, k) == iAssign_) {
              p[iCount++] = arr(i, j, k, iVar);
            }
  }
}

//==========================================================
void Pic::report_load_balance() {
  // This function report the min, max, and average of the local memory usage,
  // blocks, cells and particles among all the MPIs.

  std::string nameFunc = "Pic::monitor";
  timing_func(nameFunc);

  const int iMem_ = 0, iNCell_ = 1, iNBlk_ = 2, iNParts_ = 3, nLocal = 4;
  float localInfo[nLocal];

  int nProc = ParallelDescriptor::NProcs(), root = 0;
  float* globalInfo;
  if (ParallelDescriptor::MyProc() == root) {
    globalInfo = new float[nLocal * nProc];
  }

  std::vector<int> rc(nProc, nLocal), disp(nProc, 0);
  for (int i = 0; i < nProc; i++) {
    disp[i] = i * nLocal;
  }

  localInfo[iMem_] = (float)read_mem_usage();

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    localInfo[iNBlk_] = (float)centerB[iLevTest].local_size();
    localInfo[iNCell_] =
        (float)get_local_node_or_cell_number(centerB[iLevTest]);
  }

  {
    long npart = 0;
    for (auto& part : parts) {
      npart += part->TotalNumberOfParticles(false, true);
    }
    localInfo[iNParts_] = (float)npart;
  }

  ParallelDescriptor::Gatherv(localInfo, nLocal, globalInfo, rc, disp, root);

  if (ParallelDescriptor::MyProc() == root) {
    float maxVal[nLocal] = { 0 };
    float minVal[nLocal] = { 1e10, 1e10, 1e10, 1e10 };
    float avgVal[nLocal] = { 0 };
    int maxLoc[nLocal] = { 0 };
    for (int iProc = 0; iProc < nProc; iProc++)
      for (int iType = 0; iType < nLocal; iType++) {
        const float val = globalInfo[disp[iProc] + iType];
        if (val > maxVal[iType]) {
          maxVal[iType] = val;
          maxLoc[iType] = iProc;
        }

        if (val < minVal[iType])
          minVal[iType] = val;

        avgVal[iType] += val;
      }

    for (int iType = 0; iType < nLocal; iType++) {
      avgVal[iType] /= nProc;
    }

    const int nw = 9;
    AllPrint()
        << "==============" << printPrefix
        << "Load balance report==================\n"
        << "|    Value     |    Min   |   Avg    |    Max   |where(max)|\n"

        << "| Memory(MB)   |" << std::setw(nw) << (int)minVal[iMem_] << " |"
        << std::setw(nw) << (int)avgVal[iMem_] << " |" << std::setw(nw)
        << (int)maxVal[iMem_] << " |" << std::setw(nw) << maxLoc[iMem_]
        << " |\n"

        << "|# of Blocks   |" << std::setw(nw) << (int)minVal[iNBlk_] << " |"
        << std::setw(nw) << (int)avgVal[iNBlk_] << " |" << std::setw(nw)
        << (int)maxVal[iNBlk_] << " |" << std::setw(nw) << maxLoc[iNBlk_]
        << " |\n"

        << "|# of Cells    |" << std::setw(nw) << (int)minVal[iNCell_] << " |"
        << std::setw(nw) << (int)avgVal[iNCell_] << " |" << std::setw(nw)
        << (int)maxVal[iNCell_] << " |" << std::setw(nw) << maxLoc[iNCell_]
        << " |\n"

        << "|# of particles|" << std::setw(nw) << minVal[iNParts_] << " |"
        << std::setw(nw) << (int)avgVal[iNParts_] << " |" << std::setw(nw)
        << maxVal[iNParts_] << " |" << std::setw(nw) << maxLoc[iNParts_]
        << " |\n"

        << "============================================================"
        << std::endl;
  }

  if (ParallelDescriptor::MyProc() == root) {
    delete[] globalInfo;
  }
}

void Pic::charge_exchange() {
  if (!stateOH || !sourcePT2OH || !source)
    return;

  source->set_node_fluid_to_zero();

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->charge_exchange(tc->get_dt(), stateOH.get(), sourcePT2OH.get(),
                              source.get());
  }

  // 'source' is applied to generate new particles every step, so sum_boundary()
  // is called here to correct boundary nodes. Boundary nodes of 'sourcePT2OH'
  // should be corrected just before PT->OH coupling, instead of here.
  source->sum_boundary();
  source->convert_moment_to_velocity(true);

  // fill_source_particles();
}