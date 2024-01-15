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
  } else if (command == "#PARTICLEBOXBOUNDARY") {
    int iSpecies;
    std::string lo, hi;
    param.read_var("iSpecies", iSpecies);
    if (iSpecies >= pBCs.size()) {
      amrex::Abort("Error: wrong input or too may particle species.");
    }

    for (int i = 0; i < nDim; i++) {
      param.read_var("particleBoxBoundaryLo", lo);
      param.read_var("particleBoxBoundaryHi", hi);
      pBCs[iSpecies].lo[i] = pBCs[iSpecies].num_type(lo);
      pBCs[iSpecies].hi[i] = pBCs[iSpecies].num_type(hi);
    }
  } else if (command == "#MIDDLEPOINTSOURCE") {
    param.read_var("middlePointSource", middlePointSource);
  } else if (command == "#RANDOMPARTICLESLOCATION") {
    param.read_var("isParticleLocationRandom", isParticleLocationRandom);
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
    if (nDim == 3)
      param.read_var("npcelz", nPartPerCell[iz_]);
  } else if (command == "#SOURCEPARTICLES") {
    param.read_var("npcelx", nSourcePPC[ix_]);
    param.read_var("npcely", nSourcePPC[iy_]);
    param.read_var("npcelz", nSourcePPC[iz_]);
  } else if (command == "#KINETICSOURCE") {
    param.read_var("kineticSource", kineticSource);
  } else if (command == "#ELECTRON") {
    param.read_var("qom", qomEl);
  } else if (command == "#DISCRETIZE" || command == "#DISCRETIZATION") {
    param.read_var("theta", fsolver.theta);
    param.read_var("coefDiff", fsolver.coefDiff);
  } else if (command == "#SPLITJ") {
    param.read_var("splitJ", doSplitJ);
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
      param.read_var("maxWeightRatio", maxWeightRatio);
    }
  } else if (command == "#FASTMERGE") {
    param.read_var("fastMerge", fastMerge);
    if (fastMerge) {
      param.read_var("nMergeOld", nMergeOld);
      param.read_var("nMergeNew", nMergeNew);
      param.read_var("nMergeTry", nMergeTry);
      param.read_var("mergeRatioMax", mergeRatioMax);
    }
  } else if (command == "#ADAPTIVESOURCEPPC") {
    param.read_var("adaptiveSourcePPC", adaptiveSourcePPC);
  } else if (command == "#MERGELIGHT") {
    param.read_var("mergeLight", mergeLight);
    if (mergeLight) {
      param.read_var("mergePartRatioMax", mergePartRatioMax);
    }
  } else if (command == "#VACUUM") {
    param.read_var("vacuum", vacuum);
  } else if (command == "#PARTICLELEVRATIO") {
    param.read_var("particleLevRatio", pLevRatio);
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

    update_grid_status();
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
void Pic::distribute_arrays(const Vector<BoxArray>& cGridsOld) {

  // The last one is the sum of all species.
  if (nodePlasma.empty())
    nodePlasma.resize(nSpecies + 1);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    distribute_FabArray(centerB[iLev], cGrids[iLev], DistributionMap(iLev), 3,
                        nGst);
    distribute_FabArray(nodeB[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst);
    distribute_FabArray(nodeE[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst);
    distribute_FabArray(nodeEth[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst);

    distribute_FabArray(centerNetChargeOld[iLev], cGrids[iLev],
                        DistributionMap(iLev), 1, nGst);
    distribute_FabArray(centerNetChargeN[iLev], cGrids[iLev],
                        DistributionMap(iLev), 1, nGst);
    distribute_FabArray(centerNetChargeNew[iLev], cGrids[iLev],
                        DistributionMap(iLev), 1, nGst);

    distribute_FabArray(centerDivE[iLev], cGrids[iLev], DistributionMap(iLev),
                        1, nGst);

    distribute_FabArray(centerPhi[iLev], cGrids[iLev], DistributionMap(iLev), 1,
                        nGst);

    bool doMoveData = false;
    if (!useExplicitPIC) {
      distribute_FabArray(nodeMM[iLev], nGrids[iLev], DistributionMap(iLev), 1,
                          1, doMoveData);
    }

    distribute_FabArray(E0[iLev], nGrids[iLev], DistributionMap(iLev), 3, nGst,
                        doMoveData);

    distribute_FabArray(U0[iLev], nGrids[iLev], DistributionMap(iLev), 3, nGst,
                        doMoveData);

    distribute_FabArray(centerMM[iLev], cGrids[iLev], DistributionMap(iLev), 1,
                        nGst, doMoveData);

    distribute_FabArray(jHat[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst, doMoveData);

    for (auto& pl : nodePlasma) {
      if (pl.empty())
        pl.resize(n_lev_max());
      distribute_FabArray(pl[iLev], nGrids[iLev], DistributionMap(iLev),
                          nMoments, nGst, doMoveData);
    }
  }

  distribute_grid_arrays(cGridsOld);
}

//==========================================================
void Pic::pre_regrid() {
  if (!parts.empty()) {
    for (int i = 0; i < nSpecies; i++) {
      // Label the particles outside the OLD PIC region. It should be called
      // before active region is updated.
      parts[i]->label_particles_outside_active_region();
    }
  }
}

void Pic::post_regrid() {

  distribute_arrays(cGridsOld);

  {
    bool doMoveData = false;

    iTot = nSpecies;
    if (plasmaEnergy.empty()) {
      plasmaEnergy.resize(nSpecies + 1);
    }

    distribute_FabArray(nodeSmoothCoef, nGrids[0], DistributionMap(0), 1, nGst,
                        doMoveData);
  }
  //===========Move field data around end====================

  //--------------particles-----------------------------------
  if (parts.empty()) {
    solveEM = true;
    for (int i = 0; i < nSpecies; i++) {
      auto ptr = std::unique_ptr<PicParticles >(
          new PicParticles(this, fi, tc, i, fi->get_species_charge(i),
                          fi->get_species_mass(i), nPartPerCell, testCase));

      // If contains neutrals, assume this is OH-PT coupling, and do not solve
      // for EM fields.
      if (ptr->is_neutral())
        solveEM = false;

      //----- Set parameters------------
      if (particleMergeThreshold >= 0) {
        ptr->set_merge_threshold(particleMergeThreshold);
      }

      if (particleMergeBinBuffer >= 0) {
        ptr->set_merge_velocity_bin_buffer(particleMergeBinBuffer);
      }

      ptr->set_fast_merge(fastMerge, nMergeOld, nMergeNew, nMergeTry,
                          mergeRatioMax, mergeLight, mergePartRatioMax);

      const Real vacuumNO = vacuum * cProtonMassSI * 1e6 * fi->get_Si2NoRho();
      ptr->set_vacuum(vacuumNO);

      ptr->particle_lev_ratio(pLevRatio);

      ptr->set_bc(pBCs[i]);
      //----------------------------------

      parts.push_back(std::move(ptr));

      auto ptrSource = std::unique_ptr<PicParticles >(
          new PicParticles(this, fi, tc, i, fi->get_species_charge(i),
                          fi->get_species_mass(i), nPartPerCell, testCase));

      sourceParts.push_back(std::move(ptrSource));
    }
  } else {
    for (int i = 0; i < nSpecies; i++) {
      // Label the particles outside the NEW PIC region.
      parts[i]->label_particles_outside_active_region_general();

      parts[i]->redistribute_particles();
    }
  }
  //--------------particles-----------------------------------

  // This part does not really work for multi-level.
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    int n = get_local_node_or_cell_number(nodeE[iLev]);
    eSolver.init(n, nDim, nDim, matvec_E_solver);

    n = get_local_node_or_cell_number(centerDivE[iLev]);
    divESolver.init(n, 1, nDim, matvec_divE_accurate);
  }
}

//==========================================================
void Pic::fill_new_node_E() {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(nodeE[iLev]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = nodeE[iLev][mfi];
      const Box& box = mfi.validbox();
      const Array4<Real>& arrE = fab.array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      const auto& status = nodeStatus[iLev][mfi].array();

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_new(status(i, j, k))) {
              arrE(i, j, k, ix_) = fi->get_ex(mfi, i, j, k, iLev);
              arrE(i, j, k, iy_) = fi->get_ey(mfi, i, j, k, iLev);
              arrE(i, j, k, iz_) = fi->get_ez(mfi, i, j, k, iLev);
            }
          }
    }
  }
}

//==========================================================
void Pic::fill_new_node_B() {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(nodeB[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrB = nodeB[iLev][mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);
      const auto& status = nodeStatus[iLev][mfi].array();
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_new(status(i, j, k))) {
              arrB(i, j, k, ix_) = fi->get_bx(mfi, i, j, k, iLev);
              arrB(i, j, k, iy_) = fi->get_by(mfi, i, j, k, iLev);
              arrB(i, j, k, iz_) = fi->get_bz(mfi, i, j, k, iLev);
            }
          }
    }
  }
}

//==========================================================
void Pic::fill_new_center_B() {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& centerArr = centerB[iLev][mfi].array();
      const auto& nodeArr = nodeB[iLev][mfi].array();
      const auto& status = cellStatus[iLev][mfi].array();

      amrex::ParallelFor(
          box, centerB[iLev].nComp(),
          [&] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) noexcept {
            if (bit::is_new(status(i, j, k))) {
              centerArr(i, j, k, iVar) = 0;

              Box subBox(IntVect{ AMREX_D_DECL(i, j, k) },
                         IntVect{ AMREX_D_DECL(i + 1, j + 1, k + 1) });
              amrex::ParallelFor(subBox, [&] AMREX_GPU_DEVICE(int ii, int jj,
                                                              int kk) noexcept {
                const Real coef = (nDim == 2 ? 0.25 : 0.125);
                centerArr(i, j, k, iVar) += coef * nodeArr(ii, jj, kk, iVar);
              });
            }
          });
    }
  }
}

//==========================================================
void Pic::fill_E_B_fields() {
  if (!solveEM)
    return;

  fill_new_node_E();
  fill_new_node_B();
  fill_new_center_B();

  //-----Coarse (iLev=0) grid boundary/internal ghost cells are filled----

  nodeE[0].FillBoundary(Geom(0).periodicity());
  nodeB[0].FillBoundary(Geom(0).periodicity());
  centerB[0].FillBoundary(Geom(0).periodicity());
  apply_BC(nodeStatus[0], nodeB[0], 0, nDim, &Pic::get_node_B, 0);
  apply_BC(nodeStatus[0], nodeE[0], 0, nDim, &Pic::get_node_E, 0);
  apply_BC(cellStatus[0], centerB[0], 0, centerB[0].nComp(), &Pic::get_center_B,
           0);

  //-----Fine (iLev>0) grid boundary/internal ghost cells are filled----
  for (int iLev = 1; iLev <= finest_level; iLev++) {
    nodeE[iLev].FillBoundary();
    nodeB[iLev].FillBoundary();
    centerB[iLev].FillBoundary();

    fill_fine_lev_bny_from_coarse(nodeE[iLev - 1], nodeE[iLev], 0,
                                  nodeE[iLev - 1].nComp(), ref_ratio[iLev - 1],
                                  Geom(iLev - 1), Geom(iLev), node_status(iLev),
                                  amrex::node_bilinear_interp);

    fill_fine_lev_bny_from_coarse(nodeB[iLev - 1], nodeB[iLev], 0,
                                  nodeB[iLev - 1].nComp(), ref_ratio[iLev - 1],
                                  Geom(iLev - 1), Geom(iLev), node_status(iLev),
                                  amrex::node_bilinear_interp);

    fill_fine_lev_bny_from_coarse(
        centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
        ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), cell_status(iLev),
        amrex::cell_bilinear_interp);
  }
}

//==========================================================
void Pic::fill_particles() {
  inject_particles_for_new_cells();
  inject_particles_for_boundary_cells();
}

void Pic::fill_source_particles() {
  if (kineticSource)
    return;

  bool doSelectRegion = false;
#ifdef _PT_COMPONENT_
  doSelectRegion = (nSpecies == 4);
#endif
  for (int i = 0; i < nSpecies; i++) {
    parts[i]->add_particles_source(source, stateOH, tc->get_dt(), nSourcePPC,
                                   doSelectRegion, adaptiveSourcePPC);
  }
}

//==========================================================
void Pic::update_part_loc_to_half_stage() {
  std::string nameFunc = "Pic::update_part_loc_to_half_stage";

  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (int i = 0; i < nSpecies; i++) {
      parts[i]->update_position_to_half_stage(nodeEth[iLev], nodeB[iLev],
                                              tc->get_dt());
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
      if (maxWeightRatio > 1)
        parts[i]->limit_weight(maxWeightRatio, parts[i]->is_neutral());
      parts[i]->split(reSamplingLowLimit, parts[i]->is_neutral());
      parts[i]->merge(reSamplingHighLimit);
    }
  }
}

//==========================================================
void Pic::particle_mover() {
  std::string nameFunc = "Pic::mover";

  timing_func(nameFunc);

  // if (useExplicitPIC) {

  // MultiFab tmpE(nGrids[iLev], DistributionMap(iLev), 3, nGst);
  // // nodeE/nodeEth is at t_n/t_{n+1}, tmpE is at t_{n+0.5}
  // MultiFab::LinComb(tmpE, 0.5, nodeEth[iLev], 0, 0.5, nodeE[iLev], 0, 0,
  //                   nodeE[iLev].nComp(), nodeE[iLev].nGrow());
  // for (int i = 0; i < nSpecies; i++) {
  //   parts[i]->mover(tmpE, nodeB[iLev], iLev, tc->get_dt(),
  //                   tc->get_next_dt());
  // }

  // } else {

  Real dt = tc->get_dt();
  Real dtnext = tc->get_next_dt();

  if (middlePointSource) {
    dt *= 0.5;
    dtnext *= 0.5;
  }

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->mover(nodeEth, nodeB, dt, dtnext);
  }

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->redistribute_particles();
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

  for (int iLev = 0; iLev < n_lev(); iLev++) {

    jHat[iLev].setVal(0.0);

    if (!useExplicitPIC) {
      const RealMM mm0(0.0);
      nodeMM[iLev].setVal(mm0);
    }

    for (int i = 0; i < nSpecies; i++) {
      if (useExplicitPIC) {
        parts[i]->calc_jhat(jHat[iLev], nodeB[iLev], tc->get_dt());
      } else {
        parts[i]->calc_mass_matrix(nodeMM[iLev], jHat[iLev], nodeB[iLev],
                                   U0[iLev], tc->get_dt(), iLev);
      }
    }
    Real invVol = 1;
    for (int i = 0; i < nDim; i++) {
      invVol *= Geom(iLev).InvCellSize(i);
    }

    jHat[iLev].mult(invVol, 0, jHat[iLev].nComp(), jHat[iLev].nGrow());
    jHat[iLev].SumBoundary(Geom(iLev).periodicity());

    if (!useExplicitPIC) {
      nodeMM[iLev].SumBoundary(Geom(iLev).periodicity());
      nodeMM[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }

  for (int iLev = n_lev() - 2; iLev >= 0; iLev--) {
    sum_two_lev_interface_node(jHat[iLev], jHat[iLev + 1], 0,
                               jHat[iLev].nComp(), ref_ratio[iLev], Geom(iLev),
                               Geom(iLev + 1), node_status(iLev + 1));
  }

  for (int iLev = n_lev() - 2; iLev >= 0; iLev--) {
    sum_two_lev_interface_node(
        nodeMM[iLev], nodeMM[iLev + 1], 0, nodeMM[iLev].nComp(),
        ref_ratio[iLev], Geom(iLev), Geom(iLev + 1), node_status(iLev + 1));
  }

  // WARNING: interp_from_coarse_to_fine_for_domain_edge might be needed here
}

//==========================================================
void Pic::sum_moments(bool updateDt) {
  std::string nameFunc = "Pic::sum_moments";

  if (isGridEmpty)
    return;

  timing_func(nameFunc);

  plasmaEnergy[iTot] = 0;
  for (int i = 0; i < nSpecies; i++) {
    Real energy = parts[i]->sum_moments(nodePlasma[i], nodeB, tc->get_dt());
    plasmaEnergy[i] = energy;
    plasmaEnergy[iTot] += energy;
  }

  if (updateDt) {
    const auto& dx = Geom(0).CellSize();
    Real minDx = 1e99;
    for (int iDim = 0; iDim < nDim; iDim++) {
      if (minDx > dx[iDim])
        minDx = dx[iDim];
    }

    Real uMax = 0;
    if (tc->get_cfl() > 0 || doReport) {
      for (int i = 0; i < nSpecies; i++) {
        const int iLev = 0;
        Real uMaxSpecies =
            parts[i]->calc_max_thermal_velocity(nodePlasma[i][iLev]);
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

      if (PicParticles::particlePosition == NonStaggered) {
        tc->set_dt(dt);
      }
    }

    if (doReport) {
      if (PicParticles::particlePosition == Staggered) {
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

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    nodePlasma[nSpecies][iLev].setVal(0.0);
  }

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->convert_to_fluid_moments(nodePlasma[i]);
  }

  for (int i = 0; i < nSpecies; i++) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      // Index of 'nSpecies' represents the sum of all species.
      MultiFab::Add(nodePlasma[nSpecies][iLev], nodePlasma[i][iLev], 0, 0,
                    nMoments, nGst);
    }
  }

  isMomentsUpdated = true;
}

//==========================================================
void Pic::calc_cost_per_cell(BalanceStrategy balanceStrategy) {
  if (!isMomentsUpdated && balanceStrategy == BalanceStrategy::Particle) {
    sum_moments(false);
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    switch (balanceStrategy) {
      case BalanceStrategy::Cell: {
        cellCost[iLev].setVal(1.0);
        break;
      }
      case BalanceStrategy::Particle: {
        average_node_to_cellcenter(
            cellCost[iLev], 0, nodePlasma[nSpecies][iLev], iNum_,
            cellCost[iLev].nComp(), cellCost[iLev].nGrow());
        break;
      }
      case BalanceStrategy::Hybrid: {
        break;
      }
      default:
        break;
    }

    for (MFIter mfi(cellCost[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto lo = lbound(box);
      const auto hi = ubound(box);
      const Array4<Real>& cost = cellCost[iLev][mfi].array();
      const Array4<int const> status = cellStatus[iLev][mfi].array();

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_refined(status(i, j, k))) {
              cost(i, j, k) = 0;
            } else if (bit::is_domain_edge(status(i, j, k))) {
              // When calculating cost for each cell, the ghost cells are
              // excluded. However, ghost cells also take time to update (e.g.
              // fill boundary, launch and update boundary particles...).
              // Therefore, the cost of ghost cells is added to the cost of the
              // corresponding valid cells. The factor of 2 is just a guess.
              cost(i, j, k) *= 2;
            }

            if (balanceStrategy == BalanceStrategy::Particle) {
              // 1. The cells have been refined also allocated and use memory.
              // 2. It looks like these cells need calculations when
              // interpolating between levels.
              // 3. The number 10 is chosen by experience.
              cost(i, j, k) += 10;
            }
          }
    }
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
    // inside divE_correct_particle_position(). redistribute_particles() deletes
    // these particles. In order to get correct moments, re-inject particles in
    // the ghost cells.
    parts[i]->redistribute_particles();
  }
  inject_particles_for_boundary_cells();

  sum_to_center(false);
}

//==========================================================
void Pic::divE_correct_particle_position() {
  std::string nameFunc = "Pic::correct_position";

  timing_func(nameFunc);

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->divE_correct_position(centerPhi[0]);
  }
}

//==========================================================
void Pic::calculate_phi(LinearSolver& solver) {
  std::string nameFunc = "Pic::calculate_phi";

  timing_func(nameFunc);

  const int iLev = 0;
  MultiFab residual(cGrids[iLev], DistributionMap(iLev), 1, nGst);

  solver.reset(get_local_node_or_cell_number(centerDivE[iLev]));
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    div_node_to_center(nodeE[iLev], residual, Geom(iLev).InvCellSize());
  }

  Real coef = 1;
  if (PicParticles::particlePosition == Staggered) {
    coef = 1.0 / rhoTheta;
  }

  MultiFab::LinComb(residual, coef, residual, 0, -fourPI * coef,
                    centerNetChargeN[iLev], 0, 0, residual.nComp(),
                    residual.nGrow());

  convert_3d_to_1d(residual, solver.rhs, iLev);

  BL_PROFILE_VAR("Pic::phi_iterate", solve);
  solver.solve(doReport);
  BL_PROFILE_VAR_STOP(solve);

  convert_1d_to_3d(solver.xLeft, centerPhi[iLev], iLev);
  centerPhi[iLev].FillBoundary(Geom(iLev).periodicity());
}

//==========================================================
void Pic::divE_accurate_matvec(const double* vecIn, double* vecOut) {
  std::string nameFunc = "Pic::divE_matvec";
  timing_func(nameFunc);

  const int iLev = 0;

  zero_array(vecOut, divESolver.get_nSolve());

  MultiFab inMF(cGrids[iLev], DistributionMap(iLev), 1, nGst);

  convert_1d_to_3d(vecIn, inMF, iLev);
  inMF.FillBoundary(0, 1, IntVect(1), Geom(iLev).periodicity());

  MultiFab outMF(cGrids[iLev], DistributionMap(iLev), 1, nGst);
  outMF.setVal(0.0);

  for (MFIter mfi(inMF); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const Array4<Real>& lArr = outMF[mfi].array();
    const Array4<Real const>& rArr = inMF[mfi].array();
    const Array4<RealCMM>& mmArr = centerMM[iLev][mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i)

          for (int i2 = i - 1; i2 <= i + 1; i2++)
            for (int j2 = j - 1; j2 <= j + 1; j2++)
              for (int k2 = k - 1; k2 <= k + 1; k2++) {
                const int gp = (i2 - i + 1) * 9 + (j2 - j + 1) * 3 + k2 - k + 1;
                lArr(i, j, k) += rArr(i2, j2, k2) * mmArr(i, j, k)[gp];
              }
  }
  outMF.mult(fourPI * fourPI);
  convert_3d_to_1d(outMF, vecOut, iLev);
}

//==========================================================
void Pic::sum_to_center(bool isBeforeCorrection) {
  std::string nameFunc = "Pic::sum_to_center";

  timing_func(nameFunc);

  int iLev = 0;

  centerNetChargeNew[iLev].setVal(0.0);

  const RealCMM mm0(0.0);
  centerMM[iLev].setVal(mm0);

  bool doNetChargeOnly = !isBeforeCorrection;

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->sum_to_center(centerNetChargeNew[iLev], centerMM[iLev],
                            doNetChargeOnly);
  }

  if (!doNetChargeOnly) {
    centerMM[iLev].SumBoundary(Geom(iLev).periodicity());
  }

  centerNetChargeNew[iLev].SumBoundary(Geom(iLev).periodicity());

  apply_BC(cellStatus[iLev], centerNetChargeNew[iLev], 0,
           centerNetChargeNew[iLev].nComp(), &Pic::get_zero, iLev);

  if (PicParticles::particlePosition == NonStaggered) {
    MultiFab::Copy(centerNetChargeN[iLev], centerNetChargeNew[iLev], 0, 0,
                   centerNetChargeN[iLev].nComp(),
                   centerNetChargeN[iLev].nGrow());
  } else {
    MultiFab::LinComb(
        centerNetChargeN[iLev], 1 - rhoTheta, centerNetChargeOld[iLev], 0,
        rhoTheta, centerNetChargeNew[iLev], 0, 0,
        centerNetChargeN[iLev].nComp(), centerNetChargeN[iLev].nGrow());

    if (!isBeforeCorrection) {
      MultiFab::Copy(centerNetChargeOld[iLev], centerNetChargeNew[iLev], 0, 0,
                     centerNetChargeOld[iLev].nComp(),
                     centerNetChargeOld[iLev].nGrow());
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

  if (PicParticles::particlePosition == NonStaggered) {
    update_part_loc_to_half_stage();
  }

  calc_mass_matrix();

  update_E();

  particle_mover();

  // Calling re_sampling after particle mover so that all the particles outside
  // the domain have been deleted.
  re_sampling();

  charge_exchange();

  if (middlePointSource) {
    particle_mover();
  }

  if (source) {
    fill_source_particles();
  }

  inject_particles_for_boundary_cells();

  isMomentsUpdated = false;

  update_B();

  if (doCorrectDivE) {
    divE_correction();
  }

  tc->set_dt(tc->get_next_dt());

#ifdef _PC_COMPONENT_
  //  For PT simulations, moments are only useful for output. So, there is no
  //  need to call sum_moments() for every step.
  sum_moments(true);
#endif

  if (doReport) {
    Real tEnd = second();
    Real nPoint = activeRegion.d_numPts();
    int nProc = ParallelDescriptor::NProcs();
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
void Pic::update_U0_E0() {
  std::string nameFunc = "Pic::update_U0_E0";
  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    U0[iLev].setVal(0.0);
    E0[iLev].setVal(0.0);
    for (MFIter mfi(U0[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const Array4<Real>& arrU = U0[iLev][mfi].array();
      const Array4<const Real>& arrMoments =
          nodePlasma[nSpecies][iLev][mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            const Real rho = arrMoments(i, j, k, iRho_);
            if (rho > 1e-99) {
              const Real invRho = 1. / rho;
              for (int iu = iUx_; iu <= iUz_; iu++)
                arrU(i, j, k, iu - iUx_) = arrMoments(i, j, k, iu) * invRho;
            }
          }
    }

    U0[iLev].FillBoundary(Geom(iLev).periodicity());

    const int nSmooth = 5;

    for (int i = 0; i < nSmooth; i++)
      smooth_multifab(U0[iLev], true, 0.5);

    // MultiFab::Copy(tempNode3, nodeB, 0, 0, nodeB.nComp(), nodeB.nGrow());
    // tempNode3.FillBoundary(Geom(0).periodicity());
    // for (int i = 0; i < 5; i++)
    //   smooth_multifab(tempNode3, true, 0.5);

    for (MFIter mfi(U0[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrU = U0[iLev][mfi].array();
      const Array4<Real>& arrE = E0[iLev][mfi].array();
      const Array4<Real>& arrB = nodeB[iLev][mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            const Real& bx = arrB(i, j, k, ix_);
            const Real& by = arrB(i, j, k, iy_);
            const Real& bz = arrB(i, j, k, iz_);

            const Real& ux = arrU(i, j, k, ix_);
            const Real& uy = arrU(i, j, k, iy_);
            const Real& uz = arrU(i, j, k, iz_);

            arrE(i, j, k, ix_) = -uy * bz + uz * by;
            arrE(i, j, k, iy_) = -uz * bx + ux * bz;
            arrE(i, j, k, iz_) = -ux * by + uy * bx;
          }
    }

    E0[iLev].FillBoundary(Geom(iLev).periodicity());

    for (int i = 0; i < nSmooth; i++)
      smooth_multifab(E0[iLev], true, 0.5);

    //
    // print_MultiFab(nodeU0, "nodeU0", 1);
    // print_MultiFab(nodeE0, "nodeE0", 1);
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

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    MultiFab::Copy(nodeEth[iLev], nodeE[iLev], 0, 0, nodeE[iLev].nComp(),
                   nodeE[iLev].nGrow());
    apply_BC(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
             &Pic::get_center_B, iLev);
  }
  const Real dt = tc->get_dt();
  RealVect dt2dx;
  for (int i = 0; i < nDim; i++) {
    dt2dx[i] = dt * Geom(0).InvCellSize(i);
  }
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    curl_center_to_node(centerB[iLev], nodeE[iLev], dt2dx.begin());
    MultiFab::Saxpy(nodeE[iLev], -fourPI * dt, jHat[iLev], 0, 0,
                    nodeE[iLev].nComp(), nodeE[iLev].nGrow());

    MultiFab::Add(nodeE[iLev], nodeEth[iLev], 0, 0, nodeE[iLev].nComp(),
                  nodeE[iLev].nGrow());

    nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_BC(nodeStatus[iLev], nodeE[iLev], 0, nDim, &Pic::get_node_E, iLev);
  }
}

//==========================================================
void Pic::update_E_impl() {
  std::string nameFunc = "Pic::update_E_impl";

  timing_func(nameFunc);

  if (doSplitJ)
    update_U0_E0();

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    eSolver.reset(get_local_node_or_cell_number(nodeE[iLev]));

    update_E_rhs(eSolver.rhs, iLev);

    convert_3d_to_1d(nodeE[iLev], eSolver.xLeft, iLev);

    update_E_matvec(eSolver.xLeft, eSolver.matvec, iLev, false);

    for (int i = 0; i < eSolver.get_nSolve(); i++) {
      eSolver.rhs[i] -= eSolver.matvec[i];
      eSolver.xLeft[i] = 0;
    }

    if (doReport)
      Print() << "\n-------" << printPrefix << " E solver ------------------"
              << std::endl;

    BL_PROFILE_VAR("Pic::E_iterate", eSolver);
    eSolver.solve(iLev, doReport);
    BL_PROFILE_VAR_STOP(eSolver);

    nodeEth[iLev].setVal(0.0);
    convert_1d_to_3d(eSolver.xLeft, nodeEth[iLev], iLev);
    nodeEth[iLev].SumBoundary(Geom(iLev).periodicity());
    nodeEth[iLev].FillBoundary(Geom(iLev).periodicity());
    MultiFab::Add(nodeEth[iLev], nodeE[iLev], 0, 0, nodeEth[iLev].nComp(),
                  nGst);

    MultiFab::LinComb(nodeE[iLev], -(1.0 - fsolver.theta) / fsolver.theta,
                      nodeE[iLev], 0, 1. / fsolver.theta, nodeEth[iLev], 0, 0,
                      nodeE[iLev].nComp(), nGst);

    if (iLev == 0) {

      apply_BC(nodeStatus[iLev], nodeE[iLev], 0, nDim, &Pic::get_node_E, iLev);
      apply_BC(nodeStatus[iLev], nodeEth[iLev], 0, nDim, &Pic::get_node_E,
               iLev);

    } else {

      fill_fine_lev_bny_from_coarse(
          nodeE[iLev - 1], nodeE[iLev], 0, nodeE[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          amrex::node_bilinear_interp);

      fill_fine_lev_bny_from_coarse(
          nodeEth[iLev - 1], nodeEth[iLev], 0, nodeEth[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          amrex::node_bilinear_interp);
    }

    if (doSmoothE) {
      calc_smooth_coef();
      smooth_E(nodeEth[iLev]);
      smooth_E(nodeE[iLev]);
    }
    div_node_to_center(nodeE[iLev], centerDivE[iLev], Geom(iLev).InvCellSize());
  }
}

//==========================================================
void Pic::update_E_matvec(const double* vecIn, double* vecOut, int iLev,
                          const bool useZeroBC) {
  std::string nameFunc = "Pic::E_matvec";
  timing_func(nameFunc);

  zero_array(vecOut, eSolver.get_nSolve());

  MultiFab vecMF(nGrids[iLev], DistributionMap(iLev), 3, nGst);
  vecMF.setVal(0.0);

  MultiFab matvecMF(nGrids[iLev], DistributionMap(iLev), 3, 1);
  matvecMF.setVal(0.0);

  MultiFab tempCenter3(cGrids[iLev], DistributionMap(iLev), 3, nGst);

  MultiFab tempNode3(nGrids[iLev], DistributionMap(iLev), 3, nGst);
  tempNode3.setVal(0.0);

  MultiFab tempCenter1(cGrids[iLev], DistributionMap(iLev), 1, nGst);

  convert_1d_to_3d(vecIn, vecMF, iLev);

  // The right side edges should be filled in.
  vecMF.SumBoundary(Geom(iLev).periodicity());

  // M*E needs ghost cell information.
  vecMF.FillBoundary(Geom(iLev).periodicity());

  if (useZeroBC) {
    // The boundary nodes would not be filled in by convert_1d_3d. So, there
    // is not need to apply zero boundary conditions again here.
  } else {
    // Even after apply_BC(), the outmost layer node E is still
    // unknow. See FluidInterface::calc_current for detailed explaniation.
    if (iLev == 0) {
      apply_BC(nodeStatus[iLev], vecMF, 0, nDim, &Pic::get_node_E, iLev);
    } else {
      fill_fine_lev_bny_from_coarse(
          nodeE[iLev - 1], vecMF, 0, centerB[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          amrex::node_bilinear_interp);
    }
  }

  lap_node_to_node(vecMF, matvecMF, DistributionMap(iLev), Geom(iLev),
                   cellStatus[iLev]);

  Real delt2 = pow(fsolver.theta * tc->get_dt(), 2);
  matvecMF.mult(-delt2);

  { // grad(divE)
    div_node_to_center(vecMF, centerDivE[iLev], Geom(iLev).InvCellSize());

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

      div_center_to_center(tempCenter3, tempCenter1, Geom(iLev).InvCellSize());

      tempCenter1.FillBoundary(0, 1, IntVect(1), Geom(iLev).periodicity());

      // 1) The outmost boundary layer of tempCenter3 is not accurate.
      // 2) The 2 outmost boundary layers (all ghosts if there are 2 ghost
      // cells) of tempCenter1 are not accurate
      apply_BC(cellStatus[iLev], tempCenter1, 0, tempCenter1.nComp(),
               &Pic::get_zero, iLev);

      MultiFab::LinComb(centerDivE[iLev], 1 - fsolver.coefDiff,
                        centerDivE[iLev], 0, fsolver.coefDiff, tempCenter1, 0,
                        0, 1, 1);
    }

    grad_center_to_node(centerDivE[iLev], tempNode3, Geom(iLev).InvCellSize());

    tempNode3.mult(delt2);
    MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(),
                  matvecMF.nGrow());
  }

  tempNode3.setVal(0);
  update_E_M_dot_E(vecMF, tempNode3, iLev);

  MultiFab::Add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(), 0);

  MultiFab::Add(matvecMF, vecMF, 0, 0, matvecMF.nComp(), 0);

  convert_3d_to_1d(matvecMF, vecOut, iLev);
}

//==========================================================
void Pic::update_E_M_dot_E(const MultiFab& inMF, MultiFab& outMF, int iLev) {
  std::string nameFunc = "Pic::update_E_M_dot_E";
  timing_func(nameFunc);

  outMF.setVal(0.0);
  Real c0 = fourPI * fsolver.theta * tc->get_dt();
  for (MFIter mfi(outMF); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const Array4<Real const>& inArr = inMF[mfi].array();
    const Array4<Real>& ourArr = outMF[mfi].array();
    const Array4<RealMM>& mmArr = nodeMM[iLev][mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          auto& data0 = mmArr(i, j, k);
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
void Pic::update_E_rhs(double* rhs, int iLev) {
  std::string nameFunc = "Pic::update_E_rhs";
  timing_func(nameFunc);

  MultiFab tempNode(nGrids[iLev], DistributionMap(iLev), 3, nGst);
  tempNode.setVal(0.0);
  MultiFab temp2Node(nGrids[iLev], DistributionMap(iLev), 3, nGst);
  temp2Node.setVal(0.0);

  apply_BC(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
           &Pic::get_center_B, iLev);
  apply_BC(nodeStatus[iLev], nodeB[iLev], 0, nodeB[iLev].nComp(),
           &Pic::get_node_B, iLev);

  const Real* invDx = Geom(iLev).InvCellSize();

  curl_center_to_node(centerB[iLev], tempNode, invDx);

  MultiFab::Saxpy(temp2Node, -fourPI, jHat[iLev], 0, 0, temp2Node.nComp(),
                  temp2Node.nGrow());

  MultiFab::Add(temp2Node, tempNode, 0, 0, tempNode.nComp(), temp2Node.nGrow());

  temp2Node.mult(fsolver.theta * tc->get_dt());
  MultiFab::Add(temp2Node, nodeE[iLev], 0, 0, nodeE[iLev].nComp(),
                temp2Node.nGrow());

  tempNode.setVal(0.0);
  update_E_M_dot_E(E0[iLev], tempNode, iLev);
  MultiFab::Add(temp2Node, tempNode, 0, 0, tempNode.nComp(), tempNode.nGrow());

  convert_3d_to_1d(temp2Node, rhs, iLev);
}

//==========================================================
void Pic::update_B() {
  if (!solveEM)
    return;

  std::string nameFunc = "Pic::update_B";
  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    MultiFab dB(cGrids[iLev], DistributionMap(iLev), 3, nGst);
    curl_node_to_center(nodeEth[iLev], dB, Geom(iLev).InvCellSize());
    MultiFab::Saxpy(centerB[iLev], -tc->get_dt(), dB, 0, 0,
                    centerB[iLev].nComp(), centerB[iLev].nGrow());
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_BC(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
               &Pic::get_center_B, iLev);

    } else {
      fill_fine_lev_bny_from_coarse(
          centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          amrex::cell_bilinear_interp);
    }

    average_center_to_node(centerB[iLev], nodeB[iLev]);
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {

      apply_BC(nodeStatus[iLev], nodeB[iLev], 0, nodeB[iLev].nComp(),
               &Pic::get_node_B, iLev);

    } else {

      fill_fine_lev_bny_from_coarse(
          nodeB[iLev - 1], nodeB[iLev], 0, nodeB[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          amrex::node_bilinear_interp);
    }
  }
}

//==========================================================
void Pic::calc_smooth_coef() {
  std::string nameFunc = "Pic::calc_smooth_coef";
  timing_func(nameFunc);

  nodeSmoothCoef.setVal(coefWeakSmooth);

  if (fabs(coefStrongSmooth - coefWeakSmooth) < 1e-3)
    return;

  int iLev = 0;
  Real gamma = 5. / 3;
  for (MFIter mfi(nodePlasma[nSpecies][iLev]); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodePlasma[nSpecies][iLev][mfi];
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
void Pic::smooth_multifab(MultiFab& mf, bool useFixedCoef, double coefIn) {
  std::string nameFunc = "Pic::smooth_multifab";
  timing_func(nameFunc);

  MultiFab mfOld(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrow());

  auto smooth_dir = [&](int iDir) {
    int dIdx[3] = { 0, 0, 0 };
    dIdx[iDir] = 1;

    MultiFab::Copy(mfOld, mf, 0, 0, mf.nComp(), mf.nGrow());

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.validbox();

      Array4<Real> const& arrE = mf[mfi].array();
      Array4<Real> const& arrTmp = mfOld[mfi].array();
      Array4<Real> const& arrCoef = nodeSmoothCoef[mfi].array();

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
        Array4<Real> const& arr = mf[mfi].array();

        const Array4<const int>& statusArr = status[mfi].array();

        const auto lo = IntVect(bx.loVect());
        const auto hi = IntVect(bx.hiVect());

        for (int k = lo[iz_] + 1; k <= hi[iz_] - 1; k++)
          for (int j = lo[iy_] + 1; j <= hi[iy_] - 1; j++)
            for (int i = lo[ix_] + 1; i <= hi[ix_] - 1; i++)
              if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
                bool isNeiFound = false;

                // Find the neighboring physical cell
                for (int kk = -1; kk <= 1; kk++)
                  for (int jj = -1; jj <= 1; jj++)
                    for (int ii = -1; ii <= 1; ii++) {
                      if (!isNeiFound && !bit::is_lev_boundary(statusArr(
                                             i + ii, j + jj, k + kk, 0))) {
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
        Array4<Real> const& arr = mf[mfi].array();

        const Array4<const int>& statusArr = status[mfi].array();

        const auto lo = IntVect(bx.loVect());
        const auto hi = IntVect(bx.hiVect());

        for (int iVar = iStart; iVar < iStart + nComp; iVar++)
          for (int k = lo[iz_] + 1; k <= hi[iz_] - 1; k++)
            for (int j = lo[iy_]; j <= hi[iy_]; j++)
              for (int i = lo[ix_]; i <= hi[ix_]; i++)
                if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
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
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(nodeE[iLev]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = nodeE[iLev][mfi];
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

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = centerB[iLev][mfi];
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
void Pic::convert_1d_to_3d(const double* const p, MultiFab& MF, int iLev) {
  std::string nameFunc = "Pic::convert_1d_to_3d";
  timing_func(nameFunc);

  bool isCenter = MF.ixType().cellCentered();

  MF.setVal(0.0);

  int iCount = 0;
  for (MFIter mfi(MF, doTiling); mfi.isValid(); ++mfi) {
    const Box& box = mfi.tilebox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);
    const Array4<Real>& arr = MF[mfi].array();

    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    int iMin = lo.x, jMin = lo.y, kMin = lo.z;

    const auto& nodeArr = nodeStatus[iLev][mfi].array();
    for (int iVar = 0; iVar < MF.nComp(); iVar++)
      for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
          for (int i = iMin; i <= iMax; ++i)
            if (isCenter || bit::is_owner(nodeArr(i, j, k))) {
              arr(i, j, k, iVar) = p[iCount++];
            }
  }
}

//==========================================================
void Pic::convert_3d_to_1d(const MultiFab& MF, double* const p, int iLev) {
  std::string nameFunc = "Pic::convert_3d_to_1d";
  timing_func(nameFunc);

  bool isCenter = MF.ixType().cellCentered();

  int iCount = 0;
  for (MFIter mfi(MF, doTiling); mfi.isValid(); ++mfi) {
    const Box& box = mfi.tilebox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);
    const Array4<Real const>& arr = MF[mfi].array();

    // Avoid double counting the share edges.
    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    int iMin = lo.x, jMin = lo.y, kMin = lo.z;

    const auto& nodeArr = nodeStatus[iLev][mfi].array();
    for (int iVar = 0; iVar < MF.nComp(); iVar++)
      for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
          for (int i = iMin; i <= iMax; ++i)
            if (isCenter || bit::is_owner(nodeArr(i, j, k))) {
              p[iCount++] = arr(i, j, k, iVar);
            }
  }
}

//==========================================================
void Pic::report_load_balance(bool doReportSummary, bool doReportDetail) {
  // This function report the min, max, and average of the local memory usage,
  // blocks, cells and particles among all the MPIs.
  if (!doReportSummary && !doReportDetail)
    return;

  std::string nameFunc = "Pic::monitor";
  timing_func(nameFunc);

  int iNBlk_ = 0, iNCell_ = 1, iNParts_ = 2, iMem_ = 3 * (n_lev() + 1),
      nLocal = iMem_ + 1;

  Vector<float> localInfo(nLocal, 0);

  int nProc = ParallelDescriptor::NProcs();

  Vector<int> rc(nProc, nLocal), disp(nProc, 0);
  for (int i = 0; i < nProc; i++) {
    disp[i] = i * nLocal;
  }

  localInfo[iMem_] = (float)read_mem_usage();

  const int iBt = n_lev() * 3 + iNBlk_;
  const int iCt = n_lev() * 3 + iNCell_;
  const int iPt = n_lev() * 3 + iNParts_;
  localInfo[iBt] = 0;
  localInfo[iCt] = 0;
  localInfo[iPt] = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const int iB = iLev * 3 + iNBlk_;
    const int iC = iLev * 3 + iNCell_;
    const int iP = iLev * 3 + iNParts_;

    localInfo[iB] = (float)centerB[iLev].local_size();
    localInfo[iC] = (float)get_local_node_or_cell_number(centerB[iLev]);

    localInfo[iP] = 0;
    for (auto& part : parts) {
      localInfo[iP] += (float)part->NumberOfParticlesAtLevel(iLev, false, true);
    }

    localInfo[iBt] += localInfo[iB];
    localInfo[iCt] += localInfo[iC];
    localInfo[iPt] += localInfo[iP];
  }

  Vector<float> globalInfo;
  if (ParallelDescriptor::IOProcessor()) {
    globalInfo.resize(nLocal * nProc);
  }

  int iop = ParallelDescriptor::IOProcessorNumber();
  ParallelDescriptor::Gatherv(localInfo.dataPtr(), nLocal, globalInfo.data(),
                              rc, disp, iop);

  if (ParallelDescriptor::IOProcessor()) {

    if (doReportSummary) {

      Vector<float> maxVal(nLocal, 0);
      Vector<float> minVal(nLocal, 1e10);
      Vector<float> avgVal(nLocal, 0);
      Vector<int> maxLoc(nLocal, 0);

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

      printf("\n===============================Load balance "
             "report=============================\n");
      printf(
          "|     Value          |      Min      |     Avg      |      Max     "
          "|where(max)|\n");

      Vector<std::string> varType = {
        "|Blocks # of",
        "|Cells  # of",
        "|Parts  # of",
        "|Memory(MB)          |",
      };

      for (int iLev = 0; iLev <= n_lev(); iLev++) {
        for (int i = iNBlk_; i <= iNParts_; i++) {
          int idx = iLev * 3 + i;
          if (iLev < n_lev()) {
            printf("%s lev  %d %s %13.1f |%13.1f |%13.1f | %9d|\n",
                   varType[i].c_str(), iLev, " |", minVal[idx], avgVal[idx],
                   maxVal[idx], maxLoc[idx]);
          } else {
            printf("%s all levs| %13.1f |%13.1f |%13.1f | %9d|\n",
                   varType[i].c_str(), minVal[idx], avgVal[idx], maxVal[idx],
                   maxLoc[idx]);
          }
        }
        printf(
            "----------------------------------------------------------------"
            "---------------\n");
      }
      printf("%s %13.1f |%13.1f |%13.1f | %9d|\n", varType[3].c_str(),
             minVal[iMem_], avgVal[iMem_], maxVal[iMem_], maxLoc[iMem_]);

      printf("================================================================"
             "===============\n\n");
    }

    if (doReportDetail) {
      printf("\n");
      printf("=======================Work load of each MPI "
             "rank====================");
      for (int iLev = 1; iLev <= n_lev(); iLev++) {
        printf("=============================================");
      }
      printf("\n");

      printf("rank    |   Memory(MB) ");
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        printf("| Blocks lev %d |  Cells lev %d |  Parts lev %d ", iLev, iLev,
               iLev);
      }
      printf("| Blocks all   |  Cells all   |  Parts all   |\n");

      for (int rank = 0; rank < nProc; rank++) {
        float* info = globalInfo.data() + rank * nLocal;
        printf("%6d  |%13.1f ", rank, info[iMem_]);
        for (int iLev = 0; iLev <= n_lev(); iLev++) {
          printf("|%13.1f |%13.1f |%13.1f ", info[iLev * 3 + iNBlk_],
                 info[iLev * 3 + iNCell_], info[iLev * 3 + iNParts_]);
        }
        printf("|\n");
      }

      printf("================================================================="
             "====");
      for (int iLev = 1; iLev <= n_lev(); iLev++) {
        printf("=============================================");
      }
      printf("\n\n");
    }
  }
}

void Pic::charge_exchange() {
  timing_func("Pic::charge_exchange");

  if (!stateOH || !sourcePT2OH || !source)
    return;

  if (!kineticSource)
    source->set_node_fluid_to_zero();

  bool doSelectRegion = false;
#ifdef _PT_COMPONENT_
  doSelectRegion = (nSpecies == 4);
#endif

  int sourcePPC = 1;
  for (int i = 0; i < nDim; i++) {
    sourcePPC *= nSourcePPC[i];
  }

  for (int i = 0; i < nSpecies; i++) {
    parts[i]->charge_exchange(tc->get_dt(), stateOH, sourcePT2OH, source,
                              kineticSource, sourceParts, doSelectRegion,
                              sourcePPC);
  }

  if (kineticSource) {
    for (int i = 0; i < nSpecies; i++) {
      parts[i]->add_source_particles(sourceParts[i], nSourcePPC,
                                     adaptiveSourcePPC);
      sourceParts[i]->clearParticles();
    }

  } else {
    // 'source' is applied to generate new particles every step, so
    // sum_boundary() is called here to correct boundary nodes. Boundary nodes
    // of 'sourcePT2OH' should be corrected just before PT->OH coupling, instead
    // of here.
    source->sum_boundary();

#ifdef _PT_COMPONENT_
    bool doRegionSplit = (nSpecies == 4);
    if (doRegionSplit) {
      source->sum_to_single_source();
    }
#endif

    source->convert_moment_to_velocity(true, false);
  }
}
