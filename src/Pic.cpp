#include <algorithm>
#include <math.h>

#include <AMReX_Algorithm.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include "GridUtility.h"
#include "LinearSolver.h"
#include "Pic.h"
#include "SimDomains.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
void Pic::read_param(const std::string& command, ReadParam& param) {

  if (command == "#PIC") {
    param.read_var("usePIC", usePIC);
  } else if (command == "#SOLVEEM") {
    param.read_var("solveEM", solveEM);
  } else if (command == "#PARTMODE") {
    std::string s;
    param.read_var("partMode", s);
    if (s == "SEP")
      pMode = PartMode::SEP;
    else if (s == "PIC")
      pMode = PartMode::PIC;
    else
      Abort("Error: wrong input for partMode.");
  } else if (command == "#PARTICLEBOXBOUNDARY") {
    int iSpecies;
    std::string lo, hi;
    param.read_var("iSpecies", iSpecies);
    if (iSpecies >= pBCs.size()) {
      Abort("Error: wrong input or too may particle species.");
    }

    for (int i = 0; i < nDim; ++i) {
      param.read_var("particleBoxBoundaryLo", lo);
      param.read_var("particleBoxBoundaryHi", hi);
      pBCs[iSpecies].lo[i] = pBCs[iSpecies].num_type(lo);
      pBCs[iSpecies].hi[i] = pBCs[iSpecies].num_type(hi);
    }
  } else if (command == "#BFIELDBOXBOUNDARY") {
    std::string lo, hi;
    for (int i = 0; i < nDim; ++i) {
      param.read_var("BoxBoundaryLo", lo);
      param.read_var("BoxBoundaryHi", hi);
      bcBField.lo[i] = bcBField.num_type(lo);
      bcBField.hi[i] = bcBField.num_type(hi);
    }
  } else if (command == "#RANDOMPARTICLESLOCATION") {
    param.read_var("isParticleLocationRandom", isParticleLocationRandom);
  } else if (command == "#CONSTANTPPV") {
    param.read_var("isPPVconstant", isPPVconstant);
  } else if (command == "#PRESPLITTING") {
    param.read_var("doPreSplitting", doPreSplitting);
  } else if (command == "#OVERRIDEPRESSUREANISOTROPY") {
    param.read_var("doOverridePressureAnisotropy",
                   doOverridePressureAnisotropy);
    if (doOverridePressureAnisotropy) {
      int nspec = 2;
      param.read_var("nSpecies", nspec);
      for (int i = 0; i < nspec; ++i) {
        amrex::Real ratio = 1.0;
        param.read_var("initialAnisotropyRatio", ratio);
        initialAnisotropyRatios.push_back(ratio);
      }
    }
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
  } else if (command == "#COMOVING") {
    param.read_var("solveFieldInCoMov", solveFieldInCoMov);
    param.read_var("solvePartInCoMov", solvePartInCoMov);
    param.read_var("nSmoothBackGroundU", nSmoothBackGroundU);
  } else if (command == "#UPWINDE") {
    param.read_var("useUpwindE", useUpwindE);
    param.read_var("limiterThetaE", limiterThetaE);
    param.read_var("cMaxE", cMaxE);
  } else if (command == "#SMOOTHE") {
    param.read_var("doSmoothE", doSmoothE);
    if (doSmoothE) {
      param.read_var("nSmoothE", nSmoothE);
    }
  } else if (command == "#SMOOTHJ") {
    param.read_var("doSmoothJ", doSmoothJ);
    if (doSmoothJ) {
      param.read_var("nSmoothJ", nSmoothJ);
      param.read_var("coefSmoothJ", coefSmoothJ);
    }
  } else if (command == "#UPWINDB") {
    param.read_var("useUpwindB", useUpwindB);
    param.read_var("theta", limiterTheta);
    param.read_var("Isotropy", smoothBIso);
    if (useUpwindB) {
      useHyperbolicCleaning = true;
    }
  } else if (command == "#DIVB") {
    param.read_var("useHyperbolicCleaning", useHyperbolicCleaning);
    if (useHyperbolicCleaning) {
      param.read_var("hypDecay", hypDecay);
    }
  } else if (command == "#RESAMPLING") {
    param.read_var("doReSampling", doReSampling);
    if (doReSampling) {
      param.read_var("reSamplingLowLimit", reSamplingLowLimit);
      param.read_var("reSamplingHighLimit", reSamplingHighLimit);
      param.read_var("maxWeightRatio", maxWeightRatio);
    }
  } else if (command == "#FASTMERGE") {
    param.read_var("fastMerge", pInfo.fastMerge);
    if (pInfo.fastMerge) {
      param.read_var("nMergeOld", pInfo.nPartCombine);
      param.read_var("nMergeNew", pInfo.nPartNew);
      param.read_var("nMergeTry", pInfo.nMergeTry);
      param.read_var("mergeRatioMax", pInfo.mergeRatioMax);
    }
  } else if (command == "#ADAPTIVESOURCEPPC") {
    param.read_var("adaptiveSourcePPC", adaptiveSourcePPC);
  } else if (command == "#MERGELIGHT") {
    param.read_var("mergeLight", pInfo.mergeLight);
    if (pInfo.mergeLight) {
      param.read_var("mergePartRatioMax", pInfo.mergePartRatioMax);
    }
  } else if (command == "#VACUUM") {
    param.read_var("vacuum", pInfo.vacuumIO);
  } else if (command == "#PARTICLELEVRATIO") {
    param.read_var("particleLevRatio", pInfo.pLevRatio);
  } else if (command == "#OHION") {
    // The units are assumed to be:
    // r: AU
    // rho: amu/cc
    // T: K
    // U: km/s
    param.read_var("rAnalytic", ionOH.rAnalytic);
    param.read_var("doGetFromOH", ionOH.doGetFromOH);

    if (!ionOH.doGetFromOH) {
      param.read_var("rCutoff", ionOH.rCutoff);
      param.read_var("swRho", ionOH.swRho);
      param.read_var("swT", ionOH.swT);
      param.read_var("swU", ionOH.swU);
    }
  } else if (command == "#SUPID") {
    int n = 0;
    param.read_var("nSpecies", n);
    for (int i = 0; i < n; ++i) {
      int supid;
      param.read_var("supid", supid);
      supIDs.push_back(supid);
    }
  } else if (command == "#TESTCASE") {
    std::string testcase;
    param.read_var("testCase", testcase);
    if (testcase == "TwoStream") {
      testCase = TwoStream;
    } else if (testcase == "tophat") {
      testCase = TopHat;
      nPartPerCell = IntVect::Zero;
    } else if (testcase == "lightwave") {
      testCase = LightWave;
      nPartPerCell = IntVect::Zero;
    }
  } else if (command == "#SELECTPARTICLE") {
    param.read_var("doSelectParticle", doSelectParticle);
    if (doSelectParticle) {
      param.read_var("selectParticleInputFile", selectParticleInputFile);
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

  if (isPPVconstant || doPreSplitting) {
    SetTargetPPC(2);
    isTargetPPCDefined = true;
    for (int i = 0; i < nSpecies; i++) {
      parts[i]->set_is_target_ppc_defined(isTargetPPCDefined);
    }
  }
  if (initEM)
    fill_E_B_fields();

  if (testCase == LightWave) {
    fill_lightwaves(48.0);
  }

  if (usePIC) {
    fill_particles();
    sum_moments(true);
    if (finest_level == 0) {
      sum_to_center(false);
    } else if (doCorrectDivE) {
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        sum_to_center_new(false, iLev);
      }
    }
  }

  doNeedFillNewCell = false;
}

//==========================================================
void Pic::distribute_arrays(const Vector<BoxArray>& cGridsOld) {

  // The last one is the sum of all species.
  if (nodePlasma.empty()) {
    nodePlasma.resize(nSpecies + 1);
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (reportParticleQuality) {
      distribute_FabArray(particleQuality[iLev], cGrids[iLev],
                          DistributionMap(iLev), 18, 0);
    }
    distribute_FabArray(targetPPC[iLev], cGrids[iLev], DistributionMap(iLev), 1,
                        nGst);
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

    distribute_FabArray(divB[iLev], cGrids[iLev], DistributionMap(iLev), 3,
                        nGst, doMoveData);
    distribute_FabArray(hypPhi[iLev], cGrids[iLev], DistributionMap(iLev), 3,
                        nGst, doMoveData);
    distribute_FabArray(dBdt[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst, doMoveData);

    distribute_FabArray(eBg[iLev], nGrids[iLev], DistributionMap(iLev), 3, nGst,
                        doMoveData);

    distribute_FabArray(uBg[iLev], nGrids[iLev], DistributionMap(iLev), 3, nGst,
                        doMoveData);

    distribute_FabArray(mMach[iLev], nGrids[iLev], DistributionMap(iLev), 1,
                        nGst, doMoveData);

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
    for (int i = 0; i < nSpecies; ++i) {
      // Label the particles outside the OLD PIC region. It should be called
      // before active region is updated.
      parts[i]->label_particles_outside_active_region();
    }
  }
}

void Pic::post_regrid() {

  distribute_arrays(cGridsOld);

  {
    iTot = nSpecies;
    if (plasmaEnergy.empty()) {
      plasmaEnergy.resize(nSpecies + 1);
    }
  }
  //===========Move field data around end====================

  //--------------particles-----------------------------------
  if (parts.empty()) {
    for (int i = 0; i < nSpecies; ++i) {
      auto ptr = std::unique_ptr<PicParticles>(new PicParticles(
          this, fi, tc, i, fi->get_species_charge(i), fi->get_species_mass(i),
          nPartPerCell, pMode, testCase));

      //----- Set parameters------------
      ptr->set_info(pInfo);

      ptr->set_ion_fluid(ionOH);

      ptr->set_bc(pBCs[i]);

      if (!supIDs.empty())
        ptr->set_sup_id(supIDs[i]);
      //----------------------------------

      parts.push_back(std::move(ptr));

      auto ptrSource = std::unique_ptr<PicParticles>(new PicParticles(
          this, fi, tc, i, fi->get_species_charge(i), fi->get_species_mass(i),
          nPartPerCell, pMode, testCase));

      sourceParts.push_back(std::move(ptrSource));
    }
  } else {
    for (int i = 0; i < nSpecies; ++i) {
      // Label the particles outside the NEW PIC region.
      parts[i]->label_particles_outside_active_region_general();

      parts[i]->redistribute_particles();
    }
  }
  //--------------particles-----------------------------------

  // This part does not really work for multi-level.
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    int n = get_local_node_or_cell_number(nodeE[iLev]);
    eSolver.init(n, nDim3, nDim, matvec_E_solver);

    n = get_local_node_or_cell_number(centerDivE[iLev]);
    divESolver.init(n, 1, nDim, matvec_divE_accurate);
  }
}

//==========================================================
void Pic::fill_new_node_E() {
  {
    Real xL = 0, xR = 0;
    if (testCase == TopHat) {
      xL = 0.75 * Geom(0).ProbLo()[ix_] + 0.25 * Geom(0).ProbHi()[ix_];
      xR = 0.75 * Geom(0).ProbHi()[ix_] + 0.25 * Geom(0).ProbLo()[ix_];
    }

    int iLev = 0;
    for (MFIter mfi(nodeE[iLev]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = nodeE[iLev][mfi];
      const Box& box = mfi.validbox();
      const Array4<Real>& arrE = fab.array();
      const auto& status = nodeStatus[iLev][mfi].array();

      ParallelFor(box, [&](int i, int j, int k) {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };
        if (bit::is_new(status(ijk))) {
          if (testCase == TopHat) {
            const Real x =
                Geom(iLev).CellCenter(i, ix_) - 0.5 * Geom(iLev).CellSize(ix_);
            if (x > xL && x < xR) {
              arrE(ijk, iy_) = 1;
            }
          } else {
            arrE(ijk, ix_) = fi->get_ex(mfi, ijk, iLev);
            arrE(ijk, iy_) = fi->get_ey(mfi, ijk, iLev);
            arrE(ijk, iz_) = fi->get_ez(mfi, ijk, iLev);
          }
        }
      });
    }
  }
  if (finest_level > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_new_from_coarse(
          nodeE[iLev - 1], nodeE[iLev], 0, nodeE[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          node_bilinear_interp);
    }
  }
}

//==========================================================
void Pic::fill_new_node_B() {
  {
    Real xL = 0, xR = 0;
    if (testCase == TopHat) {
      xL = 0.75 * Geom(0).ProbLo()[ix_] + 0.25 * Geom(0).ProbHi()[ix_];
      xR = 0.75 * Geom(0).ProbHi()[ix_] + 0.25 * Geom(0).ProbLo()[ix_];
    }

    int iLev = 0;
    for (MFIter mfi(nodeB[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrB = nodeB[iLev][mfi].array();
      const auto& status = nodeStatus[iLev][mfi].array();

      ParallelFor(box, [&](int i, int j, int k) {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };
        if (bit::is_new(status(ijk))) {
          if (testCase == TopHat) {
            const Real x =
                Geom(iLev).CellCenter(i, ix_) - 0.5 * Geom(iLev).CellSize(ix_);
            if (x > xL && x < xR) {
              arrB(ijk, iz_) = 1;
            }
          } else {
            arrB(ijk, ix_) = fi->get_bx(mfi, ijk, iLev);
            arrB(ijk, iy_) = fi->get_by(mfi, ijk, iLev);
            arrB(ijk, iz_) = fi->get_bz(mfi, ijk, iLev);
          }
        }
      });
    }
  }

  if (finest_level > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_new_from_coarse(
          nodeB[iLev - 1], nodeB[iLev], 0, nodeB[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          node_bilinear_interp);
    }
  }
}

//==========================================================
void Pic::fill_new_center_B() {
  {
    int iLev = 0;
    for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& centerArr = centerB[iLev][mfi].array();
      const auto& nodeArr = nodeB[iLev][mfi].array();
      const auto& status = cellStatus[iLev][mfi].array();

      ParallelFor(
          box, centerB[iLev].nComp(), [&](int i, int j, int k, int iVar) {
            IntVect ijk = { AMREX_D_DECL(i, j, k) };

            if (bit::is_new(status(ijk))) {
              centerArr(ijk, iVar) = 0;

              Box subBox(ijk, ijk + 1);
              ParallelFor(subBox, [&](int ii, int jj, int kk) {
                const Real coef = (nDim == 2 ? 0.25 : 0.125);
                centerArr(ijk, iVar) += coef * nodeArr(ii, jj, kk, iVar);
              });
            }
          });
    }
  }
  if (finest_level > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_new_from_coarse(
          centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          cell_bilinear_interp);
    }
  }
}

//==========================================================
void Pic::fill_E_B_fields() {
  fill_new_node_E();
  fill_new_node_B();
  fill_new_center_B();

  //-----Coarse (iLev=0) grid boundary/internal ghost cells are filled----

  nodeE[0].FillBoundary(Geom(0).periodicity());
  nodeB[0].FillBoundary(Geom(0).periodicity());
  centerB[0].FillBoundary(Geom(0).periodicity());
  apply_BC(nodeStatus[0], nodeB[0], 0, nDim3, &Pic::get_node_B, 0, &bcBField);
  apply_BC(nodeStatus[0], nodeE[0], 0, nDim3, &Pic::get_node_E, 0);
  apply_BC(cellStatus[0], centerB[0], 0, centerB[0].nComp(), &Pic::get_center_B,
           0, &bcBField);

  //-----Fine (iLev>0) grid boundary/internal ghost cells are filled----
  for (int iLev = 1; iLev <= finest_level; iLev++) {
    nodeE[iLev].FillBoundary();
    nodeB[iLev].FillBoundary();
    centerB[iLev].FillBoundary();

    fill_fine_lev_bny_from_coarse(nodeE[iLev - 1], nodeE[iLev], 0,
                                  nodeE[iLev - 1].nComp(), ref_ratio[iLev - 1],
                                  Geom(iLev - 1), Geom(iLev), node_status(iLev),
                                  node_bilinear_interp);

    fill_fine_lev_bny_from_coarse(nodeB[iLev - 1], nodeB[iLev], 0,
                                  nodeB[iLev - 1].nComp(), ref_ratio[iLev - 1],
                                  Geom(iLev - 1), Geom(iLev), node_status(iLev),
                                  node_bilinear_interp);

    fill_fine_lev_bny_from_coarse(
        centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
        ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), cell_status(iLev),
        cell_bilinear_interp);
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
  for (int i = 0; i < nSpecies; ++i) {
    parts[i]->add_particles_source(source, stateOH, tc->get_dt(), nSourcePPC,
                                   doSelectRegion, adaptiveSourcePPC);
  }
}

//==========================================================
void Pic::update_part_loc_to_half_stage() {
  std::string nameFunc = "Pic::update_part_loc_to_half_stage";

  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (int i = 0; i < nSpecies; ++i) {
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
    for (int i = 0; i < nSpecies; ++i) {
      if (!doPreSplitting) {
        if (maxWeightRatio > 1) {
          parts[i]->limit_weight(maxWeightRatio, parts[i]->is_neutral());
        }
        parts[i]->split(reSamplingLowLimit, parts[i]->is_neutral());
        parts[i]->merge(reSamplingHighLimit);
      } else {
        if (maxWeightRatio > 1) {
          parts[i]->limit_weight_new(maxWeightRatio, parts[i]->is_neutral());
        }
        parts[i]->split_new(reSamplingLowLimit, parts[i]->is_neutral());
        parts[i]->merge_new(reSamplingHighLimit);
      }
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
  // for (int i = 0; i < nSpecies; ++i) {
  //   parts[i]->mover(tmpE, nodeB[iLev], iLev, tc->get_dt(),
  //                   tc->get_next_dt());
  // }

  // } else {

  Real dt = tc->get_dt();
  Real dtnext = tc->get_next_dt();

  for (int i = 0; i < nSpecies; ++i) {
    parts[i]->mover(nodeEth, nodeB, eBg, uBg, dt, dtnext, solvePartInCoMov);
  }

  for (int i = 0; i < nSpecies; ++i) {
    parts[i]->redistribute_particles();
  }
}

//==========================================================
void Pic::calc_mass_matrix() {
  std::string nameFunc = "Pic::calc_mass_matrix";

  if (isGridEmpty)
    return;

  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {

    jHat[iLev].setVal(0.0);

    if (!useExplicitPIC) {
      const RealMM mm0(0.0);
      nodeMM[iLev].setVal(mm0);
    }

    for (int i = 0; i < nSpecies; ++i) {
      if (useExplicitPIC) {
        parts[i]->calc_jhat(jHat[iLev], nodeB[iLev], tc->get_dt());
      } else {
        parts[i]->calc_mass_matrix(nodeMM[iLev], jHat[iLev], nodeB[iLev],
                                   uBg[iLev], tc->get_dt(), iLev,
                                   solveFieldInCoMov);
      }
    }
    Real invVol = 1;
    for (int i = 0; i < nDim; ++i) {
      invVol *= Geom(iLev).InvCellSize(i);
    }

    jHat[iLev].mult(invVol, 0, jHat[iLev].nComp(), jHat[iLev].nGrow());
    jHat[iLev].SumBoundary(Geom(iLev).periodicity());
    jHat[iLev].FillBoundary(Geom(iLev).periodicity());

    if (doSmoothJ) {
      for (int icount = 0; icount < nSmoothJ; icount++) {
        smooth_multifab(jHat[iLev], iLev, icount % 2 + 1, coefSmoothJ);
      }
    }

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
void Pic::calc_mass_matrix_amr() {
  std::string nameFunc = "Pic::calc_mass_matrix";

  if (isGridEmpty)
    return;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    nodeMM[iLev].setVal(0.0);
    jHat[iLev].setVal(0.0);
  }
  if (skipMassMatrix)
    return;

  timing_func(nameFunc);
  //////////////////////////////////////////////////////////////////////
  amrex::Vector<amrex::Vector<amrex::MultiFab> > jhc;
  amrex::Vector<amrex::MultiFab> jhf;
  amrex::Vector<amrex::Vector<UMultiFab<RealMM> > > nmmc;
  amrex::Vector<UMultiFab<RealMM> > nmmf;
  jhc.resize(n_lev());
  jhf.resize(n_lev());
  nmmc.resize(n_lev());
  nmmf.resize(n_lev());
  for (int iLev = 1; iLev < n_lev(); iLev++) {
    jhc[iLev].resize(iLev);
    nmmc[iLev].resize(iLev);
  }
  for (int iLev = 1; iLev < n_lev(); iLev++) {
    BoxArray bac = nodeB[iLev].boxArray();
    for (int i = iLev - 1; i >= 0; i--) {
      bac.coarsen(ref_ratio[iLev]);
      jhc[iLev][i].define(bac, nodeB[iLev].DistributionMap(), 3, 0);
      nmmc[iLev][i].define(bac, nodeB[iLev].DistributionMap(),
                           nodeMM[iLev].nComp(), 0);
      jhc[iLev][i].setVal(0.0);
      nmmc[iLev][i].setVal(0.0);
    }
  }
  for (int iLev = 0; iLev < n_lev() - 1; iLev++) {
    BoxArray baf = nodeB[iLev].boxArray();
    baf.refine(ref_ratio[iLev]);
    jhf[iLev].define(baf, nodeB[iLev].DistributionMap(), 3, 0);
    nmmf[iLev].define(baf, nodeB[iLev].DistributionMap(), nodeMM[iLev].nComp(),
                      0);
    jhf[iLev].setVal(0.0);
    nmmf[iLev].setVal(0.0);
  }
  //////////////////////////////////////////////////////////////////////
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->calc_mass_matrix_amr(nodeMM[iLev], nmmc, nmmf, jHat[iLev], jhc,
                                     jhf, nodeB[iLev], uBg[iLev], tc->get_dt(),
                                     iLev, solveFieldInCoMov, cellStatus);
    }
  }
  //////////////////////////////////////////////////////////////////////
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    jHat[iLev].SumBoundary(Geom(iLev).periodicity());
    nodeMM[iLev].SumBoundary(Geom(iLev).periodicity());
  }
  Real invVol[finest_level + 1];
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    invVol[iLev] = 1.0;
    for (int i = 0; i < nDim; ++i) {
      invVol[iLev] *= Geom(iLev).InvCellSize(i);
    }
  }

  for (int iLev = finest_level - 1; iLev >= 0; iLev--) {
    for (int i = finest_level; i > iLev; i--) {
      jHat[iLev].ParallelAdd(jhc[i][iLev]);
      nmmc[i][iLev].mult(invVol[iLev] / invVol[i]);
      nodeMM[iLev].ParallelAdd(nmmc[i][iLev]);
    }
  }
  for (int iLev = finest_level; iLev > 0; iLev--) {
    jHat[iLev].ParallelAdd(jhf[iLev - 1]);
    nmmf[iLev - 1].mult(invVol[iLev] / invVol[iLev - 1]);
    nodeMM[iLev].ParallelAdd(nmmf[iLev - 1]);
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    Real invVol = 1;
    for (int i = 0; i < nDim; ++i) {
      invVol *= Geom(iLev).InvCellSize(i);
    }
    jHat[iLev].mult(invVol, 0, jHat[iLev].nComp(), jHat[iLev].nGrow());
    jHat[iLev].FillBoundary(Geom(iLev).periodicity());
    nodeMM[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  //////////// Fill empty nodeMM elements////////////////
  // for (int iLev = 0; iLev < n_lev(); iLev++) {
  //   for (MFIter mfi(nodeMM[iLev]); mfi.isValid(); ++mfi) {
  //     // Finalize the mass matrix calculation.
  //     const Box box = mfi.validbox();
  //     const auto lo = lbound(box);
  //     const auto hi = ubound(box);

  //     Array4<RealMM> const& mmArr = nodeMM[iLev][mfi].array();

  //     // We only need the mass matrix on the physical nodes. But the first
  //     // layer
  //     // of the ghost nodes may contributes to the physical nodes below
  //     (ghost
  //     // node constributes as a sender). So, we need the '-1' and '+1' staff.
  //     const int iMin = lo.x - 1, jMin = lo.y - 1,
  //               kMin = nDim > 2 ? lo.z - 1 : 0;
  //     const int iMax = hi.x + 1, jMax = hi.y + 1,
  //               kMax = nDim > 2 ? hi.z + 1 : 0;

  //     int gps, gpr; // gp_send, gp_receive
  //     for (int k1 = kMin; k1 <= kMax; k1++)
  //       for (int j1 = jMin; j1 <= jMax; j1++)
  //         for (int i1 = iMin; i1 <= iMax; i1++) {
  //           const int kp = 2;
  //           const int kr = nDim > 2 ? k1 + kp - 1 : 0;
  //           if (kr > kMax || kr < kMin)
  //             continue;
  //           auto& datas0 = mmArr(i1, j1, k1);
  //           for (int jp = 0; jp < 3; jp++) {
  //             const int jr = j1 + jp - 1;
  //             if (jr > jMax || jr < jMin)
  //               continue;
  //             const int jpr = 2 - jp;
  //             for (int ip = 0; ip < 3; ip++) {
  //               const int ir = i1 + ip - 1;
  //               if (ir > iMax || ir < iMin)
  //                 continue;
  //               const int ipr = 2 - ip;
  //               gpr = jpr * 3 + ipr;
  //               gps = 18 + jp * 3 + ip; // gps = kp*9+jp*3+kp

  //               Real* const datar = &(mmArr(ir, jr, kr)[gpr * 9]);
  //               const Real* const datas = &(datas0[gps * 9]);
  //               for (int idx = 0; idx < 9; idx++) {
  //                 datar[idx] = datas[idx];
  //               } // idx
  //             } // kp
  //           } // jp
  //         } // k1
  //   }
  // }
  ///////////////////////////////////////////////////////////////
}
//==========================================================
void Pic::calc_mass_matrix_new() {
  std::string nameFunc = "Pic::calc_mass_matrix";
  if (isGridEmpty)
    return;
  timing_func(nameFunc);

  //////////////////////////////////////////////////////////////////////
  amrex::Vector<amrex::MultiFab> jhc;
  amrex::Vector<amrex::MultiFab> jhf;
  amrex::Vector<UMultiFab<RealMM> > nmmc;
  amrex::Vector<UMultiFab<RealMM> > nmmf;
  jhc.resize(n_lev());
  jhf.resize(n_lev());
  nmmc.resize(n_lev());
  nmmf.resize(n_lev());
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (iLev > 0) {
      BoxArray bac = nodeB[iLev].boxArray();
      bac.coarsen(ref_ratio[iLev]);
      jhc[iLev].define(bac, nodeB[iLev].DistributionMap(), 3, 0);
      nmmc[iLev].define(bac, nodeB[iLev].DistributionMap(),
                        nodeMM[iLev].nComp(), 0);
    }
    if (iLev < finest_level) {
      BoxArray baf = nodeB[iLev].boxArray();
      baf.refine(ref_ratio[iLev]);
      jhf[iLev].define(baf, nodeB[iLev].DistributionMap(), 3, 0);

      nmmf[iLev].define(baf, nodeB[iLev].DistributionMap(),
                        nodeMM[iLev].nComp(), 0);
    }
  }
  //////////////////////////////////////////////////////////////////////
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    nodeMM[iLev].setVal(0.0);
    jHat[iLev].setVal(0.0);
    jhc[iLev].setVal(0.0);
    nmmc[iLev].setVal(0.0);
    jhf[iLev].setVal(0.0);
    nmmf[iLev].setVal(0.0);
  }
  //////////////////////////////////////////////////////////////////////
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->calc_mass_matrix_new_optimized(
          nodeMM, nmmc[iLev], nmmf[iLev], jHat, jhc[iLev], jhf[iLev], nodeB,
          uBg, tc->get_dt(), iLev, solveFieldInCoMov, nodeStatus, cellStatus);
    }
  }
  //////////////////////////////////////////////////////////////////////

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    jHat[iLev].SumBoundary(Geom(iLev).periodicity());
    nodeMM[iLev].SumBoundary(Geom(iLev).periodicity());
    if (iLev > 0) {
      jhc[iLev].SumBoundary();
      nmmc[iLev].SumBoundary();
    }
    if (iLev < finest_level) {
      jhf[iLev].SumBoundary();
      nmmf[iLev].SumBoundary();
    }
  }

  amrex::MultiFab tmp;
  UMultiFab<RealMM> tmpMM;
  for (int iLev = 0; iLev < n_lev(); iLev++) {

    if (iLev > 0) {
      tmp.define(jHat[iLev - 1].boxArray(), jHat[iLev - 1].DistributionMap(),
                 jHat[iLev - 1].nComp(), 0);
      tmpMM.define(nodeMM[iLev - 1].boxArray(),
                   nodeMM[iLev - 1].DistributionMap(), nodeMM[iLev - 1].nComp(),
                   0);
      tmp.setVal(0.0);
      tmpMM.setVal(0.0);
      tmp.ParallelCopy(jhc[iLev]);
      tmpMM.ParallelCopy(nmmc[iLev]);
      MultiFab::Add(jHat[iLev - 1], tmp, 0, 0, 3, 0);
      amrex::Add(nodeMM[iLev - 1], tmpMM, 0, 0, nodeMM[iLev - 1].nComp(), 0);
    }

    if (iLev < finest_level) {
      tmp.define(jHat[iLev + 1].boxArray(), jHat[iLev + 1].DistributionMap(),
                 jHat[iLev + 1].nComp(), 0);
      tmpMM.define(nodeMM[iLev + 1].boxArray(),
                   nodeMM[iLev + 1].DistributionMap(), nodeMM[iLev + 1].nComp(),
                   0);

      tmp.setVal(0.0);
      tmpMM.setVal(0.0);

      tmp.ParallelCopy(jhf[iLev]);
      tmpMM.ParallelCopy(nmmf[iLev]);

      MultiFab::Add(jHat[iLev + 1], tmp, 0, 0, 3, 0);
      amrex::Add(nodeMM[iLev + 1], tmpMM, 0, 0, nodeMM[iLev + 1].nComp(), 0);
    }
  }

  //////////////////////////////////////////////////////////////////////

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    Real invVol = 1;
    for (int i = 0; i < nDim; ++i) {
      invVol *= Geom(iLev).InvCellSize(i);
    }

    jHat[iLev].mult(invVol, 0, jHat[iLev].nComp(), jHat[iLev].nGrow());

    if (!useExplicitPIC) {
      jHat[iLev].FillBoundary(Geom(iLev).periodicity());
      nodeMM[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }
  // for (int iLev = 0; iLev < finest_level; iLev++) {
  //   fill_fine_lev_edge_from_coarse(jHat[iLev], jHat[iLev + 1], 0,
  //                                  jHat[iLev].nComp(), ref_ratio[iLev],
  //                                  Geom(iLev), Geom(iLev + 1),
  //                                  node_status(iLev + 1),
  //                                  node_bilinear_interp);
  // }
  // average_down_nodal(jHat[1], jHat[0],ref_ratio[0]);
}
//==========================================================
void Pic::sum_moments(bool updateDt) {
  std::string nameFunc = "Pic::sum_moments";
  bool isPT = false;
#ifdef _PT_COMPONENT_
  isPT = true;
#endif
  if (isGridEmpty)
    return;

  timing_func(nameFunc);

  plasmaEnergy[iTot] = 0;
  for (int i = 0; i < nSpecies; ++i) {
    Real energy = 0.0;
    if (finest_level == 0 || isPT) {
      energy = parts[i]->sum_moments(nodePlasma[i], nodeB, tc->get_dt());
    } else {
      energy = parts[i]->sum_moments_new(nodePlasma[i], nodeB, tc->get_dt(),
                                         nodeStatus);
    }
    plasmaEnergy[i] = energy;
    plasmaEnergy[iTot] += energy;
  }

  if (updateDt) {
    amrex::Vector<amrex::Real> uMax(n_lev());
    amrex::Vector<amrex::Real> dxMin(n_lev());
    amrex::Vector<amrex::Real> dtMax(n_lev());
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      const auto& dx = Geom(iLev).CellSize();
      dxMin[iLev] = min(AMREX_D_DECL(dx[ix_], dx[iy_], dx[iz_]));

      if (tc->get_cfl() > 0 || doReport) {
        uMax[iLev] = 0.0;
        for (int i = 0; i < nSpecies; ++i) {
          Real uMaxSpecies =
              parts[i]->calc_max_thermal_velocity(nodePlasma[i][iLev]);
          ParallelDescriptor::ReduceRealMax(uMaxSpecies);

          if (doReport) {
            Print() << printPrefix << std::setprecision(5) << "lev " << iLev
                    << " Species " << i << ": max(uth) = " << uMaxSpecies
                    << std::endl;
          }

          if (uMaxSpecies > uMax[iLev]) {
            uMax[iLev] = uMaxSpecies;
          }

          if (testCase == TopHat || testCase == LightWave) {
            uMax[iLev] = 1.0;
          }
        }

        dtMax[iLev] = dxMin[iLev] / uMax[iLev];
      }
    }

    if (tc->get_cfl() > 0) {
      Real dt0 = *std::min_element(dtMax.begin(), dtMax.end());
      Real dt = tc->get_cfl() * dt0;
      tc->set_next_dt(dt);

      if (tc->get_dt() < 0) {
        tc->set_dt(dt);
      }
    }

    if (doReport) {
      Print() << printPrefix << std::setprecision(5)
              << "dt = " << tc->get_dt_si()
              << " dtNext = " << tc->get_next_dt_si() << std::endl;

      for (int iLev = 0; iLev < n_lev(); iLev++) {
        Print() << printPrefix << std::setprecision(5) << "iLev = " << iLev
                << " : CFL(dtNext) = " << tc->get_next_dt() / dtMax[iLev]
                << std::endl;
      }
    }
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    nodePlasma[nSpecies][iLev].setVal(0.0);
  }

  for (int i = 0; i < nSpecies; ++i) {
    parts[i]->convert_to_fluid_moments(nodePlasma[i]);
  }

  for (int i = 0; i < nSpecies; ++i) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      // Index of 'nSpecies' represents the sum of all species.
      MultiFab::Add(nodePlasma[nSpecies][iLev], nodePlasma[i][iLev], 0, 0,
                    nMoments, nGst);
    }
  }

  calc_mach_number();

  isMomentsUpdated = true;
}

//==========================================================
// Ma = u/vth
void Pic::calc_mach_number() {
  for (int iLev = 0; iLev < n_lev(); iLev++) {

    for (MFIter mfi(nodePlasma[nSpecies][iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const Array4<Real>& moments = nodePlasma[nSpecies][iLev][mfi].array();
      const Array4<Real>& mach = mMach[iLev][mfi].array();

      ParallelFor(box, [&](int i, int j, int k) {
        Real rho = moments(i, j, k, iRho_);
        if (rho <= 0) {
          mach(i, j, k) = 0;
          return;
        }

        Real u = moments(i, j, k, iUx_) / rho;
        Real v = moments(i, j, k, iUy_) / rho;
        Real w = moments(i, j, k, iUz_) / rho;
        Real uBulk = sqrt(u * u + v * v + w * w);

        Real p = (moments(i, j, k, iPxx_) + moments(i, j, k, iPyy_) +
                  moments(i, j, k, iPzz_)) /
                 3.0;
        Real vth = sqrt(gamma0 * p / rho);

        mach(i, j, k) = uBulk / max(vth, 1e-99);
      });
    }
  }
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

      const Array4<Real>& cost = cellCost[iLev][mfi].array();
      const Array4<int const> status = cellStatus[iLev][mfi].array();

      ParallelFor(box, [&](int i, int j, int k) {
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
      });
    }
  }
}

//==========================================================
void Pic::divE_correction() {
  std::string nameFunc = "Pic::divE_correction";

  timing_func(nameFunc);

  for (int iIter = 0; iIter < nDivECorrection; iIter++) {

    sum_to_center(true);

    if (doReport)
      Print() << "\n-----" << printPrefix << " div(E) correction at iter "
              << iIter << "----------" << std::endl;

    calculate_phi(divESolver, 0);

    divE_correct_particle_position();
  }

  for (int i = 0; i < nSpecies; ++i) {
    // The particles outside the simulation domain is marked for deletion
    // inside divE_correct_particle_position(). redistribute_particles()
    // deletes these particles. In order to get correct moments, re-inject
    // particles in the ghost cells.
    parts[i]->redistribute_particles();
  }

  inject_particles_for_boundary_cells();

  sum_to_center(false);
}

//==========================================================
void Pic::divE_correct_particle_position() {
  std::string nameFunc = "Pic::correct_position";

  timing_func(nameFunc);
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->divE_correct_position(centerPhi, iLev);
    }
  }
}

//==========================================================
void Pic::calculate_phi(LinearSolver& solver, int iLev) {
  std::string nameFunc = "Pic::calculate_phi";

  timing_func(nameFunc);

  {
    MultiFab residual(cGrids[iLev], DistributionMap(iLev), 1, nGst);

    solver.reset(get_local_node_or_cell_number(centerDivE[iLev]));
    // div_node_to_center(nodeE[iLev], residual, Geom(iLev).InvCellSize());
    MultiFab::Copy(residual, centerDivE[iLev], 0, 0, 1, nGst);
    Real coef = 1.0 / rhoTheta;

    MultiFab::LinComb(residual, coef, residual, 0, -fourPI * coef,
                      centerNetChargeN[iLev], 0, 0, residual.nComp(),
                      residual.nGrow());
    if (finest_level > 0) {
      skip_cells_divE_correction(residual, cellStatus[iLev], iLev);
    }

    convert_3d_to_1d(residual, solver.rhs, iLev);

    BL_PROFILE_VAR("Pic::phi_iterate", solve);
    solver.solve(iLev, doReport);
    BL_PROFILE_VAR_STOP(solve);

    convert_1d_to_3d(solver.xLeft, centerPhi[iLev], iLev);
    centerPhi[iLev].FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
void Pic::divE_accurate_matvec(const double* vecIn, double* vecOut, int iLev) {
  std::string nameFunc = "Pic::divE_matvec";
  timing_func(nameFunc);

  // const int iLev = 0;
  zero_array(vecOut, divESolver.get_nSolve());

  MultiFab inMF(cGrids[iLev], DistributionMap(iLev), 1, nGst);

  convert_1d_to_3d(vecIn, inMF, iLev);
  inMF.FillBoundary(0, 1, IntVect(1), Geom(iLev).periodicity());

  MultiFab outMF(cGrids[iLev], DistributionMap(iLev), 1, nGst);
  outMF.setVal(0.0);

  for (MFIter mfi(inMF); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();

    const Array4<Real>& lArr = outMF[mfi].array();
    const Array4<Real const>& rArr = inMF[mfi].array();
    const Array4<RealCMM>& mmArr = centerMM[iLev][mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      IntVect ijk = { AMREX_D_DECL(i, j, k) };
      Box subBox(ijk - 1, ijk + 1);

      ParallelFor(subBox, [&](int i2, int j2, int k2) {
        const int gp = (i2 - i + 1) * 9 + (j2 - j + 1) * 3 + k2 - k + 1;
        lArr(i, j, k) += rArr(i2, j2, k2) * mmArr(i, j, k)[gp];
      });
    });
  }
  outMF.mult(fourPI * fourPI);
  convert_3d_to_1d(outMF, vecOut, iLev);
}

//==========================================================
void Pic::sum_to_center(bool isBeforeCorrection) {
  std::string nameFunc = "Pic::sum_to_center";

  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    centerNetChargeNew[iLev].setVal(0.0);

    const RealCMM mm0(0.0);
    centerMM[iLev].setVal(mm0);

    bool doNetChargeOnly = !isBeforeCorrection;

    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->sum_to_center(centerNetChargeNew[iLev], centerMM[iLev],
                              doNetChargeOnly, iLev);
    }

    if (!doNetChargeOnly) {
      centerMM[iLev].SumBoundary(Geom(iLev).periodicity());
    }

    centerNetChargeNew[iLev].SumBoundary(Geom(iLev).periodicity());

    if (iLev == 0) {
      apply_BC(cellStatus[iLev], centerNetChargeNew[iLev], 0,
               centerNetChargeNew[iLev].nComp(), &Pic::get_zero, iLev);
    }

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
void Pic::sum_to_center_new(bool isBeforeCorrection, int iLev) {

  std::string nameFunc = "Pic::sum_to_center";
  timing_func(nameFunc);

  bool doNetChargeOnly = !isBeforeCorrection;

  centerNetChargeNew[iLev].setVal(0.0);
  const RealCMM mm0(0.0);
  centerMM[iLev].setVal(mm0);

  MultiFab jf;
  MultiFab jc;
  int fLev = iLev + 1;
  int cLev = iLev - 1;
  if (iLev == 0) {
    cLev = iLev;
  }
  if (iLev == finest_level) {
    fLev = iLev;
  }
  {
    BoxArray bac = centerB[cLev].boxArray();
    bac.refine(ref_ratio[iLev]);
    jc.define(bac, centerB[cLev].DistributionMap(), 1, 1);
    jc.setVal(0.0);
    BoxArray baf = centerB[fLev].boxArray();
    baf.coarsen(ref_ratio[iLev]);
    baf.grow(1);
    jf.define(baf, centerB[fLev].DistributionMap(), 1, 1);
    jf.setVal(0.0);
  }
  for (int i = 0; i < nSpecies; ++i) {
    parts[i]->sum_to_center_new(centerNetChargeNew[iLev], jc, jf,
                                centerMM[iLev], doNetChargeOnly, iLev);
  }

  if (!doNetChargeOnly) {
    centerMM[iLev].SumBoundary(Geom(iLev).periodicity());
  }

  centerNetChargeNew[iLev].SumBoundary(Geom(iLev).periodicity());
  jc.SumBoundary();
  jf.SumBoundary();
  amrex::MultiFab tmp;
  tmp.define(centerB[iLev].boxArray(), centerB[iLev].DistributionMap(), 1, 0);
  tmp.setVal(0.0);
  tmp.ParallelCopy(jc);
  MultiFab::Add(centerNetChargeNew[iLev], tmp, 0, 0, 1, 0);
  tmp.setVal(0.0);
  tmp.ParallelCopy(jf);
  MultiFab::Add(centerNetChargeNew[iLev], tmp, 0, 0, 1, 0);

  if (iLev == 0) {
    apply_BC(cellStatus[iLev], centerNetChargeNew[iLev], 0,
             centerNetChargeNew[iLev].nComp(), &Pic::get_zero, iLev);
  }

  MultiFab::LinComb(
      centerNetChargeN[iLev], 1 - rhoTheta, centerNetChargeOld[iLev], 0,
      rhoTheta, centerNetChargeNew[iLev], 0, 0, centerNetChargeN[iLev].nComp(),
      centerNetChargeN[iLev].nGrow());

  if (!isBeforeCorrection) {
    MultiFab::Copy(centerNetChargeOld[iLev], centerNetChargeNew[iLev], 0, 0,
                   centerNetChargeOld[iLev].nComp(),
                   centerNetChargeOld[iLev].nGrow());
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

  if (reportParticleQuality) {
    if (tc->get_cycle() % 20 == 0) {
      WriteParticleQualityToParaView();
    }
  }

  if (solveFieldInCoMov || solvePartInCoMov || useUpwindB ||
      (useUpwindE && cMaxE < 0)) {
    update_U0_E0();
  }

  if (solveEM) {
    if (finest_level == 0) {
      calc_mass_matrix();
    } else {
      calc_mass_matrix_amr();
    }
  }

  if (solveEM) {
    update_E();
  }

  particle_mover();

  // Calling re_sampling after particle mover so that all the particles
  // outside the domain have been deleted.
  re_sampling();

  charge_exchange();

  if (source) {
    fill_source_particles();
  }

  inject_particles_for_boundary_cells();

  isMomentsUpdated = false;

  if (solveEM) {
    if (projectDownEmFields) {
      project_down_E();
    }
    update_B();
  }

  // Only to be turned on if DivE error needs to be visulaized when DivE
  // cleaning is not turned on

  // for (int i = 0; i < 2; i++) {
  //   sum_moments(true);
  //   sum_moments(false);
  // }

  if (solveEM && doCorrectDivE) {
    if (finest_level == 0) {
      divE_correction();
    } else {
      amr_divE_correction();
    }
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
    uBg[iLev].setVal(0.0);
    eBg[iLev].setVal(0.0);
    for (MFIter mfi(uBg[iLev]); mfi.isValid(); ++mfi) {
      const Array4<Real>& arrU = uBg[iLev][mfi].array();
      const Array4<const Real>& arrMoments =
          nodePlasma[nSpecies][iLev][mfi].array();

      const Array4<const int>& status = nodeStatus[iLev][mfi].array();

      // Fill in the physical nodes
      ParallelFor(mfi.validbox(), [&](int i, int j, int k) {
        const Real rho = arrMoments(i, j, k, iRho_);
        if (rho > 0) {
          const Real invRho = 1. / rho;
          for (int iu = iUx_; iu <= iUz_; iu++)
            arrU(i, j, k, iu - iUx_) = arrMoments(i, j, k, iu) * invRho;
        }
      });

      // Fill in ghost nodes
      ParallelFor(mfi.fabbox(), [&](int i, int j, int k) {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };
        if (bit::is_domain_boundary(status(ijk))) {
          const int iFluid = 0;
          for (int iDir = 0; iDir < nDim3; iDir++) {
            arrU(i, j, k, iDir) =
                get_node_fluid_u(mfi, ijk, iDir, iLev, iFluid);
          }
        }
      });
    }

    uBg[iLev].FillBoundary(Geom(iLev).periodicity());

    for (int i = 0; i < nSmoothBackGroundU; ++i)
      smooth_multifab(uBg[iLev], iLev, i % 2 + 1);

    for (MFIter mfi(uBg[iLev]); mfi.isValid(); ++mfi) {
      const Array4<Real>& arrU = uBg[iLev][mfi].array();
      const Array4<Real>& arrE = eBg[iLev][mfi].array();
      const Array4<Real>& arrB = nodeB[iLev][mfi].array();

      const Array4<const int>& status = nodeStatus[iLev][mfi].array();

      // Fill in the physical nodes
      ParallelFor(mfi.validbox(), [&](int i, int j, int k) {
        const Real& bx = arrB(i, j, k, ix_);
        const Real& by = arrB(i, j, k, iy_);
        const Real& bz = arrB(i, j, k, iz_);

        const Real& ux = arrU(i, j, k, ix_);
        const Real& uy = arrU(i, j, k, iy_);
        const Real& uz = arrU(i, j, k, iz_);

        arrE(i, j, k, ix_) = -uy * bz + uz * by;
        arrE(i, j, k, iy_) = -uz * bx + ux * bz;
        arrE(i, j, k, iz_) = -ux * by + uy * bx;
      });

      // Fill in boundary nodes
      ParallelFor(mfi.fabbox(), [&](int i, int j, int k) {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };
        if (bit::is_domain_boundary(status(ijk))) {
          arrE(i, j, k, ix_) = get_node_E(mfi, ijk, ix_, iLev);
          arrE(i, j, k, iy_) = get_node_E(mfi, ijk, iy_, iLev);
          arrE(i, j, k, iz_) = get_node_E(mfi, ijk, iz_, iLev);
        }
      });
    }

    eBg[iLev].FillBoundary(Geom(iLev).periodicity());

    // for (int i = 0; i < nSmoothE; ++i)
    //   smooth_multifab(eBg[iLev], iLev, i % 2 + 1);
  }
}

//==========================================================
void Pic::update_E() {
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
             &Pic::get_center_B, iLev, &bcBField);
  }
  const Real dt = tc->get_dt();
  RealVect dt2dx;
  for (int i = 0; i < nDim; ++i) {
    dt2dx[i] = dt * Geom(0).InvCellSize(i);
  }
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    curl_center_to_node(centerB[iLev], nodeE[iLev], dt2dx.begin());
    MultiFab::Saxpy(nodeE[iLev], -fourPI * dt, jHat[iLev], 0, 0,
                    nodeE[iLev].nComp(), nodeE[iLev].nGrow());

    MultiFab::Add(nodeE[iLev], nodeEth[iLev], 0, 0, nodeE[iLev].nComp(),
                  nodeE[iLev].nGrow());

    nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_BC(nodeStatus[iLev], nodeE[iLev], 0, nDim3, &Pic::get_node_E, iLev);
  }
}

//==========================================================
void Pic::update_E_impl() {
  std::string nameFunc = "Pic::update_E_impl";

  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    eSolver.reset(get_local_node_or_cell_number(nodeE[iLev]));

    update_E_rhs(eSolver.rhs, iLev);

    convert_3d_to_1d(nodeE[iLev], eSolver.xLeft, iLev);

    update_E_matvec(eSolver.xLeft, eSolver.matvec, iLev, false);

    for (int i = 0; i < eSolver.get_nSolve(); ++i) {
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

    if (doSmoothE) {
      smooth_E(nodeEth[iLev], iLev);
    }

    MultiFab::Add(nodeEth[iLev], nodeE[iLev], 0, 0, nodeEth[iLev].nComp(),
                  nGst);

    MultiFab::LinComb(nodeE[iLev], -(1.0 - fsolver.theta) / fsolver.theta,
                      nodeE[iLev], 0, 1. / fsolver.theta, nodeEth[iLev], 0, 0,
                      nodeE[iLev].nComp(), nGst);

    if (iLev == 0) {

      apply_BC(nodeStatus[iLev], nodeE[iLev], 0, nDim3, &Pic::get_node_E, iLev);
      apply_BC(nodeStatus[iLev], nodeEth[iLev], 0, nDim3, &Pic::get_node_E,
               iLev);

    } else {

      fill_fine_lev_bny_from_coarse(
          nodeE[iLev - 1], nodeE[iLev], 0, nodeE[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          node_bilinear_interp);

      fill_fine_lev_bny_from_coarse(
          nodeEth[iLev - 1], nodeEth[iLev], 0, nodeEth[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          node_bilinear_interp);
    }

    if (doSmoothE) {
      smooth_E(nodeEth[iLev], iLev);
      smooth_E(nodeE[iLev], iLev);
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

  if (isFake2D) {
    // Make sure there is no variation in the z-direction.
    Periodicity period(IntVect(AMREX_D_DECL(0, 0, 1)));
    vecMF.FillBoundary(period);
  }

  if (useZeroBC) {
    // The boundary nodes would not be filled in by convert_1d_3d. So, there
    // is not need to apply zero boundary conditions again here.
  } else {
    // Even after apply_BC(), the outmost layer node E is still
    // unknow. See FluidInterface::calc_current for detailed explaniation.
    if (iLev == 0) {
      apply_BC(nodeStatus[iLev], vecMF, 0, nDim3, &Pic::get_node_E, iLev);
    } else {
      fill_fine_lev_bny_from_coarse(
          nodeEth[iLev - 1], vecMF, 0, nodeEth[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          node_bilinear_interp);
    }
  }

  lap_node_to_node(vecMF, matvecMF, DistributionMap(iLev), Geom(iLev));

  Real delt2 = pow(fsolver.theta * tc->get_dt(), 2);
  matvecMF.mult(-delt2);

  if (useUpwindE) {
    // Explicit scheme: add the LF artificial viscosity term to the rhs
    // vis_{i+0.5} = c_max/2*(E_i+1 - E_i)
    // E_i += dt/dx*(vis_{i+0.5} - vis_{i-0.5}) = 0.5*c_max*dt*dx*lap(E_i)
    // For implicit scheme, we add it to the lhs, so the sign changes.

    const Real dx = Geom(iLev).CellSize()[0];
    const Real coe1 = -0.5 * fsolver.theta * tc->get_dt() / dx;

    for (MFIter mfi(vecMF); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrE = vecMF[mfi].array();
      const Array4<Real>& res = matvecMF[mfi].array();
      const Array4<Real>& arrU = uBg[iLev][mfi].array();

      ParallelFor(box, vecMF.nComp(), [&](int i, int j, int k, int iVar) {
        for (int iDir = 0; iDir < nDim; iDir++) {
          Real dii[nDim3] = { 0, 0, 0 };
          dii[iDir] = 1;

          Real cR = limiter_theta(
              limiterThetaE,
              arrE(i - dii[ix_], j - dii[iy_], k - dii[iz_], iVar),
              arrE(i, j, k, iVar),
              arrE(i + dii[ix_], j + dii[iy_], k + dii[iz_], iVar));

          Real cL = limiter_theta(
              limiterThetaE,
              arrE(i - 2 * dii[ix_], j - 2 * dii[iy_], k - 2 * dii[iz_], iVar),
              arrE(i - dii[ix_], j - dii[iy_], k - dii[iz_], iVar),
              arrE(i, j, k, iVar));

          Real ur = cMaxE, ul = cMaxE;

          if (cMaxE < 0) {
            ul = fabs(0.5 *
                      (arrU(i - dii[ix_], j - dii[iy_], k - dii[iz_], iDir) +
                       arrU(i, j, k, iDir)));
            ur = fabs(0.5 *
                      (arrU(i, j, k, iDir) +
                       arrU(i + dii[ix_], j + dii[iy_], k + dii[iz_], iDir)));
          }

          Real dE = cR * ur *
                        (arrE(i + dii[ix_], j + dii[iy_], k + dii[iz_], iVar) -
                         arrE(i, j, k, iVar)) -
                    cL * ul *
                        (arrE(i, j, k, iVar) -
                         arrE(i - dii[ix_], j - dii[iy_], k - dii[iz_], iVar));

          res(i, j, k, iVar) += coe1 * dE;
        }
      });
    }
  }

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

    const Array4<Real const>& inArr = inMF[mfi].array();
    const Array4<Real>& outArr = outMF[mfi].array();
    const Array4<RealMM>& mmArr = nodeMM[iLev][mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      IntVect ijk = { AMREX_D_DECL(i, j, k) };

      auto& data0 = mmArr(ijk);

      Box subBox(ijk - 1, ijk + 1);

      ParallelFor(subBox, [&](int i2, int j2, int k2) {
        const int gp = (k2 - k + 1) * 9 + (j2 - j + 1) * 3 + i2 - i + 1;
        const int idx0 = gp * 9;

        Real* const M_I = &(data0[idx0]);

        const double& vctX = inArr(i2, j2, k2, ix_); // vectX[i2][j2][k2];
        const double& vctY = inArr(i2, j2, k2, iy_);
        const double& vctZ = inArr(i2, j2, k2, iz_);
        outArr(i, j, k, ix_) +=
            (vctX * M_I[0] + vctY * M_I[1] + vctZ * M_I[2]) * c0;
        outArr(i, j, k, iy_) +=
            (vctX * M_I[3] + vctY * M_I[4] + vctZ * M_I[5]) * c0;
        outArr(i, j, k, iz_) +=
            (vctX * M_I[6] + vctY * M_I[7] + vctZ * M_I[8]) * c0;
      });
    });
  }

  // if (doSmoothJ) {
  //   for (int icount = 0; icount < nSmoothJ; icount++) {
  //     smooth_multifab(outMF, iLev, icount % 2 + 1);
  //   }
  // }
}

//==========================================================
void Pic::update_E_rhs(double* rhs, int iLev) {
  std::string nameFunc = "Pic::update_E_rhs";
  timing_func(nameFunc);

  MultiFab tempNode(nGrids[iLev], DistributionMap(iLev), 3, nGst);
  tempNode.setVal(0.0);
  MultiFab temp2Node(nGrids[iLev], DistributionMap(iLev), 3, nGst);
  temp2Node.setVal(0.0);

  if (iLev == 0) {
    apply_BC(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
             &Pic::get_center_B, iLev, &bcBField);
    apply_BC(nodeStatus[iLev], nodeB[iLev], 0, nodeB[iLev].nComp(),
             &Pic::get_node_B, iLev, &bcBField);
  } else {
    fill_fine_lev_bny_from_coarse(
        centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
        ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), cell_status(iLev),
        cell_bilinear_interp);

    fill_fine_lev_bny_from_coarse(nodeB[iLev - 1], nodeB[iLev], 0,
                                  nodeB[iLev - 1].nComp(), ref_ratio[iLev - 1],
                                  Geom(iLev - 1), Geom(iLev), node_status(iLev),
                                  node_bilinear_interp);
  }
  const Real* invDx = Geom(iLev).InvCellSize();

  curl_center_to_node(centerB[iLev], tempNode, invDx);

  MultiFab::Saxpy(temp2Node, -fourPI, jHat[iLev], 0, 0, temp2Node.nComp(),
                  temp2Node.nGrow());

  MultiFab::Add(temp2Node, tempNode, 0, 0, tempNode.nComp(), temp2Node.nGrow());

  temp2Node.mult(fsolver.theta * tc->get_dt());
  MultiFab::Add(temp2Node, nodeE[iLev], 0, 0, nodeE[iLev].nComp(),
                temp2Node.nGrow());

  if (solveFieldInCoMov) {
    tempNode.setVal(0.0);
    update_E_M_dot_E(eBg[iLev], tempNode, iLev);
    MultiFab::Add(temp2Node, tempNode, 0, 0, tempNode.nComp(),
                  tempNode.nGrow());
  }

  convert_3d_to_1d(temp2Node, rhs, iLev);
}

//==========================================================
void Pic::update_B() {
  std::string nameFunc = "Pic::update_B";
  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    MultiFab dB(cGrids[iLev], DistributionMap(iLev), 3, nGst);
    curl_node_to_center(nodeEth[iLev], dB, Geom(iLev).InvCellSize());

    MultiFab::Saxpy(centerB[iLev], -tc->get_dt(), dB, 0, 0,
                    centerB[iLev].nComp(), centerB[iLev].nGrow());

    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
  }
  if (projectDownEmFields && finest_level > 0) {
    for (int iLev = finest_level; iLev > 0; iLev--) {
      average_down(centerB[iLev], centerB[iLev - 1], 0, 3, ref_ratio[0]);
    }
  }
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_BC(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
               &Pic::get_center_B, iLev, &bcBField);

    } else {
      fill_fine_lev_bny_from_coarse(
          centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          cell_bilinear_interp);
    }
    MultiFab::Copy(dBdt[iLev], nodeB[iLev], 0, 0, dBdt[iLev].nComp(),
                   dBdt[iLev].nGrow());

    if (useUpwindB) {
      div_node_to_center(nodeB[iLev], divB[iLev], Geom(iLev).InvCellSize());
      smooth_B(iLev);
    }
    average_center_to_node(centerB[iLev], nodeB[iLev]);
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());

    const Real invDt = 1. / tc->get_dt();
    // dBdt = (B^{n+1} - B^n)/dt;
    MultiFab::LinComb(dBdt[iLev], -invDt, dBdt[iLev], 0, invDt, nodeB[iLev], 0,
                      0, dBdt[iLev].nComp(), dBdt[iLev].nGrow());

    if (iLev == 0) {

      apply_BC(nodeStatus[iLev], nodeB[iLev], 0, nodeB[iLev].nComp(),
               &Pic::get_node_B, iLev, &bcBField);

    } else {

      fill_fine_lev_bny_from_coarse(
          nodeB[iLev - 1], nodeB[iLev], 0, nodeB[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          node_bilinear_interp);
    }
  }
}

//==========================================================
void Pic::solve_hyp_phi(int iLev) {
  std::string nameFunc = "Pic::solve_hyp_phi";
  timing_func(nameFunc);

  // divB error propagation speed
  Real ch = 0.8 * Geom(iLev).CellSize()[ix_] / tc->get_dt();

  Real coef = -tc->get_dt() * pow(ch, 2);
  for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
    Box box = mfi.validbox();

    const Array4<Real>& divBArr = divB[iLev][mfi].array();
    const Array4<Real>& phiArr = hypPhi[iLev][mfi].array();

    ParallelFor(box, [&](int i, int j, int k) {
      IntVect ijk = { AMREX_D_DECL(i, j, k) };
      phiArr(ijk) += coef * divBArr(ijk);
      phiArr(ijk) *= (1 - hypDecay);
    });
  }

  hypPhi[iLev].FillBoundary(Geom(iLev).periodicity());

  apply_BC(cellStatus[iLev], hypPhi[iLev], 0, hypPhi[iLev].nComp(), nullptr,
           iLev);
}
//==========================================================
void Pic::smooth_B(int iLev) {
  std::string nameFunc = "Pic::smooth_B";
  timing_func(nameFunc);

  MultiFab gradPhi(cGrids[iLev], DistributionMap(iLev), nDim3, 0);
  gradPhi.setVal(0.0);

  if (useHyperbolicCleaning) {
    solve_hyp_phi(iLev);

    MultiFab gradPhiNode(nGrids[iLev], DistributionMap(iLev), nDim3, 0);
    gradPhiNode.setVal(0.0);

    grad_center_to_node(hypPhi[iLev], gradPhiNode, Geom(iLev).InvCellSize());

    average_node_to_cellcenter(gradPhi, 0, gradPhiNode, 0, nDim3,
                               gradPhi.nGrow());
  }

  MultiFab centerDB(cGrids[iLev], DistributionMap(iLev), nDim3, nGst);
  centerDB.setVal(0.0);

  Real coef[nDim3];
  for (int i = 0; i < nDim3; ++i) {
    coef[i] = 0.5 * tc->get_dt() * Geom(iLev).InvCellSize()[i];
  }

  for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
    Box box = mfi.validbox();

    const Array4<Real>& cB = centerB[iLev][mfi].array();
    const Array4<Real>& nB = nodeB[iLev][mfi].array();
    const Array4<Real const>& nU = uBg[iLev][mfi].array();
    const Array4<Real>& dB = centerDB[mfi].array();
    const Array4<Real>& gradPhiArr = gradPhi[mfi].array();
    const auto& status = cellStatus[iLev][mfi].array();
    const Array4<Real const>& moments = nodePlasma[nSpecies][iLev][mfi].array();

    // Get the face along the direction iDir for the cell (i,j,k) for the iVar
    // component
    auto get_face = [&](int iDir, int i, int j, int k, int iVar,
                        Array4<Real const> const& arr, Real& l, Real& r) {
      if (testCase == TopHat) {
        l = 1;
        r = 1;
        return;
      }

      int kp1 = nDim > 2 ? k + 1 : k;
      if (iDir == ix_) {
        l = 0.25 * (arr(i, j, k, iVar) + arr(i, j + 1, k, iVar) +
                    arr(i, j, kp1, iVar) + arr(i, j + 1, kp1, iVar));
        r = 0.25 * (arr(i + 1, j, k, iVar) + arr(i + 1, j + 1, k, iVar) +
                    arr(i + 1, j, kp1, iVar) + arr(i + 1, j + 1, kp1, iVar));
      } else if (iDir == iy_) {
        l = 0.25 * (arr(i, j, k, iVar) + arr(i + 1, j, k, iVar) +
                    arr(i, j, kp1, iVar) + arr(i + 1, j, kp1, iVar));

        r = 0.25 * (arr(i, j + 1, k, iVar) + arr(i + 1, j + 1, k, iVar) +
                    arr(i, j + 1, kp1, iVar) + arr(i + 1, j + 1, kp1, iVar));

      } else if (iDir == iz_) {
        l = 0.25 * (arr(i, j, k, iVar) + arr(i, j + 1, k, iVar) +
                    arr(i + 1, j, k, iVar) + arr(i + 1, j + 1, k, iVar));

        r = 0.25 * (arr(i, j, kp1, iVar) + arr(i, j + 1, kp1, iVar) +
                    arr(i + 1, j, kp1, iVar) + arr(i + 1, j + 1, kp1, iVar));
      }
    };

    auto get_alfven = [&](int iDir, int i, int j, int k, Real& lAlfven,
                          Real& rAlfven) {
      // Left and right B
      Real lB[nDim3] = { 0, 0, 0 }, rB[nDim3] = { 0, 0, 0 }, lBt, rBt;
      Real lRho, rRho;
      for (int ivar = 0; ivar < nDim3; ivar++) {
        get_face(iDir, i, j, k, ivar, nB, lB[ivar], rB[ivar]);
      }
      lBt = l2_norm(lB, nDim3);
      rBt = l2_norm(rB, nDim3);

      get_face(iDir, i, j, k, iRho_, moments, lRho, rRho);
      lAlfven = lBt / sqrt(lRho);
      rAlfven = rBt / sqrt(rRho);
    };

    ParallelFor(box, [&](int i, int j, int k) {
      bool doDiffusion;
      Real lu[nDim3] = { 0, 0, 0 }, ru[nDim3] = { 0, 0, 0 }, lumin, rumin;
      Real ul, ur;

      IntVect ijk{ AMREX_D_DECL(i, j, k) };

      Real lAlfven = 0, rAlfven = 0;

      // Flux along  x
      for (int iDir = 0; iDir < nDim; iDir++) {
        get_face(iDir, i, j, k, iDir, nU, lu[iDir], ru[iDir]);
      }
      lumin = smoothBIso * max(fabs(lu[ix_]), fabs(lu[iy_]), fabs(lu[iz_]));
      rumin = smoothBIso * max(fabs(ru[ix_]), fabs(ru[iy_]), fabs(ru[iz_]));

      get_alfven(ix_, i, j, k, lAlfven, rAlfven);

      ul = lu[ix_];
      ur = ru[ix_];
      doDiffusion = true;
      if ((ul > 0 && bit::is_domain_boundary(status(i - 1, j, k))) ||
          (ur < 0 && bit::is_domain_boundary(status(i + 1, j, k)))) {
        doDiffusion = false;
      }

      if (smoothBIso < 0) {
        ul = fabs(ul) + lAlfven;
        ur = fabs(ur) + rAlfven;
      } else {
        ul = max(fabs(ul), lumin);
        ur = max(fabs(ur), rumin);
      }

      if (doDiffusion)
        for (int iVar = 0; iVar < nDim3; iVar++) {
          Real cR = limiter_theta(limiterTheta, cB(i - 1, j, k, iVar),
                                  cB(i, j, k, iVar), cB(i + 1, j, k, iVar));
          Real cL = limiter_theta(limiterTheta, cB(i - 2, j, k, iVar),
                                  cB(i - 1, j, k, iVar), cB(i, j, k, iVar));
          ul = min(ul, 0.5 / coef[ix_]);
          ur = min(ur, 0.5 / coef[ix_]);
          dB(i, j, k, iVar) +=
              (cR * ur * (cB(i + 1, j, k, iVar) - cB(i, j, k, iVar)) -
               cL * ul * (cB(i, j, k, iVar) - cB(i - 1, j, k, iVar))) *
              coef[ix_];
        }

      // Flux along y
      ul = lu[iy_];
      ur = ru[iy_];
      doDiffusion = true;
      if ((ul > 0 && bit::is_domain_boundary(status(i, j - 1, k))) ||
          (ur < 0 && bit::is_domain_boundary(status(i, j + 1, k)))) {
        doDiffusion = false;
      }

      get_alfven(iy_, i, j, k, lAlfven, rAlfven);

      if (smoothBIso < 0) {
        ul = fabs(ul) + lAlfven;
        ur = fabs(ur) + rAlfven;
      } else {
        ul = max(fabs(ul), lumin);
        ur = max(fabs(ur), rumin);
      }

      if (doDiffusion)
        for (int iVar = 0; iVar < nDim3; iVar++) {
          Real cR = limiter_theta(limiterTheta, cB(i, j - 1, k, iVar),
                                  cB(i, j, k, iVar), cB(i, j + 1, k, iVar));
          Real cL = limiter_theta(limiterTheta, cB(i, j - 2, k, iVar),
                                  cB(i, j - 1, k, iVar), cB(i, j, k, iVar));
          ul = min(ul, 0.5 / coef[iy_]);
          ur = min(ur, 0.5 / coef[iy_]);

          dB(i, j, k, iVar) +=
              (cR * ur * (cB(i, j + 1, k, iVar) - cB(i, j, k, iVar)) -
               cL * ul * (cB(i, j, k, iVar) - cB(i, j - 1, k, iVar))) *
              coef[iy_];
        }

      if (nDim > 2 && !isFake2D) {

        // Flux along z
        ul = lu[iz_];
        ur = ru[iz_];

        doDiffusion = true;
        if ((ul > 0 && bit::is_domain_boundary(status(i, j, k - 1))) ||
            (ur < 0 && bit::is_domain_boundary(status(i, j, k + 1)))) {
          doDiffusion = false;
        }
        get_alfven(iz_, i, j, k, lAlfven, rAlfven);

        if (smoothBIso < 0) {
          ul = fabs(ul) + lAlfven;
          ur = fabs(ur) + rAlfven;
        } else {
          ul = max(fabs(ul), lumin);
          ur = max(fabs(ur), rumin);
        }

        if (doDiffusion)
          for (int iVar = 0; iVar < nDim3; iVar++) {
            Real cR = limiter_theta(limiterTheta, cB(i, j, k - 1, iVar),
                                    cB(i, j, k, iVar), cB(i, j, k + 1, iVar));
            Real cL = limiter_theta(limiterTheta, cB(i, j, k - 2, iVar),
                                    cB(i, j, k - 1, iVar), cB(i, j, k, iVar));
            ul = min(ul, 0.5 / coef[iz_]);
            ur = min(ur, 0.5 / coef[iz_]);

            dB(i, j, k, iVar) +=
                (cR * ur * (cB(i, j, k + 1, iVar) - cB(i, j, k, iVar)) -
                 cL * ul * (cB(i, j, k, iVar) - cB(i, j, k - 1, iVar))) *
                coef[iz_];
          }
      }

      if (useHyperbolicCleaning) {
        for (int iVar = 0; iVar < nDim3; iVar++) {
          dB(ijk, iVar) += -tc->get_dt() * gradPhiArr(ijk, iVar);
        }
      }
    });
  }

  MultiFab::Add(centerB[iLev], centerDB, 0, 0, nDim3, 0);

  centerB[iLev].FillBoundary(Geom(iLev).periodicity());
}

//==========================================================
void Pic::smooth_multifab(MultiFab& mf, int iLev, int di, Real coef) {
  std::string nameFunc = "Pic::smooth_multifab";
  timing_func(nameFunc);

  MultiFab mfOld(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrow());

  auto smooth_dir = [&](int iDir) {
    int dIdx[3] = { 0, 0, 0 };
    dIdx[iDir] = di;

    MultiFab::Copy(mfOld, mf, 0, 0, mf.nComp(), mf.nGrow());

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();

      Array4<Real> const& arrE = mf[mfi].array();
      Array4<Real> const& arrTmp = mfOld[mfi].array();

      ParallelFor(box, mf.nComp(), [&](int i, int j, int k, int iVar) {
        const Real weightSelf = 1 - coef;
        const Real WeightNei = coef / 2.0;

        const Real neiSum =
            arrTmp(i - dIdx[ix_], j - dIdx[iy_], k - dIdx[iz_], iVar) +
            arrTmp(i + dIdx[ix_], j + dIdx[iy_], k + dIdx[iz_], iVar);

        arrE(i, j, k, iVar) =
            weightSelf * arrTmp(i, j, k, iVar) + WeightNei * neiSum;
      });
    }

    mf.FillBoundary(Geom(iLev).periodicity());
  };

  smooth_dir(ix_);
  if (nDim > 1)
    smooth_dir(iy_);
  if (nDim > 2 && !isFake2D)
    smooth_dir(iz_);
}

//==========================================================
void Pic::smooth_E(MultiFab& mfE, int iLev) {
  std::string nameFunc = "Pic::smooth_E";
  timing_func(nameFunc);

  for (int icount = 0; icount < nSmoothE; icount++) {
    smooth_multifab(mfE, iLev, icount % 2 + 1);
  }
}
//==========================================================
void Pic::project_down_E() {
  if (finest_level > 0) {
    for (int iLev = finest_level; iLev > 0; iLev--) {
      amrex::MultiFab tmp(nGrids[iLev], DistributionMap(iLev), 3, 0);
      tmp.setVal(0.0);
      for (MFIter mfi(tmp); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real>& arrE = nodeE[iLev][mfi].array();
        const Array4<Real>& arrTmp = tmp[mfi].array();
        ParallelFor(box, [&](int i, int j, int k) {
          for (int iVar = 0; iVar < 3; iVar++) {
            if (nDim == 3) {
              arrTmp(i, j, k, iVar) =
                  0.5 * arrE(i, j, k, iVar) +
                  (1 / 12.0) *
                      (arrE(i + 1, j, k, iVar) + arrE(i - 1, j, k, iVar) +
                       arrE(i, j + 1, k, iVar) + arrE(i, j - 1, k, iVar) +
                       arrE(i, j, k + 1, iVar) + arrE(i, j, k - 1, iVar));
            } else {
              arrTmp(i, j, k, iVar) =
                  0.5 * arrE(i, j, k, iVar) +
                  (1 / 8.0) *
                      (arrE(i + 1, j, k, iVar) + arrE(i - 1, j, k, iVar) +
                       arrE(i, j + 1, k, iVar) + arrE(i, j - 1, k, iVar));
            }
          }
        });
      }
      fill_fine_lev_edge_from_coarse(
          nodeE[iLev - 1], tmp, 0, nodeE[iLev].nComp(), ref_ratio[iLev],
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
      average_down_nodal(tmp, nodeE[iLev - 1], ref_ratio[iLev - 1]);
    }
    for (int iLev = 0; iLev <= finest_level; iLev++) {
      nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }
}
//==========================================================
void Pic::apply_BC(const iMultiFab& status, MultiFab& mf, const int iStart,
                   const int nComp, GETVALUE func, const int iLev,
                   const BC* bc) {
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
  if (nDim > 2 &&
      Geom(iLev).Domain().bigEnd(iz_) == Geom(iLev).Domain().smallEnd(iz_)) {
    ba.grow(iz_, ngrow[iz_]);
  }

  if (bc != nullptr) {
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bxFab = mfi.fabbox();
      Box bxValid = mfi.validbox();

      //! if there are cells not in the valid + periodic grown box
      //! we need to fill them here
      if (!ba.contains(bxFab)) {
        Array4<Real> const& arr = mf[mfi].array();

        const Array4<const int>& statusArr = status[mfi].array();

        ParallelFor(bxFab, [&](int i, int j, int k) {
          if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {

            int ip, jp, kp;
            bool useFloat = use_float(i, j, k, ip, jp, kp, *bc, bxValid);

            if (useFloat) {
              for (int iVar = iStart; iVar < iStart + nComp; iVar++) {
                arr(i, j, k, iVar) = arr(ip, jp, kp, iVar);
              }
            } else {
              for (int iVar = iStart; iVar < iStart + nComp; iVar++) {
                arr(i, j, k, iVar) = (this->*func)(
                    mfi, IntVect{ AMREX_D_DECL(i, j, k) }, iVar - iStart, iLev);
              }
            }
          }
        });
      }
    }

    return;
  }

  if (useFloatBC) {
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& bxFab = mfi.fabbox();
      const Box& bxValid = mfi.validbox();

      //! if there are cells not in the valid + periodic grown box
      //! we need to fill them here
      if (!ba.contains(bxFab)) {
        Array4<Real> const& arr = mf[mfi].array();

        const Array4<const int>& statusArr = status[mfi].array();

        Box box = bxValid;
        box.grow(1);

        ParallelFor(box, [&](int i, int j, int k) {
          if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
            bool isNeiFound = false;

            // Find the neighboring physical cell
            Box subBox(IntVect(-1), IntVect(1));
            ParallelFor(subBox, [&](int ii, int jj, int kk) {
              if (!isNeiFound &&
                  !bit::is_lev_boundary(statusArr(i + ii, j + jj, k + kk, 0))) {
                isNeiFound = true;
                for (int iVar = iStart; iVar < iStart + nComp; iVar++) {
                  arr(i, j, k, iVar) = arr(i + ii, j + jj, k + kk, iVar);
                }
              }
            });
          }
        });
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

        auto lo = IntVect(bx.loVect());
        auto hi = IntVect(bx.hiVect());
        if (nDim > 2) {
          lo[iz_]++;
          hi[iz_]--;
        }

        Box box0(lo, hi);

        ParallelFor(box0, nComp, [&](int i, int j, int k, int iVar) {
          if (bit::is_lev_boundary(statusArr(i, j, k, 0))) {
            arr(i, j, k, iStart + iVar) = (this->*func)(
                mfi, IntVect{ AMREX_D_DECL(i, j, k) }, iVar, iLev);
          }
        });
      }
    }
  }
}

//==========================================================
Real Pic::calc_E_field_energy() {
  Real sum = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = nodeE[iLev][mfi];
      const auto& status = cell_status(iLev)[mfi].array();
      Box box = mfi.validbox();
      const Array4<Real>& arr = fab.array();

      Real sumLoc = 0;
      ParallelFor(box, [&](int i, int j, int k) {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };

        if (!bit::is_refined(status(ijk))) {
          Box subBox(ijk, ijk + 1);
          ParallelFor(subBox, [&](int ii, int jj, int kk) {
            IntVect ijk0 = { AMREX_D_DECL(ii, jj, kk) };
            sumLoc += pow(arr(ijk0, ix_), 2) + pow(arr(ijk0, iy_), 2) +
                      pow(arr(ijk0, iz_), 2);
          });
        }
      });

      Real avg = (nDim == 3) ? 0.125 : 0.25;

      sum += sumLoc * 0.5 * avg * get_cell_volume(iLev) / fourPI;
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
      const auto& status = cell_status(iLev)[mfi].array();

      const Box& box = mfi.validbox();
      const Array4<Real>& arr = fab.array();

      Real sumLoc = 0;
      ParallelFor(box, [&](int i, int j, int k) {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };
        if (!bit::is_refined(status(ijk))) {
          sumLoc += pow(arr(i, j, k, ix_), 2) + pow(arr(i, j, k, iy_), 2) +
                    pow(arr(i, j, k, iz_), 2);
        }
      });

      sum += sumLoc * get_cell_volume(iLev) * 0.5 / fourPI;
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

    const Array4<Real>& arr = MF[mfi].array();

    const auto& nodeArr = nodeStatus[iLev][mfi].array();

    ParallelFor(box, MF.nComp(), [&](int i, int j, int k, int iVar) {
      if (isCenter || bit::is_owner(nodeArr(i, j, k))) {
        arr(i, j, k, iVar) = p[iCount++];
      }
    });
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

    const Array4<Real const>& arr = MF[mfi].array();

    const auto& nodeArr = nodeStatus[iLev][mfi].array();

    ParallelFor(box, MF.nComp(), [&](int i, int j, int k, int iVar) {
      if (isCenter || bit::is_owner(nodeArr(i, j, k))) {
        p[iCount++] = arr(i, j, k, iVar);
      }
    });
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
  for (int i = 0; i < nProc; ++i) {
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
      printf("|     Value          |      Min      |     Avg      |      Max "
             "    "
             "|where(max)|\n");

      Vector<std::string> varType = {
        "|Blocks # of",
        "|Cells  # of",
        "|Parts  # of",
        "|Memory(MB)          |",
      };

      for (int iLev = 0; iLev <= n_lev(); iLev++) {
        for (int i = iNBlk_; i <= iNParts_; ++i) {
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

  for (int i = 0; i < nSpecies; ++i) {
    parts[i]->charge_exchange(tc->get_dt(), stateOH, sourcePT2OH, source,
                              kineticSource, sourceParts, doSelectRegion,
                              product(nSourcePPC));
  }

  if (kineticSource) {
    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->add_source_particles(sourceParts[i], nSourcePPC,
                                     adaptiveSourcePPC);
      sourceParts[i]->clearParticles();
    }

  } else {
    // 'source' is applied to generate new particles every step, so
    // sum_boundary() is called here to correct boundary nodes. Boundary nodes
    // of 'sourcePT2OH' should be corrected just before PT->OH coupling,
    // instead of here.
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

void Pic::fill_lightwaves(amrex::Real wavelength, int EorB, amrex::Real time,
                          int lev) {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    nodeE[iLev].setVal(0.0);
    nodeB[iLev].setVal(0.0);
    centerB[iLev].setVal(0.0);
    if (lev != -1 && iLev != lev)
      continue;
    for (MFIter mfi(nodeE[iLev]); mfi.isValid(); ++mfi) {
      FArrayBox& fab = nodeE[iLev][mfi];
      FArrayBox& fab2 = nodeB[iLev][mfi];

      const Box& box = mfi.fabbox();
      const Array4<Real>& arrE = fab.array();
      const Array4<Real>& arrB = fab2.array();
      const auto& prob_lo = geom[iLev].ProbLo();
      const auto& dx = geom[iLev].CellSize();
      ParallelFor(box, [&](int i, int j, int k) {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };
        if (EorB == -1 || EorB == 0) {
          arrE(ijk, ix_) =
              -0.8 * sin((2.0 * (dPI) *
                          ((prob_lo[0] + dx[0] * i) * 0.6 +
                           (prob_lo[1] + dx[1] * j) * 0.8 - time)) /
                         wavelength);
          arrE(ijk, iy_) = 0.6 * sin((2.0 * (dPI) *
                                      ((prob_lo[0] + dx[0] * i) * 0.6 +
                                       (prob_lo[1] + dx[1] * j) * 0.8 - time)) /
                                     wavelength);
          arrE(ijk, iz_) = -cos((2.0 * (dPI) *
                                 ((prob_lo[0] + dx[0] * i) * 0.6 +
                                  (prob_lo[1] + dx[1] * j) * 0.8 - time)) /
                                wavelength);
        }
        if (EorB == -1 || EorB == 1) {
          arrB(ijk, ix_) =
              -0.8 * cos((2.0 * (dPI) *
                          ((prob_lo[0] + dx[0] * i) * 0.6 +
                           (prob_lo[1] + dx[1] * j) * 0.8 - time)) /
                         wavelength);
          arrB(ijk, iy_) = 0.6 * cos((2.0 * (dPI) *
                                      ((prob_lo[0] + dx[0] * i) * 0.6 +
                                       (prob_lo[1] + dx[1] * j) * 0.8 - time)) /
                                     wavelength);
          arrB(ijk, iz_) = sin((2.0 * (dPI) *
                                ((prob_lo[0] + dx[0] * i) * 0.6 +
                                 (prob_lo[1] + dx[1] * j) * 0.8 - time)) /
                               wavelength);
        }
      });
    }

    nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {

      FArrayBox& fab = centerB[iLev][mfi];

      const Box& box = mfi.fabbox();
      const Array4<Real>& arrcB = fab.array();
      const auto& prob_lo = geom[iLev].ProbLo();
      const auto& dx = geom[iLev].CellSize();
      ParallelFor(box, [&](int i, int j, int k) {
        IntVect ijk = { AMREX_D_DECL(i, j, k) };
        if (EorB == -1 || EorB == 1) {
          arrcB(ijk, ix_) =
              -0.8 * cos((2.0 * (dPI) *
                          (((prob_lo[0] + dx[0] * (i + 0.5)) * 0.6 +
                            (prob_lo[1] + dx[1] * (j + 0.5)) * 0.8 - time))) /
                         wavelength);
          arrcB(ijk, iy_) =
              0.6 * cos((2.0 * (dPI) *
                         (((prob_lo[0] + dx[0] * (i + 0.5)) * 0.6 +
                           (prob_lo[1] + dx[1] * (j + 0.5)) * 0.8 - time))) /
                        wavelength);
          arrcB(ijk, iz_) =
              sin((2.0 * (dPI) *
                   (((prob_lo[0] + dx[0] * (i + 0.5)) * 0.6 +
                     (prob_lo[1] + dx[1] * (j + 0.5)) * 0.8 - time))) /
                  wavelength);
        }
      });
    }

    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
  }
}

void Pic::amr_divE_correction() {
  std::string nameFunc = "Pic::divE_correction";

  timing_func(nameFunc);

  for (int iIter = 0; iIter < nDivECorrection; iIter++) {
    for (int iLev = finest_level; iLev >= 0; iLev--) {
      sum_to_center_new(true, iLev);
      skip_cells_divE_correction(centerMM[iLev], cell_status(iLev), iLev);
      calculate_phi(divESolver, iLev);
      for (int i = 0; i < nSpecies; ++i) {
        parts[i]->divE_correct_position(centerPhi, iLev);
      }
      if (finest_level > 0) {
        for (int i = 0; i < nSpecies; ++i) {
          parts[i]->Redistribute();
        }
      }
    }
  }

  inject_particles_for_boundary_cells();
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    sum_to_center_new(false, iLev);
    if (iLev > 0) {
      skip_cells_divE_correction(centerNetChargeN[iLev], cell_status(iLev),
                                 iLev);
      skip_cells_divE_correction(centerDivE[iLev], cell_status(iLev), iLev);
    }
  }
}
