#include <algorithm>
#include <cctype>
#include <limits>
#include <math.h>
#include <vector>

#include <AMReX_Algorithm.H>
#include <AMReX_CArena.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_MultiFabUtil.H>

#if defined(__linux__)
#include <malloc.h>
#endif

#include "GridUtility.h"
#include "LinearSolver.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

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

  if (pInfo.isPPVconstant || pInfo.doPreSplitting) {
    SetTargetPPC(2);
    isTargetPPCDefined = true;
    for (int i = 0; i < nSpecies; i++) {
      parts[i]->set_is_target_ppc_defined(isTargetPPCDefined);
    }
  }
  if (initEM) {
    fill_E_B_fields();
  }
  init_regional_fields();

  // Every registered InitialCondition plug-in seeds its fields through the
  // narrow PicICFields facade (LightWave, HybridWave, ConvectionWave, ...). The
  // hybrid-wave velocity kick and all per-particle modifications are applied
  // inside fill_particles() via the plugin.
  if (ic_) {
    PicICFields icf = ic_fields();
    ic_->set_fields(icf);
  }

  if (usePIC) {
    // Macroparticle seeding (and any per-particle modifications such as the
    // beam bulk override or the hybrid-wave Alfven velocity kick) is routed
    // through the InitialCondition plugin during fill_particles().
    fill_particles();
    sum_moments(true);
    if (useHybridPIC) {
      seed_first_hybrid_step();
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        average_center_to_node(centerB[iLev], nodeB[iLev]);
        nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
        if (iLev == 0) {
          apply_field_bc(nodeStatus[iLev], nodeB[iLev], 0, 3, &Pic::get_node_B,
                         iLev, true);
        }
      }
      if (useElectronPressureEq && !peStateInitialized_) {
        init_electron_pressure();
      }
      compute_ambipolar_E();
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        assemble_ohm_E(centerB[iLev], total_center_B(centerB[iLev], iLev),
                       nodeE[iLev], iLev, 0.0);
        nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
        if (iLev == 0) {
          apply_field_bc(nodeStatus[iLev], nodeE[iLev], 0, nDim3,
                         &Pic::get_node_E, iLev, false);
        }
      }
    }
    // div(E)-correction fields are full-PIC only.
    if (!useHybridPIC) {
      if (finest_level == 0) {
        sum_to_center(false);
      } else if (doCorrectDivE) {
        for (int iLev = 0; iLev < n_lev(); iLev++) {
          sum_to_center_amr(false, iLev);
        }
      }
    }

    if (need_divB()) {
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        compute_divB(iLev);
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

  // The relaxation workspace is needed only on coarse levels with an active
  // fine level. The active-level loop below cannot release a workspace on a
  // level that disappeared during coarsening, so clear those slots explicitly.
  if (useHybridPIC && !relaxTarget.empty()) {
    for (int iLev = finest_level; iLev < n_lev_max(); ++iLev) {
      relaxTarget[iLev].clear();
      relaxResidual[iLev].clear();
      relaxCurl[iLev].clear();
      relaxCurlT[iLev].clear();
      relaxDir[iLev].clear();
    }
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
    if (use_intrinsic_B()) {
      // B0 is analytic: there is nothing to copy on a regrid, the array is
      // simply rebuilt and refilled by fill_intrinsic_B().
      distribute_FabArray(nodeB0[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                          nGst, false);
      distribute_FabArray(centerB0[iLev], cGrids[iLev], DistributionMap(iLev),
                          3, nGst, false);
      // The Ohm's-law scratch is rebuilt on demand, and must not survive a
      // regrid with a stale BoxArray.
      centerBtotal[iLev].clear();
    }
    distribute_FabArray(nodeE[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst);
    // Keep old theta-field values where grids overlap and initialize nodes
    // introduced by AMR before the next field solve.
    distribute_FabArray(nodeEth[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst, true, 0.0);

    bool doMoveData = false;
    // div(E)/div(B) correction and implicit E-solver arrays (full-PIC only).
    if (!useHybridPIC) {
      distribute_FabArray(centerNetChargeOld[iLev], cGrids[iLev],
                          DistributionMap(iLev), 1, nGst);
      distribute_FabArray(centerNetChargeN[iLev], cGrids[iLev],
                          DistributionMap(iLev), 1, nGst);
      distribute_FabArray(centerNetChargeNew[iLev], cGrids[iLev],
                          DistributionMap(iLev), 1, nGst);
      distribute_FabArray(centerDivE[iLev], cGrids[iLev], DistributionMap(iLev),
                          1, nGst);
      distribute_FabArray(centerPhi[iLev], cGrids[iLev], DistributionMap(iLev),
                          1, nGst);

      // A diagnostic: read as zero, not as uninitialised, when not cleaned.
      distribute_FabArray(divB[iLev], cGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData, 0.0);
      distribute_FabArray(hypPhi[iLev], cGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);

      distribute_FabArray(centerDB[iLev], cGrids[iLev], DistributionMap(iLev),
                          nDim3, nGst, doMoveData);

      if (projectDownEmFields) {
        distribute_FabArray(projectScratchMF[iLev], nGrids[iLev],
                            DistributionMap(iLev), 3, 0, doMoveData);
      }

      if (!useExplicitPIC && !useHybridPIC) {
        distribute_FabArray(nodeMM[iLev], nGrids[iLev], DistributionMap(iLev),
                            1, 1, doMoveData);
        if (nodeMM_comm_data.size() != n_lev()) {
          nodeMM_comm_data.resize(n_lev());
        }
        nodeMM_comm_data[iLev].is_initialized = false;
        distribute_FabArray(solverVecMF[iLev], nGrids[iLev],
                            DistributionMap(iLev), 3, nGst, doMoveData);
        distribute_FabArray(solverMatvecMF[iLev], nGrids[iLev],
                            DistributionMap(iLev), 3, 1, doMoveData);
        distribute_FabArray(solverTempNode3[iLev], nGrids[iLev],
                            DistributionMap(iLev), 3, nGst, doMoveData);
        distribute_FabArray(solverCenterLapMF[iLev], cGrids[iLev],
                            DistributionMap(iLev), 3, 1, doMoveData);
        if (fsolver.coefDiff > 0) {
          distribute_FabArray(solverTempCenter3[iLev], cGrids[iLev],
                              DistributionMap(iLev), 3, nGst, doMoveData);
          distribute_FabArray(solverTempCenter1[iLev], cGrids[iLev],
                              DistributionMap(iLev), 1, nGst, doMoveData);
        }
        distribute_FabArray(solverRhsNode1[iLev], nGrids[iLev],
                            DistributionMap(iLev), 3, nGst, doMoveData);
        distribute_FabArray(solverRhsNode2[iLev], nGrids[iLev],
                            DistributionMap(iLev), 3, nGst, doMoveData);
      }
      if (useExplicitPIC) {
        distribute_FabArray(nodeCurlB[iLev], nGrids[iLev],
                            DistributionMap(iLev), 3, nGst, doMoveData);
        distribute_FabArray(nodeJ[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                            nGst, doMoveData);
        distribute_FabArray(centerBstar[iLev], cGrids[iLev],
                            DistributionMap(iLev), 3, nGst, doMoveData);
        distribute_FabArray(centerBstage[iLev], cGrids[iLev],
                            DistributionMap(iLev), 3, nGst, doMoveData);
      }
    }
    if (useHybridPIC || useExplicitPIC) {
      // Previous-step moments for the Ohm's-law / explicit interpolation.
      if (nodePlasmaPrev.empty()) {
        nodePlasmaPrev.resize(nSpecies + 1);
      }
    }
    if (useHybridPIC) {
      // Hyper-resistivity scratch: centerLapB = Laplacian(B); nodeHyperE
      // node-centered.
      distribute_FabArray(centerLapB[iLev], cGrids[iLev], DistributionMap(iLev),
                          3, nGst, doMoveData);
      distribute_FabArray(nodeHyperE[iLev], nGrids[iLev], DistributionMap(iLev),
                          3, nGst, doMoveData);

      // RK4 / ssprk3 shared intermediate solver scratch.
      distribute_FabArray(centerBstage[iLev], cGrids[iLev],
                          DistributionMap(iLev), 3, nGst, doMoveData);
      for (int kk = 0; kk < 4; ++kk)
        distribute_FabArray(kStage[iLev][kk], cGrids[iLev],
                            DistributionMap(iLev), 3, nGst, doMoveData);

      // rk3/rk4 persistent scratch: centerBstart = B_n; centerBstar =
      // (trial+B_n)/2.
      distribute_FabArray(centerBstart[iLev], cGrids[iLev],
                          DistributionMap(iLev), 3, nGst, doMoveData);
      distribute_FabArray(centerBstar[iLev], cGrids[iLev],
                          DistributionMap(iLev), 3, nGst, doMoveData);

      // Workspace for the divergence-free relaxation of the covered coarse B.
      // Allocate it only for coarse levels with an active fine level. Using
      // n_lev_max() here would allocate the AMR-only buffers in single-level
      // runs that merely allow AMR, which needlessly regresses memory use.
      if (finest_level > 0 && iLev < finest_level) {
        distribute_FabArray(relaxTarget[iLev], cGrids[iLev],
                            DistributionMap(iLev), nDim3, 0, doMoveData, 0.0);
        distribute_FabArray(relaxResidual[iLev], cGrids[iLev],
                            DistributionMap(iLev), nDim3, 1, doMoveData, 0.0);
        distribute_FabArray(relaxCurl[iLev], cGrids[iLev],
                            DistributionMap(iLev), nDim3, 0, doMoveData, 0.0);
        distribute_FabArray(relaxCurlT[iLev], nGrids[iLev],
                            DistributionMap(iLev), nDim3, 0, doMoveData, 0.0);
        distribute_FabArray(relaxDir[iLev], nGrids[iLev], DistributionMap(iLev),
                            nDim3, 0, doMoveData, 0.0);
      }

      // Staggered hybrid solver fields.
      distribute_FabArray(nodeEstage[iLev], nGrids[iLev], DistributionMap(iLev),
                          3, nGst, doMoveData);
      distribute_FabArray(nodeJ[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);
      distribute_FabArray(nodeBstage[iLev], nGrids[iLev], DistributionMap(iLev),
                          3, nGst, doMoveData);
      distribute_FabArray(centerPe[iLev], cGrids[iLev], DistributionMap(iLev),
                          1, nGst, doMoveData);
      distribute_FabArray(nodeEambi[iLev], nGrids[iLev], DistributionMap(iLev),
                          3, nGst, doMoveData);
      distribute_FabArray(nodeRhoTemp[iLev], nGrids[iLev],
                          DistributionMap(iLev), 1, nGst, doMoveData);
      if (useElectronPressureEq) {
        // The evolved field must survive a regrid (doCopy = true); the
        // scratch arrays are rebuilt and recomputed every step.
        distribute_FabArray(centerPeState[iLev], cGrids[iLev],
                            DistributionMap(iLev), 1, nGst, true);
        distribute_FabArray(centerPeRho[iLev], cGrids[iLev],
                            DistributionMap(iLev), 1, nGst, doMoveData);
        distribute_FabArray(centerPeTe[iLev], cGrids[iLev],
                            DistributionMap(iLev), 1, nGst, doMoveData);
        distribute_FabArray(nodePeVec[iLev], nGrids[iLev],
                            DistributionMap(iLev), 3, nGst, doMoveData);
        distribute_FabArray(nodePeRho[iLev], nGrids[iLev],
                            DistributionMap(iLev), 1, nGst, doMoveData);
        distribute_FabArray(nodePeAux[iLev], nGrids[iLev],
                            DistributionMap(iLev), 4, nGst, doMoveData);
      }
      if (hasRegionalResistivity_) {
        distribute_FabArray(nodeEtaRegional[iLev], nGrids[iLev],
                            DistributionMap(iLev), 1, nGst, false);
      }
      if (hasRegionalHyper_) {
        distribute_FabArray(nodeEtaHyperRegional[iLev], nGrids[iLev],
                            DistributionMap(iLev), 1, nGst, false);
      }
    }
    if (useHybridPIC || useExplicitPIC) {
      // Node-grid previous-step moments (J^{n-1/2}), slim layout.
      for (auto& pl : nodePlasmaPrev) {
        if (pl.empty())
          pl.resize(n_lev_max());
        distribute_FabArray(pl[iLev], nGrids[iLev], DistributionMap(iLev),
                            nHybridMomentsComps, nGst, doMoveData);
      }
    }
    distribute_FabArray(dBdt[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst, doMoveData);

    // mMach: node grid for both full-PIC and hybrid.
    distribute_FabArray(mMach[iLev], nGrids[iLev], DistributionMap(iLev), 1,
                        nGst, doMoveData);

    // Co-moving frame fields (eBg/uBg), div(E) mass matrix (centerMM), implicit
    // E current (jHat), and node-centered moments (nodePlasma): full-PIC only.
    if (!useHybridPIC) {
      distribute_FabArray(eBg[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);

      distribute_FabArray(uBg[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);

      distribute_FabArray(centerMM[iLev], cGrids[iLev], DistributionMap(iLev),
                          1, nGst, doMoveData);

      distribute_FabArray(divEInMF[iLev], cGrids[iLev], DistributionMap(iLev),
                          1, 1, false);

      distribute_FabArray(divEOutMF[iLev], cGrids[iLev], DistributionMap(iLev),
                          1, 0, false);

      distribute_FabArray(jHat[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);
    }

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

  // B0 lives on the same grids as B1 and is analytic, so the new boxes are
  // filled by re-evaluating the model rather than by interpolating the old
  // grid. This is the one place a regrid could silently leave B0 at zero.
  fill_intrinsic_B();
  init_regional_fields();

  // Boxes created by the regrid have no electron-pressure history.
  fill_new_electron_pressure();

  {
    iTot = nSpecies;
    if (plasmaEnergy.empty()) {
      plasmaEnergy.resize(nSpecies + 1);
    }
  }
  //===========Move field data around end====================

  //--------------particles-----------------------------------
  if (parts.empty()) {
    // Let the plugin apply any particle-count override (e.g. LightWave /
    // TopHat force zero macroparticles) after #PARTICLES has been fully parsed
    // so it always wins.
    if (ic_)
      ic_->apply_particle_override(pInfo);

    for (int i = 0; i < nSpecies; ++i) {
      auto ptr = std::make_unique<PicParticles>(
          this, fi, tc, i, fi->get_species_charge(i), fi->get_species_mass(i),
          pInfo, pMode, ic_.get());

      parts.push_back(std::move(ptr));

      auto ptrSource = std::make_unique<PicParticles>(
          this, fi, tc, i, fi->get_species_charge(i), fi->get_species_mass(i),
          pInfo, pMode, ic_.get());

      sourceParts.push_back(std::move(ptrSource));
    }

    if (waveBC.active) {
      for (auto& p : parts) {
        p->waveVelocityKick = [this](const Real* pos, Real t, Real& dvx,
                                     Real& dvy, Real& dvz) {
          wave_velocity_kick(pos, t, dvx, dvy, dvz);
        };
      }
    }
  } else {
    for (int i = 0; i < nSpecies; ++i) {
      // Label the particles outside the NEW PIC region.
      parts[i]->label_particles_outside_active_region_general();

      parts[i]->redistribute_particles();
      if (i < sourceParts.size() && sourceParts[i]) {
        // sourceParts is emptied after charge exchange. Redistribute also
        // resizes its AMReX level storage; otherwise get_particle_tile() can
        // access a missing level when refinement adds a finer grid.
        sourceParts[i]->redistribute_particles();
      }
    }
  }

  // Propagate field wave boundaries to each species to drive wave velocity
  // kicks.
  for (int d = 0; d < nDim; ++d) {
    const bool isWaveLo = (bcField.face(d, 0) == FieldBC::wave);
    const bool isWaveHi = (bcField.face(d, 1) == FieldBC::wave);
    for (auto& p : parts) {
      p->set_wave_face(d, 0, isWaveLo);
      p->set_wave_face(d, 1, isWaveHi);
    }
  }
  //--------------particles-----------------------------------

  // This part does not really work for multi-level.
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    int n = get_local_node_or_cell_number(nodeE[iLev]);
    eSolver.init(n, nDim3, nDim, matvec_E_solver);

    // divESolver uses the full-PIC-only centerDivE array.
    if (!useHybridPIC) {
      n = get_local_node_or_cell_number(centerDivE[iLev]);
      divESolver.init(n, 1, nDim, matvec_divE_accurate);
    }
  }
}

//==========================================================
void Pic::fill_new_node_E() {
  {
    Real xL = 0, xR = 0;
    if (ic_ && ic_->is_tophat()) {
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
          if (ic_ && ic_->is_tophat()) {
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

  // The initial/new electric field has to satisfy the inner-body condition
  // (see #BODYBOUNDARY) as well, otherwise the body would start with the
  // ambient field.
  if (useBody) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      apply_body_E_bc(nodeE[iLev], iLev);
    }
  }
}

//==========================================================
void Pic::fill_new_node_B() {
  {
    Real xL = 0, xR = 0;
    if (ic_ && ic_->is_tophat()) {
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
          if (ic_ && ic_->is_tophat()) {
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
    auto& cellInterp = *get_cell_interp();
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_new_from_coarse(centerB[iLev - 1], centerB[iLev], 0,
                                    centerB[iLev - 1].nComp(),
                                    ref_ratio[iLev - 1], Geom(iLev - 1),
                                    Geom(iLev), cell_status(iLev), cellInterp);
    }
  }
}

//==========================================================
void Pic::fill_E_B_fields() {
  fill_new_node_E();
  fill_new_node_B();
  fill_new_center_B();

  fill_intrinsic_B();

  //-----Coarse (iLev=0) grid boundary/internal ghost cells are filled----

  nodeE[0].FillBoundary(Geom(0).periodicity());
  nodeB[0].FillBoundary(Geom(0).periodicity());
  centerB[0].FillBoundary(Geom(0).periodicity());
  // NOTE: apply_field_bc() also applies the wave hard source.
  apply_field_bc(nodeStatus[0], nodeB[0], 0, nDim3, &Pic::get_node_B, 0, true);
  apply_field_bc(nodeStatus[0], nodeE[0], 0, nDim3, &Pic::get_node_E, 0, false);
  apply_field_bc(cellStatus[0], centerB[0], 0, centerB[0].nComp(),
                 &Pic::get_center_B, 0, true);

  //-----Fine (iLev>0) grid boundary/internal ghost cells are filled----
  auto& cellInterp = *get_cell_interp();
  for (int iLev = 1; iLev <= finest_level; iLev++) {
    nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());

    fill_fine_lev_bny_from_coarse(nodeE[iLev - 1], nodeE[iLev], 0,
                                  nodeE[iLev - 1].nComp(), ref_ratio[iLev - 1],
                                  Geom(iLev - 1), Geom(iLev), node_status(iLev),
                                  node_bilinear_interp);

    fill_fine_lev_bny_from_coarse(nodeB[iLev - 1], nodeB[iLev], 0,
                                  nodeB[iLev - 1].nComp(), ref_ratio[iLev - 1],
                                  Geom(iLev - 1), Geom(iLev), node_status(iLev),
                                  node_bilinear_interp);

    fill_fine_lev_bny_from_coarse(centerB[iLev - 1], centerB[iLev], 0,
                                  centerB[iLev - 1].nComp(),
                                  ref_ratio[iLev - 1], Geom(iLev - 1),
                                  Geom(iLev), cell_status(iLev), cellInterp);
  }

  // Seed the theta field from E on all levels after newly refined E nodes have
  // been filled, so the first field solve on the new mesh has valid data.
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::Copy(nodeEth[iLev], nodeE[iLev], 0, 0, nodeE[iLev].nComp(),
                   nodeE[iLev].nGrow());
  }

  // In hybrid PIC, sync nodeB from centerB at initialization.
  if (useHybridPIC) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      average_center_to_node(centerB[iLev], nodeB[iLev]);
      nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
      if (iLev == 0) {
        apply_field_bc(nodeStatus[iLev], nodeB[iLev], 0, 3, &Pic::get_node_B,
                       iLev, true);
      }
    }
  }

  if (is_body_conducting()) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      project_body_B(centerB[iLev], iLev);
      project_body_B(nodeB[iLev], iLev);
      // A perfect conductor shields its cavity: the initial uniform B is not
      // tangential and leaks into the surface nodes through the nodal average.
      mask_body_interior(centerB[iLev], cellStatus[iLev]);
      mask_body_interior(nodeB[iLev], nodeStatus[iLev]);
      centerB[iLev].FillBoundary(Geom(iLev).periodicity());
      nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }
}

//==========================================================
// Fill the frozen intrinsic field B0 from the analytic model, at
// initialization and after every regrid. Dipole-only evaluations run in
// parallel (ParallelFor), while spherical harmonic evaluations use host
// scratch.
//==========================================================
void Pic::fill_intrinsic_B() {
  if (!use_intrinsic_B())
    return;

  const IntrinsicBField& field = *intrinsicB_;
  const auto p = field.get_params();
  const int activeDim = get_dim();
  const bool onlyDipole = !field.use_crustal();

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (nodeB0[iLev].empty())
      distribute_FabArray(nodeB0[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                          nGst, false);
    if (centerB0[iLev].empty())
      distribute_FabArray(centerB0[iLev], cGrids[iLev], DistributionMap(iLev),
                          3, nGst, false);

    const auto plo = Geom(iLev).ProbLo();
    const auto dx = Geom(iLev).CellSize();

    if (onlyDipole) {
      // Dipole evaluation is scratch-free and device-callable: run in
      // ParallelFor.
      for (MFIter mfi(nodeB0[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.fabbox();
        const Array4<Real>& arr = nodeB0[iLev][mfi].array();
        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          Real b[3] = { 0.0, 0.0, 0.0 };
          const Real x = plo[ix_] + i * dx[ix_];
          const Real y = (activeDim > 1) ? (plo[iy_] + j * dx[iy_]) : 0.0;
#if (AMREX_SPACEDIM > 2)
          const Real z =
              (activeDim > 2) ? (plo[iz_] + k * dx[iz_]) : p.center[2];
#else
          const Real z = p.center[2];
#endif
          eval_dipole_b(p, x, y, z, b);
          arr(i, j, k, ix_) = b[0];
          arr(i, j, k, iy_) = b[1];
          arr(i, j, k, iz_) = b[2];
        });
      }

      for (MFIter mfi(centerB0[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.fabbox();
        const Array4<Real>& arr = centerB0[iLev][mfi].array();
        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          Real b[3] = { 0.0, 0.0, 0.0 };
          const Real x = plo[ix_] + (i + 0.5) * dx[ix_];
          const Real y =
              (activeDim > 1) ? (plo[iy_] + (j + 0.5) * dx[iy_]) : 0.0;
#if (AMREX_SPACEDIM > 2)
          const Real z =
              (activeDim > 2) ? (plo[iz_] + (k + 0.5) * dx[iz_]) : p.center[2];
#else
          const Real z = p.center[2];
#endif
          eval_dipole_b(p, x, y, z, b);
          arr(i, j, k, ix_) = b[0];
          arr(i, j, k, iy_) = b[1];
          arr(i, j, k, iz_) = b[2];
        });
      }
    } else {
      // Crustal evaluation uses host scratch workspace.
      for (MFIter mfi(nodeB0[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.fabbox();
        const Array4<Real>& arr = nodeB0[iLev][mfi].array();
        // Under a 2D build index 2 of a Box is out of range, so the z loop
        // collapses to a single pass and z comes from the model center.
#if (AMREX_SPACEDIM > 2)
        const int kLo = box.smallEnd(2);
        const int kHi = box.bigEnd(2);
#else
        const int kLo = 0;
        const int kHi = 0;
#endif
        for (int k = kLo; k <= kHi; ++k) {
#if (AMREX_SPACEDIM > 2)
          const Real z =
              (activeDim > 2) ? (plo[iz_] + k * dx[iz_]) : p.center[2];
#else
          const Real z = p.center[2];
#endif
          for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
            const Real y = (activeDim > 1) ? (plo[iy_] + j * dx[iy_]) : 0.0;
            for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
              Real b[3] = { 0.0, 0.0, 0.0 };
              field.eval(plo[ix_] + i * dx[ix_], y, z, b);
              arr(i, j, k, ix_) = b[0];
              arr(i, j, k, iy_) = b[1];
              arr(i, j, k, iz_) = b[2];
            }
          }
        }
      }

      for (MFIter mfi(centerB0[iLev]); mfi.isValid(); ++mfi) {
        const Box& box = mfi.fabbox();
        const Array4<Real>& arr = centerB0[iLev][mfi].array();
#if (AMREX_SPACEDIM > 2)
        const int kLo = box.smallEnd(2);
        const int kHi = box.bigEnd(2);
#else
        const int kLo = 0;
        const int kHi = 0;
#endif
        for (int k = kLo; k <= kHi; ++k) {
#if (AMREX_SPACEDIM > 2)
          const Real z =
              (activeDim > 2) ? (plo[iz_] + (k + 0.5) * dx[iz_]) : p.center[2];
#else
          const Real z = p.center[2];
#endif
          for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
            const Real y =
                (activeDim > 1) ? (plo[iy_] + (j + 0.5) * dx[iy_]) : 0.0;
            for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
              Real b[3] = { 0.0, 0.0, 0.0 };
              field.eval(plo[ix_] + (i + 0.5) * dx[ix_], y, z, b);
              arr(i, j, k, ix_) = b[0];
              arr(i, j, k, iy_) = b[1];
              arr(i, j, k, iz_) = b[2];
            }
          }
        }
      }
    }
  }
}

//==========================================================
amrex::MultiFab& Pic::total_center_B(amrex::MultiFab& src, const int iLev) {
  if (!use_intrinsic_B() || centerB0[iLev].empty())
    return src;

  amrex::MultiFab& dst = centerBtotal[iLev];
  const int nG = std::max(src.nGrow(), nGst);
  if (dst.empty() || dst.nGrow() < src.nGrow())
    distribute_FabArray(dst, cGrids[iLev], DistributionMap(iLev), nDim3, nG,
                        false);

  const int nGrow = std::min(dst.nGrow(), centerB0[iLev].nGrow());
  amrex::MultiFab::LinComb(dst, 1.0, src, 0, 1.0, centerB0[iLev], 0, 0, nDim3,
                           nGrow);
  return dst;
}

//==========================================================
void Pic::add_intrinsic_B(amrex::MultiFab& dst, const int iLev) {
  if (!use_intrinsic_B())
    return;

  const amrex::MultiFab& src =
      dst.ixType().cellCentered() ? centerB0[iLev] : nodeB0[iLev];
  if (src.empty())
    return;

  const int nComp = std::min(dst.nComp(), src.nComp());
  const int nGrow = std::min(dst.nGrow(), src.nGrow());
  amrex::MultiFab::Add(dst, src, 0, 0, nComp, nGrow);
}

//==========================================================
void Pic::fill_particles() {
  inject_particles_for_new_cells();
  inject_particles_for_boundary_cells();
  inject_body_ionosphere();
}

void Pic::fill_source_particles() {
  if (kineticSource)
    return;

  bool doSelectRegion = false;
#ifdef _PT_COMPONENT_
  doSelectRegion = (nSpecies == 4);
#endif

  if (source) {
    for (int i : kineticSpecies_) {
      parts[i]->add_particles_source(source, stateOH, tc->get_dt(), nSourcePPC,
                                     doSelectRegion, adaptiveSourcePPC,
                                     sourceMinWeightRatio);
    }
  }
}

//==========================================================
void Pic::inject_body_ionosphere() {
  const auto& bi = bodyIono;
  if (!bi.use || !useBody)
    return;

  const Real s0 = bi.solarDir[0], s1 = bi.solarDir[1], s2 = bi.solarDir[2];
  const Real floor = bi.cosSzaFloor;
  auto sza_factor = [=](const Real* xyz, Real p) {
    const Real r =
        std::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
    const Real cosSza =
        (r > 0) ? (xyz[0] * s0 + xyz[1] * s1 + xyz[2] * s2) / r : 1.0;
    return std::pow(std::max(cosSza, floor), p);
  };

  for (int i : kineticSpecies_) {
    const bool isElectron = parts[i]->get_charge() < 0;
    std::function<Real(const Real*)> nDensAt;
    if (isElectron) {
      // Quasi-neutral: n_e = sum(Z_i n_i) of the injected ions.
      nDensAt = [&, sza_factor](const Real* xyz) {
        Real n = 0;
        for (int iS = 0; iS < nSpecies; iS++)
          if (fi->get_species_charge(iS) > 0 && bi.n0[iS] > 0)
            n += fi->get_species_charge(iS) * bi.n0[iS] *
                 sza_factor(xyz, bi.szaExponent[iS]);
        return n;
      };
    } else {
      if (bi.n0[i] <= 0)
        continue;
      const Real n0 = bi.n0[i], p = bi.szaExponent[i];
      nDensAt = [=](const Real* xyz) { return n0 * sza_factor(xyz, p); };
    }
    parts[i]->inject_particles_at_body(nDensAt, bi.vth[i]);
  }
}

//==========================================================
bool Pic::shell_weights_up_to_date(const ShellWeights& sw) {
  bool upToDate = (int)sw.weight.size() == n_lev();
  for (int iLev = 0; upToDate && iLev < n_lev(); iLev++)
    upToDate = sw.grids[iLev] == cell_status(iLev).boxArray() &&
               sw.dms[iLev] == cell_status(iLev).DistributionMap();
  return upToDate;
}

//==========================================================
void Pic::build_shell_weights(
    ShellWeights& sw, Real radiusSI,
    const std::function<Real(const Real*, Real)>& profile,
    const std::string& name) {
  timing_func("Pic::build_shell_weights");

  sw.weight.clear();
  sw.weight.resize(n_lev());
  sw.grids.resize(n_lev());
  sw.dms.resize(n_lev());

  const Real R = radiusSI * fi->get_Si2NoL();
  const Real* c = bodyCenter;
  // Area of the sphere in a cell: the volume of a thin shell |r - R| < h/2
  // inside the cell, divided by h, sampled on a stratified nSub^3 grid.
  constexpr int nSub = 24;
  constexpr Real hRel = 0.2;

  Real total = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const auto& cs = cell_status(iLev);
    sw.grids[iLev] = cs.boxArray();
    sw.dms[iLev] = cs.DistributionMap();
    sw.weight[iLev].define(cs.boxArray(), cs.DistributionMap(), 1, 0);
    sw.weight[iLev].setVal(0.0);

    const auto& dxl = Geom(iLev).CellSize();
    const auto& plol = Geom(iLev).ProbLo();
    const Real dxMax = std::max({ AMREX_D_DECL(dxl[0], dxl[1], dxl[2]) });
    const Real h = hRel * dxMax;
    const Real dV = AMREX_D_TERM(dxl[0], *dxl[1], *dxl[2]) /
                    (nSub * nSub * (nDim > 2 ? nSub : 1));

    for (MFIter mfi(sw.weight[iLev]); mfi.isValid(); ++mfi) {
      const auto& status = cs[mfi].const_array();
      const auto& w = sw.weight[iLev][mfi].array();
      const Box& box = mfi.validbox();
      const auto lo = lbound(box);
      const auto hi = ubound(box);
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            const int st = status(i, j, k);
            if (bit::is_body(st) || bit::is_refined(st))
              continue;
            const int ijk[3] = { i, j, k };
            Real xc[3] = { 0, 0, 0 };
            for (int iDim = 0; iDim < nDim; iDim++)
              xc[iDim] = (ijk[iDim] + 0.5) * dxl[iDim] + plol[iDim] - c[iDim];
            const Real rc =
                std::sqrt(xc[0] * xc[0] + xc[1] * xc[1] + xc[2] * xc[2]);
            if (std::abs(rc - R) > 0.9 * dxMax + h)
              continue;

            int nHit = 0;
            const int nk = nDim > 2 ? nSub : 1;
            for (int kk = 0; kk < nk; kk++)
              for (int jj = 0; jj < nSub; jj++)
                for (int ii = 0; ii < nSub; ii++) {
                  const int sub[3] = { ii, jj, kk };
                  Real x[3] = { 0, 0, 0 };
                  for (int iDim = 0; iDim < nDim; iDim++)
                    x[iDim] = xc[iDim] +
                              ((sub[iDim] + 0.5) / nSub - 0.5) * dxl[iDim];
                  const Real r =
                      std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
                  if (std::abs(r - R) < 0.5 * h)
                    nHit++;
                }
            if (nHit == 0)
              continue;
            const Real area = nHit * dV / h;

            w(i, j, k) = area * profile(xc, rc);
            total += w(i, j, k);
          }
    }
  }

  ParallelDescriptor::ReduceRealSum(total);
  if (total <= 0) {
    Print() << "Warning: " << name
            << ": the sphere crosses no plasma cell; nothing is injected.\n";
    sw.weight.clear();
    return;
  }
  for (int iLev = 0; iLev < n_lev(); iLev++)
    sw.weight[iLev].mult(1.0 / total);
  Print() << "  " << name
          << " weights rebuilt, sphere area in plasma cells / 4 pi R^2 "
             "(profile-weighted) = "
          << total / (4 * M_PI * R * R) << "\n";
}

//==========================================================
void Pic::build_emission_weights() {
  auto& ie = ionoEmission;
  if (shell_weights_up_to_date(ie.shell))
    return;

  const Real* sd = ie.solarDir;
  // Emission profile a(direction from the center).
  const bool photochem = ie.profile == "photochemistry";
  Real phiSubsolar = 0;
  if (photochem) {
    if (source == nullptr)
      Abort("Error: #IONOSPHEREEMISSIONPROFILE photochemistry needs the "
            "source terms (#SOURCE).");
    phiSubsolar = source->column_photoproduction_si(sd, 0.0, ie.radiusSI);
    if (phiSubsolar <= 0)
      Abort("Error: #IONOSPHEREEMISSIONPROFILE photochemistry: no "
            "photoionization below the emission sphere at the subsolar "
            "point.");
  }
  auto profile_factor = [&](const Real* x, Real rc) {
    const Real cosSza =
        rc > 0 ? (x[0] * sd[0] + x[1] * sd[1] + x[2] * sd[2]) / rc : 1.0;
    if (!photochem)
      // rhybrid ionoCosSzaDayConstantNight profile.
      return cosSza > 0 ? ie.noon + (ie.night - ie.noon) * (1 - cosSza)
                        : ie.night;
    const Real dir[3] = { x[0] / rc, x[1] / rc, x[2] / rc };
    const Real phi =
        source->column_photoproduction_si(dir, 0.0, ie.radiusSI) / phiSubsolar;
    return ie.night +
           (ie.noon - ie.night) * std::pow(std::max(phi, 0.0), ie.exponent);
  };
  build_shell_weights(ie.shell, ie.radiusSI, profile_factor,
                      "#IONOSPHEREEMISSION");
}

//==========================================================
void Pic::inject_ionosphere_emission() {
  const auto& ie = ionoEmission;
  if (!ie.use || !usePIC || isGridEmpty)
    return;
  timing_func("Pic::inject_ionosphere_emission");

  build_emission_weights();
  if ((int)ie.shell.weight.size() != n_lev())
    return;

  // rate [#/s] * dt [s] * weight / V [m^3] -> n [/m^3] -> code units.
  const Real dtSI = tc->get_dt_si();
  const Real l3 = std::pow(fi->get_No2SiL(), 3);
  const Real nFactor = cProtonMassSI * fi->get_Si2NoRho();

  // Charge of all the ions injected in one step [#/s].
  Real ionCharge = 0;
  for (int iS = 0; iS < nSpecies; iS++)
    if (fi->get_species_charge(iS) > 0)
      ionCharge += fi->get_species_charge(iS) * ie.rateSI[iS];

  for (int i : kineticSpecies_) {
    const bool isElectron = parts[i]->get_charge() < 0;
    const Real rate = isElectron ? ionCharge : ie.rateSI[i];
    if (rate <= 0)
      continue;

    Vector<MultiFab> nAdd(n_lev());
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      const Real* dxl = Geom(iLev).CellSize();
      const Real volSI = AMREX_D_TERM(dxl[0], *dxl[1], *dxl[2]) * l3;
      const auto& w = ie.shell.weight[iLev];
      nAdd[iLev].define(w.boxArray(), w.DistributionMap(), 1, 0);
      MultiFab::Copy(nAdd[iLev], w, 0, 0, 1, 0);
      nAdd[iLev].mult(rate * dtSI / volSI * nFactor);
    }
    parts[i]->inject_emission_shell(nAdd, ie.vth[i], ie.nPPC, bodyCenter);
  }
}

//==========================================================
void Pic::update_floating_body() {
  auto& fb = floatingBody;
  if (!fb.use || !usePIC || isGridEmpty)
    return;
  timing_func("Pic::update_floating_body");

  // Global cumulative tallies: absorbed ion charge (> 0), absorbed electron
  // charge (< 0), electrons that reached the sphere, reflected electrons.
  Vector<Real> t(5, 0.0);
  for (int i : kineticSpecies_) {
    const Real q = parts[i]->get_body_absorb_charge();
    if (parts[i]->get_charge() > 0)
      t[0] += q;
    else
      t[1] += q;
    t[2] += parts[i]->get_body_hit_weight();
    t[3] += parts[i]->get_body_hit_energy();
    t[4] += parts[i]->get_body_reflect_count();
  }
  ParallelDescriptor::ReduceRealSum(t.data(), t.size());
  const Real dQi = t[0] - fb.qIon, dQe = -(t[1] - fb.qEle);
  const Real dHitW = t[2] - fb.hitW, dHitE = t[3] - fb.hitE;
  fb.qIon = t[0];
  fb.qEle = t[1];
  fb.hitW = t[2];
  fb.hitE = t[3];
  fb.nReflect = t[4];

  if (fb.type == "sheath") {
    // Integral control toward zero net absorbed current, in units of the
    // mean normal energy per charge of the electrons reaching the sphere
    // (~kTe/e), with the relaxation time tau.
    if (dQi + dQe > 0 && dHitW > 0) {
      const Real f = (dQe - dQi) / (dQe + dQi);
      fb.phi = std::max(0.0, fb.phi + tc->get_dt_si() / fb.tauSI *
                                          (dHitE / dHitW) * f);
    }
    for (int i : kineticSpecies_)
      if (parts[i]->get_charge() < 0)
        parts[i]->set_body_sheath(true, fb.phi);
    return;
  }

  // reemit: return the net electron charge absorbed in this step.
  const Real deficit = dQe - dQi;
  if (deficit <= 0)
    return;
  if (!shell_weights_up_to_date(fb.shell))
    build_shell_weights(
        fb.shell, fb.radiusSI, [](const Real*, Real) { return 1.0; },
        "#FLOATINGBODY reemit");
  if ((int)fb.shell.weight.size() != n_lev())
    return;
  for (int i : kineticSpecies_) {
    if (parts[i]->get_charge() >= 0)
      continue;
    Vector<MultiFab> nAdd(n_lev());
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      const auto& w = fb.shell.weight[iLev];
      nAdd[iLev].define(w.boxArray(), w.DistributionMap(), 1, 0);
      MultiFab::Copy(nAdd[iLev], w, 0, 0, 1, 0);
      // Charge deficit * weight / volume -> number density of electrons.
      nAdd[iLev].mult(deficit / get_cell_volume(iLev) /
                      std::abs(parts[i]->get_charge()));
    }
    parts[i]->inject_emission_shell(nAdd, fb.vth, fb.nPPC, bodyCenter);
    fb.qReemit += deficit;
    break; // The first electron species carries the current.
  }
}

//==========================================================
void Pic::fill_plasma_state(FluidInterface& dst) {
  timing_func("Pic::fill_plasma_state");

  if (isGridEmpty || !usePIC)
    return;

  const int nS = dst.get_nS();
  const int iE = useHybridPIC ? iElectron_ : -1;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    MultiFab& stateMF = dst.get_node_fluid(iLev);
    AMREX_ALWAYS_ASSERT(stateMF.boxArray() == nodePlasma[0][iLev].boxArray() &&
                        stateMF.DistributionMap() ==
                            nodePlasma[0][iLev].DistributionMap());

    // Hybrid fluid electron: n_e = sum(Z_i n_i) from the ion moments.
    MultiFab nodeNe, nodePe;
    if (iE >= 0) {
      nodeNe.define(nodePlasma[nSpecies][iLev].boxArray(),
                    nodePlasma[nSpecies][iLev].DistributionMap(), 1, 0);
      nodePe.define(nodeNe.boxArray(), nodeNe.DistributionMap(), 1, 0);
      nodeNe.setVal(0.0);
      for (int s : kineticSpecies_) {
        const Real zom = fi->get_QoQi(s) / fi->get_MoMi(s);
        MultiFab::Saxpy(nodeNe, zom, nodePlasma[s][iLev], iRho_, 0, 1, 0);
      }
      compute_node_electron_pressure(iLev, nodeNe, nodePe);
    }

    for (MFIter mfi(stateMF); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& state = stateMF[mfi].array();

      for (int s = 0; s < nS && s < nSpecies; s++) {
        const int iRho = dst.get_iRho(s), iUx = dst.get_iUx(s),
                  iUy = dst.get_iUy(s), iUz = dst.get_iUz(s),
                  iP = dst.get_iP(s);

        if (s == iE) {
          const Real massE = fi->get_MoMi(s);
          const Array4<Real const>& ne = nodeNe[mfi].const_array();
          const Array4<Real const>& pe = nodePe[mfi].const_array();
          // Charge-weighted ion velocity, the electron velocity without the
          // Hall current. Accumulated below from the per-species momenta.
          ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            state(i, j, k, iRho) = massE * ne(i, j, k);
            state(i, j, k, iUx) = 0.0;
            state(i, j, k, iUy) = 0.0;
            state(i, j, k, iUz) = 0.0;
            state(i, j, k, iP) = (ne(i, j, k) > 0.0) ? pe(i, j, k) : 0.0;
          });
          for (int si : kineticSpecies_) {
            const Real zom = fi->get_QoQi(si) / fi->get_MoMi(si);
            const Array4<Real const>& mom = nodePlasma[si][iLev][mfi].array();
            ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
              const Real n = ne(i, j, k);
              if (n > 0.0) {
                state(i, j, k, iUx) += zom * mom(i, j, k, iMx_) / n;
                state(i, j, k, iUy) += zom * mom(i, j, k, iMy_) / n;
                state(i, j, k, iUz) += zom * mom(i, j, k, iMz_) / n;
              }
            });
          }
          continue;
        }

        const Array4<Real const>& mom = nodePlasma[s][iLev][mfi].array();
        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const Real rho = mom(i, j, k, iRho_);
          if (rho > 0.0) {
            const Real invRho = 1.0 / rho;
            state(i, j, k, iRho) = rho;
            state(i, j, k, iUx) = mom(i, j, k, iMx_) * invRho;
            state(i, j, k, iUy) = mom(i, j, k, iMy_) * invRho;
            state(i, j, k, iUz) = mom(i, j, k, iMz_) * invRho;
            // nodePlasma holds the thermal pressure tensor (see
            // convert_to_fluid_moments); the scalar pressure is its trace/3.
            // With few particles at a node, the second moment minus rho*u^2
            // can come out negative: clamp it.
            state(i, j, k, iP) =
                amrex::max(0.0, (mom(i, j, k, iPxx_) + mom(i, j, k, iPyy_) +
                                 mom(i, j, k, iPzz_)) /
                                    3.0);
          } else {
            state(i, j, k, iRho) = 0.0;
            state(i, j, k, iUx) = 0.0;
            state(i, j, k, iUy) = 0.0;
            state(i, j, k, iUz) = 0.0;
            state(i, j, k, iP) = 0.0;
          }
        });
      }
    }

    stateMF.FillBoundary(Geom(iLev).periodicity());

    // The body nodes are not plasma: with #BODYIONOSPHERE they carry the
    // ghost layer, whose ions are re-injected every step and must not react.
    // Otherwise their moments are diluted by the empty body cells (the other
    // surface nodes were rescaled in sum_moments()): rescale them as well, so
    // that the state is the density of the plasma next to the surface.
    if (useBody) {
      if (bodyIono.use)
        mask_body(stateMF, node_status(iLev));
      else {
        // Only the densities and the pressures: the velocities are not
        // diluted by the empty body cells.
        for (int s = 0; s < dst.get_nS() && s < nSpecies; s++) {
          for (const int iVar : { dst.get_iRho(s), dst.get_iP(s) }) {
            MultiFab alias(stateMF, amrex::make_alias, iVar, 1);
            scale_by_body_outside_fraction(alias, cellStatus[iLev],
                                           nodeStatus[iLev],
                                           Geom(iLev).Domain(), nDim, true);
          }
        }
      }
    }
  }
}

//==========================================================
void Pic::mask_body_source(SourceInterface& src) {
  if (!useBody || isGridEmpty)
    return;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (bodyIono.use) {
      // The body nodes carry the re-injected ghost layer: no reaction there.
      mask_body(src.get_node_fluid(iLev), node_status(iLev));
      if (src.has_loss_array(iLev))
        mask_body(src.get_node_loss_fluid(iLev), node_status(iLev));
    } else {
      // The production of a node is the average over its dual cell outside
      // the body (UserSource), but only the plasma cells around the node
      // receive particles: scale it up by the inverse outside fraction. The
      // loss frequencies are intensive and stay as they are.
      scale_by_body_outside_fraction(src.get_node_fluid(iLev), cellStatus[iLev],
                                     nodeStatus[iLev], Geom(iLev).Domain(),
                                     nDim, false);
    }
  }
}

//==========================================================
void Pic::update_part_loc_to_half_stage() {
  std::string nameFunc = "Pic::update_part_loc_to_half_stage";

  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (int i = 0; i < nSpecies; ++i) {
      if (useHybridPIC && parts[i]->get_charge() < 0)
        continue;
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
      if (maxWeightRatio > 1) {
        // With #CONSTANTPPV (or pre-splitting) the reference weight is that
        // of the level target. The cell average would not do: the light
        // source particles of a cell next to the body lower it, and the
        // ordinary particles are then split again and again.
        parts[i]->limit_weight(maxWeightRatio, parts[i]->is_neutral(),
                               pInfo.doPreSplitting || pInfo.isPPVconstant);
      }
      parts[i]->split(reSamplingLowLimit, parts[i]->is_neutral(),
                      pInfo.doPreSplitting);
      // With #CONSTANTPPV the merging goal is the target of every level.
      parts[i]->merge(reSamplingHighLimit,
                      pInfo.doPreSplitting || pInfo.isPPVconstant);
    }
  }
}

//==========================================================
void Pic::particle_mover(int pushStage) {
  std::string nameFunc = "Pic::mover";

  timing_func(nameFunc);

  Real dt = tc->get_dt();
  Real dtnext = tc->get_next_dt();

  const Vector<MultiFab>& nodeEpush =
      (useHybridPIC || useExplicitPIC) ? nodeE : nodeEth;

  for (int i : kineticSpecies_) {
    parts[i]->mover(nodeEpush, nodeB, nodeB0, eBg, uBg, dt, dtnext, pushStage);
  }

  for (int i : kineticSpecies_) {
    parts[i]->redistribute_particles();
  }
}

namespace {

template <typename Array4Type>
inline void pack_node_mm_box(const Array4Type& src, const amrex::Box& box,
                             RealMM* __restrict__ buf) {
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);
  const int len_x = hi.x - lo.x + 1;
  const std::size_t row_bytes =
      static_cast<std::size_t>(len_x) * sizeof(RealMM);
  std::size_t p = 0;
  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      std::memcpy(&buf[p], &src(lo.x, j, k), row_bytes);
      p += len_x;
    }
  }
}

inline void accumulate_node_mm_box(const RealMM* __restrict__ buf,
                                   const amrex::Box& box,
                                   const amrex::Array4<RealMM>& dst) {
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);
  const int len_x = hi.x - lo.x + 1;
  const std::size_t n_doubles = static_cast<std::size_t>(len_x) * nMMComponents;
  std::size_t p = 0;
  for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
      amrex::Real* __restrict__ dst_ptr = dst(lo.x, j, k).data;
      const amrex::Real* __restrict__ src_ptr = buf[p].data;
      for (std::size_t d = 0; d < n_doubles; ++d) {
        dst_ptr[d] += src_ptr[d];
      }
      p += len_x;
    }
  }
}

} // anonymous namespace

//==========================================================
void Pic::init_boundary_node_mm_comm(int iLev) {
  if (nodeMM_comm_data.size() != n_lev()) {
    nodeMM_comm_data.resize(n_lev());
  }

  auto& cd = nodeMM_comm_data[iLev];
  if (cd.is_initialized && cd.bdkey == nodeMM[iLev].getBDKey()) {
    return;
  }

  cd.loc_tags.clear();
  cd.local_buf.clear();
  cd.box_to_loc_tags.clear();
  cd.sends.clear();
  cd.recvs.clear();

  const auto& TheFB =
      nodeMM[iLev].getFB(amrex::IntVect(0), Geom(iLev).periodicity(), false,
                         false, false, nodeMM[iLev].nGrowVect());

  if (TheFB.m_LocTags) {
    const auto& loc_tags = *TheFB.m_LocTags;
    const int n_loc = static_cast<int>(loc_tags.size());
    cd.loc_tags.resize(n_loc);

    std::size_t total_loc_pts = 0;
    for (int i = 0; i < n_loc; ++i) {
      const auto& tag = loc_tags[i];
      cd.loc_tags[i] = { tag.srcIndex, tag.dstIndex, tag.sbox, tag.dbox,
                         total_loc_pts };
      total_loc_pts += tag.sbox.numPts();
    }
    cd.local_buf.resize(total_loc_pts);

    const int nLocalBoxes = nodeMM[iLev].size();
    cd.box_to_loc_tags.resize(nLocalBoxes);
    for (int i = 0; i < n_loc; ++i) {
      int localDst = nodeMM[iLev].localindex(cd.loc_tags[i].dstIndex);
      if (localDst >= 0 && localDst < nLocalBoxes) {
        cd.box_to_loc_tags[localDst].push_back(i);
      }
    }
  }

#ifdef BL_USE_MPI
  if (TheFB.m_SndTags) {
    cd.sends.reserve(TheFB.m_SndTags->size());
    for (const auto& kv : *TheFB.m_SndTags) {
      int dst_rank = kv.first;
      const auto& tags = kv.second;
      if (tags.empty())
        continue;
      NodeMMCommData::PeerComm peer;
      peer.rank = dst_rank;
      peer.tags.reserve(tags.size());
      std::size_t offset = 0;
      for (const auto& tag : tags) {
        peer.tags.push_back({ tag.srcIndex, tag.sbox, offset });
        offset += tag.sbox.numPts();
      }
      peer.totalPts = offset;
      peer.buf.resize(offset);
      cd.sends.push_back(std::move(peer));
    }
  }

  if (TheFB.m_RcvTags) {
    cd.recvs.reserve(TheFB.m_RcvTags->size());
    for (const auto& kv : *TheFB.m_RcvTags) {
      int src_rank = kv.first;
      const auto& tags = kv.second;
      if (tags.empty())
        continue;
      NodeMMCommData::PeerComm peer;
      peer.rank = src_rank;
      peer.tags.reserve(tags.size());
      std::size_t offset = 0;
      for (const auto& tag : tags) {
        peer.tags.push_back({ tag.dstIndex, tag.dbox, offset });
        offset += tag.dbox.numPts();
      }
      peer.totalPts = offset;
      peer.buf.resize(offset);
      cd.recvs.push_back(std::move(peer));
    }
  }

  cd.recv_reqs.assign(cd.recvs.size(), MPI_REQUEST_NULL);
  cd.send_reqs.assign(cd.sends.size(), MPI_REQUEST_NULL);
#endif

  cd.bdkey = nodeMM[iLev].getBDKey();
  cd.is_initialized = true;
}

//==========================================================
void Pic::sum_boundary_node_mm(int iLev) {
  BL_PROFILE("Pic::nodeMM_SumBoundary");

  init_boundary_node_mm_comm(iLev);
  auto& cd = nodeMM_comm_data[iLev];

#ifdef BL_USE_MPI
  const int seq_num = amrex::ParallelDescriptor::SeqNum();
  const MPI_Comm comm = amrex::ParallelDescriptor::Communicator();
  const int n_recvs = static_cast<int>(cd.recvs.size());
  const int n_sends = static_cast<int>(cd.sends.size());

  cd.recv_reqs.assign(n_recvs, MPI_REQUEST_NULL);
  cd.send_reqs.assign(n_sends, MPI_REQUEST_NULL);

  // 1. Post non-blocking receives
  for (int r = 0; r < n_recvs; ++r) {
    auto& peer = cd.recvs[r];
    if (peer.totalPts > 0) {
      MPI_Irecv(reinterpret_cast<void*>(peer.buf.data()),
                peer.totalPts * sizeof(RealMM), MPI_BYTE, peer.rank, seq_num,
                comm, &cd.recv_reqs[r]);
    }
  }

  // 2. Pack and post non-blocking sends
  for (int s = 0; s < n_sends; ++s) {
    auto& peer = cd.sends[s];
    for (const auto& rt : peer.tags) {
      const auto src_arr = nodeMM[iLev].array(rt.boxIndex);
      pack_node_mm_box(src_arr, rt.box, &peer.buf[rt.bufOffset]);
    }
    if (peer.totalPts > 0) {
      MPI_Isend(reinterpret_cast<const void*>(peer.buf.data()),
                peer.totalPts * sizeof(RealMM), MPI_BYTE, peer.rank, seq_num,
                comm, &cd.send_reqs[s]);
    }
  }
#endif

  // 3. Snapshot local tags into local_buf (Phase 1)
  const int n_loc = static_cast<int>(cd.loc_tags.size());
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n_loc; ++i) {
    const auto& tag = cd.loc_tags[i];
    const auto src_arr = nodeMM[iLev].array(tag.srcIndex);
    pack_node_mm_box(src_arr, tag.sbox, &cd.local_buf[tag.bufOffset]);
  }

  // 4. Accumulate local tags into destination boxes (Phase 2)
  const int nLocalBoxes = nodeMM[iLev].size();
#ifdef AMREX_USE_OMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int localDst = 0; localDst < nLocalBoxes; ++localDst) {
    const auto& tag_indices = cd.box_to_loc_tags[localDst];
    if (tag_indices.empty())
      continue;

    auto dst_arr = nodeMM[iLev].atLocalIdx(localDst).array();
    for (int tag_idx : tag_indices) {
      const auto& tag = cd.loc_tags[tag_idx];
      accumulate_node_mm_box(&cd.local_buf[tag.bufOffset], tag.dbox, dst_arr);
    }
  }

#ifdef BL_USE_MPI
  // 5. Wait for MPI receives and unpack/accumulate
  if (n_recvs > 0) {
    MPI_Waitall(n_recvs, cd.recv_reqs.data(), MPI_STATUSES_IGNORE);
    for (int r = 0; r < n_recvs; ++r) {
      const auto& peer = cd.recvs[r];
      for (const auto& rt : peer.tags) {
        auto dst_arr = nodeMM[iLev].array(rt.boxIndex);
        accumulate_node_mm_box(&peer.buf[rt.bufOffset], rt.box, dst_arr);
      }
    }
  }

  // 6. Wait for MPI sends to finish
  if (n_sends > 0) {
    MPI_Waitall(n_sends, cd.send_reqs.data(), MPI_STATUSES_IGNORE);
  }
#endif
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

    const MultiFab* nodeB0Lev =
        (use_intrinsic_B() && !nodeB0[iLev].empty()) ? &nodeB0[iLev] : nullptr;

    for (int i = 0; i < nSpecies; ++i) {
      if (useExplicitPIC) {
        parts[i]->calc_jhat(jHat[iLev], nodeB[iLev], nodeB0Lev, tc->get_dt());
      } else {
        parts[i]->calc_mass_matrix(nodeMM[iLev], jHat[iLev], nodeB[iLev],
                                   nodeB0Lev, uBg[iLev], tc->get_dt(), iLev,
                                   solveFieldInCoMov);
      }
    }

    if (nSpecies > 0) {
      parts[0]->apply_jhat_mirror(jHat[iLev], iLev);
    }

    Real invVol = 1;
    for (int i = 0; i < nDim; ++i) {
      invVol *= Geom(iLev).InvCellSize(i);
    }

    jHat[iLev].mult(invVol, 0, jHat[iLev].nComp(), jHat[iLev].nGrow());
    jHat[iLev].SumBoundary(Geom(iLev).periodicity());
    jHat[iLev].FillBoundary(Geom(iLev).periodicity());

    // Undo the geometric dilution of the current at the staircase surface of
    // the inner body: the cells inside the body are empty, so a CIC node on
    // the surface only collects part of the current of the plasma there.
    if (useBody && !bodyIono.use)
      rescale_body_surface_nodes(jHat[iLev], cellStatus[iLev], nodeStatus[iLev],
                                 Geom(iLev).Domain(), nDim);

    if (doSmoothJ) {
      for (int icount = 0; icount < nSmoothJ; icount++) {
        smooth_multifab(jHat[iLev], iLev, icount % 2 + 1, coefSmoothJ);
      }
    }

    if (!useExplicitPIC) {
      sum_boundary_node_mm(iLev);

      // Same correction as jHat above: the mass matrix of the plasma next to
      // the body surface is diluted by the empty body cells.
      if (useBody && !bodyIono.use)
        rescale_body_surface_nodes(nodeMM[iLev], cellStatus[iLev],
                                   nodeStatus[iLev], Geom(iLev).Domain(), nDim);
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
  amrex::Vector<amrex::Vector<NodeMMFab> > nmmc;
  amrex::Vector<NodeMMFab> nmmf;
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
      const MultiFab* nodeB0Lev = (use_intrinsic_B() && !nodeB0[iLev].empty())
                                      ? &nodeB0[iLev]
                                      : nullptr;
      parts[i]->calc_mass_matrix_amr(nodeMM[iLev], nmmc, nmmf, jHat[iLev], jhc,
                                     jhf, nodeB[iLev], nodeB0Lev, uBg[iLev],
                                     tc->get_dt(), iLev, solveFieldInCoMov,
                                     cellStatus);
    }
  }
  //////////////////////////////////////////////////////////////////////
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    jHat[iLev].SumBoundary(Geom(iLev).periodicity());
    sum_boundary_node_mm(iLev);
  }
  Vector<Real> invVol(n_lev());
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    invVol[iLev] = 1.0;
    for (int i = 0; i < nDim; ++i) {
      invVol[iLev] *= Geom(iLev).InvCellSize(i);
    }
  }

  for (int iLev = finest_level - 1; iLev >= 0; iLev--) {
    for (int i = finest_level; i > iLev; i--) {
      jHat[iLev].ParallelAdd(jhc[i][iLev]);
      nmmc[i][iLev] *= (invVol[iLev] / invVol[i]);
      nodeMM[iLev].ParallelAdd(nmmc[i][iLev]);
    }
  }
  for (int iLev = finest_level; iLev > 0; iLev--) {
    jHat[iLev].ParallelAdd(jhf[iLev - 1]);
    nmmf[iLev - 1] *= (invVol[iLev] / invVol[iLev - 1]);
    nodeMM[iLev].ParallelAdd(nmmf[iLev - 1]);
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    Real invVol = 1;
    for (int i = 0; i < nDim; ++i) {
      invVol *= Geom(iLev).InvCellSize(i);
    }
    jHat[iLev].mult(invVol, 0, jHat[iLev].nComp(), jHat[iLev].nGrow());
    jHat[iLev].FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
void Pic::sum_moments(bool updateDt) {
  std::string nameFunc = "Pic::sum_moments";
  if (isGridEmpty)
    return;

  timing_func(nameFunc);

  plasmaEnergy[iTot] = 0;
  for (int i = 0; i < nSpecies; ++i) {
    Real energy = parts[i]->sum_moments(nodePlasma[i], nodeB, tc->get_dt());
    plasmaEnergy[i] = energy;
    plasmaEnergy[iTot] += energy;

    // Undo the geometric dilution of the moments at the staircase surface of
    // the inner body (see #BODY) before they are converted to fluid moments:
    // the raw moments are linear in the particle weight, so rescaling them
    // leaves the bulk velocity and the pressure unchanged while the density
    // becomes the density of the plasma that is really there.
    if (useBody && !bodyIono.use) {
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        rescale_body_surface_nodes(nodePlasma[i][iLev], cellStatus[iLev],
                                   nodeStatus[iLev], Geom(iLev).Domain(), nDim);
      }
    }
  }

  auto convert_and_sum_moments = [&]() {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      nodePlasma[nSpecies][iLev].setVal(0.0);
    }

    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->convert_to_fluid_moments(nodePlasma[i]);
    }

    for (int i : kineticSpecies_) {
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        // nodePlasma[nSpecies] holds the sum of all kinetic species, i.e.
        // all species in full PIC and the ions in hybrid PIC.
        MultiFab::Add(nodePlasma[nSpecies][iLev], nodePlasma[i][iLev], 0, 0,
                      nMoments, nGst);
      }
    }

    for (int iLev = 0; iLev < n_lev(); iLev++) {
      nodePlasma[nSpecies][iLev].FillBoundary(Geom(iLev).periodicity());
    }

    if (finest_level > 0) {
      for (int iLev = 1; iLev < n_lev(); iLev++) {
        fill_fine_lev_bny_from_coarse(
            nodePlasma[nSpecies][iLev - 1], nodePlasma[nSpecies][iLev], 0,
            nodePlasma[nSpecies][iLev].nComp(), ref_ratio[iLev - 1],
            Geom(iLev - 1), Geom(iLev), node_status(iLev),
            node_bilinear_interp);

        if (useHybridPIC || useExplicitPIC) {
          fill_fine_lev_bny_from_coarse(nodePlasmaPrev[nSpecies][iLev - 1],
                                        nodePlasmaPrev[nSpecies][iLev], 0,
                                        nodePlasmaPrev[nSpecies][iLev].nComp(),
                                        ref_ratio[iLev - 1], Geom(iLev - 1),
                                        Geom(iLev), node_status(iLev),
                                        node_bilinear_interp);
        }
      }
    }

    calc_mach_number();
  };

  if (useHybridPIC) {
    convert_and_sum_moments();
  }

  if (updateDt) {
    Vector<Real> uMax(n_lev(), 0.0);
    Vector<Real> dxMin(n_lev());
    Vector<Real> dtMax(n_lev());
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      const auto& dx = Geom(iLev).CellSize();
      dxMin[iLev] = min(AMREX_D_DECL(dx[ix_], dx[iy_], dx[iz_]));

      // Only compute thermal velocity if CFL is active or if detailed report is
      // enabled. Avoids expensive particle traversal and MPI reductions for
      // fixed-dt compact runs.
      bool needThermalSpeed =
          (tc->get_cfl() > 0) || (doReport && !domainParameters.doCompact);

      if (needThermalSpeed) {
        if (fixedUMax >= 0) {
          uMax[iLev] = fixedUMax;
        } else {
          Vector<Real> uMaxSpecies(nSpecies, 0.0);
          for (int i = 0; i < nSpecies; ++i) {
            amrex::MultiFab& momMF = nodePlasma[i][iLev];
            uMaxSpecies[i] = parts[i]->calc_max_thermal_velocity(momMF);
          }
          // Reduce all species in a single MPI reduction instead of nSpecies
          // calls
          ParallelDescriptor::ReduceRealMax(uMaxSpecies.data(), nSpecies);

          for (int i = 0; i < nSpecies; ++i) {
            if (doReport && !domainParameters.doCompact) {
              Print() << printPrefix << std::setprecision(5) << "lev " << iLev
                      << " Species " << i << ": max(uth) = " << uMaxSpecies[i]
                      << std::endl;
            }

            if (uMaxSpecies[i] > uMax[iLev]) {
              uMax[iLev] = uMaxSpecies[i];
            }
          }
        }
      }

      dtMax[iLev] = (uMax[iLev] > 0.0) ? dxMin[iLev] / uMax[iLev]
                                       : std::numeric_limits<Real>::max();
    }

    if (useHybridPIC && (tc->get_cfl() > 0 || isAutoSubcycle)) {
      Real dtNext = tc->get_next_dt();
      int nSubNext = nBSubcycle;
      calc_hybrid_dt_and_subcycle(dtNext, nSubNext, doReport);
      if (tc->get_cfl() > 0) {
        tc->set_next_dt(dtNext);
        if (tc->get_dt() < 0) {
          tc->set_dt(dtNext);
        }
      }
      if (isAutoSubcycle) {
        nBSubcycle = nSubNext;
        nSubcycle = nSubNext;
      }
    } else if (tc->get_cfl() > 0) {
      Real dt0 = *std::min_element(dtMax.begin(), dtMax.end());
      Real dt = tc->get_cfl() * dt0;
      tc->set_next_dt(dt);

      if (tc->get_dt() < 0) {
        tc->set_dt(dt);
      }
    }

    maxCFL = 0.0;
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (dtMax[iLev] > 0.0 && dtMax[iLev] < std::numeric_limits<Real>::max()) {
        Real cflLev = tc->get_next_dt() / dtMax[iLev];
        if (cflLev > maxCFL)
          maxCFL = cflLev;
      }
    }

    if (doReport && !domainParameters.doCompact) {
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

  if (!useHybridPIC) {
    convert_and_sum_moments();
  }

  isMomentsUpdated = true;
}

//==========================================================
// Ma = u/vth
void Pic::calc_mach_number() {
  for (int iLev = 0; iLev < n_lev(); iLev++) {

    const auto& momentsMF = nodePlasma[nSpecies][iLev];
    for (MFIter mfi(momentsMF); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const Array4<const Real>& moments = momentsMF[mfi].array();
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
void Pic::calc_cost_per_cell() {
  const BalanceStrategy balanceStrategy = domainParameters.balanceStrategy;
  const int cellWeight = domainParameters.cellWeight;
  if (!isMomentsUpdated && balanceStrategy == BalanceStrategy::Particle) {
    sum_moments(false);
  }
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (balanceStrategy == BalanceStrategy::Cell) {
      cellCost[iLev].setVal(1.0);
    } else {
      average_node_to_cellcenter(cellCost[iLev], 0, nodePlasma[nSpecies][iLev],
                                 iNum_, cellCost[iLev].nComp(),
                                 cellCost[iLev].nGrow());
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

        if (balanceStrategy == BalanceStrategy::Particle ||
            balanceStrategy == BalanceStrategy::Hybrid) {
          // 1. The cells have been refined also allocated and use memory.
          // 2. It looks like these cells need calculations when
          // interpolating between levels.
          // 3. The number 10 is chosen by experience.
          cost(i, j, k) += cellWeight;
        }
      });
    }
  }
}

void Pic::advance_particles_and_sources(int pushStage) {
  particle_mover(pushStage);

  // Floating body (#FLOATINGBODY): update the sheath potential from the
  // currents of this push, or re-emit the net absorbed electron charge.
  update_floating_body();

  // Calling re_sampling after particle mover so that all the particles
  // outside the domain have been deleted.
  re_sampling();

  charge_exchange();

  // Apply chemical loss (recombination, etc.) by reducing particle
  // weights.  Must come before fill_source_particles() so that loss
  // and source are applied in the correct order within one step.
  if (source && source->use_loss_source()) {
    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->apply_loss(source, tc->get_dt());
    }
  }

  if (source) {
    fill_source_particles();
  }

  inject_particles_for_boundary_cells();

  // Ghost layer of the ionospheric body boundary (#BODYIONOSPHERE). The push
  // emptied it (the particles left in the body are absorbed); refilling it
  // here, like the ghost cells of the domain boundary, keeps the ionospheric
  // plasma of the boundary in the moments, the current and the mass matrix,
  // so that the field solver sees its charge (the ambipolar field).
  inject_body_ionosphere();

  // Ionospheric emission shell (#IONOSPHEREEMISSION).
  inject_ionosphere_emission();

  // Inject incoming flux at physical inflow faces before the moment deposit.
  // Kept outside inject_particles_for_boundary_cells to avoid t=0
  // pre-injection.
  if (usePIC) {
    const Real dt = tc->get_dt();
    for (int i : kineticSpecies_) {
      parts[i]->inject_flux_at_inflow_faces(dt);
      parts[i]->inject_flux_at_open_faces(dt, nodePlasma[i][0]);
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

  if (reportParticleQuality) {
    if (tc->get_cycle() % 20 == 0) {
      WriteParticleQualityToParaView();
    }
  }

  // Co-moving frame solver is full-PIC only.
  if (!useHybridPIC &&
      (solveFieldInCoMov || useUpwindB || (useUpwindE && cMaxE <= 0))) {
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

  isMomentsUpdated = false;

  if (useExplicitPIC) {
    advance_explicit_full_pic();
  } else if (solveEM) {
    advance_particles_and_sources(0);
    if (projectDownEmFields) {
      project_down_E();
    }
    update_B();
  } else if (useHybridPIC) {
    if (!usePredictorCorrector) {
      // Standard single-push hybrid path (1st-order coupled).
      if (isFirstHybridStep) {
        seed_first_hybrid_step();
      }

      advance_particles_and_sources(0);

      // Deposit the fresh ion moments (J^{n+1/2}) AFTER the particle push. The
      // previous deposit (J^{n-1/2}) is first saved into nodePlasmaPrev so the
      // Ohm's law can time-interpolate between the two (hstep scheme).
      save_current_moments_to_prev();
      sum_moments(false);
      smooth_moments();
      update_Pe_hybrid();
      update_B_hybrid();
      isFirstHybridStep = false;
    } else {
      // 2nd-order Predictor-Corrector path.
      if (isFirstHybridStep) {
        seed_first_hybrid_step();
      }

      // 1. Backup base particle state at t^n
      for (int i : kineticSpecies_) {
        parts[i]->backup_particle_state();
      }

      // 2. Backup base fields B^n and E^n
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        MultiFab::Copy(centerBstar[iLev], centerB[iLev], 0, 0, nDim3, nGst);
        MultiFab::Copy(nodeEth[iLev], nodeE[iLev], 0, 0, nDim3, nGst);
      }

      // 3. Predictor push using (E^n, B^n) [pushStage = 1]
      advance_particles_and_sources(1);

      // 4. Deposit trial moments & advance trial field to t^(n+1)
      save_current_moments_to_prev();
      sum_moments(false);
      smooth_moments();
      update_Pe_hybrid();
      update_B_hybrid();

      // 5. Construct time-centered fields at t^(n+1/2):
      //    nodeE = 0.5 * (E^n + E*^(n+1))
      //    nodeB = 0.5 * (B^n + B*^(n+1))
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        MultiFab::LinComb(nodeE[iLev], 0.5, nodeE[iLev], 0, 0.5, nodeEth[iLev],
                          0, 0, nDim3, nGst);
        MultiFab::LinComb(centerBstage[iLev], 0.5, centerB[iLev], 0, 0.5,
                          centerBstar[iLev], 0, 0, nDim3, nGst);
        average_center_to_node(centerBstage[iLev], nodeB[iLev]);
        nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
        if (iLev == 0) {
          apply_field_bc(nodeStatus[iLev], nodeB[iLev], 0, 3, &Pic::get_node_B,
                         iLev, true);
        }
      }

      // 6. Restore base particle state at t^n
      for (int i : kineticSpecies_) {
        parts[i]->restore_particle_state();
      }

      // 7. Corrector push using time-centered fields (nodeE, nodeB) [pushStage
      // = 2]
      advance_particles_and_sources(2);

      // 8. Restore B^n and advance to final B^(n+1) using accepted moments
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        MultiFab::Copy(centerB[iLev], centerBstar[iLev], 0, 0, nDim3, nGst);
      }
      sum_moments(false);
      smooth_moments();
      update_Pe_hybrid();
      update_B_hybrid();
      isFirstHybridStep = false;
    }
  } else {
    // Pure particle tracking (PT component or static-field particle runs).
    advance_particles_and_sources(0);
  }

  if ((solveEM || useExplicitPIC) && doCorrectDivE) {
    if (finest_level == 0) {
      divE_correction();
    } else {
      amr_divE_correction();
    }
  }

#ifdef _PC_COMPONENT_
  if (useHybridPIC) {
    // For hybrid PIC, accepted moments were already deposited and smoothed in
    // advance_hybrid(). We update hybrid dt and subcycle directly without a
    // redundant sum_moments() particle pass that would overwrite smoothed
    // moments.
    if (tc->get_cfl() > 0 || isAutoSubcycle) {
      Real dtNext = tc->get_next_dt();
      int nSubNext = nBSubcycle;
      calc_hybrid_dt_and_subcycle(dtNext, nSubNext, doReport);
      if (tc->get_cfl() > 0) {
        tc->set_next_dt(dtNext);
        tc->set_dt(dtNext);
      }
      if (isAutoSubcycle) {
        nBSubcycle = nSubNext;
        nSubcycle = nSubNext;
      }
    }
  } else {
    tc->set_dt(tc->get_next_dt());
    //  For PT simulations, moments are only useful for output. So, there is no
    //  need to call sum_moments() for every step.
    sum_moments(true);
  }
#endif

#ifdef _PT_COMPONENT_
  if (maxExchangeRatio > maxExchangeRatioLimit) {
    Real dtNew = tc->get_dt() * maxExchangeRatioLimit / maxExchangeRatio;
    tc->set_dt(dtNew);
    tc->set_next_dt(dtNew);
    Print() << printPrefix << " maxExchangeRatio = " << maxExchangeRatio
            << " maxExchangeRatioLimit = " << maxExchangeRatioLimit
            << " dt is reduced to " << tc->get_dt_si() << std::endl;
  } else {
    if (tc->get_dt() < tc->get_dt_max()) {
      // Increase dt if allowed.
      Real dtnow = tc->get_dt();
      Real dtNew = min(dtnow * maxExchangeRatioLimit / maxExchangeRatio,
                       tc->get_dt_max());

      if (dtNew > dtnow * (1 + 1e-6)) {
        tc->set_dt(dtNew);
        tc->set_next_dt(dtNew);
        Print() << printPrefix << " maxExchangeRatio = " << maxExchangeRatio
                << " maxExchangeRatioLimit = " << maxExchangeRatioLimit
                << " dt is increased to " << tc->get_dt_si() << std::endl;
      }
    }
  }
#endif

  if (doReport) {
    Real tEnd = second();
    Real nPoint = activeRegion.d_numPts();
    int nProc = ParallelDescriptor::NProcs();
    // The unit of the speed is (cell per processor per second)
    Real speed = nPoint / nProc / (tEnd - tStart);

    // speedNorm is a value obtained from tests.
    Real speedNorm = 1000;
    Real normSpeed = speed / speedNorm;

    if (domainParameters.doCompact) {
      Print() << "==== " << printPrefix << "Cycle " << tc->get_cycle()
              << " | t = " << std::setprecision(6) << tc->get_time_si()
              << " (s) | dt = " << std::setprecision(5) << tc->get_dt_si()
              << " (s)";
      if (maxCFL > 0.0) {
        Print() << " | CFL = " << std::setprecision(4) << maxCFL;
      }
      Print() << " | Speed = " << std::setprecision(2) << std::fixed
              << normSpeed << " ====" << std::endl;
      Print() << std::defaultfloat;
    } else {
      Print() << printPrefix
              << "Normalized PIC simulation speed = " << normSpeed
              << " (performance is good if the value >> 1 and bad if <<1 )"
              << std::endl;
    }
  }

  // Periodic load balance report at dnReportLB interval (default: 100 steps)
  bool doReportLB = (domainParameters.dnReportLB > 0 &&
                     tc->get_cycle() % domainParameters.dnReportLB == 0);
  if (doReportLB) {
    report_load_balance();
  }

  if (dnMemory > 0 && tc->get_cycle() % dnMemory == 0) {
    Print() << printPrefix << "Load balance before freeing memory:\n";
    report_load_balance();
    Print() << printPrefix << "Freeing memory...\n";
    free_memory();
    Print() << printPrefix << "Load balance after freeing memory:\n";
    report_load_balance();
  }
}

//==========================================================
void Pic::free_memory() {
  std::string nameFunc = "Pic::free_memory";
  timing_func(nameFunc);

  if (auto* p = dynamic_cast<amrex::CArena*>(amrex::The_Arena())) {
    p->freeUnused();
  }
  if (auto* p = dynamic_cast<amrex::CArena*>(amrex::The_Pinned_Arena())) {
    p->freeUnused();
  }

  amrex::FabArrayBase::flushTileArrayCache();

  amrex::Arena::PrintUsage();

  for (int i = 0; i < nSpecies; ++i) {
    amrex::Print() << "[Particles " << i << " Before ShrinkToFit]: ";
    parts[i]->PrintCapacity();
    parts[i]->ShrinkToFit();
    amrex::Print() << "[Particles " << i << " After ShrinkToFit]: ";
    parts[i]->PrintCapacity();
  }

#if defined(__linux__)
  // Force glibc to return free pages on the heap to the OS.
  // This is necessary because glibc often holds onto freed blocks internally.
  malloc_trim(0);
#endif
}

//==========================================================
void Pic::update_U0_E0() {
  std::string nameFunc = "Pic::update_U0_E0";
  timing_func(nameFunc);

  // Full-PIC only: eBg/uBg and nodePlasma are not allocated for hybrid.
  if (useHybridPIC)
    return;

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
          for (int iu = iUx_; iu <= iUz_; iu++) {
            // A (near-)denormal density overflows 1/rho.
            const Real u = arrMoments(i, j, k, iu) * invRho;
            arrU(i, j, k, iu - iUx_) = std::isfinite(u) ? u : 0.0;
          }
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

  maxExchangeRatio = 0;
  for (int i = 0; i < nSpecies; ++i) {
    Real rate = 0;
    parts[i]->charge_exchange(tc->get_dt(), stateOH, sourcePT2OH, source,
                              kineticSource, sourceParts, doSelectRegion,
                              product(nSourcePPC), rate);
    if (rate > maxExchangeRatio)
      maxExchangeRatio = rate;
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

//==========================================================
void Pic::init_regional_fields() {
  if (!useHybridPIC || (!hasRegionalResistivity_ && !hasRegionalHyper_))
    return;
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    init_regional_fields(iLev);
  }
}

//==========================================================
void Pic::init_regional_fields(int iLev) {
  if (!useHybridPIC || (!hasRegionalResistivity_ && !hasRegionalHyper_))
    return;
  fill_regional_resistivity_field(iLev);
  fill_regional_hyper_field(iLev);
}

//==========================================================
void Pic::fill_regional_resistivity_field(int iLev) {
  if (!hasRegionalResistivity_ || nodeEtaRegional.empty() ||
      !nodeEtaRegional[iLev].isDefined())
    return;

  nodeEtaRegional[iLev].setVal(0.0);
  const auto dx = Geom(iLev).CellSizeArray();
  const auto probLo = Geom(iLev).ProbLoArray();

  for (const auto& cfg : regionalResistivityConfigs) {
    if (cfg.etaCode <= 0.0)
      continue;
    Regions reg(regionShapes, cfg.regionStr);

    for (MFIter mfi(nodeEtaRegional[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.growntilebox();
      const Array4<Real>& arr = nodeEtaRegional[iLev][mfi].array();
      const auto lo = box.smallEnd();
      const auto hi = box.bigEnd();

#if (AMREX_SPACEDIM > 2)
      for (int k = lo[2]; k <= hi[2]; ++k) {
#else
      int k = 0;
      {
#endif
        for (int j = lo[1]; j <= hi[1]; ++j) {
          for (int i = lo[0]; i <= hi[0]; ++i) {
            Real xyz[nDim];
            xyz[0] = probLo[0] + i * dx[0];
            xyz[1] = probLo[1] + j * dx[1];
#if (AMREX_SPACEDIM > 2)
            xyz[2] = probLo[2] + k * dx[2];
#endif
            if (reg.is_inside(xyz)) {
              arr(i, j, k) = amrex::max(arr(i, j, k), cfg.etaCode);
            }
          }
        }
      }
    }
  }
  nodeEtaRegional[iLev].FillBoundary(Geom(iLev).periodicity());
}

//==========================================================
void Pic::fill_regional_hyper_field(int iLev) {
  if (!hasRegionalHyper_ || nodeEtaHyperRegional.empty() ||
      !nodeEtaHyperRegional[iLev].isDefined())
    return;

  nodeEtaHyperRegional[iLev].setVal(0.0);
  const auto dx = Geom(iLev).CellSizeArray();
  const auto probLo = Geom(iLev).ProbLoArray();

  for (const auto& cfg : regionalHyperResistivityConfigs) {
    const Real etaVal =
        (iLev < static_cast<int>(cfg.etaLev.size())) ? cfg.etaLev[iLev] : 0.0;
    if (etaVal <= 0.0)
      continue;
    Regions reg(regionShapes, cfg.regionStr);

    for (MFIter mfi(nodeEtaHyperRegional[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.growntilebox();
      const Array4<Real>& arr = nodeEtaHyperRegional[iLev][mfi].array();
      const auto lo = box.smallEnd();
      const auto hi = box.bigEnd();

#if (AMREX_SPACEDIM > 2)
      for (int k = lo[2]; k <= hi[2]; ++k) {
#else
      int k = 0;
      {
#endif
        for (int j = lo[1]; j <= hi[1]; ++j) {
          for (int i = lo[0]; i <= hi[0]; ++i) {
            Real xyz[nDim];
            xyz[0] = probLo[0] + i * dx[0];
            xyz[1] = probLo[1] + j * dx[1];
#if (AMREX_SPACEDIM > 2)
            xyz[2] = probLo[2] + k * dx[2];
#endif
            if (reg.is_inside(xyz)) {
              arr(i, j, k) = amrex::max(arr(i, j, k), etaVal);
            }
          }
        }
      }
    }
  }
  nodeEtaHyperRegional[iLev].FillBoundary(Geom(iLev).periodicity());
}

//==========================================================
void Pic::update_regional_hyper_grid_mode(Real dt) {
  if (!useHybridPIC || !hasRegionalHyper_ || dt <= 0.0)
    return;

  bool changed = false;
  const int iFinest = n_lev() - 1;
  const auto dxFine = Geom(iFinest).CellSizeArray();
  Real dxMinFine = dxFine[0];
  for (int d = 1; d < nDim; ++d)
    dxMinFine = amrex::min(dxMinFine, dxFine[d]);

  for (auto& cfg : regionalHyperResistivityConfigs) {
    if (cfg.mode == "grid" && cfg.ch > 0.0) {
      const Real etaHyper = fourPI * cfg.ch * std::pow(dxMinFine, 4) / dt;
      if (cfg.etaLev.size() != static_cast<amrex::Long>(n_lev_max()) ||
          std::abs(cfg.etaLev[0] - etaHyper) > 1e-14 * etaHyper) {
        cfg.etaLev.assign(n_lev_max(), etaHyper);
        changed = true;
      }
    }
  }

  if (changed) {
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      fill_regional_hyper_field(iLev);
    }
  }
}
