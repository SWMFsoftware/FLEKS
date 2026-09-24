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
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        average_center_to_node(centerB[iLev], nodeB[iLev]);
        nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
        if (iLev == 0) {
          apply_field_bc(nodeStatus[iLev], nodeB[iLev], 0, 3, &Pic::get_node_B,
                         iLev, true);
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

      distribute_FabArray(divB[iLev], cGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);
      distribute_FabArray(hypPhi[iLev], cGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);

      distribute_FabArray(centerDB[iLev], cGrids[iLev], DistributionMap(iLev),
                          nDim3, nGst, doMoveData);

      if (projectDownEmFields) {
        distribute_FabArray(projectScratchMF[iLev], nGrids[iLev],
                            DistributionMap(iLev), 3, 0, doMoveData);
      }

      if (!useExplicitPIC) {
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
    }
    if (useHybridPIC) {
      // Previous-step ion moments for the Ohm's-law interpolation.
      if (nodePlasmaPrev.empty()) {
        nodePlasmaPrev.resize(nSpecies + 1);
      }
      // Hyper-resistivity scratch: centerLapB = Laplacian(B); nodeHyperE
      // node-centred.
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

      // Hybrid-only node-grid previous-step moments (J^{n-1/2}), slim layout.
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
    // E current (jHat), and node-centred moments (nodePlasma): full-PIC only.
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

    fill_fine_lev_bny_from_coarse(centerB[iLev - 1], centerB[iLev], 0,
                                  centerB[iLev - 1].nComp(),
                                  ref_ratio[iLev - 1], Geom(iLev - 1),
                                  Geom(iLev), cell_status(iLev), cellInterp);
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

  if (source) {
    for (int i : kineticSpecies_) {
      parts[i]->add_particles_source(source, stateOH, tc->get_dt(), nSourcePPC,
                                     doSelectRegion, adaptiveSourcePPC);
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
        parts[i]->limit_weight(maxWeightRatio, parts[i]->is_neutral(),
                               pInfo.doPreSplitting);
      }
      parts[i]->split(reSamplingLowLimit, parts[i]->is_neutral(),
                      pInfo.doPreSplitting);
      parts[i]->merge(reSamplingHighLimit, pInfo.doPreSplitting);
    }
  }
}

//==========================================================
void Pic::particle_mover() {
  std::string nameFunc = "Pic::mover";

  timing_func(nameFunc);

  Real dt = tc->get_dt();
  Real dtnext = tc->get_next_dt();

  const Vector<MultiFab>& nodeEpush = useHybridPIC ? nodeE : nodeEth;

  for (int i : kineticSpecies_) {
    parts[i]->mover(nodeEpush, nodeB, eBg, uBg, dt, dtnext);
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

    for (int i = 0; i < nSpecies; ++i) {
      if (useExplicitPIC) {
        parts[i]->calc_jhat(jHat[iLev], nodeB[iLev], tc->get_dt());
      } else {
        parts[i]->calc_mass_matrix(nodeMM[iLev], jHat[iLev], nodeB[iLev],
                                   uBg[iLev], tc->get_dt(), iLev,
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

    if (doSmoothJ) {
      for (int icount = 0; icount < nSmoothJ; icount++) {
        smooth_multifab(jHat[iLev], iLev, icount % 2 + 1, coefSmoothJ);
      }
    }

    if (!useExplicitPIC) {
      sum_boundary_node_mm(iLev);
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
      parts[i]->calc_mass_matrix_amr(nodeMM[iLev], nmmc, nmmf, jHat[iLev], jhc,
                                     jhf, nodeB[iLev], uBg[iLev], tc->get_dt(),
                                     iLev, solveFieldInCoMov, cellStatus);
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
  }

  if (updateDt) {
    Vector<Real> uMax(n_lev(), 0.0);
    Vector<Real> dxMin(n_lev());
    Vector<Real> dtMax(n_lev());
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      const auto& dx = Geom(iLev).CellSize();
      dxMin[iLev] = min(AMREX_D_DECL(dx[ix_], dx[iy_], dx[iz_]));

      // Only compute thermal velocity if CFL is active or if detailed report is enabled.
      // Avoids expensive particle traversal and MPI reductions for fixed-dt compact runs.
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
          // Reduce all species in a single MPI reduction instead of nSpecies calls
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

    if (tc->get_cfl() > 0) {
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

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    nodePlasma[nSpecies][iLev].setVal(0.0);
  }

  for (int i = 0; i < nSpecies; ++i) {
    parts[i]->convert_to_fluid_moments(nodePlasma[i]);
  }

  for (int i : kineticSpecies_) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      // nodePlasma[nSpecies] holds the sum of all ion species.
      // kineticSpecies_ excludes the (implicit fluid) electron.
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
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);

      if (useHybridPIC) {
        fill_fine_lev_bny_from_coarse(
            nodePlasmaPrev[nSpecies][iLev - 1], nodePlasmaPrev[nSpecies][iLev],
            0, nodePlasmaPrev[nSpecies][iLev].nComp(), ref_ratio[iLev - 1],
            Geom(iLev - 1), Geom(iLev), node_status(iLev),
            node_bilinear_interp);
      }
    }
  }

  calc_mach_number();

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

  // Hybrid path: the particle Boris push happens BEFORE the moment deposit and
  // the Ohm's-law E computation. The push uses the nodeEth and B^n computed at
  // the end of the previous step's B update.
  if (useHybridPIC && isFirstHybridStep) {
    seed_first_hybrid_step();
  }

  particle_mover();

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

  // Inject incoming flux at physical inflow faces before the moment deposit.
  // Kept outside inject_particles_for_boundary_cells to avoid t=0
  // pre-injection.
  if (usePIC) {
    const Real dt = tc->get_dt();
    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->inject_flux_at_inflow_faces(dt);
    }
  }

  isMomentsUpdated = false;

  if (solveEM) {
    if (projectDownEmFields) {
      project_down_E();
    }
    update_B();
  } else if (useHybridPIC) {
    // Deposit the fresh ion moments (J^{n+1/2}) AFTER the particle push. The
    // previous deposit (J^{n-1/2}) is first saved into nodePlasmaPrev so the
    // Ohm's law can time-interpolate between the two (hstep scheme).
    save_current_moments_to_prev();
    sum_moments(false);
    smooth_moments();
    update_B_hybrid();
    isFirstHybridStep = false;
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
