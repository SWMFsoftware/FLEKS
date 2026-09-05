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
  // Per-species deposit targets; last entry = sum of all species.
  if (centerPlasma.empty()) {
    centerPlasma.resize(nSpecies + 1);
  }
  if (centerPlasmaPrev.empty()) {
    centerPlasmaPrev.resize(nSpecies + 1);
  }
  if (centerPlasmaSum.empty()) {
    centerPlasmaSum.resize(nSpecies + 1);
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

      if (!useExplicitPIC) {
        distribute_FabArray(nodeMM[iLev], nGrids[iLev], DistributionMap(iLev),
                            1, 1, doMoveData);
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

      // Time-averaged B scratch (cell + node), used when useAvgFieldB is set.
      distribute_FabArray(centerBavg[iLev], cGrids[iLev], DistributionMap(iLev),
                          3, nGst, doMoveData);
      distribute_FabArray(nodeBavg[iLev], nGrids[iLev], DistributionMap(iLev),
                          3, nGst, doMoveData);

      // rk3/rk4 persistent scratch: centerBstart = B_n; centerBstar =
      // (trial+B_n)/2.
      distribute_FabArray(centerBstart[iLev], cGrids[iLev],
                          DistributionMap(iLev), 3, nGst, doMoveData);
      distribute_FabArray(centerBstar[iLev], cGrids[iLev],
                          DistributionMap(iLev), 3, nGst, doMoveData);

      // Cell-centred hybrid solver fields.
      distribute_FabArray(centerEhybrid[iLev], cGrids[iLev],
                          DistributionMap(iLev), 3, nGst, doMoveData);
      distribute_FabArray(centerJ[iLev], cGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);
      distribute_FabArray(centerEstage[iLev], cGrids[iLev],
                          DistributionMap(iLev), 3, nGst, doMoveData);
      distribute_FabArray(centerHyperE[iLev], cGrids[iLev],
                          DistributionMap(iLev), 3, nGst, doMoveData);
      for (auto& pl : centerPlasmaSum) {
        if (pl.empty())
          pl.resize(n_lev_max());
        distribute_FabArray(pl[iLev], cGrids[iLev], DistributionMap(iLev),
                            nMoments, nGst, doMoveData);
      }
      for (auto& pl : centerPlasma) {
        if (pl.empty())
          pl.resize(n_lev_max());
        distribute_FabArray(pl[iLev], cGrids[iLev], DistributionMap(iLev),
                            nMoments, nGst, doMoveData);
      }
      for (auto& pl : centerPlasmaPrev) {
        if (pl.empty())
          pl.resize(n_lev_max());
        // Ohm's law reads only rho + 3 momentum, so stored slim (like
        // nodePlasmaPrev).
        distribute_FabArray(pl[iLev], cGrids[iLev], DistributionMap(iLev),
                            nHybridMomentsComps, nGst, doMoveData);
      }
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

    // mMach: node grid for full-PIC, cell grid for hybrid.
    if (useHybridPIC) {
      distribute_FabArray(mMach[iLev], cGrids[iLev], DistributionMap(iLev), 1,
                          nGst, doMoveData);
    } else {
      distribute_FabArray(mMach[iLev], nGrids[iLev], DistributionMap(iLev), 1,
                          nGst, doMoveData);
    }

    // Co-moving frame fields (eBg/uBg), div(E) mass matrix (centerMM), implicit
    // E current (jHat), and node-centred moments (nodePlasma): full-PIC only.
    if (!useHybridPIC) {
      distribute_FabArray(eBg[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);

      distribute_FabArray(uBg[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);

      distribute_FabArray(centerMM[iLev], cGrids[iLev], DistributionMap(iLev),
                          1, nGst, doMoveData);

      distribute_FabArray(jHat[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                          nGst, doMoveData);

      for (auto& pl : nodePlasma) {
        if (pl.empty())
          pl.resize(n_lev_max());
        distribute_FabArray(pl[iLev], nGrids[iLev], DistributionMap(iLev),
                            nMoments, nGst, doMoveData);
      }
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

  // Initial-condition / restart E is node-centred (nodeE). centerEhybrid is
  // seeded from it by averaging the node values to the cell centres, which
  // plays the role of E0 for the very first hybrid particle Boris push.
  if (useHybridPIC) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      average_node_to_center(nodeE[iLev], centerEhybrid[iLev]);
      centerEhybrid[iLev].FillBoundary(Geom(iLev).periodicity());
      // Match full-PIC nodeE: the same face type closes the ghost ring of
      // the initial cell-centred E as well.
      apply_field_bc(cellStatus[iLev], centerEhybrid[iLev], 0,
                     centerEhybrid[iLev].nComp(), &Pic::get_center_E, iLev,
                     false);
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
      // Use the time-averaged B in the Boris half-stage position push when
      // enabled (falls back to the instantaneous B before the first average is
      // initialised).
      const auto& nodeBhalf =
          (useAvgFieldB && isBavgInit) ? nodeBavg[iLev] : nodeB[iLev];
      parts[i]->update_position_to_half_stage(nodeEth[iLev], nodeBhalf,
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
      if (!pInfo.doPreSplitting) {
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

  // Time-averaged B when enabled.
  const Vector<MultiFab>& nodeBpush =
      (useAvgFieldB && isBavgInit) ? nodeBavg : nodeB;
  const Vector<MultiFab>& centerBpush =
      (useAvgFieldB && isBavgInit) ? centerBavg : centerB;
  const Vector<MultiFab>& nodeEpush = nodeEth;
  if (useHybridPIC) {
    for (int i : kineticSpecies_) {
      parts[i]->mover_cell_centered(centerEhybrid, centerBpush, eBg, uBg, dt,
                                    dtnext);
    }
  } else {
    for (int i : kineticSpecies_) {
      parts[i]->mover(nodeEpush, nodeBpush, eBg, uBg, dt, dtnext);
    }
  }

  for (int i : kineticSpecies_) {
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

    if (useExplicitPIC && nSpecies > 0) {
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
void Pic::sum_moments(bool updateDt) {
  std::string nameFunc = "Pic::sum_moments";
  if (isGridEmpty)
    return;

  timing_func(nameFunc);

  plasmaEnergy[iTot] = 0;
  for (int i = 0; i < nSpecies; ++i) {
    Real energy = 0.0;
    if (useHybridPIC) {
      // Cell-centred moment deposit into centerPlasma[i].
      energy = parts[i]->sum_moments_cell_centered(centerPlasma[i]);
    } else {
      energy = parts[i]->sum_moments(nodePlasma[i], nodeB, tc->get_dt());
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
          amrex::MultiFab& momMF =
              useHybridPIC ? centerPlasma[i][iLev] : nodePlasma[i][iLev];
          Real uMaxSpecies = parts[i]->calc_max_thermal_velocity(momMF);
          ParallelDescriptor::ReduceRealMax(uMaxSpecies);

          if (doReport) {
            Print() << printPrefix << std::setprecision(5) << "lev " << iLev
                    << " Species " << i << ": max(uth) = " << uMaxSpecies
                    << std::endl;
          }

          if (uMaxSpecies > uMax[iLev]) {
            uMax[iLev] = uMaxSpecies;
          }
        }

        // Generic override of the CFL signal speed (e.g. the old TopHat
        // option used a fixed value of 1.0). A negative fixedUMax keeps the
        // particle-thermal-velocity estimate.
        if (fixedUMax >= 0) {
          uMax[iLev] = fixedUMax;
        }

        dtMax[iLev] = (uMax[iLev] > 0.0) ? dxMin[iLev] / uMax[iLev]
                                         : std::numeric_limits<Real>::max();
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

  if (useHybridPIC) {
    // Cell-centred hybrid moments: sum the per-species deposits into
    // centerPlasmaSum and sync the nodePlasma output mirror (once per step,
    // so the plot / restart / tracker path that reads nodePlasma sees correct
    // data -- the node-sync bridge of the hybrid solver).
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      centerPlasmaSum[nSpecies][iLev].setVal(0.0);
    }

    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->convert_to_fluid_moments(centerPlasma[i]);
    }

    for (int i : kineticSpecies_) {
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        // centerPlasmaSum[nSpecies] holds the sum of all kinetic-ion species.
        // kineticSpecies_ excludes the (implicit fluid) electron.
        MultiFab::Add(centerPlasmaSum[nSpecies][iLev], centerPlasma[i][iLev], 0,
                      0, nMoments, nGst);
      }
    }

    // Fill ghost cells for centerPlasmaSum so that assemble_ohm_E's 2-dx
    // stencil reads valid data at MPI boundaries and coarse-fine interfaces.
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      centerPlasmaSum[nSpecies][iLev].FillBoundary(Geom(iLev).periodicity());
    }

    // Mirror ion moments into the physical-wall ghost cells for smooth
    // pressure-gradient / Hall stencils at a wall.
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      apply_centerPlasma_BC(cell_status(iLev), centerPlasmaSum[nSpecies][iLev],
                            iLev);
      apply_centerPlasma_BC(cell_status(iLev), centerPlasmaPrev[nSpecies][iLev],
                            iLev);
    }

    // Fill coarse-fine interface ghost cells for centerPlasmaSum and
    // centerPlasmaPrev on fine levels from the coarse level. Without this,
    // assemble_ohm_E's pressure-gradient stencil reads zero/stale ghost cells
    // at the coarse-fine interface, causing incorrect E fields and runaway
    // particle heating.
    if (finest_level > 0 && useHybridPIC) {
      auto& cellInterp = *get_cell_interp();
      for (int iLev = 1; iLev < n_lev(); iLev++) {
        fill_fine_lev_bny_from_coarse(
            centerPlasmaSum[nSpecies][iLev - 1],
            centerPlasmaSum[nSpecies][iLev], 0,
            centerPlasmaSum[nSpecies][iLev].nComp(), ref_ratio[iLev - 1],
            Geom(iLev - 1), Geom(iLev), cell_status(iLev), cellInterp);

        fill_fine_lev_bny_from_coarse(
            centerPlasmaPrev[nSpecies][iLev - 1],
            centerPlasmaPrev[nSpecies][iLev], 0,
            centerPlasmaPrev[nSpecies][iLev].nComp(), ref_ratio[iLev - 1],
            Geom(iLev - 1), Geom(iLev), cell_status(iLev), cellInterp);
      }
    }
  } else {
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
  }

  if (!useHybridPIC) {
    calc_mach_number();
  }

  isMomentsUpdated = true;
}

//==========================================================
// Ma = u/vth
void Pic::calc_mach_number() {
  for (int iLev = 0; iLev < n_lev(); iLev++) {

    // Hybrid: Mach number from the live cell-centred summed ion moments
    // (centerPlasmaSum[nSpecies]). Full-PIC: from the node-centred
    // nodePlasma[nSpecies]. mMach is allocated on the matching grid.
    const auto& momentsMF = useHybridPIC ? centerPlasmaSum[nSpecies][iLev]
                                         : nodePlasma[nSpecies][iLev];
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
      // Balance by particles or hybrid. Hybrid: cellCost and the summed ion
      // moments are both cell-centred, so copy iNum_ directly (no cell->node->
      // cell roundtrip and no need to materialize the deferred nodePlasma
      // mirror). Full-PIC: average the node-centred nodePlasma particle count
      // to the cell grid.
      if (useHybridPIC) {
        MultiFab::Copy(cellCost[iLev], centerPlasmaSum[nSpecies][iLev], iNum_,
                       0, cellCost[iLev].nComp(), cellCost[iLev].nGrow());
      } else {
        average_node_to_cellcenter(
            cellCost[iLev], 0, nodePlasma[nSpecies][iLev], iNum_,
            cellCost[iLev].nComp(), cellCost[iLev].nGrow());
      }
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
    Print() << printPrefix
            << "Normalized PIC simulation speed = " << speed / speedNorm
            << " (performance is good if the value >> 1 and bad if <<1 )"
            << std::endl;

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
