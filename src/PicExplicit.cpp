#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
// Seed nodePlasmaPrev on the first explicit step, where there is no previous
// deposit: initialize it from the current deposit so the time interpolation
// evaluates to the initial state for that single step.
void Pic::seed_first_explicit_step() {
  std::string nameFunc = "Pic::seed_first_explicit_step";
  timing_func(nameFunc);

  save_current_moments_to_prev();
  isFirstExplicitStep = false;
}

//==========================================================
// Advance B by -0.5 * dtSub * curl(E) across all levels and average to nodeB.
void Pic::advance_B_half_step(Real dtSub, Real tSub) {
  std::string nameFunc = "Pic::advance_B_half_step";
  timing_func(nameFunc);

  const Real tEval = (tSub >= 0.0) ? tSub : (tc ? tc->get_time() : 0.0);

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    if (centerDB[iLev].empty()) {
      distribute_FabArray(centerDB[iLev], cGrids[iLev], DistributionMap(iLev),
                          nDim3, nGst);
    }
    curl_node_to_center(nodeE[iLev], centerDB[iLev], Geom(iLev).InvCellSize());

    if (is_body_interior_frozen()) {
      mask_body_interior(centerDB[iLev], cellStatus[iLev]);
    }

    MultiFab::Saxpy(centerB[iLev], -0.5 * dtSub, centerDB[iLev], 0, 0, nDim3,
                    centerB[iLev].nGrow());
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  if (projectDownEmFields && finest_level > 0) {
    for (int iLev = finest_level; iLev > 0; --iLev) {
      average_down(centerB[iLev], centerB[iLev - 1], 0, 3, ref_ratio[0]);
    }
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
                     &Pic::get_center_B, iLev, true, 0.5 * dtSub, tEval);
    } else {
      fill_fine_lev_bny_from_coarse(
          centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          *get_cell_interp());
    }

    average_center_to_node(centerB[iLev], nodeB[iLev]);
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(nodeStatus[iLev], nodeB[iLev], 0, nodeB[iLev].nComp(),
                     &Pic::get_node_B, iLev, true, 0.5 * dtSub, tEval);
    } else {
      fill_fine_lev_bny_from_coarse(
          nodeB[iLev - 1], nodeB[iLev], 0, nodeB[iLev].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          node_bilinear_interp);
    }
  }
}

//==========================================================
// Advance nodeE on level iLev by dtSub * (curl(B) - 4*pi*J).
void Pic::advance_E_step(Real dtSub, const MultiFab& currentMF, int iLev,
                         Real tSub) {
  std::string nameFunc = "Pic::advance_E_step";
  timing_func(nameFunc);

  const Real tEval = (tSub >= 0.0) ? tSub : (tc ? tc->get_time() : 0.0);

  curl_center_to_node(centerB[iLev], nodeCurlB[iLev], Geom(iLev).InvCellSize());

  // E += dtSub * curl(B)
  MultiFab::Saxpy(nodeE[iLev], dtSub, nodeCurlB[iLev], 0, 0, nDim3,
                  nodeE[iLev].nGrow());

  // E -= 4*pi * dtSub * J
  MultiFab::Saxpy(nodeE[iLev], -fourPI * dtSub, currentMF, 0, 0, nDim3,
                  nodeE[iLev].nGrow());

  nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(nodeStatus[iLev], nodeE[iLev], 0, nDim3, &Pic::get_node_E,
                   iLev, false, dtSub, tEval);
  } else {
    fill_fine_lev_bny_from_coarse(
        nodeE[iLev - 1], nodeE[iLev], 0, nDim3, ref_ratio[iLev - 1],
        Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
  }

  apply_body_E_bc(nodeE[iLev], iLev);
}

//==========================================================
// Subcycled explicit Maxwell advance for B and E.
void Pic::update_EB_explicit_leapfrog(Real dt, int nSub) {
  std::string nameFunc = "Pic::update_EB_explicit_leapfrog";
  timing_func(nameFunc);

  const Real subDt = dt / nSub;
  const Real t0 = tc ? tc->get_time() : 0.0;

  for (int iSub = 0; iSub < nSub; ++iSub) {
    // 1. Advance B half-step: B^m -> B^{m+1/2}
    const Real tB_half = t0 + (iSub + 0.5) * subDt;
    advance_B_half_step(subDt, tB_half);

    // 2. Advance E full-step: E^m -> E^{m+1} using B^{m+1/2} and interpolated
    // current
    const Real tE_full = t0 + (iSub + 1.0) * subDt;
    const Real hstep = (nSub == 1) ? 0.5 : (iSub + 0.5) / nSub;
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      if (nSub == 1) {
        MultiFab::LinComb(nodeJ[iLev], 0.5, nodePlasmaPrev[nSpecies][iLev],
                          iMx_, 0.5, nodePlasma[nSpecies][iLev], iMx_, 0, nDim3,
                          nGst);
      } else {
        const Real wPrev = 1.0 - hstep;
        const Real wCurr = hstep;
        MultiFab::LinComb(nodeJ[iLev], wPrev, nodePlasmaPrev[nSpecies][iLev],
                          iMx_, wCurr, nodePlasma[nSpecies][iLev], iMx_, 0,
                          nDim3, nGst);
      }
      advance_E_step(subDt, nodeJ[iLev], iLev, tE_full);
    }

    if (projectDownEmFields && finest_level > 0) {
      project_down_E();
    }

    // 3. Advance B second half-step: B^{m+1/2} -> B^{m+1}
    const Real tB_full = t0 + (iSub + 1.0) * subDt;
    advance_B_half_step(subDt, tB_full);
  }
}

//==========================================================
// 2nd-order Predictor-Corrector path for explicit full PIC.
void Pic::advance_explicit_pc(Real dt) {
  std::string nameFunc = "Pic::advance_explicit_pc";
  timing_func(nameFunc);

  if (isFirstExplicitStep) {
    seed_first_explicit_step();
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

  // 4. Deposit trial moments & advance trial fields to t^(n+1)
  save_current_moments_to_prev();
  sum_moments(false);
  if (doSmoothJ || doSmoothMoments) {
    smooth_moments();
  }
  update_EB_explicit_leapfrog(dt, nEMSubcycle);

  // 5. Construct time-centered fields at t^(n+1/2)
  const Real t0 = tc ? tc->get_time() : 0.0;
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::LinComb(nodeE[iLev], 0.5, nodeE[iLev], 0, 0.5, nodeEth[iLev], 0,
                      0, nDim3, nGst);
    MultiFab::LinComb(centerBstage[iLev], 0.5, centerB[iLev], 0, 0.5,
                      centerBstar[iLev], 0, 0, nDim3, nGst);
    average_center_to_node(centerBstage[iLev], nodeB[iLev]);
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(nodeStatus[iLev], nodeB[iLev], 0, 3, &Pic::get_node_B,
                     iLev, true, -1.0, t0 + 0.5 * dt);
    } else {
      fill_fine_lev_bny_from_coarse(
          nodeB[iLev - 1], nodeB[iLev], 0, 3, ref_ratio[iLev - 1],
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    }
  }

  // 6. Restore base particle state at t^n
  for (int i : kineticSpecies_) {
    parts[i]->restore_particle_state();
  }

  // 7. Corrector push using time-centered fields (nodeE, nodeB) [pushStage = 2]
  advance_particles_and_sources(2);

  // 8. Restore B^n and advance to final (B^{n+1}, E^{n+1}) using accepted
  // moments
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::Copy(centerB[iLev], centerBstar[iLev], 0, 0, nDim3, nGst);
    MultiFab::Copy(nodeE[iLev], nodeEth[iLev], 0, 0, nDim3, nGst);
  }
  sum_moments(false);
  if (doSmoothJ || doSmoothMoments) {
    smooth_moments();
  }
  update_EB_explicit_leapfrog(dt, nSubcycle);
}

//==========================================================
// Main advance for explicit full PIC.
void Pic::advance_explicit_full_pic() {
  std::string nameFunc = "Pic::advance_explicit_full_pic";
  timing_func(nameFunc);

  const Real dt = tc->get_dt();
  const Real t0 = tc ? tc->get_time() : 0.0;

  const bool usePC =
      (fieldIntegrator == "rk4-pc" || fieldIntegrator == "pc-rk4" ||
       fieldIntegrator == "ssprk3-pc" || fieldIntegrator == "pc-ssprk3" ||
       fieldIntegrator == "pc" || fieldIntegrator == "rk2");

  if (usePC) {
    advance_explicit_pc(dt);
  } else {
    // Standard leapfrog / Strang splitting
    if (isFirstExplicitStep) {
      seed_first_explicit_step();
    }

    if (nSubcycle == 1) {
      // 1. Advance B half-step: B^n -> B^{n+1/2}
      advance_B_half_step(dt, t0 + 0.5 * dt);

      // 2. Boris push particles using (E^n, B^{n+1/2})
      advance_particles_and_sources(0);

      // 3. Deposit moments at t^{n+1}
      save_current_moments_to_prev();
      sum_moments(false);
      if (doSmoothJ || doSmoothMoments) {
        smooth_moments();
      }

      // 4. Advance E full-step using time-centered current J^{n+1/2} = 0.5 *
      // (J^n + J^{n+1})
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        MultiFab::LinComb(nodeJ[iLev], 0.5, nodePlasmaPrev[nSpecies][iLev],
                          iMx_, 0.5, nodePlasma[nSpecies][iLev], iMx_, 0, nDim3,
                          nGst);
        advance_E_step(dt, nodeJ[iLev], iLev, t0 + dt);
      }

      if (projectDownEmFields && finest_level > 0) {
        project_down_E();
      }

      // 5. Advance B second half-step: B^{n+1/2} -> B^{n+1} using E^{n+1}
      advance_B_half_step(dt, t0 + dt);
    } else {
      // Subcycled leapfrog: push particles first with (E^n, B^n)
      advance_particles_and_sources(0);

      save_current_moments_to_prev();
      sum_moments(false);
      if (doSmoothJ || doSmoothMoments) {
        smooth_moments();
      }

      update_EB_explicit_leapfrog(dt, nSubcycle);
    }
  }

  // Div(E) correction (Gauss's law cleaning) if configured
  if (doCorrectDivE) {
    if (finest_level == 0) {
      divE_correction();
    } else {
      amr_divE_correction();
    }
  }
}
