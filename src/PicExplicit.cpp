#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
// Deposit the charge current of the initial particle state (J^0) on the
// first predictor-corrector step, where there is no previous deposit.
void Pic::seed_first_explicit_step() {
  std::string nameFunc = "Pic::seed_first_explicit_step";
  timing_func(nameFunc);

  deposit_explicit_current(nodeJexpl);
  isFirstExplicitStep = false;
}

//==========================================================
// The moments deposit q_p v in the momentum slots, and
// convert_to_fluid_moments() turns them into the mass flux (rho u)_s. Maxwell's
// equations need the charge current, so each species enters with q_s/m_s.
void Pic::deposit_explicit_current(Vector<MultiFab>& J) {
  std::string nameFunc = "Pic::deposit_explicit_current";
  timing_func(nameFunc);

  sum_moments(false);

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    J[iLev].setVal(0.0);
    for (int i : kineticSpecies_) {
      const Real qom = parts[i]->get_charge() / parts[i]->get_mass();
      MultiFab::Saxpy(J[iLev], qom, nodePlasma[i][iLev], iMx_, 0, nDim3,
                      get_n_ghost());
    }
    J[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  for (int iLev = 1; iLev < n_lev(); ++iLev) {
    fill_fine_lev_bny_from_coarse(
        J[iLev - 1], J[iLev], 0, nDim3, refRatio(iLev - 1), Geom(iLev - 1),
        Geom(iLev), node_status(iLev), node_bilinear_interp);
  }

  smooth_explicit_current(J);
}

//==========================================================
// #SMOOTHJ smoothing of the explicit current, with the same passes, walls
// and body treatment that smooth_moments() applies to the moments.
void Pic::smooth_explicit_current(Vector<MultiFab>& J) {
  const bool shouldSmooth =
      (doSmoothJ && nSmoothJ > 0) || (doSmoothMoments && nSmoothMoments > 0);
  if (!shouldSmooth)
    return;

  std::string nameFunc = "Pic::smooth_explicit_current";
  timing_func(nameFunc);

  const int nSmooth = doSmoothJ ? nSmoothJ : nSmoothMoments;
  const Real coefSmooth = doSmoothJ ? coefSmoothJ : coefSmoothMoments;

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    for (int icount = 0; icount < nSmooth; ++icount) {
      smooth_multifab(J[iLev], iLev, 1, coefSmooth);
      if (iLev > 0) {
        fill_fine_lev_bny_from_coarse(
            J[iLev - 1], J[iLev], 0, nDim3, refRatio(iLev - 1), Geom(iLev - 1),
            Geom(iLev), node_status(iLev), node_bilinear_interp);
      }
    }
    J[iLev].FillBoundary(Geom(iLev).periodicity());

    // No normal current through a conducting wall.
    if (!Geom(iLev).isAllPeriodic()) {
      const BoundaryBounds bnd(Geom(iLev), J[iLev].boxArray().ixType(),
                               &bcField);
      for (MFIter mfi(J[iLev]); mfi.isValid(); ++mfi) {
        const Box& bxValid = mfi.validbox();
        Array4<Real> const& arr = J[iLev][mfi].array();
        const Dim3 vLo = bxValid.smallEnd().dim3();
        const Dim3 vHi = bxValid.bigEnd().dim3();
        ParallelFor(bxValid, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const int ijk[3] = { i, j, k };
          const int vLoArr[3] = { vLo.x, vLo.y, vLo.z };
          const int vHiArr[3] = { vHi.x, vHi.y, vHi.z };
          for (int d = 0; d < nDim; ++d) {
            if (!bnd.isNode[d])
              continue;
            const bool onWall = ((bnd.bcLo[d] == FieldBC::conducting) &&
                                 (ijk[d] == bnd.loBnd[d])) ||
                                ((bnd.bcHi[d] == FieldBC::conducting) &&
                                 (ijk[d] == bnd.hiBnd[d]));
            if (!onWall)
              continue;
            bool inValid = true;
            for (int od = 0; od < nDim; ++od) {
              if (od != d && (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
                inValid = false;
                break;
              }
            }
            if (inValid)
              arr(i, j, k, d) = 0.0;
          }
        });
      }
    }

    if (use_body()) {
      mask_body(J[iLev], node_status(iLev));
    }
  }

  for (int iLev = 1; iLev < n_lev(); ++iLev) {
    fill_fine_lev_bny_from_coarse(
        J[iLev - 1], J[iLev], 0, nDim3, refRatio(iLev - 1), Geom(iLev - 1),
        Geom(iLev), node_status(iLev), node_bilinear_interp);
  }
}

//==========================================================
// Advance B by -0.5 * dtSub * curl(E) across all levels and average to nodeB.
void Pic::advance_B_half_step(Real dtSub, Real tSub) {
  std::string nameFunc = "Pic::advance_B_half_step";
  timing_func(nameFunc);

  const Real tEval = (tSub >= 0.0) ? tSub : (tc ? tc->get_time() : 0.0);

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    if (centerDB[iLev].empty()) {
      distribute_FabArray(centerDB[iLev], boxArray(iLev), DistributionMap(iLev),
                          nDim3, get_n_ghost());
    }
    curl_node_to_center(nodeE[iLev], centerDB[iLev], Geom(iLev).InvCellSize());

    if (is_body_interior_frozen()) {
      mask_body_interior(centerDB[iLev], cell_status(iLev));
    }

    MultiFab::Saxpy(centerB[iLev], -0.5 * dtSub, centerDB[iLev], 0, 0, nDim3,
                    centerB[iLev].nGrow());
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  if (projectDownEmFields && finestLevel() > 0) {
    for (int iLev = finestLevel(); iLev > 0; --iLev) {
      average_down(centerB[iLev], centerB[iLev - 1], 0, 3, refRatio(0));
    }
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(cell_status(iLev), centerB[iLev], 0, centerB[iLev].nComp(),
                     &Pic::get_center_B, iLev, true, 0.5 * dtSub, tEval);
    } else {
      fill_fine_lev_bny_from_coarse(
          centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
          refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          *get_cell_interp());
    }

    average_center_to_node(centerB[iLev], nodeB[iLev]);
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(node_status(iLev), nodeB[iLev], 0, nodeB[iLev].nComp(),
                     &Pic::get_node_B, iLev, true, 0.5 * dtSub, tEval);
    } else {
      fill_fine_lev_bny_from_coarse(nodeB[iLev - 1], nodeB[iLev], 0,
                                    nodeB[iLev].nComp(), refRatio(iLev - 1),
                                    Geom(iLev - 1), Geom(iLev),
                                    node_status(iLev), node_bilinear_interp);
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
    apply_field_bc(node_status(iLev), nodeE[iLev], 0, nDim3, &Pic::get_node_E,
                   iLev, false, dtSub, tEval);
  } else {
    fill_fine_lev_bny_from_coarse(
        nodeE[iLev - 1], nodeE[iLev], 0, nDim3, refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
  }

  apply_body_E_bc(nodeE[iLev], iLev);
}

//==========================================================
// Subcycled explicit Maxwell advance for B and E.
void Pic::update_EB_explicit_leapfrog(Real dt, int nSub,
                                      const Vector<MultiFab>& J0,
                                      const Vector<MultiFab>& J1,
                                      bool interpolateInTime) {
  std::string nameFunc = "Pic::update_EB_explicit_leapfrog";
  timing_func(nameFunc);

  const Real subDt = dt / nSub;
  const Real t0 = tc ? tc->get_time() : 0.0;

  for (int iSub = 0; iSub < nSub; ++iSub) {
    // 1. Advance B half-step: B^m -> B^{m+1/2}
    const Real tB_half = t0 + (iSub + 0.5) * subDt;
    advance_B_half_step(subDt, tB_half);

    // 2. Advance E full-step: E^m -> E^{m+1} using B^{m+1/2} and the current
    // at the sub-step mid time.
    const Real tE_full = t0 + (iSub + 1.0) * subDt;
    const Real w1 = interpolateInTime ? (iSub + 0.5) / nSub : 0.5;
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      MultiFab::LinComb(nodeJ[iLev], 1.0 - w1, J0[iLev], 0, w1, J1[iLev], 0, 0,
                        nDim3, get_n_ghost());
      advance_E_step(subDt, nodeJ[iLev], iLev, tE_full);
    }

    if (projectDownEmFields && finestLevel() > 0) {
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
    MultiFab::Copy(centerBstar[iLev], centerB[iLev], 0, 0, nDim3,
                   get_n_ghost());
    MultiFab::Copy(nodeEth[iLev], nodeE[iLev], 0, 0, nDim3, get_n_ghost());
  }

  // 3. Predictor push using (E^n, B^n) [pushStage = 1]
  advance_particles_and_sources(1);

  // 4. Deposit the trial current J^(n+1) & advance trial fields to t^(n+1).
  // Velocities live at whole steps here, so J^n and J^(n+1) are interpolated
  // in time.
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::Copy(nodeJexplPrev[iLev], nodeJexpl[iLev], 0, 0, nDim3,
                   get_n_ghost());
  }
  deposit_explicit_current(nodeJexpl);
  update_EB_explicit_leapfrog(dt, nEMSubcycle, nodeJexplPrev, nodeJexpl, true);

  // 5. Construct time-centered fields at t^(n+1/2)
  const Real t0 = tc ? tc->get_time() : 0.0;
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::LinComb(nodeE[iLev], 0.5, nodeE[iLev], 0, 0.5, nodeEth[iLev], 0,
                      0, nDim3, get_n_ghost());
    MultiFab::LinComb(centerBstage[iLev], 0.5, centerB[iLev], 0, 0.5,
                      centerBstar[iLev], 0, 0, nDim3, get_n_ghost());
    average_center_to_node(centerBstage[iLev], nodeB[iLev]);
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(node_status(iLev), nodeB[iLev], 0, 3, &Pic::get_node_B,
                     iLev, true, -1.0, t0 + 0.5 * dt);
    } else {
      fill_fine_lev_bny_from_coarse(
          nodeB[iLev - 1], nodeB[iLev], 0, 3, refRatio(iLev - 1),
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
    MultiFab::Copy(centerB[iLev], centerBstar[iLev], 0, 0, nDim3,
                   get_n_ghost());
    MultiFab::Copy(nodeE[iLev], nodeEth[iLev], 0, 0, nDim3, get_n_ghost());
  }
  deposit_explicit_current(nodeJexpl);
  update_EB_explicit_leapfrog(dt, nSubcycle, nodeJexplPrev, nodeJexpl, true);
}

//==========================================================
// Main advance for explicit full PIC.
void Pic::advance_explicit_full_pic() {
  std::string nameFunc = "Pic::advance_explicit_full_pic";
  timing_func(nameFunc);

  const Real dt = tc->get_dt();

  const bool usePC =
      (fieldIntegrator == "rk4-pc" || fieldIntegrator == "pc-rk4" ||
       fieldIntegrator == "ssprk3-pc" || fieldIntegrator == "pc-ssprk3" ||
       fieldIntegrator == "pc" || fieldIntegrator == "rk2");

  if (usePC) {
    advance_explicit_pc(dt);
  } else {
    // Leapfrog: x^n, v^{n-1/2}, E^n, B^n -> x^{n+1}, v^{n+1/2}, E^{n+1},
    // B^{n+1}.
    // 1. Kick v^{n-1/2} -> v^{n+1/2} with (E^n, B^n); positions stay at x^n.
    particle_mover(PicParticles::kPushKickOnly);

    // 2. Deposit q v^{n+1/2} at x^n, drift to x^{n+1}, deposit again. Their
    // average is the time-centered current J^{n+1/2}.
    deposit_explicit_current(nodeJexplPrev);
    particle_drift();
    finish_particle_push();
    deposit_explicit_current(nodeJexpl);

    // 3. B half-step, E full step with J^{n+1/2}, B half-step (subcycled).
    update_EB_explicit_leapfrog(dt, nSubcycle, nodeJexplPrev, nodeJexpl, false);
    isFirstExplicitStep = false;
  }

  // Div(E) correction (Gauss's law cleaning) if configured
  if (doCorrectDivE) {
    if (finestLevel() == 0) {
      divE_correction();
    } else {
      amr_divE_correction();
    }
  }
}
