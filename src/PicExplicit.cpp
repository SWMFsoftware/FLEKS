#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
// Push the kinetic species and deposit the charge-conserving (Esirkepov)
// current of the push, J = sum_s q_s dx/dt, into J. trialOnly: deposit the
// current of the push without moving the particles.
void Pic::push_and_deposit_explicit_current(int pushStage, Vector<MultiFab>& J,
                                            bool trialOnly) {
  std::string nameFunc = "Pic::push_and_deposit_explicit_current";
  timing_func(nameFunc);

  const Real dt = tc->get_dt();
  const Real dtNext = tc->get_next_dt();

  for (int iLev = 0; iLev < n_lev(); ++iLev)
    J[iLev].setVal(0.0);

  for (int i : kineticSpecies_) {
    parts[i]->mover(nodeE, nodeB, nodeB0, eBg, uBg, dt, dtNext, pushStage,
                    &nodeJ, trialOnly);
    for (int iLev = 0; iLev < n_lev(); ++iLev)
      MultiFab::Add(J[iLev], nodeJ[iLev], 0, 0, nDim3, get_n_ghost());
  }

  if (!trialOnly) {
    for (int i : kineticSpecies_)
      parts[i]->redistribute_particles();
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    // Same correction of the diluted surface nodes as for the moments.
    if (use_body())
      rescale_body_surface_nodes(J[iLev], cell_status(iLev), node_status(iLev),
                                 Geom(iLev).Domain(), nDim);
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
// Subcycled explicit Maxwell advance for B and E with the current J of the
// step (time-centered at t^{n+1/2}) in every sub-step.
void Pic::update_EB_explicit_leapfrog(Real dt, int nSub,
                                      const Vector<MultiFab>& J) {
  std::string nameFunc = "Pic::update_EB_explicit_leapfrog";
  timing_func(nameFunc);

  const Real subDt = dt / nSub;
  const Real t0 = tc ? tc->get_time() : 0.0;

  for (int iSub = 0; iSub < nSub; ++iSub) {
    // 1. Advance B half-step: B^m -> B^{m+1/2}
    const Real tB_half = t0 + (iSub + 0.5) * subDt;
    advance_B_half_step(subDt, tB_half);

    // 2. Advance E full-step: E^m -> E^{m+1} using B^{m+1/2} and J.
    const Real tE_full = t0 + (iSub + 1.0) * subDt;
    for (int iLev = 0; iLev < n_lev(); ++iLev)
      advance_E_step(subDt, J[iLev], iLev, tE_full);

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

  // 1. Backup base fields B^n and E^n
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::Copy(centerBstar[iLev], centerB[iLev], 0, 0, nDim3,
                   get_n_ghost());
    MultiFab::Copy(nodeEth[iLev], nodeE[iLev], 0, 0, nDim3, get_n_ghost());
  }

  // 2. Predictor: trial push with (E^n, B^n) that only deposits its current
  // J* = q (x* - x^n)/dt; the particles stay at t^n.
  push_and_deposit_explicit_current(PicParticles::kPushPredictor, nodeJexpl,
                                    true);

  // 3. Advance trial fields to t^(n+1) with J*.
  update_EB_explicit_leapfrog(dt, nEMSubcycle, nodeJexpl);

  // 4. Construct time-centered fields at t^(n+1/2)
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

  // 5. Corrector push from t^n with the time-centered fields, depositing the
  // current J^(n+1/2) of the accepted displacement.
  push_and_deposit_explicit_current(PicParticles::kPushCorrector, nodeJexpl);
  finish_particle_push();

  // 6. Restore B^n, E^n and advance to (B^{n+1}, E^{n+1}) with J^(n+1/2).
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::Copy(centerB[iLev], centerBstar[iLev], 0, 0, nDim3,
                   get_n_ghost());
    MultiFab::Copy(nodeE[iLev], nodeEth[iLev], 0, 0, nDim3, get_n_ghost());
  }
  update_EB_explicit_leapfrog(dt, nSubcycle, nodeJexpl);
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
    // 1. Kick v^{n-1/2} -> v^{n+1/2} with (E^n, B^n) and drift to x^{n+1},
    // depositing J^{n+1/2} = q (x^{n+1} - x^n)/dt in the same pass.
    push_and_deposit_explicit_current(PicParticles::kPushLeapfrog, nodeJexpl);
    finish_particle_push();

    // 2. B half-step, E full step with J^{n+1/2}, B half-step (subcycled).
    update_EB_explicit_leapfrog(dt, nSubcycle, nodeJexpl);
  }

  // A div(E) correction, when on, follows in Pic::update(), as for the
  // semi-implicit solver.
}
