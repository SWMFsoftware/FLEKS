#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
void Pic::assemble_ohm_E(const MultiFab& centerBin,
                         const MultiFab& centerBtimeAvg, MultiFab& Eout,
                         int iLev, amrex::Real hstep) {
  BL_PROFILE("Pic::assemble_ohm_E");

  const auto dx = Geom(iLev).CellSizeArray();
  const Real dxInv = 1.0 / (2.0 * dx[0]);
  const Real dyInv = 1.0 / (2.0 * dx[1]);
  const Real dzInv = (nDim > 2) ? 1.0 / (2.0 * dx[2]) : 0.0;

  // Cell-centred current J = curl(B)/(4*pi) from the trial B (2*dx central
  // difference, zero at the Nyquist wavenumber).
  const bool needJ =
      (etaResistivity > 0 || useHallTerm || etaHyperLev[iLev] > 0);
  if (needJ) {
    curl_center_to_center(centerBin, centerJ[iLev], Geom(iLev).InvCellSize());
  }

  for (MFIter mfi(Eout); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrE = Eout[mfi].array();
    const Array4<Real const>& arrB = centerBtimeAvg[mfi].array();
    const Array4<Real const>& moments =
        centerPlasmaSum[nSpecies][iLev][mfi].array();
    const Array4<Real const>& momentsPrev =
        centerPlasmaPrev[nSpecies][iLev][mfi].array();
    const Array4<Real const> arrJ =
        needJ ? centerJ[iLev][mfi].array() : Array4<Real const>();

    // Moment time-interpolation weights: X(hstep) =
    // (0.5-hstep)*X^{n-1/2} + (0.5+hstep)*X^{n+1/2}.
    const Real wPrev = 0.5 - hstep;
    const Real wCur = 0.5 + hstep;

    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real rhoPrev = momentsPrev(i, j, k, iRho_);
      const Real rhoCur = moments(i, j, k, iRho_);
      const Real rho = wPrev * rhoPrev + wCur * rhoCur;
      const Real mx =
          wPrev * momentsPrev(i, j, k, iUx_) + wCur * moments(i, j, k, iUx_);
      const Real my =
          wPrev * momentsPrev(i, j, k, iUy_) + wCur * moments(i, j, k, iUy_);
      const Real mz =
          wPrev * momentsPrev(i, j, k, iUz_) + wCur * moments(i, j, k, iUz_);
      Real ui = 0, vi = 0, wi = 0;

      if (rho > 0) {
        ui = mx / rho;
        vi = my / rho;
        wi = mz / rho;
      }

      // Interpolated density at an arbitrary cell (same hstep weights), used
      // for the electron-pressure gradient closure.
      auto rho_at = [=](int ii, int jj, int kk) AMREX_GPU_DEVICE {
        return wPrev * momentsPrev(ii, jj, kk, iRho_) +
               wCur * moments(ii, jj, kk, iRho_);
      };

      Real bx = arrB(i, j, k, ix_);
      Real by = arrB(i, j, k, iy_);
      Real bz = arrB(i, j, k, iz_);

      // Convection term: E = -U_i x B
      Real ex = -(vi * bz - wi * by);
      Real ey = -(wi * bx - ui * bz);
      Real ez = -(ui * by - vi * bx);

      // J = curl(B)/(4*pi) (CGS)
      Real jx = 0.0, jy = 0.0, jz = 0.0;
      if (needJ) {
        jx = arrJ(i, j, k, ix_) / fourPI;
        jy = arrJ(i, j, k, iy_) / fourPI;
        jz = arrJ(i, j, k, iz_) / fourPI;
      }

      // eta * J
      if (etaResistivity > 0) {
        ex += etaResistivity * jx;
        ey += etaResistivity * jy;
        ez += etaResistivity * jz;
      }

      // Electron-pressure-gradient and Hall terms. The floor caps 1/rho; the
      // pressure closure itself uses the true rho. Cells with rho == 0 are
      // left inert.
      if (rho > 0) {
        const Real invRhoEff = 1.0 / amrex::max(rho, rhoMinOhm);

        // Electron pressure gradient
        Real dPe_dx = 0.0, dPe_dy = 0.0, dPe_dz = 0.0;
        if (electronTemperature > 0) {
          if (electronGamma == 1.0) {
            // Isothermal: grad(Pe) = Te * grad(rho)
            dPe_dx = electronTemperature *
                     (rho_at(i + 1, j, k) - rho_at(i - 1, j, k)) * dxInv;
            dPe_dy = electronTemperature *
                     (rho_at(i, j + 1, k) - rho_at(i, j - 1, k)) * dyInv;
            dPe_dz = (nDim > 2)
                         ? electronTemperature *
                               (rho_at(i, j, k + 1) - rho_at(i, j, k - 1)) *
                               dzInv
                         : 0.0;
          } else {
            // Adiabatic: Pe = P0 * (rho / rho0)^gamma
            Real p0 = electronDensity0 * electronTemperature;
            Real invRho0 = 1.0 / electronDensity0;

            auto calc_Pe = [=](Real r) {
              return (r > 0) ? p0 * std::pow(r * invRho0, electronGamma) : 0.0;
            };

            dPe_dx =
                (calc_Pe(rho_at(i + 1, j, k)) - calc_Pe(rho_at(i - 1, j, k))) *
                dxInv;
            dPe_dy =
                (calc_Pe(rho_at(i, j + 1, k)) - calc_Pe(rho_at(i, j - 1, k))) *
                dyInv;
            dPe_dz = (nDim > 2) ? (calc_Pe(rho_at(i, j, k + 1)) -
                                   calc_Pe(rho_at(i, j, k - 1))) *
                                      dzInv
                                : 0.0;
          }

          ex -= dPe_dx * invRhoEff;
          ey -= dPe_dy * invRhoEff;
          ez -= dPe_dz * invRhoEff;
        }

        // Hall term: (J x B) / rho_q
        if (useHallTerm) {
          Real hall_x = (jy * bz - jz * by) * invRhoEff;
          Real hall_y = (jz * bx - jx * bz) * invRhoEff;
          Real hall_z = (jx * by - jy * bx) * invRhoEff;

          ex += hall_x;
          ey += hall_y;
          ez += hall_z;
        }
      }

      arrE(i, j, k, ix_) = ex;
      arrE(i, j, k, iy_) = ey;
      arrE(i, j, k, iz_) = ez;
    });
  }

  // Hyper-resistivity: E -= (eta_h / 4*pi) * curl(nabla^2 B) (see
  // tests/hyper_resistivity). Applied before the final BC pass so ghost cells
  // include the hyper-resistive term.
  if (etaHyperLev[iLev] > 0) {
    lap_center_to_center(centerBin, centerLapB[iLev], Geom(iLev).InvCellSize());
    centerLapB[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_field_bc(cellStatus[iLev], centerLapB[iLev], 0,
                   centerLapB[iLev].nComp(), &Pic::get_center_B, iLev, true);

    curl_center_to_center(centerLapB[iLev], centerHyperE[iLev],
                          Geom(iLev).InvCellSize());
    centerHyperE[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_field_bc(cellStatus[iLev], centerHyperE[iLev], 0,
                   centerHyperE[iLev].nComp(), &Pic::get_center_E, iLev, false);

    const Real f = etaHyperLev[iLev] / fourPI;
    MultiFab::Saxpy(Eout, -f, centerHyperE[iLev], 0, 0, 3, 0);
  }

  Eout.FillBoundary(Geom(iLev).periodicity());
  apply_field_bc(cellStatus[iLev], Eout, 0, nDim3, &Pic::get_center_E, iLev,
                 false);
}

//==========================================================
void Pic::smooth_moments() {
  std::string nameFunc = "Pic::smooth_moments";
  timing_func(nameFunc);

  if (!doSmoothMoments || nSmoothMoments <= 0)
    return;

  // Smooth the total ion moments the Ohm's law reads. Hybrid-only.
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab& moments = useHybridPIC ? centerPlasmaSum[nSpecies][iLev]
                                     : nodePlasma[nSpecies][iLev];
    moments.FillBoundary(Geom(iLev).periodicity());
    for (int icount = 0; icount < nSmoothMoments; ++icount) {
      smooth_multifab(moments, iLev, 1, coefSmoothMoments);
    }
    moments.FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
// Copy the current summed moment deposit into centerPlasmaPrev (J^{n-1/2})
// before a fresh deposit (J^{n+1/2}), so assemble_ohm_E can time-interpolate
// the two at the magnetic sub-step fraction hstep.
void Pic::save_current_moments_to_prev() {
  // Hybrid-only. Copy the rho + 3 momentum components of the current summed
  // deposit into centerPlasmaPrev so the Ohm's law can time-interpolate the
  // previous and current moments. Re-fill ghosts so grad(Pe) at box boundaries
  // reads valid values.
  std::string nameFunc = "Pic::save_current_moments_to_prev";
  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::Copy(centerPlasmaPrev[nSpecies][iLev],
                   centerPlasmaSum[nSpecies][iLev], 0, 0, nHybridMomentsComps,
                   centerPlasmaSum[nSpecies][iLev].nGrow());
    centerPlasmaPrev[nSpecies][iLev].FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
// Seed centerPlasmaPrev on the first hybrid step, where there is no previous
// deposit: initialise it from the current deposit so the time interpolation
// degrades to a plain average for that single step. Hybrid-only.
void Pic::seed_first_hybrid_step() {
  std::string nameFunc = "Pic::seed_first_hybrid_step";
  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::Copy(centerPlasmaPrev[nSpecies][iLev],
                   centerPlasmaSum[nSpecies][iLev], 0, 0, nHybridMomentsComps,
                   centerPlasmaSum[nSpecies][iLev].nGrow());
    centerPlasmaPrev[nSpecies][iLev].FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
void Pic::project_centerB_to_nodeB(int iLev) {
  centerB[iLev].FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
                   &Pic::get_center_B, iLev, true);
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
                   &Pic::get_node_B, iLev, true);
  } else {
    fill_fine_lev_bny_from_coarse(nodeB[iLev - 1], nodeB[iLev], 0,
                                  nodeB[iLev - 1].nComp(), ref_ratio[iLev - 1],
                                  Geom(iLev - 1), Geom(iLev), node_status(iLev),
                                  node_bilinear_interp);
  }
}

// BCs for the cell-centred B (the cell-centred part of
// project_centerB_to_nodeB), called at the end of each sub-step.
void Pic::apply_centerB_BC(int iLev) { apply_centerB_BC(iLev, centerB[iLev]); }

void Pic::apply_centerB_BC(int iLev, amrex::MultiFab& mfB) {
  mfB.FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(cellStatus[iLev], mfB, 0, mfB.nComp(), &Pic::get_center_B,
                   iLev, true);
  } else {
    fill_fine_lev_bny_from_coarse(
        centerB[iLev - 1], mfB, 0, mfB.nComp(), ref_ratio[iLev - 1],
        Geom(iLev - 1), Geom(iLev), cell_status(iLev), *get_cell_interp());
  }
}

//==========================================================
void Pic::project_centerB_to_nodeB_scratch(amrex::MultiFab& centerIn,
                                           amrex::MultiFab& nodeOut, int iLev) {
  // Same projection as project_centerB_to_nodeB on caller-owned scratch fields.
  centerIn.FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(cellStatus[iLev], centerIn, 0, centerIn.nComp(),
                   &Pic::get_center_B, iLev, true);
  } else {
    fill_fine_lev_bny_from_coarse(
        centerB[iLev - 1], centerIn, 0, centerIn.nComp(), ref_ratio[iLev - 1],
        Geom(iLev - 1), Geom(iLev), cell_status(iLev), *get_cell_interp());
  }
  average_center_to_node(centerIn, nodeOut);
  nodeOut.FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(nodeStatus[iLev], nodeOut, 0, nodeOut.nComp(),
                   &Pic::get_node_B, iLev, true);
  } else {
    fill_fine_lev_bny_from_coarse(
        nodeB[iLev - 1], nodeOut, 0, nodeOut.nComp(), ref_ratio[iLev - 1],
        Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
  }
}

//==========================================================
void Pic::update_B_hybrid() {
  std::string nameFunc = "Pic::update_B_hybrid";
  timing_func(nameFunc);

  const Real dt = tc->get_dt();
  const Real subDt = dt / nBSubcycle;

  // Grid-mode hyper-resistivity: uniform eta_h based on finest dx to keep
  // diffusion stable across levels.
  if (etaHyperMode == "grid" && etaHyperCh > 0) {
    const int iFinest = n_lev() - 1;
    const auto dxFine = Geom(iFinest).CellSizeArray();
    Real dxMinFine = dxFine[0];
    for (int d = 1; d < nDim; ++d)
      dxMinFine = amrex::min(dxMinFine, dxFine[d]);

    const Real etaHyper = fourPI * etaHyperCh * std::pow(dxMinFine, 4) / dt;
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      etaHyperLev[iLev] = etaHyper;
    }
  }

  // CFL stability check for explicit resistive and hyper-resistive diffusion.
  const Real cflLimit = useRK4 ? 2.785 : 2.513;
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    if (etaResistivity <= 0 && etaHyperLev[iLev] <= 0)
      continue;

    const auto dx = Geom(iLev).CellSizeArray();
    const Box& domBox = Geom(iLev).Domain();
    Real sMax = 0.0, lMax = 0.0;
    for (int iDim = 0; iDim < nDim; ++iDim) {
      const int nCellDim = domBox.length(iDim);
      if (nCellDim < 2)
        continue;
      const Real invDx2 = 1.0 / (dx[iDim] * dx[iDim]);
      lMax += 4.0 * invDx2;
      if (nCellDim >= 3)
        sMax += invDx2;
    }
    const Real cflEta = (etaResistivity / fourPI) * sMax * subDt;
    const Real cflHyper = (etaHyperLev[iLev] / fourPI) * sMax * lMax * subDt;
    if (cflEta > cflLimit)
      amrex::Print()
          << "  [CFL warning] resistivity: eta*kmax^2*dt_sub/(4pi) = " << cflEta
          << " (> " << cflLimit << ", explicit diffusion may be unstable)\n";
    if (cflHyper > cflLimit)
      amrex::Print()
          << "  [CFL warning] hyper-resistivity: eta_h*kmax^4*dt_sub/(4pi) = "
          << cflHyper << " (> " << cflLimit
          << ", explicit 4th-order diffusion may be unstable)\n";
  }

  const Real invSubcycle = 1.0 / static_cast<Real>(nBSubcycle);

  for (int subStep = 0; subStep < nBSubcycle; ++subStep) {
    // Moment time-interpolation weights hstep for RK stages within the
    // sub-step.
    const Real g = static_cast<Real>(subStep) * invSubcycle;
    const Real hstepStart = g;
    const Real hstepHalf = g + 0.5 * invSubcycle;
    const Real hstepEnd = g + invSubcycle;

    if (useRK4) {
      // Classical RK4 on dB/dt = -curl(E_Ohm(B)) across all levels.
      const Real dtSixth = -subDt / 6.0;
      const Real dtThird = -subDt / 3.0;

      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        // Stage 1: k1 = curl(E(B^n))
        assemble_ohm_E(centerB[iLev], centerB[iLev], centerEstage[iLev], iLev,
                       hstepStart);
        curl_center_to_center(centerEstage[iLev], kStage[iLev][0],
                              Geom(iLev).InvCellSize());

        // Stage 2: B2 = B^n - 0.5 dt k1; evaluate E at (B2 + B^n)/2
        MultiFab::Copy(centerBstage[iLev], centerB[iLev], 0, 0, 3, nGst);
        MultiFab::Saxpy(centerBstage[iLev], -0.5 * subDt, kStage[iLev][0], 0, 0,
                        3, nGst);
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerB[iLev], 0, 0, 3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstage[iLev], centerBstar[iLev],
                       centerEstage[iLev], iLev, hstepHalf);
        curl_center_to_center(centerEstage[iLev], kStage[iLev][1],
                              Geom(iLev).InvCellSize());

        // Stage 3: B3 = B^n - 0.5 dt k2; evaluate E at (B3 + B^n)/2
        MultiFab::Copy(centerBstage[iLev], centerB[iLev], 0, 0, 3, nGst);
        MultiFab::Saxpy(centerBstage[iLev], -0.5 * subDt, kStage[iLev][1], 0, 0,
                        3, nGst);
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerB[iLev], 0, 0, 3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstage[iLev], centerBstar[iLev],
                       centerEstage[iLev], iLev, hstepHalf);
        curl_center_to_center(centerEstage[iLev], kStage[iLev][2],
                              Geom(iLev).InvCellSize());

        // Stage 4: B4 = B^n - dt k3; evaluate E at (B4 + B^n)/2
        MultiFab::Copy(centerBstage[iLev], centerB[iLev], 0, 0, 3, nGst);
        MultiFab::Saxpy(centerBstage[iLev], -subDt, kStage[iLev][2], 0, 0, 3,
                        nGst);
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerB[iLev], 0, 0, 3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstage[iLev], centerBstar[iLev],
                       centerEstage[iLev], iLev, hstepEnd);
        curl_center_to_center(centerEstage[iLev], kStage[iLev][3],
                              Geom(iLev).InvCellSize());

        // Accumulate RK4: B^{n+1} = B^n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
        MultiFab::Saxpy(centerB[iLev], dtSixth, kStage[iLev][0], 0, 0, 3, nGst);
        MultiFab::Saxpy(centerB[iLev], dtThird, kStage[iLev][1], 0, 0, 3, nGst);
        MultiFab::Saxpy(centerB[iLev], dtThird, kStage[iLev][2], 0, 0, 3, nGst);
        MultiFab::Saxpy(centerB[iLev], dtSixth, kStage[iLev][3], 0, 0, 3, nGst);
        centerB[iLev].FillBoundary(Geom(iLev).periodicity());

        apply_centerB_BC(iLev);
      }
      continue;
    }

    if (fieldIntegrator == "ssprk3") {
      // Strong-stability-preserving RK3 with time-centred E evaluation.
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        MultiFab::Copy(centerBstart[iLev], centerB[iLev], 0, 0, 3, nGst);

        // Stage 1: B1 = B_n - subDt * curl(E(B_n))
        assemble_ohm_E(centerB[iLev], centerB[iLev], centerEstage[iLev], iLev,
                       hstepStart);
        curl_center_to_center(centerEstage[iLev], kStage[iLev][0],
                              Geom(iLev).InvCellSize());
        MultiFab::Copy(centerBstage[iLev], centerB[iLev], 0, 0, 3, nGst);
        MultiFab::Saxpy(centerBstage[iLev], -subDt, kStage[iLev][0], 0, 0, 3,
                        nGst);

        // Stage 2: B2 = (3/4)*B_n + (1/4)*(B1 - subDt * curl(E(avgB2)))
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerBstart[iLev], 0, 0, 3, nGst);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstar[iLev], centerBstar[iLev], centerEstage[iLev],
                       iLev, hstepEnd);
        curl_center_to_center(centerEstage[iLev], kStage[iLev][1],
                              Geom(iLev).InvCellSize());
        MultiFab::LinComb(centerBstage[iLev], 0.25, centerBstage[iLev], 0, 0.75,
                          centerBstart[iLev], 0, 0, 3, nGst);
        MultiFab::Saxpy(centerBstage[iLev], -0.25 * subDt, kStage[iLev][1], 0,
                        0, 3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);

        // Stage 3: B^{n+1} = (1/3)*B_n + (2/3)*(B2 - subDt * curl(E(avgB3)))
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerBstart[iLev], 0, 0, 3, nGst);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstar[iLev], centerBstar[iLev], centerEstage[iLev],
                       iLev, hstepHalf);
        curl_center_to_center(centerEstage[iLev], kStage[iLev][2],
                              Geom(iLev).InvCellSize());
        MultiFab::LinComb(centerB[iLev], 2.0 / 3.0, centerBstage[iLev], 0,
                          1.0 / 3.0, centerBstart[iLev], 0, 0, 3, nGst);
        MultiFab::Saxpy(centerB[iLev], (-2.0 / 3.0) * subDt, kStage[iLev][2], 0,
                        0, 3, nGst);
        centerB[iLev].FillBoundary(Geom(iLev).periodicity());

        apply_centerB_BC(iLev);
      }
      continue;
    }
  }

  if (projectDownEmFields && finest_level > 0) {
    for (int iLev = finest_level; iLev > 0; iLev--) {
      average_down(centerB[iLev], centerB[iLev - 1], 0, 3, ref_ratio[0]);
    }
  }

  auto& cellInterp = *get_cell_interp();
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
                     &Pic::get_center_B, iLev, true);
    } else {
      fill_fine_lev_bny_from_coarse(centerB[iLev - 1], centerB[iLev], 0,
                                    centerB[iLev - 1].nComp(),
                                    ref_ratio[iLev - 1], Geom(iLev - 1),
                                    Geom(iLev), cell_status(iLev), cellInterp);
    }
  }

  // Running time-averaged B used in Ohm's law and the particle push.
  if (useAvgFieldB) {
    const Real alpha = (nAvgFieldB > 1) ? (1.0 - 1.0 / nAvgFieldB) : 0.0;
    const bool syncNodeBavg = !useHybridPIC;
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (!isBavgInit) {
        MultiFab::Copy(centerBavg[iLev], centerB[iLev], 0, 0, 3,
                       centerBavg[iLev].nGrow());
        if (syncNodeBavg) {
          MultiFab::Copy(nodeBavg[iLev], nodeB[iLev], 0, 0, 3,
                         nodeBavg[iLev].nGrow());
        }
        isBavgInit = true;
      } else {
        centerBavg[iLev].mult(alpha);
        MultiFab::Saxpy(centerBavg[iLev], 1.0 - alpha, centerB[iLev], 0, 0, 3,
                        centerBavg[iLev].nGrow());
        if (syncNodeBavg) {
          nodeBavg[iLev].mult(alpha);
          MultiFab::Saxpy(nodeBavg[iLev], 1.0 - alpha, nodeB[iLev], 0, 0, 3,
                          nodeBavg[iLev].nGrow());
        }
      }
      centerBavg[iLev].FillBoundary(Geom(iLev).periodicity());
      if (syncNodeBavg) {
        nodeBavg[iLev].FillBoundary(Geom(iLev).periodicity());
      }
      if (iLev == 0) {
        apply_field_bc(cellStatus[iLev], centerBavg[iLev], 0, 3,
                       &Pic::get_center_B, iLev, true);
      }
    }
  }

  // Fill coarse-fine interface ghost cells for centerBavg.
  if (useAvgFieldB && isBavgInit && finest_level > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(centerBavg[iLev - 1], centerBavg[iLev], 0,
                                    3, ref_ratio[iLev - 1], Geom(iLev - 1),
                                    Geom(iLev), cell_status(iLev), cellInterp);
    }
  }

  // Evaluate integer-step E^{n+1} into centerEhybrid for the next particle
  // push.
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    bool hasParticles = false;
    for (int i : kineticSpecies_) {
      if (parts[i]->NumberOfParticlesAtLevel(iLev, true, true) > 0) {
        hasParticles = true;
        break;
      }
    }
    if (!hasParticles) {
      continue;
    }
    const auto& cBin =
        (useAvgFieldB && isBavgInit) ? centerBavg[iLev] : centerB[iLev];
    assemble_ohm_E(cBin, cBin, centerEhybrid[iLev], iLev, 1.0);
  }

  // BCs for cell-centred E (read by the hybrid Boris push).
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    centerEhybrid[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_field_bc(cellStatus[iLev], centerEhybrid[iLev], 0,
                   centerEhybrid[iLev].nComp(), &Pic::get_center_E, iLev,
                   false);
  }

  // Fill coarse-fine interface ghost cells for centerEhybrid.
  if (finest_level > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(centerEhybrid[iLev - 1],
                                    centerEhybrid[iLev], 0, 3,
                                    ref_ratio[iLev - 1], Geom(iLev - 1),
                                    Geom(iLev), cell_status(iLev), cellInterp);
    }
  }

  // Suppress grid-scale (odd-even) E component.
  if (doSmoothE) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      centerEhybrid[iLev].FillBoundary(Geom(iLev).periodicity());
      smooth_E(centerEhybrid[iLev], iLev);
    }
  }
}

//==========================================================
