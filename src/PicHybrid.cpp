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
                         int iLev, Real hstep) {
  BL_PROFILE("Pic::assemble_ohm_E");

  // Nodal total current J = curl(B)/(4*pi) from trial B (compact 1*dx stencil
  // from cell centres to nodes). Only needed for physical resistivity and Hall.
  const bool needJ = (etaResistivity > 0 || useHallTerm);
  if (needJ) {
    curl_center_to_node(centerBin, nodeJ[iLev], Geom(iLev).InvCellSize());
    nodeJ[iLev].FillBoundary(Geom(iLev).periodicity());
    average_node_to_center(nodeJ[iLev], centerJ[iLev]);
    centerJ[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  // Magnetic field interpolated from cell centres to nodes for vector cross products.
  average_center_to_node(centerBtimeAvg, nodeBstage[iLev]);
  nodeBstage[iLev].FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(nodeStatus[iLev], nodeBstage[iLev], 0, 3, &Pic::get_node_B,
                   iLev, true);
  }

  // Moment time-interpolation weights: X(hstep) =
  // (0.5-hstep)*X^{n-1/2} + (0.5+hstep)*X^{n+1/2}.
  const Real wPrev = 0.5 - hstep;
  const Real wCur = 0.5 + hstep;
  const Real invFourPI = 1.0 / fourPI;

  // Electron pressure gradient at nodes:
  // If electron temperature > 0, interpolate density in time at nodes, average to
  // cell centres, evaluate EOS Pe at cell centres, and take grad_center_to_node(Pe).
  if (electronTemperature > 0) {
    for (MFIter mfi(nodeRhoTemp[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrRho = nodeRhoTemp[iLev][mfi].array();
      const Array4<Real const>& moments =
          nodePlasma[nSpecies][iLev][mfi].array();
      const Array4<Real const>& momentsPrev =
          nodePlasmaPrev[nSpecies][iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        arrRho(i, j, k) = wPrev * momentsPrev(i, j, k, iRho_) +
                          wCur * moments(i, j, k, iRho_);
      });
    }
    nodeRhoTemp[iLev].FillBoundary(Geom(iLev).periodicity());

    average_node_to_center(nodeRhoTemp[iLev], centerPe[iLev]);
    centerPe[iLev].FillBoundary(Geom(iLev).periodicity());

    const Real p0 = electronDensity0 * electronTemperature;
    const Real invRho0 =
        (electronDensity0 > 0.0) ? (1.0 / electronDensity0) : 0.0;
    const Real gamma = electronGamma;
    const Real Te = electronTemperature;

    for (MFIter mfi(centerPe[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& arrPe = centerPe[iLev][mfi].array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        Real r = arrPe(i, j, k);
        if (gamma == 1.0) {
          arrPe(i, j, k) = Te * r;
        } else {
          arrPe(i, j, k) =
              (r > 0) ? p0 * std::pow(r * invRho0, gamma) : 0.0;
        }
      });
    }
    centerPe[iLev].FillBoundary(Geom(iLev).periodicity());

    // Zero-gradient (Neumann) BC across non-periodic domain boundaries.
    if (!Geom(iLev).isAllPeriodic() && centerPe[iLev].nGrow() > 0) {
      const Box& dom = Geom(iLev).Domain();
      for (MFIter mfi(centerPe[iLev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& arr = centerPe[iLev][mfi].array();
        for (int iDim = 0; iDim < nDim; ++iDim) {
          if (Geom(iLev).isPeriodic(iDim))
            continue;
          if (bx.smallEnd(iDim) == dom.smallEnd(iDim)) {
            IntVect lo = bx.smallEnd();
            IntVect hi = bx.bigEnd();
            lo[iDim] = dom.smallEnd(iDim) - 1;
            hi[iDim] = dom.smallEnd(iDim) - 1;
            Box ghostBox(lo, hi);
            ParallelFor(ghostBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
              IntVect src{ AMREX_D_DECL(i, j, k) };
              src[iDim] = dom.smallEnd(iDim);
              arr(i, j, k) = arr(src);
            });
          }
          if (bx.bigEnd(iDim) == dom.bigEnd(iDim)) {
            IntVect lo = bx.smallEnd();
            IntVect hi = bx.bigEnd();
            lo[iDim] = dom.bigEnd(iDim) + 1;
            hi[iDim] = dom.bigEnd(iDim) + 1;
            Box ghostBox(lo, hi);
            ParallelFor(ghostBox, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
              IntVect src{ AMREX_D_DECL(i, j, k) };
              src[iDim] = dom.bigEnd(iDim);
              arr(i, j, k) = arr(src);
            });
          }
        }
      }
    }

    grad_center_to_node(centerPe[iLev], nodeGradPe[iLev],
                        Geom(iLev).InvCellSize());
    nodeGradPe[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  for (MFIter mfi(Eout); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrE = Eout[mfi].array();
    const Array4<Real const>& arrB = nodeBstage[iLev][mfi].array();
    const Array4<Real const>& moments =
        nodePlasma[nSpecies][iLev][mfi].array();
    const Array4<Real const>& momentsPrev =
        nodePlasmaPrev[nSpecies][iLev][mfi].array();
    const Array4<Real const> arrJ =
        needJ ? nodeJ[iLev][mfi].array() : Array4<Real const>();
    const Array4<Real const> arrGradPe =
        (electronTemperature > 0) ? nodeGradPe[iLev][mfi].array()
                                  : Array4<Real const>();

    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real rhoPrev = momentsPrev(i, j, k, iRho_);
      const Real rhoCur = moments(i, j, k, iRho_);
      const Real rho = wPrev * rhoPrev + wCur * rhoCur;
      const Real mx =
          wPrev * momentsPrev(i, j, k, iMx_) + wCur * moments(i, j, k, iMx_);
      const Real my =
          wPrev * momentsPrev(i, j, k, iMy_) + wCur * moments(i, j, k, iMy_);
      const Real mz =
          wPrev * momentsPrev(i, j, k, iMz_) + wCur * moments(i, j, k, iMz_);
      Real ui = 0, vi = 0, wi = 0;

      if (rho > 0) {
        ui = mx / rho;
        vi = my / rho;
        wi = mz / rho;
      }

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
        jx = arrJ(i, j, k, ix_) * invFourPI;
        jy = arrJ(i, j, k, iy_) * invFourPI;
        jz = arrJ(i, j, k, iz_) * invFourPI;
      }

      // eta * J
      if (etaResistivity > 0) {
        ex += etaResistivity * jx;
        ey += etaResistivity * jy;
        ez += etaResistivity * jz;
      }

      // Electron-pressure-gradient and Hall terms. The floor caps 1/rho.
      if (rho > 0) {
        const Real invRhoEff = 1.0 / amrex::max(rho, rhoMinOhm);

        // Electron pressure gradient
        if (electronTemperature > 0) {
          ex -= arrGradPe(i, j, k, ix_) * invRhoEff;
          ey -= arrGradPe(i, j, k, iy_) * invRhoEff;
          ez -= arrGradPe(i, j, k, iz_) * invRhoEff;
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

  // Hyper-resistivity: E -= (eta_h / 4*pi) * curl(nabla^2 B).
  // centerLapB = Laplacian(centerBin); nodeHyperE = curl_center_to_node(centerLapB).
  if (etaHyperLev[iLev] > 0) {
    lap_center_to_center(centerBin, centerLapB[iLev], Geom(iLev).InvCellSize());
    centerLapB[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_field_bc(cellStatus[iLev], centerLapB[iLev], 0,
                   centerLapB[iLev].nComp(), &Pic::get_center_B, iLev, true);

    curl_center_to_node(centerLapB[iLev], nodeHyperE[iLev],
                        Geom(iLev).InvCellSize());
    nodeHyperE[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_field_bc(nodeStatus[iLev], nodeHyperE[iLev], 0,
                   nodeHyperE[iLev].nComp(), &Pic::get_node_E, iLev, false);

    const Real f = etaHyperLev[iLev] / fourPI;
    MultiFab::Saxpy(Eout, -f, nodeHyperE[iLev], 0, 0, nDim3, 0);
  }

  Eout.FillBoundary(Geom(iLev).periodicity());
  apply_field_bc(nodeStatus[iLev], Eout, 0, nDim3, &Pic::get_node_E, iLev,
                 false);
}

//==========================================================
void Pic::smooth_moments() {
  std::string nameFunc = "Pic::smooth_moments";
  timing_func(nameFunc);

  if (!doSmoothMoments || nSmoothMoments <= 0)
    return;

  // Smooth the total ion moments on the node grid.
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab& moments = nodePlasma[nSpecies][iLev];
    moments.FillBoundary(Geom(iLev).periodicity());
    for (int icount = 0; icount < nSmoothMoments; ++icount) {
      smooth_multifab(moments, iLev, 1, coefSmoothMoments);
    }
    moments.FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
// Copy the current summed moment deposit into nodePlasmaPrev (J^{n-1/2})
// before a fresh deposit (J^{n+1/2}), so assemble_ohm_E can time-interpolate
// the two at the magnetic sub-step fraction hstep.
void Pic::save_current_moments_to_prev() {
  std::string nameFunc = "Pic::save_current_moments_to_prev";
  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    MultiFab::Copy(nodePlasmaPrev[nSpecies][iLev],
                   nodePlasma[nSpecies][iLev], 0, 0, nHybridMomentsComps,
                   nodePlasma[nSpecies][iLev].nGrow());
    nodePlasmaPrev[nSpecies][iLev].FillBoundary(Geom(iLev).periodicity());
  }
}

//==========================================================
// Seed nodePlasmaPrev on the first hybrid step, where there is no previous
// deposit: initialise it from the current deposit so the time interpolation
// degrades to a plain average for that single step.
void Pic::seed_first_hybrid_step() {
  std::string nameFunc = "Pic::seed_first_hybrid_step";
  timing_func(nameFunc);

  save_current_moments_to_prev();
}

//==========================================================
void Pic::project_centerB_to_nodeB(int iLev) {
  project_centerB_to_nodeB_scratch(centerB[iLev], nodeB[iLev], iLev);
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
  apply_centerB_BC(iLev, centerIn);
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
        assemble_ohm_E(centerB[iLev], centerB[iLev], nodeEstage[iLev], iLev,
                       hstepStart);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][0],
                            Geom(iLev).InvCellSize());

        // Stage 2: B2 = B^n - 0.5 dt k1; evaluate E at (B2 + B^n)/2
        MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0,
                          -0.5 * subDt, kStage[iLev][0], 0, 0, nDim3, nGst);
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerB[iLev], 0, 0, nDim3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstage[iLev], centerBstar[iLev],
                       nodeEstage[iLev], iLev, hstepHalf);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][1],
                            Geom(iLev).InvCellSize());

        // Stage 3: B3 = B^n - 0.5 dt k2; evaluate E at (B3 + B^n)/2
        MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0,
                          -0.5 * subDt, kStage[iLev][1], 0, 0, nDim3, nGst);
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerB[iLev], 0, 0, nDim3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstage[iLev], centerBstar[iLev],
                       nodeEstage[iLev], iLev, hstepHalf);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][2],
                            Geom(iLev).InvCellSize());

        // Stage 4: B4 = B^n - dt k3; evaluate E at (B4 + B^n)/2
        MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0, -subDt,
                          kStage[iLev][2], 0, 0, nDim3, nGst);
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerB[iLev], 0, 0, nDim3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstage[iLev], centerBstar[iLev],
                       nodeEstage[iLev], iLev, hstepEnd);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][3],
                            Geom(iLev).InvCellSize());

        // Accumulate RK4: B^{n+1} = B^n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
        MultiFab::Saxpy(centerB[iLev], dtSixth, kStage[iLev][0], 0, 0, nDim3,
                        nGst);
        MultiFab::Saxpy(centerB[iLev], dtThird, kStage[iLev][1], 0, 0, nDim3,
                        nGst);
        MultiFab::Saxpy(centerB[iLev], dtThird, kStage[iLev][2], 0, 0, nDim3,
                        nGst);
        MultiFab::Saxpy(centerB[iLev], dtSixth, kStage[iLev][3], 0, 0, nDim3,
                        nGst);

        apply_centerB_BC(iLev);
      }
      continue;
    }

    if (fieldIntegrator == "ssprk3") {
      // Strong-stability-preserving RK3 with time-centred E evaluation.
      for (int iLev = 0; iLev < n_lev(); ++iLev) {
        MultiFab::Copy(centerBstart[iLev], centerB[iLev], 0, 0, nDim3, nGst);

        // Stage 1: B1 = B_n - subDt * curl(E(B_n))
        assemble_ohm_E(centerB[iLev], centerB[iLev], nodeEstage[iLev], iLev,
                       hstepStart);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][0],
                            Geom(iLev).InvCellSize());
        MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0, -subDt,
                          kStage[iLev][0], 0, 0, nDim3, nGst);

        // Stage 2: B2 = (3/4)*B_n + (1/4)*(B1 - subDt * curl(E(avgB2)))
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerBstart[iLev], 0, 0, nDim3, nGst);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstar[iLev], centerBstar[iLev], nodeEstage[iLev],
                       iLev, hstepEnd);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][1],
                            Geom(iLev).InvCellSize());
        MultiFab::LinComb(centerBstage[iLev], 0.25, centerBstage[iLev], 0, 0.75,
                          centerBstart[iLev], 0, 0, nDim3, nGst);
        MultiFab::Saxpy(centerBstage[iLev], -0.25 * subDt, kStage[iLev][1], 0,
                        0, nDim3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);

        // Stage 3: B^{n+1} = (1/3)*B_n + (2/3)*(B2 - subDt * curl(E(avgB3)))
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerBstart[iLev], 0, 0, nDim3, nGst);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstar[iLev], centerBstar[iLev], nodeEstage[iLev],
                       iLev, hstepHalf);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][2],
                            Geom(iLev).InvCellSize());
        MultiFab::LinComb(centerB[iLev], 2.0 / 3.0, centerBstage[iLev], 0,
                          1.0 / 3.0, centerBstart[iLev], 0, 0, nDim3, nGst);
        MultiFab::Saxpy(centerB[iLev], (-2.0 / 3.0) * subDt, kStage[iLev][2], 0,
                        0, nDim3, nGst);

        apply_centerB_BC(iLev);
      }
      continue;
    }
  }

  if (projectDownEmFields && finest_level > 0) {
    for (int iLev = finest_level; iLev > 0; iLev--) {
      average_down(centerB[iLev], centerB[iLev - 1], 0, nDim3, ref_ratio[0]);
    }
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    apply_centerB_BC(iLev);
  }

  // Update nodeB from centerB on all levels.
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    average_center_to_node(centerB[iLev], nodeB[iLev]);
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(nodeStatus[iLev], nodeB[iLev], 0, nDim3, &Pic::get_node_B,
                     iLev, true);
    }
  }

  // Running time-averaged B used in Ohm's law and the particle push.
  if (useAvgFieldB) {
    const Real alpha = (nAvgFieldB > 1) ? (1.0 - 1.0 / nAvgFieldB) : 0.0;
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (!isBavgInit) {
        MultiFab::Copy(centerBavg[iLev], centerB[iLev], 0, 0, nDim3,
                       centerBavg[iLev].nGrow());
        MultiFab::Copy(nodeBavg[iLev], nodeB[iLev], 0, 0, nDim3,
                       nodeBavg[iLev].nGrow());
        isBavgInit = true;
      } else {
        centerBavg[iLev].mult(alpha);
        MultiFab::Saxpy(centerBavg[iLev], 1.0 - alpha, centerB[iLev], 0, 0,
                        nDim3, centerBavg[iLev].nGrow());
        average_center_to_node(centerBavg[iLev], nodeBavg[iLev]);
      }
      centerBavg[iLev].FillBoundary(Geom(iLev).periodicity());
      nodeBavg[iLev].FillBoundary(Geom(iLev).periodicity());
      if (iLev == 0) {
        apply_field_bc(cellStatus[iLev], centerBavg[iLev], 0, nDim3,
                       &Pic::get_center_B, iLev, true);
        apply_field_bc(nodeStatus[iLev], nodeBavg[iLev], 0, nDim3,
                       &Pic::get_node_B, iLev, true);
      }
    }
  }

  // Fill coarse-fine interface ghost cells for centerBavg and nodeBavg.
  if (useAvgFieldB && isBavgInit && finest_level > 0) {
    auto& cellInterp = *get_cell_interp();
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(centerBavg[iLev - 1], centerBavg[iLev], 0,
                                    nDim3, ref_ratio[iLev - 1], Geom(iLev - 1),
                                    Geom(iLev), cell_status(iLev), cellInterp);
      fill_fine_lev_bny_from_coarse(nodeBavg[iLev - 1], nodeBavg[iLev], 0,
                                    nDim3, ref_ratio[iLev - 1], Geom(iLev - 1),
                                    Geom(iLev), node_status(iLev),
                                    node_bilinear_interp);
    }
  }

  // Fill coarse-fine interface ghost cells for nodeB.
  if (finest_level > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(nodeB[iLev - 1], nodeB[iLev], 0, nDim3,
                                    ref_ratio[iLev - 1], Geom(iLev - 1),
                                    Geom(iLev), node_status(iLev),
                                    node_bilinear_interp);
    }
  }

  // Evaluate E^{n+1} into nodeE for the next push.
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const auto& cBin =
        (useAvgFieldB && isBavgInit) ? centerBavg[iLev] : centerB[iLev];
    assemble_ohm_E(cBin, cBin, nodeE[iLev], iLev, 1.0);
  }

  // Fill coarse-fine interface ghost cells for nodeE.
  if (finest_level > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(nodeE[iLev - 1], nodeE[iLev], 0, nDim3,
                                    ref_ratio[iLev - 1], Geom(iLev - 1),
                                    Geom(iLev), node_status(iLev),
                                    node_bilinear_interp);
    }
  }

  // Suppress grid-scale E component if enabled.
  if (doSmoothE) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
      smooth_E(nodeE[iLev], iLev);
    }
  }

  // Sync centerEhybrid for backwards compatibility with any remaining legacy readers.
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    average_node_to_center(nodeE[iLev], centerEhybrid[iLev]);
    centerEhybrid[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_field_bc(cellStatus[iLev], centerEhybrid[iLev], 0, nDim3,
                   &Pic::get_center_E, iLev, false);
  }
}
