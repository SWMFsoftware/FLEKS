#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_PhysBCFunct.H>

#include "GridUtility.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
void Pic::assemble_ohm_E(const MultiFab& centerBin,
                         const MultiFab& centerBtimeAvg, MultiFab& Eout,
                         int iLev, Real hstep, bool includeAmbi) {
  BL_PROFILE("Pic::assemble_ohm_E");

  // Nodal total current J = curl(B)/(4*pi) from trial B (compact 1*dx stencil
  // from cell centres to nodes). Only needed for physical resistivity and Hall.
  const bool needJ = (etaResistivity > 0 || useHallTerm);
  if (needJ) {
    curl_center_to_node(centerBin, nodeJ[iLev], Geom(iLev).InvCellSize());
    nodeJ[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  // Magnetic field interpolated from cell centres to nodes for vector cross
  // products.
  average_center_to_node(centerBtimeAvg, nodeBstage[iLev]);
  nodeBstage[iLev].FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(nodeStatus[iLev], nodeBstage[iLev], 0, 3, &Pic::get_node_B,
                   iLev, true);
  }

  // Moment time-interpolation weights:
  // X(hstep) = (0.5-hstep)*X^{n-1/2} + (0.5+hstep)*X^{n+1/2}.
  const Real wPrev = 0.5 - hstep;
  const Real wCur = 0.5 + hstep;
  const Real invFourPI = 1.0 / fourPI;

  for (MFIter mfi(Eout); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrE = Eout[mfi].array();
    const Array4<Real const>& arrB = nodeBstage[iLev][mfi].array();
    const Array4<Real const>& moments = nodePlasma[nSpecies][iLev][mfi].array();
    const Array4<Real const>& momentsPrev =
        nodePlasmaPrev[nSpecies][iLev][mfi].array();
    const Array4<Real const> arrJ =
        needJ ? nodeJ[iLev][mfi].array() : Array4<Real const>();
    const Array4<Real const> arrEambi = (electronTemperature > 0)
                                            ? nodeEambi[iLev][mfi].array()
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

      // J = curl(B)/(4*pi)
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

      // Ambipolar electric field: E_ambi = -grad(p_e)/(e*n_e)
      // Analytically curl(E_ambi) == 0 for isothermal/polytropic electrons.
      if (includeAmbi && electronTemperature > 0) {
        ex += arrEambi(i, j, k, ix_);
        ey += arrEambi(i, j, k, iy_);
        ez += arrEambi(i, j, k, iz_);
      }

      // Hall term: (J x B) / rho_q
      if (rho > 0 && useHallTerm) {
        const Real invRhoEff = 1.0 / amrex::max(rho, rhoMinOhm);
        Real hall_x = (jy * bz - jz * by) * invRhoEff;
        Real hall_y = (jz * bx - jx * bz) * invRhoEff;
        Real hall_z = (jx * by - jy * bx) * invRhoEff;

        ex += hall_x;
        ey += hall_y;
        ez += hall_z;
      }

      arrE(i, j, k, ix_) = ex;
      arrE(i, j, k, iy_) = ey;
      arrE(i, j, k, iz_) = ez;
    });
  }

  // Hyper-resistivity: E -= (eta_h / 4*pi) * curl(nabla^2 B).
  // centerLapB = Laplacian(centerBin);
  // nodeHyperE = curl_center_to_node(centerLapB).
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
// Ambipolar electric field: E_ambi = -grad(p_e) / (e * n_e).
// Precomputed once per PIC timestep outside the magnetic subcycling
// steps, avoiding redundant gradient and EOS operations during the
// high-frequency whistler integration.
void Pic::compute_ambipolar_E() {
  std::string nameFunc = "Pic::compute_ambipolar_E";
  timing_func(nameFunc);

  if (electronTemperature <= 0) {
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      nodeEambi[iLev].setVal(0.0);
    }
    return;
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    compute_ambipolar_E(iLev);
  }

  if (finest_level > 0) {
    for (int iLev = 1; iLev < n_lev(); ++iLev) {
      fill_fine_lev_bny_from_coarse(
          nodeEambi[iLev - 1], nodeEambi[iLev], 0, nDim3, ref_ratio[iLev - 1],
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    }
  }
}

//==========================================================
void Pic::compute_ambipolar_E(int iLev) {
  if (electronTemperature <= 0) {
    nodeEambi[iLev].setVal(0.0);
    return;
  }

  // Copy nodal ion density to nodeRhoTemp and fill periodic boundaries
  for (MFIter mfi(nodeRhoTemp[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrRho = nodeRhoTemp[iLev][mfi].array();
    const Array4<Real const>& moments = nodePlasma[nSpecies][iLev][mfi].array();
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      arrRho(i, j, k) = moments(i, j, k, iRho_);
    });
  }
  nodeRhoTemp[iLev].FillBoundary(Geom(iLev).periodicity());

  // Average nodal density to cell centres
  average_node_to_center(nodeRhoTemp[iLev], centerPe[iLev]);
  centerPe[iLev].FillBoundary(Geom(iLev).periodicity());

  // Evaluate electron pressure Pe at cell centres via EOS
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
        arrPe(i, j, k) = (r > 0) ? p0 * std::pow(r * invRho0, gamma) : 0.0;
      }
    });
  }
  centerPe[iLev].FillBoundary(Geom(iLev).periodicity());

  // Zero-gradient (Neumann / foextrap) BC across non-periodic domain boundaries
  if (!Geom(iLev).isAllPeriodic() && centerPe[iLev].nGrow() > 0) {
    Vector<BCRec> bcr(1);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      if (Geom(iLev).isPeriodic(d)) {
        bcr[0].setLo(d, BCType::int_dir);
        bcr[0].setHi(d, BCType::int_dir);
      } else {
        bcr[0].setLo(d, BCType::foextrap);
        bcr[0].setHi(d, BCType::foextrap);
      }
    }
    GpuBndryFuncFab<FabFillNoOp> bfunc(FabFillNoOp{});
    PhysBCFunct<GpuBndryFuncFab<FabFillNoOp> > physbcf(Geom(iLev), bcr, bfunc);
    physbcf(centerPe[iLev], 0, 1, centerPe[iLev].nGrowVect(), 0.0, 0);
    centerPe[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  if (isFake2D) {
    for (amrex::MFIter mfi(centerPe[iLev]); mfi.isValid(); ++mfi) {
      const auto& vbox = mfi.validbox();
      const auto& fbox = mfi.fabbox();
      auto arr = centerPe[iLev][mfi].array();
      const int klo = vbox.smallEnd(2);
      const int khi = vbox.bigEnd(2);
      amrex::ParallelFor(fbox,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                           const int k_src = std::clamp(k, klo, khi);
                           if (k != k_src) {
                             arr(i, j, k) = arr(i, j, k_src);
                           }
                         });
    }
  }

  // Compute grad_center_to_node(Pe) directly into nodeEambi
  grad_center_to_node(centerPe[iLev], nodeEambi[iLev],
                      Geom(iLev).InvCellSize());

  // Scale in-place: E_ambi = -grad(Pe) / max(rho, rhoMinOhm)
  for (MFIter mfi(nodeEambi[iLev]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrEambi = nodeEambi[iLev][mfi].array();
    const Array4<Real const>& moments = nodePlasma[nSpecies][iLev][mfi].array();

    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      const Real rho = moments(i, j, k, iRho_);
      if (rho > 0) {
        const Real invRhoEff = -1.0 / amrex::max(rho, rhoMinOhm);
        arrEambi(i, j, k, ix_) *= invRhoEff;
        arrEambi(i, j, k, iy_) *= invRhoEff;
        arrEambi(i, j, k, iz_) *= invRhoEff;
      } else {
        arrEambi(i, j, k, ix_) = 0.0;
        arrEambi(i, j, k, iy_) = 0.0;
        arrEambi(i, j, k, iz_) = 0.0;
      }
    });
  }

  // At physical conducting/reflecting walls, the normal ambipolar field
  // vanishes analytically (dPe/dn = 0)
  if (!Geom(iLev).isAllPeriodic()) {
    const BoundaryBounds bnd(Geom(iLev), nodeEambi[iLev].boxArray().ixType(),
                             &bcField);
    for (MFIter mfi(nodeEambi[iLev]); mfi.isValid(); ++mfi) {
      const Box& bxValid = mfi.validbox();
      Array4<Real> const& arr = nodeEambi[iLev][mfi].array();
      const Dim3 vLo = bxValid.smallEnd().dim3();
      const Dim3 vHi = bxValid.bigEnd().dim3();
      ParallelFor(bxValid, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const int ijk[3] = { i, j, k };
        const int vLoArr[3] = { vLo.x, vLo.y, vLo.z };
        const int vHiArr[3] = { vHi.x, vHi.y, vHi.z };
        for (int d = 0; d < nDim; ++d) {
          if (!bnd.isNode[d])
            continue;
          const bool onLoWall =
              (bnd.bcLo[d] == FieldBC::conducting) && (ijk[d] == bnd.loBnd[d]);
          const bool onHiWall =
              (bnd.bcHi[d] == FieldBC::conducting) && (ijk[d] == bnd.hiBnd[d]);
          if (onLoWall || onHiWall) {
            bool inValid = true;
            for (int od = 0; od < nDim; ++od) {
              if (od != d && (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
                inValid = false;
                break;
              }
            }
            if (inValid) {
              arr(i, j, k, d) = 0.0;
            }
          }
        }
      });
    }
  }

  if (isFake2D) {
    for (amrex::MFIter mfi(nodeEambi[iLev]); mfi.isValid(); ++mfi) {
      const auto& vbox = mfi.validbox();
      const auto& fbox = mfi.fabbox();
      auto arr = nodeEambi[iLev][mfi].array();
      const int klo = vbox.smallEnd(2);
      const int khi = vbox.bigEnd(2);
      amrex::ParallelFor(fbox,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                           const int k_src = std::clamp(k, klo, khi);
                           if (k != k_src) {
                             for (int n = 0; n < nDim3; ++n) {
                               arr(i, j, k, n) = arr(i, j, k_src, n);
                             }
                           }
                         });
    }
  }

  nodeEambi[iLev].FillBoundary(Geom(iLev).periodicity());
  apply_field_bc(nodeStatus[iLev], nodeEambi[iLev], 0, nDim3, &Pic::get_node_E,
                 iLev, false);
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

    if (!Geom(iLev).isAllPeriodic()) {
      const BoundaryBounds bnd(Geom(iLev), moments.boxArray().ixType(),
                               &bcField);
      for (MFIter mfi(moments); mfi.isValid(); ++mfi) {
        const Box& bxValid = mfi.validbox();
        Array4<Real> const& arr = moments[mfi].array();
        const Dim3 vLo = bxValid.smallEnd().dim3();
        const Dim3 vHi = bxValid.bigEnd().dim3();
        ParallelFor(bxValid, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
          const int ijk[3] = { i, j, k };
          const int vLoArr[3] = { vLo.x, vLo.y, vLo.z };
          const int vHiArr[3] = { vHi.x, vHi.y, vHi.z };
          for (int d = 0; d < nDim; ++d) {
            if (!bnd.isNode[d])
              continue;
            const bool onLoWall = (bnd.bcLo[d] == FieldBC::conducting) &&
                                  (ijk[d] == bnd.loBnd[d]);
            const bool onHiWall = (bnd.bcHi[d] == FieldBC::conducting) &&
                                  (ijk[d] == bnd.hiBnd[d]);
            if (onLoWall || onHiWall) {
              bool inValid = true;
              for (int od = 0; od < nDim; ++od) {
                if (od != d && (ijk[od] < vLoArr[od] || ijk[od] > vHiArr[od])) {
                  inValid = false;
                  break;
                }
              }
              if (inValid) {
                if (d == 0) {
                  arr(i, j, k, iMx_) = 0.0;
                  arr(i, j, k, iPxy_) = 0.0;
                  arr(i, j, k, iPxz_) = 0.0;
                } else if (d == 1) {
                  arr(i, j, k, iMy_) = 0.0;
                  arr(i, j, k, iPxy_) = 0.0;
                  arr(i, j, k, iPyz_) = 0.0;
                } else if (d == 2) {
                  arr(i, j, k, iMz_) = 0.0;
                  arr(i, j, k, iPxz_) = 0.0;
                  arr(i, j, k, iPyz_) = 0.0;
                }
              }
            }
          }
        });
      }
    }
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
    MultiFab::Copy(nodePlasmaPrev[nSpecies][iLev], nodePlasma[nSpecies][iLev],
                   0, 0, nHybridMomentsComps,
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
// BCs for the cell-centred B, applied to the RK trial states and to the new
// state at the end of each B sub-step.
void Pic::apply_centerB_BC(int iLev) { apply_centerB_BC(iLev, centerB[iLev]); }

void Pic::apply_centerB_BC(int iLev, amrex::MultiFab& mfB) {
  mfB.FillBoundary(Geom(iLev).periodicity());
  if (iLev == 0) {
    apply_field_bc(cellStatus[iLev], mfB, 0, mfB.nComp(), &Pic::get_center_B,
                   iLev, true);
  } else {
    MultiFab& coarseB = (&mfB == &centerBstage[iLev])  ? centerBstage[iLev - 1]
                        : (&mfB == &centerBstar[iLev]) ? centerBstar[iLev - 1]
                                                       : centerB[iLev - 1];
    fill_fine_lev_bny_from_coarse(
        coarseB, mfB, 0, mfB.nComp(), ref_ratio[iLev - 1], Geom(iLev - 1),
        Geom(iLev), cell_status(iLev), *get_cell_interp());
    apply_field_bc(cellStatus[iLev], mfB, 0, mfB.nComp(), &Pic::get_center_B,
                   iLev, true);
  }

  if (isFake2D) {
    const int nComp = mfB.nComp();
    for (amrex::MFIter mfi(mfB); mfi.isValid(); ++mfi) {
      const auto& vbox = mfi.validbox();
      const auto& fbox = mfi.fabbox();
      auto arr = mfB[mfi].array();
      const int klo = vbox.smallEnd(2);
      const int khi = vbox.bigEnd(2);
      amrex::ParallelFor(fbox,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                           const int k_src = std::clamp(k, klo, khi);
                           if (k != k_src) {
                             for (int n = 0; n < nComp; ++n) {
                               arr(i, j, k, n) = arr(i, j, k_src, n);
                             }
                           }
                         });
    }
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
  // CFL stability check for explicit resistive, hyper-resistive, and whistler
  // terms.
  const Real cflLimit = useRK4 ? 2.785 : 2.513;
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    if (etaResistivity <= 0 && etaHyperLev[iLev] <= 0 && !useHallTerm)
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

    if (useHallTerm) {
      Real bMaxLev = 0.0;
      for (int d = 0; d < 3; ++d) {
        bMaxLev = amrex::max(bMaxLev, centerB[iLev].norm0(d, 0, false));
      }
      const Real rhoMinLev = nodePlasma[nSpecies][iLev].min(iRho_, 0, false);
      const Real rhoEff = amrex::max(rhoMinLev, rhoMinOhm);
      const Real cflWhistler = (sMax * bMaxLev / (fourPI * rhoEff)) * subDt;
      if (cflWhistler > cflLimit) {
        const int nSubReq =
            static_cast<int>(std::ceil(nBSubcycle * (cflWhistler / cflLimit)));
        amrex::Print()
            << "  [CFL warning] whistler wave: "
               "(kmax^2*Bmax)/(4pi*rho_eff)*dt_sub = "
            << cflWhistler << " (> " << cflLimit
            << ", whistler wave is linearly unstable in RK4! rho_min="
            << rhoMinLev << ", rho_eff=" << rhoEff << ", B_max=" << bMaxLev
            << "). Recommend nBSubcycle >= " << nSubReq
            << " or increasing rhoMinOhm.\n";
      }
    }
  }

  // Precalculate the ambipolar electric field E_ambi = -grad(p_e)/(e*n_e)
  // once per PIC timestep outside the magnetic subcycling steps.
  compute_ambipolar_E();

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
                       hstepStart, false);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][0],
                            Geom(iLev).InvCellSize());

        // Stage 2: B2 = B^n - 0.5 dt k1; evaluate E at (B2 + B^n)/2
        MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0,
                          -0.5 * subDt, kStage[iLev][0], 0, 0, nDim3, nGst);
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerB[iLev], 0, 0, nDim3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstage[iLev], centerBstar[iLev], nodeEstage[iLev],
                       iLev, hstepHalf, false);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][1],
                            Geom(iLev).InvCellSize());

        // Stage 3: B3 = B^n - 0.5 dt k2; evaluate E at (B3 + B^n)/2
        MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0,
                          -0.5 * subDt, kStage[iLev][1], 0, 0, nDim3, nGst);
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerB[iLev], 0, 0, nDim3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstage[iLev], centerBstar[iLev], nodeEstage[iLev],
                       iLev, hstepHalf, false);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][2],
                            Geom(iLev).InvCellSize());

        // Stage 4: B4 = B^n - dt k3; evaluate E at (B4 + B^n)/2
        MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0, -subDt,
                          kStage[iLev][2], 0, 0, nDim3, nGst);
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerB[iLev], 0, 0, nDim3, nGst);
        apply_centerB_BC(iLev, centerBstage[iLev]);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstage[iLev], centerBstar[iLev], nodeEstage[iLev],
                       iLev, hstepEnd, false);
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
                       hstepStart, false);
        curl_node_to_center(nodeEstage[iLev], kStage[iLev][0],
                            Geom(iLev).InvCellSize());
        MultiFab::LinComb(centerBstage[iLev], 1.0, centerB[iLev], 0, -subDt,
                          kStage[iLev][0], 0, 0, nDim3, nGst);

        // Stage 2: B2 = (3/4)*B_n + (1/4)*(B1 - subDt * curl(E(avgB2)))
        MultiFab::LinComb(centerBstar[iLev], 0.5, centerBstage[iLev], 0, 0.5,
                          centerBstart[iLev], 0, 0, nDim3, nGst);
        apply_centerB_BC(iLev, centerBstar[iLev]);
        assemble_ohm_E(centerBstar[iLev], centerBstar[iLev], nodeEstage[iLev],
                       iLev, hstepEnd, false);
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
                       iLev, hstepHalf, false);
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

  // Fill coarse-fine interface ghost cells for nodeB.
  if (finest_level > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(
          nodeB[iLev - 1], nodeB[iLev], 0, nDim3, ref_ratio[iLev - 1],
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    }
  }

  // Evaluate E^{n+1} into nodeE for the next push.
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    assemble_ohm_E(centerB[iLev], centerB[iLev], nodeE[iLev], iLev, 1.0);
  }

  // Fill coarse-fine interface ghost cells for nodeE.
  if (finest_level > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(
          nodeE[iLev - 1], nodeE[iLev], 0, nDim3, ref_ratio[iLev - 1],
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
    }
  }

  // Suppress grid-scale E component if enabled.
  if (doSmoothE) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
      smooth_E(nodeE[iLev], iLev);
    }
  }
}

//==========================================================
void Pic::calc_hybrid_dt_and_subcycle(Real& dtNext, int& nSubNext,
                                      bool doReport) {
  std::string nameFunc = "Pic::calc_hybrid_dt_and_subcycle";
  timing_func(nameFunc);

  const Real cflLimit = useRK4 ? 2.785 : 2.513;
  const Real fieldSafety = 0.85;

  Real dtMacroMin = std::numeric_limits<Real>::max();
  Real dtSubMin = std::numeric_limits<Real>::max();
  std::string localReason = "ion kinetics";

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    const auto dx = Geom(iLev).CellSizeArray();
    Real dxMin = dx[0];
    for (int d = 1; d < nDim; ++d) {
      dxMin = amrex::min(dxMin, dx[d]);
    }

    // Ion particle motion
    Real uMaxLev = 0.0;
    if (fixedUMax >= 0) {
      uMaxLev = fixedUMax;
    } else {
      Vector<Real> uMaxSpecies(nSpecies, 0.0);
      for (int i = 0; i < nSpecies; ++i) {
        amrex::MultiFab& momMF = nodePlasma[i][iLev];
        uMaxSpecies[i] = parts[i]->calc_max_thermal_velocity(momMF);
      }
      ParallelDescriptor::ReduceRealMax(uMaxSpecies.data(), nSpecies);
      for (int i = 0; i < nSpecies; ++i) {
        if (uMaxSpecies[i] > uMaxLev) {
          uMaxLev = uMaxSpecies[i];
        }
      }
    }

    // Maximum magnetic field on this level
    Real bMaxLev = 0.0;
    for (int d = 0; d < 3; ++d) {
      bMaxLev = amrex::max(bMaxLev, centerB[iLev].norm0(d, 0, false));
    }

    // 1. Kinetic advection limit
    Real dtPart =
        (uMaxLev > 0.0) ? (dxMin / uMaxLev) : std::numeric_limits<Real>::max();

    // 2. Ion gyrofrequency limit (Omega_ci * dt <= thetaGyro)
    Real dtGyro = std::numeric_limits<Real>::max();
    if (bMaxLev > 0.0) {
      Real maxQom = 0.0;
      for (int i = 0; i < nSpecies; ++i) {
        Real qom = std::abs(parts[i]->get_charge() / parts[i]->get_mass());
        maxQom = amrex::max(maxQom, qom);
      }
      if (maxQom > 0.0) {
        const Real thetaGyro = 0.35; // radians (~18 steps per gyroperiod)
        dtGyro = thetaGyro / (maxQom * bMaxLev);
      }
    }

    Real dtMacroLev = amrex::min(dtPart, dtGyro);
    if (dtMacroLev < dtMacroMin) {
      dtMacroMin = dtMacroLev;
      if (dtGyro < dtPart) {
        localReason = "ion gyro-frequency";
      } else {
        localReason = "ion particle CFL";
      }
    }

    // 3. Field subcycling limits
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

    if (sMax > 0.0) {
      if (useHallTerm && bMaxLev > 0.0) {
        const Real rhoMinLev = nodePlasma[nSpecies][iLev].min(iRho_, 0, false);
        const Real rhoEff = amrex::max(
            rhoMinLev, amrex::max(rhoMinOhm, static_cast<Real>(1e-10)));
        const Real omegaWhistler = (sMax * bMaxLev) / (fourPI * rhoEff);
        if (omegaWhistler > 0.0) {
          Real dtWhistler = (cflLimit * fieldSafety) / omegaWhistler;
          if (dtWhistler < dtSubMin) {
            dtSubMin = dtWhistler;
            localReason = "whistler wave";
          }
        }
      }

      if (etaResistivity > 0.0) {
        const Real omegaEta = (etaResistivity / fourPI) * sMax;
        if (omegaEta > 0.0) {
          Real dtEta = (cflLimit * fieldSafety) / omegaEta;
          if (dtEta < dtSubMin) {
            dtSubMin = dtEta;
            localReason = "resistive diffusion";
          }
        }
      }

      if (etaHyperLev.size() > iLev && etaHyperLev[iLev] > 0.0 && lMax > 0.0) {
        const Real omegaHyper = (etaHyperLev[iLev] / fourPI) * sMax * lMax;
        if (omegaHyper > 0.0) {
          Real dtHyper = (cflLimit * fieldSafety) / omegaHyper;
          if (dtHyper < dtSubMin) {
            dtSubMin = dtHyper;
            localReason = "hyper-resistivity";
          }
        }
      }
    }
  }

  const Real userCFL = (tc->get_cfl() > 0.0) ? tc->get_cfl() : 0.2;
  const Real dtMacroTarget = userCFL * dtMacroMin;
  dtLimitingReason = localReason;

  if (tc->get_cfl() <= 0.0) {
    // Fixed macro dt mode
    Real dtFixed = tc->get_dt();
    if (dtFixed <= 0.0) {
      dtFixed = tc->get_next_dt();
    }
    dtNext = dtFixed;

    if (isAutoSubcycle && dtSubMin < std::numeric_limits<Real>::max()) {
      int nSubReq = static_cast<int>(std::ceil(dtFixed / dtSubMin));
      nSubReq = std::clamp(nSubReq, nBSubcycleMin, nBSubcycleMax);

      if (nSubReq > nBSubcycle) {
        nSubNext = nSubReq;
        subcycleHoldCount = 0;
      } else if (nSubReq < nBSubcycle) {
        if (dtFixed / (nBSubcycle - 1) <= 0.85 * dtSubMin) {
          subcycleHoldCount++;
          if (subcycleHoldCount >= 3) {
            nSubNext = nBSubcycle - 1;
            subcycleHoldCount = 0;
          } else {
            nSubNext = nBSubcycle;
          }
        } else {
          subcycleHoldCount = 0;
          nSubNext = nBSubcycle;
        }
      } else {
        nSubNext = nBSubcycle;
        subcycleHoldCount = 0;
      }
    } else {
      nSubNext = nBSubcycle;
    }
  } else {
    // Adaptive macro dt mode
    if (!isAutoSubcycle) {
      Real dtFieldLimit = nBSubcycle * dtSubMin;
      Real dtTarget = std::min(dtMacroTarget, dtFieldLimit);
      if (dtTarget == dtFieldLimit &&
          dtSubMin < std::numeric_limits<Real>::max()) {
        dtLimitingReason = "field advance (fixed subcycles)";
      }

      Real dtOld = tc->get_dt();
      if (dtOld > 0.0) {
        dtNext = std::min(dtTarget, static_cast<Real>(1.10) * dtOld);
        dtNext = std::max(dtNext, static_cast<Real>(0.50) * dtOld);
      } else {
        dtNext = dtTarget;
      }
      nSubNext = nBSubcycle;
    } else {
      int nIdeal = (dtSubMin < std::numeric_limits<Real>::max())
                       ? static_cast<int>(std::ceil(dtMacroTarget / dtSubMin))
                       : 1;
      Real dtTarget = dtMacroTarget;

      if (nIdeal > nBSubcycleMax &&
          dtSubMin < std::numeric_limits<Real>::max()) {
        nSubNext = nBSubcycleMax;
        dtTarget = nBSubcycleMax * dtSubMin;
        dtLimitingReason = "whistler / field (capped at max subcycles)";
      } else {
        int nTarget = std::clamp(nIdeal, nBSubcycleMin, nBSubcycleMax);
        if (nTarget > nBSubcycle) {
          nSubNext = nTarget;
          subcycleHoldCount = 0;
        } else if (nTarget < nBSubcycle) {
          if (dtTarget / (nBSubcycle - 1) <= 0.85 * dtSubMin) {
            subcycleHoldCount++;
            if (subcycleHoldCount >= 3) {
              nSubNext = nBSubcycle - 1;
              subcycleHoldCount = 0;
            } else {
              nSubNext = nBSubcycle;
            }
          } else {
            subcycleHoldCount = 0;
            nSubNext = nBSubcycle;
          }
        } else {
          nSubNext = nBSubcycle;
          subcycleHoldCount = 0;
        }
      }

      Real dtOld = tc->get_dt();
      if (dtOld > 0.0) {
        dtNext = std::min(dtTarget, static_cast<Real>(1.10) * dtOld);
        dtNext = std::max(dtNext, static_cast<Real>(0.50) * dtOld);
      } else {
        dtNext = dtTarget;
      }
    }
  }

  if (doReport && !domainParameters.doCompact) {
    amrex::Print() << printPrefix << "[Hybrid Adaptive Step] dt = " << dtNext
                   << ", nSub = " << nSubNext
                   << " (limited by: " << dtLimitingReason << ")\n";
  }
}
