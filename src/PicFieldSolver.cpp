#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_Algorithm.H>
#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "LinearSolver.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

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
    apply_field_bc(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
                   &Pic::get_center_B, iLev, true);
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
    apply_field_bc(nodeStatus[iLev], nodeE[iLev], 0, nDim3, &Pic::get_node_E,
                   iLev, false);
  }
}

//==========================================================
void Pic::update_E_impl() {
  std::string nameFunc = "Pic::update_E_impl";

  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    eSolver.reset(get_local_node_or_cell_number(nodeE[iLev]));

    update_E_rhs(eSolver.rhs, iLev);

    // Both solver modes compute one correction in eSolver.xLeft.
    if (fsolver.mode == FieldSolverMode::NewtonKrylov) {
      solve_E_newton_krylov(iLev);
    } else {
      solve_E_gmres(iLev);
    }

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

      // NOTE: the wave hard source is applied inside apply_field_bc().
      apply_field_bc(nodeStatus[iLev], nodeE[iLev], 0, nDim3, &Pic::get_node_E,
                     iLev, false);
      apply_field_bc(nodeStatus[iLev], nodeEth[iLev], 0, nDim3,
                     &Pic::get_node_E, iLev, false);
    }

    if (doSmoothE) {
      smooth_E(nodeEth[iLev], iLev);
      smooth_E(nodeE[iLev], iLev);
    }
    div_node_to_center(nodeE[iLev], centerDivE[iLev], Geom(iLev).InvCellSize());
  }
}

//==========================================================
void Pic::solve_E_gmres(int iLev) {
  convert_3d_to_1d(nodeE[iLev], eSolver.xLeft, iLev);

  update_E_matvec(eSolver.xLeft, eSolver.matvec, iLev, false);

  // Original linear solve: A * delta = rhs - A(E_old).
  for (int i = 0; i < eSolver.get_nSolve(); ++i) {
    eSolver.rhs[i] -= eSolver.matvec[i];
    eSolver.xLeft[i] = 0;
  }

  if (doReport)
    Print() << "\n-------" << printPrefix
            << " GMRES E solver ------------------" << std::endl;

  BL_PROFILE_VAR("Pic::E_iterate", eSolver);
  eSolver.solve(iLev, doReport);
  BL_PROFILE_VAR_STOP(eSolver);
}

//==========================================================
void Pic::solve_E_newton_krylov(int iLev) {
  const int nSolve = eSolver.get_nSolve();

  std::vector<double> base(nSolve);
  std::vector<double> baseMatvec(nSolve);
  std::vector<double> work(nSolve);

  convert_3d_to_1d(nodeE[iLev], base.data(), iLev);

  update_E_matvec(base.data(), baseMatvec.data(), iLev, false);

  // One Newton step: J(E_old) * delta = rhs - F(E_old).
  for (int i = 0; i < nSolve; ++i) {
    eSolver.rhs[i] -= baseMatvec[i];
    eSolver.xLeft[i] = 0;
  }

  if (doReport)
    Print() << "\n-------" << printPrefix
            << " Newton-Krylov E solver ------------------" << std::endl;

  BL_PROFILE_VAR("Pic::E_iterate", eSolver);
  const MPI_Comm iComm = ParallelDescriptor::Communicator();
  const double normBase =
      sqrt(dot_product_mpi(base.data(), base.data(), nSolve, iComm));

  // GMRES only sees this linearized Jacobian-vector product.
  auto jacobianFreeMatvec = [this, &base, &baseMatvec, &work, nSolve, normBase,
                             iComm](const double* vecIn, double* vecOut,
                                    const int iLevIn) {
    const double normDirection =
        sqrt(dot_product_mpi(vecIn, vecIn, nSolve, iComm));
    const double epsilon =
        fleks_jfnk::finite_difference_epsilon(normBase, normDirection);

    auto nonlinearMatvec = [this](const double* in, double* out,
                                  const int iLevMatvec) {
      update_E_matvec(in, out, iLevMatvec, false);
    };

    fleks_jfnk::jacobian_free_matvec(nonlinearMatvec, base.data(),
                                     baseMatvec.data(), vecIn, vecOut,
                                     work.data(), nSolve, iLevIn, epsilon);
  };

  eSolver.solve(jacobianFreeMatvec, iLev, doReport);
  BL_PROFILE_VAR_STOP(eSolver);
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
    // Even after apply_field_bc(), the outmost layer node E is still
    // unknow. See FluidInterface::calc_current for detailed explaniation.
    if (iLev == 0) {
      apply_field_bc(nodeStatus[iLev], vecMF, 0, nDim3, &Pic::get_node_E, iLev,
                     false);
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
      const Array4<Real>& arrE0 = nodeE[iLev][mfi].array();
      const Array4<Real>& limitE = fsolver.useLaggedLimiter ? arrE0 : arrE;
      const Array4<Real>& res = matvecMF[mfi].array();
      const Array4<Real>& arrU = uBg[iLev][mfi].array();

      ParallelFor(box, vecMF.nComp(), [&](int i, int j, int k, int iVar) {
        for (int iDir = 0; iDir < nDim; iDir++) {
          Real dii[nDim3] = { 0, 0, 0 };
          dii[iDir] = 1;

          Real cR = limiter_theta(
              limiterThetaE,
              limitE(i - dii[ix_], j - dii[iy_], k - dii[iz_], iVar),
              limitE(i, j, k, iVar),
              limitE(i + dii[ix_], j + dii[iy_], k + dii[iz_], iVar));

          Real cL = limiter_theta(
              limiterThetaE,
              limitE(i - 2 * dii[ix_], j - 2 * dii[iy_], k - 2 * dii[iz_],
                     iVar),
              limitE(i - dii[ix_], j - dii[iy_], k - dii[iz_], iVar),
              limitE(i, j, k, iVar));

          Real ur = cMaxE, ul = cMaxE;

          if (cMaxE <= 0) {
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
    apply_field_bc(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
                   &Pic::get_center_B, iLev, true);
    apply_field_bc(nodeStatus[iLev], nodeB[iLev], 0, nodeB[iLev].nComp(),
                   &Pic::get_node_B, iLev, true);
  } else {
    fill_fine_lev_bny_from_coarse(
        centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
        ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), cell_status(iLev),
        *get_cell_interp());

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
      apply_field_bc(cellStatus[iLev], centerB[iLev], 0, centerB[iLev].nComp(),
                     &Pic::get_center_B, iLev, true);

    } else {
      fill_fine_lev_bny_from_coarse(
          centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          *get_cell_interp());
    }
    MultiFab::Copy(dBdt[iLev], nodeB[iLev], 0, 0, dBdt[iLev].nComp(),
                   dBdt[iLev].nGrow());

    if (useHyperbolicCleaning) {
      div_node_to_center(nodeB[iLev], divB[iLev], Geom(iLev).InvCellSize());
    }

    if (useUpwindB || useHyperbolicCleaning) {
      correct_B(iLev);
    }

    average_center_to_node(centerB[iLev], nodeB[iLev]);
    nodeB[iLev].FillBoundary(Geom(iLev).periodicity());

    const Real invDt = 1. / tc->get_dt();
    // dBdt = (B^{n+1} - B^n)/dt;
    MultiFab::LinComb(dBdt[iLev], -invDt, dBdt[iLev], 0, invDt, nodeB[iLev], 0,
                      0, dBdt[iLev].nComp(), dBdt[iLev].nGrow());

    if (iLev == 0) {
      apply_field_bc(nodeStatus[iLev], nodeB[iLev], 0, nodeB[iLev].nComp(),
                     &Pic::get_node_B, iLev, true);
    } else {
      fill_fine_lev_bny_from_coarse(
          nodeB[iLev - 1], nodeB[iLev], 0, nodeB[iLev - 1].nComp(),
          ref_ratio[iLev - 1], Geom(iLev - 1), Geom(iLev), node_status(iLev),
          node_bilinear_interp);
    }
  }
}

//==========================================================
// Generalized Ohm's law
//   E = -U_i x B + eta J + (J x B)/rho_q - grad(Pe)/rho_q
// on cell-centred fields at an arbitrary (off-member) B state. `centerBin` is
// the trial B used for J = curl(B)/(4*pi); `centerBtimeAvg` is the
// time-averaged B used for the Hall/convection factor. Ion moments are
// time-interpolated between centerPlasmaPrev (J^{n-1/2}) and centerPlasmaSum
// (J^{n+1/2}) at hstep. Particle weights are initialized with dt = 1, so iRho_
// is the true charge density rho_q (the Hall / pressure terms divide by it
// directly).
//==========================================================

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
void Pic::correct_B(int iLev) {
  std::string nameFunc = "Pic::correct_B";
  timing_func(nameFunc);

  if (!useUpwindB && !useHyperbolicCleaning) {
    return;
  }

  MultiFab centerDB(cGrids[iLev], DistributionMap(iLev), nDim3, nGst);
  centerDB.setVal(0.0);

  if (useUpwindB) {
    Real coef[nDim3];
    for (int i = 0; i < nDim3; ++i) {
      coef[i] = 0.5 * tc->get_dt() * Geom(iLev).InvCellSize()[i];
    }

    for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
      Box box = mfi.validbox();

      const Array4<Real>& cB = centerB[iLev][mfi].array();
      const Array4<Real const>& nU = uBg[iLev][mfi].array();
      const Array4<Real>& dB = centerDB[mfi].array();
      const auto& status = cellStatus[iLev][mfi].array();

      // Get the face along the direction iDir for the cell (i,j,k) for the iVar
      // component
      auto get_face = [&](int iDir, int i, int j, int k, int iVar,
                          Array4<Real const> const& arr, Real& l, Real& r) {
        // Generic fixed-upwind-velocity override (e.g. the old TopHat
        // "bypass_limiter" used a constant speed of 1.0). Zero keeps the
        // normal plasma-background-velocity reconstruction below.
        if (fixedUpwindVel > 0) {
          l = fixedUpwindVel;
          r = fixedUpwindVel;
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

      ParallelFor(box, [&](int i, int j, int k) {
        bool doDiffusion;
        Real lu[nDim3] = { 0, 0, 0 }, ru[nDim3] = { 0, 0, 0 };
        Real ul, ur;

        IntVect ijk{ AMREX_D_DECL(i, j, k) };

        // Flux along  x
        for (int iDir = 0; iDir < nDim; iDir++) {
          get_face(iDir, i, j, k, iDir, nU, lu[iDir], ru[iDir]);
        }

        ul = lu[ix_];
        ur = ru[ix_];
        doDiffusion = true;
        if ((ul > 0 && bit::is_domain_boundary(status(i - 1, j, k))) ||
            (ur < 0 && bit::is_domain_boundary(status(i + 1, j, k)))) {
          doDiffusion = false;
        }

        if (doDiffusion) {
          ul = fabs(ul);
          ur = fabs(ur);
          for (int iVar = 0; iVar < nDim3; iVar++) {
            Real cR = limiter_theta(limiterThetaB, cB(i - 1, j, k, iVar),
                                    cB(i, j, k, iVar), cB(i + 1, j, k, iVar));
            Real cL = limiter_theta(limiterThetaB, cB(i - 2, j, k, iVar),
                                    cB(i - 1, j, k, iVar), cB(i, j, k, iVar));
            ul = min(ul, 0.5 / coef[ix_]);
            ur = min(ur, 0.5 / coef[ix_]);
            dB(i, j, k, iVar) +=
                (cR * ur * (cB(i + 1, j, k, iVar) - cB(i, j, k, iVar)) -
                 cL * ul * (cB(i, j, k, iVar) - cB(i - 1, j, k, iVar))) *
                coef[ix_];
          }
        }

        // Flux along y
        ul = lu[iy_];
        ur = ru[iy_];
        doDiffusion = true;
        if ((ul > 0 && bit::is_domain_boundary(status(i, j - 1, k))) ||
            (ur < 0 && bit::is_domain_boundary(status(i, j + 1, k)))) {
          doDiffusion = false;
        }

        if (doDiffusion) {
          ul = fabs(ul);
          ur = fabs(ur);

          for (int iVar = 0; iVar < nDim3; iVar++) {
            Real cR = limiter_theta(limiterThetaB, cB(i, j - 1, k, iVar),
                                    cB(i, j, k, iVar), cB(i, j + 1, k, iVar));
            Real cL = limiter_theta(limiterThetaB, cB(i, j - 2, k, iVar),
                                    cB(i, j - 1, k, iVar), cB(i, j, k, iVar));
            ul = min(ul, 0.5 / coef[iy_]);
            ur = min(ur, 0.5 / coef[iy_]);

            dB(i, j, k, iVar) +=
                (cR * ur * (cB(i, j + 1, k, iVar) - cB(i, j, k, iVar)) -
                 cL * ul * (cB(i, j, k, iVar) - cB(i, j - 1, k, iVar))) *
                coef[iy_];
          }
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

          if (doDiffusion) {
            ul = fabs(ul);
            ur = fabs(ur);

            for (int iVar = 0; iVar < nDim3; iVar++) {
              Real cR = limiter_theta(limiterThetaB, cB(i, j, k - 1, iVar),
                                      cB(i, j, k, iVar), cB(i, j, k + 1, iVar));
              Real cL = limiter_theta(limiterThetaB, cB(i, j, k - 2, iVar),
                                      cB(i, j, k - 1, iVar), cB(i, j, k, iVar));
              ul = min(ul, 0.5 / coef[iz_]);
              ur = min(ur, 0.5 / coef[iz_]);

              dB(i, j, k, iVar) +=
                  (cR * ur * (cB(i, j, k + 1, iVar) - cB(i, j, k, iVar)) -
                   cL * ul * (cB(i, j, k, iVar) - cB(i, j, k - 1, iVar))) *
                  coef[iz_];
            }
          }
        }
      });
    } // end MFIter
  } // end useUpwindB

  if (useHyperbolicCleaning) {
    MultiFab gradPhi(cGrids[iLev], DistributionMap(iLev), nDim3, 0);
    gradPhi.setVal(0.0);

    solve_hyp_phi(iLev);

    MultiFab gradPhiNode(nGrids[iLev], DistributionMap(iLev), nDim3, 0);
    gradPhiNode.setVal(0.0);

    grad_center_to_node(hypPhi[iLev], gradPhiNode, Geom(iLev).InvCellSize());

    average_node_to_cellcenter(gradPhi, 0, gradPhiNode, 0, nDim3,
                               gradPhi.nGrow());

    for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
      Box box = mfi.validbox();

      const Array4<Real>& dB = centerDB[mfi].array();
      const Array4<Real>& gradPhiArr = gradPhi[mfi].array();

      ParallelFor(box, [&](int i, int j, int k) {
        IntVect ijk{ AMREX_D_DECL(i, j, k) };
        for (int iVar = 0; iVar < nDim3; iVar++) {
          dB(ijk, iVar) += -tc->get_dt() * gradPhiArr(ijk, iVar);
        }
      });
    } // end MFIter
  } // end useHyperbolicCleaning

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
// Apply boundary conditions to electromagnetic fields (E or B).
// Visits each boundary face and applies the corresponding BC configured in
// bcField.
//==========================================================

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
