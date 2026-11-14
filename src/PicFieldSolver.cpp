#include <algorithm>
#include <cmath>
#include <vector>

#include <AMReX_Algorithm.H>
#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "LinearSolver.h"
#include "Pic.h"
#include "PicUpwind.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
void Pic::update_E() {
  if (useExplicitPIC) {
    update_E_expl();
  } else {
    update_E_impl();
  }
  if (finestLevel() > 0) {
    if (useHybridPIC) {
      sync_emf_at_interface(nodeE);
      fill_ghost_emf_from_coarse(nodeE);
      sync_emf_at_interface(nodeEth);
      fill_ghost_emf_from_coarse(nodeEth);
    } else if (useCurlUpwindB) {
      fill_interface_and_ghost_emf_from_coarse(nodeE);
      fill_interface_and_ghost_emf_from_coarse(nodeEth);
      for (int iLev = 1; iLev < n_lev(); ++iLev) {
        div_node_to_center(nodeE[iLev], centerDivE[iLev],
                           Geom(iLev).InvCellSize());
      }
    }
  }
}

//==========================================================
void Pic::update_E_expl() {
  std::string nameFunc = "Pic::update_E_expl";

  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    MultiFab::Copy(nodeEth[iLev], nodeE[iLev], 0, 0, nodeE[iLev].nComp(),
                   nodeE[iLev].nGrow());
    if (useHybridPIC || useCurlUpwindB) {
      apply_centerB_BC(iLev);
    } else {
      centerB[iLev].FillBoundary(Geom(iLev).periodicity());
      if (iLev == 0) {
        apply_field_bc(cell_status(iLev), centerB[iLev], 0,
                       centerB[iLev].nComp(), &Pic::get_center_B, iLev, true);
      } else {
        fill_fine_lev_bny_from_coarse(
            centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
            refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), cell_status(iLev),
            *get_cell_interp());
      }
    }
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
    apply_field_bc(node_status(iLev), nodeE[iLev], 0, nDim3, &Pic::get_node_E,
                   iLev, false);

    // Electric field boundary on the inner body (see #BODYBOUNDARY).
    apply_body_E_bc(nodeE[iLev], iLev);
  }
}

//==========================================================
// The tangential electric field on a conducting wall node is a Dirichlet
// unknown of the implicit E solve (E_t = 0): its matrix row is the identity
// and its right-hand side is E^n, so the GMRES correction there is zero.
AMREX_GPU_HOST_DEVICE inline bool is_conducting_wall_Et(
    const BoundaryBounds& bnd, const int ijk[3], int c) {
  for (int d = 0; d < nDim; ++d) {
    if (!bnd.isNode[d] || c == d)
      continue;
    if ((bnd.bcLo[d] == FieldBC::conducting && ijk[d] == bnd.loBnd[d]) ||
        (bnd.bcHi[d] == FieldBC::conducting && ijk[d] == bnd.hiBnd[d]))
      return true;
  }
  return false;
}

//==========================================================
void Pic::update_E_impl() {
  std::string nameFunc = "Pic::update_E_impl";

  timing_func(nameFunc);

  if (hasAbsorbBC_ && !isAbsorbCflWarned_)
    warn_absorb_cfl();

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (fsolver.solverBackend == ESolverBackend::AmrexGMRES) {
      solve_E_gmres_amrex(iLev);
    } else {
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
    }

    if (doSmoothE) {
      smooth_E(nodeEth[iLev], iLev);
    }

    MultiFab::Add(nodeEth[iLev], nodeE[iLev], 0, 0, nodeEth[iLev].nComp(),
                  get_n_ghost());

    MultiFab::LinComb(nodeE[iLev], -(1.0 - fsolver.theta) / fsolver.theta,
                      nodeE[iLev], 0, 1. / fsolver.theta, nodeEth[iLev], 0, 0,
                      nodeE[iLev].nComp(), get_n_ghost());

    if (iLev == 0) {

      // NOTE: the wave hard source is applied inside apply_field_bc().
      apply_field_bc(node_status(iLev), nodeE[iLev], 0, nDim3, &Pic::get_node_E,
                     iLev, false);
      apply_field_bc(node_status(iLev), nodeEth[iLev], 0, nDim3,
                     &Pic::get_node_E, iLev, false);
    } else {
      fill_fine_lev_bny_from_coarse(
          nodeE[iLev - 1], nodeE[iLev], 0, nDim3, refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
      nodeE[iLev].FillBoundary(Geom(iLev).periodicity());
      apply_field_bc(node_status(iLev), nodeE[iLev], 0, nDim3, &Pic::get_node_E,
                     iLev, false);

      fill_fine_lev_bny_from_coarse(
          nodeEth[iLev - 1], nodeEth[iLev], 0, nDim3, refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
      nodeEth[iLev].FillBoundary(Geom(iLev).periodicity());
      apply_field_bc(node_status(iLev), nodeEth[iLev], 0, nDim3,
                     &Pic::get_node_E, iLev, false);
    }

    if (doSmoothE) {
      smooth_E(nodeEth[iLev], iLev);
      smooth_E(nodeE[iLev], iLev);
    }

    // Electric field boundary on the inner body (see #BODYBOUNDARY). For
    // 'linetied' the body nodes are not part of the linear system, so nodeEth
    // is zero there by construction and nodeE is zeroed here, after the last
    // operation that can write into the body (the smoothing above). For
    // 'conducting' both are projected onto the radial direction, which also
    // removes the round-off of the Krylov solve.
    apply_body_E_bc(nodeE[iLev], iLev);
    apply_body_E_bc(nodeEth[iLev], iLev);

    div_node_to_center(nodeE[iLev], centerDivE[iLev], Geom(iLev).InvCellSize());
  }
}

//==========================================================
void Pic::solve_E_gmres(int iLev) {
  convert_3d_to_1d(nodeE[iLev], eSolver.xLeft, iLev);

  update_E_matvec(eSolver.xLeft, eSolver.matvec, iLev, false);

  // Original linear solve: A * delta = rhs - A(E_old).
  Real* const r = eSolver.rhs;
  const Real* const m = eSolver.matvec;
  Real* const x = eSolver.xLeft;
  parallel_for_host_1d(eSolver.get_nSolve(), [=](int i) {
    r[i] -= m[i];
    x[i] = 0;
  });

  if (doReport)
    Print() << "\n-------" << printPrefix
            << " GMRES E solver ------------------" << std::endl;

  BL_PROFILE_VAR("Pic::E_iterate", eSolver);
  eSolver.solve(iLev, doReport);
  BL_PROFILE_VAR_STOP(eSolver);
}

//==========================================================
void Pic::sync_nodal_solver_vec(MultiFab& vecMF, int iLev) {
  const bool skipAllBody = is_body_linetied();
  const bool skipInterior = is_body_conducting();
  const auto& statusMF = node_status(iLev);

  for (MFIter mfi(vecMF, doTiling); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const Array4<Real> arr = vecMF[mfi].array();
    const Array4<int const> nodeArr = statusMF[mfi].const_array();
    parallel_for_omp(bx, vecMF.nComp(),
                     [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                       const int st = nodeArr(i, j, k);
                       if (!bit::is_owner(st) ||
                           (skipAllBody && bit::is_body(st)) ||
                           (skipInterior && bit::is_body_interior(st))) {
                         arr(i, j, k, c) = 0.0;
                       }
                     });
  }

  vecMF.SumBoundary(0, vecMF.nComp(), IntVect(0), IntVect(0),
                    Geom(iLev).periodicity());
}

//==========================================================
void Pic::solve_E_gmres_amrex(int iLev) {
  if (static_cast<int>(amrexESolver.size()) < n_lev_max()) {
    amrexESolver.resize(n_lev_max());
  }
  if (!amrexESolver[iLev]) {
    amrexESolver[iLev] = std::make_unique<AmrexLinearSolver>(this, iLev);
  }
  amrexESolver[iLev]->define(this, iLev, fsolver.precondType,
                             fsolver.precondIters, fsolver.smootherWeight,
                             fsolver.mode);

  if (solverRhsNode2[iLev].empty()) {
    distribute_FabArray(solverRhsNode1[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
    distribute_FabArray(solverRhsNode2[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
  }
  if (solverMatvecMF[iLev].empty()) {
    distribute_FabArray(solverVecMF[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
    distribute_FabArray(solverMatvecMF[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, 1);
    distribute_FabArray(solverTempNode3[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
    distribute_FabArray(solverCenterLapMF[iLev], boxArray(iLev),
                        DistributionMap(iLev), 3, 1);
  }

  MultiFab& rhsMF = solverRhsNode2[iLev];
  update_E_rhs(rhsMF, iLev);

  MultiFab& matvecMF = solverMatvecMF[iLev];
  update_E_matvec(nodeE[iLev], matvecMF, iLev, false);

  if (fsolver.mode == FieldSolverMode::NewtonKrylov) {
    amrexESolver[iLev]->set_newton_base(nodeE[iLev], matvecMF);
  }

  MultiFab::Saxpy(rhsMF, -1.0, matvecMF, 0, 0, 3, 0);

  if (doReport) {
    Print() << "\n-------" << printPrefix
            << (fsolver.mode == FieldSolverMode::NewtonKrylov
                    ? " Newton-Krylov E solver ------------------"
                    : " GMRES E solver ------------------")
            << std::endl;
    Print() << "Before 0 nMatVec, relative tol: " << eSolver.get_nIter() << " "
            << eSolver.get_tol() << std::endl;
  }

  const Real res0 = std::sqrt(
      mf_dot(rhsMF, 0, rhsMF, 0, 3, 0, &amrexESolver[iLev]->dot_mask()));
  if (doReport) {
    Print() << "initial residual norm: " << res0 << std::endl;
  }

#ifdef GMRES_VERBOSE
  const int verbose = doReport ? 2 : 0;
#else
  const int verbose = 0;
#endif

  int nMatvec = 0;
  Real resFinal = 0.0;
  BL_PROFILE_VAR("Pic::E_iterate", amrex_gmres);
  amrexESolver[iLev]->solve(nodeEth[iLev], rhsMF, eSolver.get_tol(), 0.0,
                            eSolver.get_nIter(), fsolver.gmresRestartLength,
                            verbose, nMatvec, resFinal);
  BL_PROFILE_VAR_STOP(amrex_gmres);

  const Real relErr = (res0 > 0.0) ? (resFinal / res0) : 0.0;
  eSolver.set_stats(res0, nMatvec, relErr);

  if (doReport) {
    Print() << "After nMatVec, relative error = " << nMatvec << " " << relErr
            << std::endl;
  }

  nodeEth[iLev].SumBoundary(Geom(iLev).periodicity());
  nodeEth[iLev].FillBoundary(Geom(iLev).periodicity());
}

//==========================================================
void Pic::solve_E_newton_krylov(int iLev) {
  const int nSolve = eSolver.get_nSolve();

  std::vector<Real> base(nSolve);
  std::vector<Real> baseMatvec(nSolve);
  std::vector<Real> work(nSolve);

  convert_3d_to_1d(nodeE[iLev], base.data(), iLev);

  update_E_matvec(base.data(), baseMatvec.data(), iLev, false);

  // One Newton step: J(E_old) * delta = rhs - F(E_old).
  Real* const r = eSolver.rhs;
  const Real* const bm = baseMatvec.data();
  Real* const x = eSolver.xLeft;
  parallel_for_host_1d(nSolve, [=](int i) {
    r[i] -= bm[i];
    x[i] = 0;
  });

  if (doReport)
    Print() << "\n-------" << printPrefix
            << " Newton-Krylov E solver ------------------" << std::endl;

  BL_PROFILE_VAR("Pic::E_iterate", eSolver);
  const MPI_Comm iComm = ParallelDescriptor::Communicator();
  const Real normBase =
      sqrt(dot_product_mpi(base.data(), base.data(), nSolve, iComm));

  // GMRES only sees this linearized Jacobian-vector product.
  auto jacobianFreeMatvec = [this, &base, &baseMatvec, &work, nSolve, normBase,
                             iComm](const Real* vecIn, Real* vecOut,
                                    const int iLevIn) {
    const Real normDirection =
        sqrt(dot_product_mpi(vecIn, vecIn, nSolve, iComm));
    const Real epsilon =
        fleks_jfnk::finite_difference_epsilon(normBase, normDirection);

    auto nonlinearMatvec = [this](const Real* in, Real* out,
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
void Pic::update_E_matvec(const Real* vecIn, Real* vecOut, int iLev,
                          const bool useZeroBC) {
  if (solverVecMF[iLev].empty()) {
    distribute_FabArray(solverVecMF[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
    distribute_FabArray(solverMatvecMF[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, 1);
    distribute_FabArray(solverTempNode3[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
    distribute_FabArray(solverCenterLapMF[iLev], boxArray(iLev),
                        DistributionMap(iLev), 3, 1);
  }

  MultiFab& vecMF = solverVecMF[iLev];
  MultiFab& matvecMF = solverMatvecMF[iLev];

  convert_1d_to_3d(vecIn, vecMF, iLev);

  // The right side edges should be filled in.
  vecMF.SumBoundary(Geom(iLev).periodicity());

  update_E_matvec(vecMF, matvecMF, iLev, useZeroBC);

  convert_3d_to_1d(matvecMF, vecOut, iLev);
}

//==========================================================
void Pic::update_E_matvec(const MultiFab& vecIn, MultiFab& matvecMF, int iLev,
                          const bool useZeroBC) {
  std::string nameFunc = "Pic::E_matvec";
  timing_func(nameFunc);

  if (solverVecMF[iLev].empty()) {
    distribute_FabArray(solverVecMF[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
    distribute_FabArray(solverMatvecMF[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, 1);
    distribute_FabArray(solverTempNode3[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
    distribute_FabArray(solverCenterLapMF[iLev], boxArray(iLev),
                        DistributionMap(iLev), 3, 1);
  }

  const bool hasWallEt = hasConductingBC_;
  if (hasWallEt && solverWallInMF[iLev].empty())
    distribute_FabArray(solverWallInMF[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, 0);

  MultiFab& vecMF = solverVecMF[iLev];
  if (&vecIn != &vecMF) {
    vecMF.setVal(0.0);
    MultiFab::Copy(vecMF, vecIn, 0, 0, 3, 0);
    if (vecIn.nGrow() == 0) {
      vecMF.SumBoundary(0, 3, IntVect(0), IntVect(0), Geom(iLev).periodicity());
    } else {
      sync_nodal_solver_vec(vecMF, iLev);
    }
  }

  MultiFab& tempNode3 = solverTempNode3[iLev];

  // M*E needs ghost cell information.
  vecMF.FillBoundary(Geom(iLev).periodicity());

  if (is_fake_2d()) {
    // Make sure there is no variation in the z-direction.
    Periodicity period(IntVect(AMREX_D_DECL(0, 0, 1)));
    vecMF.FillBoundary(period);
  }

  if (hasWallEt)
    MultiFab::Copy(solverWallInMF[iLev], vecMF, 0, 0, 3, 0);

  const BoundaryBounds bnd(Geom(iLev), vecMF.boxArray().ixType(), &bcField);
  bool hasDirichlet = false;
  for (int d = 0; d < nDim; ++d) {
    if (bnd.isNode[d] && (FieldBC::is_dirichlet(bnd.bcLo[d]) ||
                          FieldBC::is_dirichlet(bnd.bcHi[d]))) {
      hasDirichlet = true;
    }
  }

  if (useZeroBC) {
    // For Dirichlet physical boundaries (inflow, fixed, wave), the correction
    // in the Krylov subspace must be identically zero on the boundary nodes.
    if (hasDirichlet) {
      for (MFIter mfi(vecMF); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const Array4<Real> arrVec = vecMF[mfi].array();
        parallel_for_omp(bx, vecMF.nComp(),
                         [=] AMREX_GPU_DEVICE(int i, int j, int k, int c)
                             noexcept {
                               const int ijk[3] = { i, j, k };
                               for (int d = 0; d < nDim; ++d) {
                                 if (!bnd.isNode[d])
                                   continue;
                                 if ((FieldBC::is_dirichlet(bnd.bcLo[d]) &&
                                      ijk[d] == bnd.loBnd[d]) ||
                                     (FieldBC::is_dirichlet(bnd.bcHi[d]) &&
                                      ijk[d] == bnd.hiBnd[d])) {
                                   arrVec(i, j, k, c) = 0.0;
                                 }
                               }
                             });
      }
    }

    // The correction obeys the linear, homogeneous part of the field BC:
    // E_t = 0 and mirrored ghosts on conducting walls, zero-gradient ghosts
    // on outflow faces. Without it GMRES solved with zero ghosts, a different
    // boundary than the one imposed on E afterwards, and the solver went
    // unstable above c*dt/dx ~ 2-3.
    apply_field_bc_homogeneous(node_status(iLev), vecMF, iLev);
  } else {
    // Even after apply_field_bc(), the outmost layer node E is still
    // unknow. See FluidInterface::calc_current for detailed explaniation.
    if (iLev == 0) {
      apply_field_bc(node_status(iLev), vecMF, 0, nDim3, &Pic::get_node_E, iLev,
                     false);
    } else {
      fill_fine_lev_bny_from_coarse(
          nodeEth[iLev - 1], vecMF, 0, nodeEth[iLev - 1].nComp(),
          refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), node_status(iLev),
          node_bilinear_interp);
      vecMF.FillBoundary(Geom(iLev).periodicity());
      apply_field_bc(node_status(iLev), vecMF, 0, nDim3, &Pic::get_node_E, iLev,
                     false);
    }
  }

  lap_node_to_node(vecMF, matvecMF, DistributionMap(iLev), Geom(iLev),
                   &solverCenterLapMF[iLev]);

  Real delt2 = pow(fsolver.theta * tc->get_dt(), 2);
  mf_mult(matvecMF, -delt2, 0, matvecMF.nComp(), 0);

  if (useUpwindE) {
    // Explicit scheme: add the LF artificial viscosity term to the rhs
    // vis_{i+0.5} = c_max/2*(E_i+1 - E_i)
    // E_i += dt/dx*(vis_{i+0.5} - vis_{i-0.5}) = 0.5*c_max*dt*dx*lap(E_i)
    // For implicit scheme, we add it to the lhs, so the sign changes.

    const Real dx = Geom(iLev).CellSize()[0];
    const Real coe1 = -0.5 * fsolver.theta * tc->get_dt() / dx;
    const Real limThetaE = limiterThetaE;
    const Real cMaxVal = cMaxE;

    for (MFIter mfi(vecMF); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real const> arrE = vecMF[mfi].const_array();
      const Array4<Real const> arrE0 = nodeE[iLev][mfi].const_array();
      const Array4<Real const> limitE = fsolver.useLaggedLimiter ? arrE0 : arrE;
      const Array4<Real> res = matvecMF[mfi].array();
      const Array4<Real const> arrU = uBg[iLev][mfi].const_array();

      parallel_for_omp(
          box, vecMF.nComp(),
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
            for (int iDir = 0; iDir < nDim; iDir++) {
              int dii[nDim3] = { 0, 0, 0 };
              dii[iDir] = 1;

              Real cR = limiter_theta(
                  limThetaE,
                  limitE(i - dii[ix_], j - dii[iy_], k - dii[iz_], iVar),
                  limitE(i, j, k, iVar),
                  limitE(i + dii[ix_], j + dii[iy_], k + dii[iz_], iVar));

              Real cL = limiter_theta(
                  limThetaE,
                  limitE(i - 2 * dii[ix_], j - 2 * dii[iy_], k - 2 * dii[iz_],
                         iVar),
                  limitE(i - dii[ix_], j - dii[iy_], k - dii[iz_], iVar),
                  limitE(i, j, k, iVar));

              Real ur = cMaxVal, ul = cMaxVal;

              if (cMaxVal <= 0) {
                ul = fabs(0.5 * (arrU(i - dii[ix_], j - dii[iy_], k - dii[iz_],
                                      iDir) +
                                 arrU(i, j, k, iDir)));
                ur = fabs(0.5 * (arrU(i, j, k, iDir) +
                                 arrU(i + dii[ix_], j + dii[iy_], k + dii[iz_],
                                      iDir)));
              }

              Real dE =
                  cR * ur *
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
      if (solverTempCenter3[iLev].empty()) {
        distribute_FabArray(solverTempCenter3[iLev], boxArray(iLev),
                            DistributionMap(iLev), 3, get_n_ghost());
        distribute_FabArray(solverTempCenter1[iLev], boxArray(iLev),
                            DistributionMap(iLev), 1, get_n_ghost());
      }
      MultiFab& tempCenter3 = solverTempCenter3[iLev];
      MultiFab& tempCenter1 = solverTempCenter1[iLev];
      tempCenter3.setVal(0.0);
      tempCenter1.setVal(0.0);

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
      apply_BC(cell_status(iLev), tempCenter1, 0, tempCenter1.nComp(),
               &Pic::get_zero, iLev);

      MultiFab::LinComb(centerDivE[iLev], 1 - fsolver.coefDiff,
                        centerDivE[iLev], 0, fsolver.coefDiff, tempCenter1, 0,
                        0, 1, 1);
    }

    grad_center_to_node(centerDivE[iLev], tempNode3, Geom(iLev).InvCellSize());

    mf_saxpy(matvecMF, delt2, tempNode3, 0, 0, matvecMF.nComp(), 0);
  }

  update_E_M_dot_E(vecMF, tempNode3, iLev);

  mf_add(matvecMF, tempNode3, 0, 0, matvecMF.nComp(), 0);

  mf_add(matvecMF, vecMF, 0, 0, matvecMF.nComp(), 0);

  // Identity rows of the tangential E on conducting walls, from the input
  // before the wall BC zeroed it.
  if (hasWallEt) {
    for (MFIter mfi(matvecMF); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.validbox();
      const Array4<Real> arrOut = matvecMF[mfi].array();
      const Array4<Real const> arrIn = solverWallInMF[iLev][mfi].const_array();
      parallel_for_omp(bx, matvecMF.nComp(),
                       [=] AMREX_GPU_DEVICE(int i, int j, int k, int c)
                           noexcept {
                             const int ijk[3] = { i, j, k };
                             if (is_conducting_wall_Et(bnd, ijk, c))
                               arrOut(i, j, k, c) = arrIn(i, j, k, c);
                           });
    }
  }

  // Dirichlet physical boundaries (inflow, fixed, wave): enforce A_ii = 1, A_ij
  // = 0 by copying vecMF to matvecMF on boundary nodes.
  if (hasDirichlet) {
    for (MFIter mfi(matvecMF); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.validbox();
      const Array4<Real> arrOut = matvecMF[mfi].array();
      const Array4<Real const> arrIn = vecMF[mfi].const_array();
      parallel_for_omp(bx, matvecMF.nComp(),
                       [=] AMREX_GPU_DEVICE(int i, int j, int k, int c)
                           noexcept {
                             const int ijk[3] = { i, j, k };
                             for (int d = 0; d < nDim; ++d) {
                               if (!bnd.isNode[d])
                                 continue;
                               if ((FieldBC::is_dirichlet(bnd.bcLo[d]) &&
                                    ijk[d] == bnd.loBnd[d]) ||
                                   (FieldBC::is_dirichlet(bnd.bcHi[d]) &&
                                    ijk[d] == bnd.hiBnd[d])) {
                                 arrOut(i, j, k, c) = arrIn(i, j, k, c);
                               }
                             }
                           });
    }
  }

  // 'conducting' body: keep the Krylov vectors inside the constrained space
  // by removing the tangential electric field on the body nodes. Together with
  // the same projection of the right-hand side and of the solution, this makes
  // E_t = 0 an exact constraint of the linear system.
  if (is_body_conducting())
    project_body_E(matvecMF, iLev);
}

//==========================================================
void Pic::update_E_M_dot_E(const MultiFab& inMF, MultiFab& outMF, int iLev) {
  std::string nameFunc = "Pic::update_E_M_dot_E";
  timing_func(nameFunc);

  Real c0 = fourPI * fsolver.theta * tc->get_dt();
  parallel_for_mf_omp(
      outMF, [](const MFIter& mfi) { return mfi.validbox(); },
      [&](const MFIter& mfi) {
        const Array4<Real const> inArr = inMF[mfi].const_array();
        const Array4<Real> outArr = outMF[mfi].array();
        const Array4<RealMM const> mmArr = nodeMM[iLev][mfi].const_array();
        return [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          IntVect ijk = { AMREX_D_DECL(i, j, k) };

          const auto& data0 = mmArr(ijk);
          Real sumX = 0.0;
          Real sumY = 0.0;
          Real sumZ = 0.0;

#if (AMREX_SPACEDIM == 2)
          for (int j2 = j - 1; j2 <= j + 1; ++j2) {
            for (int i2 = i - 1; i2 <= i + 1; ++i2) {
              const int gp = (j2 - j + 1) * 3 + i2 - i + 1;
              const int idx0 = gp * 9;

              const Real* AMREX_RESTRICT const M_I = &(data0[idx0]);

              const Real vctX = inArr(i2, j2, 0, ix_);
              const Real vctY = inArr(i2, j2, 0, iy_);
              const Real vctZ = inArr(i2, j2, 0, iz_);
              sumX += vctX * M_I[0] + vctY * M_I[1] + vctZ * M_I[2];
              sumY += vctX * M_I[3] + vctY * M_I[4] + vctZ * M_I[5];
              sumZ += vctX * M_I[6] + vctY * M_I[7] + vctZ * M_I[8];
            }
          }
#else
          const int kMin2 = nDim > 2 ? k - 1 : k;
          const int kMax2 = nDim > 2 ? k + 1 : k;

          for (int k2 = kMin2; k2 <= kMax2; ++k2) {
            for (int j2 = j - 1; j2 <= j + 1; ++j2) {
              for (int i2 = i - 1; i2 <= i + 1; ++i2) {
                const int gp = (k2 - k + 1) * 9 + (j2 - j + 1) * 3 + i2 - i + 1;
                const int idx0 = gp * 9;

                const Real* AMREX_RESTRICT const M_I = &(data0[idx0]);

                const Real vctX = inArr(i2, j2, k2, ix_);
                const Real vctY = inArr(i2, j2, k2, iy_);
                const Real vctZ = inArr(i2, j2, k2, iz_);
                sumX += vctX * M_I[0] + vctY * M_I[1] + vctZ * M_I[2];
                sumY += vctX * M_I[3] + vctY * M_I[4] + vctZ * M_I[5];
                sumZ += vctX * M_I[6] + vctY * M_I[7] + vctZ * M_I[8];
              }
            }
          }
#endif
          outArr(i, j, k, ix_) = sumX * c0;
          outArr(i, j, k, iy_) = sumY * c0;
          outArr(i, j, k, iz_) = sumZ * c0;
        };
      });

  // if (doSmoothJ) {
  //   for (int icount = 0; icount < nSmoothJ; icount++) {
  //     smooth_multifab(outMF, iLev, icount % 2 + 1);
  //   }
  // }
}

//==========================================================
void Pic::update_E_rhs(Real* rhs, int iLev) {
  if (solverRhsNode1[iLev].empty()) {
    distribute_FabArray(solverRhsNode1[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
    distribute_FabArray(solverRhsNode2[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
  }

  MultiFab& temp2Node = solverRhsNode2[iLev];

  update_E_rhs(temp2Node, iLev);

  convert_3d_to_1d(temp2Node, rhs, iLev);
}

//==========================================================
void Pic::update_E_rhs(MultiFab& temp2Node, int iLev) {
  std::string nameFunc = "Pic::update_E_rhs";
  timing_func(nameFunc);

  if (solverRhsNode1[iLev].empty()) {
    distribute_FabArray(solverRhsNode1[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
    distribute_FabArray(solverRhsNode2[iLev], node_box_array(iLev),
                        DistributionMap(iLev), 3, get_n_ghost());
  }

  MultiFab& tempNode = solverRhsNode1[iLev];
  tempNode.setVal(0.0);
  temp2Node.setVal(0.0);

  if (useHybridPIC || useCurlUpwindB) {
    apply_centerB_BC(iLev);
  } else {
    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
    if (iLev == 0) {
      apply_field_bc(cell_status(iLev), centerB[iLev], 0, centerB[iLev].nComp(),
                     &Pic::get_center_B, iLev, true);
    } else {
      fill_fine_lev_bny_from_coarse(
          centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
          refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          *get_cell_interp());
    }
  }

  nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
  // The same base-versus-fine distinction applies to both B layouts used by
  // the implicit field solve after a regrid.
  if (iLev == 0) {
    apply_field_bc(node_status(iLev), nodeB[iLev], 0, nodeB[iLev].nComp(),
                   &Pic::get_node_B, iLev, true);
  } else {
    fill_fine_lev_bny_from_coarse(nodeB[iLev - 1], nodeB[iLev], 0,
                                  nodeB[iLev - 1].nComp(), refRatio(iLev - 1),
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

  // Dirichlet physical boundaries (inflow, fixed, wave): the right-hand side
  // for the total field must equal the boundary value that update_E_matvec()
  // puts into A(E_old), so that the initial residual rhs - A(E_old) is exactly
  // zero on Dirichlet boundary nodes. The Krylov corrections are zero there, so
  // any residual left on these nodes is never removed and grows. On level 0
  // that value is the prescribed state at the solve time (including the #WAVEBC
  // source), not the nodeE of the previous step.
  const BoundaryBounds bnd(Geom(iLev), temp2Node.boxArray().ixType(), &bcField);
  bool hasDirichlet = false;
  for (int d = 0; d < nDim; ++d) {
    if (bnd.isNode[d] && (FieldBC::is_dirichlet(bnd.bcLo[d]) ||
                          FieldBC::is_dirichlet(bnd.bcHi[d]))) {
      hasDirichlet = true;
    }
  }

  if (hasDirichlet) {
    const MultiFab* dirichletE = &nodeE[iLev];
    if (hasInflowBC_) {
      MultiFab::Copy(tempNode, nodeE[iLev], 0, 0, nodeE[iLev].nComp(),
                     tempNode.nGrow());
      apply_field_bc(node_status(iLev), tempNode, 0, nDim3, &Pic::get_node_E,
                     iLev, false);
      dirichletE = &tempNode;
    }

    for (MFIter mfi(temp2Node); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.validbox();
      const Array4<Real> arrRHS = temp2Node[mfi].array();
      const Array4<Real const> arrE = (*dirichletE)[mfi].const_array();
      parallel_for_omp(bx, temp2Node.nComp(),
                       [=] AMREX_GPU_DEVICE(int i, int j, int k, int c)
                           noexcept {
                             const int ijk[3] = { i, j, k };
                             for (int d = 0; d < nDim; ++d) {
                               if (!bnd.isNode[d])
                                 continue;
                               if ((FieldBC::is_dirichlet(bnd.bcLo[d]) &&
                                    ijk[d] == bnd.loBnd[d]) ||
                                   (FieldBC::is_dirichlet(bnd.bcHi[d]) &&
                                    ijk[d] == bnd.hiBnd[d])) {
                                 arrRHS(i, j, k, c) = arrE(i, j, k, c);
                               }
                             }
                           });
    }
  }

  // Tangential E on conducting walls: rhs = E^n (zero after the wall BC).
  if (hasConductingBC_) {
    for (MFIter mfi(temp2Node); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.validbox();
      const Array4<Real> arrRHS = temp2Node[mfi].array();
      const Array4<Real const> arrE = nodeE[iLev][mfi].const_array();
      parallel_for_omp(bx, temp2Node.nComp(),
                       [=] AMREX_GPU_DEVICE(int i, int j, int k, int c)
                           noexcept {
                             const int ijk[3] = { i, j, k };
                             if (is_conducting_wall_Et(bnd, ijk, c))
                               arrRHS(i, j, k, c) = arrE(i, j, k, c);
                           });
    }
  }

  // 'conducting' body: project the right-hand side onto the radial direction
  // so that the Krylov space stays inside the constrained space.
  if (is_body_conducting())
    project_body_E(temp2Node, iLev);
}

//==========================================================
void Pic::build_E_preconditioner(MultiFab& invDiagMF, int iLev,
                                 EPrecondType type) {
  std::string nameFunc = "Pic::build_E_precond";
  timing_func(nameFunc);

  const int nCompInv = (type == EPrecondType::Jacobi) ? 3 : 9;
  const BoxArray& ba = node_box_array(iLev);
  const DistributionMapping& dm = DistributionMap(iLev);
  if (invDiagMF.empty() || invDiagMF.boxArray() != ba ||
      invDiagMF.DistributionMap() != dm || invDiagMF.nComp() != nCompInv) {
    invDiagMF.define(ba, dm, nCompInv, 0);
  }

  const Real dtTheta = fsolver.theta * tc->get_dt();
  const Real delt2 = dtTheta * dtTheta;
  const Real c0 = fourPI * dtTheta;
  const Real* invDx = Geom(iLev).InvCellSize();
  const Real wDim = (nDim > 2) ? 0.5 : 1.0;

  Real lapDiag = 0.0;
  for (int d = 0; d < nDim; ++d) {
    lapDiag += invDx[d] * invDx[d];
  }

  GpuArray<Real, 3> emDiag = { 0.0, 0.0, 0.0 };
  const Real c2cWeight = (nDim > 2) ? (1.0 / 18.0) : (1.0 / 6.0);
  const Real gradDivCoeff =
      (1.0 - fsolver.coefDiff) * wDim + fsolver.coefDiff * c2cWeight;
  for (int a = 0; a < 3; ++a) {
    const Real gradDivDiag =
        (a < nDim) ? gradDivCoeff * invDx[a] * invDx[a] : 0.0;
    emDiag[a] = delt2 * (wDim * lapDiag - gradDivDiag);
  }

  const Real dx0 = Geom(iLev).CellSize()[0];
  const Real coeUpwind = useUpwindE ? (0.5 * dtTheta / dx0) : 0.0;
  const Real cMaxVal = cMaxE;
  const bool doUpwindE = useUpwindE;

  const bool skipAllBody = is_body_linetied();
  const bool skipInterior = is_body_conducting();
  const bool hasWallEt = hasConductingBC_;
  const BoundaryBounds bnd(Geom(iLev), ba.ixType(), &bcField);
  const auto& statusMF = node_status(iLev);

  for (MFIter mfi(invDiagMF, doTiling); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const Array4<Real> invArr = invDiagMF[mfi].array();
    const Array4<RealMM const> mmArr = nodeMM[iLev][mfi].const_array();
    const Array4<int const> statusArr = statusMF[mfi].const_array();
    const Array4<Real const> arrU = uBg[iLev][mfi].const_array();

    parallel_for_omp(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int st = statusArr(i, j, k);
      const bool isExcludedBody = (skipAllBody && bit::is_body(st)) ||
                                  (skipInterior && bit::is_body_interior(st));
      const int ijk[3] = { i, j, k };
      bool isDirichletNode = false;
      for (int d = 0; d < nDim; ++d) {
        if (!bnd.isNode[d])
          continue;
        if ((FieldBC::is_dirichlet(bnd.bcLo[d]) && ijk[d] == bnd.loBnd[d]) ||
            (FieldBC::is_dirichlet(bnd.bcHi[d]) && ijk[d] == bnd.hiBnd[d])) {
          isDirichletNode = true;
        }
      }

      if (isExcludedBody || isDirichletNode) {
        if (nCompInv == 3) {
          invArr(i, j, k, 0) = 1.0;
          invArr(i, j, k, 1) = 1.0;
          invArr(i, j, k, 2) = 1.0;
        } else {
          for (int m = 0; m < 9; ++m)
            invArr(i, j, k, m) = 0.0;
          invArr(i, j, k, 0) = 1.0;
          invArr(i, j, k, 4) = 1.0;
          invArr(i, j, k, 8) = 1.0;
        }
        return;
      }

      Real upwindDiag = 0.0;
      if (doUpwindE) {
        for (int iDir = 0; iDir < nDim; ++iDir) {
          int dii[3] = { 0, 0, 0 };
          dii[iDir] = 1;
          Real ur = cMaxVal;
          Real ul = cMaxVal;
          if (cMaxVal <= 0.0) {
            ul =
                std::abs(0.5 * (arrU(i - dii[0], j - dii[1], k - dii[2], iDir) +
                                arrU(i, j, k, iDir)));
            ur = std::abs(0.5 *
                          (arrU(i, j, k, iDir) +
                           arrU(i + dii[0], j + dii[1], k + dii[2], iDir)));
          }
          upwindDiag += coeUpwind * (ur + ul);
        }
      }

#if (AMREX_SPACEDIM == 2)
      constexpr int idx0 = 4 * 9;
#else
      constexpr int idx0 = 13 * 9;
#endif
      const Real* const M_I = &(mmArr(i, j, k)[idx0]);

      Real D[3][3];
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          D[a][b] = c0 * M_I[a * 3 + b];
        }
        D[a][a] += 1.0 + emDiag[a] + upwindDiag;
      }

      if (hasWallEt) {
        for (int c = 0; c < 3; ++c) {
          if (is_conducting_wall_Et(bnd, ijk, c)) {
            D[c][0] = 0.0;
            D[c][1] = 0.0;
            D[c][2] = 0.0;
            D[0][c] = 0.0;
            D[1][c] = 0.0;
            D[2][c] = 0.0;
            D[c][c] = 1.0;
          }
        }
      }

      if (nCompInv == 3) {
        for (int c = 0; c < 3; ++c) {
          invArr(i, j, k, c) =
              (std::abs(D[c][c]) > 1e-30) ? (1.0 / D[c][c]) : 1.0;
        }
      } else {
        const Real c00 = D[1][1] * D[2][2] - D[1][2] * D[2][1];
        const Real c01 = -(D[1][0] * D[2][2] - D[1][2] * D[2][0]);
        const Real c02 = D[1][0] * D[2][1] - D[1][1] * D[2][0];

        const Real c10 = -(D[0][1] * D[2][2] - D[0][2] * D[2][1]);
        const Real c11 = D[0][0] * D[2][2] - D[0][2] * D[2][0];
        const Real c12 = -(D[0][0] * D[2][1] - D[0][1] * D[2][0]);

        const Real c20 = D[0][1] * D[1][2] - D[0][2] * D[1][1];
        const Real c21 = -(D[0][0] * D[1][2] - D[0][2] * D[1][0]);
        const Real c22 = D[0][0] * D[1][1] - D[0][1] * D[1][0];

        const Real det = D[0][0] * c00 + D[0][1] * c01 + D[0][2] * c02;
        if (std::abs(det) > 1e-30) {
          const Real invDet = 1.0 / det;
          invArr(i, j, k, 0) = c00 * invDet;
          invArr(i, j, k, 1) = c10 * invDet;
          invArr(i, j, k, 2) = c20 * invDet;
          invArr(i, j, k, 3) = c01 * invDet;
          invArr(i, j, k, 4) = c11 * invDet;
          invArr(i, j, k, 5) = c21 * invDet;
          invArr(i, j, k, 6) = c02 * invDet;
          invArr(i, j, k, 7) = c12 * invDet;
          invArr(i, j, k, 8) = c22 * invDet;
        } else {
          for (int m = 0; m < 9; ++m)
            invArr(i, j, k, m) = 0.0;
          invArr(i, j, k, 0) = 1.0;
          invArr(i, j, k, 4) = 1.0;
          invArr(i, j, k, 8) = 1.0;
        }
      }
    });
  }
}

//==========================================================
const Pic::NodeSystemLayout& Pic::node_system_layout(const MultiFab& MF,
                                                     int iLev) {
  if (static_cast<int>(nodeLayout.size()) < n_lev_max())
    nodeLayout.resize(n_lev_max());

  // A 'linetied' body drops all its nodes from the linear system (E = 0),
  // a 'conducting' body only drops the interior ones (E = 0 inside, E_t = 0
  // on the surface), and an 'insulating' body keeps all of them.
  const bool skipAllBody = is_body_linetied();
  const bool skipInterior = is_body_conducting();

  // The selection depends only on the grid (node ownership) and on the body,
  // which is static, so the layout is rebuilt only after a regrid.
  NodeSystemLayout& layout = nodeLayout[iLev];
  if (layout.boxes.size() == static_cast<size_t>(MF.local_size()) &&
      layout.ba == MF.boxArray() && layout.dm == MF.DistributionMap() &&
      layout.skipAllBody == skipAllBody &&
      layout.skipInterior == skipInterior) {
    return layout;
  }

  layout.ba = MF.boxArray();
  layout.dm = MF.DistributionMap();
  layout.skipAllBody = skipAllBody;
  layout.skipInterior = skipInterior;
  layout.boxes.assign(MF.local_size(), {});

  HostFabCache statusCache;
  for (MFIter mfi(MF); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<int const> nodeArr = statusCache.array(node_status(iLev), mfi);
    const Dim3 lo = lbound(box);
    const Dim3 hi = ubound(box);

    auto is_unknown = [&](int i, int j, int k) {
      const int st = nodeArr(i, j, k);
      return bit::is_owner(st) && !(skipAllBody && bit::is_body(st)) &&
             !(skipInterior && bit::is_body_interior(st));
    };

    NodeSystemLayout::BoxLayout& boxLayout = layout.boxes[mfi.LocalIndex()];
    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        int i = lo.x;
        while (i <= hi.x) {
          while (i <= hi.x && !is_unknown(i, j, k))
            ++i;
          if (i > hi.x)
            break;
          const int iBeg = i;
          while (i <= hi.x && is_unknown(i, j, k))
            ++i;
          const int iEnd = i - 1;
          boxLayout.runs.push_back({ iBeg, iEnd, j, k, boxLayout.nUnknown });
          boxLayout.nUnknown += iEnd - iBeg + 1;
        }
      }
    }
  }

  return layout;
}

//==========================================================
void Pic::convert_1d_to_3d(const Real* const p, MultiFab& MF, int iLev) {
  std::string nameFunc = "Pic::convert_1d_to_3d";
  timing_func(nameFunc);

  MF.setVal(0.0);
  HostStaging<FArrayBox> stageMF(MF);

  const int nComp = MF.nComp();
  int offset = 0;
  if (MF.ixType().cellCentered()) {
    for (MFIter mfi(MF); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      const int ni = hi.x - lo.x + 1;
      const int nj = hi.y - lo.y + 1;
      const int nk = hi.z - lo.z + 1;
      const int nUnknown = ni * nj * nk;
      const Array4<Real> arr = MF[mfi].array();
      const Real* const pBox = p + offset;
      const int totalRows = nComp * nk * nj;

      auto copyRow = [&](int row) {
        const int iVar = row / (nk * nj);
        const int rem = row - iVar * (nk * nj);
        const int k = lo.z + rem / nj;
        const int j = lo.y + rem % nj;
        Real* AMREX_RESTRICT dst = arr.ptr(lo.x, j, k, iVar);
        const Real* AMREX_RESTRICT src =
            pBox + iVar * nUnknown + ((k - lo.z) * nj + (j - lo.y)) * ni;
        fleks_omp_detail::eval_simd_1d(0, ni,
                                       [&](int ii) { dst[ii] = src[ii]; });
      };

      if (nComp * nUnknown >= kMinOmp1DSize) {
        parallel_for_host_items(totalRows, copyRow, false);
      } else {
        for (int row = 0; row < totalRows; ++row)
          copyRow(row);
      }
      offset += nComp * nUnknown;
    }
    return;
  }

  const NodeSystemLayout& layout = node_system_layout(MF, iLev);
  for (MFIter mfi(MF); mfi.isValid(); ++mfi) {
    const NodeSystemLayout::BoxLayout& boxLayout =
        layout.boxes[mfi.LocalIndex()];
    const int nUnknown = boxLayout.nUnknown;
    const int nRuns = static_cast<int>(boxLayout.runs.size());
    const Array4<Real> arr = MF[mfi].array();
    const Real* const pBox = p + offset;
    const int totalTasks = nComp * nRuns;

    auto copyRun = [&](int task) {
      const int iVar = task / nRuns;
      const NodeSystemLayout::XRun& run = boxLayout.runs[task - iVar * nRuns];
      const int len = run.iEnd - run.iBeg + 1;
      Real* AMREX_RESTRICT dst = arr.ptr(run.iBeg, run.j, run.k, iVar);
      const Real* AMREX_RESTRICT src = pBox + iVar * nUnknown + run.vecOffset;
      fleks_omp_detail::eval_simd_1d(0, len,
                                     [&](int ii) { dst[ii] = src[ii]; });
    };

    if (nComp * nUnknown >= kMinOmp1DSize) {
      parallel_for_host_items(totalTasks, copyRun, false);
    } else {
      for (int task = 0; task < totalTasks; ++task)
        copyRun(task);
    }
    offset += nComp * nUnknown;
  }
}

//==========================================================
void Pic::convert_3d_to_1d(const MultiFab& MF, Real* const p, int iLev) {
  std::string nameFunc = "Pic::convert_3d_to_1d";
  timing_func(nameFunc);

  HostFabCache mfCache;
  const int nComp = MF.nComp();
  int offset = 0;

  if (MF.ixType().cellCentered()) {
    for (MFIter mfi(MF); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Dim3 lo = lbound(box);
      const Dim3 hi = ubound(box);
      const int ni = hi.x - lo.x + 1;
      const int nj = hi.y - lo.y + 1;
      const int nk = hi.z - lo.z + 1;
      const int nUnknown = ni * nj * nk;
      const Array4<Real const> arr = mfCache.array(MF, mfi);
      Real* const pBox = p + offset;
      const int totalRows = nComp * nk * nj;

      auto copyRow = [&](int row) {
        const int iVar = row / (nk * nj);
        const int rem = row - iVar * (nk * nj);
        const int k = lo.z + rem / nj;
        const int j = lo.y + rem % nj;
        const Real* AMREX_RESTRICT src = arr.ptr(lo.x, j, k, iVar);
        Real* AMREX_RESTRICT dst =
            pBox + iVar * nUnknown + ((k - lo.z) * nj + (j - lo.y)) * ni;
        fleks_omp_detail::eval_simd_1d(0, ni,
                                       [&](int ii) { dst[ii] = src[ii]; });
      };

      if (nComp * nUnknown >= kMinOmp1DSize) {
        parallel_for_host_items(totalRows, copyRow, false);
      } else {
        for (int row = 0; row < totalRows; ++row)
          copyRow(row);
      }
      offset += nComp * nUnknown;
    }
    return;
  }

  const NodeSystemLayout& layout = node_system_layout(MF, iLev);
  for (MFIter mfi(MF); mfi.isValid(); ++mfi) {
    const NodeSystemLayout::BoxLayout& boxLayout =
        layout.boxes[mfi.LocalIndex()];
    const int nUnknown = boxLayout.nUnknown;
    const int nRuns = static_cast<int>(boxLayout.runs.size());
    const Array4<Real const> arr = mfCache.array(MF, mfi);
    Real* const pBox = p + offset;
    const int totalTasks = nComp * nRuns;

    auto copyRun = [&](int task) {
      const int iVar = task / nRuns;
      const NodeSystemLayout::XRun& run = boxLayout.runs[task - iVar * nRuns];
      const int len = run.iEnd - run.iBeg + 1;
      const Real* AMREX_RESTRICT src = arr.ptr(run.iBeg, run.j, run.k, iVar);
      Real* AMREX_RESTRICT dst = pBox + iVar * nUnknown + run.vecOffset;
      fleks_omp_detail::eval_simd_1d(0, len,
                                     [&](int ii) { dst[ii] = src[ii]; });
    };

    if (nComp * nUnknown >= kMinOmp1DSize) {
      parallel_for_host_items(totalTasks, copyRun, false);
    } else {
      for (int task = 0; task < totalTasks; ++task)
        copyRun(task);
    }
    offset += nComp * nUnknown;
  }

  // Explicitly zero trailing slots up to the solver dimension to prevent
  // uninitialized or stale memory from polluting GMRES Arnoldi inner products.
  const int nSolve = eSolver.get_nSolve();
  if (offset < nSolve) {
    Real* const pTail = p + offset;
    parallel_for_host_1d(nSolve - offset, [=](int i) { pTail[i] = 0.0; });
  }
}

//==========================================================
void Pic::ensure_divB(int iLev) {
  // Full PIC allocates divB in distribute_arrays(); hybrid PIC has no such call
  // path, so it is created on demand and only when a deck asks for it.
  if (static_cast<int>(divB.size()) != n_lev_max()) {
    divB.resize(n_lev_max());
  }
  if (divB[iLev].empty() || divB[iLev].boxArray() != boxArray(iLev) ||
      divB[iLev].DistributionMap() != DistributionMap(iLev)) {
    // div_node_to_center() only fills component 0, so one component is enough.
    distribute_FabArray(divB[iLev], boxArray(iLev), DistributionMap(iLev), 1,
                        get_n_ghost(), false, 0.0);
  }
}

//==========================================================
void Pic::ensure_centerDivB(int iLev) {
  // Diagnostic only, so always on demand.
  if (static_cast<int>(centerDivB.size()) != n_lev_max()) {
    centerDivB.resize(n_lev_max());
  }
  if (centerDivB[iLev].empty() ||
      centerDivB[iLev].boxArray() != boxArray(iLev) ||
      centerDivB[iLev].DistributionMap() != DistributionMap(iLev)) {
    // div_center_to_center() only fills component 0.
    distribute_FabArray(centerDivB[iLev], boxArray(iLev), DistributionMap(iLev),
                        1, get_n_ghost(), false, 0.0);
  }
}

//==========================================================
void Pic::ensure_hypPhi(int iLev) {
  // On demand: distribute_arrays() does not provision it for hybrid PIC.
  if (static_cast<int>(hypPhi.size()) != n_lev_max()) {
    hypPhi.resize(n_lev_max());
  }
  if (hypPhi[iLev].empty() || hypPhi[iLev].boxArray() != boxArray(iLev) ||
      hypPhi[iLev].DistributionMap() != DistributionMap(iLev)) {
    distribute_FabArray(hypPhi[iLev], boxArray(iLev), DistributionMap(iLev),
                        nDim3, get_n_ghost(), false, 0.0);
  }
}

//==========================================================
void Pic::compute_divB(int iLev) {
  std::string nameFunc = "Pic::compute_divB";
  timing_func(nameFunc);

  // div of the nodal field, the quantity the cleaning consumes.
  ensure_divB(iLev);
  div_node_to_center(nodeB[iLev], divB[iLev], Geom(iLev).InvCellSize());

  // The same for the cell-centred field the update advances; 'divBc'.
  ensure_centerDivB(iLev);
  div_center_to_center(centerB[iLev], centerDivB[iLev],
                       Geom(iLev).InvCellSize());
}

//==========================================================
void Pic::update_B() {
  std::string nameFunc = "Pic::update_B";
  timing_func(nameFunc);

  if (finestLevel() > 0) {
    if (useHybridPIC) {
      sync_emf_at_interface(nodeEth);
      fill_ghost_emf_from_coarse(nodeEth);
      sync_emf_at_interface(nodeE);
      fill_ghost_emf_from_coarse(nodeE);
    } else if (useCurlUpwindB) {
      fill_interface_and_ghost_emf_from_coarse(nodeEth);
      fill_interface_and_ghost_emf_from_coarse(nodeE);
    }
  }

  // When curl-form upwind is active in the full-PIC solver, compute the upwind
  // dissipation as an effective EMF on nodes: E_upwind = sum_d (a_d/2) e_d x
  // jump(B). Faraday's curl update then advances B strictly divergence-free.
  const bool doCurlUpwind = useUpwindB && useCurlUpwindB && !useHybridPIC;
  if (doCurlUpwind) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (nodeEmf[iLev].empty()) {
        distribute_FabArray(nodeEmf[iLev], node_box_array(iLev),
                            DistributionMap(iLev), 3, get_n_ghost(), false,
                            0.0);
      }
      MultiFab::Copy(nodeEmf[iLev], nodeEth[iLev], 0, 0, nDim3,
                     nodeEth[iLev].nGrow());
      add_curl_upwind_E(nodeEmf[iLev], iLev);
      nodeEmf[iLev].OverrideSync(Geom(iLev).periodicity());
      nodeEmf[iLev].FillBoundary(Geom(iLev).periodicity());
      apply_field_bc(node_status(iLev), nodeEmf[iLev], 0, nDim3,
                     &Pic::get_node_E, iLev, false);
      if (use_body()) {
        apply_body_E_bc(nodeEmf[iLev], iLev);
      }
    }
    if (finestLevel() > 0) {
      fill_interface_and_ghost_emf_from_coarse(nodeEmf);
    }
  }

  Vector<MultiFab>& faradayEmf = doCurlUpwind ? nodeEmf : nodeEth;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (centerDB[iLev].empty()) {
      distribute_FabArray(centerDB[iLev], boxArray(iLev), DistributionMap(iLev),
                          nDim3, get_n_ghost());
    }
    MultiFab& dB = centerDB[iLev];
    curl_node_to_center(faradayEmf[iLev], dB, Geom(iLev).InvCellSize());

    // The interior of the body is a cavity: Faraday's law does not change the
    // magnetic field there, so B keeps its initial value.
    if (is_body_interior_frozen())
      mask_body_interior(dB, cell_status(iLev));

    MultiFab::Saxpy(centerB[iLev], -tc->get_dt(), dB, 0, 0,
                    centerB[iLev].nComp(), centerB[iLev].nGrow());

    centerB[iLev].FillBoundary(Geom(iLev).periodicity());
  }

  // Divergence-free relaxation of covered coarse B toward fine average
  if ((useHybridPIC || useCurlUpwindB) && finestLevel() > 0) {
    for (int iLev = finestLevel() - 1; iLev >= 0; iLev--) {
      relax_covered_B_to_fine(iLev);
    }
  } else if (projectDownEmFields && finestLevel() > 0) {
    for (int iLev = finestLevel(); iLev > 0; iLev--) {
      average_down(centerB[iLev], centerB[iLev - 1], 0, 3, refRatio(0));
    }
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (useHybridPIC || useCurlUpwindB) {
      apply_centerB_BC(iLev);
    } else {
      centerB[iLev].FillBoundary(Geom(iLev).periodicity());
      if (iLev == 0) {
        apply_field_bc(cell_status(iLev), centerB[iLev], 0,
                       centerB[iLev].nComp(), &Pic::get_center_B, iLev, true);
      } else {
        fill_fine_lev_bny_from_coarse(
            centerB[iLev - 1], centerB[iLev], 0, centerB[iLev - 1].nComp(),
            refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), cell_status(iLev),
            *get_cell_interp());
      }
    }

    MultiFab::Copy(dBdt[iLev], nodeB[iLev], 0, 0, dBdt[iLev].nComp(),
                   dBdt[iLev].nGrow());

    if (need_divB()) {
      // The cleaning used to be fed the PREVIOUS step's nodal field, nodeB
      // being rebuilt only further down. Average the current centerB in first;
      // the average below repeats it so that nodeB also carries the correction.
      average_center_to_node(centerB[iLev], nodeB[iLev]);
      nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
      div_node_to_center(nodeB[iLev], divB[iLev], Geom(iLev).InvCellSize());
      // Also refresh 'divBc' (see compute_divB).
      ensure_centerDivB(iLev);
      div_center_to_center(centerB[iLev], centerDivB[iLev],
                           Geom(iLev).InvCellSize());
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
      apply_field_bc(node_status(iLev), nodeB[iLev], 0, nodeB[iLev].nComp(),
                     &Pic::get_node_B, iLev, true);
    } else {
      fill_fine_lev_bny_from_coarse(nodeB[iLev - 1], nodeB[iLev], 0,
                                    nodeB[iLev - 1].nComp(), refRatio(iLev - 1),
                                    Geom(iLev - 1), Geom(iLev),
                                    node_status(iLev), node_bilinear_interp);
    }

    // 'conducting' body: B_r = 0, the tangential B carries the surface
    // current. 'linetied' and 'insulating' leave B untouched.
    if (is_body_conducting()) {
      project_body_B(centerB[iLev], iLev);
      project_body_B(nodeB[iLev], iLev);
      nodeB[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }

  if (need_divB() && finestLevel() > 0) {
    for (int iLev = 1; iLev < n_lev(); iLev++) {
      fill_fine_lev_bny_from_coarse(
          divB[iLev - 1], divB[iLev], 0, divB[iLev].nComp(), refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), cell_status(iLev), cell_bilinear_interp);
      fill_fine_lev_bny_from_coarse(
          centerDivB[iLev - 1], centerDivB[iLev], 0, centerDivB[iLev].nComp(),
          refRatio(iLev - 1), Geom(iLev - 1), Geom(iLev), cell_status(iLev),
          cell_bilinear_interp);
    }
  }

  if (need_divB() && n_lev() > 1 && doReport) {
    report_divB_amr();
  }
}

//==========================================================
void Pic::solve_hyp_phi(int iLev) {
  std::string nameFunc = "Pic::solve_hyp_phi";
  timing_func(nameFunc);

  // divB error propagation speed
  Real ch = 0.8 * Geom(iLev).CellSize()[ix_] / tc->get_dt();

  Real coef = -tc->get_dt() * (ch * ch);
  const Real decayFactor = 1.0 - hypDecay;
  for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
    Box box = mfi.validbox();

    const Array4<Real const> divBArr = divB[iLev][mfi].const_array();
    const Array4<Real> phiArr = hypPhi[iLev][mfi].array();

    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      IntVect ijk = { AMREX_D_DECL(i, j, k) };
      phiArr(ijk) = (phiArr(ijk) + coef * divBArr(ijk)) * decayFactor;
    });
  }

  hypPhi[iLev].FillBoundary(Geom(iLev).periodicity());

  if (iLev > 0) {
    fill_fine_lev_bny_from_coarse(
        hypPhi[iLev - 1], hypPhi[iLev], 0, 1, refRatio(iLev - 1),
        Geom(iLev - 1), Geom(iLev), cell_status(iLev), *get_cell_interp());
    const Real scale =
        pow(Geom(iLev).CellSize()[ix_] / Geom(iLev - 1).CellSize()[ix_], 2);
    for (MFIter mfi(hypPhi[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.fabbox();
      const auto phi = hypPhi[iLev][mfi].array();
      const auto st = cell_status(iLev)[mfi].const_array();
      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (bit::is_lev_boundary(st(i, j, k))) {
          phi(i, j, k) *= scale;
        }
      });
    }
  }

  apply_BC(cell_status(iLev), hypPhi[iLev], 0, hypPhi[iLev].nComp(), nullptr,
           iLev, nullptr, iLev > 0);
}

//==========================================================
void Pic::add_curl_upwind_E(MultiFab& nodeEmfLev, int iLev) {
  std::string nameFunc = "Pic::add_curl_upwind_E";
  timing_func(nameFunc);

  if (!useUpwindB || !useCurlUpwindB)
    return;

  GpuArray<int, nDim3> upwindDir = { 0, 0, 0 };
  GpuArray<Real, nDim3> invDxUp = { 0.0, 0.0, 0.0 };
  for (int iDim = 0; iDim < nDim; ++iDim) {
    upwindDir[iDim] = (Geom(iLev).Domain().length(iDim) > 1) ? 1 : 0;
    invDxUp[iDim] = Geom(iLev).InvCellSize(iDim);
  }

  // Map limiterThetaB to limType: 0 = upwind1, 1 = minmod, 3 = mc
  int limTypeUp = 1;
  if (limiterThetaB <= 0.0) {
    limTypeUp = 0;
  } else if (limiterThetaB <= 1.0) {
    limTypeUp = 1;
  } else {
    limTypeUp = 3;
  }

  const Real sensorFloorUp = 1.0;
  const Real fixedSpeedUp = fixedUpwindVel;
  const Real dt = tc->get_dt();
  const bool isFake2D = is_fake_2d();

  for (MFIter mfi(nodeEmfLev); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& arrE = nodeEmfLev[mfi].array();
    const Array4<Real const>& arrBin = centerB[iLev][mfi].const_array();
    const Array4<Real const>& arrU = uBg[iLev][mfi].const_array();

    FArrayBox faceJump[nDim3];
    Array4<Real const> arrFace[nDim3];
    FArrayBox slopeFab;
    for (int iDir = 0; iDir < nDim; ++iDir) {
      if (!upwindDir[iDir])
        continue;
      Box faceBox(box.smallEnd(), box.bigEnd());
      for (int t = 0; t < nDim; ++t) {
        if (t != iDir)
          faceBox.growLo(t, 1);
      }
      Box slopeBox = faceBox;
      slopeBox.growLo(iDir, 1);
      slopeFab.resize(slopeBox, nDim3, The_Async_Arena());
      faceJump[iDir].resize(faceBox, nDim3, The_Async_Arena());
      const Array4<Real>& arrS = slopeFab.array();
      const Array4<Real const>& arrSc = slopeFab.const_array();
      const Array4<Real>& arrF = faceJump[iDir].array();

      if (iDir == ix_) {
        parallel_for_omp(slopeBox, nDim3,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                           arrS(i, j, k, iVar) = upwind_cell_slope<ix_>(
                               arrBin, i, j, k, iVar, limTypeUp);
                         });
        parallel_for_omp(faceBox, nDim3,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                           arrF(i, j, k, iVar) = upwind_face_jump<ix_>(
                               arrBin, arrSc, i, j, k, iVar, sensorFloorUp);
                         });
      } else if (iDir == iy_) {
        parallel_for_omp(slopeBox, nDim3,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                           arrS(i, j, k, iVar) = upwind_cell_slope<iy_>(
                               arrBin, i, j, k, iVar, limTypeUp);
                         });
        parallel_for_omp(faceBox, nDim3,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                           arrF(i, j, k, iVar) = upwind_face_jump<iy_>(
                               arrBin, arrSc, i, j, k, iVar, sensorFloorUp);
                         });
      } else {
        parallel_for_omp(slopeBox, nDim3,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                           arrS(i, j, k, iVar) = upwind_cell_slope<iz_>(
                               arrBin, i, j, k, iVar, limTypeUp);
                         });
        parallel_for_omp(faceBox, nDim3,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
                           arrF(i, j, k, iVar) = upwind_face_jump<iz_>(
                               arrBin, arrSc, i, j, k, iVar, sensorFloorUp);
                         });
      }
      arrFace[iDir] = faceJump[iDir].const_array();
    }

    const Array4<Real const> arrFaceX = arrFace[ix_];
    const Array4<Real const> arrFaceY = arrFace[iy_];
    const Array4<Real const> arrFaceZ = arrFace[iz_];

    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
      Real jump[nDim3];
      for (int iDir = 0; iDir < nDim; ++iDir) {
        if (!upwindDir[iDir])
          continue;
        if (iDir == iz_ && isFake2D)
          continue;

        if (iDir == ix_) {
          upwind_node_jump<ix_>(arrFaceX, i, j, k, jump);
        } else if (iDir == iy_) {
          upwind_node_jump<iy_>(arrFaceY, i, j, k, jump);
        } else {
          upwind_node_jump<iz_>(arrFaceZ, i, j, k, jump);
        }

        Real a = 0.0;
        if (fixedSpeedUp > 0.0) {
          a = fixedSpeedUp;
        } else {
          a = std::abs(arrU(i, j, k, iDir));
        }

        const Real maxA = 1.0 / (dt * invDxUp[iDir]);
        a = amrex::min(a, maxA);
        const Real halfA = 0.5 * a;

        if (iDir == ix_) {
          arrE(i, j, k, iy_) -= halfA * jump[iz_];
          arrE(i, j, k, iz_) += halfA * jump[iy_];
        } else if (iDir == iy_) {
          arrE(i, j, k, ix_) += halfA * jump[iz_];
          arrE(i, j, k, iz_) -= halfA * jump[ix_];
        } else {
          arrE(i, j, k, ix_) -= halfA * jump[iy_];
          arrE(i, j, k, iy_) += halfA * jump[ix_];
        }
      }
    });
  }
}

//==========================================================
void Pic::correct_B(int iLev) {
  std::string nameFunc = "Pic::correct_B";
  timing_func(nameFunc);

  if (!useUpwindB && !useHyperbolicCleaning) {
    return;
  }

  // The upwind diffusion needs uBg, which only the full-PIC solver allocates,
  // so #UPWINDB stays a full-PIC-only option. When curl-form upwind is active,
  // upwind dissipation is already applied inside Faraday's law.
  const bool doUpwind = useUpwindB && !useHybridPIC && !useCurlUpwindB;
  if (!doUpwind && !useHyperbolicCleaning) {
    return;
  }

  if (useHyperbolicCleaning) {
    ensure_hypPhi(iLev);
  }

  if (centerDB[iLev].empty()) {
    distribute_FabArray(centerDB[iLev], boxArray(iLev), DistributionMap(iLev),
                        nDim3, get_n_ghost());
  }
  MultiFab& cDB = centerDB[iLev];
  cDB.setVal(0.0);

  if (doUpwind) {
    const bool isFake2D = is_fake_2d();
    const Real limThetaB = limiterThetaB;
    const Real fixedVel = fixedUpwindVel;
    GpuArray<Real, nDim3> coef = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < nDim; ++i) {
      coef[i] = 0.5 * tc->get_dt() * Geom(iLev).InvCellSize()[i];
    }

    for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
      Box box = mfi.validbox();

      const Array4<Real const> cB = centerB[iLev][mfi].const_array();
      const Array4<Real const> nU = uBg[iLev][mfi].const_array();
      const Array4<Real> dB = cDB[mfi].array();
      const Array4<int const> status = cell_status(iLev)[mfi].const_array();

      // Get the face along the direction iDir for the cell (i,j,k) for the iVar
      // component
      auto get_face = [=] AMREX_GPU_DEVICE(
                          int iDir, int i, int j, int k, int iVar,
                          Array4<Real const> const& arr, Real& l, Real& r) {
        // Generic fixed-upwind-velocity override (e.g. the old TopHat
        // "bypass_limiter" used a constant speed of 1.0). Zero keeps the
        // normal plasma-background-velocity reconstruction below.
        if (fixedVel > 0) {
          l = fixedVel;
          r = fixedVel;
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

      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        Real lu[nDim3] = { 0, 0, 0 }, ru[nDim3] = { 0, 0, 0 };
        Real ul, ur;

        auto is_face_blocked = [&](int ni, int nj, int nk, int, int, int) {
          const int nStat = status(ni, nj, nk);
          if (bit::is_domain_boundary(nStat))
            return true;
          return false;
        };

        // Flux along  x
        for (int iDir = 0; iDir < nDim; iDir++) {
          get_face(iDir, i, j, k, iDir, nU, lu[iDir], ru[iDir]);
        }

        ul = fabs(lu[ix_]);
        ur = fabs(ru[ix_]);
        ul = min(ul, 0.5 / coef[ix_]);
        ur = min(ur, 0.5 / coef[ix_]);
        const bool blockLx = is_face_blocked(i - 1, j, k, i, j, k);
        const bool blockRx = is_face_blocked(i + 1, j, k, i, j, k);

        for (int iVar = 0; iVar < nDim3; iVar++) {
          Real fluxR = 0.0;
          if (!blockRx) {
            Real cR = limiter_theta(limThetaB, cB(i - 1, j, k, iVar),
                                    cB(i, j, k, iVar), cB(i + 1, j, k, iVar));
            fluxR = cR * ur * (cB(i + 1, j, k, iVar) - cB(i, j, k, iVar));
          }
          Real fluxL = 0.0;
          if (!blockLx) {
            Real cL = limiter_theta(limThetaB, cB(i - 2, j, k, iVar),
                                    cB(i - 1, j, k, iVar), cB(i, j, k, iVar));
            fluxL = cL * ul * (cB(i, j, k, iVar) - cB(i - 1, j, k, iVar));
          }
          dB(i, j, k, iVar) += (fluxR - fluxL) * coef[ix_];
        }

        // Flux along y
        ul = fabs(lu[iy_]);
        ur = fabs(ru[iy_]);
        ul = min(ul, 0.5 / coef[iy_]);
        ur = min(ur, 0.5 / coef[iy_]);
        const bool blockLy = is_face_blocked(i, j - 1, k, i, j, k);
        const bool blockRy = is_face_blocked(i, j + 1, k, i, j, k);

        for (int iVar = 0; iVar < nDim3; iVar++) {
          Real fluxR = 0.0;
          if (!blockRy) {
            Real cR = limiter_theta(limThetaB, cB(i, j - 1, k, iVar),
                                    cB(i, j, k, iVar), cB(i, j + 1, k, iVar));
            fluxR = cR * ur * (cB(i, j + 1, k, iVar) - cB(i, j, k, iVar));
          }
          Real fluxL = 0.0;
          if (!blockLy) {
            Real cL = limiter_theta(limThetaB, cB(i, j - 2, k, iVar),
                                    cB(i, j - 1, k, iVar), cB(i, j, k, iVar));
            fluxL = cL * ul * (cB(i, j, k, iVar) - cB(i, j - 1, k, iVar));
          }
          dB(i, j, k, iVar) += (fluxR - fluxL) * coef[iy_];
        }

        if (nDim > 2 && !isFake2D) {
          // Flux along z
          ul = fabs(lu[iz_]);
          ur = fabs(ru[iz_]);
          ul = min(ul, 0.5 / coef[iz_]);
          ur = min(ur, 0.5 / coef[iz_]);
          const bool blockLz = is_face_blocked(i, j, k - 1, i, j, k);
          const bool blockRz = is_face_blocked(i, j, k + 1, i, j, k);

          for (int iVar = 0; iVar < nDim3; iVar++) {
            Real fluxR = 0.0;
            if (!blockRz) {
              Real cR = limiter_theta(limThetaB, cB(i, j, k - 1, iVar),
                                      cB(i, j, k, iVar), cB(i, j, k + 1, iVar));
              fluxR = cR * ur * (cB(i, j, k + 1, iVar) - cB(i, j, k, iVar));
            }
            Real fluxL = 0.0;
            if (!blockLz) {
              Real cL = limiter_theta(limThetaB, cB(i, j, k - 2, iVar),
                                      cB(i, j, k - 1, iVar), cB(i, j, k, iVar));
              fluxL = cL * ul * (cB(i, j, k, iVar) - cB(i, j, k - 1, iVar));
            }
            dB(i, j, k, iVar) += (fluxR - fluxL) * coef[iz_];
          }
        }
      });
    } // end MFIter
  } // end useUpwindB

  if (useHyperbolicCleaning) {
    if (solverCenterLapMF[iLev].empty()) {
      distribute_FabArray(solverCenterLapMF[iLev], boxArray(iLev),
                          DistributionMap(iLev), 3, 1);
      distribute_FabArray(solverTempNode3[iLev], node_box_array(iLev),
                          DistributionMap(iLev), 3, get_n_ghost());
    }
    MultiFab& gradPhi = solverCenterLapMF[iLev];
    solve_hyp_phi(iLev);

    MultiFab& gradPhiNode = solverTempNode3[iLev];
    grad_center_to_node(hypPhi[iLev], gradPhiNode, Geom(iLev).InvCellSize());

    average_node_to_cellcenter(gradPhi, 0, gradPhiNode, 0, nDim3,
                               gradPhi.nGrow());

    const Real dt = tc->get_dt();
    for (MFIter mfi(centerB[iLev]); mfi.isValid(); ++mfi) {
      Box box = mfi.validbox();

      const Array4<Real> dB = cDB[mfi].array();
      const Array4<Real const> gradPhiArr = gradPhi[mfi].const_array();

      parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        IntVect ijk{ AMREX_D_DECL(i, j, k) };
        for (int iVar = 0; iVar < nDim3; iVar++) {
          dB(ijk, iVar) += -dt * gradPhiArr(ijk, iVar);
        }
      });
    } // end MFIter
  } // end useHyperbolicCleaning

  // The div(B) cleaning must not reach into the frozen interior of the body.
  if (is_body_interior_frozen())
    mask_body_interior(cDB, cell_status(iLev));

  MultiFab::Add(centerB[iLev], cDB, 0, 0, nDim3, 0);

  centerB[iLev].FillBoundary(Geom(iLev).periodicity());
}

//==========================================================
void Pic::smooth_multifab(MultiFab& mf, int iLev, int di, Real coef) {
  std::string nameFunc = "Pic::smooth_multifab";
  timing_func(nameFunc);

  if (smoothScratchMF[iLev].empty() ||
      smoothScratchMF[iLev].boxArray() != mf.boxArray() ||
      smoothScratchMF[iLev].DistributionMap() != mf.DistributionMap() ||
      smoothScratchMF[iLev].nComp() != mf.nComp() ||
      smoothScratchMF[iLev].nGrow() != mf.nGrow()) {
    smoothScratchMF[iLev].define(mf.boxArray(), mf.DistributionMap(),
                                 mf.nComp(), mf.nGrow());
  }

  MultiFab& mfOld = smoothScratchMF[iLev];

  auto smooth_dir = [&](int iDir) {
    const int dIdx0 = (iDir == ix_) ? di : 0;
    const int dIdx1 = (iDir == iy_) ? di : 0;
    const int dIdx2 = (iDir == iz_) ? di : 0;

    MultiFab::Copy(mfOld, mf, 0, 0, mf.nComp(), mf.nGrow());

    const bool isPeriodicDir = Geom(iLev).isPeriodic(iDir);
    const BoundaryBounds bnd(Geom(iLev), mf.boxArray().ixType());
    const int domLo = bnd.loBnd[iDir];
    const int domHi = bnd.hiBnd[iDir];
    const Real weightSelf = 1 - coef;
    const Real WeightNei = coef / 2.0;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();

      Array4<Real> const& arrE = mf[mfi].array();
      Array4<Real const> const& arrTmp = mfOld[mfi].const_array();

      if (isPeriodicDir) {
        parallel_for_omp(
            box, mf.nComp(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
              const Real neiLo = arrTmp(i - dIdx0, j - dIdx1, k - dIdx2, iVar);
              const Real neiHi = arrTmp(i + dIdx0, j + dIdx1, k + dIdx2, iVar);
              arrE(i, j, k, iVar) = weightSelf * arrTmp(i, j, k, iVar) +
                                    WeightNei * (neiLo + neiHi);
            });
      } else {
        parallel_for_omp(
            box, mf.nComp(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int iVar) {
              Real neiLo = arrTmp(i - dIdx0, j - dIdx1, k - dIdx2, iVar);
              Real neiHi = arrTmp(i + dIdx0, j + dIdx1, k + dIdx2, iVar);
              const int idx = (iDir == ix_) ? i : ((iDir == iy_) ? j : k);
              if (idx <= domLo) {
                neiLo = neiHi;
              } else if (idx >= domHi) {
                neiHi = neiLo;
              }
              arrE(i, j, k, iVar) = weightSelf * arrTmp(i, j, k, iVar) +
                                    WeightNei * (neiLo + neiHi);
            });
      }
    }

    mf.FillBoundary(Geom(iLev).periodicity());
  };

  smooth_dir(ix_);
  if (nDim > 1)
    smooth_dir(iy_);
  if (nDim > 2 && !is_fake_2d())
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
  if (finestLevel() > 0) {
    for (int iLev = finestLevel(); iLev > 0; iLev--) {
      if (projectScratchMF[iLev].empty() ||
          projectScratchMF[iLev].boxArray() != node_box_array(iLev) ||
          projectScratchMF[iLev].DistributionMap() != DistributionMap(iLev)) {
        distribute_FabArray(projectScratchMF[iLev], node_box_array(iLev),
                            DistributionMap(iLev), 3, 0);
      }
      amrex::MultiFab& tmp = projectScratchMF[iLev];
      tmp.setVal(0.0);
      for (MFIter mfi(tmp); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<Real const> arrE = nodeE[iLev][mfi].const_array();
        const Array4<Real> arrTmp = tmp[mfi].array();
        parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
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
          nodeE[iLev - 1], tmp, 0, nodeE[iLev].nComp(), refRatio(iLev - 1),
          Geom(iLev - 1), Geom(iLev), node_status(iLev), node_bilinear_interp);
      average_down_nodal(tmp, nodeE[iLev - 1], refRatio(iLev - 1));
    }
    for (int iLev = 0; iLev <= finestLevel(); iLev++) {
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
      const Array4<int const> status = cell_status(iLev)[mfi].const_array();
      Box box = mfi.validbox();
      const Array4<Real const> arr = nodeE[iLev][mfi].const_array();

      const Dim3 lo = lbound(box);
      const Dim3 len = length(box);
      const int nNode = nDim > 2 ? 8 : 4;
      const Real sumLoc = deterministic_sum(
          static_cast<int>(box.numPts()) * nNode,
          [=] AMREX_GPU_DEVICE(int l) -> Real {
            const int c = l / nNode, m = l - c * nNode;
            const int k = lo.z + c / (len.x * len.y);
            const int j = lo.y + (c / len.x) % len.y;
            const int i = lo.x + c % len.x;
            if (bit::is_refined(status(i, j, k)))
              return 0.0;
            const int ii = i + (m & 1);
            const int jj = j + ((m >> 1) & 1);
            const int kk = k + (m >> 2);
            return arr(ii, jj, kk, ix_) * arr(ii, jj, kk, ix_) +
                   arr(ii, jj, kk, iy_) * arr(ii, jj, kk, iy_) +
                   arr(ii, jj, kk, iz_) * arr(ii, jj, kk, iz_);
          });

      Real avg = (nDim == 3) ? 0.125 : 0.25;

      sum += sumLoc * 0.5 * avg * get_cell_volume(iLev) / fourPI;
    }
  }
  ParallelDescriptor::ReduceRealSum(sum,
                                    ParallelDescriptor::IOProcessorNumber());

  if (!ParallelDescriptor::IOProcessor())
    sum = 0;

  return sum;
}

//==========================================================
Real Pic::calc_B_field_energy() {
  Real sum = 0;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    // The magnetic energy is that of the *total* field, so the frozen
    // intrinsic field B0 is folded in here when it is configured.
    const bool hasB0 = use_intrinsic_B() && !centerB0[iLev].empty();
    MultiFab* const centerBt =
        hasB0 ? &total_center_B(centerB[iLev], iLev) : nullptr;
    MultiFab& mfB = hasB0 ? *centerBt : centerB[iLev];

    for (MFIter mfi(mfB); mfi.isValid(); ++mfi) {
      const Array4<int const> status = cell_status(iLev)[mfi].const_array();

      const Box& box = mfi.validbox();
      const Array4<Real const> arr = mfB[mfi].const_array();
      const Dim3 lo = lbound(box);
      const Dim3 len = length(box);

      const Real sumLoc = deterministic_sum(
          static_cast<int>(box.numPts()), [=] AMREX_GPU_DEVICE(int c) -> Real {
            const int k = lo.z + c / (len.x * len.y);
            const int j = lo.y + (c / len.x) % len.y;
            const int i = lo.x + c % len.x;
            if (bit::is_refined(status(i, j, k)))
              return 0.0;
            return arr(i, j, k, ix_) * arr(i, j, k, ix_) +
                   arr(i, j, k, iy_) * arr(i, j, k, iy_) +
                   arr(i, j, k, iz_) * arr(i, j, k, iz_);
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
