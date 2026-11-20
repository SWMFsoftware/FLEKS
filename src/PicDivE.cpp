#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "LinearSolver.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
void Pic::divE_correction() {
  std::string nameFunc = "Pic::divE_correction";

  timing_func(nameFunc);

#ifdef GMRES_VERBOSE
  const bool useCompactDivELog = false;
#else
  const bool useCompactDivELog = doCompact;
#endif

  struct DivELogEntry {
    Real res0;
    int nMatVec;
    Real relErr;
  };
  std::vector<DivELogEntry> logEntries;
  if (doReport && useCompactDivELog) {
    logEntries.reserve(nDivECorrection);
  }

  for (int iIter = 0; iIter < nDivECorrection; iIter++) {

    sum_to_center(true);

    if (doReport && !useCompactDivELog)
      Print() << "\n-----" << printPrefix << " div(E) correction at iter "
              << iIter << "----------" << std::endl;

    calculate_phi(divESolver, 0, doReport && !useCompactDivELog);

    if (doReport && useCompactDivELog) {
      logEntries.push_back({ divESolver.get_initial_residual(),
                             divESolver.get_actual_nIter(),
                             divESolver.get_achieved_error() });
    }

    divE_correct_particle_position();
  }

  if (doReport && useCompactDivELog && !logEntries.empty()) {
    int totalMatVec = 0;
    for (const auto& entry : logEntries) {
      totalMatVec += entry.nMatVec;
    }
    const Real resInitial = logEntries.front().res0;
    const Real resFinal = logEntries.back().res0 * logEntries.back().relErr;

    if (logEntries.size() == 1) {
      Print() << printPrefix << "div(E): res " << std::setprecision(5)
              << resInitial << " -> " << resFinal << " (" << totalMatVec
              << " mv)" << std::defaultfloat << std::setprecision(6)
              << std::endl;
    } else {
      Print() << printPrefix << "div(E) [" << logEntries.size()
              << " iters]: res " << std::setprecision(5) << resInitial << " -> "
              << resFinal << " (total " << totalMatVec << " mv)"
              << std::defaultfloat << std::setprecision(6) << std::endl;
    }
  }

  for (int i = 0; i < nSpecies; ++i) {
    // The particles outside the simulation domain is marked for deletion
    // inside divE_correct_particle_position(). redistribute_particles()
    // deletes these particles. In order to get correct moments, re-inject
    // particles in the ghost cells.
    parts[i]->redistribute_particles();
  }

  inject_particles_for_boundary_cells();

  sum_to_center(false);
}

//==========================================================
void Pic::divE_correct_particle_position() {
  std::string nameFunc = "Pic::correct_position";

  timing_func(nameFunc);
  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->divE_correct_position(centerPhi, iLev);
    }
  }
}

//==========================================================
void Pic::calculate_phi(LinearSolver& solver, int iLev, bool reportSolver) {
  std::string nameFunc = "Pic::calculate_phi";

  timing_func(nameFunc);

  {
    // Reuse divEOutMF[iLev] (boxArray, dmap, nComp=1, nGrow=0) as residual
    // scratch before solver.solve() uses it inside divE_accurate_matvec.
    MultiFab& residual = divEOutMF[iLev];

    Real coef = 1.0 / rhoTheta;

    MultiFab::LinComb(residual, coef, centerDivE[iLev], 0, -fourPI * coef,
                      centerNetChargeN[iLev], 0, 0, residual.nComp(),
                      residual.nGrow());
    if (finestLevel() > 0) {
      skip_cells_divE_correction(residual, cell_status(iLev), iLev,
                                 bcField.has(FieldBC::coupled));
    }

    // The div(E) equation cannot be satisfied inside a 'linetied' or
    // 'conducting' body, where E is constrained, so those cells are dropped;
    // for an 'insulating' body the charge inside is a real source and is kept.
    if (use_body() && bodyFieldBC != BodyFieldBC::insulating) {
      mask_body(residual, cell_status(iLev));
    }

    if (iLev < static_cast<int>(centerMM.size()) && !centerMM[iLev].empty()) {
#if (AMREX_SPACEDIM == 2)
      constexpr int centerGp = 4;
#else
      constexpr int centerGp = 13;
#endif
      for (MFIter mfi(residual); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const Array4<Real>& resArr = residual[mfi].array();
        const Array4<RealCMM const>& mmArr = centerMM[iLev][mfi].const_array();
        parallel_for_omp(bx,
                         [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                           if (std::abs(mmArr(i, j, k)[centerGp]) == 0.0) {
                             resArr(i, j, k) = 0.0;
                           }
                         });
      }

      // On pure-Neumann refined levels (no Dirichlet coupling to refined,
      // body, or domain-boundary cells), A * 1 == 0 identically, so any
      // nonzero mean residual lies in ker(A) and cannot be reduced by
      // in-level particle shifts. Projecting out the active-cell mean makes
      // the linear system compatible so GMRES converges cleanly.
      if (iLev > 0) {
        const bool maskBody =
            use_body() && bodyFieldBC != BodyFieldBC::insulating;
        const auto& statusMF = cell_status(iLev);
        amrex::ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum>
            reduceOp;
        amrex::ReduceData<Real, Real, Real, Real> reduceData(reduceOp);
        using ReduceTuple = typename decltype(reduceData)::Type;
        for (MFIter mfi(residual); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          const Array4<Real const> resArr = residual[mfi].const_array();
          const Array4<RealCMM const> mmArr = centerMM[iLev][mfi].const_array();
          const Array4<int const> stArr = statusMF[mfi].const_array();
          reduceOp.eval(
              bx, reduceData,
              [=] AMREX_GPU_DEVICE(int i, int j, int k)
                  noexcept -> ReduceTuple {
                    const int st = stArr(i, j, k);
                    if (bit::is_refined(st) || bit::is_lev_boundary(st) ||
                        bit::is_domain_boundary(st) ||
                        (maskBody && bit::is_body(st))) {
                      return { 0.0, 0.0, 0.0, 0.0 };
                    }
                    const Real diag = std::abs(mmArr(i, j, k)[centerGp]);
                    if (diag == 0.0) {
                      return { 0.0, 0.0, 0.0, 0.0 };
                    }
                    Real activeRowSum = 0.0;
#if (AMREX_SPACEDIM == 2)
                    for (int dj = -1; dj <= 1; ++dj) {
                      for (int di = -1; di <= 1; ++di) {
                        const int nst = stArr(i + di, j + dj, 0);
                        if (bit::is_refined(nst) || bit::is_lev_boundary(nst) ||
                            bit::is_domain_boundary(nst) ||
                            (maskBody && bit::is_body(nst))) {
                          continue;
                        }
                        const int gp = (di + 1) * 3 + (dj + 1);
                        activeRowSum += mmArr(i, j, 0)[gp];
                      }
                    }
#else
                    for (int dk = -1; dk <= 1; ++dk) {
                      for (int dj = -1; dj <= 1; ++dj) {
                        for (int di = -1; di <= 1; ++di) {
                          const int nst = stArr(i + di, j + dj, k + dk);
                          if (bit::is_refined(nst) ||
                              bit::is_lev_boundary(nst) ||
                              bit::is_domain_boundary(nst) ||
                              (maskBody && bit::is_body(nst))) {
                            continue;
                          }
                          const int gp = (di + 1) * 9 + (dj + 1) * 3 + (dk + 1);
                          activeRowSum += mmArr(i, j, k)[gp];
                        }
                      }
                    }
#endif
                    return { std::abs(activeRowSum), diag, resArr(i, j, k),
                             Real(1.0) };
                  });
        }
        auto hv = reduceData.value();
        Real vals[4] = { amrex::get<0>(hv), amrex::get<1>(hv),
                         amrex::get<2>(hv), amrex::get<3>(hv) };
        ParallelDescriptor::ReduceRealSum(vals, 4);
        const Real sumActiveRow = vals[0];
        const Real sumDiag = vals[1];
        const Real sumRes = vals[2];
        const Real nActive = vals[3];
        if (nActive > 0.0 && sumDiag > 0.0 &&
            sumActiveRow <= Real(1e-12) * sumDiag) {
          const Real meanRes = sumRes / nActive;
          for (MFIter mfi(residual); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const Array4<Real> resArr = residual[mfi].array();
            const Array4<RealCMM const> mmArr =
                centerMM[iLev][mfi].const_array();
            const Array4<int const> stArr = statusMF[mfi].const_array();
            parallel_for_omp(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  const int st = stArr(i, j, k);
                  if (bit::is_refined(st) || bit::is_lev_boundary(st) ||
                      bit::is_domain_boundary(st) ||
                      (maskBody && bit::is_body(st)) ||
                      std::abs(mmArr(i, j, k)[centerGp]) == 0.0) {
                    return;
                  }
                  resArr(i, j, k) -= meanRes;
                });
          }
        }
      }
    }

    if (fsolver.divESolverBackend == ESolverBackend::AmrexGMRES) {
      if (static_cast<int>(amrexDivESolver.size()) < n_lev_max()) {
        amrexDivESolver.resize(n_lev_max());
      }
      if (!amrexDivESolver[iLev]) {
        amrexDivESolver[iLev] = std::make_unique<AmrexDivESolver>(this, iLev);
      }
      amrexDivESolver[iLev]->define(this, iLev, fsolver.divEPrecondType);

      const Real res0 = std::sqrt(mf_dot(residual, 0, residual, 0, 1, 0,
                                         &amrexDivESolver[iLev]->dot_mask()));
      if (reportSolver) {
        Print() << "Before 0 nMatVec, relative tol: " << solver.get_nIter()
                << " " << solver.get_tol() << std::endl;
        Print() << "initial residual norm: " << res0 << std::endl;
      }

#ifdef GMRES_VERBOSE
      const int verbose = reportSolver ? 2 : 0;
#else
      const int verbose = 0;
#endif

      int nMatvec = 0;
      Real resFinal = 0.0;
      BL_PROFILE_VAR("Pic::phi_iterate", solve);
      amrexDivESolver[iLev]->solve(
          centerPhi[iLev], residual, solver.get_tol(), 0.0, solver.get_nIter(),
          fsolver.gmresRestartLength, verbose, nMatvec, resFinal);
      BL_PROFILE_VAR_STOP(solve);

      const Real relErr = (res0 > 0.0) ? (resFinal / res0) : 0.0;
      solver.set_stats(res0, nMatvec, relErr);

      if (reportSolver) {
        Print() << "After nMatVec, relative error = " << nMatvec << " "
                << relErr << std::endl;
      }
      centerPhi[iLev].FillBoundary(Geom(iLev).periodicity());
    } else {
      solver.reset(get_local_node_or_cell_number(centerDivE[iLev]));

      convert_3d_to_1d(residual, solver.rhs, iLev);

      BL_PROFILE_VAR("Pic::phi_iterate", solve);
      solver.solve(iLev, reportSolver);
      BL_PROFILE_VAR_STOP(solve);

      centerPhi[iLev].setVal(0.0);
      convert_1d_to_3d(solver.xLeft, centerPhi[iLev], iLev);
      centerPhi[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }
}

//==========================================================
void Pic::divE_accurate_matvec_impl(MultiFab& inMF, MultiFab& outMF, int iLev) {
  std::string nameFunc = "Pic::divE_matvec";
  timing_func(nameFunc);

  inMF.FillBoundary(0, 1, IntVect(1), Geom(iLev).periodicity());

  const Real factor = fourPI * fourPI;
#if AMREX_SPACEDIM > 1
  constexpr int jMin = -1, jMax = 1;
#else
  constexpr int jMin = 0, jMax = 0;
#endif
#if AMREX_SPACEDIM > 2
  constexpr int kMin = -1, kMax = 1;
#endif

  for (MFIter mfi(inMF); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();

    const Array4<Real> lArr = outMF[mfi].array();
    const Array4<Real const> rArr = inMF[mfi].const_array();
    const Array4<RealCMM const> mmArr = centerMM[iLev][mfi].const_array();

    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      Real sum = 0.0;
#if (AMREX_SPACEDIM == 2)
      for (int dj = jMin; dj <= jMax; ++dj) {
        for (int di = -1; di <= 1; ++di) {
          const int gp = (di + 1) * 3 + (dj + 1);
          sum += rArr(i + di, j + dj, 0) * mmArr(i, j, 0)[gp];
        }
      }
#else
      for (int dk = kMin; dk <= kMax; ++dk) {
        for (int dj = jMin; dj <= jMax; ++dj) {
          for (int di = -1; di <= 1; ++di) {
            const int gp = (di + 1) * 9 + (dj + 1) * 3 + (dk + 1);
            sum += rArr(i + di, j + dj, k + dk) * mmArr(i, j, k)[gp];
          }
        }
      }
#endif
      lArr(i, j, k) = sum * factor;
    });
  }
}

//==========================================================
void Pic::divE_accurate_matvec(const Real* vecIn, Real* vecOut, int iLev) {
  MultiFab& inMF = divEInMF[iLev];
  MultiFab& outMF = divEOutMF[iLev];

  inMF.setVal(0.0);
  convert_1d_to_3d(vecIn, inMF, iLev);
  divE_accurate_matvec_impl(inMF, outMF, iLev);
  convert_3d_to_1d(outMF, vecOut, iLev);
}

//==========================================================
void Pic::divE_accurate_matvec(const MultiFab& inMF, MultiFab& outMF,
                               int iLev) {
  divEInMF[iLev].setVal(0.0);
  MultiFab::Copy(divEInMF[iLev], inMF, 0, 0, 1, 0);
  divE_accurate_matvec_impl(divEInMF[iLev], outMF, iLev);
}

//==========================================================
void Pic::build_divE_preconditioner(MultiFab& invDiagMF,
                                    const iMultiFab& dotMask, int iLev) const {
  const BoxArray& ba = boxArray(iLev);
  const DistributionMapping& dm = DistributionMap(iLev);
  if (invDiagMF.empty() || invDiagMF.boxArray() != ba ||
      invDiagMF.DistributionMap() != dm) {
    invDiagMF.define(ba, dm, 1, 0);
  }

  const Real factor = fourPI * fourPI;
#if (AMREX_SPACEDIM == 2)
  constexpr int centerGp = 4;
#else
  constexpr int centerGp = 13;
#endif

  for (MFIter mfi(invDiagMF); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real> invArr = invDiagMF[mfi].array();
    const Array4<RealCMM const> mmArr = centerMM[iLev][mfi].const_array();
    const Array4<int const> maskArr = dotMask[mfi].const_array();

    parallel_for_omp(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (maskArr(i, j, k) == 0) {
        invArr(i, j, k) = 0.0;
        return;
      }
      const Real diag = mmArr(i, j, k)[centerGp] * factor;
      invArr(i, j, k) = (std::abs(diag) > 0.0) ? (1.0 / diag) : 0.0;
    });
  }
}

//==========================================================
void Pic::sum_to_center(bool isBeforeCorrection) {
  std::string nameFunc = "Pic::sum_to_center";

  timing_func(nameFunc);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    centerNetChargeNew[iLev].setVal(0.0);

    bool doNetChargeOnly = !isBeforeCorrection;
    if (!doNetChargeOnly) {
      const RealCMM mm0(0.0);
      centerMM[iLev].setVal(mm0);
    }

    for (int i = 0; i < nSpecies; ++i) {
      parts[i]->sum_to_center(centerNetChargeNew[iLev], centerMM[iLev],
                              doNetChargeOnly, iLev);
    }
    if (!doNetChargeOnly) {
      centerMM[iLev].SumBoundary(Geom(iLev).periodicity());
    }

    centerNetChargeNew[iLev].SumBoundary(Geom(iLev).periodicity());

    if (iLev == 0) {
      apply_BC(cell_status(iLev), centerNetChargeNew[iLev], 0,
               centerNetChargeNew[iLev].nComp(), &Pic::get_zero, iLev);
    }

    MultiFab::LinComb(
        centerNetChargeN[iLev], 1 - rhoTheta, centerNetChargeOld[iLev], 0,
        rhoTheta, centerNetChargeNew[iLev], 0, 0,
        centerNetChargeN[iLev].nComp(), centerNetChargeN[iLev].nGrow());

    if (!isBeforeCorrection) {
      MultiFab::Copy(centerNetChargeOld[iLev], centerNetChargeNew[iLev], 0, 0,
                     centerNetChargeOld[iLev].nComp(),
                     centerNetChargeOld[iLev].nGrow());
    }
  }
}

//==========================================================
void Pic::sum_to_center_amr(bool isBeforeCorrection, int iLev) {

  std::string nameFunc = "Pic::sum_to_center";
  timing_func(nameFunc);

  bool doNetChargeOnly = !isBeforeCorrection;

  centerNetChargeNew[iLev].setVal(0.0);
  if (!doNetChargeOnly) {
    const RealCMM mm0(0.0);
    centerMM[iLev].setVal(mm0);
  }

  if (static_cast<int>(sumToCenterJc.size()) != n_lev()) {
    sumToCenterJc.resize(n_lev());
    sumToCenterJf.resize(n_lev());
    sumToCenterCmmc.resize(n_lev());
  }
  MultiFab& jc = sumToCenterJc[iLev];
  MultiFab& jf = sumToCenterJf[iLev];
  CenterMMFab& cmmc = sumToCenterCmmc[iLev];
  const bool hasCoarser = (iLev > 0);
  const bool hasFiner = (iLev < finestLevel());
  if (hasCoarser) {
    if (jc.empty()) {
      BoxArray bac = centerB[iLev - 1].boxArray();
      bac.refine(refRatio(iLev - 1));
      jc.define(bac, centerB[iLev - 1].DistributionMap(), 1, 1);
    }
    jc.setVal(0.0);
    if (!doNetChargeOnly) {
      if (cmmc.empty()) {
        cmmc.define(jc.boxArray(), jc.DistributionMap(), 1, 1);
      }
      cmmc.setVal(RealCMM(0.0));
    }
  }
  if (hasFiner) {
    if (jf.empty()) {
      BoxArray baf = centerB[iLev + 1].boxArray();
      baf.coarsen(refRatio(iLev));
      jf.define(baf, centerB[iLev + 1].DistributionMap(), 1, 1);
    }
    jf.setVal(0.0);
  }
  CenterMMFab* cmmcPtr = (hasCoarser && !doNetChargeOnly) ? &cmmc : nullptr;
  for (int i = 0; i < nSpecies; ++i) {
    parts[i]->sum_to_center_amr(centerNetChargeNew[iLev], jc, jf,
                                centerMM[iLev], cmmcPtr, doNetChargeOnly, iLev);
  }

  if (!doNetChargeOnly) {
    centerMM[iLev].SumBoundary(Geom(iLev).periodicity());
    if (hasCoarser) {
      centerMM[iLev].ParallelAdd(cmmc, 0, 0, 1, IntVect(1), IntVect(0),
                                 Geom(iLev).periodicity());
    }
  }

  centerNetChargeNew[iLev].SumBoundary(Geom(iLev).periodicity());
  if (hasCoarser) {
    centerNetChargeNew[iLev].ParallelAdd(jc, 0, 0, 1, IntVect(1), IntVect(0),
                                         Geom(iLev).periodicity());
  }
  if (hasFiner) {
    centerNetChargeNew[iLev].ParallelAdd(jf, 0, 0, 1, IntVect(1), IntVect(0),
                                         Geom(iLev).periodicity());
  }

  if (iLev == 0) {
    apply_BC(cell_status(iLev), centerNetChargeNew[iLev], 0,
             centerNetChargeNew[iLev].nComp(), &Pic::get_zero, iLev);
  }

  MultiFab::LinComb(
      centerNetChargeN[iLev], 1 - rhoTheta, centerNetChargeOld[iLev], 0,
      rhoTheta, centerNetChargeNew[iLev], 0, 0, centerNetChargeN[iLev].nComp(),
      centerNetChargeN[iLev].nGrow());

  if (!isBeforeCorrection) {
    MultiFab::Copy(centerNetChargeOld[iLev], centerNetChargeNew[iLev], 0, 0,
                   centerNetChargeOld[iLev].nComp(),
                   centerNetChargeOld[iLev].nGrow());
  }
}

//==========================================================
void Pic::amr_divE_correction() {
  std::string nameFunc = "Pic::divE_correction";

  timing_func(nameFunc);

#ifdef GMRES_VERBOSE
  const bool useCompactDivELog = false;
#else
  const bool useCompactDivELog = doCompact;
#endif

  struct DivELogEntry {
    int iter;
    int lev;
    Real res0;
    int nMatVec;
    Real relErr;
  };
  std::vector<DivELogEntry> logEntries;
  if (doReport && useCompactDivELog) {
    logEntries.reserve(nDivECorrection * (finestLevel() + 1));
  }

  for (int iIter = 0; iIter < nDivECorrection; iIter++) {
    for (int iLev = finestLevel(); iLev >= 0; iLev--) {
      sum_to_center_amr(true, iLev);
      skip_cells_divE_correction(centerMM[iLev], cell_status(iLev), iLev,
                                 bcField.has(FieldBC::coupled));

      if (doReport && !useCompactDivELog)
        Print() << "\n-----" << printPrefix << " div(E) correction at iter "
                << iIter << ", lev " << iLev << "----------" << std::endl;

      calculate_phi(divESolver, iLev, doReport && !useCompactDivELog);

      if (doReport && useCompactDivELog) {
        logEntries.push_back({ iIter, iLev, divESolver.get_initial_residual(),
                               divESolver.get_actual_nIter(),
                               divESolver.get_achieved_error() });
      }

      const MultiFab* phiCoarsePtr = nullptr;
      if (iLev > 0) {
        if (static_cast<int>(divEPhiCoarseView.size()) != n_lev()) {
          divEPhiCoarseView.resize(n_lev());
        }
        MultiFab& phiCoarse = divEPhiCoarseView[iLev];
        if (phiCoarse.empty()) {
          phiCoarse.define(sumToCenterJc[iLev].boxArray(),
                           sumToCenterJc[iLev].DistributionMap(), 1, 1);
        }
        phiCoarse.setVal(0.0);
        phiCoarse.ParallelCopy(centerPhi[iLev], 0, 0, 1, IntVect(0), IntVect(1),
                               Geom(iLev).periodicity());
        phiCoarsePtr = &phiCoarse;
      }

      for (int i = 0; i < nSpecies; ++i) {
        parts[i]->divE_correct_position(centerPhi, iLev, phiCoarsePtr);
      }
      if (finestLevel() > 0) {
        for (int i = 0; i < nSpecies; ++i) {
          parts[i]->redistribute_particles();
        }
      }
    }
  }

  if (doReport && useCompactDivELog && !logEntries.empty()) {
    int totalMatVec = 0;
    for (const auto& entry : logEntries) {
      totalMatVec += entry.nMatVec;
    }

    if (finestLevel() == 0) {
      const Real resInitial = logEntries.front().res0;
      const Real resFinal = logEntries.back().res0 * logEntries.back().relErr;
      if (logEntries.size() == 1) {
        Print() << printPrefix << "div(E): res " << std::setprecision(5)
                << resInitial << " -> " << resFinal << " (" << totalMatVec
                << " mv)" << std::defaultfloat << std::setprecision(6)
                << std::endl;
      } else {
        Print() << printPrefix << "div(E) [" << logEntries.size()
                << " iters]: res " << std::setprecision(5) << resInitial
                << " -> " << resFinal << " (total " << totalMatVec << " mv)"
                << std::defaultfloat << std::setprecision(6) << std::endl;
      }
    } else {
      const int nLev = finestLevel() + 1;
      std::vector<Real> levResInit(nLev, 0.0);
      std::vector<Real> levResFinal(nLev, 0.0);
      std::vector<int> levMatVec(nLev, 0);
      std::vector<bool> levSeen(nLev, false);
      for (const auto& entry : logEntries) {
        if (entry.lev < 0 || entry.lev >= nLev)
          continue;
        if (!levSeen[entry.lev]) {
          levResInit[entry.lev] = entry.res0;
          levSeen[entry.lev] = true;
        }
        levResFinal[entry.lev] = entry.res0 * entry.relErr;
        levMatVec[entry.lev] += entry.nMatVec;
      }
      Real worstResInit = 0.0;
      Real worstResFinal = 0.0;
      for (int l = 0; l < nLev; ++l) {
        worstResInit = std::max(worstResInit, levResInit[l]);
        worstResFinal = std::max(worstResFinal, levResFinal[l]);
      }

      std::ostringstream oss;
      oss << printPrefix << "div(E) [" << nDivECorrection << " iters, " << nLev
          << " levs]: res " << std::setprecision(5) << worstResInit << " -> "
          << worstResFinal << " (total " << totalMatVec << " mv) [";
      for (int l = 0; l < nLev; ++l) {
        if (l > 0)
          oss << ", ";
        oss << "L" << l << ": " << levResInit[l] << "->" << levResFinal[l]
            << "/" << levMatVec[l];
      }
      oss << "]";
      Print() << oss.str() << std::endl;
    }
  }

  inject_particles_for_boundary_cells();
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    sum_to_center_amr(false, iLev);
    if (iLev > 0) {
      const bool skipLevEdge0 = bcField.has(FieldBC::coupled);
      skip_cells_divE_correction(centerNetChargeN[iLev], cell_status(iLev),
                                 iLev, skipLevEdge0);
      skip_cells_divE_correction(centerDivE[iLev], cell_status(iLev), iLev,
                                 skipLevEdge0);
    }
  }
}
