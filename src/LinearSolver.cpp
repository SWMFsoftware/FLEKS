#include "LinearSolver.h"
#include "SimDomains.h"

using namespace amrex;

void matvec_E_solver(const Real *vecIn, Real *vecOut, int iLev) {
  fleksDomains(fleksDomains.selected())
      .pic->update_E_matvec(vecIn, vecOut, iLev);
}

void matvec_divE_accurate(const Real *vecIn, Real *vecOut, int iLev) {
  fleksDomains(fleksDomains.selected())
      .pic->divE_accurate_matvec(vecIn, vecOut, iLev);
}

void LinearSolver::solve(
    std::function<void(const Real *, Real *, const int)> matvec, int iLev,
    bool doReport) {
  const bool doTest = (doReport && ParallelDescriptor::MyProc() == 0);
  const MPI_Comm iComm = ParallelDescriptor::Communicator();
  const StopType typeStop = REL;
  const bool useInitialGuess = false;

  // Initialize solution vector to zero
  Real *const x = xLeft;
  parallel_for_1d(nSolve, [=] AMREX_GPU_DEVICE(int i) { x[i] = 0.0; });

  // Initialize stopping conditions. Solver will return actual values.
  int nMatvec = nIter;
  Real error = tol;

  // Label the tolerance/error explicitly as relative or absolute
  const char *typeError = typeStop == REL ? "relative" : "absolute";

  if (doTest)
    std::cout << "Before " << KrylovType::GMRES << " nMatVec, " << typeError
              << " tol: " << nMatvec << " " << error << std::endl;

  gmres(matvec, iLev, rhs, xLeft, useInitialGuess, nSolve, nIter, error,
        typeStop, nMatvec, doTest, iComm, &gmresWork, &initialResidual);

  if (doTest)
    std::cout << "After nMatVec, " << typeError << " error = " << nMatvec << " "
              << error << std::endl;

  actualMatvec = nMatvec;
  achievedError = error;
}

/*
Initially written by Youcef Saad (May 23, 1985)
Revised by Henk A. van der Vorst and Mike Botchev (Oct 1996)
Rewritten into F90 and parallelized by Gabor Toth (May 2002)
Moved into LinearSolver.f90 for SWMF by Gabor Toth (Dec 2006)
Rewritten into C++ by Hongyang Zhou (Oct 2020)
*/
int gmres(std::function<void(const Real *, Real *, const int)> matvec,
          int iLev,                // Func for matrix std::vector multiplication
          const Real *rhs,         // Right hand side std::vector
          Real *sol,               // Initial guess / solution std::vector
          const bool isInit,       // true if sol contains initial guess
          const int n,             // Number of unknowns
          const int nKrylov,       // Size of krylov subspace
          Real &tol,               // Required / achieved norm of residual
          const StopType typeStop, // Determine stopping criterion
          int &nIter,              // Maximum/actual number of iterations
          const bool doTest,       // Write debug info if true
          MPI_Comm iComm,          // MPI communicator
          GmresWorkspace *work, Real *initialResidual) {
  int info = -1;
  // gives reason for returning:
  //    abs(info)=  0: solution found satisfying given tolerance.
  //                2: no convergence within maximum number of iterations.
  //                3: initial guess satisfies the stopping criterion.
  //   sign(info)=  +: residual decreased
  //                -: residual did not reduce

#ifndef GMRES_VERBOSE
  // Concise mode: only the caller's summary lines are printed.
  (void)doTest; // Avoid the unused-parameter warning.
#endif

#ifdef GMRES_VERBOSE
  if (doTest)
    std::cout << "GMRES tol,nIter:" << tol << " " << nIter << std::endl;
#endif

  int nKrylov1 = nKrylov + 1;
  const size_t krylovSize = static_cast<size_t>(n) * (nKrylov + 2);
  const size_t hhSize = static_cast<size_t>(nKrylov1) * nKrylov;

  std::vector<Real> cFallback;
  std::vector<Real> sFallback;
  std::vector<Real> rsFallback;
  std::vector<Real> krylovFallback;
  std::vector<Real> hhFallback;

  std::vector<Real> &c = work ? work->c : cFallback;
  std::vector<Real> &s = work ? work->s : sFallback;
  std::vector<Real> &rs = work ? work->rs : rsFallback;
  std::vector<Real> &krylovBuffer = work ? work->krylovBuffer : krylovFallback;
  std::vector<Real> &hhBuffer = work ? work->hhBuffer : hhFallback;

  if (c.size() < static_cast<size_t>(nKrylov))
    c.resize(nKrylov);
  if (s.size() < static_cast<size_t>(nKrylov))
    s.resize(nKrylov);
  if (rs.size() < static_cast<size_t>(nKrylov1))
    rs.resize(nKrylov1);
  if (krylovBuffer.size() < krylovSize)
    krylovBuffer.resize(krylovSize);
  if (hhBuffer.size() < hhSize)
    hhBuffer.resize(hhSize);

  auto *Krylov_II = krylovBuffer.data();
  auto *hh = hhBuffer.data();

  constexpr Real epsmac = (sizeof(Real) == 8) ? Real(1e-16) : Real(1e-8);

  Real Tol1, ro, ro0;
  int its = 0;

  do // Restart loop
  {
    // Compute initial residual std::vector
    // Krylov_II[1]:=A*sol
    if (isInit || its > 0) {
      matvec(sol, Krylov_II, iLev);
      parallel_for_1d(n, [=] AMREX_GPU_DEVICE(int i) {
        Krylov_II[i] = rhs[i] - Krylov_II[i];
      });
    } else {
      // Save a matvec when starting from zero initial condition
      parallel_for_1d(n,
                      [=] AMREX_GPU_DEVICE(int i) { Krylov_II[i] = rhs[i]; });
    }
    //-------------------------------------------------------------
    ro = sqrt(dot_product_mpi(Krylov_II, Krylov_II, n, iComm));
    if (ro == 0.0) {
      if (its == 0) {
        info = 3;
      } else {
        info = 0;
      }

      tol = ro;
      nIter = its;
      if (initialResidual)
        *initialResidual = 0.0;
      return info;
    }

    // Set Tol1 for stopping criterion
    if (its == 0) {
      ro0 = ro;
      if (initialResidual)
        *initialResidual = ro0;
      if (doTest)
        std::cout << "initial residual norm: " << ro0 << std::endl;
      if (typeStop == ABS) {
        Tol1 = tol;
        if (ro <= Tol1) { // Quit if accurate enough
          info = 3;
          tol = ro;
          nIter = its;
#ifdef GMRES_VERBOSE
          if (doTest)
            std::cout << "GMRES: nothing to do. info = " << info;
#endif
          return info;
        }
      } else {
        Tol1 = tol * ro;
      }
    }

    auto coef = 1.0 / ro;
    parallel_for_1d(n, [=] AMREX_GPU_DEVICE(int i) { Krylov_II[i] *= coef; });

    // Initialize 1st term of RHS of Hessenberg system
    rs[0] = ro;
    int i = 0;
    do // KRYLOVLOOP
    {
      int i1;
      its += 1;
      i1 = i + 1;

      // Krylov_II[i1]:=A*Krylov_II[i]
      matvec(&Krylov_II[i * n], &Krylov_II[i1 * n], iLev);

      // Modified Gram-Schmidt
      Real *AMREX_RESTRICT v_i1 = &Krylov_II[i1 * n];
      for (int j = 0; j <= i; ++j) {
        const Real *AMREX_RESTRICT v_j = &Krylov_II[j * n];
        Real t = dot_product_mpi(v_j, v_i1, n, iComm);
        hh[i * nKrylov1 + j] = t;
        Real *const vNew = &Krylov_II[i1 * n];
        const Real *const vOld = &Krylov_II[j * n];
        parallel_for_1d(
            n, [=] AMREX_GPU_DEVICE(int k) { vNew[k] -= t * vOld[k]; });
      }
      Real cDot = sqrt(dot_product_mpi(v_i1, v_i1, n, iComm));
      hh[i * nKrylov1 + i1] = cDot;
      if (cDot != 0.0) {
        cDot = 1.0 / cDot;
        Real *const vNew = &Krylov_II[i1 * n];
        parallel_for_1d(n, [=] AMREX_GPU_DEVICE(int k) { vNew[k] *= cDot; });
      }
      // Done with modified Gram-Schmidt and Arnoldi step.
      // Update factorization of hh.
      // Perform previous transformations on i-th column of h.
      for (int k = 1; k <= i; ++k) {
        int k1 = k - 1;
        Real t = hh[i * nKrylov1 + k1];
        hh[i * nKrylov1 + k1] = c[k1] * t + s[k1] * hh[i * nKrylov1 + k];
        hh[i * nKrylov1 + k] = -s[k1] * t + c[k1] * hh[i * nKrylov1 + k];
      }
      Real g = sqrt(hh[i * nKrylov1 + i] * hh[i * nKrylov1 + i] +
                    hh[i * nKrylov1 + i1] * hh[i * nKrylov1 + i1]);
      if (g == 0.0)
        g = epsmac;
      // Determine next plane rotation
      c[i] = hh[i * nKrylov1 + i] / g;
      s[i] = hh[i * nKrylov1 + i1] / g;
      rs[i1] = -s[i] * rs[i];
      rs[i] *= c[i];
      // Determine residual norm and test for convergence
      hh[i * nKrylov1 + i] =
          c[i] * hh[i * nKrylov1 + i] + s[i] * hh[i * nKrylov1 + i1];
      ro = std::abs(rs[i1]);
#ifdef GMRES_VERBOSE
      if (doTest) {
        if (typeStop == REL) {
          std::cout << its << " matvecs, "
                    << " ||rn||/||r0|| = " << ro / ro0 << std::endl;
        } else if (typeStop == ABS) {
          std::cout << its << " matvecs, "
                    << " ||rn|| = " << ro << std::endl;
        }
      }
#endif
      i += 1;
    } while (i < nKrylov && (ro > Tol1));

    // Compute solution. First solve upper triangular system.
    // rs := hh(1:i,1:i) ^-1 * rs
    for (int j = i - 1; j >= 0; --j) {
      if (rs[j] != 0.0) {
        rs[j] /= hh[j * nKrylov1 + j];
        for (int k = j - 1; k >= 0; k--) {
          rs[k] -= rs[j] * hh[j * nKrylov1 + k];
        }
      }
    }

    // Done with back substitution.
    // Form linear combination to get solution.
    for (int j = 0; j < i; ++j) {
      const Real rsj = rs[j];
      const Real *const v = &Krylov_II[j * n];
      parallel_for_1d(n, [=] AMREX_GPU_DEVICE(int k) { sol[k] += rsj * v[k]; });
    }
  } while (ro > Tol1 && its < nIter);

  // eps == tolerance for stopping criterion.
  // Process is stopped as soon as ( ||.|| is the euclidean norm):
  // || current residual||/||initial residual|| <= eps
  // on OUTPUT: actual achieved norm residual (if iabs!=0)
  // or achieved relative residual norm reduction.

  nIter = its;
  tol = tol / Tol1 * ro; // (relative) tolerance achieved

  if (ro < Tol1) {
    info = 0;
  } else if (ro < ro0) {
    info = 2;
  } else {
    info = -2;
  }

  return info;
}
