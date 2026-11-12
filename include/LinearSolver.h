#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_REAL.H>

#include "GpuBackend.h"

typedef void (*MATVEC)(const amrex::Real *vecIn, amrex::Real *vecOut, int n);

enum KrylovType { GMRES, BICGSTAB, CG };
enum class FieldSolverMode { GMRES, NewtonKrylov };
enum StopType { // (||.|| denotes the 2-norm):
  REL,          // relative stopping crit.:||res|| <= Tol*||res0||
  ABS,          // absolute stopping crit.:||res|| <= Tol
  MAX           // maximum  stopping crit.: max(abs(res)) <= Tol
};

namespace fleks_jfnk {

inline amrex::Real finite_difference_epsilon(amrex::Real normBase,
                                             amrex::Real normDirection) {
  if (normDirection <= 0.0) {
    return 0.0;
  }

  // Scale the finite-difference step with the base state and perturbation
  // norms.
  return std::sqrt(std::numeric_limits<amrex::Real>::epsilon()) *
         (1.0 + normBase) / normDirection;
}

// Vector operations of the Jacobian-free matvec. They are not templates:
// nvcc does not allow device lambdas in a function template instantiated with
// a local lambda type (the nonlinear matvec below).
inline void set_zero(amrex::Real *out, const int n) {
  parallel_for_1d(n, [=] AMREX_GPU_DEVICE(int i) { out[i] = 0.0; });
}

inline void perturb(const amrex::Real *base, const amrex::Real *direction,
                    amrex::Real epsilon, amrex::Real *work, const int n) {
  parallel_for_1d(n, [=] AMREX_GPU_DEVICE(int i) {
    work[i] = base[i] + epsilon * direction[i];
  });
}

inline void difference_quotient(const amrex::Real *baseMatvec,
                                amrex::Real invEpsilon, amrex::Real *out,
                                const int n) {
  parallel_for_1d(n, [=] AMREX_GPU_DEVICE(int i) {
    out[i] = (out[i] - baseMatvec[i]) * invEpsilon;
  });
}

template <typename NonlinearMatvec>
void jacobian_free_matvec(NonlinearMatvec nonlinearMatvec,
                          const amrex::Real *base,
                          const amrex::Real *baseMatvec,
                          const amrex::Real *direction, amrex::Real *out,
                          amrex::Real *work, const int n, const int iLev,
                          const amrex::Real epsilon) {
  if (epsilon == 0.0) {
    set_zero(out, n);
    return;
  }

  // Evaluate (F(base + epsilon * direction) - F(base)) / epsilon.
  perturb(base, direction, epsilon, work, n);

  nonlinearMatvec(work, out, iLev);

  const amrex::Real invEpsilon = 1.0 / epsilon;
  difference_quotient(baseMatvec, invEpsilon, out, n);
}

} // namespace fleks_jfnk

struct GmresWorkspace {
  std::vector<amrex::Real> c;
  std::vector<amrex::Real> s;
  std::vector<amrex::Real> rs;
  std::vector<amrex::Real> krylovBuffer;
  std::vector<amrex::Real> hhBuffer;
};

void matvec_E_solver(const amrex::Real *vecIn, amrex::Real *vecOut, int iLev);
void matvec_divE_accurate(const amrex::Real *vecIn, amrex::Real *vecOut,
                          int iLev);

class LinearSolver {
  int nGrid;
  int nVar;
  int nSolve;
  int nDim;
  amrex::Real tol;
  int nIter;
  MATVEC fMatvec;
  std::vector<amrex::Real> rhsBuffer;
  std::vector<amrex::Real> xLeftBuffer;
  std::vector<amrex::Real> matvecBuffer;
  GmresWorkspace gmresWork;
  int actualMatvec = 0;
  amrex::Real initialResidual = 0.0;
  amrex::Real achievedError = 0.0;

public:
  amrex::Real *rhs;
  amrex::Real *xLeft;
  amrex::Real *matvec;

  LinearSolver()
      : nGrid(0),
        nVar(0),
        nSolve(0),
        nDim(0),
        tol(1),
        nIter(0),
        fMatvec(nullptr),
        rhs(nullptr),
        xLeft(nullptr),
        matvec(nullptr) {}

  LinearSolver(const LinearSolver &) = delete;
  LinearSolver &operator=(const LinearSolver &) = delete;

  void de_alloc() {
    rhsBuffer.clear();
    xLeftBuffer.clear();
    matvecBuffer.clear();
    gmresWork.c.clear();
    gmresWork.s.clear();
    gmresWork.rs.clear();
    gmresWork.krylovBuffer.clear();
    gmresWork.hhBuffer.clear();
    actualMatvec = 0;
    initialResidual = 0.0;
    achievedError = 0.0;
    rhs = nullptr;
    xLeft = nullptr;
    matvec = nullptr;
  }

  void reset(int in) {
    const int nSolveNew = in * nVar;
    if (in != nGrid || nSolveNew != nSolve) {
      nGrid = in;
      de_alloc();
      nSolve = nSolveNew;
      if (nSolve > 0) {
        rhsBuffer.assign(nSolve, 0.0);
        xLeftBuffer.assign(nSolve, 0.0);
        matvecBuffer.assign(nSolve, 0.0);
        rhs = rhsBuffer.data();
        xLeft = xLeftBuffer.data();
        matvec = matvecBuffer.data();
      }
    }

    amrex::Real *const r = rhs;
    amrex::Real *const x = xLeft;
    amrex::Real *const m = matvec;
    parallel_for_1d(nSolve, [=] AMREX_GPU_DEVICE(int i) {
      r[i] = 0;
      x[i] = 0;
      m[i] = 0;
    });
    actualMatvec = 0;
    initialResidual = 0.0;
    achievedError = 0.0;
  }

  void set_tol(amrex::Real in) { tol = in; }
  void set_nIter(int in) { nIter = in; }

  amrex::Real get_tol() const { return tol; }
  int get_nIter() const { return nIter; }

  void set_stats(amrex::Real res0, int nMatvecIn, amrex::Real relErr) {
    initialResidual = res0;
    actualMatvec = nMatvecIn;
    achievedError = relErr;
  }

  void init(int nGridIn, int nVarIn, int nDimIn, MATVEC fIn) {
    nVar = nVarIn;
    nDim = nDimIn;
    fMatvec = fIn;
    reset(nGridIn);
  }

  void solve(int iLev, bool doReport = true) { solve(fMatvec, iLev, doReport); }

  void solve(
      std::function<void(const amrex::Real *, amrex::Real *, const int)> matvec,
      int iLev, bool doReport = true);

  int get_nSolve() const { return nSolve; }
  int get_actual_nIter() const { return actualMatvec; }
  amrex::Real get_initial_residual() const { return initialResidual; }
  amrex::Real get_achieved_error() const { return achievedError; }
};

int gmres(
    std::function<void(const amrex::Real *, amrex::Real *, const int)> matvec,
    int iLev,                       // Func for matrix vector multiplication
    const amrex::Real *rhs,         // Right hand side vector
    amrex::Real *sol,               // Initial guess / solution vector
    const bool isInit,              // true if Sol contains initial guess
    const int n,                    // Number of unknowns
    const int nKrylov,              // Size of krylov subspace
    amrex::Real &tol,               // Required / achieved residual
    const StopType typeStop,        // Determine stopping criterion
    int &nIter,                     // Maximum/actual number of iterations
    const bool doTest,              // Write debug info if true
    MPI_Comm iComm = MPI_COMM_SELF, // MPI communicator
    GmresWorkspace *work = nullptr, amrex::Real *initialResidual = nullptr);

#endif
