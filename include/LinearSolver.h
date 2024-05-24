#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_

#include <AMReX_ParallelDescriptor.H>

typedef void (*MATVEC)(const double *vecIn, double *vecOut, int n);

enum PrecondType {
  NONE,  // No preconditioner
  MBILU, // Gustaffson modification of diagonal blocks
  DILU,  // LU for diagonal, keep off-diagonal blocks
  BILU,  // LU for diagonal, premultiply U with D^-1
};
enum KrylovType { GMRES, BICGSTAB, CG };
enum StopType { // (||.|| denotes the 2-norm):
  REL,          // relative stopping crit.:||res|| <= Tol*||res0||
  ABS,          // absolute stopping crit.:||res|| <= Tol
  MAX           // maximum  stopping crit.: max(abs(res)) <= Tol
};
enum PrecondSideType { LEFT };

struct LinearSolverParam {
  bool doPrecond;                  // Do preconditioning
  PrecondSideType typePrecondSide; // Precondition left, right, symmetric
  PrecondType typePrecond;         // Preconditioner type
  KrylovType typeKrylov;           // Krylov solver type
  StopType typeStop;               // Stopping criterion type
  double errorMax;                 // Tolerance for solver
  int maxMatvec;                   // Maximum number of iterations
  int nKrylovVector;               // Number of vectors for GMRES
  bool useInitialGuess;            // Non-zero initial guess
  double error;                    // Actual accuracy achieved
  int nMatvec;                     // Actual number of iterations
};

struct Block {
  int nDim;       // Number of spatial dimensions
  int nVar;       // Number of impl. variables/cell
  int nI, nJ, nK; // Number of cells in a block
  int nBlock;     // Number of impl. blocks on current proc
  double *precondMatrix_II;
};

void matvec_E_solver(const double *vecIn, double *vecOut, int iLev);
void matvec_divE_accurate(const double *vecIn, double *vecOut, int iLev);

void linear_solver_gmres(double tolerance, int nIteration, int nVarSolve,
                         int nDim, int nGrid, double *rhs, double *xLeft,
                         MATVEC fMatvec, int iLev, bool doReport = true);

// hyzhou: eventually we should use this and merge the above into this class!
class LinearSolver {
  int nGrid;
  int nVar;
  int nSolve;
  int nDim;
  double tol;
  int nIter;
  MATVEC fMatvec;

public:
  double *rhs;
  double *xLeft;
  double *matvec;

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

  ~LinearSolver() { de_alloc(); }

  void de_alloc() {
    if (rhs != nullptr) {
      delete[] rhs;
      delete[] xLeft;
      delete[] matvec;
      rhs = nullptr;
      xLeft = nullptr;
      matvec = nullptr;
    }
  }

  void reset(int in) {
    if (in != nGrid) {
      nGrid = in;
      de_alloc();
      nSolve = nGrid * nVar;
      if (nSolve > 0) {
        rhs = new double[nSolve];
        xLeft = new double[nSolve];
        matvec = new double[nSolve];
      }
    }

    for (int i = 0; i < nSolve; ++i) {
      rhs[i] = 0;
      xLeft[i] = 0;
      matvec[i] = 0;
    }
  }

  void set_tol(amrex::Real in) { tol = in; }
  void set_nIter(int in) { nIter = in; }

  void init(int nGridIn, int nVarIn, int nDimIn, MATVEC fIn) {
    nVar = nVarIn;
    nDim = nDimIn;
    fMatvec = fIn;
    reset(nGridIn);
  }

  void solve(int iLev, bool doReport = true) {
    linear_solver_gmres(tol, nIter, nVar, nDim, nGrid, rhs, xLeft, fMatvec,
                        iLev, doReport);
  }

  int get_nSolve() const { return nSolve; }
};

void linear_solver_wrapper_hy(
    std::function<void(const double *, double *, const int)> matvec, int iLev,
    const KrylovType solverType, const double tolerance, const int nIteration,
    const int nVar, const int nDim, const int nI, const int nJ, const int nK,
    const int nBlock, MPI_Comm iComm, double *Rhs_I, double *x_I,
    const PrecondType TypePrecond, double *precond_matrix, const int lTest);

int gmres(std::function<void(const double *, double *, const int)> matvec,
          int iLev,                // Func for matrix vector multiplication
          const double *rhs,       // Right hand side vector
          double *sol,             // Initial guess / solution vector
          const bool isInit,       // true if Sol contains initial guess
          const int n,             // Number of unknowns
          const int nKrylov,       // Size of krylov subspace
          double &tol,             // Required / achieved residual
          const StopType typeStop, // Determine stopping criterion
          int &nIter,              // Maximum/actual number of iterations
          const bool doTest,       // Write debug info if true
          MPI_Comm iComm);         // MPI communicator

#endif
