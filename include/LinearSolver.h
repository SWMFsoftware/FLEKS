#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_

#include <AMReX_ParallelDescriptor.H>

#include "linear_solver_wrapper_c.h"

typedef void (*MATVEC)(double* vecIn, double* vecOut, int n);

void matvec_E_solver(double* vecIn, double* vecOut, int n);
void matvec_divE_accurate(double* vecIn, double* vecOut, int n);

void linear_solver_gmres(double tolerance, int nIteration, int nVarSolve,
                         int nDim, int nGrid, double* rhs, double* xLeft,
                         MATVEC fMatvec);

class LinearSolver {
  int nGrid;
  int nVar;
  int nSolve;
  int nDim;
  double tol;
  int nIter;
  MATVEC fMatvec;

public:
  double* rhs;
  double* xLeft;
  double* matvec;

  LinearSolver()
      : rhs(nullptr),
        xLeft(nullptr),
        matvec(nullptr),
        fMatvec(nullptr),        
        nGrid(0),
        nVar(0),
        nSolve(0),
        nDim(0),
        tol(1),
        nIter(0) {}

  ~LinearSolver() { de_alloc(); }

  void de_alloc() {
    if (rhs != nullptr) {
      delete[] rhs;
      delete[] xLeft;
      delete[] matvec;
      rhs = nullptr;
      xLeft = nullptr;
      fMatvec = nullptr; 
    }
  }

  void reset(int in) {
    if (in != nGrid) {
      nGrid = in;
      de_alloc();
      nSolve = nGrid * nVar;
      if(nSolve>0){
	rhs = new double[nSolve];
	xLeft = new double[nSolve];
	matvec = new double[nSolve];
      }
    }

    for (int i = 0; i < nSolve; i++) {
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

  void solve() {
    linear_solver_gmres(tol, nIter, nVar, nDim, nGrid, rhs, xLeft, fMatvec);
  }

  int get_nSolve() const { return nSolve; }
};

#endif
