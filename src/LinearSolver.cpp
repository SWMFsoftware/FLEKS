
#include "SWMFDomains.h"
#include "linear_solver_wrapper_c.h"

void matvec_E_solver(double *vecIn, double *vecOut, int n) {
  FLEKSs->pic.update_E_matvec(vecIn, vecOut);
}

void matvec_divE_accurate(double *vecIn, double *vecOut, int n) {
  FLEKSs->pic.divE_accurate_matvec(vecIn, vecOut);
}

void linear_solver_gmres(double tolerance, int nIteration, int nVarSolve,
                         int nDim, int nGrid, double *rhs, double *xLeft,
                         MATVEC fMatvec) {

  int nJ = 1, nK = 1, nBlock = 1;
  MPI_Fint iComm = MPI_Comm_c2f(amrex::ParallelDescriptor::Communicator());
  double precond_matrix_II[1][1];
  precond_matrix_II[0][0] = 0;
  // parameter to choose preconditioner types
  // 0:No precondition; 1: BILU; 2:DILU;
  //[-1,0): MBILU;
  double PrecondParam = 0;
  int lTest = amrex::ParallelDescriptor::MyProc() == 0;
  linear_solver_matvec_c = fMatvec;
  linear_solver_wrapper("GMRES", &tolerance, &nIteration, &nVarSolve, &nDim,
                        &nGrid, &nJ, &nK, &nBlock, &iComm, rhs, xLeft,
                        &PrecondParam, precond_matrix_II[0], &lTest);
}
