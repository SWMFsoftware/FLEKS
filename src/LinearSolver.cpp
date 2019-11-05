
#include "linear_solver_wrapper_c.h"
#include "../../srcInterface/multi_ipic3d_domain.h"

void matvec_E_solver(double *vecIn, double *vecOut, int n) {
  MPICs->update_E_matvec(vecIn, vecOut);
}

void matvec_divE_accurate(double *vecIn, double *vecOut, int n) {
  MPICs->divE_accurate_matvec(vecIn, vecOut);
}


