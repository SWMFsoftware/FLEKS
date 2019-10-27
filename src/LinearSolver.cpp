
#include "linear_solver_wrapper_c.h"
#include "../../srcInterface/multi_ipic3d_domain.h"

void pic_matvec(double *vecIn, double *vecOut, int n) {
  MPICs->update_E_matvec(vecIn, vecOut);
}
