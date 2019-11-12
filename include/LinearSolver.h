#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_

#include "linear_solver_wrapper_c.h"

typedef void (*MATVEC)(double *vecIn, double *vecOut, int n);

void matvec_E_solver(double *vecIn, double *vecOut, int n);
void matvec_divE_accurate(double *vecIn, double *vecOut, int n);


#endif
