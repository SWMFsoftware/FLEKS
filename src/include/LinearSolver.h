#ifndef _LINEARSOLVER_H_
#define _LINEARSOLVER_H_

#include "linear_solver_wrapper_c.h"

typedef void (*MATVEC)(double *, double *, int);
void pic_matvec(double *vecIn, double *vecOut, int n);
// void iPIC3D_PoissonImage(double *vecIn, double *vecOut, int n);
// void iPIC3D_matvec_particle_correction(double *vecIn, double *vecOut, int n);

#endif
