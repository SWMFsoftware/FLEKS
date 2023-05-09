#ifndef _SWMFINTERFACE_H_
#define _SWMFINTERFACE_H_

#include <mpi.h>
#include <sstream>

#ifdef _PT_COMPONENT_
extern "C" {
void OH_get_charge_exchange_wrapper(double *rhoIon, double *cs2Ion,
                                    double uIon_D[3], double *rhoNeu,
                                    double *cs2Neu, double uNeu_D[3],
                                    double sourceIon_V[5],
                                    double sourceNeu_V[5]);

void OH_get_charge_exchange_region(int *iRegion, double *r, double *rhoDim,
                                   double *u2Dim, double *tempDim,
                                   double *mach2);
}
#else
void OH_get_charge_exchange_wrapper(double *rhoIon, double *cs2Ion,
                                    double uIon_D[3], double *rhoNeu,
                                    double *cs2Neu, double uNeu_D[3],
                                    double sourceIon_V[5],
                                    double sourceNeu_V[5]) {}

void OH_get_charge_exchange_region(int *iRegion, double *r, double *rhoDim,
                                   double *u2Dim, double *tempDim,
                                   double *mach2) {}
#endif

#endif
