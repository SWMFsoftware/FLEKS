#ifndef _SWMFINTERFACE_H_
#define _SWMFINTERFACE_H_

#include <sstream>

#include <mpi.h>

#ifdef _PT_COMPONENT_
extern "C" {
void OH_get_charge_exchange_wrapper(double *rhoIon, double *cs2Ion,
                                    double uIon_D[3], double *rhoNeu,
                                    double *cs2Neu, double uNeu_D[3],
                                    double sourceIon_V[5],
                                    double sourceNeu_V[5]);

void OH_get_charge_exchange_region(int *iRegion, double *r, double *rhoDim,
                                   double *u2Dim, double *uSW2Dim,
                                   double *tempDim, double *tempPu2Dim,
                                   double *mach2, double *machPUI2,
                                   double *machSW2, double *levHP);

void OH_get_solar_wind(double *x, double *y, double *z, double *numDen,
                       double *ur, double *temp, double b[3]);
}

inline void OH_get_charge_exchange_wrapper(float *rhoIon, float *cs2Ion,
                                           float uIon_D[3], float *rhoNeu,
                                           float *cs2Neu, float uNeu_D[3],
                                           float sourceIon_V[5],
                                           float sourceNeu_V[5]) {
  double rhoIonD = *rhoIon, cs2IonD = *cs2Ion;
  double uIonD[3] = { uIon_D[0], uIon_D[1], uIon_D[2] };
  double rhoNeuD = *rhoNeu, cs2NeuD = *cs2Neu;
  double uNeuD[3] = { uNeu_D[0], uNeu_D[1], uNeu_D[2] };
  double sourceIonD[5] = { 0 }, sourceNeuD[5] = { 0 };
  OH_get_charge_exchange_wrapper(&rhoIonD, &cs2IonD, uIonD, &rhoNeuD, &cs2NeuD,
                                 uNeuD, sourceIonD, sourceNeuD);
  for (int i = 0; i < 5; ++i) {
    sourceIon_V[i] = static_cast<float>(sourceIonD[i]);
    sourceNeu_V[i] = static_cast<float>(sourceNeuD[i]);
  }
}

inline void OH_get_charge_exchange_region(int *iRegion, float *r, float *rhoDim,
                                          float *u2Dim, float *uSW2Dim,
                                          float *tempDim, float *tempPu2Dim,
                                          float *mach2, float *machPUI2,
                                          float *machSW2, float *levHP) {
  double rD = *r, rhoDimD = *rhoDim, u2DimD = *u2Dim, uSW2DimD = *uSW2Dim;
  double tempDimD = *tempDim, tempPu2DimD = *tempPu2Dim;
  double mach2D = *mach2, machPUI2D = *machPUI2, machSW2D = *machSW2;
  double levHPD = *levHP;
  OH_get_charge_exchange_region(iRegion, &rD, &rhoDimD, &u2DimD, &uSW2DimD,
                                &tempDimD, &tempPu2DimD, &mach2D, &machPUI2D,
                                &machSW2D, &levHPD);
}

inline void OH_get_solar_wind(float *x, float *y, float *z, float *numDen,
                              float *ur, float *temp, float b[3]) {
  double xD = *x, yD = *y, zD = *z;
  double numDenD = 0, urD = 0, tempD = 0, bD[3] = { 0 };
  OH_get_solar_wind(&xD, &yD, &zD, &numDenD, &urD, &tempD, bD);
  *numDen = static_cast<float>(numDenD);
  *ur = static_cast<float>(urD);
  *temp = static_cast<float>(tempD);
  for (int i = 0; i < 3; ++i)
    b[i] = static_cast<float>(bD[i]);
}
#else
template <typename RealT>
inline void OH_get_charge_exchange_wrapper(RealT *rhoIon, RealT *cs2Ion,
                                           RealT uIon_D[3], RealT *rhoNeu,
                                           RealT *cs2Neu, RealT uNeu_D[3],
                                           RealT sourceIon_V[5],
                                           RealT sourceNeu_V[5]) {}

template <typename RealT>
inline void OH_get_charge_exchange_region(int *iRegion, RealT *r, RealT *rhoDim,
                                          RealT *u2Dim, RealT *uSW2Dim,
                                          RealT *tempDim, RealT *tempPu2Dim,
                                          RealT *mach2, RealT *machPUI2,
                                          RealT *machSW2, RealT *levHP) {}

template <typename RealT>
inline void OH_get_solar_wind(RealT *x, RealT *y, RealT *z, RealT *numDen,
                              RealT *ur, RealT *temp, RealT b[3]) {}

#endif

#endif
