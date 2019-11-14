#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <AMReX.H>

const int nVirGst = 1;

constexpr static int nDim = AMREX_SPACEDIM, nDimMax = 3;
constexpr static int ix_ = 0, iy_ = 1, iz_ = 2;

constexpr static int nMoments = 13;
constexpr static int iRho_ = 0;
constexpr static int iUx_ = 1, iUy_ = 2, iUz_ = 3;
constexpr static int iMx_ = 1, iMy_ = 2, iMz_ = 3;
constexpr static int iPxx_ = 4, iPyy_ = 5, iPzz_ = 6, iPxy_ = 7, iPxz_ = 8,
                     iPyz_ = 9, iJhx_ = 10, iJhy_ = 11, iJhz_ = 12;

constexpr static bool doTiling=false; 



constexpr static double fourPI = 3.14159265359*4; 

#endif
