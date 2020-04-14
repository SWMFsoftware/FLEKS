#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <AMReX.H>

// FLEKS is always 3D. But it can be fake 2D with one cell in the z-direction.
constexpr static int nDim = 3;
constexpr static int ix_ = 0, iy_ = 1, iz_ = 2;

constexpr static int nMoments = 14;
constexpr static int iRho_ = 0;
constexpr static int iUx_ = 1, iUy_ = 2, iUz_ = 3;
constexpr static int iMx_ = 1, iMy_ = 2, iMz_ = 3;
constexpr static int iPxx_ = 4, iPyy_ = 5, iPzz_ = 6, iPxy_ = 7, iPxz_ = 8,
                     iPyz_ = 9, iJhx_ = 10, iJhy_ = 11, iJhz_ = 12, iNum_ = 13;

constexpr static bool doTiling = false;

constexpr static double fourPI = 3.14159265359 * 4;

// Integers to label the the status of a cell.
constexpr static int iBoundary_ = 0, iOnNew_ = 1, iOnOld_ = 2;
#endif
