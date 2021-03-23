#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#include <AMReX.H>

// FLEKS is always 3D. But it can be fake 2D with one cell in the z-direction.
constexpr static int nDim = 3;
constexpr static int ix_ = 0, iy_ = 1, iz_ = 2;

constexpr static int nMoments = 11;
constexpr static int iRho_ = 0;
constexpr static int iUx_ = 1, iUy_ = 2, iUz_ = 3;
constexpr static int iMx_ = 1, iMy_ = 2, iMz_ = 3;
constexpr static int iPxx_ = 4, iPyy_ = 5, iPzz_ = 6, iPxy_ = 7, iPxz_ = 8,
                     iPyz_ = 9, iNum_ = 10;

constexpr static bool doTiling = false;

constexpr static double fourPI = 3.14159265359 * 4;

// Integers to label the the status of a cell.
constexpr static int iBoundary_ = 0, iOnNew_ = 1, iOnOld_ = 2,
                     iAddPTParticle_ = 3;

constexpr static int nPicPartReal = 4;

constexpr static int nPTRecord = 100;
constexpr static int ptRecordSize = 7;
constexpr static int nPTPartReal = nPicPartReal + ptRecordSize * nPTRecord,
                     nPTPartInt = 1;

enum ParticleStaggering { Staggered, NonStaggered };

enum TestCase { RegularSimulation, TwoStream };
#endif
