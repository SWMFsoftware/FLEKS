#include "Utility.h"
#include "Constants.h"

using namespace amrex;
void curl_center_to_node(const MultiFab& centerMF, MultiFab& nodeMF,
                         const Real* invDx) {
  Real compZDY, compYDZ;

  Real compXDZ, compZDX;
  Real compYDX, compXDY;

  for (MFIter mfi(nodeMF); mfi.isValid(); ++mfi) // Loop over grids
  {
    const Box& box = mfi.validbox();
    const Array4<Real>& nodeArr = nodeMF[mfi].array();
    const Array4<Real const>& centerArr = centerMF[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {

          compZDY =
              .25 * (centerArr(i, j, k, iz_) - centerArr(i, j - 1, k, iz_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i, j, k - 1, iz_) -
                   centerArr(i, j - 1, k - 1, iz_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i - 1, j, k, iz_) -
                   centerArr(i - 1, j - 1, k, iz_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i - 1, j, k - 1, iz_) -
                   centerArr(i - 1, j - 1, k - 1, iz_)) *
                  invDx[iy_];
          compYDZ =
              .25 * (centerArr(i, j, k, iy_) - centerArr(i, j, k - 1, iy_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i - 1, j, k, iy_) -
                   centerArr(i - 1, j, k - 1, iy_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i, j - 1, k, iy_) -
                   centerArr(i, j - 1, k - 1, iy_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i - 1, j - 1, k, iy_) -
                   centerArr(i - 1, j - 1, k - 1, iy_)) *
                  invDx[iz_];
          // curl - Y
          compXDZ =
              .25 * (centerArr(i, j, k, ix_) - centerArr(i, j, k - 1, ix_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i - 1, j, k, ix_) -
                   centerArr(i - 1, j, k - 1, ix_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i, j - 1, k, ix_) -
                   centerArr(i, j - 1, k - 1, ix_)) *
                  invDx[iz_] +
              .25 *
                  (centerArr(i - 1, j - 1, k, ix_) -
                   centerArr(i - 1, j - 1, k - 1, ix_)) *
                  invDx[iz_];
          compZDX =
              .25 * (centerArr(i, j, k, iz_) - centerArr(i - 1, j, k, iz_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j, k - 1, iz_) -
                   centerArr(i - 1, j, k - 1, iz_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j - 1, k, iz_) -
                   centerArr(i - 1, j - 1, k, iz_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j - 1, k - 1, iz_) -
                   centerArr(i - 1, j - 1, k - 1, iz_)) *
                  invDx[ix_];
          // curl - Z
          compYDX =
              .25 * (centerArr(i, j, k, iy_) - centerArr(i - 1, j, k, iy_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j, k - 1, iy_) -
                   centerArr(i - 1, j, k - 1, iy_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j - 1, k, iy_) -
                   centerArr(i - 1, j - 1, k, iy_)) *
                  invDx[ix_] +
              .25 *
                  (centerArr(i, j - 1, k - 1, iy_) -
                   centerArr(i - 1, j - 1, k - 1, iy_)) *
                  invDx[ix_];
          compXDY =
              .25 * (centerArr(i, j, k, ix_) - centerArr(i, j - 1, k, ix_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i, j, k - 1, ix_) -
                   centerArr(i, j - 1, k - 1, ix_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i - 1, j, k, ix_) -
                   centerArr(i - 1, j - 1, k, ix_)) *
                  invDx[iy_] +
              .25 *
                  (centerArr(i - 1, j, k - 1, ix_) -
                   centerArr(i - 1, j - 1, k - 1, ix_)) *
                  invDx[iy_];

          nodeArr(i, j, k, ix_) = compZDY - compYDZ;
          nodeArr(i, j, k, iy_) = compXDZ - compZDX;
          nodeArr(i, j, k, iz_) = compYDX - compXDY;

          Print() << " i = " << i << " j = " << j << " k = " << k
                  << " nodeX = " << nodeArr(i, j, k, ix_)
                  << " nodeY = " << nodeArr(i, j, k, iy_)
                  << " nodeZ = " << nodeArr(i, j, k, iz_) << std::endl;
        }
  }
}
