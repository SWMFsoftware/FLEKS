#include <AMReX_BC_TYPES.H>
#include <AMReX_PhysBCFunct.H>
#include <Constants.h>

#include "BC.h"

using namespace amrex;

void apply_float_boundary(const iMultiFab& status, MultiFab& mf,
                          const Geometry& gm, const int iStart, const int nComp,
                          const int nshift) {

  Abort("Error: function apply_float_boundary has not been implemented!");

  //----The commented implementation is problematic!!------------------

  // if (gm.isAllPeriodic())
  //   return;
  // if (mf.nGrow() == 0)
  //   return;

  // BoxArray ba = mf.boxArray();
  // for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
  //   const Box& bx = mfi.fabbox();

  //   //! if there are cells not in the valid + periodic grown box
  //   //! we need to fill them here
  //   if (!ba.contains(bx)) {
  //     const auto lo = IntVect(bx.loVect());
  //     const auto hi = IntVect(bx.hiVect());

  //     IntVect mid = (lo + hi) / 2;

  //     Vector<BCRec> bcr(1);
  //     for (int iDim = 0; iDim < nDim; iDim++) {
  //       auto idxLo = mid;
  //       idxLo[iDim] = lo[iDim];
  //       bcr[0].setLo(iDim,
  //                    ba.contains(idxLo) ? BCType::int_dir :
  //                    BCType::foextrap);

  //       auto idxHi = mid;
  //       idxHi[iDim] = hi[iDim];
  //       bcr[0].setHi(iDim,
  //                    ba.contains(idxHi) ? BCType::int_dir :
  //                    BCType::foextrap);
  //     }

  //     Array4<Real> const& arr = mf[mfi].array();
  //     const auto& statusArr = status[mfi].array();

  //     // Include ghost cells.
  //     int iMin = lo[ix_], iMax = hi[ix_];
  //     int jMin = lo[iy_], jMax = hi[iy_];
  //     int kMin = lo[iz_], kMax = hi[iz_];

  //     IntVect nGst = mf.nGrowVect();

  //     // x left
  //     // if (bcr[0].lo(ix_) == BCType::foextrap)
  //     for (int iVar = iStart; iVar < nComp; iVar++)
  //       for (int k = kMin; k <= kMax; ++k)
  //         for (int j = jMin; j <= jMax; ++j)
  //           for (int i = iMin; i <= iMin + nGst[ix_] - 1 + nshift; ++i) {
  //             if (statusArr(i, j, k) == iBoundary_)
  //               arr(i, j, k, iVar) = arr(iMin + nGst[ix_] + nshift, j, k,
  //               iVar);
  //           }

  //     // x right
  //     // if (bcr[0].hi(ix_) == BCType::foextrap)
  //     for (int iVar = iStart; iVar < nComp; iVar++)
  //       for (int k = kMin; k <= kMax; ++k)
  //         for (int j = jMin; j <= jMax; ++j)
  //           for (int i = iMax - nGst[ix_] + 1 - nshift; i <= iMax; ++i) {
  //             if (statusArr(i, j, k) == iBoundary_)
  //               arr(i, j, k, iVar) = arr(iMax - nGst[ix_] - nshift, j, k,
  //               iVar);
  //           }

  //     // y left
  //     // if (bcr[0].lo(iy_) == BCType::foextrap)
  //     for (int iVar = iStart; iVar < nComp; iVar++)
  //       for (int k = kMin; k <= kMax; ++k)
  //         for (int j = jMin; j <= jMin + nGst[iy_] - 1 + nshift; ++j)
  //           for (int i = iMin; i <= iMax; ++i) {
  //             if (statusArr(i, j, k) == iBoundary_)
  //               arr(i, j, k, iVar) = arr(i, jMin + nGst[iy_] + nshift, k,
  //               iVar);
  //           }

  //     // y right
  //     // if (bcr[0].hi(iy_) == BCType::foextrap)
  //     for (int iVar = iStart; iVar < nComp; iVar++)
  //       for (int k = kMin; k <= kMax; ++k)
  //         for (int j = jMax - nGst[iy_] + 1 - nshift; j <= jMax; ++j)
  //           for (int i = iMin; i <= iMax; ++i) {
  //             if (statusArr(i, j, k) == iBoundary_)
  //               arr(i, j, k, iVar) = arr(i, jMax - nGst[iy_] - nshift, k,
  //               iVar);
  //           }

  //     // z left
  //     // if (bcr[0].lo(iz_) == BCType::foextrap)
  //     for (int iVar = iStart; iVar < nComp; iVar++)
  //       for (int k = kMin; k <= kMin + nGst[iz_] - 1 + nshift; ++k)
  //         for (int j = jMin; j <= jMax; ++j)
  //           for (int i = iMin; i <= iMax; ++i) {
  //             if (statusArr(i, j, k) == iBoundary_)
  //               arr(i, j, k, iVar) = arr(i, j, kMin + nGst[iz_] + nshift,
  //               iVar);
  //           }

  //     // z right
  //     // if (bcr[0].hi(iz_) == BCType::foextrap)
  //     for (int iVar = iStart; iVar < nComp; iVar++)
  //       for (int k = kMax - nGst[iz_] + 1 - nshift; k <= kMax; ++k)
  //         for (int j = jMin; j <= jMax; ++j)
  //           for (int i = iMin; i <= iMax; ++i) {
  //             if (statusArr(i, j, k) == iBoundary_)
  //               arr(i, j, k, iVar) = arr(i, j, kMax - nGst[iz_] - nshift,
  //               iVar);
  //           }
  //   }
  // }
}
