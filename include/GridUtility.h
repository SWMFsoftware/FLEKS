#ifndef _GRIDUTILITY_H_
#define _GRIDUTILITY_H_

#include <AMReX_DistributionMapping.H>
#include <AMReX_FabArray.H>

#include "Constants.h"
#include "GridInfo.h"

template <class FAB>
void redistribute_FabArray(amrex::FabArray<FAB>& fa,
                           const amrex::DistributionMapping& dm,
                           bool doRedistribute = true) {
  // Assume 'dm' is the new dm.

  amrex::FabArray<FAB> tmp;
  tmp.define(fa.boxArray(), dm, fa.nComp(), fa.nGrow());
  if (doRedistribute) {
    tmp.Redistribute(fa, 0, 0, tmp.nComp(), tmp.nGrowVect());
  } else {
    tmp.setVal(0);
  }
  fa = std::move(tmp);
}

template <class FAB>
void distribute_FabArray(amrex::FabArray<FAB>& fa, amrex::BoxArray baNew,
                         const amrex::DistributionMapping& dm, int nComp,
                         int nGst, bool doCopy = true) {
  // Assume 'dm' is the new dm.
  amrex::FabArray<FAB> tmp;

  if (!baNew.empty()) {
    tmp.define(baNew, dm, nComp, nGst);
    tmp.setVal(-7777);
    if (doCopy && !fa.empty()) {
      tmp.ParallelCopy(fa, 0, 0, fa.nComp(), nGst, nGst);
      // tmp.Redistribute(fa, 0, 0, tmp.nComp(), tmp.nGrowVect());
    } else {
      tmp.setVal(0);
    }
  }

  fa = std::move(tmp);
}

template <class FAB>
void distribute_FabArray(amrex::FabArray<FAB>& fa, amrex::BoxArray baNew,
                         const amrex::DistributionMapping& dm,
                         bool doCopy = true) {

  distribute_FabArray(fa, baNew, dm, fa.nComp(), fa.nGrow(), doCopy);
}

template <class T>
void ChangeMultiFabDM(T& mf, amrex::DistributionMapping dm)

{
  T tmpMF(mf.boxArray(), dm, mf.nComp(), mf.nGrowVect());
  tmpMF.ParallelCopy(mf);
  mf = std::move(tmpMF);
}

// This function is called recursively to combine active patahces into larger
// boxes. The number of boxes should be minimized.
inline void get_boxlist_from_region(amrex::BoxList& bl, GridInfo& gridInfo,
                                    amrex::IntVect imin, amrex::IntVect imax) {

  amrex::IntVect patchSize;
  for (int i = 0; i < nDim; i++) {
    patchSize[i] = gridInfo.get_patch_size(i);
  }

  amrex::IntVect patchLen = imax - imin + 1;

  for (int i = 0; i < nDim; i++) {
    if ((patchLen[i] % patchSize[i]) != 0)
      amrex::Abort("Error: the sub-region range is wrong!");
    patchLen[i] /= patchSize[i];
  }

  if (patchLen == 1) {
    if (gridInfo.get_status(imin[ix_], imin[iy_], imin[iz_]) ==
        GridInfo::iPicOn_) {
      bl.push_back(amrex::Box(imin, imax, { AMREX_D_DECL(0, 0, 0) }));
    }
    return;
  }

  amrex::IntVect dhalf;
  for (int i = 0; i < nDim; i++) {
    dhalf[i] = ceil(patchLen[i] / 2.0) * patchSize[i];
  }

  amrex::BoxList blLoc;
  for (int i = imin[ix_]; i <= imax[ix_]; i += dhalf[ix_])
    for (int j = imin[iy_]; j <= imax[iy_]; j += dhalf[iy_])
      for (int k = imin[iz_]; k <= imax[iz_]; k += dhalf[iz_]) {
        amrex::IntVect iminsub = { AMREX_D_DECL(i, j, k) };
        amrex::IntVect imaxsub = iminsub + dhalf - 1;

        for (int i = 0; i < nDim; i++) {
          if (imaxsub[i] > imax[i])
            imaxsub[i] = imax[i];
        }
        get_boxlist_from_region(blLoc, gridInfo, iminsub, imaxsub);
      }

  // The number '3' is chosen based on numerical experiments.
  for (int i = 0; i < 3; i++)
    blLoc.simplify();

  bl.join(blLoc);

  for (int i = 0; i < 3; i++)
    bl.simplify();
}

#endif
