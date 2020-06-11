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

  if (!dm.empty()) {
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

inline void add_cells_to_BoxArray(amrex::BoxArray& ba,
                                  const amrex::Vector<amrex::IntVect>& cells) {
  amrex::BoxList bl(ba);
  for (const amrex::IntVect& cell : cells) {
    bl.push_back(amrex::Box(cell, cell, bl.ixType()));
  }
  bl.simplify();
  bl.simplify();
  bl.simplify();
  ba.define(bl);
}

inline void add_boxes_to_BoxArray(amrex::BoxArray& ba,
                                  const amrex::Vector<amrex::Box>& boxes) {
  amrex::BoxList bl(ba);
  for (const amrex::Box& bx : boxes) {
    bl.push_back(bx);
  }

  // nSimplify is chosen based on tests. It seems 3 is a reasonable choice.
  const int nSimplify = 3;
  for (int i = 0; i < nSimplify; i++) {
    bl.simplify();
  }

  ba.define(bl);
}

inline void get_boxlist_from_region(amrex::BoxList& bl, GridInfo& gridInfo,
                                    amrex::IntVect imin, amrex::IntVect imax) {

  amrex::IntVect patchSize;
  for (int i = 0; i < 3; i++) {
    patchSize[i] = gridInfo.get_patch_size(i);
  }

  amrex::IntVect patchLen = imax - imin + 1;

  for (int i = 0; i < 3; i++) {
    if ((patchLen[i] % patchSize[i]) != 0)
      amrex::Abort("Error: the sub-region range is wrong!");
    patchLen[i] /= patchSize[i];
  }
  // amrex::Print() << "imin = " << imin << " imax = " << imax
  //                << " patchLen = " << patchLen << " patchenlen == 1"
  //                << ((patchLen == 1) ? 'T' : 'F') << std::endl;

  if (patchLen == 1) {
    if (gridInfo.get_status(imin[ix_], imin[iy_], imin[iz_]) ==
        GridInfo::iPicOn_) {
      bl.push_back(amrex::Box(imin, imax, { 0, 0, 0 }));
    }
    return;
  }

  amrex::IntVect dhalf;
  for (int i = 0; i < 3; i++) {
    dhalf[i] = ceil(patchLen[i] / 2.0) * patchSize[i];
  }

  amrex::BoxList blLoc;
  for (int i = imin[ix_]; i <= imax[ix_]; i += dhalf[ix_])
    for (int j = imin[iy_]; j <= imax[iy_]; j += dhalf[iy_])
      for (int k = imin[iz_]; k <= imax[iz_]; k += dhalf[iz_]) {
        amrex::IntVect iminsub = { i, j, k };
        amrex::IntVect imaxsub = iminsub + dhalf - 1;

        for (int i = 0; i < 3; i++) {
          if (imaxsub[i] > imax[i])
            imaxsub[i] = imax[i];
        }
        get_boxlist_from_region(blLoc, gridInfo, iminsub, imaxsub);
      }
  for (int i = 0; i < 3; i++)
    blLoc.simplify();

  bl.join(blLoc);

  for (int i = 0; i < 3; i++)
    bl.simplify();
}

#endif