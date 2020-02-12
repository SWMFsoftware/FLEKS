#ifndef _GRIDUTILITY_H_
#define _GRIDUTILITY_H_

#include <AMReX_DistributionMapping.H>
#include <AMReX_FabArray.H>

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
void redistribute_FabArray(amrex::FabArray<FAB>& fa, amrex::BoxArray baNew,
                           const amrex::DistributionMapping& dm,
                           bool doRedistribute = true) {
  // Assume 'dm' is the new dm.
  amrex::FabArray<FAB> tmp;
  tmp.define(baNew, dm, fa.nComp(), fa.nGrow());
  if (doRedistribute) {
    tmp.ParallelCopy(fa, 0, 0, fa.nComp(),0, 0);
    // tmp.Redistribute(fa, 0, 0, tmp.nComp(), tmp.nGrowVect());
  } else {
    tmp.setVal(0);
  }
  fa = std::move(tmp);
}

#endif