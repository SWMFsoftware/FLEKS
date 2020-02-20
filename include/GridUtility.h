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
void distribute_FabArray(amrex::FabArray<FAB>& fa, amrex::BoxArray baNew,
                         const amrex::DistributionMapping& dm, int nComp,
                         int nGst, bool doCopy = true) {
  // Assume 'dm' is the new dm.
  amrex::FabArray<FAB> tmp;

  tmp.define(baNew, dm, nComp, nGst);
  tmp.setVal(-7777);
  if (doCopy && !fa.empty()) {
    tmp.ParallelCopy(fa, 0, 0, fa.nComp(), nGst, nGst);
    // tmp.Redistribute(fa, 0, 0, tmp.nComp(), tmp.nGrowVect());
  } else {
    tmp.setVal(0);
  }
  fa = std::move(tmp);

  // fa.FillBoundary();
}

template <class FAB>
void distribute_FabArray(amrex::FabArray<FAB>& fa, amrex::BoxArray baNew,
                         const amrex::DistributionMapping& dm,
                         bool doCopy = true) {

  distribute_FabArray(fa, baNew, dm, fa.nComp(), fa.nGrow(), doCopy);

  // Assume 'dm' is the new dm.
  // amrex::FabArray<FAB> tmp;
  // tmp.define(baNew, dm, fa.nComp(), fa.nGrow());
  // if (doRedistribute) {
  //   tmp.ParallelCopy(fa, 0, 0, fa.nComp(), 0, 0);
  //   // tmp.Redistribute(fa, 0, 0, tmp.nComp(), tmp.nGrowVect());
  // } else {
  //   tmp.setVal(0);
  // }
  // fa = std::move(tmp);
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
  bl.simplify();
  bl.simplify();
  bl.simplify();
  ba.define(bl);
}

#endif