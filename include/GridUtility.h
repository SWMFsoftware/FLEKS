#ifndef _GRIDUTILITY_H_
#define _GRIDUTILITY_H_

#include <AMReX_DistributionMapping.H>
#include <AMReX_FabArray.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PhysBCFunct.H>

#include "Bit.h"
#include "Constants.h"
#include "GridInfo.h"
#include "UInterp.h"
#include "Utility.h"

template <class FAB> class PhysBCFunctNoOpFab {
public:
  void operator()(amrex::FabArray<FAB>& /*mf*/, int /*dcomp*/, int /*ncomp*/,
                  amrex::IntVect const& /*nghost*/, amrex::Real /*time*/,
                  int /*bccomp*/) {}
};

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

// Interpolate from coarse lev to fine lev
template <class FAB, class Interp>
void interp_from_coarse_to_fine(amrex::FabArray<FAB>& coarse,
                                amrex::FabArray<FAB>& fine, const int iStart,
                                const int nComp, const amrex::IntVect ratio,
                                const amrex::Geometry& cgeom,
                                const amrex::Geometry& fgeom, Interp* mapper,
                                int nGst = 0) {
  BL_PROFILE("interp_from_coarse_to_fine");
  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  int lo_bc[3] = { amrex::BCType::int_dir, amrex::BCType::int_dir,
                   amrex::BCType::int_dir };
  int hi_bc[3] = { amrex::BCType::int_dir, amrex::BCType::int_dir,
                   amrex::BCType::int_dir };

  amrex::Vector<amrex::BCRec> bcs(nComp, amrex::BCRec(lo_bc, hi_bc));

  PhysBCFunctNoOpFab<FAB> cphysbc, fphysbc;

  // See definition in
  // AMREX/InstallDir/include/AMReX_FillPatchUtil_I.H
  amrex::InterpFromCoarseLevel(f, amrex::IntVect(nGst), 0.0, c, 0, 0, nComp,
                               cgeom, fgeom, cphysbc, 0, fphysbc, 0, ratio,
                               mapper, bcs, 0);
}

// Sum from fine level to coarse level for nodes at the interface of two levels.
template <class FAB>
void sum_fine_to_coarse_lev_bny_node(amrex::FabArray<FAB>& coarse,
                                     amrex::FabArray<FAB>& fine,
                                     const int iStart, const int nComp,
                                     const amrex::IntVect ratio) {
  BL_PROFILE("sum_fine_to_coarse_lev_bny_node");
  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ctmp(c.boxArray(), c.DistributionMap(), nComp, 0);
  ctmp.setVal(0.0);

  // // f and c/ctmp have different distribution maps. Inside
  // average_down_nodal(), parallel copy is called to copy coarsen(f) to ctmp.
  amrex::average_down_nodal(f, ctmp, ratio);

  // This is the Add() declared in AMReX_FabArray.H, and it works for any FAB.
  amrex::Add(c, ctmp, 0, 0, nComp, 0);
}

// Sum from coarse level to fine level for nodes at the boundary of two levels.
template <class FAB>
void fill_fine_lev_bny_from_coarse(amrex::FabArray<FAB>& coarse,
                                   amrex::FabArray<FAB>& fine, const int iStart,
                                   const int nComp, const amrex::IntVect ratio,
                                   const amrex::Geometry& cgeom,
                                   const amrex::Geometry& fgeom,
                                   const amrex::iMultiFab& fstatus,
                                   amrex::Interpolater& mapper) {
  BL_PROFILE("fill_fine_lev_bny_from_coarse");

  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ftmp(f.boxArray(), f.DistributionMap(), nComp,
                            fine.nGrow());
  ftmp.setVal(0.0);

  interp_from_coarse_to_fine(c, ftmp, 0, nComp, ratio, cgeom, fgeom, &mapper,
                             f.nGrow());

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    FAB& fab = f[mfi];
    const auto& box = mfi.fabbox();
    const auto& data = fab.array();

    const auto& statusArr = fstatus[mfi].array();
    const auto& tmp = ftmp[mfi].array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = 0; iVar < f.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_lev_boundary(statusArr(i, j, k))) {
              data(i, j, k, iVar) = tmp(i, j, k, iVar);
            }
          }
  }
}

template <class FAB>
void fill_fine_lev_from_coarse(amrex::FabArray<FAB>& coarse,
                               amrex::FabArray<FAB>& fine, const int iStart,
                               const int nComp, const amrex::IntVect ratio,
                               const amrex::Geometry& cgeom,
                               const amrex::Geometry& fgeom,
                               const amrex::iMultiFab& fstatus,
                               amrex::Interpolater& mapper) {
  BL_PROFILE("fill_fine_lev_from_coarse");

  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ftmp(f.boxArray(), f.DistributionMap(), nComp,
                            fine.nGrow());
  ftmp.setVal(0.0);

  interp_from_coarse_to_fine(c, ftmp, 0, nComp, ratio, cgeom, fgeom, &mapper,
                             f.nGrow());

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    FAB& fab = f[mfi];
    const auto& box = mfi.fabbox();
    const auto& data = fab.array();

    const auto& statusArr = fstatus[mfi].array();
    const auto& tmp = ftmp[mfi].array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = 0; iVar < f.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {

            data(i, j, k, iVar) = tmp(i, j, k, iVar);
          }
  }
  fine.FillBoundary();
}

// Sum from coarse level to fine level for nodes at the boundary of two levels.
template <class FAB>
void sum_coarse_to_fine_lev_bny_node(
    amrex::FabArray<FAB>& coarse, amrex::FabArray<FAB>& fine, const int iStart,
    const int nComp, const amrex::IntVect ratio, const amrex::Geometry& cgeom,
    const amrex::Geometry& fgeom, const amrex::iMultiFab& fstatus) {
  BL_PROFILE("sum_coarse_to_fine_lev_bny_node");

  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ftmp(f.boxArray(), f.DistributionMap(), nComp, 0);
  ftmp.setVal(0.0);
  amrex::UNodeBilinear<typename FAB::value_type> mapper;
  interp_from_coarse_to_fine(c, ftmp, 0, nComp, ratio, cgeom, fgeom, &mapper);

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    FAB& fab = f[mfi];
    const auto& box = mfi.validbox();
    const auto& data = fab.array();

    const auto& statusArr = fstatus[mfi].array();
    const auto& tmp = ftmp[mfi].array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = 0; iVar < f.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_lev_edge(statusArr(i, j, k))) {
              data(i, j, k, iVar) = tmp(i, j, k, iVar);
            }
          }
  }
}

// Sum from coarse level to fine level for domain boundary edge nodes
template <class FAB>
void interp_from_coarse_to_fine_for_domain_edge(
    amrex::FabArray<FAB>& coarse, amrex::FabArray<FAB>& fine, const int iStart,
    const int nComp, const amrex::IntVect ratio, const amrex::Geometry& cgeom,
    const amrex::Geometry& fgeom, const amrex::iMultiFab& fstatus) {
  BL_PROFILE("interp_from_coarse_to_fine_for_domain_edge");

  amrex::FabArray<FAB> f(fine, amrex::make_alias, iStart, nComp);
  amrex::FabArray<FAB> c(coarse, amrex::make_alias, iStart, nComp);

  amrex::FabArray<FAB> ftmp(f.boxArray(), f.DistributionMap(), nComp, 0);
  ftmp.setVal(0.0);

  interp_from_coarse_to_fine(c, ftmp, 0, nComp, ratio, cgeom, fgeom,
                             &amrex::node_bilinear_interp);

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    FAB& fab = f[mfi];
    const auto& box = mfi.validbox();
    const auto& data = fab.array();

    const auto& statusArr = fstatus[mfi].array();
    const auto& tmp = ftmp[mfi].array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int iVar = 0; iVar < f.nComp(); iVar++)
      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i) {
            if (bit::is_domain_edge(statusArr(i, j, k))) {
              data(i, j, k, iVar) = tmp(i, j, k, iVar);
            }
          }
  }
}

template <class FAB>
void sum_two_lev_interface_node(amrex::FabArray<FAB>& coarse,
                                amrex::FabArray<FAB>& fine, int iStart,
                                const int nComp, const amrex::IntVect ratio,
                                const amrex::Geometry& cgeom,
                                const amrex::Geometry& fgeom,
                                const amrex::iMultiFab& fstatus) {
  BL_PROFILE("sum_two_lev_interface_node");
  sum_fine_to_coarse_lev_bny_node(coarse, fine, iStart, nComp, ratio);
  sum_coarse_to_fine_lev_bny_node(coarse, fine, iStart, nComp, ratio, cgeom,
                                  fgeom, fstatus);
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
    if (gridInfo.get_status(AMREX_D_DECL(imin[ix_], imin[iy_], imin[iz_])) ==
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
  AMREX_D_TERM(for (int i = imin[ix_]; i <= imax[ix_]; i += dhalf[ix_]),
               for (int j = imin[iy_]; j <= imax[iy_]; j += dhalf[iy_]),
               for (int k = imin[iz_]; k <= imax[iz_]; k += dhalf[iz_])) {
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
