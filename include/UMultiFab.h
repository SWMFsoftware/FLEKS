#ifndef _UMULTIFAB_H_
#define _UMULTIFAB_H_

#include "Array1D.h"
#include <AMReX_BLassert.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_Periodicity.H>

namespace amrex {

class UMFInterp_NodeBilinear : public Interpolater {

public:
  // UMFInterp_NodeBilinear () {}
  ~UMFInterp_NodeBilinear() override {}

  Box CoarseBox(const Box& fine, int ratio) override {
    Box b = amrex::coarsen(fine, ratio);

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      if (b.length(i) < 2) {
        b.growHi(i, 1);
      }
    }

    return b;
  }

  Box CoarseBox(const Box& fine, const IntVect& ratio) override {
    Box b = amrex::coarsen(fine, ratio);

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      if (b.length(i) < 2) {
        b.growHi(i, 1);
      }
    }

    return b;
  }

  void interp(const FArrayBox& crse, int crse_comp, FArrayBox& fine,
              int fine_comp, int ncomp, const Box& fine_region,
              const IntVect& ratio, const Geometry& crse_geom,
              const Geometry& fine_geom, Vector<BCRec> const& bcr,
              int actual_comp, int actual_state, RunOn gpu_or_cpu) override {}

  void interp(const BaseFab<RealMM>& crse, int crse_comp, BaseFab<RealMM>& fine,
              int fine_comp, int ncomp, const Box& fine_region,
              const IntVect& ratio, const Geometry& crse_geom,
              const Geometry& fine_geom, Vector<BCRec> const& bcr,
              int actual_comp, int actual_state, RunOn gpu_or_cpu) {
    ncomp = 243;
    Array4<RealMM const> const& crsearr = crse.const_array();
    Array4<RealMM> const& finearr = fine.array();
    RunOn runon = (gpu_or_cpu == RunOn::Gpu) ? RunOn::Gpu : RunOn::Cpu;
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(
        runon, fine_region, ncomp, i, j, k, n, {
          umf_nodebilin_interp(i, j, k, n, finearr, fine_comp, crsearr,
                               crse_comp, ratio);
        });
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void umf_nodebilin_interp(
      int i, int j, int k, int n, Array4<RealMM> const& fine, int fcomp,
      Array4<RealMM const> const& crse, int ccomp,
      IntVect const& ratio) noexcept {

    int ic = amrex::coarsen(i, ratio[0]);
    int jc = amrex::coarsen(j, ratio[1]);
    int kc = amrex::coarsen(k, ratio[2]);
    int ioff = i - ic * ratio[0];
    int joff = j - jc * ratio[1];
    int koff = k - kc * ratio[2];
    Real rxinv = Real(1.0) / Real(ratio[0]);
    Real ryinv = Real(1.0) / Real(ratio[1]);
    Real rzinv = Real(1.0) / Real(ratio[2]);
    if (ioff != 0 && joff != 0 && koff != 0) {
      // Fine node at center of cell
      fine(i, j, k).data[n + fcomp] =
          rxinv * ryinv * rzinv *
          (crse(ic, jc, kc).data[n + ccomp] * (ratio[0] - ioff) *
               (ratio[1] - joff) * (ratio[2] - koff) +
           crse(ic + 1, jc, kc).data[n + ccomp] * (ioff) * (ratio[1] - joff) *
               (ratio[2] - koff) +
           crse(ic, jc + 1, kc).data[n + ccomp] * (ratio[0] - ioff) * (joff) *
               (ratio[2] - koff) +
           crse(ic + 1, jc + 1, kc).data[n + ccomp] * (ioff) * (joff) *
               (ratio[2] - koff) +
           crse(ic, jc, kc + 1).data[n + ccomp] * (ratio[0] - ioff) *
               (ratio[1] - joff) * (koff) +
           crse(ic + 1, jc, kc + 1).data[n + ccomp] * (ioff) *
               (ratio[1] - joff) * (koff) +
           crse(ic, jc + 1, kc + 1).data[n + ccomp] * (ratio[0] - ioff) *
               (joff) * (koff) +
           crse(ic + 1, jc + 1, kc + 1).data[n + ccomp] * (ioff) * (joff) *
               (koff));
    } else if (joff != 0 && koff != 0) {
      // Node on a Y-Z face
      fine(i, j, k).data[n + fcomp] =
          ryinv * rzinv *
          (crse(ic, jc, kc).data[n + ccomp] * (ratio[1] - joff) *
               (ratio[2] - koff) +
           crse(ic, jc + 1, kc).data[n + ccomp] * (joff) * (ratio[2] - koff) +
           crse(ic, jc, kc + 1).data[n + ccomp] * (ratio[1] - joff) * (koff) +
           crse(ic, jc + 1, kc + 1).data[n + ccomp] * (joff) * (koff));
    } else if (ioff != 0 && koff != 0) {
      // Node on a Z-X face
      fine(i, j, k).data[n + fcomp] =
          rxinv * rzinv *
          (crse(ic, jc, kc).data[n + ccomp] * (ratio[0] - ioff) *
               (ratio[2] - koff) +
           crse(ic + 1, jc, kc).data[n + ccomp] * (ioff) * (ratio[2] - koff) +
           crse(ic, jc, kc + 1).data[n + ccomp] * (ratio[0] - ioff) * (koff) +
           crse(ic + 1, jc, kc + 1).data[n + ccomp] * (ioff) * (koff));
    } else if (ioff != 0 && joff != 0) {
      // Node on a X-Y face
      fine(i, j, k).data[n + fcomp] =
          rxinv * ryinv *
          (crse(ic, jc, kc).data[n + ccomp] * (ratio[0] - ioff) *
               (ratio[1] - joff) +
           crse(ic + 1, jc, kc).data[n + ccomp] * (ioff) * (ratio[1] - joff) +
           crse(ic, jc + 1, kc).data[n + ccomp] * (ratio[0] - ioff) * (joff) +
           crse(ic + 1, jc + 1, kc).data[n + ccomp] * (ioff) * (joff));
    } else if (ioff != 0) {
      // Node on X line
      fine(i, j, k).data[n + fcomp] =
          rxinv * ((ratio[0] - ioff) * crse(ic, jc, kc).data[n + ccomp] +
                   (ioff)*crse(ic + 1, jc, kc).data[n + ccomp]);
    } else if (joff != 0) {
      // Node on Y line
      fine(i, j, k).data[n + fcomp] =
          ryinv * ((ratio[1] - joff) * crse(ic, jc, kc).data[n + ccomp] +
                   (joff)*crse(ic, jc + 1, kc).data[n + ccomp]);
    } else if (koff != 0) {
      // Node on Z line
      fine(i, j, k).data[n + fcomp] =
          rzinv * ((ratio[2] - koff) * crse(ic, jc, kc).data[n + ccomp] +
                   (koff)*crse(ic, jc, kc + 1).data[n + ccomp]);
    } else {
      // Node coincident with coarse node
      fine(i, j, k).data[n + fcomp] = crse(ic, jc, kc).data[n + ccomp];
    }
  }
};

template <class T> class UMultiFab : public FabArray<BaseFab<T> > {
public:
  using FabArray<BaseFab<T> >::n_grow;
  using FabArray<BaseFab<T> >::n_comp;
  using FabArray<BaseFab<T> >::boxArray;
  using FabArray<BaseFab<T> >::DistributionMap;
  using FabArray<BaseFab<T> >::Factory;

  UMultiFab() noexcept {}

  void Copy(UMultiFab& dst, const UMultiFab& src, int srccomp, int dstcomp,
            int numcomp, const IntVect& nghost) {
    BL_ASSERT(dst.distributionMap == src.distributionMap);
    BL_ASSERT(dst.nGrowVect().allGE(nghost));

    BL_PROFILE("UMultiFab::Copy()");

    amrex::Copy(dst, src, srccomp, dstcomp, numcomp, nghost);
  }

  void SumBoundary(int scomp, int ncomp, IntVect const& nghost,
                   const Periodicity& period = Periodicity::NonPeriodic()) {
    BL_PROFILE("UMultiFab::SumBoundary()");

    if (n_grow == IntVect::TheZeroVector() &&
        boxArray().ixType().cellCentered())
      return;

    if (boxArray().ixType().cellCentered()) {
      // Self copy is safe only for cell-centered MultiFab
      this->ParallelCopy(*this, scomp, scomp, ncomp, n_grow, nghost, period,
                         FabArrayBase::ADD);
    } else {
      UMultiFab tmp;
      tmp.define(boxArray(), DistributionMap(), ncomp, n_grow, MFInfo(),
                 Factory());
      UMultiFab::Copy(tmp, *this, scomp, 0, ncomp, n_grow);
      this->setVal(0.0, scomp, ncomp, nghost);
      this->ParallelCopy(tmp, 0, scomp, ncomp, n_grow, nghost, period,
                         FabArrayBase::ADD);
    }
  }
  //---------------------------------------------------------

  void SumBoundary(const Periodicity& period = Periodicity::NonPeriodic()) {
    SumBoundary(0, n_comp, IntVect(0), period);
  }
  //---------------------------------------------------------

  void SumBoundary(int scomp, int ncomp,
                   const Periodicity& period = Periodicity::NonPeriodic()) {
    SumBoundary(scomp, ncomp, IntVect(0), period);
  }
  //---------------------------------------------------------
};
} // namespace amrex
#endif
