#ifndef _UINTERP_H_
#define _UINTERP_H_

namespace amrex {

template <class T> class UInterp : public InterpBase {
public:
  virtual ~UInterp() = default;

  // Copied from AMReX_Interpolater.cpp
  Box CoarseBox(const Box& fine, int ratio) override {
    Box b = amrex::coarsen(fine, ratio);
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      if (b.length(i) < 2) {
        b.growHi(i, 1);
      }
    }
    return b;
  }

  // Copied from AMReX_Interpolater.cpp
  Box CoarseBox(const Box& fine, const IntVect& ratio) override {
    Box b = amrex::coarsen(fine, ratio);

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      if (b.length(i) < 2) {
        b.growHi(i, 1);
      }
    }

    return b;
  }

  virtual void interp(const BaseFab<T>& crse, int crse_comp, BaseFab<T>& fine,
                      int fine_comp, int ncomp, const Box& fine_region,
                      const IntVect& ratio, const Geometry& crse_geom,
                      const Geometry& fine_geom, Vector<BCRec> const& bcr,
                      int actual_comp, int actual_state, RunOn gpu_or_cpu) = 0;
};

template <class T> class UNodeBilinear : public UInterp<T> {
public:
  ~UNodeBilinear() override {}

  void interp(const BaseFab<T>& crse, int crse_comp, BaseFab<T>& fine,
              int fine_comp, int ncomp, const Box& fine_region,
              const IntVect& ratio, const Geometry& crse_geom,
              const Geometry& fine_geom, Vector<BCRec> const& bcr,
              int actual_comp, int actual_state, RunOn gpu_or_cpu) override {
    Array4<T const> const& crsearr = crse.const_array();
    Array4<T> const& finearr = fine.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(
        gpu_or_cpu, fine_region, ncomp, i, j, k, n, {
          node_bilinear_interp(i, j, k, n, finearr, fine_comp, crsearr,
                               crse_comp, ratio);
        });
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void node_bilinear_interp(
      int i, int j, int k, int n, Array4<T> const& fine, int fcomp,
      Array4<T const> const& crse, int ccomp, IntVect const& ratio) noexcept {

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
      fine(i, j, k, n + fcomp) =
          rxinv * ryinv * rzinv *
          (crse(ic, jc, kc, n + ccomp) * (ratio[0] - ioff) * (ratio[1] - joff) *
               (ratio[2] - koff) +
           crse(ic + 1, jc, kc, n + ccomp) * (ioff) * (ratio[1] - joff) *
               (ratio[2] - koff) +
           crse(ic, jc + 1, kc, n + ccomp) * (ratio[0] - ioff) * (joff) *
               (ratio[2] - koff) +
           crse(ic + 1, jc + 1, kc, n + ccomp) * (ioff) * (joff) *
               (ratio[2] - koff) +
           crse(ic, jc, kc + 1, n + ccomp) * (ratio[0] - ioff) *
               (ratio[1] - joff) * (koff) +
           crse(ic + 1, jc, kc + 1, n + ccomp) * (ioff) * (ratio[1] - joff) *
               (koff) +
           crse(ic, jc + 1, kc + 1, n + ccomp) * (ratio[0] - ioff) * (joff) *
               (koff) +
           crse(ic + 1, jc + 1, kc + 1, n + ccomp) * (ioff) * (joff) * (koff));
    } else if (joff != 0 && koff != 0) {
      // Node on a Y-Z face
      fine(i, j, k, n + fcomp) =
          ryinv * rzinv *
          (crse(ic, jc, kc, n + ccomp) * (ratio[1] - joff) * (ratio[2] - koff) +
           crse(ic, jc + 1, kc, n + ccomp) * (joff) * (ratio[2] - koff) +
           crse(ic, jc, kc + 1, n + ccomp) * (ratio[1] - joff) * (koff) +
           crse(ic, jc + 1, kc + 1, n + ccomp) * (joff) * (koff));
    } else if (ioff != 0 && koff != 0) {
      // Node on a Z-X face
      fine(i, j, k, n + fcomp) =
          rxinv * rzinv *
          (crse(ic, jc, kc, n + ccomp) * (ratio[0] - ioff) * (ratio[2] - koff) +
           crse(ic + 1, jc, kc, n + ccomp) * (ioff) * (ratio[2] - koff) +
           crse(ic, jc, kc + 1, n + ccomp) * (ratio[0] - ioff) * (koff) +
           crse(ic + 1, jc, kc + 1, n + ccomp) * (ioff) * (koff));
    } else if (ioff != 0 && joff != 0) {
      // Node on a X-Y face
      fine(i, j, k, n + fcomp) =
          rxinv * ryinv *
          (crse(ic, jc, kc, n + ccomp) * (ratio[0] - ioff) * (ratio[1] - joff) +
           crse(ic + 1, jc, kc, n + ccomp) * (ioff) * (ratio[1] - joff) +
           crse(ic, jc + 1, kc, n + ccomp) * (ratio[0] - ioff) * (joff) +
           crse(ic + 1, jc + 1, kc, n + ccomp) * (ioff) * (joff));
    } else if (ioff != 0) {
      // Node on X line
      fine(i, j, k, n + fcomp) =
          rxinv * ((ratio[0] - ioff) * crse(ic, jc, kc, n + ccomp) +
                   (ioff)*crse(ic + 1, jc, kc, n + ccomp));
    } else if (joff != 0) {
      // Node on Y line
      fine(i, j, k, n + fcomp) =
          ryinv * ((ratio[1] - joff) * crse(ic, jc, kc, n + ccomp) +
                   (joff)*crse(ic, jc + 1, kc, n + ccomp));
    } else if (koff != 0) {
      // Node on Z line
      fine(i, j, k, n + fcomp) =
          rzinv * ((ratio[2] - koff) * crse(ic, jc, kc, n + ccomp) +
                   (koff)*crse(ic, jc, kc + 1, n + ccomp));
    } else {
      // Node coincident with coarse node
      fine(i, j, k, n + fcomp) = crse(ic, jc, kc, n + ccomp);
    }
  }
};

} // namespace amrex

#endif