#ifndef DEMO_FIELD_KERNELS_H
#define DEMO_FIELD_KERNELS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_BaseFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

// Same 27 neighboring nodes x 3x3 tensor layout as FLEKS RealMM. Production
// Arr1D lacks device annotations; keep this GPU-enabled counterpart local.
struct DemoRealMM {
  amrex::Real data[243];
  AMREX_GPU_HOST_DEVICE DemoRealMM(amrex::Real value = 0) noexcept {
    for (int n = 0; n < 243; ++n)
      data[n] = value;
  }
  AMREX_GPU_HOST_DEVICE amrex::Real& operator[](int n) noexcept {
    return data[n];
  }
  AMREX_GPU_HOST_DEVICE const amrex::Real& operator[](int n) const noexcept {
    return data[n];
  }
  AMREX_GPU_HOST_DEVICE DemoRealMM& operator+=(
      const DemoRealMM& other) noexcept {
    for (int n = 0; n < 243; ++n)
      data[n] += other[n];
    return *this;
  }
  AMREX_GPU_HOST_DEVICE DemoRealMM& operator*=(amrex::Real scale) noexcept {
    for (int n = 0; n < 243; ++n)
      data[n] *= scale;
    return *this;
  }
};
static_assert(sizeof(DemoRealMM) == 243 * sizeof(amrex::Real),
              "Mass matrix layout must match production");
using DemoNodeMMFab = amrex::FabArray<amrex::BaseFab<DemoRealMM> >;

namespace demo_field {
constexpr amrex::Real pi = 3.1415926535897932384626433832795;
constexpr amrex::Real fourPi = 4 * pi;
// Production Pic::rhoTheta weights old/new charge in div(E) correction. It is
// fixed independently of the field-solver theta.
constexpr amrex::Real rhoTheta = 0.51;

// FLEKS dual-grid derivative: difference in direction d, averaged over the
// four transverse corners. toCenter selects node->center or center->node.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real derivative(
    const amrex::Array4<const amrex::Real>& a, int i, int j, int k, int c,
    int d, const amrex::GpuArray<amrex::Real, 3>& invDx,
    bool toCenter) noexcept {
  amrex::Real result = 0;
  for (int p = 0; p < 2; ++p)
    for (int q = 0; q < 2; ++q) {
      int lo[3] = { i, j, k };
      int e = (d + 1) % 3, f = (d + 2) % 3;
      lo[e] += toCenter ? p : -p;
      lo[f] += toCenter ? q : -q;
      int hi[3] = { lo[0], lo[1], lo[2] };
      if (toCenter)
        ++hi[d];
      else
        --lo[d];
      result += a(hi[0], hi[1], hi[2], c) - a(lo[0], lo[1], lo[2], c);
    }
  return .25 * invDx[d] * result;
}
// Production div_center_to_center averages central differences over a 3x3
// transverse stencil. This smoothing is part of coefDiff, not optional here.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real center_divergence(
    const amrex::Array4<const amrex::Real>& a, int i, int j, int k,
    const amrex::GpuArray<amrex::Real, 3>& invDx) noexcept {
  amrex::Real result = 0;
  for (int d = 0; d < 3; ++d) {
    amrex::Real difference = 0;
    for (int p = -1; p <= 1; ++p)
      for (int q = -1; q <= 1; ++q) {
        int lo[3] = { i, j, k };
        const int e = (d + 1) % 3, f = (d + 2) % 3;
        lo[e] += p;
        lo[f] += q;
        int hi[3] = { lo[0], lo[1], lo[2] };
        --lo[d];
        ++hi[d];
        difference += a(hi[0], hi[1], hi[2], d) - a(lo[0], lo[1], lo[2], d);
      }
    result += .5 * invDx[d] * difference;
  }
  return result / 9;
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real curl(
    const amrex::Array4<const amrex::Real>& a, int i, int j, int k, int c,
    const amrex::GpuArray<amrex::Real, 3>& invDx, bool toCenter) noexcept {
  const int d = (c + 1) % 3, e = (c + 2) % 3;
  return derivative(a, i, j, k, e, d, invDx, toCenter) -
         derivative(a, i, j, k, d, e, invDx, toCenter);
}
} // namespace demo_field
#endif
