#ifndef DEMO_PARTICLE_KERNELS_H
#define DEMO_PARTICLE_KERNELS_H

#include <cmath>

#include <AMReX_Array4.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace demo_particle {

// Bit zero selects x, bit one y, bit two z. A center shape differs from a
// nodal shape only by its half-cell shift; both use the same CIC convention.
struct Shape {
  int cell[3];
  amrex::Real weight[8];
  amrex::Real derivative[8][3];
};

template <class Particle>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Shape
shape(const Particle& p, const amrex::GpuArray<amrex::Real, 3>& lo,
      const amrex::GpuArray<amrex::Real, 3>& invDx,
      amrex::Real shift = 0) noexcept {
  Shape s;
  amrex::Real f[3];
  for (int d = 0; d < 3; ++d) {
    const amrex::Real u = (p.pos(d) - lo[d]) * invDx[d] - shift;
    s.cell[d] = static_cast<int>(std::floor(u));
    f[d] = u - s.cell[d];
  }
  for (int a = 0; a < 8; ++a) {
    amrex::Real w[3];
    for (int d = 0; d < 3; ++d)
      w[d] = (a & (1 << d)) ? f[d] : 1 - f[d];
    s.weight[a] = w[0] * w[1] * w[2];
    for (int d = 0; d < 3; ++d)
      s.derivative[a][d] = ((a & (1 << d)) ? 1 : -1) * invDx[d] *
                           w[(d + 1) % 3] * w[(d + 2) % 3];
  }
  return s;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void index(const Shape& s, int a,
                                                    int& i, int& j,
                                                    int& k) noexcept {
  i = s.cell[0] + (a & 1);
  j = s.cell[1] + ((a >> 1) & 1);
  k = s.cell[2] + ((a >> 2) & 1);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gather(
    const amrex::Array4<const amrex::Real>& field, const Shape& s,
    amrex::Real (&value)[3]) noexcept {
  for (int c = 0; c < 3; ++c)
    value[c] = 0;
  for (int a = 0; a < 8; ++a) {
    int i, j, k;
    index(s, a, i, j, k);
    for (int c = 0; c < 3; ++c)
      value[c] += s.weight[a] * field(i, j, k, c);
  }
}

// FLEKS Utility.h boris_push_nonrelativistic: alpha =
// (I + cross(omega) + omega outer omega)/(1+omega^2).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void alpha(
    const amrex::Real (&b)[3], amrex::Real qdt2m,
    amrex::Real (&matrix)[9]) noexcept {
  const amrex::Real x = qdt2m * b[0];
  const amrex::Real y = qdt2m * b[1];
  const amrex::Real z = qdt2m * b[2];
  const amrex::Real r = 1 / (1 + x * x + y * y + z * z);
  matrix[0] = (1 + x * x) * r;
  matrix[1] = (z + x * y) * r;
  matrix[2] = (-y + x * z) * r;
  matrix[3] = (-z + x * y) * r;
  matrix[4] = (1 + y * y) * r;
  matrix[5] = (x + y * z) * r;
  matrix[6] = (y + x * z) * r;
  matrix[7] = (-x + y * z) * r;
  matrix[8] = (1 + z * z) * r;
}

template <class Particle>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void push(
    Particle& p, const amrex::Real (&e)[3], const amrex::Real (&b)[3],
    amrex::Real qom, amrex::Real dt, const amrex::GpuArray<amrex::Real, 3>& lo,
    const amrex::GpuArray<amrex::Real, 3>& length, bool averageDrift) noexcept {
  amrex::Real a[9];
  alpha(b, qom * dt / 2, a);
  amrex::Real u[3], v[3];
  for (int d = 0; d < 3; ++d)
    u[d] = p.rdata(d) + qom * dt * e[d] / 2;
  for (int d = 0; d < 3; ++d) {
    v[d] = -p.rdata(d);
    for (int c = 0; c < 3; ++c)
      v[d] += 2 * a[3 * d + c] * u[c];
  }
  for (int d = 0; d < 3; ++d) {
    const amrex::Real drift = averageDrift ? (p.rdata(d) + v[d]) / 2 : v[d];
    p.pos(d) += dt * drift;
    p.pos(d) -= length[d] * std::floor((p.pos(d) - lo[d]) / length[d]);
    p.rdata(d) = v[d];
  }
}

} // namespace demo_particle
#endif
