#ifndef _INTRINSICBFIELD_H_
#define _INTRINSICBFIELD_H_

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <AMReX_GpuContainers.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

// GPU-safe POD geometry descriptors for external current primitives (#COIL and
// #TOROIDALFIELD), all stored in code units after convert_units().
struct CoilLoopParam {
  amrex::Real center[3] = { 0.0, 0.0, 0.0 };
  amrex::Real axis[3] = { 0.0, 0.0, 1.0 };
  amrex::Real radius = 0.0;
  amrex::Real coreRadius = 0.0;
  // Nominal (mu0 * N * I / 4*pi) converted to code_B * code_L.
  amrex::Real kI0 = 0.0;
  int groupIdx = 0;
};

struct WireSegmentParam {
  amrex::Real x1[3] = { 0.0, 0.0, 0.0 };
  amrex::Real dir[3] = { 1.0, 0.0, 0.0 };
  amrex::Real length = 0.0;
  amrex::Real coreRadius = 0.0;
  amrex::Real kI0 = 0.0;
  int groupIdx = 0;
};

struct Filament2DParam {
  amrex::Real x0 = 0.0;
  amrex::Real y0 = 0.0;
  amrex::Real coreRadius = 0.0;
  amrex::Real rRef2 = 1.0;
  amrex::Real kI0 = 0.0;
  int groupIdx = 0;
};

struct ToroidalFieldParam {
  amrex::Real center[3] = { 0.0, 0.0, 0.0 };
  amrex::Real axis[3] = { 0.0, 0.0, 1.0 };
  amrex::Real rRef = 1.0;
  amrex::Real rCore = 0.0;
  // Nominal (B_ref * R_ref) in code_B * code_L.
  amrex::Real kTor0 = 0.0;
  int groupIdx = 0;
};

// Intrinsic / external magnetic field B0 (planetary #DIPOLE, #CRUSTALFIELD,
// external #COIL, and #TOROIDALFIELD) plus optional inductive electric field
// E0(t) = -dA0/dt when #COILWAVEFORM drives time-varying coil currents.
// B0 and E0 are background fields excluded from Faraday evolution of B1 and
// div(B) cleaning; total fields are B = B1 + B0 and (in non-hybrid solvers)
// E = E1 + E0, while hybrid Ohm's law yields total E directly and subtracts E0
// in Faraday's law for B1.
struct B0Params {
  // Center of both the dipole and the spherical harmonic expansion, code
  // units.
  amrex::Real center[3] = { 0.0, 0.0, 0.0 };

  // Dipole moment m in code units: B = [3 (m.r^) r^ - m] / r^3, so
  // |m| = Bc * rRef^3 for the equatorial strength Bc at the reference radius.
  amrex::Real moment[3] = { 0.0, 0.0, 0.0 };
  bool useDipole = false;

  // Crustal coefficients, flat layout coef[n * nMax + m], already rescaled to
  // code units. Null when the crustal field is off.
  const amrex::Real *coefC = nullptr;
  const amrex::Real *coefD = nullptr;

  // Number of harmonic degrees (n = 0..nMax - 1); 0 disables the crustal part.
  int nMax = 0;

  // External coil and toroidal field arrays (managed device memory).
  const CoilLoopParam *loops = nullptr;
  int nLoops = 0;

  const WireSegmentParam *segments = nullptr;
  int nSegments = 0;

  const Filament2DParam *filaments = nullptr;
  int nFilaments = 0;

  const ToroidalFieldParam *toroidals = nullptr;
  int nToroidals = 0;

  // Per-group time-dependent factors at the active evaluation time:
  // groupScaleB[g] = f_g(t), groupScaleE[g] = -(df_g/dt_code)(t).
  const amrex::Real *groupScaleB = nullptr;
  const amrex::Real *groupScaleE = nullptr;
  int nGroups = 0;
};

// Dipole part of B0 at the code-unit position (x, y, z), added into b[].
// Device callable and scratch free, so it may run inside a ParallelFor.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_dipole_b(const B0Params &p,
                                                            amrex::Real x,
                                                            amrex::Real y,
                                                            amrex::Real z,
                                                            amrex::Real b[3]) {
  const amrex::Real xc = x - p.center[0];
  const amrex::Real yc = y - p.center[1];
  const amrex::Real zc = z - p.center[2];

  const amrex::Real r2 = xc * xc + yc * yc + zc * zc;
  if (r2 <= amrex::Real(1.0e-12))
    return;

  const amrex::Real term =
      3.0 * (p.moment[0] * xc + p.moment[1] * yc + p.moment[2] * zc) / r2;
  const amrex::Real invR3 = 1.0 / (r2 * std::sqrt(r2));

  b[0] += (term * xc - p.moment[0]) * invR3;
  b[1] += (term * yc - p.moment[1]) * invR3;
  b[2] += (term * zc - p.moment[2]) * invR3;
}

// Complete elliptic integrals K(m) and E(m) via the arithmetic-geometric mean.
// Parameter convention: m = k^2 in [0, 1).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_elliptic_ke(
    amrex::Real m, amrex::Real &kVal, amrex::Real &eVal) {
  constexpr amrex::Real cHalfPi = 1.57079632679489661923;
  const amrex::Real mClamped =
      (m < amrex::Real(0.0))
          ? amrex::Real(0.0)
          : ((m > amrex::Real(1.0 - 1.0e-15)) ? amrex::Real(1.0 - 1.0e-15) : m);

  amrex::Real a = 1.0;
  amrex::Real b = std::sqrt(amrex::Real(1.0) - mClamped);
  amrex::Real sum = mClamped;
  amrex::Real pow2 = 1.0;

  for (int iter = 0; iter < 16; ++iter) {
    const amrex::Real c = amrex::Real(0.5) * (a - b);
    const amrex::Real aNext = amrex::Real(0.5) * (a + b);
    b = std::sqrt(a * b);
    a = aNext;
    pow2 *= amrex::Real(2.0);
    sum += pow2 * c * c;
    if (std::abs(c) <= amrex::Real(1.0e-15) * a)
      break;
  }

  kVal = cHalfPi / a;
  eVal = kVal * (amrex::Real(1.0) - amrex::Real(0.5) * sum);
}

// Circular current loop: adds B0 * scaleB into b[] and A0 * scaleE into e[].
// Uses Plummer axial softening Z_eff^2 = zeta^2 + coreRadius^2, preserving
// curl(A0) == B0 identically even for coreRadius > 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_coil_loop_ba(
    const CoilLoopParam &c, amrex::Real scaleB, amrex::Real scaleE,
    amrex::Real x, amrex::Real y, amrex::Real z, amrex::Real b[3],
    amrex::Real e[3]) {
  if (scaleB == amrex::Real(0.0) && scaleE == amrex::Real(0.0))
    return;

  constexpr amrex::Real cPi = 3.14159265358979323846;

  const amrex::Real dx = x - c.center[0];
  const amrex::Real dy = y - c.center[1];
  const amrex::Real dz = z - c.center[2];

  const amrex::Real zeta = dx * c.axis[0] + dy * c.axis[1] + dz * c.axis[2];
  const amrex::Real rx = dx - zeta * c.axis[0];
  const amrex::Real ry = dy - zeta * c.axis[1];
  const amrex::Real rz = dz - zeta * c.axis[2];

  const amrex::Real rho2 = rx * rx + ry * ry + rz * rz;
  const amrex::Real rho = std::sqrt(rho2);
  const amrex::Real r2 = c.radius * c.radius;
  const amrex::Real zEff2 = zeta * zeta + c.coreRadius * c.coreRadius;
  const amrex::Real eps2 = (c.coreRadius > amrex::Real(0.0))
                               ? amrex::Real(0.0)
                               : (amrex::Real(1.0e-24) * r2);

  const amrex::Real dMinus = c.radius - rho;
  const amrex::Real dPlus = c.radius + rho;
  const amrex::Real alpha2 = dMinus * dMinus + zEff2 + eps2;
  const amrex::Real beta2 = dPlus * dPlus + zEff2 + eps2;
  const amrex::Real beta = std::sqrt(beta2);

  amrex::Real m = (amrex::Real(4.0) * c.radius * rho) / beta2;
  if (m > amrex::Real(1.0 - 1.0e-15))
    m = amrex::Real(1.0 - 1.0e-15);

  amrex::Real K = 0.0;
  amrex::Real E = 0.0;
  eval_elliptic_ke(m, K, E);

  if (scaleB != amrex::Real(0.0)) {
    const amrex::Real kB = c.kI0 * scaleB;
    const amrex::Real bZeta = (amrex::Real(2.0) * kB / beta) *
                              (K + ((r2 - rho2 - zEff2) / alpha2) * E);
    amrex::Real bRhoOverRho = 0.0;
    if (m < amrex::Real(1.0e-3)) {
      const amrex::Real beta5 = beta2 * beta2 * beta;
      bRhoOverRho =
          (amrex::Real(8.0) * cPi * kB * r2 * zeta / beta5) *
          (amrex::Real(0.375) +
           m * (amrex::Real(0.46875) + m * (amrex::Real(525.0 / 1024.0) +
                                            m * amrex::Real(2205.0 / 4096.0))));
    } else {
      bRhoOverRho = (amrex::Real(2.0) * kB * zeta / (rho2 * beta)) *
                    (-K + ((r2 + rho2 + zEff2) / alpha2) * E);
    }

    b[0] += bZeta * c.axis[0] + bRhoOverRho * rx;
    b[1] += bZeta * c.axis[1] + bRhoOverRho * ry;
    b[2] += bZeta * c.axis[2] + bRhoOverRho * rz;
  }

  if (scaleE != amrex::Real(0.0)) {
    const amrex::Real kE = c.kI0 * scaleE;
    amrex::Real aPhiOverRho = 0.0;
    if (m < amrex::Real(1.0e-3)) {
      const amrex::Real beta3 = beta2 * beta;
      aPhiOverRho =
          (amrex::Real(8.0) * cPi * kE * r2 / beta3) *
          (amrex::Real(0.125) +
           m * (amrex::Real(0.09375) + m * (amrex::Real(75.0 / 1024.0) +
                                            m * amrex::Real(245.0 / 4096.0))));
    } else {
      aPhiOverRho = (kE * beta / rho2) *
                    ((amrex::Real(2.0) - m) * K - amrex::Real(2.0) * E);
    }

    const amrex::Real cx = c.axis[1] * rz - c.axis[2] * ry;
    const amrex::Real cy = c.axis[2] * rx - c.axis[0] * rz;
    const amrex::Real cz = c.axis[0] * ry - c.axis[1] * rx;
    e[0] += aPhiOverRho * cx;
    e[1] += aPhiOverRho * cy;
    e[2] += aPhiOverRho * cz;
  }
}

// Finite straight wire segment x1 -> x2 with core regularization:
// curl(A0) == B0 identically everywhere.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_wire_segment_ba(
    const WireSegmentParam &s, amrex::Real scaleB, amrex::Real scaleE,
    amrex::Real x, amrex::Real y, amrex::Real z, amrex::Real b[3],
    amrex::Real e[3]) {
  if (scaleB == amrex::Real(0.0) && scaleE == amrex::Real(0.0))
    return;

  const amrex::Real rx1 = x - s.x1[0];
  const amrex::Real ry1 = y - s.x1[1];
  const amrex::Real rz1 = z - s.x1[2];

  const amrex::Real s1 = rx1 * s.dir[0] + ry1 * s.dir[1] + rz1 * s.dir[2];
  const amrex::Real s2 = s1 - s.length;
  const amrex::Real px = rx1 - s1 * s.dir[0];
  const amrex::Real py = ry1 - s1 * s.dir[1];
  const amrex::Real pz = rz1 - s1 * s.dir[2];

  const amrex::Real eps2 = amrex::Real(1.0e-24) * s.length * s.length;
  const amrex::Real rhoEff2 =
      px * px + py * py + pz * pz + s.coreRadius * s.coreRadius + eps2;
  const amrex::Real r1 = std::sqrt(s1 * s1 + rhoEff2);
  const amrex::Real r2 = std::sqrt(s2 * s2 + rhoEff2);

  if (scaleB != amrex::Real(0.0)) {
    const amrex::Real facB = (s.kI0 * scaleB / rhoEff2) * (s1 / r1 - s2 / r2);
    const amrex::Real cx = s.dir[1] * pz - s.dir[2] * py;
    const amrex::Real cy = s.dir[2] * px - s.dir[0] * pz;
    const amrex::Real cz = s.dir[0] * py - s.dir[1] * px;
    b[0] += facB * cx;
    b[1] += facB * cy;
    b[2] += facB * cz;
  }

  if (scaleE != amrex::Real(0.0)) {
    const amrex::Real rhoEff = std::sqrt(rhoEff2);
    const amrex::Real facE =
        (s.kI0 * scaleE) * (std::asinh(s1 / rhoEff) - std::asinh(s2 / rhoEff));
    e[0] += facE * s.dir[0];
    e[1] += facE * s.dir[1];
    e[2] += facE * s.dir[2];
  }
}

// Infinite straight wire along +z at (x0, y0) in 2D / fake-2D runs.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_filament2d_ba(
    const Filament2DParam &f, amrex::Real scaleB, amrex::Real scaleE,
    amrex::Real x, amrex::Real y, amrex::Real b[3], amrex::Real e[3]) {
  if (scaleB == amrex::Real(0.0) && scaleE == amrex::Real(0.0))
    return;

  const amrex::Real dx = x - f.x0;
  const amrex::Real dy = y - f.y0;
  const amrex::Real rEff2 = dx * dx + dy * dy + f.coreRadius * f.coreRadius +
                            amrex::Real(1.0e-24) * f.rRef2;

  if (scaleB != amrex::Real(0.0)) {
    const amrex::Real facB = amrex::Real(2.0) * f.kI0 * scaleB / rEff2;
    b[0] += -facB * dy;
    b[1] += facB * dx;
  }

  if (scaleE != amrex::Real(0.0)) {
    e[2] += -(f.kI0 * scaleE) * std::log(rEff2 / f.rRef2);
  }
}

// Vacuum 1/R_perp toroidal magnetic field around (center, axis), regularized
// with core radius rCore so div(B) == 0 identically everywhere.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_toroidal_ba(
    const ToroidalFieldParam &t, amrex::Real scaleB, amrex::Real scaleE,
    amrex::Real x, amrex::Real y, amrex::Real z, amrex::Real b[3],
    amrex::Real e[3]) {
  if (scaleB == amrex::Real(0.0) && scaleE == amrex::Real(0.0))
    return;

  const amrex::Real dx = x - t.center[0];
  const amrex::Real dy = y - t.center[1];
  const amrex::Real dz = z - t.center[2];

  const amrex::Real zeta = dx * t.axis[0] + dy * t.axis[1] + dz * t.axis[2];
  const amrex::Real rx = dx - zeta * t.axis[0];
  const amrex::Real ry = dy - zeta * t.axis[1];
  const amrex::Real rz = dz - zeta * t.axis[2];

  const amrex::Real eps2 = amrex::Real(1.0e-24) * t.rRef * t.rRef;
  const amrex::Real rhoEff2 =
      rx * rx + ry * ry + rz * rz + t.rCore * t.rCore + eps2;

  if (scaleB != amrex::Real(0.0)) {
    const amrex::Real facB = (t.kTor0 * scaleB) / rhoEff2;
    const amrex::Real cx = t.axis[1] * rz - t.axis[2] * ry;
    const amrex::Real cy = t.axis[2] * rx - t.axis[0] * rz;
    const amrex::Real cz = t.axis[0] * ry - t.axis[1] * rx;
    b[0] += facB * cx;
    b[1] += facB * cy;
    b[2] += facB * cz;
  }

  if (scaleE != amrex::Real(0.0)) {
    const amrex::Real facE = -amrex::Real(0.5) * (t.kTor0 * scaleE) *
                             std::log(rhoEff2 / (t.rRef * t.rRef));
    e[0] += facE * t.axis[0];
    e[1] += facE * t.axis[1];
    e[2] += facE * t.axis[2];
  }
}

// Add all GPU-callable intrinsic B0 contributions (dipole + loops + segments +
// 2D filaments + toroidal fields) at (x, y, z) into b[].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_intrinsic_b_dev(
    const B0Params &p, amrex::Real x, amrex::Real y, amrex::Real z,
    int activeDim, amrex::Real b[3]) {
  if (p.useDipole) {
    const amrex::Real zDip = (activeDim > 2) ? z : p.center[2];
    eval_dipole_b(p, x, y, zDip, b);
  }

  amrex::Real dummyE[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < p.nLoops; ++i) {
    const int g = p.loops[i].groupIdx;
    const amrex::Real sB = p.groupScaleB ? p.groupScaleB[g] : amrex::Real(1.0);
    eval_coil_loop_ba(p.loops[i], sB, amrex::Real(0.0), x, y, z, b, dummyE);
  }
  for (int i = 0; i < p.nSegments; ++i) {
    const int g = p.segments[i].groupIdx;
    const amrex::Real sB = p.groupScaleB ? p.groupScaleB[g] : amrex::Real(1.0);
    eval_wire_segment_ba(p.segments[i], sB, amrex::Real(0.0), x, y, z, b,
                         dummyE);
  }
  for (int i = 0; i < p.nFilaments; ++i) {
    const int g = p.filaments[i].groupIdx;
    const amrex::Real sB = p.groupScaleB ? p.groupScaleB[g] : amrex::Real(1.0);
    eval_filament2d_ba(p.filaments[i], sB, amrex::Real(0.0), x, y, b, dummyE);
  }
  for (int i = 0; i < p.nToroidals; ++i) {
    const int g = p.toroidals[i].groupIdx;
    const amrex::Real sB = p.groupScaleB ? p.groupScaleB[g] : amrex::Real(1.0);
    eval_toroidal_ba(p.toroidals[i], sB, amrex::Real(0.0), x, y, z, b, dummyE);
  }
}

// Add all GPU-callable inductive E0(t) = -dA0/dt contributions at (x, y, z)
// into e[].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_inductive_e_dev(
    const B0Params &p, amrex::Real x, amrex::Real y, amrex::Real z,
    amrex::Real e[3]) {
  amrex::Real dummyB[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < p.nLoops; ++i) {
    const int g = p.loops[i].groupIdx;
    const amrex::Real sE = p.groupScaleE ? p.groupScaleE[g] : amrex::Real(0.0);
    eval_coil_loop_ba(p.loops[i], amrex::Real(0.0), sE, x, y, z, dummyB, e);
  }
  for (int i = 0; i < p.nSegments; ++i) {
    const int g = p.segments[i].groupIdx;
    const amrex::Real sE = p.groupScaleE ? p.groupScaleE[g] : amrex::Real(0.0);
    eval_wire_segment_ba(p.segments[i], amrex::Real(0.0), sE, x, y, z, dummyB,
                         e);
  }
  for (int i = 0; i < p.nFilaments; ++i) {
    const int g = p.filaments[i].groupIdx;
    const amrex::Real sE = p.groupScaleE ? p.groupScaleE[g] : amrex::Real(0.0);
    eval_filament2d_ba(p.filaments[i], amrex::Real(0.0), sE, x, y, dummyB, e);
  }
  for (int i = 0; i < p.nToroidals; ++i) {
    const int g = p.toroidals[i].groupIdx;
    const amrex::Real sE = p.groupScaleE ? p.groupScaleE[g] : amrex::Real(0.0);
    eval_toroidal_ba(p.toroidals[i], amrex::Real(0.0), sE, x, y, z, dummyB, e);
  }
}

//==========================================================
// Everything else about B0/E0: SI inputs for #DIPOLE, #CRUSTALFIELD, #COIL,
// #TOROIDALFIELD, #COILWAVEFORM, conversion to code units, waveform time
// updates, and spherical harmonic evaluation.
//==========================================================
class IntrinsicBField {
public:
  enum class CoilKind { Loop, Solenoid, Segment, Filament2D };

  enum class WaveformKind {
    Constant,
    LinearRamp,
    CosineRamp,
    Sinusoid,
    CrowbarRLC,
    Table
  };

  struct CoilSpecSI {
    CoilKind kind = CoilKind::Loop;
    std::string group = "default";
    amrex::Real p1[3] = { 0.0, 0.0, 0.0 };
    amrex::Real p2[3] = { 0.0, 0.0, 1.0 };
    amrex::Real radius = 0.0;
    amrex::Real length = 0.0;
    int nTurns = 1;
    amrex::Real current = 0.0;
    amrex::Real coreRadius = 0.0;
  };

  struct ToroidalSpecSI {
    std::string group = "default";
    amrex::Real center[3] = { 0.0, 0.0, 0.0 };
    amrex::Real axis[3] = { 0.0, 0.0, 1.0 };
    amrex::Real bRefT = 0.0;
    amrex::Real rRefM = 1.0;
    amrex::Real rCoreM = 0.0;
  };

  struct WaveformSpecSI {
    WaveformKind kind = WaveformKind::Constant;
    amrex::Real scale = 1.0;
    amrex::Real tStart = 0.0;
    amrex::Real tEnd = 0.0;
    amrex::Real fStart = 0.0;
    amrex::Real fEnd = 1.0;
    amrex::Real fOffset = 0.0;
    amrex::Real fAmp = 0.0;
    amrex::Real freqHz = 0.0;
    amrex::Real phaseDeg = 0.0;
    amrex::Real tRise = 0.0;
    amrex::Real tauDecay = 0.0;
    bool doCrowbar = false;
    amrex::Real tauCrowbar = 0.0;
    std::string fileName;
    std::vector<amrex::Real> tableTimeSI;
    std::vector<amrex::Real> tableScale;
  };

  IntrinsicBField() = default;
  ~IntrinsicBField() = default;

  bool is_active() const {
    return useDipole_ || nMax_ > 0 || !coilSpecs_.empty() ||
           !toroidalSpecs_.empty();
  }
  bool use_dipole() const { return useDipole_; }
  bool use_crustal() const { return nMax_ > 0; }
  bool has_coils_or_toroidal() const {
    return !coilSpecs_.empty() || !toroidalSpecs_.empty();
  }
  bool is_time_dependent() const { return hasTimeDepWaveform_; }

  // Raw SI input filled by the read_param branches of #DIPOLE, #CRUSTALFIELD,
  // #COIL, #TOROIDALFIELD, and #COILWAVEFORM.
  void read_dipole(amrex::Real bEqNT, amrex::Real thetaDeg, amrex::Real phiDeg,
                   amrex::Real rRefSI);
  void read_crustal(const std::string &fileName, int nMax);
  void add_coil(const CoilSpecSI &spec);
  void add_toroidal_field(const ToroidalSpecSI &spec);
  void set_waveform(const std::string &group, const WaveformSpecSI &spec);

  // True when #DIPOLE asked for its own reference radius.
  bool has_reference_radius_si() const { return rRefSI_ > 0.0; }

  // Number of harmonic degrees requested by #CRUSTALFIELD (0 = off).
  int get_n_max() const { return nMax_; }

  const std::string &get_file_name() const { return fileName_; }

  // Turn the SI input into code units. Must be called once, after the
  // normalization is finalized. Reads any coefficient/table files on rank 0
  // and broadcasts them.
  void convert_units(amrex::Real si2NoB, amrex::Real si2NoL, amrex::Real si2NoV,
                     amrex::Real rPlanetSI, const amrex::Real *bodyCenterCode,
                     amrex::Real bodyRadiusCode, bool useBody, int nDim);

  // Evaluate waveform amplitudes at code time tCode and refresh managed
  // groupScaleB / groupScaleE arrays. Returns true when values were updated.
  bool update_time(amrex::Real tCode);
  bool is_at_time(amrex::Real tCode) const;

  // Total B0 at a code-unit position, added into b[]. Host only when crustal
  // harmonics are active (O(nMax^2) scratch per point).
  void eval(amrex::Real x, amrex::Real y, amrex::Real z,
            amrex::Real b[3]) const;

  const B0Params &get_params() const { return params_; }

  // Resolved quantities after the conversion.
  amrex::Real get_reference_radius_code() const { return rRefCode_; }
  amrex::Real get_planet_radius_code() const { return rPlanetCode_; }
  amrex::Real get_dipole_strength_code() const { return bEqCode_; }

  // Multi-line summary for the log.
  std::string describe() const;

private:
  // Read BATSRUS crustal field layout into coefC_/coefD_ and broadcast.
  void read_coefficients();
  // Read two-column ASCII waveform table (t_SI, scale) on rank 0 & broadcast.
  static void read_waveform_table(WaveformSpecSI &spec);
  // Evaluate waveform scale f(tSI) and derivative df/dtSI [1/s].
  static void eval_waveform(const WaveformSpecSI &spec, amrex::Real tSI,
                            amrex::Real &scale, amrex::Real &dScaleDtSI);
  // Fold the rPlanet^(n+2) power and the unit conversion into the coefficients
  // so that the evaluation works in code units.
  void rescale_coefficients(amrex::Real si2NoB);
  // Add Cartesian crustal field into b[] for a code-unit position.
  void eval_crustal(amrex::Real x, amrex::Real y, amrex::Real z,
                    amrex::Real b[3]) const;

  bool useDipole_ = false;
  amrex::Real bEqNT_ = 0.0;
  amrex::Real thetaDeg_ = 0.0;
  amrex::Real phiDeg_ = 0.0;
  // Negative means "not given": fall back to #BODY / #PLANETRADIUS.
  amrex::Real rRefSI_ = -1.0;

  // Number of harmonic degrees to read (the BATSRUS NNm); the highest retained
  // degree is nMax_ - 1. 0 disables the crustal field.
  int nMax_ = 0;
  std::string fileName_;

  std::vector<CoilSpecSI> coilSpecs_;
  std::vector<ToroidalSpecSI> toroidalSpecs_;
  std::map<std::string, WaveformSpecSI> waveformSpecs_;
  std::vector<std::string> groupNames_;
  std::vector<WaveformSpecSI> groupWaveforms_;

  amrex::Real rRefCode_ = 0.0;
  amrex::Real rPlanetCode_ = 1.0;
  amrex::Real bEqCode_ = 0.0;
  amrex::Real si2NoT_ = 0.0;
  int nDim_ = 3;
  bool hasTimeDepWaveform_ = false;
  bool hasEvaluatedTime_ = false;
  amrex::Real lastEvalTimeCode_ = 0.0;

  amrex::Gpu::ManagedVector<amrex::Real> coefC_, coefD_;
  amrex::Gpu::ManagedVector<CoilLoopParam> loops_;
  amrex::Gpu::ManagedVector<WireSegmentParam> segments_;
  amrex::Gpu::ManagedVector<Filament2DParam> filaments_;
  amrex::Gpu::ManagedVector<ToroidalFieldParam> toroidals_;
  amrex::Gpu::ManagedVector<amrex::Real> groupScaleB_;
  amrex::Gpu::ManagedVector<amrex::Real> groupScaleE_;
  B0Params params_;
};

#endif
