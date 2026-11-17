#ifndef _EXTERNALCOILFIELD_H_
#define _EXTERNALCOILFIELD_H_

#include <cmath>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <AMReX_Array.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "ReadParam.h"

// GPU-safe POD geometry descriptors for external current primitives (#COIL and
// #TOROIDALFIELD), stored in code units with invariants precomputed once in
// convert_units().
struct CoilLoopParam {
  amrex::Real center[3] = { 0.0, 0.0, 0.0 };
  amrex::Real axis[3] = { 0.0, 0.0, 1.0 };
  amrex::Real radius = 0.0;
  amrex::Real r2 = 0.0;
  amrex::Real fourR = 0.0;
  amrex::Real coreEps2 = 0.0;
  // Nominal (mu0 * N * I / 4*pi) converted to code_B * code_L (includes static
  // Constant waveform scale when no time-dependent waveform is active).
  amrex::Real kI0 = 0.0;
  int groupIdx = 0;
};

struct WireSegmentParam {
  amrex::Real x1[3] = { 0.0, 0.0, 0.0 };
  amrex::Real dir[3] = { 1.0, 0.0, 0.0 };
  amrex::Real length = 0.0;
  amrex::Real coreEps2 = 0.0;
  amrex::Real kI0 = 0.0;
  int groupIdx = 0;
};

struct Filament2DParam {
  amrex::Real x0 = 0.0;
  amrex::Real y0 = 0.0;
  amrex::Real coreEps2 = 0.0;
  amrex::Real invRRef2 = 1.0;
  amrex::Real kI0 = 0.0;
  int groupIdx = 0;
};

struct ToroidalFieldParam {
  amrex::Real center[3] = { 0.0, 0.0, 0.0 };
  amrex::Real axis[3] = { 0.0, 0.0, 1.0 };
  amrex::Real coreEps2 = 0.0;
  amrex::Real invRRef2 = 1.0;
  // Nominal (B_ref * R_ref) in code_B * code_L (includes static Constant scale
  // when no time-dependent waveform is active).
  amrex::Real kTor0 = 0.0;
  int groupIdx = 0;
};

// Lightweight value-captured device descriptor for all external coil/toroidal
// primitives and instantaneous per-group scales. When nGroups <=
// kMaxInlineGroups (16), scales are captured by value inside kernels, avoiding
// UVM page faults and asynchronous host/device races across RK/subcycle stages.
struct CoilFieldDeviceParams {
  static constexpr int kMaxInlineGroups = 16;

  const CoilLoopParam *loops = nullptr;
  int nLoops = 0;

  const WireSegmentParam *segments = nullptr;
  int nSegments = 0;

  const Filament2DParam *filaments = nullptr;
  int nFilaments = 0;

  const ToroidalFieldParam *toroidals = nullptr;
  int nToroidals = 0;

  amrex::GpuArray<amrex::Real, kMaxInlineGroups> inlineScaleB{};
  amrex::GpuArray<amrex::Real, kMaxInlineGroups> inlineScaleE{};
  const amrex::Real *extScaleB = nullptr;
  const amrex::Real *extScaleE = nullptr;
  int nGroups = 0;
  bool hasTimeDep = false;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real scale_b(int g) const {
    if (!hasTimeDep)
      return amrex::Real(1.0);
    return extScaleB ? extScaleB[g] : inlineScaleB[g];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real scale_e(int g) const {
    if (!hasTimeDep)
      return amrex::Real(0.0);
    return extScaleE ? extScaleE[g] : inlineScaleE[g];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool any_active() const {
    return nLoops > 0 || nSegments > 0 || nFilaments > 0 || nToroidals > 0;
  }
};

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
// Uses Plummer axial softening Z_eff^2 = zeta^2 + coreEps2, preserving
// curl(A0) == B0 identically even for coreRadius > 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_coil_loop_be(
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
  const amrex::Real zEff2 = zeta * zeta + c.coreEps2;

  const amrex::Real dMinus = c.radius - rho;
  const amrex::Real dPlus = c.radius + rho;
  const amrex::Real alpha2 = dMinus * dMinus + zEff2;
  const amrex::Real beta2 = dPlus * dPlus + zEff2;
  const amrex::Real beta = std::sqrt(beta2);

  amrex::Real m = (c.fourR * rho) / beta2;
  if (m > amrex::Real(1.0 - 1.0e-15))
    m = amrex::Real(1.0 - 1.0e-15);

  amrex::Real K = 0.0;
  amrex::Real E = 0.0;
  eval_elliptic_ke(m, K, E);

  if (scaleB != amrex::Real(0.0)) {
    const amrex::Real kB = c.kI0 * scaleB;
    const amrex::Real bZeta = (amrex::Real(2.0) * kB / beta) *
                              (K + ((c.r2 - rho2 - zEff2) / alpha2) * E);
    amrex::Real bRhoOverRho = 0.0;
    if (m < amrex::Real(1.0e-3)) {
      const amrex::Real beta5 = beta2 * beta2 * beta;
      bRhoOverRho =
          (amrex::Real(8.0) * cPi * kB * c.r2 * zeta / beta5) *
          (amrex::Real(0.375) +
           m * (amrex::Real(0.46875) + m * (amrex::Real(525.0 / 1024.0) +
                                            m * amrex::Real(2205.0 / 4096.0))));
    } else {
      bRhoOverRho = (amrex::Real(2.0) * kB * zeta / (rho2 * beta)) *
                    (-K + ((c.r2 + rho2 + zEff2) / alpha2) * E);
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
          (amrex::Real(8.0) * cPi * kE * c.r2 / beta3) *
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_wire_segment_be(
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

  const amrex::Real rhoEff2 = px * px + py * py + pz * pz + s.coreEps2;

  if (scaleB != amrex::Real(0.0)) {
    const amrex::Real r1 = std::sqrt(s1 * s1 + rhoEff2);
    const amrex::Real r2 = std::sqrt(s2 * s2 + rhoEff2);
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_filament2d_be(
    const Filament2DParam &f, amrex::Real scaleB, amrex::Real scaleE,
    amrex::Real x, amrex::Real y, amrex::Real b[3], amrex::Real e[3]) {
  if (scaleB == amrex::Real(0.0) && scaleE == amrex::Real(0.0))
    return;

  const amrex::Real dx = x - f.x0;
  const amrex::Real dy = y - f.y0;
  const amrex::Real rEff2 = dx * dx + dy * dy + f.coreEps2;

  if (scaleB != amrex::Real(0.0)) {
    const amrex::Real facB = amrex::Real(2.0) * f.kI0 * scaleB / rEff2;
    b[0] += -facB * dy;
    b[1] += facB * dx;
  }

  if (scaleE != amrex::Real(0.0)) {
    e[2] += -(f.kI0 * scaleE) * std::log(rEff2 * f.invRRef2);
  }
}

// Vacuum 1/R_perp toroidal magnetic field around (center, axis), regularized
// with core radius rCore so div(B) == 0 identically everywhere.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_toroidal_be(
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

  const amrex::Real rhoEff2 = rx * rx + ry * ry + rz * rz + t.coreEps2;

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
    const amrex::Real facE =
        -amrex::Real(0.5) * (t.kTor0 * scaleE) * std::log(rhoEff2 * t.invRRef2);
    e[0] += facE * t.axis[0];
    e[1] += facE * t.axis[1];
    e[2] += facE * t.axis[2];
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_external_coils_b_dev(
    const CoilFieldDeviceParams &p, amrex::Real x, amrex::Real y, amrex::Real z,
    amrex::Real b[3]) {
  amrex::Real dummyE[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < p.nLoops; ++i) {
    eval_coil_loop_be(p.loops[i], p.scale_b(p.loops[i].groupIdx),
                      amrex::Real(0.0), x, y, z, b, dummyE);
  }
  for (int i = 0; i < p.nSegments; ++i) {
    eval_wire_segment_be(p.segments[i], p.scale_b(p.segments[i].groupIdx),
                         amrex::Real(0.0), x, y, z, b, dummyE);
  }
  for (int i = 0; i < p.nFilaments; ++i) {
    eval_filament2d_be(p.filaments[i], p.scale_b(p.filaments[i].groupIdx),
                       amrex::Real(0.0), x, y, b, dummyE);
  }
  for (int i = 0; i < p.nToroidals; ++i) {
    eval_toroidal_be(p.toroidals[i], p.scale_b(p.toroidals[i].groupIdx),
                     amrex::Real(0.0), x, y, z, b, dummyE);
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_external_coils_e_dev(
    const CoilFieldDeviceParams &p, amrex::Real x, amrex::Real y, amrex::Real z,
    amrex::Real e[3]) {
  amrex::Real dummyB[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < p.nLoops; ++i) {
    eval_coil_loop_be(p.loops[i], amrex::Real(0.0),
                      p.scale_e(p.loops[i].groupIdx), x, y, z, dummyB, e);
  }
  for (int i = 0; i < p.nSegments; ++i) {
    eval_wire_segment_be(p.segments[i], amrex::Real(0.0),
                         p.scale_e(p.segments[i].groupIdx), x, y, z, dummyB, e);
  }
  for (int i = 0; i < p.nFilaments; ++i) {
    eval_filament2d_be(p.filaments[i], amrex::Real(0.0),
                       p.scale_e(p.filaments[i].groupIdx), x, y, dummyB, e);
  }
  for (int i = 0; i < p.nToroidals; ++i) {
    eval_toroidal_be(p.toroidals[i], amrex::Real(0.0),
                     p.scale_e(p.toroidals[i].groupIdx), x, y, z, dummyB, e);
  }
}

// Fused evaluation of both B0(t) and E0(t) = -dA0/dt in one element pass,
// sharing geometric transforms and elliptic integrals across both fields.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_external_coils_be_dev(
    const CoilFieldDeviceParams &p, amrex::Real x, amrex::Real y, amrex::Real z,
    amrex::Real b[3], amrex::Real e[3]) {
  for (int i = 0; i < p.nLoops; ++i) {
    const int g = p.loops[i].groupIdx;
    eval_coil_loop_be(p.loops[i], p.scale_b(g), p.scale_e(g), x, y, z, b, e);
  }
  for (int i = 0; i < p.nSegments; ++i) {
    const int g = p.segments[i].groupIdx;
    eval_wire_segment_be(p.segments[i], p.scale_b(g), p.scale_e(g), x, y, z, b,
                         e);
  }
  for (int i = 0; i < p.nFilaments; ++i) {
    const int g = p.filaments[i].groupIdx;
    eval_filament2d_be(p.filaments[i], p.scale_b(g), p.scale_e(g), x, y, b, e);
  }
  for (int i = 0; i < p.nToroidals; ++i) {
    const int g = p.toroidals[i].groupIdx;
    eval_toroidal_be(p.toroidals[i], p.scale_b(g), p.scale_e(g), x, y, z, b, e);
  }
}

//==========================================================
// External current coils (#COIL, #TOROIDALFIELD) and current waveforms
// (#COILWAVEFORM): parsing, unit conversion, waveform evaluation, and
// device-ready descriptor management.
//==========================================================
class ExternalCoilField {
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

  struct TimeUpdateStatus {
    bool changedB = false;
    bool changedE = false;
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
    // Precomputed constants populated during convert_units().
    amrex::Real omega = 0.0;
    amrex::Real phase0 = 0.0;
    amrex::Real gamma = 0.0;
    amrex::Real tPeak = 0.0;
    amrex::Real aNorm = 1.0;
    amrex::Real invTauCrowbar = 0.0;
  };

  ExternalCoilField() = default;
  ~ExternalCoilField() = default;

  bool is_active() const {
    return !coilSpecs_.empty() || !toroidalSpecs_.empty();
  }
  bool is_time_dependent() const { return hasTimeDepWaveform_; }
  bool has_nonzero_e() const { return anyNonzeroE_; }

  void read_coil_param(ReadParam &param);
  void read_toroidal_param(ReadParam &param);
  void read_waveform_param(ReadParam &param);

  void add_coil(const CoilSpecSI &spec);
  void add_toroidal_field(const ToroidalSpecSI &spec);
  void set_waveform(const std::string &group, const WaveformSpecSI &spec);

  void convert_units(amrex::Real si2NoB, amrex::Real si2NoL,
                     amrex::Real si2NoV);

  TimeUpdateStatus update_time(amrex::Real tCode);
  bool is_at_time(amrex::Real tCode) const;

  const CoilFieldDeviceParams &get_params() const { return params_; }

  void describe(std::ostream &os) const;

private:
  static void read_waveform_table(WaveformSpecSI &spec);
  static void precompute_waveform(WaveformSpecSI &spec);
  static void eval_waveform(const WaveformSpecSI &spec, amrex::Real tSI,
                            amrex::Real &scale, amrex::Real &dScaleDtSI);

  std::vector<CoilSpecSI> coilSpecs_;
  std::vector<ToroidalSpecSI> toroidalSpecs_;
  std::map<std::string, WaveformSpecSI> waveformSpecs_;
  std::vector<std::string> groupNames_;
  std::vector<WaveformSpecSI> groupWaveforms_;

  amrex::Real si2NoT_ = 0.0;
  bool hasTimeDepWaveform_ = false;
  bool anyNonzeroE_ = false;
  bool hasEvaluatedTime_ = false;
  amrex::Real lastEvalTimeCode_ = 0.0;

  std::vector<amrex::Real> curScaleB_;
  std::vector<amrex::Real> curScaleE_;

  amrex::Gpu::DeviceVector<CoilLoopParam> loops_;
  amrex::Gpu::DeviceVector<WireSegmentParam> segments_;
  amrex::Gpu::DeviceVector<Filament2DParam> filaments_;
  amrex::Gpu::DeviceVector<ToroidalFieldParam> toroidals_;
  amrex::Gpu::DeviceVector<amrex::Real> devScaleB_;
  amrex::Gpu::DeviceVector<amrex::Real> devScaleE_;

  CoilFieldDeviceParams params_;
};

#endif // _EXTERNALCOILFIELD_H_
