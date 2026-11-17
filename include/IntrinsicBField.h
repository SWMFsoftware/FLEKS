#ifndef _INTRINSICBFIELD_H_
#define _INTRINSICBFIELD_H_

#include <cmath>
#include <string>

#include <AMReX_GpuContainers.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "ExternalCoilField.h"

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

  CoilFieldDeviceParams coils;
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

// Add all GPU-callable intrinsic B0 contributions (dipole + external coils /
// toroidal fields) at (x, y, z) into b[].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_intrinsic_b_dev(
    const B0Params &p, amrex::Real x, amrex::Real y, amrex::Real z,
    int activeDim, amrex::Real b[3]) {
  if (p.useDipole) {
    const amrex::Real zDip = (activeDim > 2) ? z : p.center[2];
    eval_dipole_b(p, x, y, zDip, b);
  }
  eval_external_coils_b_dev(p.coils, x, y, z, b);
}

// Add all GPU-callable inductive E0(t) = -dA0/dt contributions at (x, y, z)
// into e[].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_inductive_e_dev(
    const B0Params &p, amrex::Real x, amrex::Real y, amrex::Real z,
    amrex::Real e[3]) {
  eval_external_coils_e_dev(p.coils, x, y, z, e);
}

// Fused GPU-callable evaluation of both B0 and inductive E0 at (x, y, z).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_intrinsic_be_dev(
    const B0Params &p, amrex::Real x, amrex::Real y, amrex::Real z,
    int activeDim, amrex::Real b[3], amrex::Real e[3]) {
  if (p.useDipole) {
    const amrex::Real zDip = (activeDim > 2) ? z : p.center[2];
    eval_dipole_b(p, x, y, zDip, b);
  }
  eval_external_coils_be_dev(p.coils, x, y, z, b, e);
}

//==========================================================
// Management class for intrinsic planetary (#DIPOLE, #CRUSTALFIELD) and
// composed external coil/toroidal fields (#COIL, #TOROIDALFIELD,
// #COILWAVEFORM): SI input storage, unit conversion, waveform updates, and
// spherical harmonic evaluation.
//==========================================================
class IntrinsicBField {
public:
  IntrinsicBField() = default;
  ~IntrinsicBField() = default;

  bool is_active() const {
    return useDipole_ || nMax_ > 0 || coils_.is_active();
  }
  bool use_dipole() const { return useDipole_; }
  bool use_crustal() const { return nMax_ > 0; }
  bool has_coils_or_toroidal() const { return coils_.is_active(); }
  bool is_time_dependent() const { return coils_.is_time_dependent(); }
  bool has_nonzero_inductive_E() const { return coils_.has_nonzero_e(); }

  ExternalCoilField &coils() { return coils_; }
  const ExternalCoilField &coils() const { return coils_; }

  // Raw SI input filled by #DIPOLE and #CRUSTALFIELD parameter parsing.
  void read_dipole(amrex::Real bEqNT, amrex::Real thetaDeg, amrex::Real phiDeg,
                   amrex::Real rRefSI);
  void read_crustal(const std::string &fileName, int nMax);

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

  // Evaluate waveform amplitudes at code time tCode and refresh device scale
  // descriptors. Returns independent change flags for B0 and E0.
  ExternalCoilField::TimeUpdateStatus update_time(amrex::Real tCode);
  bool is_at_time(amrex::Real tCode) const { return coils_.is_at_time(tCode); }

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

  amrex::Real rRefCode_ = 0.0;
  amrex::Real rPlanetCode_ = 1.0;
  amrex::Real bEqCode_ = 0.0;
  int nDim_ = 3;

  ExternalCoilField coils_;
  amrex::Gpu::ManagedVector<amrex::Real> coefC_, coefD_;
  B0Params params_;
};

#endif
