#ifndef _INTRINSICBFIELD_H_
#define _INTRINSICBFIELD_H_

#include <string>

#include <AMReX_GpuContainers.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

// Static intrinsic magnetic field B0: an analytic dipole and/or a crustal
// field from a spherical harmonic expansion. B0 is frozen: it never enters
// Faraday's law and it is not written to the restart files. It only takes part
// where a total magnetic field is needed, as B = B1 + B0. Keeping the two
// apart is what stops the div(B) cleaning from eroding the planetary field.
struct B0Params {
  // Center of both the dipole and the spherical harmonic expansion, code
  // units.
  amrex::Real center[3] = {0.0, 0.0, 0.0};

  // Dipole moment m in code units: B = [3 (m.r^) r^ - m] / r^3, so
  // |m| = Bc * rRef^3 for the equatorial strength Bc at the reference radius.
  amrex::Real moment[3] = {0.0, 0.0, 0.0};

  // Crustal coefficients, flat layout coef[n * nMax + m], already rescaled to
  // code units. Null when the crustal field is off.
  const amrex::Real *coefC = nullptr;
  const amrex::Real *coefD = nullptr;

  // Number of harmonic degrees (n = 0..nMax - 1); 0 disables the crustal part.
  int nMax = 0;
};

// Dipole part of B0 at the code-unit position (x, y, z), added into b[].
// Device callable and scratch free, so it may run inside a ParallelFor.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
eval_dipole_b(const B0Params &p, amrex::Real x, amrex::Real y, amrex::Real z,
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

//==========================================================
// Everything else about B0: the SI input of the two commands, the conversion
// to code units and the spherical harmonic evaluation.
//==========================================================
class IntrinsicBField {
public:
  IntrinsicBField() = default;
  ~IntrinsicBField() = default;

  bool is_active() const { return useDipole_ || nMax_ > 0; }
  bool use_dipole() const { return useDipole_; }
  bool use_crustal() const { return nMax_ > 0; }

  // Raw SI input filled by the read_param branches of #DIPOLE and
  // #CRUSTALFIELD. Nothing is converted yet: the normalization may not have
  // been read at that point.
  void read_dipole(amrex::Real bEqNT, amrex::Real thetaDeg, amrex::Real phiDeg,
                   amrex::Real rRefSI);
  void read_crustal(const std::string &fileName, int nMax);

  // True when #DIPOLE asked for its own reference radius.
  bool has_reference_radius_si() const { return rRefSI_ > 0.0; }

  // Number of harmonic degrees requested by #CRUSTALFIELD (0 = off).
  int get_n_max() const { return nMax_; }

  const std::string &get_file_name() const { return fileName_; }

  // Turn the SI input into code units. Must be called once, after the
  // normalization is finalized. Reads the coefficient file on rank 0 and
  // broadcasts it.
  void convert_units(double si2NoB, double si2NoL, double rPlanetSI,
                     const double *bodyCenterCode, double bodyRadiusCode,
                     bool useBody, int nDim);

  // Total B0 at a code-unit position, added into b[]. Host only: the crustal
  // part needs O(nMax^2) scratch per point.
  void eval(amrex::Real x, amrex::Real y, amrex::Real z,
            amrex::Real b[3]) const;

  const B0Params &get_params() const { return params_; }

  // Resolved quantities after the conversion.
  double get_reference_radius_code() const { return rRefCode_; }
  double get_planet_radius_code() const { return rPlanetCode_; }
  double get_dipole_strength_code() const { return bEqCode_; }

  // Multi-line summary for the log.
  std::string describe() const;

private:
  // Read BATSRUS crustal field layout into coefC_/coefD_ and broadcast.
  void read_coefficients();
  // Fold the rPlanet^(n+2) power and the unit conversion into the coefficients
  // so that the evaluation works in code units.
  void rescale_coefficients(double si2NoB);
  // Add Cartesian crustal field into b[] for a code-unit position.
  void eval_crustal(amrex::Real x, amrex::Real y, amrex::Real z,
                    amrex::Real b[3]) const;

  bool useDipole_ = false;
  amrex::Real bEqNT_ = 0.0;
  amrex::Real thetaDeg_ = 0.0;
  amrex::Real phiDeg_ = 0.0;
  // < 0 means "not given": fall back to #BODY / #PLANETRADIUS.
  amrex::Real rRefSI_ = -1.0;

  // Number of harmonic degrees to read (the BATSRUS NNm); the highest retained
  // degree is nMax_ - 1. 0 disables the crustal field.
  int nMax_ = 0;
  std::string fileName_;

  double rRefCode_ = 0.0;
  double rPlanetCode_ = 1.0;
  double bEqCode_ = 0.0;
  int nDim_ = 3;

  amrex::Gpu::ManagedVector<amrex::Real> coefC_, coefD_;
  B0Params params_;
};

#endif

