#ifndef _IONOSPHERE_H_
#define _IONOSPHERE_H_

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <AMReX_GpuContainers.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "Constants.h"

//==========================================================
// Self-contained Magnetosphere-Ionosphere (MI) coupling module.
//
// Solves height-integrated ionospheric current continuity in spherical
// magnetic dipole coordinates (SMG):
//
//   div_h ( Sigma . grad_h Phi ) = J_r^{outward}
//
// where E_h = -grad_h Phi, J_h = Sigma . E_h, and J_r^{outward} = J_|| *
// (B_r / |B|) is the outward radial component of the field-aligned current
// mapped from the magnetospheric extraction shell r = rFac to the ionospheric
// shell r = rIono.
//
// The 2D PDE is solved entirely in SI units (m, S, A/m^2, V) on rank 0 and
// broadcast to all ranks. Mapped electric fields E = -grad Phi are evaluated
// on the inner-body surface nodes via the exact analytic dipole mapping and
// its 2x3 Jacobian d(theta, phi)/d(x, y, z), converted to FLEKS code units.
//==========================================================

struct DipoleFrame {
  // Body / dipole center in simulation code units.
  amrex::Real center[3] = { 0.0, 0.0, 0.0 };

  // Orthonormal rotation matrix R from simulation frame to SMG frame:
  //   x_smg = R * (x_sim - center),   x_sim = center + R^T * x_smg.
  // Row 0 = e_x (sunward projected orthogonal to dipole axis),
  // Row 1 = e_y = e_z x e_x,
  // Row 2 = e_z (northern magnetic dipole axis, e_z.z >= 0).
  amrex::Real rot[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  // Sign of dipole moment along +e_z:
  //   -1.0 for Earth-like polarity (B0 enters North pole, B_z < 0 at equator),
  //   +1.0 for reversed polarity (B0 exits North pole, B_z > 0 at equator).
  amrex::Real dipoleSign = -1.0;

  // Ionosphere radius in code units and SI [m].
  amrex::Real rIonoCode = 1.0;
  amrex::Real rIonoSI = 6.481e6;

  // Conversion from SI electric field [V/m] to FLEKS code units.
  amrex::Real si2NoE = 1.0;
  // Conversion from code length to SI [m].
  amrex::Real no2SiL = 1.0;
};

//==========================================================
// Analytic dipole field-line mapping in SMG coordinates (device-callable).
//
// Maps xSmg (at radius r) along a dipole field line r = L sin^2(theta) to
// target radius rMap (same length units as xSmg).
//
// Returns iHemisphere:
//   +1 (north, z >= 0), -1 (south, z < 0), or 0 (field line does not reach
//   rMap, e.g. equatorial field line with L < rMap).
//
// If computeJac is true, fills jacSmg[2][3] = d(theta_m, phi_m)/d(x, y, z)
// in SMG coordinates (units: 1 / length_unit).
//==========================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int map_dipole_smg(
    const amrex::Real xSmg[3], amrex::Real rMap, amrex::Real xMap[3],
    amrex::Real &thetaMap, amrex::Real &phiMap, bool computeJac,
    amrex::Real jacSmg[2][3]) {
  const amrex::Real x = xSmg[0];
  const amrex::Real y = xSmg[1];
  const amrex::Real z = xSmg[2];

  const amrex::Real r2 = x * x + y * y + z * z;
  const amrex::Real r = std::sqrt(r2);
  if (r <= amrex::Real(1.0e-12) || rMap <= amrex::Real(1.0e-12)) {
    xMap[0] = xMap[1] = xMap[2] = 0.0;
    thetaMap = 0.0;
    phiMap = 0.0;
    if (computeJac) {
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 3; ++b)
          jacSmg[a][b] = 0.0;
    }
    return 0;
  }

  const int iHem = (z >= amrex::Real(0.0)) ? 1 : -1;
  const amrex::Real ratio = rMap / r;
  const amrex::Real xyRatio = ratio * std::sqrt(ratio); // (rMap / r)^(3/2)

  xMap[0] = xyRatio * x;
  xMap[1] = xyRatio * y;
  const amrex::Real rhoMap2 = xMap[0] * xMap[0] + xMap[1] * xMap[1];
  const amrex::Real rMap2 = rMap * rMap;

  if (rhoMap2 > rMap2) {
    // Field line closes below rMap (L < rMap).
    const amrex::Real rho = std::sqrt(x * x + y * y);
    if (rho > amrex::Real(0.0)) {
      xMap[0] = (rMap / rho) * x;
      xMap[1] = (rMap / rho) * y;
    } else {
      xMap[0] = rMap;
      xMap[1] = 0.0;
    }
    xMap[2] = 0.0;
    thetaMap = amrex::Real(0.5) * dPI;
    phiMap = std::atan2(xMap[1], xMap[0]);
    if (phiMap < amrex::Real(0.0))
      phiMap += amrex::Real(2.0) * dPI;
    if (computeJac) {
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 3; ++b)
          jacSmg[a][b] = 0.0;
    }
    return 0;
  }

  xMap[2] = iHem * std::sqrt(std::max(amrex::Real(0.0), rMap2 - rhoMap2));

  amrex::Real cosTh = xMap[2] / rMap;
  if (cosTh > amrex::Real(1.0))
    cosTh = amrex::Real(1.0);
  if (cosTh < amrex::Real(-1.0))
    cosTh = amrex::Real(-1.0);
  thetaMap = std::acos(cosTh);

  phiMap = std::atan2(xMap[1], xMap[0]);
  if (phiMap < amrex::Real(0.0))
    phiMap += amrex::Real(2.0) * dPI;

  if (computeJac) {
    const amrex::Real rho2 = x * x + y * y;
    const amrex::Real rho = std::sqrt(rho2);
    const amrex::Real rhoMap = std::sqrt(rhoMap2);
    const amrex::Real zMap = xMap[2];

    if (std::abs(zMap) < amrex::Real(1.0e-6) * rMap ||
        rho < amrex::Real(1.0e-8) * r) {
      // Near magnetic equator (unmapped / low-lat zero potential) or exact
      // pole axis (where dPhi/dtheta, dPhi/dphi -> 0 by regularity): zero
      // Jacobian to avoid 0/0.
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 3; ++b)
          jacSmg[a][b] = 0.0;
    } else {
      const amrex::Real zOverR = z / r;
      const amrex::Real factor =
          -(amrex::Real(0.5) - amrex::Real(1.5) * zOverR * zOverR) /
          (zMap * rho);
      jacSmg[0][0] = xMap[0] * factor;
      jacSmg[0][1] = xMap[1] * factor;
      jacSmg[0][2] = -amrex::Real(1.5) * (rhoMap / zMap) * (z / r2);

      jacSmg[1][0] = -y / rho2;
      jacSmg[1][1] = x / rho2;
      jacSmg[1][2] = 0.0;
    }
  }

  return iHem;
}

//==========================================================
// GPU-ready view of the solved ionospheric potential gradient.
//==========================================================
struct IonosphereDeviceView {
  DipoleFrame frame;
  int nTheta = 0; // Full pole-to-pole grid size (odd, 0..pi)
  int nPhi = 0;   // Periodic longitude grid size (0..2pi, last == first)
  amrex::Real dTheta = 1.0;
  amrex::Real dPhi = 1.0;
  amrex::Real latBoundaryRad = 0.0;        // Low-latitude cutoff [rad]
  const amrex::Real *dPhiDTheta = nullptr; // SI [V/rad], size nTheta * nPhi
  const amrex::Real *dPhiDPhi = nullptr;   // SI [V/rad], size nTheta * nPhi
};

//==========================================================
// Evaluate mapped ionospheric electric field E_iono = -grad(Phi) in FLEKS
// code units at simulation position (x, y, z). Device-callable.
//==========================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void eval_ionosphere_E_dev(
    const IonosphereDeviceView &v, amrex::Real xSim, amrex::Real ySim,
    amrex::Real zSim, amrex::Real eSim[3]) {
  eSim[0] = eSim[1] = eSim[2] = 0.0;
  if (!v.dPhiDTheta || !v.dPhiDPhi || v.nTheta < 3 || v.nPhi < 3)
    return;

  // Shift to body center and rotate to SMG frame.
  const amrex::Real dx[3] = {
    xSim - v.frame.center[0],
    ySim - v.frame.center[1],
    zSim - v.frame.center[2],
  };
  amrex::Real xSmg[3] = { 0.0, 0.0, 0.0 };
  for (int a = 0; a < 3; ++a) {
    xSmg[a] = v.frame.rot[a][0] * dx[0] + v.frame.rot[a][1] * dx[1] +
              v.frame.rot[a][2] * dx[2];
  }

  amrex::Real xMap[3] = { 0.0, 0.0, 0.0 };
  amrex::Real theta = 0.0, phi = 0.0;
  amrex::Real jacSmg[2][3];
  const int iHem =
      map_dipole_smg(xSmg, v.frame.rIonoCode, xMap, theta, phi, true, jacSmg);
  if (iHem == 0)
    return;

  // Below the low-latitude boundary, mapped potential is zero (or constant).
  const amrex::Real magLat = std::abs(amrex::Real(0.5) * dPI - theta);
  if (magLat <= v.latBoundaryRad)
    return;

  // Bilinear interpolation of (dPhi/dTheta, dPhi/dPhi) on the (nTheta x nPhi)
  // uniform spherical grid.
  const amrex::Real tNorm = theta / v.dTheta;
  const amrex::Real pNorm = phi / v.dPhi;

  int i0 = static_cast<int>(std::floor(tNorm));
  int j0 = static_cast<int>(std::floor(pNorm));
  if (i0 < 0)
    i0 = 0;
  if (i0 >= v.nTheta - 1)
    i0 = v.nTheta - 2;
  if (j0 < 0)
    j0 = 0;
  if (j0 >= v.nPhi - 1)
    j0 = v.nPhi - 2;

  const amrex::Real wt = tNorm - i0;
  const amrex::Real wp = pNorm - j0;

  const int idx00 = i0 * v.nPhi + j0;
  const int idx01 = i0 * v.nPhi + (j0 + 1);
  const int idx10 = (i0 + 1) * v.nPhi + j0;
  const int idx11 = (i0 + 1) * v.nPhi + (j0 + 1);

  const amrex::Real dPdT =
      (1.0 - wt) *
          ((1.0 - wp) * v.dPhiDTheta[idx00] + wp * v.dPhiDTheta[idx01]) +
      wt * ((1.0 - wp) * v.dPhiDTheta[idx10] + wp * v.dPhiDTheta[idx11]);
  const amrex::Real dPdP =
      (1.0 - wt) * ((1.0 - wp) * v.dPhiDPhi[idx00] + wp * v.dPhiDPhi[idx01]) +
      wt * ((1.0 - wp) * v.dPhiDPhi[idx10] + wp * v.dPhiDPhi[idx11]);

  // E_smg [V/m] = -(dPhi/dTheta * grad_SI(theta) + dPhi/dPhi * grad_SI(phi))
  // Since jacSmg is in [1 / code_length], divide by no2SiL to get [1 / m],
  // then multiply by si2NoE to convert [V/m] -> FLEKS code units.
  const amrex::Real scale = -v.frame.si2NoE / v.frame.no2SiL;
  amrex::Real eSmg[3];
  for (int b = 0; b < 3; ++b) {
    eSmg[b] = scale * (dPdT * jacSmg[0][b] + dPdP * jacSmg[1][b]);
  }

  // Rotate back from SMG to simulation frame: eSim = R^T * eSmg.
  for (int b = 0; b < 3; ++b) {
    eSim[b] = v.frame.rot[0][b] * eSmg[0] + v.frame.rot[1][b] * eSmg[1] +
              v.frame.rot[2][b] * eSmg[2];
  }
}

//==========================================================
// Configuration and solver class for the 2D spherical ionosphere.
//==========================================================
struct IonosphereParams {
  // Spherical grid resolution:
  // nThetaHem = number of colatitude points from pole to equator (inclusive).
  // Full pole-to-pole grid has nTheta = 2 * nThetaHem - 1 points.
  int nThetaHem = 91;
  int nPhi = 181; // Includes periodic endpoint phi = 2*pi == 0.

  // Geometry:
  amrex::Real rIonoRatio = 1.017; // rIono / rPlanet (1 R_E + ~110 km)
  amrex::Real rFacCode =
      -1.0; // FAC extraction radius in code units (<= 0: auto)
  amrex::Real latBoundaryDeg = 20.0; // Low-latitude boundary [deg]

  // Conductance model: "uniform", "solar_euv", or "auroral"
  std::string conductanceModel = "uniform";
  amrex::Real sigmaPedersen = 5.0;      // Uniform / nightside Pedersen [S]
  amrex::Real sigmaHall = 0.0;          // Uniform / nightside Hall [S]
  amrex::Real sigmaParallel = 1000.0;   // Field-aligned conductance Sigma_0 [S]
  amrex::Real f107 = 150.0;             // Solar F10.7 radio flux for EUV model
  amrex::Real auroralFacScale = 1.0e-6; // Reference upward FAC [A/m^2]
  amrex::Real auroralSigmaPMax = 15.0;  // Max auroral Pedersen enhancement [S]
  amrex::Real auroralHallRatio = 1.5;   // Auroral Sigma_H / Sigma_P ratio

  // Low-latitude boundary condition: "dirichlet" (Phi=0) or "neumann"
  // (dPhi/dTheta=0)
  std::string lowLatBC = "dirichlet";

  // Linear solver settings:
  amrex::Real solverTol = 1.0e-7;
  int solverMaxIter = 500;

  // Coupling cadence & noise control:
  int dnCouple = 1; // Solve every dnCouple steps (<=0: use dtCouple)
  amrex::Real dtCoupleSI =
      0.0;            // Solve every dtCoupleSI seconds (>0 overrides dn)
  int nSmoothFac = 2; // Spherical smoothing passes on J_r
  amrex::Real relaxAlpha =
      1.0; // Temporal blending: Phi_eff = (1-a)*old + a*new

  // Optional corotation potential:
  bool useCorotation = false;
  amrex::Real omegaPlanetSI = 7.2921159e-5; // Earth sidereal rotation [rad/s]

  // Diagnostic 2D output cadence (0 = disabled):
  int dnOutput = 0;
};

class Ionosphere {
public:
  Ionosphere() = default;
  explicit Ionosphere(const IonosphereParams &params) : params_(params) {}

  void set_params(const IonosphereParams &params) { params_ = params; }
  const IonosphereParams &get_params() const { return params_; }
  IonosphereParams &params() { return params_; }

  // Initialize grid geometry, dipole frame, and conductances.
  // dipoleMomentCode: 3-component dipole moment in code units.
  // bodyCenterCode:   3-component body center in code units.
  // rPlanetSI:        planetary radius [m].
  // si2NoL, si2NoE, si2NoB: unit conversion factors (SI -> code).
  void init(const amrex::Real dipoleMomentCode[3],
            const amrex::Real bodyCenterCode[3], amrex::Real rPlanetSI,
            amrex::Real si2NoL, amrex::Real si2NoE, amrex::Real si2NoB,
            amrex::Real bEqSI = -3.0e-5);

  bool is_initialized() const { return initialized_; }

  // Grid dimensions (full pole-to-pole).
  int n_theta() const { return nTheta_; }
  int n_phi() const { return nPhi_; }
  amrex::Real d_theta() const { return dTheta_; }
  amrex::Real d_phi() const { return dPhi_; }
  const std::vector<amrex::Real> &theta() const { return theta_; }
  const std::vector<amrex::Real> &phi() const { return phi_; }
  const DipoleFrame &frame() const { return frame_; }

  // Set outward radial FAC J_r [A/m^2] on the full (nTheta x nPhi) grid,
  // optionally smooth it, update conductances, solve both hemispheres,
  // compute gradients, and sync to GPU device memory.
  // Returns true if both hemispheres converged.
  bool solve(const std::vector<amrex::Real> &jrOutwardSI,
             bool doBroadcast = true);

  // Access solved fields (host, size nTheta * nPhi, row-major [i*nPhi + j]).
  const std::vector<amrex::Real> &get_potential() const { return phiPot_; }
  const std::vector<amrex::Real> &get_jr() const { return jr_; }
  const std::vector<amrex::Real> &get_sigma_p() const { return sigmaP_; }
  const std::vector<amrex::Real> &get_sigma_h() const { return sigmaH_; }
  const std::vector<amrex::Real> &get_dphi_dtheta() const {
    return dPhiDTheta_;
  }
  const std::vector<amrex::Real> &get_dphi_dphi() const { return dPhiDPhi_; }

  // Restore potential from restart (recomputes gradients and syncs device).
  void set_potential(const std::vector<amrex::Real> &phiIn);

  // Cross-polar-cap potentials [kV] and last solver stats.
  amrex::Real cpcp_north_kV() const { return cpcpNorthKV_; }
  amrex::Real cpcp_south_kV() const { return cpcpSouthKV_; }
  int last_iter_north() const { return lastIterNorth_; }
  int last_iter_south() const { return lastIterSouth_; }
  amrex::Real last_resid_north() const { return lastResidNorth_; }
  amrex::Real last_resid_south() const { return lastResidSouth_; }

  // Device view for GPU/OpenMP kernels in apply_body_E_ionosphere().
  IonosphereDeviceView device_view() const;

  // Write 2D spherical diagnostics file (ASCII table, rank 0 only).
  void write_diagnostics(const std::string &fileName, amrex::Real timeSI,
                         int cycle) const;

  // Build SMG rotation matrix from a dipole moment vector (public for tests).
  static void build_dipole_frame(const amrex::Real dipoleMoment[3],
                                 amrex::Real rot[3][3],
                                 amrex::Real &dipoleSign);

private:
  void compute_conductances();
  void compute_conductance_tensor();
  void smooth_jr();
  bool solve_hemisphere(bool isNorth, int &nIterOut, amrex::Real &residOut);
  void add_corotation_potential();
  void compute_gradients();
  void sync_device_gradients();

  IonosphereParams params_;
  DipoleFrame frame_;
  bool initialized_ = false;
  amrex::Real rPlanetSI_ = 6.371e6;
  amrex::Real bEqSI_ = -3.0e-5;

  int nTheta_ = 0;
  int nPhi_ = 0;
  amrex::Real dTheta_ = 0.0;
  amrex::Real dPhi_ = 0.0;
  amrex::Real latBoundaryRad_ = 0.0;

  std::vector<amrex::Real> theta_;
  std::vector<amrex::Real> phi_;
  std::vector<amrex::Real> sinTheta_;
  std::vector<amrex::Real> cosTheta_;

  // Full-grid 2D arrays (size nTheta_ * nPhi_, index = i * nPhi_ + j)
  std::vector<amrex::Real> jr_;
  std::vector<amrex::Real> sigmaP_;
  std::vector<amrex::Real> sigmaH_;
  std::vector<amrex::Real> sigmaThTh_;
  std::vector<amrex::Real> sigmaThPh_;
  std::vector<amrex::Real> sigmaPhPh_;
  std::vector<amrex::Real> dSigmaThTh_dTh_;
  std::vector<amrex::Real> dSigmaThPh_dTh_;
  std::vector<amrex::Real> dSigmaThPh_dPh_;
  std::vector<amrex::Real> dSigmaPhPh_dPh_;

  std::vector<amrex::Real> phiPot_;     // Current blended potential [V]
  std::vector<amrex::Real> phiSolved_;  // Raw solved potential [V]
  std::vector<amrex::Real> dPhiDTheta_; // [V/rad]
  std::vector<amrex::Real> dPhiDPhi_;   // [V/rad]
  bool hasPrevSolution_ = false;

  // Device copies of gradients for GPU/OpenMP body-BC kernels.
  amrex::Gpu::DeviceVector<amrex::Real> dPhiDThetaDev_;
  amrex::Gpu::DeviceVector<amrex::Real> dPhiDPhiDev_;

  amrex::Real cpcpNorthKV_ = 0.0;
  amrex::Real cpcpSouthKV_ = 0.0;
  int lastIterNorth_ = 0;
  int lastIterSouth_ = 0;
  amrex::Real lastResidNorth_ = 0.0;
  amrex::Real lastResidSouth_ = 0.0;
};

#endif
