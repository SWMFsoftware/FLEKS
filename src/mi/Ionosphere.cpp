#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include "Ionosphere.h"

using namespace amrex;

namespace {

constexpr Real cDegToRad = dPI / 180.0;

inline int idx2d(int i, int j, int nPhi) { return i * nPhi + j; }

// Five-point sparse operator on a 2D (nTh x nPh) hemisphere grid with:
//   - pole average coupling at the pole neighbor row
//   - periodic boundary condition in longitude (j = 0 .. nPh - 1)
//   - Dirichlet (Phi = 0) or Neumann (dPhi/dTheta = 0) at low-latitude edge
struct FivePointOp {
  int nTh = 0;
  int nPh = 0;
  bool isNorth = true;
  bool isNeumann = false;
  // Stencil coefficients for unknown k = j * nTh + iLoc:
  //   (A x)_k = cA*x(i,j) + cB*x(i-1,j) + cC*x(i+1,j) + cD*x(i,j-1) +
  //   cE*x(i,j+1)
  std::vector<Real> cA, cB, cC, cD, cE;

  // ILU(0) factors stored in-place for the interior tridiagonal+2-band part
  // (ordered with iLoc fastest, j slowest, matching Ridley_serial prehepta):
  std::vector<Real> dInv, eL, fU, e1L, f1U;

  void build_ilu0() {
    const int n = nTh * nPh;
    dInv.assign(n, 0.0);
    eL.assign(n, 0.0);
    fU.assign(n, 0.0);
    e1L.assign(n, 0.0);
    f1U.assign(n, 0.0);

    // Copy non-periodic interior bands (boundary wrap terms handled in matvec).
    for (int j = 0; j < nPh; ++j) {
      for (int i = 0; i < nTh; ++i) {
        const int k = j * nTh + i;
        dInv[k] = cA[k];
        eL[k] = (i > 0) ? cB[k] : 0.0;
        fU[k] = (i < nTh - 1) ? cC[k] : 0.0;
        e1L[k] = (j > 0) ? cD[k] : 0.0;
        f1U[k] = (j < nPh - 1) ? cE[k] : 0.0;
      }
    }

    // Standard ILU(0) factorization for 5-point stencil (same structure as
    // SWMF/IE prehepta with omega = 0).
    for (int k = 0; k < n; ++k) {
      Real pivot = dInv[k];
      if (k >= 1)
        pivot -= eL[k] * fU[k - 1] * dInv[k - 1];
      if (k >= nTh)
        pivot -= e1L[k] * f1U[k - nTh] * dInv[k - nTh];

      if (std::abs(pivot) < Real(1.0e-30))
        pivot = (pivot >= 0.0 ? Real(1.0e-30) : Real(-1.0e-30));
      dInv[k] = Real(1.0) / pivot;
    }
  }

  // Apply M^{-1} = U^{-1} L^{-1} to vector r, writing result to z.
  void apply_precond(const std::vector<Real> &r, std::vector<Real> &z) const {
    const int n = nTh * nPh;
    // Forward solve: L y = r  (unit lower triangular with eL, e1L scaled by
    // dInv)
    for (int k = 0; k < n; ++k) {
      Real val = r[k];
      if (k >= 1)
        val -= eL[k] * z[k - 1];
      if (k >= nTh)
        val -= e1L[k] * z[k - nTh];
      z[k] = val * dInv[k];
    }
    // Backward solve: U z = y
    for (int k = n - 1; k >= 0; --k) {
      Real val = 0.0;
      if (k + 1 < n)
        val += fU[k] * z[k + 1];
      if (k + nTh < n)
        val += f1U[k] * z[k + nTh];
      z[k] -= dInv[k] * val;
    }
  }

  // Apply full operator y = A x (including pole average and phi periodicity).
  void matvec(const std::vector<Real> &x, std::vector<Real> &y) const {
    // Compute pole average from the row adjacent to the pole.
    Real poleAvg = 0.0;
    const int iPoleNeighbor = isNorth ? 0 : (nTh - 1);
    for (int j = 0; j < nPh; ++j) {
      poleAvg += x[j * nTh + iPoleNeighbor];
    }
    poleAvg /= static_cast<Real>(nPh);

    for (int j = 0; j < nPh; ++j) {
      const int jPrev = (j == 0) ? (nPh - 1) : (j - 1);
      const int jNext = (j == nPh - 1) ? 0 : (j + 1);
      const int colOff = j * nTh;
      const int prevOff = jPrev * nTh;
      const int nextOff = jNext * nTh;

      for (int i = 0; i < nTh; ++i) {
        const int k = colOff + i;
        Real xIm1 = 0.0;
        Real xIp1 = 0.0;

        if (i > 0) {
          xIm1 = x[colOff + i - 1];
        } else {
          // Top edge of local slice: pole for North, low-lat boundary for South
          if (isNorth) {
            xIm1 = poleAvg;
          } else {
            xIm1 = isNeumann ? x[colOff + i] : Real(0.0);
          }
        }

        if (i < nTh - 1) {
          xIp1 = x[colOff + i + 1];
        } else {
          // Bottom edge of local slice: low-lat for North, pole for South
          if (isNorth) {
            xIp1 = isNeumann ? x[colOff + i] : Real(0.0);
          } else {
            xIp1 = poleAvg;
          }
        }

        y[k] = cA[k] * x[k] + cB[k] * xIm1 + cC[k] * xIp1 +
               cD[k] * x[prevOff + i] + cE[k] * x[nextOff + i];
      }
    }

    // For pure Neumann low-lat BC, pin one gauge point to remove the nullspace.
    if (isNeumann && nTh > 0 && nPh > 0) {
      const int kPin = isNorth ? (nTh - 1) : 0;
      y[kPin] = x[kPin];
    }
  }
};

Real dot_vec(const std::vector<Real> &a, const std::vector<Real> &b) {
  Real s = 0.0;
  const size_t n = a.size();
  for (size_t i = 0; i < n; ++i)
    s += a[i] * b[i];
  return s;
}

Real norm2_vec(const std::vector<Real> &a) { return std::sqrt(dot_vec(a, a)); }

// Right-preconditioned BiCGSTAB for A x = b.
bool solve_bicgstab(const FivePointOp &op, const std::vector<Real> &b,
                    std::vector<Real> &x, Real tol, int maxIter, int &nIterOut,
                    Real &residOut) {
  const size_t n = b.size();
  std::vector<Real> r(n), rHat(n), p(n, 0.0), v(n, 0.0), s(n), t(n);
  std::vector<Real> y(n), z(n);

  const Real bNorm = norm2_vec(b);
  if (bNorm <= Real(1.0e-30)) {
    std::fill(x.begin(), x.end(), 0.0);
    nIterOut = 0;
    residOut = 0.0;
    return true;
  }

  op.matvec(x, v);
  for (size_t i = 0; i < n; ++i) {
    r[i] = b[i] - v[i];
    rHat[i] = r[i];
  }

  Real rNorm = norm2_vec(r);
  residOut = rNorm / bNorm;
  if (residOut <= tol) {
    nIterOut = 0;
    return true;
  }

  Real rhoOld = 1.0;
  Real alpha = 1.0;
  Real omega = 1.0;

  for (int iter = 1; iter <= maxIter; ++iter) {
    const Real rhoNew = dot_vec(rHat, r);
    if (std::abs(rhoNew) < Real(1.0e-30)) {
      nIterOut = iter;
      residOut = rNorm / bNorm;
      return residOut <= tol;
    }

    if (iter == 1) {
      p = r;
    } else {
      const Real beta = (rhoNew / rhoOld) * (alpha / omega);
      for (size_t i = 0; i < n; ++i) {
        p[i] = r[i] + beta * (p[i] - omega * v[i]);
      }
    }

    op.apply_precond(p, y);
    op.matvec(y, v);

    const Real rHatV = dot_vec(rHat, v);
    if (std::abs(rHatV) < Real(1.0e-30)) {
      nIterOut = iter;
      residOut = rNorm / bNorm;
      return residOut <= tol;
    }
    alpha = rhoNew / rHatV;

    for (size_t i = 0; i < n; ++i)
      s[i] = r[i] - alpha * v[i];

    const Real sNorm = norm2_vec(s);
    if (sNorm / bNorm <= tol) {
      for (size_t i = 0; i < n; ++i)
        x[i] += alpha * y[i];
      nIterOut = iter;
      residOut = sNorm / bNorm;
      return true;
    }

    op.apply_precond(s, z);
    op.matvec(z, t);

    const Real tt = dot_vec(t, t);
    if (tt < Real(1.0e-30)) {
      for (size_t i = 0; i < n; ++i)
        x[i] += alpha * y[i];
      nIterOut = iter;
      residOut = sNorm / bNorm;
      return residOut <= tol;
    }
    omega = dot_vec(t, s) / tt;

    for (size_t i = 0; i < n; ++i) {
      x[i] += alpha * y[i] + omega * z[i];
      r[i] = s[i] - omega * t[i];
    }

    rNorm = norm2_vec(r);
    residOut = rNorm / bNorm;
    nIterOut = iter;
    if (residOut <= tol)
      return true;

    if (std::abs(omega) < Real(1.0e-30))
      return false;
    rhoOld = rhoNew;
  }

  return residOut <= tol;
}

} // namespace

//==========================================================
void Ionosphere::build_dipole_frame(const Real dipoleMoment[3], Real rot[3][3],
                                    Real &dipoleSign) {
  const Real mNorm = std::sqrt(dipoleMoment[0] * dipoleMoment[0] +
                               dipoleMoment[1] * dipoleMoment[1] +
                               dipoleMoment[2] * dipoleMoment[2]);
  if (mNorm <= Real(1.0e-30)) {
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        rot[a][b] = (a == b) ? Real(1.0) : Real(0.0);
    dipoleSign = -1.0;
    return;
  }

  // Unit dipole moment vector.
  Real ax[3] = {
    dipoleMoment[0] / mNorm,
    dipoleMoment[1] / mNorm,
    dipoleMoment[2] / mNorm,
  };

  // Orient +e_z toward northern hemisphere (e_z.z >= 0).
  dipoleSign = (ax[2] < Real(0.0)) ? Real(-1.0) : Real(1.0);
  Real ez[3] = { dipoleSign * ax[0], dipoleSign * ax[1], dipoleSign * ax[2] };

  // Build +e_x by projecting simulation +x (sunward) orthogonal to e_z.
  // Fall back to +y if the dipole lies almost along the x-axis.
  Real ref[3] = { 1.0, 0.0, 0.0 };
  if (std::abs(ez[0]) > Real(0.99)) {
    ref[0] = 0.0;
    ref[1] = 1.0;
    ref[2] = 0.0;
  }
  const Real dotRef = ref[0] * ez[0] + ref[1] * ez[1] + ref[2] * ez[2];
  Real ex[3] = {
    ref[0] - dotRef * ez[0],
    ref[1] - dotRef * ez[1],
    ref[2] - dotRef * ez[2],
  };
  const Real exNorm = std::sqrt(ex[0] * ex[0] + ex[1] * ex[1] + ex[2] * ex[2]);
  ex[0] /= exNorm;
  ex[1] /= exNorm;
  ex[2] /= exNorm;

  // e_y = e_z x e_x (right-handed).
  Real ey[3] = {
    ez[1] * ex[2] - ez[2] * ex[1],
    ez[2] * ex[0] - ez[0] * ex[2],
    ez[0] * ex[1] - ez[1] * ex[0],
  };

  for (int b = 0; b < 3; ++b) {
    rot[0][b] = ex[b];
    rot[1][b] = ey[b];
    rot[2][b] = ez[b];
  }
}

//==========================================================
void Ionosphere::init(const Real dipoleMomentCode[3],
                      const Real bodyCenterCode[3], Real rPlanetSI, Real si2NoL,
                      Real si2NoE, Real si2NoB, Real bEqSI) {
  (void)si2NoB;
  if (params_.nThetaHem < 5)
    Abort("Ionosphere::init: nThetaHem must be >= 5.");
  if (params_.nPhi < 5)
    Abort("Ionosphere::init: nPhi must be >= 5.");
  if (rPlanetSI <= 0.0)
    Abort("Ionosphere::init: rPlanetSI must be positive.");

  rPlanetSI_ = rPlanetSI;
  bEqSI_ = bEqSI;

  for (int d = 0; d < 3; ++d)
    frame_.center[d] = bodyCenterCode[d];

  build_dipole_frame(dipoleMomentCode, frame_.rot, frame_.dipoleSign);

  frame_.rIonoSI = params_.rIonoRatio * rPlanetSI;
  frame_.rIonoCode = frame_.rIonoSI * si2NoL;
  frame_.si2NoE = si2NoE;
  frame_.no2SiL = Real(1.0) / si2NoL;

  // Full pole-to-pole grid: nTheta = 2*nThetaHem - 1 (north pole = 0,
  // equator = nThetaHem - 1, south pole = nTheta - 1).
  nTheta_ = 2 * params_.nThetaHem - 1;
  nPhi_ = params_.nPhi;
  dTheta_ = dPI / static_cast<Real>(nTheta_ - 1);
  dPhi_ = Real(2.0) * dPI / static_cast<Real>(nPhi_ - 1);
  latBoundaryRad_ = params_.latBoundaryDeg * cDegToRad;

  theta_.resize(nTheta_);
  sinTheta_.resize(nTheta_);
  cosTheta_.resize(nTheta_);
  for (int i = 0; i < nTheta_; ++i) {
    theta_[i] = i * dTheta_;
    sinTheta_[i] = std::sin(theta_[i]);
    cosTheta_[i] = std::cos(theta_[i]);
  }

  phi_.resize(nPhi_);
  for (int j = 0; j < nPhi_; ++j) {
    phi_[j] = j * dPhi_;
  }

  const size_t nTot = static_cast<size_t>(nTheta_) * nPhi_;
  jr_.assign(nTot, 0.0);
  sigmaP_.assign(nTot, params_.sigmaPedersen);
  sigmaH_.assign(nTot, params_.sigmaHall);
  sigmaThTh_.assign(nTot, 0.0);
  sigmaThPh_.assign(nTot, 0.0);
  sigmaPhPh_.assign(nTot, 0.0);
  dSigmaThTh_dTh_.assign(nTot, 0.0);
  dSigmaThPh_dTh_.assign(nTot, 0.0);
  dSigmaThPh_dPh_.assign(nTot, 0.0);
  dSigmaPhPh_dPh_.assign(nTot, 0.0);

  phiPot_.assign(nTot, 0.0);
  phiSolved_.assign(nTot, 0.0);
  dPhiDTheta_.assign(nTot, 0.0);
  dPhiDPhi_.assign(nTot, 0.0);
  hasPrevSolution_ = false;

  compute_conductances();
  compute_conductance_tensor();
  sync_device_gradients();

  initialized_ = true;
}

//==========================================================
void Ionosphere::compute_conductances() {
  const bool useSolar = (params_.conductanceModel == "solar_euv" ||
                         params_.conductanceModel == "auroral");
  const bool useAuroral = (params_.conductanceModel == "auroral");

  // Moen & Brekke (1993) solar EUV conductance coefficients.
  const Real f107Term = std::pow(std::max(params_.f107, Real(1.0)), Real(0.49));
  const Real f107Hall = std::pow(std::max(params_.f107, Real(1.0)), Real(0.53));

  for (int i = 0; i < nTheta_; ++i) {
    const Real sn = sinTheta_[i];
    for (int j = 0; j < nPhi_; ++j) {
      const int k = idx2d(i, j, nPhi_);

      // Base nightside / uniform conductance.
      Real sp = params_.sigmaPedersen;
      Real sh = params_.sigmaHall;

      if (useSolar) {
        // In SMG coordinates, +x (phi = 0) points sunward (ignoring seasonal
        // dipole tilt in the simplest approximation; cos(chi) =
        // sin(theta)*cos(phi) rotated back to sim x-axis).
        const Real xSmg = sn * std::cos(phi_[j]);
        const Real ySmg = sn * std::sin(phi_[j]);
        const Real zSmg = cosTheta_[i];
        // Sun direction in simulation frame is +x => cos(SZA) = (R^T * r_smg).x
        const Real cosChi = frame_.rot[0][0] * xSmg + frame_.rot[1][0] * ySmg +
                            frame_.rot[2][0] * zSmg;

        if (cosChi > Real(0.0)) {
          const Real spEuv =
              f107Term * (Real(0.34) * cosChi + Real(0.93) * std::sqrt(cosChi));
          const Real shEuv =
              f107Hall * (Real(0.81) * cosChi + Real(0.54) * std::sqrt(cosChi));
          sp = std::sqrt(sp * sp + spEuv * spEuv);
          sh = std::sqrt(sh * sh + shEuv * shEuv);
        }
      }

      if (useAuroral && params_.auroralFacScale > Real(0.0)) {
        // Upward FAC (J_r^{outward} > 0) corresponds to precipitating electrons
        // enhancing auroral Pedersen and Hall conductances.
        const Real jUp = std::max(Real(0.0), jr_[k]);
        if (jUp > Real(0.0)) {
          const Real ratio = jUp / params_.auroralFacScale;
          const Real spAur =
              params_.auroralSigmaPMax * std::tanh(std::sqrt(ratio));
          const Real shAur = params_.auroralHallRatio * spAur;
          sp = std::sqrt(sp * sp + spAur * spAur);
          sh = std::sqrt(sh * sh + shAur * shAur);
        }
      }

      sigmaP_[k] = std::max(sp, Real(0.1));
      sigmaH_[k] = sh;
    }
  }
}

//==========================================================
void Ionosphere::compute_conductance_tensor() {
  // Exact dipole dip-angle conductance tensor matching SWMF/IE
  // (iono_conductances.f90:2411-2425):
  //   C = 4 * Sigma0 * cos^2(theta) + SigmaP * sin^2(theta)
  //   Sigma_th_th = Sigma0 * SigmaP * (1 + 3 cos^2(theta)) / C
  //   Sigma_th_ph = (-dipoleSign) * 2 * Sigma0 * SigmaH * cos(theta) *
  //                 sqrt(1 + 3 cos^2(theta)) / C
  //   Sigma_ph_ph = SigmaP + SigmaH^2 * sin^2(theta) / C
  // Note: for Earth (dipoleSign = -1), -dipoleSign = +1, matching
  // Ridley_serial.
  const Real s0 = std::max(params_.sigmaParallel, Real(10.0));
  const Real hallSign = -frame_.dipoleSign;

  for (int i = 0; i < nTheta_; ++i) {
    const Real sn = sinTheta_[i];
    const Real cs = cosTheta_[i];
    const Real sn2 = sn * sn;
    const Real cs2 = cs * cs;
    const Real cs3 = Real(1.0) + Real(3.0) * cs2;
    const Real cs4 = std::sqrt(cs3);

    for (int j = 0; j < nPhi_; ++j) {
      const int k = idx2d(i, j, nPhi_);
      const Real sp = sigmaP_[k];
      const Real sh = sigmaH_[k];
      const Real C = Real(4.0) * s0 * cs2 + sp * sn2;

      sigmaThTh_[k] = s0 * sp * cs3 / C;
      sigmaThPh_[k] = hallSign * Real(2.0) * s0 * sh * cs * cs4 / C;
      sigmaPhPh_[k] = sp + (sh * sh * sn2) / C;
    }
  }

  // Central differences of tensor components on interior points (periodic in
  // phi).
  const Real inv2dTh = Real(0.5) / dTheta_;
  const Real inv2dPh = Real(0.5) / dPhi_;
  const int nPhiUsed = nPhi_ - 1;

  for (int i = 0; i < nTheta_; ++i) {
    for (int j = 0; j < nPhi_; ++j) {
      const int k = idx2d(i, j, nPhi_);
      if (i > 0 && i < nTheta_ - 1) {
        const int km = idx2d(i - 1, j, nPhi_);
        const int kp = idx2d(i + 1, j, nPhi_);
        dSigmaThTh_dTh_[k] = (sigmaThTh_[kp] - sigmaThTh_[km]) * inv2dTh;
        dSigmaThPh_dTh_[k] = (sigmaThPh_[kp] - sigmaThPh_[km]) * inv2dTh;
      } else {
        dSigmaThTh_dTh_[k] = 0.0;
        dSigmaThPh_dTh_[k] = 0.0;
      }

      const int jUsed = j % nPhiUsed;
      const int jPrev = (jUsed == 0) ? (nPhiUsed - 1) : (jUsed - 1);
      const int jNext = (jUsed == nPhiUsed - 1) ? 0 : (jUsed + 1);
      const int kPrev = idx2d(i, jPrev, nPhi_);
      const int kNext = idx2d(i, jNext, nPhi_);
      dSigmaThPh_dPh_[k] = (sigmaThPh_[kNext] - sigmaThPh_[kPrev]) * inv2dPh;
      dSigmaPhPh_dPh_[k] = (sigmaPhPh_[kNext] - sigmaPhPh_[kPrev]) * inv2dPh;
    }
  }
}

//==========================================================
void Ionosphere::smooth_jr() {
  if (params_.nSmoothFac <= 0)
    return;

  const int nPhiUsed = nPhi_ - 1;
  std::vector<Real> tmp(jr_.size(), 0.0);

  for (int pass = 0; pass < params_.nSmoothFac; ++pass) {
    // 5-point binomial smoothing (0.5 center + 0.125 * 4 neighbors) on the
    // spherical grid, preserving periodicity in phi and pole averages.
    for (int i = 1; i < nTheta_ - 1; ++i) {
      for (int j = 0; j < nPhiUsed; ++j) {
        const int jPrev = (j == 0) ? (nPhiUsed - 1) : (j - 1);
        const int jNext = (j == nPhiUsed - 1) ? 0 : (j + 1);
        tmp[idx2d(i, j, nPhi_)] =
            Real(0.5) * jr_[idx2d(i, j, nPhi_)] +
            Real(0.125) *
                (jr_[idx2d(i - 1, j, nPhi_)] + jr_[idx2d(i + 1, j, nPhi_)] +
                 jr_[idx2d(i, jPrev, nPhi_)] + jr_[idx2d(i, jNext, nPhi_)]);
      }
      tmp[idx2d(i, nPhi_ - 1, nPhi_)] = tmp[idx2d(i, 0, nPhi_)];
    }

    // Poles: average of adjacent ring.
    Real northAvg = 0.0, southAvg = 0.0;
    for (int j = 0; j < nPhiUsed; ++j) {
      northAvg += tmp[idx2d(1, j, nPhi_)];
      southAvg += tmp[idx2d(nTheta_ - 2, j, nPhi_)];
    }
    northAvg /= static_cast<Real>(nPhiUsed);
    southAvg /= static_cast<Real>(nPhiUsed);
    for (int j = 0; j < nPhi_; ++j) {
      tmp[idx2d(0, j, nPhi_)] = northAvg;
      tmp[idx2d(nTheta_ - 1, j, nPhi_)] = southAvg;
    }

    jr_.swap(tmp);
  }
}

//==========================================================
bool Ionosphere::solve_hemisphere(bool isNorth, int &nIterOut, Real &residOut) {
  const int nPhiUsed = nPhi_ - 1;
  const int iEq = params_.nThetaHem - 1;

  // Determine active colatitude index range [iMin, iMax] (inclusive) strictly
  // between the pole and the low-latitude boundary. Use a small ULP margin so
  // grid points lying exactly on latBoundaryDeg are treated as boundary nodes.
  constexpr Real epsBnd = Real(1.0e-12);
  int iMin = 1, iMax = 0;
  if (isNorth) {
    iMin = 1;
    iMax = 0;
    for (int i = 1; i < iEq; ++i) {
      if (Real(0.5) * dPI - theta_[i] > latBoundaryRad_ + epsBnd)
        iMax = i;
    }
  } else {
    iMax = nTheta_ - 2;
    iMin = nTheta_ - 1;
    for (int i = nTheta_ - 2; i > iEq; --i) {
      if (theta_[i] - Real(0.5) * dPI > latBoundaryRad_ + epsBnd)
        iMin = i;
    }
  }

  const int nTh = iMax - iMin + 1;
  if (nTh <= 0) {
    nIterOut = 0;
    residOut = 0.0;
    return true;
  }

  const bool isNeumann = (params_.lowLatBC == "neumann");
  const int nX = nTh * nPhiUsed;

  FivePointOp op;
  op.nTh = nTh;
  op.nPh = nPhiUsed;
  op.isNorth = isNorth;
  op.isNeumann = isNeumann;
  op.cA.resize(nX);
  op.cB.resize(nX);
  op.cC.resize(nX);
  op.cD.resize(nX);
  op.cE.resize(nX);

  std::vector<Real> rhs(nX, 0.0);
  std::vector<Real> x(nX, 0.0);

  const Real dTh2 = dTheta_ * dTheta_;
  const Real dPh2 = dPhi_ * dPhi_;
  const Real twoDTh = Real(2.0) * dTheta_;
  const Real twoDPh = Real(2.0) * dPhi_;
  const Real rIono2 = frame_.rIonoSI * frame_.rIonoSI;

  // For Neumann BC, enforce discrete current balance (sum of Jr*dA = 0) so
  // the singular system is consistent.
  Real jrShift = 0.0;
  if (isNeumann) {
    Real num = 0.0, den = 0.0;
    for (int j = 0; j < nPhiUsed; ++j) {
      for (int i = iMin; i <= iMax; ++i) {
        const Real w = sinTheta_[i];
        num += jr_[idx2d(i, j, nPhi_)] * w;
        den += w;
      }
    }
    if (den > Real(0.0))
      jrShift = num / den;
  }

  for (int j = 0; j < nPhiUsed; ++j) {
    for (int iLoc = 0; iLoc < nTh; ++iLoc) {
      const int i = iMin + iLoc;
      const int kGrid = idx2d(i, j, nPhi_);
      const int k = j * nTh + iLoc;

      const Real sn = sinTheta_[i];
      const Real cs = cosTheta_[i];
      const Real sn2 = sn * sn;

      const Real termTh2 = sigmaThTh_[kGrid] * sn2 / dTh2;
      const Real termTh1 =
          (sigmaThTh_[kGrid] * sn * cs + dSigmaThTh_dTh_[kGrid] * sn2 -
           dSigmaThPh_dPh_[kGrid] * sn) /
          twoDTh;
      const Real termPh2 = sigmaPhPh_[kGrid] / dPh2;
      const Real termPh1 =
          (dSigmaThPh_dTh_[kGrid] * sn + dSigmaPhPh_dPh_[kGrid]) / twoDPh;

      op.cA[k] = -Real(2.0) * (termTh2 + termPh2);
      op.cB[k] = termTh2 - termTh1;
      op.cC[k] = termTh2 + termTh1;
      op.cD[k] = termPh2 - termPh1;
      op.cE[k] = termPh2 + termPh1;

      rhs[k] = (jr_[kGrid] - jrShift) * rIono2 * sn2;
      x[k] = hasPrevSolution_ ? phiSolved_[kGrid] : Real(0.0);
    }
  }

  if (isNeumann) {
    const int kPin = isNorth ? (nTh - 1) : 0;
    rhs[kPin] = 0.0;
    op.cA[kPin] = 1.0;
    op.cB[kPin] = op.cC[kPin] = op.cD[kPin] = op.cE[kPin] = 0.0;
  }

  op.build_ilu0();
  const bool converged = solve_bicgstab(
      op, rhs, x, params_.solverTol, params_.solverMaxIter, nIterOut, residOut);

  // Scatter solution into phiSolved_ (zeroing low-latitude region).
  const int iStartClear = isNorth ? 0 : iEq;
  const int iEndClear = isNorth ? iEq : (nTheta_ - 1);
  for (int i = iStartClear; i <= iEndClear; ++i) {
    for (int j = 0; j < nPhi_; ++j) {
      phiSolved_[idx2d(i, j, nPhi_)] = 0.0;
    }
  }

  Real phiMin = x[0], phiMax = x[0];
  for (int j = 0; j < nPhiUsed; ++j) {
    for (int iLoc = 0; iLoc < nTh; ++iLoc) {
      const int i = iMin + iLoc;
      const Real val = x[j * nTh + iLoc];
      phiSolved_[idx2d(i, j, nPhi_)] = val;
      if (val < phiMin)
        phiMin = val;
      if (val > phiMax)
        phiMax = val;
    }
  }

  // Pole average.
  const int iPole = isNorth ? 0 : (nTheta_ - 1);
  const int iAdj = isNorth ? 1 : (nTheta_ - 2);
  Real poleVal = 0.0;
  for (int j = 0; j < nPhiUsed; ++j)
    poleVal += phiSolved_[idx2d(iAdj, j, nPhi_)];
  poleVal /= static_cast<Real>(nPhiUsed);
  for (int j = 0; j < nPhi_; ++j)
    phiSolved_[idx2d(iPole, j, nPhi_)] = poleVal;

  // Periodic wrap in longitude.
  for (int i = iStartClear; i <= iEndClear; ++i) {
    phiSolved_[idx2d(i, nPhi_ - 1, nPhi_)] = phiSolved_[idx2d(i, 0, nPhi_)];
  }

  if (isNorth) {
    cpcpNorthKV_ = (phiMax - phiMin) * Real(1.0e-3);
  } else {
    cpcpSouthKV_ = (phiMax - phiMin) * Real(1.0e-3);
  }

  return converged;
}

//==========================================================
void Ionosphere::add_corotation_potential() {
  if (!params_.useCorotation)
    return;

  // Corotation potential on a dipole shell at r = rIonoSI:
  //   E_corot = -(Omega x r) x B = -grad(Phi_corot)
  //   Phi_corot(theta) = -Omega * B_eq * R_p^3 * sin^2(theta) / r_I
  // (with B_eq < 0 for Earth, Phi_corot > 0 toward equator, E_r inward).
  const Real coef = -params_.omegaPlanetSI * bEqSI_ *
                    (rPlanetSI_ * rPlanetSI_ * rPlanetSI_) / frame_.rIonoSI;
  for (int i = 0; i < nTheta_; ++i) {
    const Real dPhiCorot = coef * sinTheta_[i] * sinTheta_[i];
    for (int j = 0; j < nPhi_; ++j) {
      phiPot_[idx2d(i, j, nPhi_)] += dPhiCorot;
    }
  }
}

//==========================================================
void Ionosphere::compute_gradients() {
  const int nPhiUsed = nPhi_ - 1;
  const Real inv2dTh = Real(0.5) / dTheta_;
  const Real inv2dPh = Real(0.5) / dPhi_;

  // Interior central differences in theta.
  for (int i = 1; i < nTheta_ - 1; ++i) {
    for (int j = 0; j < nPhi_; ++j) {
      dPhiDTheta_[idx2d(i, j, nPhi_)] =
          (phiPot_[idx2d(i + 1, j, nPhi_)] - phiPot_[idx2d(i - 1, j, nPhi_)]) *
          inv2dTh;
    }
  }

  // Second-order one-sided differences at the poles (matching BATSRUS
  // ModIeCoupling.f90:262-272).
  for (int j = 0; j < nPhi_; ++j) {
    dPhiDTheta_[idx2d(0, j, nPhi_)] = (Real(4.0) * phiPot_[idx2d(1, j, nPhi_)] -
                                       Real(3.0) * phiPot_[idx2d(0, j, nPhi_)] -
                                       phiPot_[idx2d(2, j, nPhi_)]) *
                                      inv2dTh;
    dPhiDTheta_[idx2d(nTheta_ - 1, j, nPhi_)] =
        (Real(3.0) * phiPot_[idx2d(nTheta_ - 1, j, nPhi_)] -
         Real(4.0) * phiPot_[idx2d(nTheta_ - 2, j, nPhi_)] +
         phiPot_[idx2d(nTheta_ - 3, j, nPhi_)]) *
        inv2dTh;
  }

  // Periodic central differences in phi.
  for (int i = 0; i < nTheta_; ++i) {
    for (int j = 0; j < nPhiUsed; ++j) {
      const int jPrev = (j == 0) ? (nPhiUsed - 1) : (j - 1);
      const int jNext = (j == nPhiUsed - 1) ? 0 : (j + 1);
      dPhiDPhi_[idx2d(i, j, nPhi_)] =
          (phiPot_[idx2d(i, jNext, nPhi_)] - phiPot_[idx2d(i, jPrev, nPhi_)]) *
          inv2dPh;
    }
    dPhiDPhi_[idx2d(i, nPhi_ - 1, nPhi_)] = dPhiDPhi_[idx2d(i, 0, nPhi_)];
  }
}

//==========================================================
void Ionosphere::sync_device_gradients() {
  const size_t nTot = dPhiDTheta_.size();
  dPhiDThetaDev_.resize(nTot);
  dPhiDPhiDev_.resize(nTot);
  if (nTot > 0) {
    Gpu::copyAsync(Gpu::hostToDevice, dPhiDTheta_.begin(), dPhiDTheta_.end(),
                   dPhiDThetaDev_.begin());
    Gpu::copyAsync(Gpu::hostToDevice, dPhiDPhi_.begin(), dPhiDPhi_.end(),
                   dPhiDPhiDev_.begin());
    Gpu::streamSynchronize();
  }
}

//==========================================================
bool Ionosphere::solve(const std::vector<Real> &jrOutwardSI, bool doBroadcast) {
  if (!initialized_)
    Abort("Ionosphere::solve called before init().");

  const size_t nTot = static_cast<size_t>(nTheta_) * nPhi_;
  if (jrOutwardSI.size() != nTot)
    Abort("Ionosphere::solve: jrOutwardSI size mismatch.");

  int convFlag = 1;

  if (!doBroadcast || ParallelDescriptor::IOProcessor()) {
    jr_ = jrOutwardSI;
    smooth_jr();
    compute_conductances();
    compute_conductance_tensor();

    const bool convN = solve_hemisphere(true, lastIterNorth_, lastResidNorth_);
    const bool convS = solve_hemisphere(false, lastIterSouth_, lastResidSouth_);
    convFlag = (convN && convS) ? 1 : 0;

    // Temporal blending with previous solution.
    const Real alpha = std::clamp(params_.relaxAlpha, Real(0.0), Real(1.0));
    if (!hasPrevSolution_ || alpha >= Real(1.0)) {
      phiPot_ = phiSolved_;
    } else {
      for (size_t k = 0; k < nTot; ++k) {
        phiPot_[k] = (Real(1.0) - alpha) * phiPot_[k] + alpha * phiSolved_[k];
      }
    }
    hasPrevSolution_ = true;

    add_corotation_potential();
    compute_gradients();
  }

  if (doBroadcast && ParallelDescriptor::NProcs() > 1) {
    const int root = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::Bcast(phiPot_.data(), nTot, root);
    ParallelDescriptor::Bcast(dPhiDTheta_.data(), nTot, root);
    ParallelDescriptor::Bcast(dPhiDPhi_.data(), nTot, root);
    ParallelDescriptor::Bcast(&cpcpNorthKV_, 1, root);
    ParallelDescriptor::Bcast(&cpcpSouthKV_, 1, root);
    ParallelDescriptor::Bcast(&convFlag, 1, root);
    hasPrevSolution_ = true;
  }

  sync_device_gradients();
  return convFlag != 0;
}

//==========================================================
void Ionosphere::set_potential(const std::vector<Real> &phiIn) {
  const size_t nTot = static_cast<size_t>(nTheta_) * nPhi_;
  if (phiIn.size() != nTot)
    Abort("Ionosphere::set_potential: size mismatch.");
  phiPot_ = phiIn;
  phiSolved_ = phiIn;
  hasPrevSolution_ = true;
  compute_gradients();
  sync_device_gradients();
}

//==========================================================
IonosphereDeviceView Ionosphere::device_view() const {
  IonosphereDeviceView v;
  v.frame = frame_;
  v.nTheta = nTheta_;
  v.nPhi = nPhi_;
  v.dTheta = dTheta_;
  v.dPhi = dPhi_;
  v.latBoundaryRad = latBoundaryRad_;
  v.dPhiDTheta = dPhiDThetaDev_.data();
  v.dPhiDPhi = dPhiDPhiDev_.data();
  return v;
}

//==========================================================
void Ionosphere::write_diagnostics(const std::string &fileName, Real timeSI,
                                   int cycle) const {
  if (!ParallelDescriptor::IOProcessor())
    return;

  std::ofstream out(fileName);
  if (!out.is_open()) {
    Print() << "WARNING: Ionosphere::write_diagnostics could not open "
            << fileName << "\n";
    return;
  }

  out << "# FLEKS Ionosphere Diagnostics\n";
  out << "# cycle = " << cycle << "  time_SI = " << std::scientific
      << std::setprecision(6) << timeSI << "  cpcpN_kV = " << cpcpNorthKV_
      << "  cpcpS_kV = " << cpcpSouthKV_ << "\n";
  out << "# nTheta = " << nTheta_ << "  nPhi = " << nPhi_
      << "  rIono_m = " << frame_.rIonoSI << "\n";
  out << "# theta_deg phi_deg Jr_uAm2 SigmaP_S SigmaH_S Phi_kV Eth_mVm "
         "Eph_mVm\n";

  const Real invRI = Real(1.0) / frame_.rIonoSI;
  for (int i = 0; i < nTheta_; ++i) {
    const Real thDeg = theta_[i] / cDegToRad;
    const Real sn = sinTheta_[i];
    for (int j = 0; j < nPhi_; ++j) {
      const int k = idx2d(i, j, nPhi_);
      const Real phDeg = phi_[j] / cDegToRad;
      const Real eth = -dPhiDTheta_[k] * invRI * Real(1.0e3); // mV/m
      const Real eph = (sn > Real(1.0e-8))
                           ? (-dPhiDPhi_[k] * invRI / sn * Real(1.0e3))
                           : Real(0.0);
      out << std::fixed << std::setprecision(4) << thDeg << " " << phDeg << " "
          << std::scientific << std::setprecision(6) << jr_[k] * Real(1.0e6)
          << " " << sigmaP_[k] << " " << sigmaH_[k] << " "
          << phiPot_[k] * Real(1.0e-3) << " " << eth << " " << eph << "\n";
    }
  }
}
