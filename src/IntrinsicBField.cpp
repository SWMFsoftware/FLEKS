#include "IntrinsicBField.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include "Constants.h"

using namespace amrex;

namespace {
constexpr double cDegToRad = dPI / 180.0;
// The recurrence divides by sin(theta); below this the azimuthal component is
// set to zero instead, exactly as the reference implementation does.
constexpr double cThetaSinMin = 1.0e-6;
// nanotesla -> tesla, the unit of #DIPOLE / #CRUSTALFIELD input.
constexpr double cNTToT = 1.0e-9;
}

//==========================================================
void IntrinsicBField::read_dipole(Real bEqNT, Real thetaDeg, Real phiDeg,
                                  Real rRefSI) {
  useDipole_ = true;
  bEqNT_ = bEqNT;
  thetaDeg_ = thetaDeg;
  phiDeg_ = phiDeg;
  rRefSI_ = rRefSI;
}

//==========================================================
void IntrinsicBField::read_crustal(const std::string &fileName, int nMax) {
  fileName_ = fileName;
  nMax_ = nMax;
}

//==========================================================
// Read the spherical harmonic coefficients from the crustal field file,
// matching the latest layout in BATSRUS ModUserMars:
// - 3 header lines (skipped)
// - For each degree i = 1..n:
//     one 'n m value' triple for g(i,0) followed by a g/h pair for each m = 1..i.
//
// Only rank 0 touches the file; the result is broadcast so that every rank
// holds an identical copy.
//==========================================================
void IntrinsicBField::read_coefficients() {
  const int n = nMax_;
  const long nCoef = static_cast<long>(n) * n;

  std::vector<Real> c(nCoef, 0.0), d(nCoef, 0.0);
  int nRead = 0;

  if (ParallelDescriptor::IOProcessor()) {
    std::ifstream inFile(fileName_);
    if (!inFile.is_open()) {
      Abort("IntrinsicBField: cannot open the crustal field file '" +
            fileName_ + "'.");
    }

    std::string line;
    // Skip the three header lines.
    for (int i = 0; i < 3 && std::getline(inFile, line); ++i) {
    }
    // Then, per degree i = 1..n, one 'n m value' triple for g(i,0) followed
    // by a g/h pair for every m = 1..i.
    for (int iRow = 1; iRow <= n; ++iRow) {
      int nn = -1, mm = -1;
      Real val = 0.0;
      if (!(inFile >> nn >> mm >> val))
        break;
      // Degrees 1..n-1 feed the evaluation, which loops n = 0..n-1 exactly
      // like the reference implementation; the last row of the file is
      // therefore skipped, again like the reference.
      if (nn >= 1 && nn < n && mm == 0)
        c[nn * n + 0] = val;
      for (int m = 1; m <= iRow; ++m) {
        if (!(inFile >> nn >> mm >> val))
          break;
        if (nn >= 1 && nn < n && mm >= 0 && mm < n)
          c[nn * n + mm] = val;
        if (!(inFile >> nn >> mm >> val))
          break;
        if (nn >= 1 && nn < n && mm >= 0 && mm < n)
          d[nn * n + mm] = val;
      }
      ++nRead;
    }
    inFile.close();
  }

  ParallelDescriptor::Bcast(&nRead, 1, ParallelDescriptor::IOProcessorNumber());
  if (nRead == 0) {
    Abort("IntrinsicBField: no coefficient was read from '" + fileName_ +
          "'. Check nMax against the file layout.");
  }

  coefC_.resize(nCoef);
  coefD_.resize(nCoef);
  ParallelDescriptor::Bcast(c.data(), nCoef,
                            ParallelDescriptor::IOProcessorNumber());
  ParallelDescriptor::Bcast(d.data(), nCoef,
                            ParallelDescriptor::IOProcessorNumber());
  for (long i = 0; i < nCoef; ++i) {
    coefC_[i] = c[i];
    coefD_[i] = d[i];
  }
}

//==========================================================
// The expansion is written for r measured in planetary radii and returns nT.
// Folding rPlanet^(n+2) and the nT -> code conversion into the coefficients
// makes the evaluation a plain function of the code-unit radius returning
// code-unit field, so no unit is touched at run time.
//==========================================================
void IntrinsicBField::rescale_coefficients(double si2NoB) {
  const int n = nMax_;
  const Real unit = static_cast<Real>(cNTToT * si2NoB);

  Real power = static_cast<Real>(rPlanetCode_ * rPlanetCode_); // n+2 with n = 0
  for (int iDeg = 0; iDeg < n; ++iDeg) {
    const Real scale = power * unit;
    for (int m = 0; m < n; ++m) {
      coefC_[iDeg * n + m] *= scale;
      coefD_[iDeg * n + m] *= scale;
    }
    power *= static_cast<Real>(rPlanetCode_);
  }
}

//==========================================================
void IntrinsicBField::convert_units(double si2NoB, double si2NoL,
                                    double rPlanetSI,
                                    const double *bodyCenterCode,
                                    double bodyRadiusCode, bool useBody,
                                    int nDim) {
  nDim_ = nDim;
  rPlanetCode_ = rPlanetSI * si2NoL;
  if (!(rPlanetCode_ > 0.0)) {
    Abort("IntrinsicBField: the planetary radius is not positive; check "
          "#PLANETRADIUS and #NORMALIZATION.");
  }

  for (int i = 0; i < 3; ++i)
    params_.center[i] = useBody ? static_cast<Real>(bodyCenterCode[i]) : 0.0;

  // Dipole: the reference radius is #DIPOLE's own radius when given, else the
  // radius of the inner body, else the planetary radius.
  if (rRefSI_ > 0.0)
    rRefCode_ = rRefSI_ * si2NoL;
  else if (useBody)
    rRefCode_ = bodyRadiusCode;
  else
    rRefCode_ = rPlanetCode_;

  bEqCode_ = bEqNT_ * cNTToT * si2NoB;
  if (useDipole_) {
    // The axis is the spherical direction (theta, phi) measured from +z,
    // tipped towards -x at phi = 0 so that it agrees with the BATSRUS dipole
    // tilt (Dipole_D = Bdp * [-sin(theta), 0, cos(theta)]).
    const Real theta = thetaDeg_ * cDegToRad;
    const Real phi = phiDeg_ * cDegToRad;
    const Real sinTheta = std::sin(theta);
    const Real moment =
        static_cast<Real>(bEqCode_ * rRefCode_ * rRefCode_ * rRefCode_);

    params_.moment[0] = -moment * sinTheta * std::cos(phi);
    params_.moment[1] = -moment * sinTheta * std::sin(phi);
    params_.moment[2] = moment * std::cos(theta);
  }

  if (nMax_ > 0) {
    read_coefficients();
    rescale_coefficients(si2NoB);
    params_.coefC = coefC_.data();
    params_.coefD = coefD_.data();
    params_.nMax = nMax_;
  }
}

//==========================================================
// Schmidt semi-normalized spherical harmonics and the (Br, Btheta, Bphi) of
// the crustal part, for a position given in code units. This is a direct port
// of ModUserMars::set_mars_b0, adapted to take r in code units (the powers of
// the planetary radius are already folded into the coefficients).
//==========================================================
void IntrinsicBField::eval_crustal(Real x, Real y, Real z, Real b[3]) const {
  const int n = nMax_;
  const int NN = n - 1;
  if (NN < 0 || params_.coefC == nullptr || params_.coefD == nullptr)
    return;

  const Real xc = x - params_.center[0];
  const Real yc = y - params_.center[1];
  const Real zc = (nDim_ > 2) ? (z - params_.center[2]) : 0.0;

  const Real r = std::sqrt(xc * xc + yc * yc + zc * zc);
  if (r <= 0.0)
    return;

  // Colatitude and longitude of the point; theta is in [0, pi] so xtsin >= 0.
  const Real xtcos = zc / r;
  const Real xtsin = std::sqrt(std::max(0.0, 1.0 - xtcos * xtcos));
  const Real phi = std::atan2(yc, xc);

  // Scratch, grown once and then reused: the fill loop runs on the host, so a
  // per-point allocation would dominate the cost.
  const int stride = NN + 2;
  const std::size_t nRnm = static_cast<std::size_t>(stride) * stride;
  const std::size_t nAorn = static_cast<std::size_t>(NN) + 3;
  static thread_local std::vector<Real> work;
  const std::size_t need = nRnm + nAorn + 2 * static_cast<std::size_t>(n);
  if (work.size() < need)
    work.resize(need);

  Real *const Rnm = work.data();
  Real *const aorn = Rnm + nRnm;
  Real *const xpcos = aorn + nAorn;
  Real *const xpsin = xpcos + n;

  for (int im = 0; im <= NN; ++im) {
    xpcos[im] = std::cos(im * phi);
    xpsin[im] = std::sin(im * phi);
  }

  const Real invR = 1.0 / r;
  aorn[0] = 1.0;
  for (int k = 1; k <= NN + 2; ++k)
    aorn[k] = invR * aorn[k - 1];

  auto RN = [&](const int i, const int j) -> Real & {
    return Rnm[static_cast<std::size_t>(i) * stride + j];
  };
  for (std::size_t i = 0; i < nRnm; ++i)
    Rnm[i] = 0.0;

  RN(0, 0) = 1.0;
  RN(1, 0) = xtcos;
  for (int nn = 1; nn <= NN; ++nn) {
    if (nn == 1)
      RN(nn, nn) = xtsin * RN(nn - 1, nn - 1);
    else
      RN(nn, nn) = std::sqrt((nn - 0.5) / nn) * xtsin * RN(nn - 1, nn - 1);
    RN(nn + 1, nn) = xtcos * std::sqrt(2.0 * nn + 1.0) * RN(nn, nn);
  }
  for (int m = 0; m <= NN; ++m) {
    for (int l = m + 2; l <= NN; ++l) {
      RN(l, m) = (xtcos * (2.0 * l - 1.0) * RN(l - 1, m) -
                  RN(l - 2, m) * std::sqrt((l + m - 1.0) * (l - m - 1.0))) /
                 std::sqrt(1.0 * l * l - 1.0 * m * m);
    }
  }

  Real bsph[3] = {0.0, 0.0, 0.0};
  for (int m = 0; m <= NN; ++m) {
    for (int nn = m; nn <= NN; ++nn) {
      Real dRnm;
      if (m == 0) {
        dRnm = -std::sqrt((nn + 1.0) * nn / 2.0) * RN(nn, m + 1);
      } else if (xtsin <= cThetaSinMin) {
        dRnm = -std::sqrt((nn + m + 1.0) * (nn - m)) * RN(nn, m + 1);
      } else {
        dRnm = m * xtcos * RN(nn, m) / xtsin -
               std::sqrt((nn + m + 1.0) * (nn - m)) * RN(nn, m + 1);
      }

      const Real cc = params_.coefC[nn * n + m];
      const Real dd = params_.coefD[nn * n + m];
      const Real cd = cc * xpcos[m] + dd * xpsin[m];

      bsph[0] += (nn + 1) * aorn[nn + 2] * RN(nn, m) * cd;
      bsph[1] -= aorn[nn + 2] * dRnm * cd;
      if (xtsin > cThetaSinMin)
        bsph[2] -=
            aorn[nn + 2] * RN(nn, m) * m / xtsin * (-cc * xpsin[m] + dd * xpcos[m]);
    }
  }

  // Spherical -> Cartesian with the rows of the rotation matrix being the unit
  // vectors r^, theta^, phi^ (this is BATSRUS rot_xyz_sph).
  const Real cp = xpcos[1];
  const Real sp = xpsin[1];
  const Real st = xtsin;
  const Real ct = xtcos;

  b[0] += bsph[0] * (st * cp) + bsph[1] * (ct * cp) + bsph[2] * (-sp);
  b[1] += bsph[0] * (st * sp) + bsph[1] * (ct * sp) + bsph[2] * cp;
  b[2] += bsph[0] * ct + bsph[1] * (-st);
}

//==========================================================
void IntrinsicBField::eval(Real x, Real y, Real z, Real b[3]) const {
  const Real zEval = (nDim_ > 2) ? z : params_.center[2];
  eval_dipole_b(params_, x, y, zEval, b);

  if (nMax_ <= 0)
    return;

  eval_crustal(x, y, z, b);
}

//==========================================================
std::string IntrinsicBField::describe() const {
  std::ostringstream os;
  os << "  intrinsic magnetic field B0 (static, B = B1 + B0):\n";
  os << "    center [code]         = (" << params_.center[0] << ", "
     << params_.center[1] << ", " << params_.center[2] << ")\n";
  if (useDipole_) {
    os << "    dipole strength       = " << bEqNT_ << " [nT] at rRef\n";
    os << "    dipole tilt           = theta " << thetaDeg_ << ", phi "
       << phiDeg_ << " [deg]\n";
    os << "    rRef [code]           = " << rRefCode_ << "\n";
    os << "    |B0| at rRef [code]   = " << bEqCode_ << "\n";
  } else {
    os << "    dipole                = off\n";
  }
  if (nMax_ > 0) {
    os << "    crustal file          = " << fileName_ << "\n";
    os << "    crustal degrees       = " << nMax_ << " (n = 0.." << nMax_ - 1
       << ")\n";
  } else {
    os << "    crustal field         = off\n";
  }
  os << "    rPlanet [code]        = " << rPlanetCode_;
  return os.str();
}




