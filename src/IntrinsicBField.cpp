#include <algorithm>
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
#include "IntrinsicBField.h"

using namespace amrex;

namespace {
constexpr Real cDegToRad = dPI / 180.0;
// The recurrence divides by sin(theta); below this the azimuthal component is
// set to zero instead, exactly as the reference implementation does.
constexpr Real cThetaSinMin = 1.0e-6;
// nanotesla -> tesla, the unit of #DIPOLE / #CRUSTALFIELD input.
constexpr Real cNTToT = 1.0e-9;
// mu0 / (4 * pi) in SI [T * m / A].
constexpr Real cMu0Over4PiSI = 1.0e-7;

std::string normalize_group(const std::string &group) {
  return group.empty() ? std::string("default") : group;
}

const char *waveform_kind_name(IntrinsicBField::WaveformKind kind) {
  switch (kind) {
    case IntrinsicBField::WaveformKind::Constant:
      return "constant";
    case IntrinsicBField::WaveformKind::LinearRamp:
      return "linear_ramp";
    case IntrinsicBField::WaveformKind::CosineRamp:
      return "cosine_ramp";
    case IntrinsicBField::WaveformKind::Sinusoid:
      return "sinusoid";
    case IntrinsicBField::WaveformKind::CrowbarRLC:
      return "crowbar_rlc";
    case IntrinsicBField::WaveformKind::Table:
      return "table";
  }
  return "unknown";
}
} // namespace

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
void IntrinsicBField::add_coil(const CoilSpecSI &spec) {
  coilSpecs_.push_back(spec);
  coilSpecs_.back().group = normalize_group(spec.group);
}

//==========================================================
void IntrinsicBField::add_toroidal_field(const ToroidalSpecSI &spec) {
  toroidalSpecs_.push_back(spec);
  toroidalSpecs_.back().group = normalize_group(spec.group);
}

//==========================================================
void IntrinsicBField::set_waveform(const std::string &group,
                                   const WaveformSpecSI &spec) {
  waveformSpecs_[normalize_group(group)] = spec;
}

//==========================================================
// Read the spherical harmonic coefficients from the crustal field file,
// matching the latest layout in BATSRUS ModUserMars:
// - 3 header lines (skipped)
// - For each degree i = 1..n:
//     one 'n m value' triple for g(i,0) followed by a g/h pair for each m
//     = 1..i.
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
        const int absM_c = std::abs(mm);
        if (nn >= 1 && nn < n && absM_c >= 0 && absM_c < n)
          c[nn * n + absM_c] = val;
        if (!(inFile >> nn >> mm >> val))
          break;
        const int absM_d = std::abs(mm);
        if (nn >= 1 && nn < n && absM_d >= 0 && absM_d < n)
          d[nn * n + absM_d] = val;
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
void IntrinsicBField::read_waveform_table(WaveformSpecSI &spec) {
  std::vector<Real> times;
  std::vector<Real> scales;

  if (ParallelDescriptor::IOProcessor()) {
    std::ifstream inFile(spec.fileName);
    if (!inFile.is_open()) {
      Abort("IntrinsicBField: cannot open waveform table file '" +
            spec.fileName + "'.");
    }

    std::string line;
    while (std::getline(inFile, line)) {
      const auto hashPos = line.find('#');
      if (hashPos != std::string::npos)
        line.erase(hashPos);
      std::istringstream iss(line);
      Real tVal = 0.0;
      Real fVal = 0.0;
      if (!(iss >> tVal >> fVal))
        continue;
      if (!times.empty() && !(tVal > times.back())) {
        Abort("IntrinsicBField: waveform table '" + spec.fileName +
              "' timestamps must be strictly increasing.");
      }
      times.push_back(tVal);
      scales.push_back(fVal);
    }
    inFile.close();
  }

  int nRows = static_cast<int>(times.size());
  ParallelDescriptor::Bcast(&nRows, 1, ParallelDescriptor::IOProcessorNumber());
  if (nRows < 2) {
    Abort("IntrinsicBField: waveform table '" + spec.fileName +
          "' must contain at least 2 (t_SI, scale) entries.");
  }

  if (!ParallelDescriptor::IOProcessor()) {
    times.resize(nRows);
    scales.resize(nRows);
  }
  ParallelDescriptor::Bcast(times.data(), nRows,
                            ParallelDescriptor::IOProcessorNumber());
  ParallelDescriptor::Bcast(scales.data(), nRows,
                            ParallelDescriptor::IOProcessorNumber());

  spec.tableTimeSI = std::move(times);
  spec.tableScale = std::move(scales);
}

//==========================================================
void IntrinsicBField::eval_waveform(const WaveformSpecSI &spec, Real tSI,
                                    Real &scale, Real &dScaleDtSI) {
  switch (spec.kind) {
    case WaveformKind::Constant:
      scale = spec.scale;
      dScaleDtSI = 0.0;
      return;

    case WaveformKind::LinearRamp: {
      const Real dt = spec.tEnd - spec.tStart;
      if (tSI < spec.tStart) {
        scale = spec.fStart;
        dScaleDtSI = 0.0;
      } else if (tSI >= spec.tEnd) {
        scale = spec.fEnd;
        dScaleDtSI = 0.0;
      } else {
        const Real slope = (spec.fEnd - spec.fStart) / dt;
        scale = spec.fStart + slope * (tSI - spec.tStart);
        dScaleDtSI = slope;
      }
      return;
    }

    case WaveformKind::CosineRamp: {
      const Real dt = spec.tEnd - spec.tStart;
      if (tSI < spec.tStart) {
        scale = spec.fStart;
        dScaleDtSI = 0.0;
      } else if (tSI >= spec.tEnd) {
        scale = spec.fEnd;
        dScaleDtSI = 0.0;
      } else {
        const Real u = (tSI - spec.tStart) / dt;
        const Real df = spec.fEnd - spec.fStart;
        scale = spec.fStart + 0.5 * df * (1.0 - std::cos(dPI * u));
        dScaleDtSI = (0.5 * dPI / dt) * df * std::sin(dPI * u);
      }
      return;
    }

    case WaveformKind::Sinusoid: {
      const Real phase0 = spec.phaseDeg * cDegToRad;
      if (tSI < spec.tStart) {
        scale = spec.fOffset + spec.fAmp * std::sin(phase0);
        dScaleDtSI = 0.0;
      } else {
        const Real omega = 2.0 * dPI * spec.freqHz;
        const Real angle = omega * (tSI - spec.tStart) + phase0;
        scale = spec.fOffset + spec.fAmp * std::sin(angle);
        dScaleDtSI = omega * spec.fAmp * std::cos(angle);
      }
      return;
    }

    case WaveformKind::CrowbarRLC: {
      const Real tau = tSI - spec.tStart;
      if (tau <= 0.0) {
        scale = 0.0;
        dScaleDtSI = 0.0;
        return;
      }
      const Real omega = 0.5 * dPI / spec.tRise;
      const Real gamma = (spec.tauDecay > 0.0) ? (1.0 / spec.tauDecay) : 0.0;
      const Real tPeak =
          (gamma > 0.0) ? (std::atan(omega / gamma) / omega) : spec.tRise;
      const Real aNorm = std::exp(gamma * tPeak) / std::sin(omega * tPeak);

      if (!spec.doCrowbar || tau <= tPeak) {
        const Real env = aNorm * std::exp(-gamma * tau);
        const Real sinWt = std::sin(omega * tau);
        const Real cosWt = std::cos(omega * tau);
        scale = env * sinWt;
        dScaleDtSI = env * (omega * cosWt - gamma * sinWt);
      } else if (spec.tauCrowbar > 0.0) {
        const Real dtPost = tau - tPeak;
        scale = std::exp(-dtPost / spec.tauCrowbar);
        dScaleDtSI = -scale / spec.tauCrowbar;
      } else {
        scale = 1.0;
        dScaleDtSI = 0.0;
      }
      return;
    }

    case WaveformKind::Table: {
      const auto &ts = spec.tableTimeSI;
      const auto &fs = spec.tableScale;
      if (ts.size() < 2) {
        scale = ts.empty() ? 1.0 : fs.front();
        dScaleDtSI = 0.0;
        return;
      }
      if (tSI < ts.front()) {
        scale = fs.front();
        dScaleDtSI = 0.0;
        return;
      }
      if (tSI >= ts.back()) {
        scale = fs.back();
        dScaleDtSI = 0.0;
        return;
      }
      const auto it = std::upper_bound(ts.begin(), ts.end(), tSI);
      const std::size_t idx = static_cast<std::size_t>(it - ts.begin());
      const Real slope = (fs[idx] - fs[idx - 1]) / (ts[idx] - ts[idx - 1]);
      scale = fs[idx - 1] + slope * (tSI - ts[idx - 1]);
      dScaleDtSI = slope;
      return;
    }
  }
}

//==========================================================
void IntrinsicBField::convert_units(Real si2NoB, Real si2NoL, Real si2NoV,
                                    Real rPlanetSI, const Real *bodyCenterCode,
                                    Real bodyRadiusCode, bool useBody,
                                    int nDim) {
  nDim_ = nDim;
  si2NoT_ = (si2NoV > 0.0) ? (si2NoL / si2NoV) : 0.0;
  rPlanetCode_ = rPlanetSI * si2NoL;

  const bool needsPlanetRadius =
      nMax_ > 0 || (useDipole_ && rRefSI_ <= 0.0 && !useBody);
  if (needsPlanetRadius && !(rPlanetCode_ > 0.0)) {
    Abort("IntrinsicBField: the planetary radius is not positive; check "
          "#PLANETRADIUS and #NORMALIZATION.");
  }

  for (int i = 0; i < 3; ++i)
    params_.center[i] = useBody ? bodyCenterCode[i] : 0.0;

  // Dipole: the reference radius is #DIPOLE's own radius when given, else the
  // radius of the inner body, else the planetary radius.
  if (rRefSI_ > 0.0)
    rRefCode_ = rRefSI_ * si2NoL;
  else if (useBody)
    rRefCode_ = bodyRadiusCode;
  else
    rRefCode_ = rPlanetCode_;

  bEqCode_ = bEqNT_ * cNTToT * si2NoB;
  params_.useDipole = useDipole_;
  if (useDipole_) {
    // The axis is the spherical direction (theta, phi) measured from +z,
    // tipped towards -x at phi = 0 so that it agrees with the BATSRUS dipole
    // tilt (Dipole_D = Bdp * [-sin(theta), 0, cos(theta)]).
    const Real theta = thetaDeg_ * cDegToRad;
    const Real phi = phiDeg_ * cDegToRad;
    const Real sinTheta = std::sin(theta);
    const Real moment = bEqCode_ * rRefCode_ * rRefCode_ * rRefCode_;

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

  // Register unique coil/toroidal groups and resolve attached waveforms.
  std::map<std::string, int> groupToIdx;
  groupNames_.clear();
  groupWaveforms_.clear();

  auto resolve_group = [&](const std::string &rawGroup) -> int {
    const std::string g = normalize_group(rawGroup);
    const auto it = groupToIdx.find(g);
    if (it != groupToIdx.end())
      return it->second;

    const int idx = static_cast<int>(groupNames_.size());
    groupToIdx[g] = idx;
    groupNames_.push_back(g);

    const auto wit = waveformSpecs_.find(g);
    if (wit != waveformSpecs_.end()) {
      WaveformSpecSI wf = wit->second;
      if (wf.kind == WaveformKind::Table) {
        read_waveform_table(wf);
      }
      groupWaveforms_.push_back(std::move(wf));
    } else {
      groupWaveforms_.push_back(WaveformSpecSI{});
    }
    return idx;
  };

  std::vector<CoilLoopParam> hostLoops;
  std::vector<WireSegmentParam> hostSegments;
  std::vector<Filament2DParam> hostFilaments;
  std::vector<ToroidalFieldParam> hostToroidals;

  const Real unitKI = cMu0Over4PiSI * si2NoB * si2NoL;

  for (const auto &spec : coilSpecs_) {
    const int gIdx = resolve_group(spec.group);

    if (spec.kind == CoilKind::Loop || spec.kind == CoilKind::Solenoid) {
      const Real nMag =
          std::sqrt(spec.p2[0] * spec.p2[0] + spec.p2[1] * spec.p2[1] +
                    spec.p2[2] * spec.p2[2]);
      if (!(nMag > 0.0)) {
        Abort("IntrinsicBField: #COIL axis vector must be nonzero.");
      }
      if (!(spec.radius > 0.0)) {
        Abort("IntrinsicBField: #COIL loop/solenoid radius must be positive.");
      }
      const Real unitAxis[3] = { spec.p2[0] / nMag, spec.p2[1] / nMag,
                                 spec.p2[2] / nMag };
      const int nTurns =
          (spec.kind == CoilKind::Solenoid) ? std::max(1, spec.nTurns) : 1;
      const Real lengthCode =
          (spec.kind == CoilKind::Solenoid) ? (spec.length * si2NoL) : 0.0;

      for (int iTurn = 0; iTurn < nTurns; ++iTurn) {
        const Real sOffset =
            (nTurns == 1)
                ? 0.0
                : (-0.5 * lengthCode +
                   (static_cast<Real>(iTurn) + 0.5) * (lengthCode / nTurns));
        CoilLoopParam p;
        for (int d = 0; d < 3; ++d) {
          p.center[d] = spec.p1[d] * si2NoL + sOffset * unitAxis[d];
          p.axis[d] = unitAxis[d];
        }
        p.radius = spec.radius * si2NoL;
        p.coreRadius = spec.coreRadius * si2NoL;
        p.kI0 = unitKI * spec.current;
        p.groupIdx = gIdx;
        hostLoops.push_back(p);
      }
    } else if (spec.kind == CoilKind::Segment) {
      const Real dx[3] = { (spec.p2[0] - spec.p1[0]) * si2NoL,
                           (spec.p2[1] - spec.p1[1]) * si2NoL,
                           (spec.p2[2] - spec.p1[2]) * si2NoL };
      const Real segLen =
          std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
      if (!(segLen > 0.0)) {
        Abort("IntrinsicBField: #COIL segment endpoints must not coincide.");
      }
      WireSegmentParam p;
      for (int d = 0; d < 3; ++d) {
        p.x1[d] = spec.p1[d] * si2NoL;
        p.dir[d] = dx[d] / segLen;
      }
      p.length = segLen;
      p.coreRadius = spec.coreRadius * si2NoL;
      p.kI0 = unitKI * spec.current;
      p.groupIdx = gIdx;
      hostSegments.push_back(p);
    } else if (spec.kind == CoilKind::Filament2D) {
      Filament2DParam p;
      const Real rRefCode = (si2NoL > 0.0) ? si2NoL : 1.0;
      p.x0 = spec.p1[0] * si2NoL;
      p.y0 = spec.p1[1] * si2NoL;
      p.coreRadius = spec.coreRadius * si2NoL;
      p.rRef2 = rRefCode * rRefCode;
      p.kI0 = unitKI * spec.current;
      p.groupIdx = gIdx;
      hostFilaments.push_back(p);
    }
  }

  for (const auto &spec : toroidalSpecs_) {
    const int gIdx = resolve_group(spec.group);
    const Real nMag =
        std::sqrt(spec.axis[0] * spec.axis[0] + spec.axis[1] * spec.axis[1] +
                  spec.axis[2] * spec.axis[2]);
    if (!(nMag > 0.0)) {
      Abort("IntrinsicBField: #TOROIDALFIELD axis vector must be nonzero.");
    }
    if (!(spec.rRefM > 0.0)) {
      Abort("IntrinsicBField: #TOROIDALFIELD rRef must be positive.");
    }
    ToroidalFieldParam p;
    for (int d = 0; d < 3; ++d) {
      p.center[d] = spec.center[d] * si2NoL;
      p.axis[d] = spec.axis[d] / nMag;
    }
    p.rRef = spec.rRefM * si2NoL;
    p.rCore = spec.rCoreM * si2NoL;
    p.kTor0 = (spec.bRefT * si2NoB) * p.rRef;
    p.groupIdx = gIdx;
    hostToroidals.push_back(p);
  }

  auto assign_managed = [](auto &dst, const auto &src) {
    dst.resize(src.size());
    for (std::size_t i = 0; i < src.size(); ++i)
      dst[i] = src[i];
  };
  assign_managed(loops_, hostLoops);
  assign_managed(segments_, hostSegments);
  assign_managed(filaments_, hostFilaments);
  assign_managed(toroidals_, hostToroidals);

  const int nGroups = static_cast<int>(groupWaveforms_.size());
  groupScaleB_.resize(nGroups, 1.0);
  groupScaleE_.resize(nGroups, 0.0);

  params_.loops = loops_.empty() ? nullptr : loops_.data();
  params_.nLoops = static_cast<int>(loops_.size());
  params_.segments = segments_.empty() ? nullptr : segments_.data();
  params_.nSegments = static_cast<int>(segments_.size());
  params_.filaments = filaments_.empty() ? nullptr : filaments_.data();
  params_.nFilaments = static_cast<int>(filaments_.size());
  params_.toroidals = toroidals_.empty() ? nullptr : toroidals_.data();
  params_.nToroidals = static_cast<int>(toroidals_.size());
  params_.groupScaleB = groupScaleB_.empty() ? nullptr : groupScaleB_.data();
  params_.groupScaleE = groupScaleE_.empty() ? nullptr : groupScaleE_.data();
  params_.nGroups = nGroups;

  hasTimeDepWaveform_ = false;
  for (const auto &wf : groupWaveforms_) {
    if (wf.kind != WaveformKind::Constant) {
      hasTimeDepWaveform_ = true;
      break;
    }
  }
  if (hasTimeDepWaveform_ && !(si2NoT_ > 0.0)) {
    Abort("IntrinsicBField: time-dependent #COILWAVEFORM requires positive "
          "normalization time conversion (si2NoL / si2NoV).");
  }

  for (const auto &kv : waveformSpecs_) {
    if (groupToIdx.find(kv.first) == groupToIdx.end()) {
      Print() << "  Warning: #COILWAVEFORM for group '" << kv.first
              << "' has no matching #COIL or #TOROIDALFIELD elements.\n";
    }
  }

  hasEvaluatedTime_ = false;
  update_time(0.0);
}

//==========================================================
bool IntrinsicBField::is_at_time(Real tCode) const {
  if (!hasEvaluatedTime_)
    return false;
  return std::abs(tCode - lastEvalTimeCode_) <=
         Real(1.0e-14) * (Real(1.0) + std::abs(tCode));
}

//==========================================================
bool IntrinsicBField::update_time(Real tCode) {
  if (groupWaveforms_.empty()) {
    lastEvalTimeCode_ = tCode;
    hasEvaluatedTime_ = true;
    return false;
  }
  if (is_at_time(tCode))
    return false;

  const Real tSI = (si2NoT_ > 0.0) ? (tCode / si2NoT_) : 0.0;
  for (std::size_t g = 0; g < groupWaveforms_.size(); ++g) {
    Real s = 1.0;
    Real dsDtSI = 0.0;
    eval_waveform(groupWaveforms_[g], tSI, s, dsDtSI);
    groupScaleB_[g] = s;
    groupScaleE_[g] = (si2NoT_ > 0.0) ? (-dsDtSI / si2NoT_) : 0.0;
  }

  lastEvalTimeCode_ = tCode;
  hasEvaluatedTime_ = true;
  return true;
}

//==========================================================
void IntrinsicBField::rescale_coefficients(Real si2NoB) {
  const int n = nMax_;
  const Real unit = cNTToT * si2NoB;

  Real power = rPlanetCode_ * rPlanetCode_; // n+2 with n = 0
  for (int iDeg = 0; iDeg < n; ++iDeg) {
    const Real scale = power * unit;
    for (int m = 0; m < n; ++m) {
      coefC_[iDeg * n + m] *= scale;
      coefD_[iDeg * n + m] *= scale;
    }
    power *= rPlanetCode_;
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

  // Spherical harmonics are only physically valid in the source-free region
  // outside the planet (r >= rPlanet). Inside the inner body (r < rPlanet),
  // (rPlanet / r)^(n+2) diverges rapidly (e.g. 10^40 for nMax=60 as r -> 0).
  // Clamp evaluation radius to rPlanet to prevent unphysical numerical blowup
  // inside the dead core.
  const Real rEff = (rPlanetCode_ > 0.0) ? std::max(r, rPlanetCode_) : r;

  // Colatitude and longitude of the point; theta is in [0, pi] so xtsin >= 0.
  const Real xtcos = zc / r;
  const Real xtsin = std::sqrt(std::max(Real(0.0), Real(1.0) - xtcos * xtcos));
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

  const Real invR = 1.0 / rEff;
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

  Real bsph[3] = { 0.0, 0.0, 0.0 };
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
        bsph[2] -= aorn[nn + 2] * RN(nn, m) * m / xtsin *
                   (-cc * xpsin[m] + dd * xpcos[m]);
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
  eval_intrinsic_b_dev(params_, x, y, z, nDim_, b);

  if (nMax_ <= 0)
    return;

  eval_crustal(x, y, z, b);
}

//==========================================================
std::string IntrinsicBField::describe() const {
  std::ostringstream os;
  os << "  intrinsic / external magnetic field B0 ("
     << (hasTimeDepWaveform_ ? "time-dependent, E0 = -dA0/dt active" : "static")
     << ", B = B1 + B0):\n";
  if (useDipole_ || nMax_ > 0) {
    os << "    center [code]         = (" << params_.center[0] << ", "
       << params_.center[1] << ", " << params_.center[2] << ")\n";
  }
  if (useDipole_) {
    os << "    dipole strength       = " << bEqNT_ << " [nT] at rRef\n";
    os << "    dipole tilt           = theta " << thetaDeg_ << ", phi "
       << phiDeg_ << " [deg]\n";
    os << "    rRef [code]           = " << rRefCode_ << "\n";
    os << "    |B0| at rRef [code]   = " << bEqCode_ << "\n";
    os << "    rPlanet [code]        = " << rPlanetCode_ << "\n";
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
  if (params_.nLoops > 0 || params_.nSegments > 0 || params_.nFilaments > 0 ||
      params_.nToroidals > 0) {
    os << "    external elements     = " << params_.nLoops << " loop(s), "
       << params_.nSegments << " segment(s), " << params_.nFilaments
       << " filament2d(s), " << params_.nToroidals << " toroidal(s)\n";
    for (std::size_t g = 0; g < groupNames_.size(); ++g) {
      os << "    group '" << groupNames_[g]
         << "' waveform  = " << waveform_kind_name(groupWaveforms_[g].kind);
      if (g + 1 < groupNames_.size())
        os << "\n";
    }
  } else {
    os << "    external coils        = off";
  }
  return os.str();
}
