#include <algorithm>
#include <cmath>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include "Constants.h"
#include "ExternalCoilField.h"

using namespace amrex;

namespace {
constexpr Real cDegToRad = dPI / 180.0;
// mu0 / (4 * pi) in SI [T * m / A].
constexpr Real cMu0Over4PiSI = 1.0e-7;

std::string normalize_group(const std::string &group) {
  return group.empty() ? std::string("default") : group;
}

const char *waveform_kind_name(ExternalCoilField::WaveformKind kind) {
  switch (kind) {
    case ExternalCoilField::WaveformKind::Constant:
      return "constant";
    case ExternalCoilField::WaveformKind::LinearRamp:
      return "linear_ramp";
    case ExternalCoilField::WaveformKind::CosineRamp:
      return "cosine_ramp";
    case ExternalCoilField::WaveformKind::Sinusoid:
      return "sinusoid";
    case ExternalCoilField::WaveformKind::CrowbarRLC:
      return "crowbar_rlc";
    case ExternalCoilField::WaveformKind::Table:
      return "table";
  }
  return "unknown";
}
} // namespace

//==========================================================
void ExternalCoilField::read_coil_param(ReadParam &param) {
  CoilSpecSI spec;
  std::string type;
  param.read_var("type", type);
  param.read_var("group", spec.group);

  if (type == "loop" || type == "solenoid") {
    spec.kind = (type == "loop") ? CoilKind::Loop : CoilKind::Solenoid;
    param.read_var("centerX", spec.p1[0]);
    param.read_var("centerY", spec.p1[1]);
    param.read_var("centerZ", spec.p1[2]);
    param.read_var("axisX", spec.p2[0]);
    param.read_var("axisY", spec.p2[1]);
    param.read_var("axisZ", spec.p2[2]);
    param.read_var("radius", spec.radius);
    if (type == "solenoid") {
      param.read_var("length", spec.length);
      param.read_var("nTurns", spec.nTurns);
      if (spec.length < 0.0 || spec.nTurns < 1) {
        Abort("Error: #COIL solenoid requires length >= 0 and nTurns >= 1.");
      }
    }
    param.read_var("current", spec.current);
    param.read_var("coreRadius", spec.coreRadius);

    const Real axisMag2 = spec.p2[0] * spec.p2[0] + spec.p2[1] * spec.p2[1] +
                          spec.p2[2] * spec.p2[2];
    if (!(axisMag2 > 0.0)) {
      Abort("Error: #COIL loop/solenoid axis vector must be nonzero.");
    }
    if (!(spec.radius > 0.0)) {
      Abort("Error: #COIL loop/solenoid radius must be positive.");
    }
  } else if (type == "segment") {
    spec.kind = CoilKind::Segment;
    param.read_var("x1", spec.p1[0]);
    param.read_var("y1", spec.p1[1]);
    param.read_var("z1", spec.p1[2]);
    param.read_var("x2", spec.p2[0]);
    param.read_var("y2", spec.p2[1]);
    param.read_var("z2", spec.p2[2]);
    param.read_var("current", spec.current);
    param.read_var("coreRadius", spec.coreRadius);

    const Real dx = spec.p2[0] - spec.p1[0];
    const Real dy = spec.p2[1] - spec.p1[1];
    const Real dz = spec.p2[2] - spec.p1[2];
    if (!(dx * dx + dy * dy + dz * dz > 0.0)) {
      Abort("Error: #COIL segment endpoints must not coincide.");
    }
  } else if (type == "filament2d") {
    spec.kind = CoilKind::Filament2D;
    param.read_var("x0", spec.p1[0]);
    param.read_var("y0", spec.p1[1]);
    spec.p1[2] = 0.0;
    param.read_var("current", spec.current);
    param.read_var("coreRadius", spec.coreRadius);
  } else {
    Abort("Error: unknown #COIL type '" + type +
          "'. Expected loop, solenoid, segment, or filament2d.");
  }

  if (spec.coreRadius < 0.0) {
    Abort("Error: #COIL coreRadius must be non-negative.");
  }

  add_coil(spec);

  Print() << "  intrinsic field: coil type = " << type
          << ", group = " << spec.group << ", I = " << spec.current << " [A]\n";
}

//==========================================================
void ExternalCoilField::read_toroidal_param(ReadParam &param) {
  ToroidalSpecSI spec;
  param.read_var("group", spec.group);
  param.read_var("centerX", spec.center[0]);
  param.read_var("centerY", spec.center[1]);
  param.read_var("centerZ", spec.center[2]);
  param.read_var("axisX", spec.axis[0]);
  param.read_var("axisY", spec.axis[1]);
  param.read_var("axisZ", spec.axis[2]);
  param.read_var("bRef", spec.bRefT);
  param.read_var("rRef", spec.rRefM);
  param.read_var("rCore", spec.rCoreM);

  const Real axisMag2 = spec.axis[0] * spec.axis[0] +
                        spec.axis[1] * spec.axis[1] +
                        spec.axis[2] * spec.axis[2];
  if (!(axisMag2 > 0.0)) {
    Abort("Error: #TOROIDALFIELD axis vector must be nonzero.");
  }
  if (!(spec.rRefM > 0.0)) {
    Abort("Error: #TOROIDALFIELD rRef must be positive.");
  }
  if (spec.rCoreM < 0.0) {
    Abort("Error: #TOROIDALFIELD rCore must be non-negative.");
  }

  add_toroidal_field(spec);

  Print() << "  intrinsic field: toroidal group = " << spec.group
          << ", bRef = " << spec.bRefT << " [T] at rRef = " << spec.rRefM
          << " [m]\n";
}

//==========================================================
void ExternalCoilField::read_waveform_param(ReadParam &param) {
  std::string group;
  std::string type;
  param.read_var("group", group);
  param.read_var("type", type);

  WaveformSpecSI spec;
  if (type == "constant") {
    spec.kind = WaveformKind::Constant;
    param.read_var("scale", spec.scale);
  } else if (type == "linear_ramp" || type == "cosine_ramp") {
    spec.kind = (type == "linear_ramp") ? WaveformKind::LinearRamp
                                        : WaveformKind::CosineRamp;
    param.read_var("tStart", spec.tStart);
    param.read_var("tEnd", spec.tEnd);
    param.read_var("fStart", spec.fStart);
    param.read_var("fEnd", spec.fEnd);
    if (!(spec.tEnd > spec.tStart)) {
      Abort("Error: #COILWAVEFORM " + type + " requires tEnd > tStart.");
    }
  } else if (type == "sinusoid") {
    spec.kind = WaveformKind::Sinusoid;
    param.read_var("tStart", spec.tStart);
    param.read_var("fOffset", spec.fOffset);
    param.read_var("fAmp", spec.fAmp);
    param.read_var("freqHz", spec.freqHz);
    param.read_var("phaseDeg", spec.phaseDeg);
    if (spec.freqHz < 0.0) {
      Abort("Error: #COILWAVEFORM sinusoid requires freqHz >= 0.");
    }
  } else if (type == "crowbar_rlc") {
    spec.kind = WaveformKind::CrowbarRLC;
    param.read_var("tStart", spec.tStart);
    param.read_var("tRise", spec.tRise);
    param.read_var("tauDecay", spec.tauDecay);
    param.read_var("doCrowbar", spec.doCrowbar);
    param.read_var("tauCrowbar", spec.tauCrowbar);
    if (!(spec.tRise > 0.0)) {
      Abort("Error: #COILWAVEFORM crowbar_rlc requires tRise > 0.");
    }
    if (spec.tauDecay < 0.0 || spec.tauCrowbar < 0.0) {
      Abort("Error: #COILWAVEFORM crowbar_rlc requires non-negative decay "
            "times.");
    }
  } else if (type == "table") {
    spec.kind = WaveformKind::Table;
    param.read_var("fileName", spec.fileName);
  } else {
    Abort("Error: unknown #COILWAVEFORM type '" + type +
          "'. Expected constant, linear_ramp, cosine_ramp, sinusoid, "
          "crowbar_rlc, or table.");
  }

  set_waveform(group, spec);

  Print() << "  intrinsic field: waveform group = " << group
          << ", type = " << type << "\n";
}

//==========================================================
void ExternalCoilField::add_coil(const CoilSpecSI &spec) {
  coilSpecs_.push_back(spec);
  coilSpecs_.back().group = normalize_group(spec.group);
}

//==========================================================
void ExternalCoilField::add_toroidal_field(const ToroidalSpecSI &spec) {
  toroidalSpecs_.push_back(spec);
  toroidalSpecs_.back().group = normalize_group(spec.group);
}

//==========================================================
void ExternalCoilField::set_waveform(const std::string &group,
                                     const WaveformSpecSI &spec) {
  waveformSpecs_[normalize_group(group)] = spec;
}

//==========================================================
void ExternalCoilField::read_waveform_table(WaveformSpecSI &spec) {
  std::vector<Real> times;
  std::vector<Real> scales;

  if (ParallelDescriptor::IOProcessor()) {
    std::ifstream inFile(spec.fileName);
    if (!inFile.is_open()) {
      Abort("ExternalCoilField: cannot open waveform table file '" +
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
        Abort("ExternalCoilField: waveform table '" + spec.fileName +
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
    Abort("ExternalCoilField: waveform table '" + spec.fileName +
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
void ExternalCoilField::precompute_waveform(WaveformSpecSI &spec) {
  if (spec.kind == WaveformKind::Sinusoid) {
    spec.phase0 = spec.phaseDeg * cDegToRad;
    spec.omega = 2.0 * dPI * spec.freqHz;
  } else if (spec.kind == WaveformKind::CrowbarRLC) {
    spec.omega = 0.5 * dPI / spec.tRise;
    spec.gamma = (spec.tauDecay > 0.0) ? (1.0 / spec.tauDecay) : 0.0;
    spec.tPeak = (spec.gamma > 0.0)
                     ? (std::atan(spec.omega / spec.gamma) / spec.omega)
                     : spec.tRise;
    spec.aNorm =
        std::exp(spec.gamma * spec.tPeak) / std::sin(spec.omega * spec.tPeak);
    spec.invTauCrowbar =
        (spec.tauCrowbar > 0.0) ? (1.0 / spec.tauCrowbar) : 0.0;
  }
}

//==========================================================
void ExternalCoilField::eval_waveform(const WaveformSpecSI &spec, Real tSI,
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
      if (tSI < spec.tStart) {
        scale = spec.fOffset + spec.fAmp * std::sin(spec.phase0);
        dScaleDtSI = 0.0;
      } else {
        const Real angle = spec.omega * (tSI - spec.tStart) + spec.phase0;
        scale = spec.fOffset + spec.fAmp * std::sin(angle);
        dScaleDtSI = spec.omega * spec.fAmp * std::cos(angle);
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
      if (!spec.doCrowbar || tau <= spec.tPeak) {
        const Real env = spec.aNorm * std::exp(-spec.gamma * tau);
        const Real sinWt = std::sin(spec.omega * tau);
        const Real cosWt = std::cos(spec.omega * tau);
        scale = env * sinWt;
        dScaleDtSI = env * (spec.omega * cosWt - spec.gamma * sinWt);
      } else if (spec.invTauCrowbar > 0.0) {
        const Real dtPost = tau - spec.tPeak;
        scale = std::exp(-dtPost * spec.invTauCrowbar);
        dScaleDtSI = -scale * spec.invTauCrowbar;
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
void ExternalCoilField::convert_units(Real si2NoB, Real si2NoL, Real si2NoV) {
  si2NoT_ = (si2NoV > 0.0) ? (si2NoL / si2NoV) : 0.0;

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
      precompute_waveform(wf);
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
        Abort("ExternalCoilField: #COIL axis vector must be nonzero.");
      }
      if (!(spec.radius > 0.0)) {
        Abort(
            "ExternalCoilField: #COIL loop/solenoid radius must be positive.");
      }
      const Real unitAxis[3] = { spec.p2[0] / nMag, spec.p2[1] / nMag,
                                 spec.p2[2] / nMag };
      const int nTurns =
          (spec.kind == CoilKind::Solenoid) ? std::max(1, spec.nTurns) : 1;
      const Real lengthCode =
          (spec.kind == CoilKind::Solenoid) ? (spec.length * si2NoL) : 0.0;
      const Real radiusCode = spec.radius * si2NoL;
      const Real r2 = radiusCode * radiusCode;
      const Real coreCode = spec.coreRadius * si2NoL;
      const Real coreEps2 =
          (coreCode > 0.0) ? (coreCode * coreCode) : (Real(1.0e-24) * r2);

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
        p.radius = radiusCode;
        p.r2 = r2;
        p.fourR = Real(4.0) * radiusCode;
        p.coreEps2 = coreEps2;
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
        Abort("ExternalCoilField: #COIL segment endpoints must not coincide.");
      }
      const Real coreCode = spec.coreRadius * si2NoL;
      WireSegmentParam p;
      for (int d = 0; d < 3; ++d) {
        p.x1[d] = spec.p1[d] * si2NoL;
        p.dir[d] = dx[d] / segLen;
      }
      p.length = segLen;
      p.coreEps2 = coreCode * coreCode + Real(1.0e-24) * segLen * segLen;
      p.kI0 = unitKI * spec.current;
      p.groupIdx = gIdx;
      hostSegments.push_back(p);
    } else if (spec.kind == CoilKind::Filament2D) {
      Filament2DParam p;
      const Real rRefCode = (si2NoL > 0.0) ? si2NoL : 1.0;
      const Real rRef2 = rRefCode * rRefCode;
      const Real coreCode = spec.coreRadius * si2NoL;
      p.x0 = spec.p1[0] * si2NoL;
      p.y0 = spec.p1[1] * si2NoL;
      p.coreEps2 = coreCode * coreCode + Real(1.0e-24) * rRef2;
      p.invRRef2 = Real(1.0) / rRef2;
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
      Abort("ExternalCoilField: #TOROIDALFIELD axis vector must be nonzero.");
    }
    if (!(spec.rRefM > 0.0)) {
      Abort("ExternalCoilField: #TOROIDALFIELD rRef must be positive.");
    }
    ToroidalFieldParam p;
    for (int d = 0; d < 3; ++d) {
      p.center[d] = spec.center[d] * si2NoL;
      p.axis[d] = spec.axis[d] / nMag;
    }
    const Real rRefCode = spec.rRefM * si2NoL;
    const Real rRef2 = rRefCode * rRefCode;
    const Real rCoreCode = spec.rCoreM * si2NoL;
    p.coreEps2 = rCoreCode * rCoreCode + Real(1.0e-24) * rRef2;
    p.invRRef2 = Real(1.0) / rRef2;
    p.kTor0 = (spec.bRefT * si2NoB) * rRefCode;
    p.groupIdx = gIdx;
    hostToroidals.push_back(p);
  }

  hasTimeDepWaveform_ = false;
  for (const auto &wf : groupWaveforms_) {
    if (wf.kind != WaveformKind::Constant) {
      hasTimeDepWaveform_ = true;
      break;
    }
  }
  if (hasTimeDepWaveform_ && !(si2NoT_ > 0.0)) {
    Abort("ExternalCoilField: time-dependent #COILWAVEFORM requires positive "
          "normalization time conversion (si2NoL / si2NoV).");
  }

  // Fold static Constant waveform scale into kI0/kTor0 so static runs skip
  // per-group scale lookups entirely.
  if (!hasTimeDepWaveform_) {
    for (auto &p : hostLoops)
      p.kI0 *= groupWaveforms_[p.groupIdx].scale;
    for (auto &p : hostSegments)
      p.kI0 *= groupWaveforms_[p.groupIdx].scale;
    for (auto &p : hostFilaments)
      p.kI0 *= groupWaveforms_[p.groupIdx].scale;
    for (auto &p : hostToroidals)
      p.kTor0 *= groupWaveforms_[p.groupIdx].scale;
  }

  auto copy_to_device = [](auto &dst, const auto &src) {
    dst.resize(src.size());
    if (!src.empty()) {
      Gpu::copy(Gpu::hostToDevice, src.begin(), src.end(), dst.begin());
    }
  };
  copy_to_device(loops_, hostLoops);
  copy_to_device(segments_, hostSegments);
  copy_to_device(filaments_, hostFilaments);
  copy_to_device(toroidals_, hostToroidals);

  const int nGroups = static_cast<int>(groupWaveforms_.size());
  params_.loops = loops_.empty() ? nullptr : loops_.data();
  params_.nLoops = static_cast<int>(loops_.size());
  params_.segments = segments_.empty() ? nullptr : segments_.data();
  params_.nSegments = static_cast<int>(segments_.size());
  params_.filaments = filaments_.empty() ? nullptr : filaments_.data();
  params_.nFilaments = static_cast<int>(filaments_.size());
  params_.toroidals = toroidals_.empty() ? nullptr : toroidals_.data();
  params_.nToroidals = static_cast<int>(toroidals_.size());
  params_.nGroups = nGroups;
  params_.hasTimeDep = hasTimeDepWaveform_;
  params_.inlineScaleB.fill(Real(1.0));
  params_.inlineScaleE.fill(Real(0.0));
  params_.extScaleB = nullptr;
  params_.extScaleE = nullptr;

  if (!hasTimeDepWaveform_) {
    curScaleB_.clear();
    curScaleE_.clear();
    devScaleB_.clear();
    devScaleE_.clear();
    anyNonzeroE_ = false;
    hasEvaluatedTime_ = true;
    lastEvalTimeCode_ = 0.0;
  } else {
    curScaleB_.assign(nGroups, Real(1.0));
    curScaleE_.assign(nGroups, Real(0.0));
    if (nGroups > CoilFieldDeviceParams::kMaxInlineGroups) {
      devScaleB_.resize(nGroups);
      devScaleE_.resize(nGroups);
      params_.extScaleB = devScaleB_.data();
      params_.extScaleE = devScaleE_.data();
    } else {
      devScaleB_.clear();
      devScaleE_.clear();
    }
    anyNonzeroE_ = false;
    hasEvaluatedTime_ = false;
    update_time(0.0);
  }

  for (const auto &kv : waveformSpecs_) {
    if (groupToIdx.find(kv.first) == groupToIdx.end()) {
      Print() << "  Warning: #COILWAVEFORM for group '" << kv.first
              << "' has no matching #COIL or #TOROIDALFIELD elements.\n";
    }
  }
}

//==========================================================
bool ExternalCoilField::is_at_time(Real tCode) const {
  if (!hasEvaluatedTime_)
    return false;
  return std::abs(tCode - lastEvalTimeCode_) <=
         Real(1.0e-14) * (Real(1.0) + std::abs(tCode));
}

//==========================================================
ExternalCoilField::TimeUpdateStatus ExternalCoilField::update_time(Real tCode) {
  if (!hasTimeDepWaveform_ || groupWaveforms_.empty()) {
    lastEvalTimeCode_ = tCode;
    hasEvaluatedTime_ = true;
    return {};
  }
  if (hasEvaluatedTime_ && is_at_time(tCode)) {
    return {};
  }

  const bool firstEval = !hasEvaluatedTime_;
  bool changedB = firstEval;
  bool changedE = firstEval;
  bool anyNonzeroE = false;

  const Real tSI = (si2NoT_ > 0.0) ? (tCode / si2NoT_) : 0.0;
  const int nGroups = static_cast<int>(groupWaveforms_.size());
  for (int g = 0; g < nGroups; ++g) {
    Real s = 1.0;
    Real dsDtSI = 0.0;
    eval_waveform(groupWaveforms_[g], tSI, s, dsDtSI);
    const Real sE = (si2NoT_ > 0.0) ? (-dsDtSI / si2NoT_) : 0.0;

    if (!firstEval) {
      if (s != curScaleB_[g])
        changedB = true;
      if (sE != curScaleE_[g])
        changedE = true;
    }
    curScaleB_[g] = s;
    curScaleE_[g] = sE;
    if (sE != Real(0.0))
      anyNonzeroE = true;

    if (nGroups <= CoilFieldDeviceParams::kMaxInlineGroups) {
      params_.inlineScaleB[g] = s;
      params_.inlineScaleE[g] = sE;
    }
  }

  if (nGroups > CoilFieldDeviceParams::kMaxInlineGroups) {
    if (changedB) {
      Gpu::copy(Gpu::hostToDevice, curScaleB_.begin(), curScaleB_.end(),
                devScaleB_.begin());
    }
    if (changedE) {
      Gpu::copy(Gpu::hostToDevice, curScaleE_.begin(), curScaleE_.end(),
                devScaleE_.begin());
    }
  }

  anyNonzeroE_ = anyNonzeroE;
  lastEvalTimeCode_ = tCode;
  hasEvaluatedTime_ = true;
  return { changedB, changedE };
}

//==========================================================
void ExternalCoilField::describe(std::ostream &os) const {
  if (params_.any_active()) {
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
}
