#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include "Constants.h"
#include "FluidInterface.h"
#include "WaveBC.h"

namespace {
inline constexpr amrex::Real cTwoPi = 2.0 * dPI;

const MonoWave sMonoWave;
const WavePacket sWavePacket;
const SpectralWave sSpectralWave;

const WaveProfile* get_profile_evaluator(int profile) {
  switch (profile) {
    case WaveComponent::kPacket:
      return &sWavePacket;
    case WaveComponent::kSpectral:
      return &sSpectralWave;
    case WaveComponent::kMono:
    default:
      return &sMonoWave;
  }
}

inline amrex::Real si2no(const FluidInterface* fi, int iField) {
  if (!fi)
    return 1.0;
  switch (iField) {
    case 0: // B (T)
      return fi->get_No2SiB() > 0.0 ? (1.0 / fi->get_No2SiB()) : 1.0;
    case 1: // E (V/m)
      return fi->get_No2SiE() > 0.0 ? (1.0 / fi->get_No2SiE()) : 1.0;
    case 2: // velocity (m/s)
      return fi->get_Si2NoV() > 0.0 ? fi->get_Si2NoV() : 1.0;
    case 3: // density (kg/m^3)
      return fi->get_Si2NoRho() > 0.0 ? fi->get_Si2NoRho() : 1.0;
    default:
      return 1.0;
  }
}

inline void cross3(const amrex::Real a[3], const amrex::Real b[3],
                   amrex::Real out[3]) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

inline amrex::Real dot3(const amrex::Real a[3], const amrex::Real b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline amrex::Real norm3(const amrex::Real a[3]) {
  return std::sqrt(dot3(a, a));
}
} // namespace

amrex::Real MonoWave::value(const WaveComponent& c, amrex::Real t,
                            const amrex::Real* pos) const {
  if (c.amplitude == 0.0)
    return 0.0;
  const amrex::Real kdotx =
      c.k_vec[0] * pos[0] + c.k_vec[1] * pos[1] + c.k_vec[2] * pos[2];
  const amrex::Real arg = kdotx - c.frequency * t + c.phase;
  return c.amplitude * std::sin(arg);
}

amrex::Real WavePacket::value(const WaveComponent& c, amrex::Real t,
                              const amrex::Real* pos) const {
  const amrex::Real carrier = MonoWave::value(c, t, pos);
  if (carrier == 0.0)
    return 0.0;
  const amrex::Real tw = (c.tWidth > 0.0) ? c.tWidth : 1.0;
  const amrex::Real tau = (t - c.tCenter) / tw;
  const amrex::Real g = std::exp(-tau * tau);
  return carrier * g;
}

amrex::Real SpectralWave::value(const WaveComponent& c, amrex::Real t,
                                const amrex::Real* pos) const {
  if (c.comps.empty())
    return 0.0;
  amrex::Real sum = 0.0;
  for (const auto& sub : c.comps)
    sum += sMonoWave.value(sub, t, pos);
  return sum;
}

std::unique_ptr<WaveProfile> WaveRegistry::create(const WaveComponent& c) {
  switch (c.profile) {
    case WaveComponent::kPacket:
      return std::make_unique<WavePacket>();
    case WaveComponent::kSpectral:
      return std::make_unique<SpectralWave>();
    case WaveComponent::kCustom:
      break; // handled via c.custom; no dedicated evaluator
    case WaveComponent::kMono:
    default:
      break;
  }
  return std::make_unique<MonoWave>();
}

amrex::Real WaveWindow::weight(const amrex::Real* pos) const {
  if (!active)
    return 1.0;
  const amrex::Real x = pos[direction];
  if (x < lo || x > hi)
    return 0.0;
  if (taperLo > 0.0 && (x - lo) < taperLo)
    return 0.5 * (1.0 - std::cos(dPI * (x - lo) / taperLo));
  if (taperHi > 0.0 && (hi - x) < taperHi)
    return 0.5 * (1.0 - std::cos(dPI * (hi - x) / taperHi));
  return 1.0;
}

void WaveWindow::weight_and_grad(const amrex::Real* pos, amrex::Real& w,
                                 amrex::Real gradW[3]) const {
  gradW[0] = 0.0;
  gradW[1] = 0.0;
  gradW[2] = 0.0;
  if (!active) {
    w = 1.0;
    return;
  }
  const amrex::Real x = pos[direction];
  if (x < lo || x > hi) {
    w = 0.0;
    return;
  }
  if (taperLo > 0.0 && (x - lo) < taperLo) {
    const amrex::Real arg = dPI * (x - lo) / taperLo;
    w = 0.5 * (1.0 - std::cos(arg));
    gradW[direction] = 0.5 * (dPI / taperLo) * std::sin(arg);
    return;
  }
  if (taperHi > 0.0 && (hi - x) < taperHi) {
    const amrex::Real arg = dPI * (hi - x) / taperHi;
    w = 0.5 * (1.0 - std::cos(arg));
    gradW[direction] = -0.5 * (dPI / taperHi) * std::sin(arg);
    return;
  }
  w = 1.0;
}

void WaveBoundaryManager::read_window(ReadParam& param) {
  param.read_var("direction", window.direction);
  if (window.direction < 0 || window.direction >= 3)
    amrex::Abort("WAVEBCWINDOW: direction must be 0 (x), 1 (y), or 2 (z).");
  param.read_var("min", window.lo);
  param.read_var("max", window.hi);
  param.read_var("taperWidth", window.taper);
  if (window.hi <= window.lo)
    amrex::Abort("WAVEBCWINDOW: max must be larger than min.");
  if (window.taper < 0.0)
    amrex::Abort("WAVEBCWINDOW: taperWidth must be non-negative.");
  window.taperLo = window.taper;
  window.taperHi = window.taper;
  window.active = true;
  window.isCodeUnits = false;
}

void WaveBoundaryManager::read_turbulent_inflow(ReadParam& param) {
  active = true;
  turb.active = true;
  turb.initialized = false;

  param.read_var("nMode", turb.nMode);
  param.read_var("lambdaMin", turb.lambdaMin);
  param.read_var("lambdaMax", turb.lambdaMax);
  param.read_var("spectralIndex", turb.spectralIndex);
  param.read_var("dBrmsOverB0", turb.dBrmsOverB0);
  param.read_var("seed", turb.seed);
  param.read_var("rampTime", turb.rampTimeSI);

  if (turb.nMode <= 0)
    amrex::Abort("TURBULENTINFLOW: nMode must be > 0.");
  if (turb.lambdaMin <= 0.0)
    amrex::Abort("TURBULENTINFLOW: lambdaMin must be > 0.");
  if (turb.lambdaMax < turb.lambdaMin)
    amrex::Abort("TURBULENTINFLOW: lambdaMax must be >= lambdaMin.");
  if (turb.dBrmsOverB0 < 0.0)
    amrex::Abort("TURBULENTINFLOW: dBrmsOverB0 must be >= 0.");
  if (turb.rampTimeSI < 0.0)
    amrex::Abort("TURBULENTINFLOW: rampTime must be >= 0.");

  // Optional trailing parameters
  bool progress = true;
  while (progress) {
    progress = false;
    if (param.read_optional("polarization", turb.polarization) ||
        param.read_optional("family", turb.polarization)) {
      for (char& ch : turb.polarization)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
      if (turb.polarization != "out_of_plane" &&
          turb.polarization != "in_plane" && turb.polarization != "alfven" &&
          turb.polarization != "isotropic") {
        amrex::Abort("TURBULENTINFLOW: polarization must be 'out_of_plane', "
                     "'in_plane', 'alfven', or 'isotropic', got '" +
                     turb.polarization + "'.");
      }
      progress = true;
    }
    if (param.read_optional("sigmaC", turb.sigmaC)) {
      if (turb.sigmaC < -1.0 || turb.sigmaC > 1.0)
        amrex::Abort("TURBULENTINFLOW: sigmaC must be in [-1, 1].");
      progress = true;
    }
    if (param.read_optional("tStart", turb.tStartSI)) {
      progress = true;
    }
    if (param.read_optional("modeFile", turb.modeFile)) {
      progress = true;
    }
  }
}

void WaveBoundaryManager::clear() {
  faces.clear();
  turb = TurbulentInflowConfig{};
  active = false;
}

WaveFace& WaveBoundaryManager::face(int direction, int side) {
  for (auto& f : faces)
    if (f.direction == direction && f.side == side)
      return f;
  faces.push_back(WaveFace{});
  faces.back().direction = direction;
  faces.back().side = side;
  return faces.back();
}

void WaveBoundaryManager::add_component(int direction, int side,
                                        const WaveComponent& c) {
  WaveComponent comp = c; // already in code units
  comp.isCodeUnits = true;
  if (comp.waveLength > 0.0) {
    const amrex::Real k = cTwoPi / comp.waveLength;
    for (int d = 0; d < 3; ++d)
      comp.k_vec[d] = k * comp.dir[d];
  } else {
    for (int d = 0; d < 3; ++d)
      comp.k_vec[d] = 0.0;
  }
  face(direction, side).comps.push_back(comp);
}

// cos^2 ramp-in, and smooth ramp-out before tEnd if rampTime > 0.
amrex::Real WaveBoundaryManager::envelope(const WaveComponent& c,
                                          amrex::Real t) const {
  if (t < c.tStart)
    return 0.0;
  if (c.tEnd > 0.0 && t > c.tEnd)
    return 0.0;
  amrex::Real ramp = 1.0;
  if (c.rampTime > 0.0) {
    if (t < c.tStart + c.rampTime) {
      amrex::Real f = (t - c.tStart) / c.rampTime;
      ramp = 0.5 * (1.0 - std::cos(dPI * f));
    } else if (c.tEnd > 0.0 && t > c.tEnd - c.rampTime) {
      amrex::Real f = (c.tEnd - t) / c.rampTime;
      ramp = 0.5 * (1.0 - std::cos(dPI * f));
    }
  }
  return ramp;
}

amrex::Real WaveBoundaryManager::turbulent_envelope(amrex::Real t) const {
  if (!turb.active || turb.dBrmsOverB0 <= 0.0 || t < turb.tStart)
    return 0.0;
  if (turb.rampTime > 0.0 && t < turb.tStart + turb.rampTime) {
    const amrex::Real f = (t - turb.tStart) / turb.rampTime;
    return 0.5 * (1.0 - std::cos(dPI * f));
  }
  return 1.0;
}

amrex::Real WaveBoundaryManager::value(const WaveComponent& c, amrex::Real t,
                                       const amrex::Real* pos) const {
  if (!active || c.amplitude == 0.0)
    return 0.0;
  const amrex::Real env = envelope(c, t);
  if (env == 0.0)
    return 0.0;
  if (c.profile == WaveComponent::kCustom && c.custom)
    return c.custom(c, t, pos) * env;
  const WaveProfile* prof = get_profile_evaluator(c.profile);
  return prof->value(c, t, pos) * env;
}

void WaveBoundaryManager::read_param(ReadParam& param) {
  active = true;
  int nFace = 0;
  param.read_var("nWaveFace", nFace);

  for (int iFace = 0; iFace < nFace; ++iFace) {
    int direction = 0;
    std::string sideStr;
    param.read_var("direction", direction);
    if (direction < 0 || direction >= 3) {
      amrex::Abort("WAVEBC: direction must be 0 (x), 1 (y), or 2 (z).");
    }
    param.read_var("side", sideStr);
    int side = 0;
    if (sideStr == "lo") {
      side = 0;
    } else if (sideStr == "hi") {
      side = 1;
    } else {
      amrex::Abort("WAVEBC: side must be 'lo' or 'hi', got '" + sideStr + "'.");
    }

    int nComp = 0;
    param.read_var("nComp", nComp);

    WaveFace& f = face(direction, side);
    for (int ic = 0; ic < nComp; ++ic) {
      WaveComponent c;
      param.read_var("iField", c.iField);
      param.read_var("iSpecies", c.iSpecies);
      param.read_var("profile", c.profile);
      param.read_var("amplitude", c.amplitude);
      param.read_var("frequency", c.frequency);
      param.read_var("waveLength", c.waveLength);
      param.read_var("phase", c.phase);
      param.read_var("dir", c.dir[0]);
      param.read_var("dir", c.dir[1]);
      param.read_var("dir", c.dir[2]);
      param.read_var("pol", c.pol[0]);
      param.read_var("pol", c.pol[1]);
      param.read_var("pol", c.pol[2]);
      param.read_var("tStart", c.tStart);
      param.read_var("tEnd", c.tEnd);
      param.read_var("rampTime", c.rampTime);
      param.read_var("tCenter", c.tCenter);
      param.read_var("tWidth", c.tWidth);

      // Normalize propagation direction vector if non-zero.
      const amrex::Real dirNorm = std::sqrt(
          c.dir[0] * c.dir[0] + c.dir[1] * c.dir[1] + c.dir[2] * c.dir[2]);
      if (dirNorm > 0.0) {
        for (int d = 0; d < 3; ++d)
          c.dir[d] /= dirNorm;
      }

      // Normalize polarization vector if non-zero.
      const amrex::Real polNorm = std::sqrt(
          c.pol[0] * c.pol[0] + c.pol[1] * c.pol[1] + c.pol[2] * c.pol[2]);
      if (polNorm > 0.0) {
        for (int d = 0; d < 3; ++d)
          c.pol[d] /= polNorm;
      }

      f.comps.push_back(c);
    }
  }
}

void WaveBoundaryManager::convert_units(const FluidInterface* fi) {
  if (window.active && !window.isCodeUnits) {
    if (siInput && fi && fi->get_Si2NoL() > 0.0) {
      window.lo *= fi->get_Si2NoL();
      window.hi *= fi->get_Si2NoL();
      window.taper *= fi->get_Si2NoL();
      window.taperLo *= fi->get_Si2NoL();
      window.taperHi *= fi->get_Si2NoL();
    }
    window.isCodeUnits = true;
  }

  if (turb.active && fi) {
    const amrex::Real si2noT = fi->get_Si2NoT();
    turb.tStart = (si2noT > 0.0) ? (turb.tStartSI * si2noT) : turb.tStartSI;
    turb.rampTime =
        (si2noT > 0.0) ? (turb.rampTimeSI * si2noT) : turb.rampTimeSI;
  }

  for (auto& f : faces) {
    for (auto& c : f.comps) {
      if (c.isCodeUnits)
        continue;

      // SI -> code conversion (unless the block was marked 'code').
      if (siInput && fi) {
        const amrex::Real f_si2no = si2no(fi, c.iField);
        c.amplitude *= f_si2no;
        if (fi->get_Si2NoL() > 0.0)
          c.frequency *= fi->get_Si2NoV() / fi->get_Si2NoL();
        if (c.waveLength > 0.0)
          c.waveLength *= fi->get_Si2NoL();

        const amrex::Real si2noT = fi->get_Si2NoT();
        if (si2noT > 0.0) {
          c.tStart *= si2noT;
          if (c.tEnd > 0.0)
            c.tEnd *= si2noT;
          c.rampTime *= si2noT;
          c.tCenter *= si2noT;
          c.tWidth *= si2noT;
        }
      }

      // Precompute wave vector k_vec = (2*pi/lambda) * dir.
      if (c.waveLength > 0.0) {
        const amrex::Real k = cTwoPi / c.waveLength;
        for (int d = 0; d < 3; ++d)
          c.k_vec[d] = k * c.dir[d];
      } else {
        for (int d = 0; d < 3; ++d)
          c.k_vec[d] = 0.0;
      }

      // Amplitude guard: reject nonlinear/unstable injection.
      if (maxAmplitude > 0.0 && std::abs(c.amplitude) > maxAmplitude) {
        amrex::Abort("WAVEBC: amplitude exceeds maxAmplitude (code units); "
                     "injection may be nonlinear/unstable.");
      }

      c.isCodeUnits = true;
    }
  }
}

void WaveBoundaryManager::ensure_turbulence_initialized(
    const FluidInterface* fi, const amrex::Geometry& geom, int nDim,
    bool isFake2D, const bool isWaveFace[6]) {
  if (!turb.active || turb.initialized)
    return;
  if (!fi)
    return;

  const amrex::Real si2noT = fi->get_Si2NoT();
  turb.tStart = (si2noT > 0.0) ? (turb.tStartSI * si2noT) : turb.tStartSI;
  turb.rampTime = (si2noT > 0.0) ? (turb.rampTimeSI * si2noT) : turb.rampTimeSI;

  const int nS = fi->get_nS();
  const bool useHybrid = fi->get_useHybridPIC();
  const bool swmfSingleFluid = fi->is_init_from_swmf();

  amrex::Real mIon = 0.0, mElec = 0.0;
  for (int iS = 0; iS < nS; ++iS) {
    const amrex::Real q = fi->get_species_charge(iS);
    const amrex::Real m = fi->get_species_mass(iS);
    if (q > 0.0 && mIon <= 0.0)
      mIon = m;
    else if (q < 0.0 && mElec <= 0.0 && !useHybrid)
      mElec = m;
  }
  if (mIon <= 0.0)
    mIon = 1.0;
  const amrex::Real mSum = mIon + mElec;
  turb.massRatioIon = mIon / mSum;
  turb.massRatioElec = mElec / mSum;

  // Helper to evaluate node state (rho, rhoQ, V, B) from FluidInterface
  auto eval_node_state = [&](const amrex::MFIter& mfi,
                             const amrex::IntVect& ijk, amrex::Real& rhoOut,
                             amrex::Real& rhoQOut, amrex::Real vOut[3],
                             amrex::Real bOut[3]) -> bool {
    rhoOut = 0.0;
    rhoQOut = 0.0;
    vOut[0] = vOut[1] = vOut[2] = 0.0;
    bOut[0] = fi->get_bx(mfi, ijk, 0);
    bOut[1] = fi->get_by(mfi, ijk, 0);
    bOut[2] = fi->get_bz(mfi, ijk, 0);

    amrex::Real rhoIon = 0.0, rhoAll = 0.0;
    amrex::Real pxIon = 0.0, pyIon = 0.0, pzIon = 0.0;
    for (int iS = 0; iS < nS; ++iS) {
      const amrex::Real q = fi->get_species_charge(iS);
      const amrex::Real m = fi->get_species_mass(iS);
      if (useHybrid && q < 0.0)
        continue;
      const amrex::Real n = fi->get_number_density(mfi, ijk, iS, 0);
      if (n <= 0.0)
        continue;
      const amrex::Real rhoS = n * m;
      rhoAll += rhoS;
      if (q > 0.0) {
        rhoIon += rhoS;
        rhoQOut += n * q;
        pxIon += rhoS * fi->get_ux(mfi, ijk, iS, 0);
        pyIon += rhoS * fi->get_uy(mfi, ijk, iS, 0);
        pzIon += rhoS * fi->get_uz(mfi, ijk, iS, 0);
      }
    }
    rhoOut = swmfSingleFluid ? rhoAll : ((rhoIon > 0.0) ? rhoIon : rhoAll);
    if (rhoIon > 0.0) {
      vOut[0] = pxIon / rhoIon;
      vOut[1] = pyIon / rhoIon;
      vOut[2] = pzIon / rhoIon;
    }
    return rhoOut > 0.0;
  };

  // Sample all active 'wave' faces in FluidInterface to identify the primary
  // upstream inflow face and measure the reference solar-wind state.
  constexpr int nSlotPerFace = 9; // count, rho, rhoQ, vx, vy, vz, bx, by, bz
  std::vector<amrex::Real> faceSums(6 * nSlotPerFace, 0.0);

  const int effDim = (isFake2D && nDim == 3) ? 2 : nDim;
  const amrex::Box domNode =
      amrex::convert(geom.Domain(), amrex::IntVect::TheNodeVector());

  if (fi->is_node_fluid_ready() && !fi->get_node_fluid(0).empty()) {
    const amrex::MultiFab& nodeFluid = fi->get_node_fluid(0);
    for (int d = 0; d < effDim; ++d) {
      for (int side = 0; side < 2; ++side) {
        const int fIdx = 2 * d + side;
        if (!isWaveFace[fIdx])
          continue;
        const int wallIdx =
            (side == 0) ? domNode.smallEnd(d) : domNode.bigEnd(d);
        amrex::Box faceBox = domNode;
        faceBox.setSmall(d, wallIdx);
        faceBox.setBig(d, wallIdx);
        if (isFake2D && nDim == 3) {
          faceBox.setSmall(2, domNode.smallEnd(2));
          faceBox.setBig(2, domNode.smallEnd(2));
        }

        for (amrex::MFIter mfi(nodeFluid); mfi.isValid(); ++mfi) {
          const amrex::Box& vbx = mfi.validbox();
          amrex::Box sub = vbx & faceBox;
          if (!sub.ok())
            continue;
          const auto lo = amrex::lbound(sub);
          const auto hi = amrex::ubound(sub);
          for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                const amrex::IntVect ijk{ AMREX_D_DECL(i, j, k) };
                // Unique nodal ownership across tangential block boundaries:
                bool uniqueOwner = true;
                for (int od = 0; od < effDim; ++od) {
                  if (od != d && ijk[od] == vbx.bigEnd(od) &&
                      ijk[od] != domNode.bigEnd(od)) {
                    uniqueOwner = false;
                    break;
                  }
                }
                if (!uniqueOwner)
                  continue;

                amrex::Real rho = 0.0, rhoQ = 0.0, v[3], b[3];
                if (eval_node_state(mfi, ijk, rho, rhoQ, v, b)) {
                  amrex::Real* s = &faceSums[fIdx * nSlotPerFace];
                  s[0] += 1.0;
                  s[1] += rho;
                  s[2] += rhoQ;
                  s[3] += v[0];
                  s[4] += v[1];
                  s[5] += v[2];
                  s[6] += b[0];
                  s[7] += b[1];
                  s[8] += b[2];
                }
              }
            }
          }
        }
      }
    }
    amrex::ParallelDescriptor::ReduceRealSum(faceSums.data(),
                                             static_cast<int>(faceSums.size()));
  }

  int bestFace = -1;
  amrex::Real bestInflowSpeed = -1.0e30;
  for (int fIdx = 0; fIdx < 2 * effDim; ++fIdx) {
    if (!isWaveFace[fIdx])
      continue;
    const amrex::Real* s = &faceSums[fIdx * nSlotPerFace];
    if (s[0] <= 0.0)
      continue;
    const int d = fIdx / 2;
    const int side = fIdx % 2;
    const amrex::Real vd = s[3 + d] / s[0];
    const amrex::Real vInward = (side == 1) ? -vd : vd;
    if (vInward > bestInflowSpeed) {
      bestInflowSpeed = vInward;
      bestFace = fIdx;
    }
  }

  if (bestFace >= 0) {
    const amrex::Real* s = &faceSums[bestFace * nSlotPerFace];
    const amrex::Real invN = 1.0 / s[0];
    turb.primaryDir = bestFace / 2;
    turb.primarySide = bestFace % 2;
    turb.rho0 = s[1] * invN;
    turb.rhoQ0 = s[2] * invN;
    turb.v0[0] = s[3] * invN;
    turb.v0[1] = s[4] * invN;
    turb.v0[2] = s[5] * invN;
    turb.b0[0] = s[6] * invN;
    turb.b0[1] = s[7] * invN;
    turb.b0[2] = s[8] * invN;
  } else {
    // Fallback to #INFLOW / #UNIFORMSTATE if nodeFluid was not populated on
    // boundary nodes
    const auto& uState = fi->get_uniform_state();
    if (nS > 0 && static_cast<int>(uState.size()) >= 5 * nS + 3 &&
        static_cast<int>(fi->inflowState.size()) >= nS) {
      amrex::Real rhoIon = 0.0, rhoQ = 0.0;
      amrex::Real px = 0.0, py = 0.0, pz = 0.0;
      for (int iS = 0; iS < nS; ++iS) {
        const amrex::Real q = fi->get_species_charge(iS);
        const amrex::Real m = fi->get_species_mass(iS);
        const amrex::Real n = fi->inflowState[iS].nDens;
        if (q > 0.0 && n > 0.0) {
          rhoIon += n * m;
          rhoQ += n * q;
          px += n * m * fi->inflowState[iS].ux;
          py += n * m * fi->inflowState[iS].uy;
          pz += n * m * fi->inflowState[iS].uz;
        }
      }
      turb.rho0 = rhoIon;
      turb.rhoQ0 = rhoQ;
      if (rhoIon > 0.0) {
        turb.v0[0] = px / rhoIon;
        turb.v0[1] = py / rhoIon;
        turb.v0[2] = pz / rhoIon;
      }
      for (int d = 0; d < 3; ++d)
        turb.b0[d] = uState[5 * nS + d] * fi->get_Si2NoB();
    }
  }

  if (turb.rho0 <= 0.0) {
    amrex::Abort("TURBULENTINFLOW: failed to determine a positive upstream "
                 "density rho0 on the 'wave' boundary faces.");
  }

  const amrex::Real invSqrtFourPiRho = 1.0 / std::sqrt(fourPI * turb.rho0);
  for (int d = 0; d < 3; ++d)
    turb.va0[d] = turb.b0[d] * invSqrtFourPiRho;

  // Automatic solar-wind windowing on side 'wave' faces: if side wave faces
  // cross a bow shock (density jump > 1.20 * rho0) and no explicit
  // #WAVEBCWINDOW was given, automatically taper wave injection upstream of the
  // flank shock crossing.
  bool hasSideWaveFace = false;
  for (int fIdx = 0; fIdx < 2 * effDim; ++fIdx) {
    if (isWaveFace[fIdx] && (fIdx / 2) != turb.primaryDir)
      hasSideWaveFace = true;
  }

  if (hasSideWaveFace && !window.active && fi->is_node_fluid_ready() &&
      !fi->get_node_fluid(0).empty()) {
    const int pDir = turb.primaryDir;
    const int pSide = turb.primarySide;
    const amrex::Real* probLo = geom.ProbLo();
    const amrex::Real* probHi = geom.ProbHi();
    const amrex::Real* dx = geom.CellSize();
    const amrex::Real rhoThresh = 1.20 * turb.rho0;

    amrex::Real xShock = (pSide == 1) ? -1.0e30 : 1.0e30;
    const amrex::MultiFab& nodeFluid = fi->get_node_fluid(0);

    for (int d = 0; d < effDim; ++d) {
      if (d == pDir)
        continue;
      for (int side = 0; side < 2; ++side) {
        const int fIdx = 2 * d + side;
        if (!isWaveFace[fIdx])
          continue;
        const int wallIdx =
            (side == 0) ? domNode.smallEnd(d) : domNode.bigEnd(d);
        amrex::Box faceBox = domNode;
        faceBox.setSmall(d, wallIdx);
        faceBox.setBig(d, wallIdx);
        if (isFake2D && nDim == 3) {
          faceBox.setSmall(2, domNode.smallEnd(2));
          faceBox.setBig(2, domNode.smallEnd(2));
        }

        for (amrex::MFIter mfi(nodeFluid); mfi.isValid(); ++mfi) {
          amrex::Box sub = mfi.validbox() & faceBox;
          if (!sub.ok())
            continue;
          const auto lo = amrex::lbound(sub);
          const auto hi = amrex::ubound(sub);
          for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                const amrex::IntVect ijk{ AMREX_D_DECL(i, j, k) };
                amrex::Real rho = 0.0, rhoQ = 0.0, v[3], b[3];
                if (eval_node_state(mfi, ijk, rho, rhoQ, v, b) &&
                    rho > rhoThresh) {
                  const amrex::Real xCoord =
                      probLo[pDir] + ijk[pDir] * dx[pDir];
                  if (pSide == 1)
                    xShock = std::max(xShock, xCoord);
                  else
                    xShock = std::min(xShock, xCoord);
                }
              }
            }
          }
        }
      }
    }

    if (pSide == 1) {
      amrex::ParallelDescriptor::ReduceRealMax(xShock);
      if (xShock > probLo[pDir] && xShock < probHi[pDir]) {
        const amrex::Real lUp = probHi[pDir] - xShock;
        window.active = true;
        window.isCodeUnits = true;
        window.direction = pDir;
        window.lo = xShock + 0.05 * lUp;
        window.hi = 1.0e30;
        window.taperLo = 0.25 * lUp;
        window.taperHi = 0.0;
        window.taper = window.taperLo;
      }
    } else {
      amrex::ParallelDescriptor::ReduceRealMin(xShock);
      if (xShock > probLo[pDir] && xShock < probHi[pDir]) {
        const amrex::Real lUp = xShock - probLo[pDir];
        window.active = true;
        window.isCodeUnits = true;
        window.direction = pDir;
        window.lo = -1.0e30;
        window.hi = xShock - 0.05 * lUp;
        window.taperLo = 0.0;
        window.taperHi = 0.25 * lUp;
        window.taper = window.taperHi;
      }
    }
  }

  // Synthesize the N modes deterministically across all ranks
  const int N = turb.nMode;
  turb.kx.resize(N, 0.0);
  turb.ky.resize(N, 0.0);
  turb.kz.resize(N, 0.0);
  turb.Bx.resize(N, 0.0);
  turb.By.resize(N, 0.0);
  turb.Bz.resize(N, 0.0);
  turb.Ax.resize(N, 0.0);
  turb.Ay.resize(N, 0.0);
  turb.Az.resize(N, 0.0);
  turb.Jx.resize(N, 0.0);
  turb.Jy.resize(N, 0.0);
  turb.Jz.resize(N, 0.0);
  turb.Vmx.resize(N, 0.0);
  turb.Vmy.resize(N, 0.0);
  turb.Vmz.resize(N, 0.0);
  turb.omega.resize(N, 0.0);
  turb.phase.resize(N, 0.0);
  turb.signProp.resize(N, 1.0);

  std::mt19937_64 rng(static_cast<std::uint64_t>(turb.seed));
  std::uniform_real_distribution<double> uniTwoPi(0.0, cTwoPi);

  std::vector<amrex::Real> kmag(N, 0.0);
  std::vector<amrex::Real> theta(N, 0.0);
  std::vector<amrex::Real> polAngle(N, 0.0);

  for (int m = 0; m < N; ++m) {
    const amrex::Real frac =
        (N > 1)
            ? (static_cast<amrex::Real>(m) / static_cast<amrex::Real>(N - 1))
            : 0.0;
    const amrex::Real lam =
        turb.lambdaMax * std::pow(turb.lambdaMin / turb.lambdaMax, frac);
    kmag[m] = cTwoPi / lam;
  }
  for (int m = 0; m < N; ++m)
    theta[m] = static_cast<amrex::Real>(uniTwoPi(rng));
  for (int m = 0; m < N; ++m)
    turb.phase[m] = static_cast<amrex::Real>(uniTwoPi(rng));
  if (turb.polarization == "isotropic") {
    for (int m = 0; m < N; ++m)
      polAngle[m] = static_cast<amrex::Real>(uniTwoPi(rng));
  }

  // Outward propagation sign in the plasma frame: follows V0 along B0
  const amrex::Real sOut = (dot3(turb.v0, turb.b0) >= 0.0) ? 1.0 : -1.0;

  // Spectral amplitude weights w_m ~ k_m^((1 - spectralIndex)/2)
  std::vector<amrex::Real> rawWeight(N, 0.0);
  amrex::Real sumSqPlus = 0.0, sumSqMinus = 0.0;
  for (int m = 0; m < N; ++m) {
    rawWeight[m] = std::pow(kmag[m], 0.5 * (1.0 - turb.spectralIndex));
    if (turb.sigmaC >= 1.0 - 1.0e-12) {
      turb.signProp[m] = sOut;
      sumSqPlus += rawWeight[m] * rawWeight[m];
    } else if (turb.sigmaC <= -1.0 + 1.0e-12) {
      turb.signProp[m] = -sOut;
      sumSqMinus += rawWeight[m] * rawWeight[m];
    } else {
      // Interleave outward and inward modes across the spectrum so both
      // propagation directions span all wavelengths, then normalize each group
      // to match the exact cross-helicity sigmaC = (E+ - E-) / (E+ + E-).
      if (m % 2 == 0) {
        turb.signProp[m] = sOut;
        sumSqPlus += rawWeight[m] * rawWeight[m];
      } else {
        turb.signProp[m] = -sOut;
        sumSqMinus += rawWeight[m] * rawWeight[m];
      }
    }
  }

  const amrex::Real b0Mag = norm3(turb.b0);
  const amrex::Real targetRms = turb.dBrmsOverB0 * b0Mag;
  const amrex::Real fPlus = 0.5 * (1.0 + turb.sigmaC);
  const amrex::Real fMinus = 0.5 * (1.0 - turb.sigmaC);
  const amrex::Real scalePlus =
      (sumSqPlus > 0.0) ? (targetRms * std::sqrt(fPlus / (0.5 * sumSqPlus)))
                        : 0.0;
  const amrex::Real scaleMinus =
      (sumSqMinus > 0.0) ? (targetRms * std::sqrt(fMinus / (0.5 * sumSqMinus)))
                         : 0.0;

  const amrex::Real* probLo = geom.ProbLo();
  const amrex::Real* probHi = geom.ProbHi();
  const bool perY = (effDim > 1) && geom.isPeriodic(1);
  const amrex::Real ly = perY ? (probHi[1] - probLo[1]) : 0.0;
  const amrex::Real dky = (ly > 0.0) ? (cTwoPi / ly) : 0.0;
  const amrex::Real invFourPi = 1.0 / fourPI;

  for (int m = 0; m < N; ++m) {
    amrex::Real kx = kmag[m] * std::cos(theta[m]);
    amrex::Real ky = kmag[m] * std::sin(theta[m]);
    amrex::Real kz = 0.0;

    if (perY && dky > 0.0) {
      int mult = static_cast<int>(std::llround(ky / dky));
      while (std::abs(mult) * dky > kmag[m]) {
        mult -= (mult > 0) ? 1 : -1;
      }
      ky = static_cast<amrex::Real>(mult) * dky;
      const amrex::Real kxAbs =
          std::sqrt(std::max(kmag[m] * kmag[m] - ky * ky, amrex::Real(0.0)));
      kx = (kx >= 0.0) ? kxAbs : -kxAbs;
    }

    const amrex::Real kVec[3] = { kx, ky, kz };
    const amrex::Real kNorm = norm3(kVec);
    const amrex::Real kHat[3] = { (kNorm > 0.0) ? (kx / kNorm) : 1.0,
                                  (kNorm > 0.0) ? (ky / kNorm) : 0.0,
                                  (kNorm > 0.0) ? (kz / kNorm) : 0.0 };

    // Determine transverse unit polarization vector ePol perpendicular to kHat
    amrex::Real ePol[3] = { 0.0, 0.0, 1.0 };
    if (turb.polarization == "out_of_plane") {
      ePol[0] = 0.0;
      ePol[1] = 0.0;
      ePol[2] = 1.0;
    } else if (turb.polarization == "in_plane") {
      ePol[0] = -kHat[1];
      ePol[1] = kHat[0];
      ePol[2] = 0.0;
    } else {
      // "alfven" (k x B0) or "isotropic"
      amrex::Real e1[3];
      cross3(kHat, turb.b0, e1);
      amrex::Real e1Len = norm3(e1);
      if (e1Len <= 1.0e-12 * b0Mag) {
        e1[0] = -kHat[1];
        e1[1] = kHat[0];
        e1[2] = 0.0;
        e1Len = norm3(e1);
      }
      for (int d = 0; d < 3; ++d)
        e1[d] /= e1Len;
      if (turb.polarization == "alfven") {
        for (int d = 0; d < 3; ++d)
          ePol[d] = e1[d];
      } else {
        amrex::Real e2[3];
        cross3(kHat, e1, e2);
        const amrex::Real ca = std::cos(polAngle[m]);
        const amrex::Real sa = std::sin(polAngle[m]);
        for (int d = 0; d < 3; ++d)
          ePol[d] = ca * e1[d] + sa * e2[d];
      }
    }

    const amrex::Real ampB =
        rawWeight[m] *
        ((turb.signProp[m] * sOut > 0.0) ? scalePlus : scaleMinus);
    const amrex::Real bVec[3] = { ampB * ePol[0], ampB * ePol[1],
                                  ampB * ePol[2] };

    amrex::Real kCrossB[3];
    cross3(kVec, bVec, kCrossB);
    const amrex::Real invKSq = (kNorm > 0.0) ? (1.0 / (kNorm * kNorm)) : 0.0;
    const amrex::Real alphaV = -turb.signProp[m] * invSqrtFourPiRho;

    turb.kx[m] = kx;
    turb.ky[m] = ky;
    turb.kz[m] = kz;
    turb.Bx[m] = bVec[0];
    turb.By[m] = bVec[1];
    turb.Bz[m] = bVec[2];
    turb.Ax[m] = kCrossB[0] * invKSq;
    turb.Ay[m] = kCrossB[1] * invKSq;
    turb.Az[m] = kCrossB[2] * invKSq;
    turb.Jx[m] = kCrossB[0] * invFourPi;
    turb.Jy[m] = kCrossB[1] * invFourPi;
    turb.Jz[m] = kCrossB[2] * invFourPi;
    turb.Vmx[m] = alphaV * bVec[0];
    turb.Vmy[m] = alphaV * bVec[1];
    turb.Vmz[m] = alphaV * bVec[2];
    turb.omega[m] =
        dot3(kVec, turb.v0) + turb.signProp[m] * dot3(kVec, turb.va0);
  }

  turb.initialized = true;

  if (amrex::ParallelDescriptor::IOProcessor()) {
    const amrex::Real no2siL = fi->get_No2SiL();
    const amrex::Real no2siV = fi->get_No2SiV();
    const amrex::Real no2siB = fi->get_No2SiB();
    const amrex::Real no2siRho = fi->get_No2SiRho();
    const amrex::Real no2siT = fi->get_No2SiT();
    const amrex::Real vaMagSI = norm3(turb.va0) * no2siV;
    const amrex::Real b0MagSI = b0Mag * no2siB;

    amrex::Print() << "  #TURBULENTINFLOW initialized:\n"
                   << "    primary inflow face: dir = " << turb.primaryDir
                   << ", side = " << (turb.primarySide == 0 ? "lo" : "hi")
                   << "\n"
                   << "    upstream B0 = (" << turb.b0[0] * no2siB * 1e9 << ", "
                   << turb.b0[1] * no2siB * 1e9 << ", "
                   << turb.b0[2] * no2siB * 1e9
                   << ") [nT] (|B0| = " << b0MagSI * 1e9 << " nT)\n"
                   << "    upstream V0 = (" << turb.v0[0] * no2siV * 1e-3
                   << ", " << turb.v0[1] * no2siV * 1e-3 << ", "
                   << turb.v0[2] * no2siV * 1e-3
                   << ") [km/s], vA = " << vaMagSI * 1e-3 << " [km/s]\n"
                   << "    nMode = " << N << ", lambda = [" << turb.lambdaMin
                   << ", " << turb.lambdaMax << "] [code] (["
                   << turb.lambdaMin * no2siL * 1e-3 << ", "
                   << turb.lambdaMax * no2siL * 1e-3 << "] km)\n"
                   << "    spectralIndex = " << turb.spectralIndex
                   << ", dBrms/B0 = " << turb.dBrmsOverB0
                   << ", sigmaC = " << turb.sigmaC
                   << ", polarization = " << turb.polarization
                   << ", seed = " << turb.seed << "\n";
    if (window.active) {
      amrex::Print() << "    spatial window: dir = " << window.direction
                     << ", lo = " << window.lo << ", hi = " << window.hi
                     << ", taperLo = " << window.taperLo
                     << ", taperHi = " << window.taperHi << " [code]\n";
    }

    if (!turb.modeFile.empty() && turb.modeFile != "none") {
      std::ofstream ofs(turb.modeFile);
      if (ofs) {
        ofs << std::setprecision(10);
        const amrex::Real nCm3 =
            (turb.rho0 * no2siRho) / (cProtonMassSI * mIon * 1.0e6);
        ofs << "{\n"
            << " \"info\": {\n"
            << "  \"n\": " << nCm3 << ",\n"
            << "  \"mi\": " << mIon << ",\n"
            << "  \"b0\": [" << turb.b0[0] * no2siB << ", "
            << turb.b0[1] * no2siB << ", " << turb.b0[2] * no2siB << "],\n"
            << "  \"v0\": [" << turb.v0[0] * no2siV << ", "
            << turb.v0[1] * no2siV << ", " << turb.v0[2] * no2siV << "],\n"
            << "  \"va\": [" << turb.va0[0] * no2siV << ", "
            << turb.va0[1] * no2siV << ", " << turb.va0[2] * no2siV << "],\n"
            << "  \"rho\": " << turb.rho0 * no2siRho << ",\n"
            << "  \"s\": " << sOut * (turb.sigmaC >= 0.0 ? 1.0 : -1.0) << ",\n"
            << "  \"sigma_c\": " << turb.sigmaC << ",\n"
            << "  \"polarization\": \"" << turb.polarization << "\",\n"
            << "  \"seed\": " << turb.seed << ",\n"
            << "  \"slope\": " << turb.spectralIndex << ",\n"
            << "  \"dbrms\": " << turb.dBrmsOverB0 << ",\n"
            << "  \"lmin\": " << turb.lambdaMin * no2siL << ",\n"
            << "  \"lmax\": " << turb.lambdaMax * no2siL << "\n"
            << " },\n"
            << " \"args\": {\n"
            << "  \"ramp\": " << turb.rampTimeSI << ",\n"
            << "  \"tstart\": " << turb.tStartSI << "\n"
            << " },\n"
            << " \"modes\": [\n";
        for (int m = 0; m < N; ++m) {
          const amrex::Real kSI[3] = { turb.kx[m] / no2siL, turb.ky[m] / no2siL,
                                       turb.kz[m] / no2siL };
          const amrex::Real kMagSI = norm3(kSI);
          const amrex::Real lamSI = (kMagSI > 0.0) ? (cTwoPi / kMagSI) : 0.0;
          const amrex::Real bModeCode[3] = { turb.Bx[m], turb.By[m],
                                             turb.Bz[m] };
          const amrex::Real ampBSI = norm3(bModeCode) * no2siB;
          const amrex::Real vModeCode[3] = { turb.Vmx[m], turb.Vmy[m],
                                             turb.Vmz[m] };
          const amrex::Real ampVSI = norm3(vModeCode) * no2siV;
          const amrex::Real signV = -turb.signProp[m];
          const amrex::Real omegaSI =
              (no2siT > 0.0) ? (turb.omega[m] / no2siT) : turb.omega[m];
          ofs << "  {\n"
              << "   \"k\": [" << kSI[0] << ", " << kSI[1] << "],\n"
              << "   \"wavelength\": " << lamSI << ",\n"
              << "   \"dir\": [" << (kMagSI > 0.0 ? kSI[0] / kMagSI : 0.0)
              << ", " << (kMagSI > 0.0 ? kSI[1] / kMagSI : 0.0) << "],\n"
              << "   \"amp_b\": " << ampBSI << ",\n"
              << "   \"amp_v\": " << ampVSI << ",\n"
              << "   \"sign_v\": " << signV << ",\n"
              << "   \"omega\": " << omegaSI << ",\n"
              << "   \"phase\": " << turb.phase[m] << "\n"
              << "  }" << (m + 1 < N ? "," : "") << "\n";
        }
        ofs << " ]\n}\n";
      }
    }
  }
}

void WaveBoundaryManager::eval_turbulent_field(const amrex::Real* pos,
                                               const amrex::Real* phase_t,
                                               amrex::Real env, int iField,
                                               amrex::Real out[3]) const {
  out[0] = 0.0;
  out[1] = 0.0;
  out[2] = 0.0;
  if (!turb.active || !turb.initialized || env <= 0.0 ||
      turb.dBrmsOverB0 <= 0.0)
    return;

  amrex::Real w = 1.0;
  amrex::Real gradW[3] = { 0.0, 0.0, 0.0 };
  window.weight_and_grad(pos, w, gradW);
  if (w == 0.0 && gradW[0] == 0.0 && gradW[1] == 0.0 && gradW[2] == 0.0)
    return;

  const int N = turb.nMode;
  const amrex::Real* AMREX_RESTRICT kx = turb.kx.data();
  const amrex::Real* AMREX_RESTRICT ky = turb.ky.data();
  const amrex::Real* AMREX_RESTRICT kz = turb.kz.data();
  const amrex::Real* AMREX_RESTRICT Bx = turb.Bx.data();
  const amrex::Real* AMREX_RESTRICT By = turb.By.data();
  const amrex::Real* AMREX_RESTRICT Bz = turb.Bz.data();
  const amrex::Real* AMREX_RESTRICT Ax = turb.Ax.data();
  const amrex::Real* AMREX_RESTRICT Ay = turb.Ay.data();
  const amrex::Real* AMREX_RESTRICT Az = turb.Az.data();

  amrex::Real dBxRaw = 0.0, dByRaw = 0.0, dBzRaw = 0.0;
  amrex::Real dAxRaw = 0.0, dAyRaw = 0.0, dAzRaw = 0.0;

  if (iField == 0) {
// Magnetic field perturbation: dB = env * (W * dB_raw + grad(W) x dA_raw)
#pragma omp simd reduction(+ : dBxRaw, dByRaw, dBzRaw, dAxRaw, dAyRaw, dAzRaw)
    for (int m = 0; m < N; ++m) {
      const amrex::Real arg =
          kx[m] * pos[0] + ky[m] * pos[1] + kz[m] * pos[2] + phase_t[m];
      const amrex::Real s = std::sin(arg);
      const amrex::Real c = std::cos(arg);
      dBxRaw += Bx[m] * s;
      dByRaw += By[m] * s;
      dBzRaw += Bz[m] * s;
      dAxRaw += Ax[m] * c;
      dAyRaw += Ay[m] * c;
      dAzRaw += Az[m] * c;
    }
    const amrex::Real dA[3] = { dAxRaw, dAyRaw, dAzRaw };
    amrex::Real gradWCrossA[3];
    cross3(gradW, dA, gradWCrossA);
    out[0] = env * (w * dBxRaw + gradWCrossA[0]);
    out[1] = env * (w * dByRaw + gradWCrossA[1]);
    out[2] = env * (w * dBzRaw + gradWCrossA[2]);
  } else if (iField == 1) {
    // Electric field perturbation from full kinetic consistency:
    //   dV_e = dV_mhd - (m_i / (m_i + m_e)) * dJ / rho_q0
    //   dE   = -(V0 x dB + dV_e x B0 + dV_e x dB)
    const amrex::Real* AMREX_RESTRICT Jx = turb.Jx.data();
    const amrex::Real* AMREX_RESTRICT Jy = turb.Jy.data();
    const amrex::Real* AMREX_RESTRICT Jz = turb.Jz.data();
    const amrex::Real* AMREX_RESTRICT Vmx = turb.Vmx.data();
    const amrex::Real* AMREX_RESTRICT Vmy = turb.Vmy.data();
    const amrex::Real* AMREX_RESTRICT Vmz = turb.Vmz.data();

    amrex::Real dJxRaw = 0.0, dJyRaw = 0.0, dJzRaw = 0.0;
    amrex::Real dVxRaw = 0.0, dVyRaw = 0.0, dVzRaw = 0.0;

#pragma omp simd reduction(+ : dBxRaw, dByRaw, dBzRaw, dAxRaw, dAyRaw, dAzRaw, \
                               dJxRaw, dJyRaw, dJzRaw, dVxRaw, dVyRaw, dVzRaw)
    for (int m = 0; m < N; ++m) {
      const amrex::Real arg =
          kx[m] * pos[0] + ky[m] * pos[1] + kz[m] * pos[2] + phase_t[m];
      const amrex::Real s = std::sin(arg);
      const amrex::Real c = std::cos(arg);
      dBxRaw += Bx[m] * s;
      dByRaw += By[m] * s;
      dBzRaw += Bz[m] * s;
      dAxRaw += Ax[m] * c;
      dAyRaw += Ay[m] * c;
      dAzRaw += Az[m] * c;
      dJxRaw += Jx[m] * c;
      dJyRaw += Jy[m] * c;
      dJzRaw += Jz[m] * c;
      dVxRaw += Vmx[m] * s;
      dVyRaw += Vmy[m] * s;
      dVzRaw += Vmz[m] * s;
    }

    const amrex::Real dA[3] = { dAxRaw, dAyRaw, dAzRaw };
    amrex::Real gradWCrossA[3];
    cross3(gradW, dA, gradWCrossA);
    const amrex::Real dB[3] = { env * (w * dBxRaw + gradWCrossA[0]),
                                env * (w * dByRaw + gradWCrossA[1]),
                                env * (w * dBzRaw + gradWCrossA[2]) };

    const amrex::Real dBRawVec[3] = { dBxRaw, dByRaw, dBzRaw };
    amrex::Real gradWCrossB[3];
    cross3(gradW, dBRawVec, gradWCrossB);
    const amrex::Real invFourPi = 1.0 / fourPI;
    const amrex::Real dJ[3] = { env * (w * dJxRaw + invFourPi * gradWCrossB[0]),
                                env * (w * dJyRaw + invFourPi * gradWCrossB[1]),
                                env *
                                    (w * dJzRaw + invFourPi * gradWCrossB[2]) };

    const amrex::Real jFactorElec =
        (turb.rhoQ0 > 0.0) ? (turb.massRatioIon / turb.rhoQ0) : 0.0;
    const amrex::Real dVe[3] = { env * w * dVxRaw - jFactorElec * dJ[0],
                                 env * w * dVyRaw - jFactorElec * dJ[1],
                                 env * w * dVzRaw - jFactorElec * dJ[2] };

    // dE = -(V0 x dB + dVe x B0 + dVe x dB) = dB x V0 + (B0 + dB) x dVe
    const amrex::Real bTot[3] = { turb.b0[0] + dB[0], turb.b0[1] + dB[1],
                                  turb.b0[2] + dB[2] };
    amrex::Real term1[3], term2[3];
    cross3(dB, turb.v0, term1);
    cross3(bTot, dVe, term2);
    out[0] = term1[0] + term2[0];
    out[1] = term1[1] + term2[1];
    out[2] = term1[2] + term2[2];
  }
}

void WaveBoundaryManager::eval_turbulent_velocity_kick(
    const amrex::Real* pos, amrex::Real t, amrex::Real charge, amrex::Real& dvx,
    amrex::Real& dvy, amrex::Real& dvz) const {
  if (!turb.active || !turb.initialized || turb.dBrmsOverB0 <= 0.0)
    return;
  const amrex::Real env = turbulent_envelope(t);
  if (env <= 0.0)
    return;

  amrex::Real w = 1.0;
  amrex::Real gradW[3] = { 0.0, 0.0, 0.0 };
  window.weight_and_grad(pos, w, gradW);
  if (w == 0.0 && gradW[0] == 0.0 && gradW[1] == 0.0 && gradW[2] == 0.0)
    return;

  const int N = turb.nMode;
  const amrex::Real* AMREX_RESTRICT kx = turb.kx.data();
  const amrex::Real* AMREX_RESTRICT ky = turb.ky.data();
  const amrex::Real* AMREX_RESTRICT kz = turb.kz.data();
  const amrex::Real* AMREX_RESTRICT Bx = turb.Bx.data();
  const amrex::Real* AMREX_RESTRICT By = turb.By.data();
  const amrex::Real* AMREX_RESTRICT Bz = turb.Bz.data();
  const amrex::Real* AMREX_RESTRICT Jx = turb.Jx.data();
  const amrex::Real* AMREX_RESTRICT Jy = turb.Jy.data();
  const amrex::Real* AMREX_RESTRICT Jz = turb.Jz.data();
  const amrex::Real* AMREX_RESTRICT Vmx = turb.Vmx.data();
  const amrex::Real* AMREX_RESTRICT Vmy = turb.Vmy.data();
  const amrex::Real* AMREX_RESTRICT Vmz = turb.Vmz.data();
  const amrex::Real* AMREX_RESTRICT omega = turb.omega.data();
  const amrex::Real* AMREX_RESTRICT phase = turb.phase.data();

  amrex::Real dBxRaw = 0.0, dByRaw = 0.0, dBzRaw = 0.0;
  amrex::Real dJxRaw = 0.0, dJyRaw = 0.0, dJzRaw = 0.0;
  amrex::Real dVxRaw = 0.0, dVyRaw = 0.0, dVzRaw = 0.0;

#pragma omp simd reduction(+ : dBxRaw, dByRaw, dBzRaw, dJxRaw, dJyRaw, dJzRaw, \
                               dVxRaw, dVyRaw, dVzRaw)
  for (int m = 0; m < N; ++m) {
    const amrex::Real arg = kx[m] * pos[0] + ky[m] * pos[1] + kz[m] * pos[2] -
                            omega[m] * t + phase[m];
    const amrex::Real s = std::sin(arg);
    const amrex::Real c = std::cos(arg);
    dBxRaw += Bx[m] * s;
    dByRaw += By[m] * s;
    dBzRaw += Bz[m] * s;
    dJxRaw += Jx[m] * c;
    dJyRaw += Jy[m] * c;
    dJzRaw += Jz[m] * c;
    dVxRaw += Vmx[m] * s;
    dVyRaw += Vmy[m] * s;
    dVzRaw += Vmz[m] * s;
  }

  const amrex::Real dBRawVec[3] = { dBxRaw, dByRaw, dBzRaw };
  amrex::Real gradWCrossB[3];
  cross3(gradW, dBRawVec, gradWCrossB);
  const amrex::Real invFourPi = 1.0 / fourPI;
  const amrex::Real dJ[3] = { env * (w * dJxRaw + invFourPi * gradWCrossB[0]),
                              env * (w * dJyRaw + invFourPi * gradWCrossB[1]),
                              env * (w * dJzRaw + invFourPi * gradWCrossB[2]) };

  amrex::Real jDriftScale = 0.0;
  if (turb.rhoQ0 > 0.0) {
    if (charge < 0.0) {
      jDriftScale = -turb.massRatioIon / turb.rhoQ0;
    } else if (charge > 0.0) {
      jDriftScale = turb.massRatioElec / turb.rhoQ0;
    }
  }

  dvx += env * w * dVxRaw + jDriftScale * dJ[0];
  dvy += env * w * dVyRaw + jDriftScale * dJ[1];
  dvz += env * w * dVzRaw + jDriftScale * dJ[2];
}
