#ifndef FLEKS_WAVEBC_H
#define FLEKS_WAVEBC_H

#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "ReadParam.h"

class FluidInterface;

// One spectral component of a wave injected at a boundary face.
struct WaveComponent {
  enum {
    kMono = 0,    // A sin(k.x - w t + phi)
    kPacket = 1,  // mono * Gaussian envelope
    kCustom = 2,  // std::function hook
    kSpectral = 3 // superposition of sub-components
  };

  int iField = 0;    // 0 = B, 1 = E, 2 = velocity, 3 = density weight
  int iSpecies = -1; // species for velocity/density; -1 = EM field

  amrex::Real amplitude = 0.0;
  amrex::Real frequency = 0.0;  // omega
  amrex::Real waveLength = 0.0; // 0 => pump (uniform along dir)
  amrex::Real phase = 0.0;
  amrex::Real dir[3] = { 1, 0, 0 };

  amrex::Real tStart = 0.0;
  amrex::Real tEnd = -1.0;    // -1 => forever
  amrex::Real rampTime = 0.0; // cos^2 ramp-in
  amrex::Real tCenter = 0.0;  // packet center (kPacket)
  amrex::Real tWidth = 1.0;   // packet width  (kPacket)

  int profile = kMono;

  std::vector<WaveComponent> comps; // sub-components (kSpectral)
  int seed = 1;                     // random-phase seed (kSpectral)
  amrex::Real spectralIndex = 0.0;  // power-law index (unused in v1)

  std::function<amrex::Real(const WaveComponent&, amrex::Real,
                            const amrex::Real*)>
      custom = nullptr; // kCustom

  amrex::Real pol[3] = { 0, 1, 0 };   // polarisation vector
  amrex::Real k_vec[3] = { 0, 0, 0 }; // precomputed (2*pi/waveLength) * dir

  bool isCodeUnits = false; // set once convert_units() has processed it
};

// Smooth spatial window (#WAVEBCWINDOW or auto-detected upstream solar-wind
// window for #TURBULENTINFLOW): 1 inside [lo, hi] along one axis, with
// independent cos^2 tapers of width taperLo / taperHi at each end, 0 outside.
struct WaveWindow {
  bool active = false;
  int direction = 0; // 0,1,2 = x,y,z
  amrex::Real lo = 0.0;
  amrex::Real hi = 0.0;
  amrex::Real taper = 0.0;
  amrex::Real taperLo = 0.0;
  amrex::Real taperHi = 0.0;
  bool isCodeUnits = false;

  amrex::Real weight(const amrex::Real* pos) const;
  void weight_and_grad(const amrex::Real* pos, amrex::Real& w,
                       amrex::Real gradW[3]) const;
};

// Broadband Alfvenic / Elsasser turbulent inflow spectrum (#TURBULENTINFLOW).
// Synthesizes a divergence-free wave spectrum from physical parameters and the
// upstream reference state in FluidInterface, evaluating dB, dE, dV_i, and dV_e
// self-consistently in a single SIMD pass.
struct TurbulentInflowConfig {
  bool active = false;
  bool initialized = false;

  int nMode = 48;
  amrex::Real lambdaMin = 0.5; // domain length units (code units)
  amrex::Real lambdaMax = 8.0; // domain length units (code units)
  amrex::Real spectralIndex = 5.0 / 3.0;
  amrex::Real dBrmsOverB0 = 0.3;
  int seed = 1;
  amrex::Real rampTimeSI = 20.0; // [s]
  amrex::Real tStartSI = 0.0;    // [s]
  amrex::Real sigmaC = 1.0;      // cross-helicity in [-1, 1]
  std::string polarization =
      "out_of_plane"; // out_of_plane | in_plane | alfven | isotropic
  std::string modeFile = "wave_modes.json";

  // Resolved code-unit reference state and timing:
  amrex::Real rampTime = 0.0;
  amrex::Real tStart = 0.0;
  amrex::Real rho0 = 0.0;          // upstream mass density
  amrex::Real rhoQ0 = 0.0;         // upstream ion charge density
  amrex::Real massRatioIon = 1.0;  // m_i / (m_i + m_e)
  amrex::Real massRatioElec = 0.0; // m_e / (m_i + m_e)
  amrex::Real v0[3] = { 0.0, 0.0, 0.0 };
  amrex::Real b0[3] = { 0.0, 0.0, 0.0 };
  amrex::Real va0[3] = { 0.0, 0.0, 0.0 };
  int primaryDir = 0;
  int primarySide = 1;

  // Flat Structure-of-Arrays (SoA) for SIMD mode summation (length nMode):
  std::vector<amrex::Real> kx, ky, kz;
  std::vector<amrex::Real> Bx, By, Bz; // B_n = amp_n * e_n
  std::vector<amrex::Real> Ax, Ay, Az; // A_n = (k_n x B_n) / |k_n|^2
  std::vector<amrex::Real> Jx, Jy, Jz; // J_n = (k_n x B_n) / (4*pi)
  std::vector<amrex::Real> Vmx, Vmy,
      Vmz; // V_mhd,n = -s_n * B_n / sqrt(4*pi*rho0)
  std::vector<amrex::Real> omega;
  std::vector<amrex::Real> phase;
  std::vector<amrex::Real> signProp;
};

// Per-face emitter stack.
struct WaveFace {
  int direction = 0; // 0,1,2 = x,y,z
  int side = 0;      // 0 = lo, 1 = hi

  std::vector<WaveComponent> comps;
};

class WaveProfile {
public:
  virtual ~WaveProfile() = default;
  virtual amrex::Real value(const WaveComponent& c, amrex::Real t,
                            const amrex::Real* pos) const = 0;
};

class MonoWave : public WaveProfile {
public:
  amrex::Real value(const WaveComponent& c, amrex::Real t,
                    const amrex::Real* pos) const override;
};

class WavePacket : public MonoWave {
public:
  amrex::Real value(const WaveComponent& c, amrex::Real t,
                    const amrex::Real* pos) const override;
};

class SpectralWave : public WaveProfile {
public:
  amrex::Real value(const WaveComponent& c, amrex::Real t,
                    const amrex::Real* pos) const override;
};

class WaveRegistry {
public:
  static std::unique_ptr<WaveProfile> create(const WaveComponent& c);
};

// Runtime-editable wave-boundary manager held by Pic.
class WaveBoundaryManager {
public:
  bool active = false;
  std::vector<WaveFace> faces; // indexed by (direction, side)

  // Parse a #WAVEBC block. Values stay in input units until convert_units(),
  // because the normalization is only final after all commands are read.
  void read_param(ReadParam& param);
  // Parse a #WAVEBCWINDOW block.
  void read_window(ReadParam& param);
  // Parse a #TURBULENTINFLOW block.
  void read_turbulent_inflow(ReadParam& param);
  // Convert newly read components (and the window) to code units and
  // precompute k_vec.
  void convert_units(const FluidInterface* fi);
  // Lazily sample upstream solar-wind state from FluidInterface, auto-detect
  // any flank bow-shock window, and synthesize the #TURBULENTINFLOW mode table.
  void ensure_turbulence_initialized(const FluidInterface* fi,
                                     const amrex::Geometry& geom, int nDim,
                                     bool isFake2D, const bool isWaveFace[6]);
  // Window weight at pos (code units); 1 when no window is set.
  amrex::Real window_weight(const amrex::Real* pos) const {
    return window.active ? window.weight(pos) : 1.0;
  }
  amrex::Real value(const WaveComponent& c, amrex::Real t,
                    const amrex::Real* pos) const;

  // Evaluate #TURBULENTINFLOW time envelope at code time t.
  amrex::Real turbulent_envelope(amrex::Real t) const;
  // Evaluate #TURBULENTINFLOW field perturbation (iField: 0 = dB, 1 = dE)
  // using precomputed phase_t[m] = phase[m] - omega[m] * t and time envelope.
  void eval_turbulent_field(const amrex::Real* pos, const amrex::Real* phase_t,
                            amrex::Real env, int iField,
                            amrex::Real out[3]) const;
  // Evaluate #TURBULENTINFLOW species velocity kick (charge > 0 => ion kick,
  // charge < 0 => electron kick including current drift -J / rho_q0).
  void eval_turbulent_velocity_kick(const amrex::Real* pos, amrex::Real t,
                                    amrex::Real charge, amrex::Real& dvx,
                                    amrex::Real& dvy, amrex::Real& dvz) const;

  WaveFace& face(int direction, int side);
  void clear();
  void add_component(int direction, int side, const WaveComponent& c);

  amrex::Real envelope(const WaveComponent& c, amrex::Real t) const;

  bool siInput = true;
  amrex::Real maxAmplitude = 0.0; // 0 => amplitude guard disabled

  WaveWindow window;
  TurbulentInflowConfig turb;
};

#endif // FLEKS_WAVEBC_H
