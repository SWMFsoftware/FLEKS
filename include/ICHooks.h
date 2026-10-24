#ifndef _IC_HOOKS_H_
#define _IC_HOOKS_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "Constants.h"

// Per-particle state handed to a particle-seeding IC (position, index,
// in/out bulk velocity and weight, thermal-seeding fields and draws).
struct ParticleICState {
  int iLev = 0;
  int iSpec = 0;
  int iCount = 0;   // particle index within the cell
  int nPerCell = 0; // macroparticles to seed in this cell
  amrex::Real x = 0.0, y = 0.0, z = 0.0;
  amrex::Real q = 0.0; // in/out: macroparticle weight
  amrex::Real uBulk = 0.0, vBulk = 0.0, wBulk = 0.0; // in/out: bulk velocity
  amrex::Real qScale = 1.0; // out: extra weight multiplier (Beam)
  amrex::Real temperature = 0.0;
  amrex::Real charge = 0.0; // species charge (code units)

  // Anisotropic (bi-Maxwellian) seeding: isotropic thermal velocity and the
  // four uniform [0,1) draws may be overwritten for custom sampling.
  amrex::Real uThermal = 0.0, vThermal = 0.0, wThermal = 0.0;
  amrex::Real rand[4] = { 0.0, 0.0, 0.0, 0.0 };
};

//==========================================================
// The parameters of the initial conditions that modify the seeded particles,
// and their per-particle hooks. Each IC class derives from its parameter
// struct, and its virtual hooks call the methods below, so the hooks are
// written once and can also run in the particle-seeding kernels (see
// ICHooks). The methods are defined in src/ic/*.cpp.
//==========================================================

struct BeamICParams {
  int iSpecies_ = -1;                    // beam species
  amrex::Real vel_[nDim3] = { 0, 0, 0 }; // beam bulk velocity
  amrex::Real ratio_ = 0.0;              // beam charge/mass fraction

  AMREX_GPU_HOST_DEVICE void hook_velocity(ParticleICState& s) const;
};

struct WaveICParams {
  enum Profile {
    LightWave,
    HybridWave,
    AlfvenPulse,
    ConvectionWave,
    IonAcousticWave,
    Generic
  };

  Profile profile_ = Generic;

  // Sub-parameters (defaults filled by apply_preset).
  bool seedE_ = false;      // seed the E field
  bool seedB_ = false;      // seed the B field
  bool oblique_ = false;    // oblique plane wave vs x-aligned kx
  bool guideField_ = false; // add uniform guide field Bx0
  bool velKick_ = false;    // matching Alfven ion velocity kick
  bool seedWeight_ = false; // sinusoidal density perturbation
  bool rightHand_ = false;  // transverse (B, u) helicity: T = right-hand

  amrex::Real dir_[nDim3] = { 1, 0, 0 }; // propagation direction (oblique)
  amrex::Real waveLength_ = 48.0;        // wavelength in cells (oblique)
  int waveMode_ = 1;                     // mode number for x-aligned kx
  amrex::Real frac_ = 0.02;              // B perturbation amplitude (B1/Bx0)
  amrex::Real pert_ = 0.1;               // density perturbation amplitude
  amrex::Real anisoTPerpOverTPar_ = 0.0; // T_perp/T_par (0 = isotropic)
  // Gaussian pulse width (0 = disabled, use global sinusoidal mode).
  // When > 0 the sinusoidal seed is multiplied by
  //   exp(-((x - xCenter_) / gaussWidth_)^2)
  // giving a spatially localised Alfven pulse.
  amrex::Real gaussWidth_ = 0.0;
  amrex::Real xCenter_ = 0.0; // Gaussian pulse center (code units)
  // Transverse velocity kick in units of the Alfvenic one (u_perp = -f * B1).
  // 1.0 is the incompressible Alfven relation u_perp = -B_perp/B0; the whistler
  // eigenmode needs -(k d_i)/(omega/Omega_i) instead (see PARAM.XML).
  amrex::Real walenFactor_ = 1.0;

  // Cached in set_fields from the domain / guide field.
  mutable amrex::Real B1_ = 0.0; // B perturbation amplitude (code units)
  mutable amrex::Real kx_ = 0.0; // x wavenumber
  mutable amrex::Real Lx_ = 0.0; // domain length in x

  // +1 for the left-hand (default) seed, -1 for the right-hand one.
  AMREX_GPU_HOST_DEVICE amrex::Real helicity() const {
    return rightHand_ ? -1.0 : 1.0;
  }

  // Effective Gaussian envelope width; <= 0 selects the global sinusoidal mode.
  AMREX_GPU_HOST_DEVICE amrex::Real gaussWidthEffective() const;

  // Transverse envelope (f, g) at x, shared by the seeded B, E and ion kick.
  AMREX_GPU_HOST_DEVICE void envelope(amrex::Real x, amrex::Real gw,
                                      amrex::Real& f, amrex::Real& g) const;

  AMREX_GPU_HOST_DEVICE void hook_weight(ParticleICState& s) const;
  AMREX_GPU_HOST_DEVICE void hook_velocity(ParticleICState& s) const;
  AMREX_GPU_HOST_DEVICE void hook_thermal_velocity(ParticleICState& s) const;
};

struct FadeevICParams {
  amrex::Real L_ = 5.0;         // sheet half-thickness / d_i
  amrex::Real eps_ = 0.4;       // island size vs sheet thickness
  int numIslands_ = 2;          // number of islands along x
  amrex::Real b0_ = 1.0;        // reference field (code units)
  amrex::Real bg_ = 0.0;        // guide-field ratio (B_guide/b0)
  amrex::Real perturb_ = -0.05; // m=1 perturbation amplitude
  amrex::Real nb_ = 0.2;        // background-to-peak density ratio

  // Cached derived quantities.
  mutable amrex::Real Lx_ = 0.0;         // domain length in x
  mutable amrex::Real Ly_ = 0.0;         // domain length in y
  mutable amrex::Real profileMax_ = 0.0; // (1+eps)/(1-eps)
  mutable amrex::Real invLx_ = 0.0;      // 1/Lx
  mutable amrex::Real invLy_ = 0.0;      // 1/Ly

  AMREX_GPU_HOST_DEVICE void hook_weight(ParticleICState& s) const;
  AMREX_GPU_HOST_DEVICE void hook_velocity(ParticleICState& s) const;
};

struct ForceFreeICParams {
  amrex::Real lambda_ = 1.0;          // current sheet thickness / d_i
  amrex::Real b0_ = 1.0;              // asymptotic in-plane field
  amrex::Real bg_ = 0.3;              // guide field
  amrex::Real perturb_ = 0.01;        // initial perturbation amplitude
  bool useUniformIonPressure_ = true; // electrons carry full current by default
  amrex::Real teOverTi_ = -1.0; // optional temperature-weighted current split

  AMREX_GPU_HOST_DEVICE void hook_velocity(ParticleICState& s) const;
};

struct GemICParams {
  // Equilibrium parameters
  amrex::Real b0_ = 1.0;
  amrex::Real b1_ = 1.0;
  amrex::Real b2_ = 1.0;
  amrex::Real tp_ = 1.0;
  amrex::Real t1_ = 1.0;
  amrex::Real t2_ = 1.0;
  amrex::Real lambda0_ = 0.5; // sheet thickness
  amrex::Real apert_ = 0.1;   // perturbation amplitude
  amrex::Real bg_ = 0.0;      // guide field along z (ratio to b0)
  amrex::Real nb_ = 0.2;      // background density ratio
  amrex::Real teOverTi_ = 0.2;

  bool useDoubleCurrentSheet_ = false;
  bool isAsymmetryReconnection_ = false;
  bool useGEMReflected_ = false;
  bool useStandardGem_ = true;
  bool useUniformPressure_ = false;
  bool useUniformIonPressure_ = false;

  amrex::Real gaussX_ = 5.0;
  amrex::Real gaussY_ = 5.0;
  amrex::Real waveLengthX_ = 0.0;
  amrex::Real waveLengthY_ = 0.0;
  // pertType == "gaussian" (GemIC keeps the string).
  bool isGaussianPert_ = false;

  // Cached derived quantities
  mutable bool geomInit_ = false;
  mutable amrex::Real Lx_ = 0.0;
  mutable amrex::Real Ly_ = 0.0;
  mutable amrex::Real Wx_ = 0.0;
  mutable amrex::Real Wy_ = 0.0;
  mutable amrex::Real ySheet_ = 0.0;
  mutable amrex::Real xT_ = 0.0;
  mutable amrex::Real xB_ = 0.0;
  mutable amrex::Real Kx_ = 0.0;
  mutable amrex::Real Ky_ = 0.0;
  mutable amrex::Real gaussXInv_ = 0.0;
  mutable amrex::Real gaussYInv_ = 0.0;
  mutable amrex::Real b0Eff_ = 1.0;
  mutable amrex::Real tpEff_ = 1.0;
  mutable amrex::Real pB_ = 1.0;

  void init_geometry(amrex::Real Lx, amrex::Real Ly) const;
  AMREX_GPU_HOST_DEVICE amrex::Real eval_Bx0(amrex::Real y) const;
  AMREX_GPU_HOST_DEVICE amrex::Real eval_Jz(amrex::Real y) const;
  AMREX_GPU_HOST_DEVICE void eval_perturbation(amrex::Real x, amrex::Real y,
                                               amrex::Real& dbx,
                                               amrex::Real& dby) const;
  AMREX_GPU_HOST_DEVICE amrex::Real eval_Tp(amrex::Real y) const;
  AMREX_GPU_HOST_DEVICE amrex::Real eval_density(amrex::Real y,
                                                 amrex::Real Bx0) const;

  AMREX_GPU_HOST_DEVICE void hook_weight(ParticleICState& s) const;
  AMREX_GPU_HOST_DEVICE void hook_velocity(ParticleICState& s) const;
  AMREX_GPU_HOST_DEVICE void hook_thermal_velocity(ParticleICState& s) const;
};

// The particle hooks of the active initial condition in a form a kernel can
// capture by value: which IC it is, which hooks it uses, and a copy of its
// parameters. See InitialCondition::particle_hooks().
struct ICHooks {
  enum Type { none, beam, wave, fadeev, forceFree, gem };
  int type = none;
  bool modifiesWeights = false;
  bool modifiesVelocities = false;
  bool modifiesThermalVelocity = false;

  BeamICParams beamParams;
  WaveICParams waveParams;
  FadeevICParams fadeevParams;
  ForceFreeICParams forceFreeParams;
  GemICParams gemParams;

  AMREX_GPU_HOST_DEVICE void modify_weight(ParticleICState& s) const {
    switch (type) {
      case wave:
        waveParams.hook_weight(s);
        break;
      case fadeev:
        fadeevParams.hook_weight(s);
        break;
      case gem:
        gemParams.hook_weight(s);
        break;
      default:
        break;
    }
  }

  AMREX_GPU_HOST_DEVICE void modify_velocity(ParticleICState& s) const {
    switch (type) {
      case beam:
        beamParams.hook_velocity(s);
        break;
      case wave:
        waveParams.hook_velocity(s);
        break;
      case fadeev:
        fadeevParams.hook_velocity(s);
        break;
      case forceFree:
        forceFreeParams.hook_velocity(s);
        break;
      case gem:
        gemParams.hook_velocity(s);
        break;
      default:
        break;
    }
  }

  AMREX_GPU_HOST_DEVICE void modify_thermal_velocity(ParticleICState& s) const {
    switch (type) {
      case wave:
        waveParams.hook_thermal_velocity(s);
        break;
      case gem:
        gemParams.hook_thermal_velocity(s);
        break;
      default:
        break;
    }
  }
};

#endif
