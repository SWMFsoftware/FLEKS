#ifndef _WAVE_IC_H_
#define _WAVE_IC_H_

#include "InitialCondition.h"

// Unified transverse / sinusoidal wave initial condition.
//
// The four former wave plugins (LightWaveIC, HybridWaveIC, IonAcousticWaveIC)
// were folded into this single class. Every wave test now configures through
// ONE #WAVEIC parameter block, so adding a new wave test costs zero C++:
// register nothing, just write a #TESTCASE waveic (or one of the legacy
// aliases) plus a #WAVEIC block in the test's PARAM.in.
//
// The four legacy names are still registered (lightwave, hybridwave,
// convectionwave, ionacousticwave); each maps to a preset that reproduces the
// exact field / particle seeding of the original plugin, so existing PARAM.in
// files keep working unchanged. All behaviour is also overridable from a
// #WAVEIC block, which is how a brand-new wave test is expressed.
//
// Presets (overridable via the matching #WAVEIC sub-parameter):
//   lightwave       : E+B oblique circularly polarized plane wave,
//                     direction (0.6, 0.8, 0), wavelength 48 (in cells).
//   hybridwave      : transverse B perturbation B1*(cos kx, sin kx) on top of
//                     the guide field Bx0, kx = 2*pi/Lx, B1 = 0.02*Bx0, with a
//                     matching Alfven ion velocity kick.
//   convectionwave   : same B perturbation as hybridwave but B1 = 0.2*Bx0 and
//   NO
//                     velocity kick (advected by a uniform bulk flow instead).
//   ionacousticwave : NO field; sinusoidal density perturbation via weight
//                     scaling 1 + pert*sin(kx*x), kx = waveMode*2*pi/Lx.
class WaveIC : public InitialCondition {
public:
  enum Profile {
    LightWave,
    HybridWave,
    ConvectionWave,
    IonAcousticWave,
    Generic
  };

  explicit WaveIC(Profile profile) : profile_(profile) {}

  std::string name() const override;

  void read_param(ReadParam& param) override;
  void set_fields(PicICFields& fields) const override;

  // The pure-EM lightwave test has no macroparticles: zero the per-cell count
  // so the electron background remains a fluid. Wave tests with kinetic ions
  // (hybridwave / convectionwave / ionacousticwave) keep their particles.
  void apply_particle_override(class ParticlesInfo& pInfo) const override;

  bool modifies_weights() const override { return profile_ == IonAcousticWave; }
  void modify_particle_weight(ParticleICState& s) const override;

  bool modifies_velocities() const override { return profile_ == HybridWave; }
  void modify_particle_velocity(ParticleICState& s) const override;

  // Anisotropic (bi-Maxwellian) thermal seeding: when anisoTPerpOverTPar_ > 0
  // the isotropic thermal velocity sampled from #UNIFORMSTATE is rescaled so
  // the component along the guide field keeps the #UNIFORMSTATE temperature
  // (interpreted as T_par) and the two perpendicular components are inflated by
  // sqrt(T_perp/T_par).  This drives the proton-cyclotron anisotropy
  // instability from an anisotropic free-energy source.
  bool modifies_thermal_velocity() const override {
    return anisoTPerpOverTPar_ > 0.0;
  }
  void modify_particle_thermal_velocity(ParticleICState& s) const override;

private:
  // Apply the preset defaults for the active profile *before* any #WAVEIC
  // sub-parameters are read, so user overrides win.
  void apply_preset();

  Profile profile_;

  // --- #WAVEIC sub-parameters (defaults filled by apply_preset). ---
  bool seedE_ = false;      // seed the E field (lightwave)
  bool seedB_ = false;      // seed the B field (hybrid/convection/lightwave)
  bool oblique_ = false;    // oblique plane wave (lightwave) vs x-aligned kx
  bool guideField_ = false; // add the uniform guide field Bx0 under the pert
  bool velKick_ = false;    // matching Alfven ion velocity kick (hybridwave)
  bool seedWeight_ = false; // sinusoidal density perturbation (ionacousticwave)

  amrex::Real dir_[nDim3] = { 1, 0, 0 }; // propagation direction (oblique case)
  amrex::Real waveLength_ = 48.0;        // wavelength in cells (oblique case)
  int waveMode_ = 1;                     // mode number for x-aligned kx
  amrex::Real frac_ = 0.02;              // B perturbation amplitude (B1/Bx0)
  amrex::Real pert_ = 0.1;               // density perturbation amplitude (IAW)
  amrex::Real anisoTPerpOverTPar_ = 0.0; // T_perp/T_par ratio (0 = isotropic)

  // --- cached in set_fields from the domain / guide field. ---
  mutable amrex::Real B1_ = 0.0; // B perturbation amplitude (code units)
  mutable amrex::Real kx_ = 0.0; // x wavenumber
  mutable amrex::Real Lx_ = 0.0; // domain length in x
};

#endif
