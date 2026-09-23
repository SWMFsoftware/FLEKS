#ifndef _WAVE_IC_H_
#define _WAVE_IC_H_

#include "InitialCondition.h"

// Wave IC via preset profiles, configurable through #WAVEIC.
//
// Presets (overridable via #WAVEIC sub-parameters):
//   lightwave       : EM oblique circularly polarized plane wave.
//   hybridwave      : transverse B perturbation B1*(cos kx, +-sin kx) on guide
//                     field Bx0 + matching Alfven velocity kick.
//   alfvenpulse     : the same transverse B seed with a Gaussian envelope
//                     (gaussWidth, xCenter) and no velocity kick; the pulse
//                     splits into +-x travelling Alfven packets.
//   ionacousticwave : no field; sinusoidal density perturbation via weight
//                     scaling 1 + pert*sin(kx*x).
//
// The presets are solver-agnostic: nothing keys off useHybridPIC/solveEM.
//
// The +- sign of the B_z (and matching u_z) perturbation is the seed helicity,
// selected with the rightHand sub-parameter:
//   rightHand = F (default, left-hand ): B_z = +B1 sin kx  -> ion-cyclotron
//   rightHand = T (right-hand, whistler): B_z = -B1 sin kx  -> whistler
// relative to a +x guide field B0 = Bx0 with the wave propagating along +B0
// (waveMode > 0).  The right-hand (whistler) sense is the one that rotates
// y -> z about +B0, i.e. the sense of the electron gyration.
class WaveIC : public InitialCondition {
public:
  enum Profile {
    LightWave,
    HybridWave,
    AlfvenPulse,
    ConvectionWave,
    IonAcousticWave,
    Generic
  };

  explicit WaveIC(Profile profile) : profile_(profile) {}

  std::string name() const override;

  void read_param(ReadParam& param) override;
  void set_fields(PicICFields& fields) const override;

  // EM lightwave has no macroparticles; others keep their kinetic ions.
  void apply_particle_override(class ParticlesInfo& pInfo) const override;

  bool modifies_weights() const override { return profile_ == IonAcousticWave; }
  void modify_particle_weight(ParticleICState& s) const override;

  bool modifies_velocities() const override {
    return profile_ == HybridWave || profile_ == AlfvenPulse;
  }
  void modify_particle_velocity(ParticleICState& s) const override;

  // Anisotropic thermal seeding: treat #UNIFORMSTATE T as T_par, inflate the
  // two perpendicular draws by sqrt(T_perp/T_par).
  bool modifies_thermal_velocity() const override {
    return anisoTPerpOverTPar_ > 0.0;
  }
  void modify_particle_thermal_velocity(ParticleICState& s) const override;

private:
  // Apply preset defaults before reading #WAVEIC sub-params, so overrides win.
  void apply_preset();

  Profile profile_;

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
  amrex::Real pert_ = 0.1;              // density perturbation amplitude
  amrex::Real anisoTPerpOverTPar_ = 0.0; // T_perp/T_par (0 = isotropic)
  // Gaussian pulse width (0 = disabled, use global sinusoidal mode).
  // When > 0 the sinusoidal seed is multiplied by
  //   exp(-((x - xCenter_) / gaussWidth_)^2)
  // giving a spatially localised Alfven pulse.
  amrex::Real gaussWidth_ = 0.0;
  amrex::Real xCenter_ = 0.0;           // Gaussian pulse centre (code units)
  // Transverse velocity kick in units of the Alfvenic one (u_perp = -f * B1).
  // 1.0 is the incompressible Alfven relation u_perp = -B_perp/B0; the whistler
  // eigenmode needs -(k d_i)/(omega/Omega_i) instead (see PARAM.XML).
  amrex::Real walenFactor_ = 1.0;

  // +1 for the left-hand (default) seed, -1 for the right-hand one.
  amrex::Real helicity() const { return rightHand_ ? -1.0 : 1.0; }

  // Effective Gaussian envelope width; <= 0 selects the global sinusoidal mode.
  amrex::Real gaussWidthEffective() const;

  // Transverse envelope (f, g) at x, shared by the seeded B, E and ion kick.
  void envelope(amrex::Real x, amrex::Real gw, amrex::Real& f,
                amrex::Real& g) const;

  // Cached in set_fields from the domain / guide field.
  mutable amrex::Real B1_ = 0.0; // B perturbation amplitude (code units)
  mutable amrex::Real kx_ = 0.0; // x wavenumber
  mutable amrex::Real Lx_ = 0.0; // domain length in x
};

#endif
