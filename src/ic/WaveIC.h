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
class WaveIC : public InitialCondition, public WaveICParams {
public:
  explicit WaveIC(Profile profile) { profile_ = profile; }

  std::string name() const override;

  void read_param(ReadParam& param) override;
  void set_fields(PicICFields& fields) const override;

  // EM lightwave has no macroparticles; others keep their kinetic ions.
  void apply_particle_override(class ParticlesInfo& pInfo) const override;

  bool modifies_weights() const override { return profile_ == IonAcousticWave; }
  void modify_particle_weight(ParticleICState& s) const override {
    hook_weight(s);
  }

  bool modifies_velocities() const override {
    return profile_ == HybridWave || profile_ == AlfvenPulse;
  }
  void modify_particle_velocity(ParticleICState& s) const override {
    hook_velocity(s);
  }

  // Anisotropic thermal seeding: treat #UNIFORMSTATE T as T_par, inflate the
  // two perpendicular draws by sqrt(T_perp/T_par).
  bool modifies_thermal_velocity() const override {
    return anisoTPerpOverTPar_ > 0.0;
  }
  void modify_particle_thermal_velocity(ParticleICState& s) const override {
    hook_thermal_velocity(s);
  }

  ICHooks particle_hooks() const override {
    ICHooks h;
    h.type = ICHooks::wave;
    h.modifiesWeights = modifies_weights();
    h.modifiesVelocities = modifies_velocities();
    h.modifiesThermalVelocity = modifies_thermal_velocity();
    h.waveParams = *this;
    return h;
  }

private:
  // Apply preset defaults before reading #WAVEIC sub-params, so overrides win.
  void apply_preset();
};

#endif
