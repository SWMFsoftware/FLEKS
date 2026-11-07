#ifndef _BEAM_IC_H_
#define _BEAM_IC_H_

#include "InitialCondition.h"

// First nBeam = npcel*ratio macroparticles per cell get the bulk velocity `vel`
// (weight scaled up); the rest are the background (weight scaled down),
// preserving total count/mass with `ratio` of the charge.
class BeamIC : public InitialCondition, public BeamICParams {
public:
  std::string name() const override { return "beam"; }

  void read_param(ReadParam& param) override;
  bool modifies_velocities() const override { return true; }
  void modify_particle_velocity(ParticleICState& s) const override {
    hook_velocity(s);
  }
  ICHooks particle_hooks() const override {
    ICHooks h;
    h.type = ICHooks::beam;
    h.modifiesWeights = false;
    h.modifiesVelocities = true;
    h.modifiesThermalVelocity = false;
    h.beamParams = *this;
    return h;
  }
};

#endif
