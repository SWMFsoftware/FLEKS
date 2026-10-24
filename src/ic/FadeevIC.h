#ifndef _FADEEV_IC_H_
#define _FADEEV_IC_H_

#include "InitialCondition.h"

// Fadeev (force-free) current-sheet reconnection equilibrium.
// See tests/reconnection/README.md for the coordinate mapping and formulas.
class FadeevIC : public InitialCondition, public FadeevICParams {
public:
  std::string name() const override { return "fadeev"; }

  void read_param(ReadParam& param) override;
  void set_fields(PicICFields& fields) const override;

  // Boost the uniform background (#UNIFORMSTATE) to the Fadeev profile.
  bool modifies_weights() const override { return true; }
  void modify_particle_weight(ParticleICState& s) const override {
    hook_weight(s);
  }

  // Seed the diamagnetic (current) drift for force balance.
  bool modifies_velocities() const override { return true; }
  void modify_particle_velocity(ParticleICState& s) const override {
    hook_velocity(s);
  }

  ICHooks particle_hooks() const override {
    ICHooks h;
    h.type = ICHooks::fadeev;
    h.modifiesWeights = true;
    h.modifiesVelocities = true;
    h.modifiesThermalVelocity = false;
    h.fadeevParams = *this;
    return h;
  }
};

#endif
