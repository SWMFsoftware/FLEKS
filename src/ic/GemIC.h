#ifndef _GEM_IC_H_
#define _GEM_IC_H_

#include <string>

#include "InitialCondition.h"

/**
 * @brief GEM Challenge and Asymmetric Magnetic Reconnection Initial Condition.
 *
 * Implements magnetic reconnection equilibria and perturbation models for:
 *   1. Classic standard GEM challenge (Harris sheet, conducting walls or
 *      periodic, sin/cos or Gaussian perturbation).
 *   2. Symmetric double current sheet (Harris equilibrium with sheets at
 *      +-0.25*WaveLengthY).
 *   3. Asymmetric double current sheet (asymmetric magnetic field B1 != B2 and
 *      temperatures T1 != T2 across periodic quadruple tanh current sheets).
 *   4. Reflected GEM current sheet.
 */
class GemIC : public InitialCondition, public GemICParams {
public:
  std::string name() const override { return "gem"; }

  void read_param(ReadParam& param) override;
  void set_fields(PicICFields& fields) const override;

  bool modifies_weights() const override { return true; }
  void modify_particle_weight(ParticleICState& s) const override {
    hook_weight(s);
  }

  bool modifies_velocities() const override { return true; }
  void modify_particle_velocity(ParticleICState& s) const override {
    hook_velocity(s);
  }

  bool modifies_thermal_velocity() const override {
    return isAsymmetryReconnection_;
  }
  void modify_particle_thermal_velocity(ParticleICState& s) const override {
    hook_thermal_velocity(s);
  }

  ICHooks particle_hooks() const override {
    ICHooks h;
    h.type = ICHooks::gem;
    h.modifiesWeights = true;
    h.modifiesVelocities = true;
    h.modifiesThermalVelocity = modifies_thermal_velocity();
    h.gemParams = *this;
    return h;
  }

private:
  std::string pertType_ = "";
};

#endif
