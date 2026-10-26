#ifndef _FORCE_FREE_IC_H_
#define _FORCE_FREE_IC_H_

#include "InitialCondition.h"

// Force-free current-sheet reconnection equilibrium (Le et al. 2016).
// In 2D (x, y) with current sheet along x and variation along y:
//   Bx(x, y) = b0 * tanh(y / lambda) - perturb * b0 * (Lx / (2*Ly)) *
//   cos(2*pi*x / Lx) * sin(pi*y / Ly) By(x, y) = perturb * b0 * sin(2*pi*x /
//   Lx) * cos(pi*y / Ly) Bz(x, y) = sqrt(bg^2 + b0^2 * sech^2(y / lambda))
// Total |B|^2 = b0^2 + bg^2 = const everywhere, so total magnetic pressure is
// constant and plasma density is strictly uniform.
class ForceFreeIC : public InitialCondition, public ForceFreeICParams {
public:
  std::string name() const override { return "forcefree"; }

  void read_param(ReadParam& param) override;
  void set_fields(PicICFields& fields) const override;

  bool modifies_velocities() const override { return true; }
  void modify_particle_velocity(ParticleICState& s) const override {
    hook_velocity(s);
  }

  ICHooks particle_hooks() const override {
    ICHooks h;
    h.type = ICHooks::forceFree;
    h.modifiesWeights = false;
    h.modifiesVelocities = true;
    h.modifiesThermalVelocity = false;
    h.forceFreeParams = *this;
    return h;
  }
};

#endif
