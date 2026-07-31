#ifndef _HYBRID_WAVE_IC_H_
#define _HYBRID_WAVE_IC_H_

#include "InitialCondition.h"

// Parallel-propagating circularly-polarized Alfven/whistler wave seed.
//
// The field part (a (B1 cos kx, B1 sin kx) transverse perturbation on the
// uniform guide field Bx0) is written in set_fields(); the matching ion
// velocity kick (-B1 cos kx, -B1 sin kx) is applied per particle in
// modify_particle_velocity() so that the initial state is a clean normal mode.
//
// Two flavours are registered from one class:
//   "hybridwave"     -> perturbVel = true,  waveFrac = 0.02
//   "convectionwave" -> perturbVel = false, waveFrac = 0.20
class HybridWaveIC : public InitialCondition {
public:
  HybridWaveIC(bool perturbVel, amrex::Real waveFrac)
    : perturbVel_(perturbVel), waveFrac_(waveFrac) {}

  std::string name() const override {
    return perturbVel_ ? "hybridwave" : "convectionwave";
  }

  void set_fields(PicICFields& fields) const override;
  bool modifies_velocities() const override { return perturbVel_; }
  void modify_particle_velocity(ParticleICState& s) const override;

private:
  bool perturbVel_;
  amrex::Real waveFrac_;

  // Cached from set_fields(): per-level perturbation amplitude and wavenumber.
  // Mutable so the const set_fields() can populate them for the (also const)
  // modify_particle_velocity() called later during fill_particles().
  mutable amrex::Vector<amrex::Real> B1_;
  mutable amrex::Vector<amrex::Real> kx_;
};

#endif
