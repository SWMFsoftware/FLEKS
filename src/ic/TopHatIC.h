#ifndef _TOP_HAT_IC_H_
#define _TOP_HAT_IC_H_

#include "InitialCondition.h"

// TopHat is a pure EM step test: it seeds a discontinuous field (a unit Ey /
// Bz step between x in [xL, xR]) in new (refined) cells, and forces zero
// macroparticles. The step seeding is a genuine IC behaviour expressed through
// the is_tophat() hook (Pic::fill_new_node_E/B performs the actual writing).
//
// The two solver-behaviour overrides it previously needed -- a fixed CFL signal
// speed (uMax = 1.0) and a fixed upwind-B velocity (1.0) -- were promoted to
// ordinary PARAM options in Phase 4: set #FIXEDUMAX 1.0 and #UPWINDB with
// fixedUpwindVel 1.0 in the test's PARAM.in. The TopHat plugin no longer needs
// to touch solver behaviour.
class TopHatIC : public InitialCondition {
public:
  std::string name() const override { return "tophat"; }

  void apply_particle_override(class ParticlesInfo& pInfo) const override;
  bool is_tophat() const override { return true; }
};

#endif
