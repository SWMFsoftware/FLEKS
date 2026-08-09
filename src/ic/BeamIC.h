#ifndef _BEAM_IC_H_
#define _BEAM_IC_H_

#include "InitialCondition.h"

// Charged-particle beam. Within each cell the first nBeam = npcel*ratio
// macroparticles become the beam (bulk velocity = vel, weight scaled up),
// the remaining npcel - nBeam become the colder background (weight scaled down)
// so that the total particle count and total mass are preserved while the beam
// carries a fraction `ratio` of the charge / mass.
class BeamIC : public InitialCondition {
public:
  std::string name() const override { return "beam"; }

  void read_param(ReadParam& param) override;
  bool modifies_velocities() const override { return true; }
  void modify_particle_velocity(ParticleICState& s) const override;

private:
  int iSpecies_ = -1;                    // species that carries the beam
  amrex::Real vel_[nDim3] = { 0, 0, 0 }; // beam bulk velocity
  amrex::Real ratio_ = 0.0;              // beam charge/mass fraction
};

#endif
