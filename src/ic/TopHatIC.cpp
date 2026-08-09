#include "Particles.h"
#include "TopHatIC.h"

// The TopHat test is a pure EM step test: no macroparticles are seeded.
void TopHatIC::apply_particle_override(ParticlesInfo& pInfo) const {
  pInfo.nPartPerCell = amrex::IntVect::Zero;
}
