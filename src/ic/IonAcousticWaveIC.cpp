#include "IonAcousticWaveIC.h"

#include <AMReX.H>

#include "ReadParam.h"

void IonAcousticWaveIC::read_param(ReadParam& param) {
  param.read_var("pert", pert_);
  param.read_var("waveMode", waveMode_);
}

// No EM field seeding; just cache the domain length so the weight hook can
// build kx without needing geom access during add_particles_cell.
void IonAcousticWaveIC::set_fields(PicICFields& fields) const {
  if (fields.n_lev() > 0) {
    Lx_ = (fields.geom(0).ProbHi())[0] - (fields.geom(0).ProbLo())[0];
  }
}

// Exactly reproduces the add_particles_cell weight scaling:
//   factor = 1 + pert * sin(kx * x),   kx = waveMode * 2*pi / Lx
//   q *= (factor > 0) ? factor : 0  -> a non-positive factor zeroes the
//   particle weight so it is skipped by the q != 0 guard.
void IonAcousticWaveIC::modify_particle_weight(ParticleICState& s) const {
  if (Lx_ <= 0.0)
    return;
  const amrex::Real kx = waveMode_ * (2.0 * dPI) / Lx_;
  const amrex::Real factor = 1.0 + pert_ * std::sin(kx * s.x);
  s.q *= (factor > 0.0) ? factor : 0.0;
}
