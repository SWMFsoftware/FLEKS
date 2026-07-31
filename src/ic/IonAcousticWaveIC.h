#ifndef _ION_ACOUSTIC_WAVE_IC_H_
#define _ION_ACOUSTIC_WAVE_IC_H_

#include "InitialCondition.h"

// Ion-acoustic-wave seed. Seeds a sinusoidal ion-density perturbation
//     n(x) = n0 * (1 + pert * sin(kx * x)),   kx = 2*pi*waveMode / Lx
// by scaling each macroparticle's weight q. Because all moments scale with q,
// the bulk velocity stays uniform while the density / charge / pressure moments
// carry the perturbation. E = B = 0 otherwise.
class IonAcousticWaveIC : public InitialCondition {
public:
  std::string name() const override { return "ionacousticwave"; }

  void read_param(ReadParam& param) override;
  void set_fields(PicICFields& fields) const override;
  bool modifies_weights() const override { return pert_ != 0.0; }
  void modify_particle_weight(ParticleICState& s) const override;

private:
  amrex::Real pert_ = 0.0;   // density perturbation amplitude (< 1)
  int waveMode_ = 1;         // number of wavelengths across the x-domain
  mutable amrex::Real Lx_ = 0.0;  // cached domain length (from set_fields)
};

#endif
