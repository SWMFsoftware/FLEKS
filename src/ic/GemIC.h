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
class GemIC : public InitialCondition {
public:
  std::string name() const override { return "gem"; }

  void read_param(ReadParam& param) override;
  void set_fields(PicICFields& fields) const override;

  bool modifies_weights() const override { return true; }
  void modify_particle_weight(ParticleICState& s) const override;

  bool modifies_velocities() const override { return true; }
  void modify_particle_velocity(ParticleICState& s) const override;

  bool modifies_thermal_velocity() const override {
    return isAsymmetryReconnection_;
  }
  void modify_particle_thermal_velocity(ParticleICState& s) const override;

private:
  void init_geometry(amrex::Real Lx, amrex::Real Ly) const;
  amrex::Real eval_Bx0(amrex::Real y) const;
  amrex::Real eval_Jz(amrex::Real y) const;
  void eval_perturbation(amrex::Real x, amrex::Real y, amrex::Real& dbx,
                         amrex::Real& dby) const;
  amrex::Real eval_Tp(amrex::Real y) const;
  amrex::Real eval_density(amrex::Real y, amrex::Real Bx0) const;

  // Equilibrium parameters
  amrex::Real b0_ = 1.0;
  amrex::Real b1_ = 1.0;
  amrex::Real b2_ = 1.0;
  amrex::Real tp_ = 1.0;
  amrex::Real t1_ = 1.0;
  amrex::Real t2_ = 1.0;
  amrex::Real lambda0_ = 0.5; // sheet thickness
  amrex::Real apert_ = 0.1;   // perturbation amplitude
  amrex::Real bg_ = 0.0;      // guide field along z (ratio to b0)
  amrex::Real nb_ = 0.2;      // background density ratio
  amrex::Real teOverTi_ = 0.2;

  bool useDoubleCurrentSheet_ = false;
  bool isAsymmetryReconnection_ = false;
  bool useGEMReflected_ = false;
  bool useStandardGem_ = true;
  bool useUniformPressure_ = false;
  bool useUniformIonPressure_ = false;

  amrex::Real gaussX_ = 5.0;
  amrex::Real gaussY_ = 5.0;
  amrex::Real waveLengthX_ = 0.0;
  amrex::Real waveLengthY_ = 0.0;
  std::string pertType_ = "";

  // Cached derived quantities
  mutable bool geomInit_ = false;
  mutable amrex::Real Lx_ = 0.0;
  mutable amrex::Real Ly_ = 0.0;
  mutable amrex::Real Wx_ = 0.0;
  mutable amrex::Real Wy_ = 0.0;
  mutable amrex::Real ySheet_ = 0.0;
  mutable amrex::Real xT_ = 0.0;
  mutable amrex::Real xB_ = 0.0;
  mutable amrex::Real Kx_ = 0.0;
  mutable amrex::Real Ky_ = 0.0;
  mutable amrex::Real gaussXInv_ = 0.0;
  mutable amrex::Real gaussYInv_ = 0.0;
  mutable amrex::Real b0Eff_ = 1.0;
  mutable amrex::Real tpEff_ = 1.0;
  mutable amrex::Real pB_ = 1.0;
};

#endif
