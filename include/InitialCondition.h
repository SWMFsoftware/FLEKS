#ifndef _INITIAL_CONDITION_H_
#define _INITIAL_CONDITION_H_

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "Constants.h"

class Pic;        // forward: facade uses public getters only, never the full type
class ReadParam;  // forward: read_param takes a reference
class ParticlesInfo;  // forward: apply_particle_override takes a reference

// Narrow per-particle state handed to a particle-seeding IC. The IC may read the
// particle's position / index and modify its bulk velocity and macroparticle
// weight. It never touches Particles internals directly.
struct ParticleICState {
  int iLev = 0;
  int iSpec = 0;
  int iCount = 0;       // particle index within the cell
  int nPerCell = 0;     // total macroparticles to seed in this cell
  amrex::Real x = 0.0, y = 0.0, z = 0.0;
  amrex::Real q = 0.0;             // in/out: macroparticle weight q
  amrex::Real uBulk = 0.0, vBulk = 0.0, wBulk = 0.0;  // in/out: bulk velocity
  amrex::Real qScale = 1.0;        // out: extra weight multiplier (Beam)
  amrex::Real temperature = 0.0;

  // --- Anisotropic (bi-Maxwellian) thermal seeding hook. ---
  // When the IC declares modifies_thermal_velocity(), the isotropic thermal
  // velocity sampled from #UNIFORMSTATE is handed back through these fields
  // together with the four uniform [0,1) draws used to sample it.  The IC may
  // overwrite uThermal/vThermal/wThermal with its own anisotropic sampling
  // (e.g. distinct T_perp vs T_par along the guide field), and the overwritten
  // values are used in place of the isotropic one before the bulk velocity is
  // added.  rand[0..3] are the raw draws so the IC can run its own Box-Muller.
  amrex::Real uThermal = 0.0, vThermal = 0.0, wThermal = 0.0;
  amrex::Real rand[4] = {0.0, 0.0, 0.0, 0.0};
};

// Registry of named InitialCondition factories. Lookup is case-insensitive.
class ICRegistry {
public:
  using Factory = std::function<std::unique_ptr<class InitialCondition>()>;

  static ICRegistry& instance();

  void register_ic(const std::string& name, Factory f);
  std::unique_ptr<class InitialCondition> create(const std::string& name) const;
  std::vector<std::string> names() const;

private:
  ICRegistry() = default;
  std::map<std::string, Factory> registry_;
};

// Built-in InitialCondition registration. Defined in src/ic/RegisterAll.cpp and
// called once (lazily) before the first registry lookup.
void register_all_initial_conditions();

// Abstract base class for every test-case initial condition.
class InitialCondition {
public:
  virtual ~InitialCondition() = default;

  // Read IC-specific sub-parameters (called right after the #TESTCASE name is
  // parsed, with the ReadParam cursor positioned at the first sub-parameter).
  virtual void read_param(ReadParam& param) { (void)param; }

  // Seed the EM fields (nodeE / nodeB / centerB) through the narrow facade.
  virtual void set_fields(class PicICFields& fields) const { (void)fields; }

  // Apply a particle-count override to the running ParticlesInfo. Called once
  // before particle construction; used by LightWave / TopHat to force zero
  // particles. Default: no-op.
  virtual void apply_particle_override(class ParticlesInfo& pInfo) const {
    (void)pInfo;
  }

  // ---- Per-particle hooks used inside Particles::add_particles_cell. ----
  // Weight hook runs BEFORE the q != 0 guard (so an IC can zero a particle's
  // weight to suppress it). Velocity hook runs AFTER the bulk velocity is
  // computed (so an IC can override / perturb the bulk velocity).
  virtual bool modifies_weights() const { return false; }
  virtual void modify_particle_weight(ParticleICState& s) const { (void)s; }
  virtual bool modifies_velocities() const { return false; }
  virtual void modify_particle_velocity(ParticleICState& s) const { (void)s; }

  // ---- Anisotropic (bi-Maxwellian) thermal-seeding hook. ----
  // Runs right after the isotropic thermal velocity is sampled from
  // #UNIFORMSTATE but BEFORE the bulk velocity is added.  An IC that declares
  // modifies_thermal_velocity() may overwrite s.uThermal/vThermal/wThermal with
  // a distribution with distinct perpendicular / parallel temperature (e.g.
  // T_perp != T_par for the proton-cyclotron anisotropy instability), using the
  // four uniform draws s.rand[0..3] for its own Box-Muller sampling.  The
  // overwritten thermal velocity is added to the bulk velocity downstream.
  virtual bool modifies_thermal_velocity() const { return false; }
  virtual void modify_particle_thermal_velocity(ParticleICState& s) const {
    (void)s;
  }

  // Field-step seeding hook for tests that seed a discontinuous field in a
  // sub-region of new (refined) cells. Kept as a genuine IC hook (it is a
  // field-seeding behaviour, not a solver override); the old solver overrides
  // uMax / upwind-velocity were promoted to #FIXEDUMAX / #UPWINDB fixedUpwindVel
  // options and the override_umax()/bypass_limiter() hooks removed.
  virtual bool is_tophat() const { return false; }

  virtual std::string name() const = 0;
};

// Narrow facade exposing only what an InitialCondition needs from the running
// Pic. Deliberately NOT a friend of Pic: it goes through public getters.
class PicICFields {
public:
  explicit PicICFields(Pic& pic) : pic_(&pic) {}

  int n_lev() const;
  const amrex::Geometry& geom(int iLev) const;

  amrex::MultiFab& node_E(int iLev);
  amrex::MultiFab& node_B(int iLev);
  amrex::MultiFab& center_B(int iLev);

  // Fill ghost cells of all three field arrays (used by field ICs).
  void fill_boundary_E_B();

private:
  Pic* pic_;
};

#endif
