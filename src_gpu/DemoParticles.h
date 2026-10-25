#ifndef DEMO_PARTICLES_H
#define DEMO_PARTICLES_H

#include <map>
#include <utility>
#include <vector>

#include <AMReX_AmrParticles.H>
#include <AMReX_DenseBins.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>

class DemoConfig;
class DemoGrid;
class DemoFields;

struct DemoParticleDiagnostics {
  amrex::Long count = 0;
  amrex::Real charge = 0;
  amrex::Real kineticEnergy = 0;
  amrex::GpuArray<amrex::Real, 3> current{ 0, 0, 0 };
  amrex::Real idSum = 0;
};

// Production FLEKS AoS layout, including its integer attribute. Runtime
// gathered fields live separately and do not change the particle ABI.
class DemoParticles : public amrex::AmrParticleContainer<4, 1> {
public:
  using Base = amrex::AmrParticleContainer<4, 1>;
  using ParticleType = Base::ParticleType;
  static constexpr int iup_ = 0, ivp_ = 1, iwp_ = 2, iqp_ = 3;

  DemoParticles(const DemoGrid& gridIn, const DemoConfig& configIn,
                amrex::Real chargeIn, amrex::Real massIn, int speciesIn);
  void initialize(amrex::Real perturbation = 0,
                  amrex::GpuArray<amrex::Real, 3> drift = { 0, 0, 0 });
  void drift_positions(amrex::Real dt);
  void gather(const amrex::MultiFab& e, const amrex::MultiFab& b);
  void push(amrex::Real dt, bool fused = false,
            const amrex::MultiFab* e = nullptr,
            const amrex::MultiFab* b = nullptr, bool averageDrift = false);
  // All scatter methods add to existing storage. Caller zeros, sums periodic
  // contributions once, and fills ghosts after ALL species have contributed.
  void deposit(amrex::MultiFab& moments);
  void accumulate_response(DemoFields& fields, amrex::Real dt);
  void accumulate_response_buffered(DemoFields& fields, amrex::Real dt);
  // GPU-optimized: counting-sort particles by cell, reduce each cell's 64
  // node-pair tensors in a thread block, then issue 576 atomics per cell
  // rather than per particle. Falls back to accumulate_response on CPUs.
  void accumulate_response_cell(DemoFields& fields, amrex::Real dt);
  // Precision-specific implementations, selected by cell_precision. Public
  // only because nvcc forbids extended lambdas in non-public members.
  template <class T>
  void accumulate_response_cell_impl(DemoFields& fields, amrex::Real dt);
  template <class T>
  void deposit_center_charge_cell_impl(amrex::MultiFab& rho,
                                       amrex::MultiFab* centerMatrix);
  int benchmark_layouts(const amrex::MultiFab& e, const amrex::MultiFab& b);
  void deposit_center_charge(amrex::MultiFab& rho,
                             amrex::MultiFab* centerMatrix = nullptr);
  // GPU cell-reduced counterparts of deposit_center_charge and deposit, with
  // the same results and caller contract. They fall back to the per-particle
  // versions on CPUs.
  void deposit_center_charge_cell(amrex::MultiFab& rho,
                                  amrex::MultiFab* centerMatrix = nullptr);
  void deposit_cell(amrex::MultiFab& moments);
  void correct_positions(const amrex::MultiFab& phi);
  void sort_by_cell();
  DemoParticleDiagnostics diagnostics() const;
  std::vector<ParticleType> snapshots(bool global = false) const;
  void set_particles(const std::vector<ParticleType>& globalParticles);
  int run_checks(DemoFields& fields);
  amrex::Real species_charge() const { return charge; }
  amrex::Real species_mass() const { return mass; }

private:
  const DemoGrid& grid;
  const DemoConfig& config;
  amrex::Real charge;
  amrex::Real mass;
  int species;
  bool gatheredValid = false;
  std::map<std::pair<int, int>, amrex::Gpu::DeviceVector<amrex::Real> >
      gathered;
  // Persistent per-tile cell bins so the counting sort does not reallocate.
  std::map<std::pair<int, int>, amrex::DenseBins<ParticleType> > cellBins;
  // Optional prefix sums scheduling bounded particle chunks per cell.
  std::map<std::pair<int, int>, amrex::Gpu::DeviceVector<int> > cellTasks;
  // Half-cell-shifted bins for the cell-centered charge deposit.
  std::map<std::pair<int, int>, amrex::DenseBins<ParticleType> > centerBins;
};
#endif
