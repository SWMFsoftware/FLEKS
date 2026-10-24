#include <AMReX_Gpu.H>
#include <AMReX_Print.H>

#include "DemoGrid.h"

DemoGrid::DemoGrid(const DemoConfig& config) {
  static_assert(AMREX_SPACEDIM == 3, "The GPU porting demo requires 3D AMReX");
  const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(config.nCell - 1));
  const amrex::RealBox physical({ 0., 0., 0. }, { 1., 1., 1. });
  const amrex::Array<int, 3> periodic{ 1, 1, 1 };
  geom.define(domain, physical, amrex::CoordSys::cartesian, periodic);
  ba.define(domain);
  ba.maxSize(config.maxGridSize);
  dm.define(ba);
  nodeBa = amrex::convert(ba, amrex::IntVect::TheNodeVector());
  // The 243-real response dominates memory. Include persistent Krylov vectors,
  // field scratch, two species, and benchmark copies before allocating them.
  double bytes = 0;
  for (int box = 0; box < ba.size(); ++box) {
    if (dm[box] != amrex::ParallelDescriptor::MyProc())
      continue;
    const double nodes = amrex::grow(nodeBa[box], 1).numPts();
    const double cells = amrex::grow(ba[box], 1).numPts();
    bytes += sizeof(amrex::Real) *
             (nodes * (243 + 3 * (config.gmresRestart + 1) + 80) +
              cells * (config.gmresRestart + 80));
    bytes += ba[box].numPts() * double(config.particlesPerCell) * 2 *
             (20 * sizeof(amrex::Real) + 32);
    if (config.mode == "benchmark")
      bytes += nodes * 243 * sizeof(amrex::Real);
  }
  amrex::ParallelDescriptor::ReduceRealMax(bytes);
  amrex::Print() << "MEMORY estimated_max_rank_bytes=" << bytes << "\n";
#ifdef AMREX_USE_GPU
  if (bytes > 0.7 * amrex::Gpu::Device::totalGlobalMem())
    amrex::Abort("Demo memory estimate exceeds 70% of device memory; reduce "
                 "n_cell, particles_per_cell, or gmres_restart");
#endif
}
