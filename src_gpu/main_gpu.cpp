// Standalone FLEKS-style GPU PIC demo. See README.md and
// IMPLEMENTATION_PLAN.md.
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "DemoConfig.h"
#include "DemoPic.h"

static void set_default_parameters() {
  amrex::ParmParse pp("amrex");
  if (!pp.contains("use_gpu_aware_mpi"))
    pp.add("use_gpu_aware_mpi", 0);
  if (!pp.contains("the_arena_init_size"))
    pp.add("the_arena_init_size", 64 * 1024 * 1024);
}

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, set_default_parameters);
  int failures = 0;
  {
    const DemoConfig config;
    DemoPic simulation(config);
    failures = simulation.run();
  }
  amrex::ParallelDescriptor::ReduceIntMax(failures);
  if (failures == 0)
    amrex::Print() << "ALL GPU PIC DEMO CHECKS PASSED\n";
  else
    amrex::Print() << failures << " GPU PIC DEMO CHECK(S) FAILED\n";
  amrex::Finalize();
  return failures == 0 ? 0 : 1;
}
