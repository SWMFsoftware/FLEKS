#include "show_git_info.h"
#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

int main(int argc, char* argv[]) {
  using namespace amrex;

  Initialize(argc, argv);
  {
    if (ParallelDescriptor::MyProc() == 0)
      print_git_info();

    // Domain domain;

    // domain.init();
  }

  Finalize();
}
