#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include "Domain.h"
#include "show_git_info.h"

int main(int argc, char* argv[]) {

  amrex::Initialize(argc, argv);
  {
    if (amrex::ParallelDescriptor::MyProc() == 0)
      print_git_info();

    // Domain domain;

    // domain.init();
  }

  amrex::Finalize();
}
