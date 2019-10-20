#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include "Domain.h"

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv);

  amrex::IntVect lo(AMREX_D_DECL(64, 64, 64));

  amrex::AllPrint() << "Hello world from AMReX version " << amrex::Version()
                    << " MyProc = " << amrex::ParallelDescriptor::MyProc()
                    << "\n";

  amrex::Print(1) << "lo = " << lo;

  amrex::Finalize();
}
