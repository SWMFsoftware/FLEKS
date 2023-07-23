#include <fstream>
#include <iostream>
#include <vector>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "Converter.h"
#include "GridUtility.h"

using namespace amrex;
int main(int argc, char* argv[]) {
  Initialize(MPI_COMM_WORLD);

  std::vector<std::string> commandLine;
  for (int i = 0; i < argc; ++i) {
    commandLine.push_back((std::string)(argv[i]));
  }
  commandLine.push_back("3d_fluid_region0_2_t00020036_n00000010_amrex");

  std::array<std::string, 2> arg = { "-h", "-help" };
  if (argc > 1 && find(arg.begin(), arg.end(), commandLine[1]) != arg.end()) {
    std::cout << " \n"
              << " This exectuable converts a *_amrex file to a tecplot ascii "
                 "*.dat file. "
                 " \n\n Usage:\n "
                 " ./Converter.exe *_amrex\n\n";
    return 0;
  }

  for (std::vector<std::string>::size_type i = 1; i < commandLine.size(); i++) {
    std::cout << commandLine[i] << std::endl;
    Converter cv(commandLine[i]);
    cv.read();
    cv.write();
  }

  Finalize();

  return 1;
}
