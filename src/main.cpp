#include "show_git_info.h"
#include <AMReX.H>
#include <AMReX_Print.H>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "Domain.h"
#include "SimDomains.h"

Domains fleksDomains;

extern "C" {
void timing_start_c(size_t* nameLen, char* name) {}
void timing_stop_c(size_t* nameLen, char* name) {}
}

int main(int argc, char* argv[]) {
  using namespace amrex;

  Initialize(argc, argv);
  {
    if (ParallelDescriptor::MyProc() == 0)
      print_git_info();

    // 1. Read PARAM.in into a string
    std::string paramString;
    if (ParallelDescriptor::MyProc() == 0) {
      std::ifstream infile("PARAM.in");
      if (infile.is_open()) {
        std::string line;
        while (std::getline(infile, line)) {
          paramString += line + "\n";
        }
        infile.close();
      } else {
        std::cerr << "Error: Could not open PARAM.in" << std::endl;
        ParallelDescriptor::Abort();
      }
    }

    // Broadcast paramString to all processors
    int paramLen = paramString.length();
    ParallelDescriptor::Bcast(&paramLen, 1, 0);
    if (ParallelDescriptor::MyProc() != 0) {
      paramString.resize(paramLen);
    }
    ParallelDescriptor::Bcast(&paramString[0], paramLen, 0);

    // Create output directories
    if (ParallelDescriptor::IOProcessor()) {
      std::filesystem::create_directories("FLEKS1/plots");
      std::filesystem::create_directories("FLEKS1/restartOUT");
    }

    // 2. Initialize Domain
    fleksDomains.add_new_domain();
    fleksDomains.select(0);
    Domain& domain = fleksDomains(0);
    domain.init(0.0, 1, paramString);

    // Turn on all cells.
    domain.receive_grid_info();

    // Create grids for all components.
    domain.regrid();

    // 3. Set Initial Conditions
    domain.set_ic();

    // 4. Run Loop
    double timeMax = 0.0;
    int maxIter = -1;

    // Extract stop criteria from PARAM.in
    {
      ReadParam reader;
      reader = paramString;
      std::string command;
      while (reader.get_next_command(command)) {
        if (command == "#STOP") {
          reader.read_var("MaxIter", maxIter);
          reader.read_var("TimeMax", timeMax);
        }
      }
    }

    while (true) {
      if (maxIter >= 0 && domain.tc->get_cycle() >= maxIter)
        break;
      if (timeMax > 0 && domain.tc->get_time_si() >= timeMax - 1e-10)
        break;

      domain.update();
    }

    amrex::Print() << "\nSimulation finished at time = "
                   << domain.tc->get_time_si() << std::endl;

    // 5. Final output
    domain.write_plots(true);
  }

  Finalize();
  return 0;
}
