#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "Domain.h"
#include "SimDomains.h"
#include "show_git_info.h"

// Normal completion clears this registry before Finalize(). Keep the empty
// registry alive on abnormal process termination to avoid static destruction
// after AMReX has finalized its allocators.
Domains& fleksDomains = *new Domains;

extern "C" {
void timing_start_c(size_t* nameLen, char* name) {}
void timing_stop_c(size_t* nameLen, char* name) {}

#ifdef _PT_COMPONENT_
void OH_get_charge_exchange_wrapper(double* /*rhoIon*/, double* /*cs2Ion*/,
                                    double /*uIon_D*/[3], double* /*rhoNeu*/,
                                    double* /*cs2Neu*/, double /*uNeu_D*/[3],
                                    double /*sourceIon_V*/[5],
                                    double /*sourceNeu_V*/[5]) {}

void OH_get_charge_exchange_region(int* iRegion, double* /*r*/,
                                   double* /*rhoDim*/, double* /*u2Dim*/,
                                   double* /*uSW2Dim*/, double* /*tempDim*/,
                                   double* /*tempPu2Dim*/, double* /*mach2*/,
                                   double* /*machPUI2*/, double* /*machSW2*/,
                                   double* /*levHP*/) {
  if (iRegion)
    *iRegion = -1;
}

void OH_get_solar_wind(double* /*x*/, double* /*y*/, double* /*z*/,
                       double* /*numDen*/, double* /*ur*/, double* /*temp*/,
                       double /*b*/[3]) {}
#endif
}

namespace {

std::string prepare_standalone_run() {
  std::string paramString;

  if (amrex::ParallelDescriptor::MyProc() == 0) {
    std::ifstream infile("PARAM.in");
    if (!infile.is_open()) {
      std::cerr << "Error: Could not open PARAM.in" << std::endl;
      amrex::ParallelDescriptor::Abort();
    }

    std::string line;
    while (std::getline(infile, line)) {
      paramString += line + "\n";
    }
  }

  int paramLen = static_cast<int>(paramString.length());
  amrex::ParallelDescriptor::Bcast(&paramLen, 1, 0);
  if (amrex::ParallelDescriptor::MyProc() != 0) {
    paramString.resize(paramLen);
  }
  if (paramLen > 0) {
    amrex::ParallelDescriptor::Bcast(paramString.data(), paramLen, 0);
  }

  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::filesystem::create_directories("PC/plots");
    std::filesystem::create_directories("PC/restartOUT");
  }

  return paramString;
}

std::vector<std::string> split_sessions(const std::string& paramString) {
  // #RUN ends one parameter session. The final commands also form a session
  // when the input has no trailing #RUN.
  std::vector<std::string> lines;
  std::istringstream stream(paramString);
  std::string line;
  while (std::getline(stream, line)) {
    lines.push_back(line + "\n");
  }

  auto has_content = [](const std::string& s) {
    return s.find_first_not_of(" \t\r\n") != std::string::npos;
  };

  std::vector<std::string> sessions;
  std::string currentSession;

  for (size_t i = 0; i < lines.size(); ++i) {
    const std::string& l = lines[i];
    auto pos = l.find_first_not_of(" \t");
    bool isRun =
        (pos != std::string::npos && l.rfind("#RUN", pos) == pos &&
         (pos + 4 >= l.size() || l[pos + 4] == ' ' || l[pos + 4] == '\t' ||
          l[pos + 4] == '\r' || l[pos + 4] == '\n' || l[pos + 4] == '#'));

    currentSession += l;

    if (isRun) {
      if (has_content(currentSession)) {
        sessions.push_back(currentSession);
      }
      currentSession.clear();
    }
  }

  if (has_content(currentSession)) {
    sessions.push_back(currentSession);
  }

  if (sessions.empty()) {
    sessions.push_back(paramString);
  }

  return sessions;
}

void read_stop_criteria(const std::string& paramString, int& maxIter,
                        double& timeMax) {
  ReadParam reader;
  reader = paramString;

  reader.set_verbose(false);

  std::string command;
  while (reader.get_next_command(command)) {
    if (command == "#STOP") {
      reader.set_verbose(amrex::ParallelDescriptor::IOProcessor());
      reader.read_var("MaxIter", maxIter);
      reader.read_var("TimeMax", timeMax);
      break;
    }
  }
}

} // namespace

int main(int argc, char* argv[]) {
  using namespace amrex;

  Initialize(argc, argv);
  {
    if (ParallelDescriptor::MyProc() == 0)
      print_git_info();

    // 1. Read PARAM.in, broadcast it, and create standalone output directories.
    std::string paramString = prepare_standalone_run();
    std::vector<std::string> sessions = split_sessions(paramString);

    // 2. Initialize Domain
    fleksDomains.add_new_domain();
    fleksDomains.select(0);
    Domain& domain = fleksDomains(0);

    domain.init(0.0, 1, sessions[0], {}, {}, {}, /*isStandalone=*/true);

    // Turn on all cells.
    domain.receive_grid_info();

    // Create grids for all components.
    domain.regrid();

    // 3. Set Initial Conditions
    domain.set_ic();

    // #STOP cycle and time limits are absolute endpoints; each later session
    // updates the same Domain before advancing it to its next endpoint.
    for (size_t iSession = 0; iSession < sessions.size(); ++iSession) {
      if (iSession > 0) {
        domain.update_param(sessions[iSession]);
      }

      int maxIter = -1;
      double timeMax = 0.0;
      read_stop_criteria(sessions[iSession], maxIter, timeMax);

      if (maxIter < 0 && timeMax <= 0.0) {
        continue;
      }

      while ((maxIter < 0 || domain.tc->get_cycle() < maxIter) &&
             (timeMax <= 0.0 || domain.tc->get_time_si() < timeMax - 1e-10)) {
        domain.update();
      }
    }

    amrex::Print() << "\nSimulation finished at time = "
                   << domain.tc->get_time_si() << std::endl;

    // 5. Final output
    domain.write_plots(true);
    fleksDomains.clear();
  }

  Finalize();
  return 0;
}
