#include "show_git_info.h"
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "Constants.h"
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
      std::filesystem::create_directories(component + "/plots");
      std::filesystem::create_directories(component + "/restartOUT");
    }

    // 2. Initialize Domain
    fleksDomains.add_new_domain();
    fleksDomains.select(0);
    Domain& domain = fleksDomains(0);
    domain.init(0.0, 1, paramString);

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

      // Species diagnostics loop
      if (domain.pic && domain.pic->has_particles()) {
        int nSpecies = domain.pic->get_nSpecies();
        for (int iS = 0; iS < nSpecies; iS++) {
          if (!domain.pic->get_particle_pointer(iS)) continue;
          double local_weight = 0;
          double local_vx = 0, local_vy = 0, local_vz = 0;
          double local_ke = 0;
          long local_macro = 0;

          for (int iLev = 0; iLev < domain.pic->n_lev(); iLev++) {
            for (amrex::ParIter<nPicPartReal, nPicPartInt> pti(*domain.pic->get_particle_pointer(iS), iLev); pti.isValid(); ++pti) {
              const auto& tile = pti.GetArrayOfStructs();
              for (const auto& p : tile) {
                if (p.id() < 0) continue;
                double q = p.rdata(PicParticles::iqp_); // weight
                double vx = p.rdata(PicParticles::iup_);
                double vy = p.rdata(PicParticles::ivp_);
                double wp = p.rdata(PicParticles::iwp_);
                local_macro++;
                local_weight += std::abs(q);
                local_vx += vx * std::abs(q);
                local_vy += vy * std::abs(q);
                local_vz += wp * std::abs(q);
                local_ke += 0.5 * (vx*vx + vy*vy + wp*wp) * std::abs(q);
              }
            }
          }

          double global_weight = local_weight;
          double global_vx = local_vx;
          double global_vy = local_vy;
          double global_vz = local_vz;
          double global_ke = local_ke;
          long global_macro = local_macro;

          amrex::ParallelDescriptor::ReduceRealSum(global_weight);
          amrex::ParallelDescriptor::ReduceRealSum(global_vx);
          amrex::ParallelDescriptor::ReduceRealSum(global_vy);
          amrex::ParallelDescriptor::ReduceRealSum(global_vz);
          amrex::ParallelDescriptor::ReduceRealSum(global_ke);
          amrex::ParallelDescriptor::ReduceLongSum(global_macro);

          if (global_weight > 0) {
            global_vx /= global_weight;
            global_vy /= global_weight;
            global_vz /= global_weight;
          }

          if (amrex::ParallelDescriptor::IOProcessor()) {
            std::cout << std::scientific << std::setprecision(6)
                      << "DIAGNOSTIC:"
                      << " Species=" << iS
                      << " Time=" << domain.tc->get_time_si()
                      << " Cycle=" << domain.tc->get_cycle()
                      << " MacroParticles=" << global_macro
                      << " PhysParticles=" << global_weight
                      << " MeanVx=" << global_vx
                      << " MeanVy=" << global_vy
                      << " MeanVz=" << global_vz
                      << " KineticEnergy=" << global_ke
                      << "\n";
          }
        }
      }

      if (domain.pic) {
        if (amrex::ParallelDescriptor::IOProcessor()) {
          double max_By = 0.0;
          double max_Bz = 0.0;
          for (int iLev = 0; iLev < domain.pic->n_lev(); iLev++) {
            max_By = std::max(max_By, domain.pic->get_nodeB()[iLev].norm0(1));
            max_Bz = std::max(max_Bz, domain.pic->get_nodeB()[iLev].norm0(2));
          }
          std::cout << std::scientific << std::setprecision(6)
                    << "DIAGNOSTIC_FIELD:"
                    << " Time=" << domain.tc->get_time_si()
                    << " Cycle=" << domain.tc->get_cycle()
                    << " MaxBy=" << max_By
                    << " MaxBz=" << max_Bz
                    << "\n";
        }
      }
    }

    amrex::Print() << "\nSimulation finished at time = "
                   << domain.tc->get_time_si() << std::endl;

    // 5. Final output
    domain.write_plots(true);
  }

  fleksDomains.clear();
  Finalize();
  return 0;
}
