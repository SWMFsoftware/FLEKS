#include "ParticleTracker.h"
#include "GridUtility.h"

using namespace amrex;

void ParticleTracker::set_ic(Pic& pic) {
  if (isGridEmpty || !usePT)
    return;

  complete_parameters();

  update_field(pic);

  for (int i = 0; i < parts.size(); i++) {
    auto& tps = parts[i];
    if (doInitFromPIC) {
      tps->read_test_particle_list(listFiles);
      tps->add_test_particles_from_pic(pic.get_particle_pointer(i));
    } else {
      tps->add_test_particles_from_fluid(tpStates);
    }
    tps->update_initial_particle_number();

    Print() << printPrefix << " initial particle # is "
            << tps->init_particle_number() << " for species " << i << std::endl;
  }

  for (int i = 0; i < parts.size(); i++) {
    // The initial state is special. Here, we do a fake update with dt=0, and
    // write the initial state to disk.
    bool doSave = true;
    auto& tps = parts[i];
    for (int iLev = 0; iLev <= finest_level; iLev++) {
      tps->move_and_save_particles(nodeE[iLev], nodeB[iLev], 0, 0,
                                   tc->get_time_si(), doSave);
    }
    tps->write_particles(tc->get_cycle());
  }
}

//==========================================================
void ParticleTracker::write_log(bool doForce, bool doCreateFile) {
  if (isGridEmpty || !usePT)
    return;

  if (doCreateFile && ParallelDescriptor::IOProcessor()) {
    std::stringstream ss;
    ss << component << "/plots/log_pt_n" << std::setfill('0') << std::setw(8)
       << tc->get_cycle() << ".log";
    logFile = ss.str();
    std::ofstream of(logFile.c_str());
    of << "time nStep ";
    for (int i = 0; i < nSpecies; i++)
      of << " mass_" << i << " moment_x_" << i << " moment_y_" << i
         << " moment_z_" << i << " energy_" << i;

    of << std::endl;
    of.close();
  }

  if (tc->ptLog.is_time_to(doForce)) {

    Vector<std::array<Real, 5> > moments;
    for (int i = 0; i < parts.size(); i++) {
      moments.push_back(parts[i]->total_moments());
    }

    if (ParallelDescriptor::IOProcessor()) {
      std::ofstream of(logFile.c_str(), std::fstream::app);
      of.precision(15);
      of << std::scientific;
      of << tc->get_time_si() << "\t" << tc->get_cycle();

      for (int i = 0; i < parts.size(); i++) {
        for (auto& m : moments[i])
          of << "\t" << m;
      }
      of << std::endl;
      of.close();
    }
  }
}

//==========================================================
void ParticleTracker::update(Pic& pic) {
  std::string funcName = "PTracker::update";
  timing_func(funcName);

  if (isGridEmpty || !usePT)
    return;

  Print() << printPrefix
          << " updating test particles. t =  " << std::setprecision(6)
          << tc->get_time_si() << " (s), cycle = " << tc->get_cycle()
          << std::endl;

  update_field(pic);

  bool doSave = savectr->is_time_to();
  for (int i = 0; i < parts.size(); i++) {
    auto& tps = parts[i];

    for (int iLev = 0; iLev <= finest_level; iLev++) {
      tps->move_and_save_particles(nodeE[iLev], nodeB[iLev], tc->get_dt(),
                                   tc->get_next_dt(), tc->get_time_si(),
                                   tc->get_cycle() % dnSave == 0);
    }

    if (doSave) {
      Print() << printPrefix << "particle number of species " << i
              << ": initial = " << tps->init_particle_number()
              << ". current = " << tps->TotalNumberOfParticles() << ". ratio = "
              << (double)tps->TotalNumberOfParticles() /
                     tps->init_particle_number()
              << std::endl;

      tps->write_particles(tc->get_cycle());

      // Refill test particles if necessary.
      if (doInitFromPIC) {
        tps->add_test_particles_from_pic(pic.get_particle_pointer(i));
      } else if (tps->TotalNumberOfParticles() <
                 0.5 * tps->init_particle_number()) {
        tps->add_test_particles_from_fluid(tpStates);
      }
    }
  }
}

void ParticleTracker::update_field(Pic& pic) {
  for (int iLev = 0; iLev <= finest_level; iLev++) {
    MultiFab::Copy(nodeE[iLev], pic.nodeE[iLev], 0, 0, nodeE[iLev].nComp(),
                   nodeE[iLev].nGrow());
    MultiFab::Copy(nodeB[iLev], pic.nodeB[iLev], 0, 0, nodeB[iLev].nComp(),
                   nodeB[iLev].nGrow());
  }
}

void ParticleTracker::post_process_param() {
  savectr = std::unique_ptr<PlotCtr>(
      new PlotCtr(tc.get(), gridID, -1, nPTRecord * dnSave));
  savectr->set_multiple(dnSave);
}

void ParticleTracker::regrid(const BoxArray& region, const Grid* const grid,
                             Pic& pic) {
  std::string nameFunc = "PT::regrid";

  timing_func(nameFunc);

  if (!usePT)
    return;

  // Why need 'isGridInitialized'? See the explanation in Domain::regrid().
  if (region == activeRegion && isGridInitialized)
    return;

  if (!parts.empty()) {
    for (int i = 0; i < nSpecies; i++) {
      // Label the particles outside the OLD PIC region. It should be called
      // before active region is updated.
      parts[i]->label_particles_outside_active_region();
    }
  }

  activeRegion = region;
  isGridEmpty = activeRegion.empty();

  if (isGridEmpty) {
    cGrids.clear();
    cGrids.push_back(amrex::BoxArray());
  } else {
    if (grid) {
      finest_level = grid->finestLevel();
      for (int iLev = 0; iLev < nLev; iLev++) {
        SetBoxArray(iLev, grid->boxArray(iLev));
        SetDistributionMap(iLev, grid->DistributionMap(iLev));
      }
    } else {
      // This method will call MakeNewLevelFromScratch() and
      // PostProcessBaseGrids()
      InitFromScratch(tc->get_time());
    }    
  }

  calc_node_grids();

  print_grid_info();

  if (nodeB.empty()) {
    nodeB.resize(nLev);
  }
  if (nodeE.empty()) {
    nodeE.resize(nLev);
  }

  for (int iLev = 0; iLev <= finest_level; iLev++) {
    distribute_FabArray(nodeE[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst, false);
    distribute_FabArray(nodeB[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst, false);
  }

  distribute_grid_arrays();

  //--------------test particles-----------------------------------
  if (parts.empty()) {
    for (int i = 0; i < nSpecies; i++) {
      auto ptr = std::unique_ptr<TestParticles>(new TestParticles(
          this, fi.get(), tc.get(), i, fi->get_species_charge(i),
          fi->get_species_mass(i), gridID));
      ptr->set_region_range(activeRegion);
      ptr->set_ppc(nTPPerCell);
      ptr->set_interval(nTPIntervalCell);
      ptr->set_particle_region(sPartRegion);
      ptr->set_relativistic(isRelativistic);
      parts.push_back(std::move(ptr));
    }
  } else {
    for (int i = 0; i < nSpecies; i++) {
      parts[i]->SetParticleBoxArray(0, cGrids[0]);
      parts[i]->set_region_range(activeRegion);
      parts[i]->SetParticleDistributionMap(0, DistributionMap(0));
      // Label the particles outside the NEW PIC region.
      parts[i]->label_particles_outside_active_region_general();
      parts[i]->Redistribute();
    }
  }
  //--------------test particles-----------------------------------

  activeRegion = activeRegion.simplified();

  isGridInitialized = true;
}

void ParticleTracker::save_restart_data() {
  if (isGridEmpty || !usePT)
    return;

  std::string restartDir = component + "/restartOUT/";

  bool doSavePlot = savectr->is_time_to(true);
  for (int iPart = 0; iPart < parts.size(); iPart++) {
    // Keep the following two lines for safety.
    parts[iPart]->label_particles_outside_active_region();
    parts[iPart]->Redistribute();

    if (doSavePlot) {
      parts[iPart]->write_particles(tc->get_cycle());
    }

    parts[iPart]->Checkpoint(restartDir, gridName + "_test_particles" +
                                             std::to_string(iPart));
  }
}

void ParticleTracker::read_restart() {
  if (!usePT)
    return;

  std::string restartDir = component + "/restartIN/";
  for (int iPart = 0; iPart < parts.size(); iPart++) {
    parts[iPart]->Restart(restartDir,
                          gridName + "_test_particles" + std::to_string(iPart));
    parts[iPart]->reset_record_counter();
    parts[iPart]->init_particle_number(initPartNumber[iPart]);
  }
  complete_parameters();
}

void ParticleTracker::complete_parameters() {
  PlotWriter& writer = savectr->writer;
  {
    // The plotCtr requires writing at the very begining, which should not
    // happen for test particle. So skip the first saving.
    savectr->is_time_to();
  }

  // Pass information to writers.
  writer.set_plotString("3d fluid test_particle real4 " + sIOUnit);
  writer.set_rank(ParallelDescriptor::MyProc());
  writer.set_nProcs(ParallelDescriptor::NProcs());
  writer.set_nDim(fi->get_fluid_dimension());
  writer.set_units(fi->get_No2SiL(), fi->get_No2SiV(), fi->get_No2SiB(),
                   fi->get_No2SiRho(), fi->get_No2SiP(), fi->get_No2SiJ(),
                   fi->get_rPlanet_SI());
  writer.set_No2NoL(fi->get_MhdNo2NoL());
  //--------------------------------------------------
  writer.init();

  for (auto& tps : parts) {
    tps->set_IO_units(writer.No2OutTable("X"), writer.No2OutTable("u"),
                      writer.No2OutTable("mass"), writer.No2OutTable("B"),
                      writer.No2OutTable("E"));
  }
}

void ParticleTracker::save_restart_header(std::ofstream& headerFile) {
  if (!usePT)
    return;

  std::string command_suffix = "_" + gridName + "\n";

  if (ParallelDescriptor::IOProcessor()) {
    headerFile << "#TESTPARTICLENUMBER" + command_suffix;
    for (auto& tp : parts) {
      headerFile << tp->init_particle_number() << "\n";
    }
    headerFile << "\n";
  }
}

void ParticleTracker::read_param(const std::string& command, ReadParam& param) {

  if (command == "#PARTICLETRACKER") {
    param.read_var("usePT", usePT);
  } else if (command == "#TPPARTICLES") {
    param.read_var("npcelx", nTPPerCell[ix_]);
    param.read_var("npcely", nTPPerCell[iy_]);
    param.read_var("npcelz", nTPPerCell[iz_]);
  } else if (command == "#TPCELLINTERVAL") {
    param.read_var("nIntervalX", nTPIntervalCell[ix_]);
    param.read_var("nIntervalY", nTPIntervalCell[iy_]);
    param.read_var("nIntervalZ", nTPIntervalCell[iz_]);
  } else if (command == "#TPREGION") {
    param.read_var("region", sPartRegion);
  } else if (command == "#TPSAVE") {
    param.read_var("IOUnit", sIOUnit);
    param.read_var("dnSave", dnSave);
  } else if (command == "#TPRELATIVISTIC") {
    param.read_var("isRelativistic", isRelativistic);
  } else if (command == "#TPSTATESI") {
    double si2noV = fi->get_Si2NoV();
    int nState;
    param.read_var("nState", nState);
    for (int i = 0; i < nState; i++) {
      Vel state;
      param.read_var("iSpecies", state.tag);
      param.read_var("vth", state.vth);
      param.read_var("vx", state.vx);
      param.read_var("vy", state.vy);
      param.read_var("vz", state.vz);
      state.vth *= si2noV;
      state.vx *= si2noV;
      state.vy *= si2noV;
      state.vz *= si2noV;
      tpStates.push_back(state);
    }
  } else if (command == "#TPINITFROMPIC") {
    param.read_var("doInitFromPIC", doInitFromPIC);
    if (doInitFromPIC) {
      int nList;
      param.read_var("nList", nList);
      for (int i = 0; i < nList; i++) {
        std::string s;
        param.read_var("list", s);
        listFiles.push_back(s);
      }
    }
  } else if (command == "#TESTPARTICLENUMBER") {
    initPartNumber.clear();
    unsigned long int num;
    for (int iPart = 0; iPart < nSpecies; iPart++) {
      param.read_var("Number", num);
      initPartNumber.push_back(num);
    }
  }
}
