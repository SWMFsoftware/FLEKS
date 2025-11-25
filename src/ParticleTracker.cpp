#include "ParticleTracker.h"
#include "GridUtility.h"

using namespace amrex;

void ParticleTracker::set_ic(Pic& pic) {
  if (isGridEmpty)
    return;

  complete_parameters();

  update_field(pic);

  for (int i = 0; i < parts.size(); ++i) {
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

  for (int i = 0; i < parts.size(); ++i) {
    // The initial state is special. Here, we do a fake update with dt=0, and
    // write the initial state to disk.
    bool doSave = true;
    auto& tps = parts[i];
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      tps->move_and_save_particles(nodeE[iLev], nodeB[iLev], 0, 0,
                                   tc->get_time_si(), doSave);
    }
    tps->write_particles(tc->get_cycle());
  }
}

//==========================================================
void ParticleTracker::write_log(bool doForce, bool doCreateFile) {
  if (isGridEmpty)
    return;

  if (doCreateFile && ParallelDescriptor::IOProcessor()) {
    std::stringstream ss;
    ss << component << "/plots/log_pt_n" << std::setfill('0') << std::setw(8)
       << tc->get_cycle() << ".log";
    logFile = ss.str();
    std::ofstream of(logFile.c_str());
    of << "time nStep";
    for (int i = 0; i < nSpecies; ++i)
      of << " mass_" << i << " moment_x_" << i << " moment_y_" << i
         << " moment_z_" << i << " energy_" << i;

    of << std::endl;
    of.close();
  }

  if (tc->ptLog.is_time_to(doForce)) {

    Vector<std::array<Real, 5> > moments;
    for (int i = 0; i < parts.size(); ++i) {
      moments.push_back(parts[i]->total_moments());
    }

    if (ParallelDescriptor::IOProcessor()) {
      std::ofstream of(logFile.c_str(), std::fstream::app);
      of.precision(15);
      of << std::scientific;
      of << tc->get_time_si() << "\t" << tc->get_cycle();

      for (int i = 0; i < parts.size(); ++i) {
        for (auto& m : moments[i])
          of << "\t" << m;
      }
      of << std::endl;
      of.close();
    }
  }
}

//==========================================================
void ParticleTracker::update(Pic& pic, bool doReport) {
  std::string funcName = "PTracker::update";
  timing_func(funcName);

  if (isGridEmpty)
    return;

  if (doReport) {
    Print() << printPrefix
            << " updating test particles. t =  " << std::setprecision(6)
            << tc->get_time_si() << " (s), cycle = " << tc->get_cycle()
            << std::endl;
  }

  update_field(pic);

  bool doSave = savectr->is_time_to();
  for (int i = 0; i < parts.size(); ++i) {
    auto& tps = parts[i];

    for (int iLev = 0; iLev < n_lev(); iLev++) {
      tps->move_and_save_particles(nodeE[iLev], nodeB[iLev], tc->get_dt(),
                                   tc->get_next_dt(), tc->get_time_si(),
                                   tc->get_cycle() % dnSave[i] == 0);
    }

    if (doSave) {
      auto nt = tps->TotalNumberOfParticles();
      auto n0 = tps->init_particle_number();
      Print() << printPrefix << "particle number of species " << i
              << ": initial = " << n0 << ". current = " << nt
              << ". ratio = " << (n0 > 0 ? (double)nt / n0 : 0.0) << std::endl;

      tps->write_particles(tc->get_cycle());

      // Refill test particles if necessary.
      if (doInitFromPIC) {
        tps->add_test_particles_from_pic(pic.get_particle_pointer(i));
      } else if (tps->TotalNumberOfParticles() <
                 launchThreshold[i] * tps->init_particle_number()) {
        tps->add_test_particles_from_fluid(tpStates);
      }
    }
  }
}

void ParticleTracker::update_field(Pic& pic) {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    MultiFab::Copy(nodeE[iLev], pic.nodeE[iLev], 0, 0, nodeE[iLev].nComp(),
                   nodeE[iLev].nGrow());
    MultiFab::Copy(nodeB[iLev], pic.nodeB[iLev], 0, 0, nodeB[iLev].nComp(),
                   nodeB[iLev].nGrow());
  }
}

void ParticleTracker::post_process_param() {
  int min_dnSave = dnSave[0];
  for (int i = 1; i < nSpecies; ++i) {
    if (dnSave[i] < min_dnSave) {
      min_dnSave = dnSave[i];
    }
  }
  savectr = std::unique_ptr<PlotCtr>(
      new PlotCtr(ParallelDescriptor::Communicator(), tc, gridID, -1,
                  nPTRecord * min_dnSave));
  savectr->set_multiple(min_dnSave);
}

void ParticleTracker::pre_regrid() {
  if (!parts.empty()) {
    for (int i = 0; i < nSpecies; ++i) {
      // Label the particles outside the OLD PIC region. It should be called
      // before active region is updated.
      parts[i]->label_particles_outside_active_region();
    }
  }
}

void ParticleTracker::post_regrid() {

  if (nodeB.empty()) {
    nodeB.resize(n_lev_max());
  }
  if (nodeE.empty()) {
    nodeE.resize(n_lev_max());
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    distribute_FabArray(nodeE[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst, false);
    distribute_FabArray(nodeB[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst, false);
  }

  distribute_grid_arrays();

  //--------------test particles-----------------------------------
  if (parts.empty()) {
    for (int i = 0; i < nSpecies; ++i) {
      auto ptr = std::unique_ptr<TestParticles>(
          new TestParticles(this, fi, tc, i, fi->get_species_charge(i),
                            fi->get_species_mass(i), gridID));
      ptr->set_ppc(nTPPerCell);
      ptr->set_interval(nTPIntervalCell);
      ptr->set_particle_region(sRegion, tpShapes);
      ptr->set_relativistic(isRelativistic);
      parts.push_back(std::move(ptr));
    }
    Print() << gridName << " pt: Number of test particle species: " << nSpecies
            << std::endl;
    Print() << gridName
            << " pt: TPSAVE parameters (nPTRecord, ptRecordSize): " << nPTRecord
            << ", " << ptRecordSize << std::endl;
    Print() << gridName << " pt: TPREGION: " << sRegion << std::endl;
  } else {
    for (int i = 0; i < nSpecies; ++i) {
      // Label the particles outside the NEW PIC region.
      parts[i]->label_particles_outside_active_region_general();
      parts[i]->redistribute_particles();
    }
  }
  //--------------test particles-----------------------------------
}

void ParticleTracker::set_tp_init_shapes(
    amrex::Vector<std::shared_ptr<Shape> >& shapes) {
  tpShapes = shapes;
}

void ParticleTracker::save_restart_data() {
  if (isGridEmpty)
    return;

  bool doSavePlot = savectr->is_time_to(true);
  for (int iPart = 0; iPart < parts.size(); iPart++) {
    // Keep the following two lines for safety.
    parts[iPart]->label_particles_outside_active_region();
    parts[iPart]->redistribute_particles();

    if (doSavePlot) {
      parts[iPart]->write_particles(tc->get_cycle());
    }

    parts[iPart]->Checkpoint(fi->get_restart_out_dir(),
                             gridName + "_test_particles" +
                                 std::to_string(iPart));
  }
}

void ParticleTracker::read_restart() {
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

  if (command == "#TPPARTICLES") {
    param.read_var("npcelx", nTPPerCell[ix_]);
    param.read_var("npcely", nTPPerCell[iy_]);
    if (nDim == 3)
      param.read_var("npcelz", nTPPerCell[iz_]);
  } else if (command == "#TPCELLINTERVAL") {
    param.read_var("nIntervalX", nTPIntervalCell[ix_]);
    param.read_var("nIntervalY", nTPIntervalCell[iy_]);
    if (nDim == 3)
      param.read_var("nIntervalZ", nTPIntervalCell[iz_]);
  } else if (command == "#TPREGION") {
    param.read_var("region", sRegion);
  } else if (command == "#TPSAVE") {
    int iSpecies;
    param.read_var("iSpecies", iSpecies);
    if (iSpecies >= nSpecies)
      amrex::Abort("Error: iSpecies is out of bound in #TPSAVE.");
    param.read_var("IOUnit", sIOUnit);
    param.read_var("dnSave", dnSave[iSpecies]);
    param.read_var("launchThreshold", launchThreshold[iSpecies]);
  } else if (command == "#TPRELATIVISTIC") {
    param.read_var("isRelativistic", isRelativistic);
  } else if (command == "#TPSTATESI") {
    double si2noV = fi->get_Si2NoV();
    int nState;
    param.read_var("nState", nState);
    for (int i = 0; i < nState; ++i) {
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
      for (int i = 0; i < nList; ++i) {
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
