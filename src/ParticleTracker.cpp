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
      tps->add_test_particles_from_fluid(cellStatus, tpStates);
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
    for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
      tps->move_and_save_particles(nodeE[iLevTest], nodeB[iLevTest], 0, 0,
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

    for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
      tps->move_and_save_particles(nodeE[iLevTest], nodeB[iLevTest], tc->get_dt(),
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
        tps->add_test_particles_from_fluid(cellStatus, tpStates);
      }
    }
  }
}

void ParticleTracker::update_field(Pic& pic) {
  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    MultiFab::Copy(nodeE[iLevTest], pic.nodeE[iLevTest], 0, 0,
                   nodeE[iLevTest].nComp(), nodeE[iLevTest].nGrow());
    MultiFab::Copy(nodeB[iLevTest], pic.nodeB[iLevTest], 0, 0,
                   nodeB[iLevTest].nComp(), nodeB[iLevTest].nGrow());
  }
}

void ParticleTracker::update_cell_status(Pic& pic) {

  if (cellStatus.empty())
    return;

  iMultiFab::Copy(cellStatus, pic.cellStatus, 0, 0, cellStatus.nComp(),
                  cellStatus.nGrow());

  for (MFIter mfi(cellStatus); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const auto& cellArrPT = cellStatus[mfi].array();
    const auto& cellArrPIC = pic.cellStatus[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {

          for (int kk = k - 1; kk <= k + 1; kk++)
            for (int jj = j - 1; jj <= j + 1; jj++)
              for (int ii = i - 1; ii <= i + 1; ii++) {
                if (cellArrPIC(ii, jj, kk) == iBoundary_)
                  cellArrPT(i, j, k) = iAddPTParticle_;
              }
        }
  }
}

void ParticleTracker::post_process_param() {
  savectr = std::unique_ptr<PlotCtr>(
      new PlotCtr(tc.get(), gridID, -1, nPTRecord * dnSave));
  savectr->set_multiple(dnSave);
}

void ParticleTracker::regrid(const BoxArray& ptRegionIn,
                             const BoxArray& centerBAIn,
                             const DistributionMapping& dmIn, Pic& pic) {
  std::string nameFunc = "PT::regrid";

  timing_func(nameFunc);

  if (!usePT)
    return;

  // Why need 'isGridInitialized'? See the explanation in Domain::regrid().
  if (centerBAIn == cGrids[0] && isGridInitialized)
    return;

  isGridEmpty = ptRegionIn.empty();

  activeRegionBA = ptRegionIn;
  baseGrid = centerBAIn;  
  
  if (baseGrid.empty()) {
    cGrids.clear();
    cGrids.push_back(amrex::BoxArray());
  } else {
    // This method will call MakeNewLevelFromScratch() and
    // PostProcessBaseGrids()
    InitFromScratch(tc->get_time());
    SetDistributionMap(0, dmIn);
  }

  calc_node_grids();

  if (nodeB.empty()) {
    nodeB.resize(max_level + 1);
  }
  if (nodeE.empty()) {
    nodeE.resize(max_level + 1);
  }

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    distribute_FabArray(nodeE[iLevTest], nGrids[0], DistributionMap(iLevTest), 3,
                        nGst, false);
    distribute_FabArray(nodeB[iLevTest], nGrids[0], DistributionMap(iLevTest), 3,
                        nGst, false);
  }
  distribute_FabArray(cellStatus, cGrids[0], DistributionMap(0), 1, nGst, false);

  update_cell_status(pic);

  //--------------test particles-----------------------------------
  if (parts.empty()) {
    for (int i = 0; i < nSpecies; i++) {
      auto ptr = std::unique_ptr<TestParticles>(new TestParticles(
          this, fi.get(), tc.get(), i, fi->get_species_charge(i),
          fi->get_species_mass(i), gridID));
      ptr->set_region_range(activeRegionBA);
      ptr->set_ppc(nTPPerCell);
      ptr->set_interval(nTPIntervalCell);
      ptr->set_particle_region(sPartRegion);
      ptr->set_relativistic(isRelativistic);
      parts.push_back(std::move(ptr));
    }
  } else {
    for (int i = 0; i < nSpecies; i++) {
      // Label the particles outside the OLD PIC region.
      parts[i]->label_particles_outside_ba();
      parts[i]->SetParticleBoxArray(0, cGrids[0]);
      parts[i]->set_region_range(activeRegionBA);
      parts[i]->SetParticleDistributionMap(0, DistributionMap(0));
      // Label the particles outside the NEW PIC region.
      parts[i]->label_particles_outside_ba_general();
      parts[i]->Redistribute();
    }
  }

  { // Copy cell Status to Particles objects.
    for (int i = 0; i < nSpecies; i++) {
      distribute_FabArray(parts[i]->cellStatus, cGrids[0], DistributionMap(0), 1,
                          nGst, false);

      if (!cellStatus.empty()) {
        iMultiFab::Copy(parts[i]->cellStatus, cellStatus, 0, 0,
                        cellStatus.nComp(), cellStatus.nGrow());
      }
    }
  }
  //--------------test particles-----------------------------------

  isGridInitialized = true;
}

void ParticleTracker::save_restart_data() {
  if (isGridEmpty || !usePT)
    return;

  std::string restartDir = component + "/restartOUT/";

  bool doSavePlot = savectr->is_time_to(true);
  for (int iPart = 0; iPart < parts.size(); iPart++) {
    // Keep the following two lines for safety.
    parts[iPart]->label_particles_outside_ba();
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
  // writer.set_iRegion(gridID);
  // writer.set_domainMin_D({ { 0, 0, 0 } });

  // writer.set_domainMax_D({ { 1, 1, 1 } });

  // const Real* dx = Geom(0).CellSize();
  // writer.set_dx_D({ { dx[ix_], dx[iy_], dx[iz_] } });
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
