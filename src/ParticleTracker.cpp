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
    tps->add_test_particles(cellStatus);
    tps->update_initial_particle_number();

    Print() << printPrefix << " initial particle # is "
            << tps->init_particle_number() << " for species " << i << std::endl;
  }
}

void ParticleTracker::update(Pic& pic) {
  std::string funcName = "PTracker::update";
  timing_func(funcName);

  if (isGridEmpty || !usePT)
    return;

  update_field(pic);
  bool doSave = savectr->is_time_to();
  for (int i = 0; i < parts.size(); i++) {
    auto& tps = parts[i];
    tps->move_and_save_particles(nodeE, nodeB, tc->get_dt(), tc->get_next_dt(),
                                 tc->get_time_si(),
                                 tc->get_cycle() % dnSave == 0);

    if (doSave) {
      Print() << printPrefix << "particle number of species " << i
              << ": initial = " << tps->init_particle_number()
              << ". current = " << tps->TotalNumberOfParticles() << ". ratio = "
              << (double)tps->TotalNumberOfParticles() /
                     tps->init_particle_number()
              << std::endl;

      tps->write_particles(tc->get_cycle());
      // Refill test particles if necessary.
      if (tps->TotalNumberOfParticles() < 0.5 * tps->init_particle_number()) {
        tps->add_test_particles(cellStatus);
      }
    }
  }
}

void ParticleTracker::update_field(Pic& pic) {
  MultiFab::Copy(nodeE, pic.nodeE, 0, 0, nodeE.nComp(), nodeE.nGrow());
  MultiFab::Copy(nodeB, pic.nodeB, 0, 0, nodeB.nComp(), nodeB.nGrow());
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

void ParticleTracker::init(std::shared_ptr<FluidInterface>& fluidIn,
                           std::shared_ptr<TimeCtr>& tcIn, int domainIDIn) {
  tc = tcIn;
  fluidInterface = fluidIn;
  nSpecies = fluidInterface->get_nS();

  domainID = domainIDIn;

  {
    std::stringstream ss;
    ss << "FLEKS" << domainID;
    domainName = ss.str();
    printPrefix = domainName + " PT: ";
  }
}

void ParticleTracker::post_process_param() {
  savectr = std::unique_ptr<PlotCtr>(
      new PlotCtr(tc.get(), domainID, -1, nPTRecord * dnSave));
}

//==========================================================
void ParticleTracker::set_geom(int nGstIn, const Geometry& geomIn) {
  if (!usePT)
    return;

  nGst = nGstIn;
  geom = geomIn;
}

void ParticleTracker::regrid(const BoxArray& ptRegionIn,
                             const BoxArray& centerBAIn,
                             const DistributionMapping& dmIn, Pic& pic) {
  std::string nameFunc = "PT::regrid";

  timing_func(nameFunc);

  if (!usePT)
    return;

  // Why need 'isGridInitialized'? See the explanation in Domain::regrid().
  if (centerBAIn == centerBA && isGridInitialized)
    return;

  isGridEmpty = ptRegionIn.empty();

  ptRegionBA = ptRegionIn;
  centerBA = centerBAIn;
  nodeBA = convert(centerBA, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });
  dm = dmIn;

  distribute_FabArray(nodeE, nodeBA, dm, 3, nGst, false);
  distribute_FabArray(nodeB, nodeBA, dm, 3, nGst, false);
  distribute_FabArray(cellStatus, centerBA, dm, 1, nGst, false);

  update_cell_status(pic);

  //--------------test particles-----------------------------------
  if (parts.empty()) {
    for (int i = 0; i < nSpecies; i++) {
      auto ptr = std::unique_ptr<TestParticles>(new TestParticles(
          ptRegionBA, geom, dm, centerBA, fluidInterface.get(), tc.get(), i,
          fluidInterface->get_species_charge(i),
          fluidInterface->get_species_mass(i), domainID));
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
      parts[i]->SetParticleBoxArray(0, centerBA);
      parts[i]->set_region_ba(ptRegionBA);
      parts[i]->SetParticleDistributionMap(0, dm);
      // Label the particles outside the NEW PIC region.
      parts[i]->label_particles_outside_ba_general();
      parts[i]->Redistribute();
    }
  }

  { // Copy cell Status to Particles objects.
    for (int i = 0; i < nSpecies; i++) {
      distribute_FabArray(parts[i]->cellStatus, centerBA, dm, 1, nGst, false);

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

  std::string restartDir = "PC/restartOUT/";

  bool doSavePlot = savectr->is_time_to(true);
  for (int iPart = 0; iPart < parts.size(); iPart++) {
    // Keep the following two lines for safety.
    parts[iPart]->label_particles_outside_ba();
    parts[iPart]->Redistribute();

    if (doSavePlot) {
      parts[iPart]->write_particles(tc->get_cycle());
    }

    parts[iPart]->Checkpoint(restartDir, domainName + "_test_particles" +
                                             std::to_string(iPart));
  }
}

void ParticleTracker::read_restart() {
  if (!usePT)
    return;

  std::string restartDir = "PC/restartIN/";
  for (int iPart = 0; iPart < parts.size(); iPart++) {
    parts[iPart]->Restart(restartDir, domainName + "_test_particles" +
                                          std::to_string(iPart));
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
  writer.set_nDim(fluidInterface->get_fluid_dimension());
  // writer.set_iRegion(domainID);
  // writer.set_domainMin_D({ { 0, 0, 0 } });

  // writer.set_domainMax_D({ { 1, 1, 1 } });

  // const Real* dx = geom.CellSize();
  // writer.set_dx_D({ { dx[ix_], dx[iy_], dx[iz_] } });
  writer.set_units(fluidInterface->get_No2SiL(), fluidInterface->get_No2SiV(),
                   fluidInterface->get_No2SiB(), fluidInterface->get_No2SiRho(),
                   fluidInterface->get_No2SiP(), fluidInterface->get_No2SiJ(),
                   fluidInterface->get_rPlanet_SI());
  writer.set_No2NoL(fluidInterface->get_MhdNo2NoL());
  //--------------------------------------------------
  writer.init();

  for (auto& tps : parts) {
    tps->set_IO_units(writer.No2OutTable("X"), writer.No2OutTable("u"),
                      writer.No2OutTable("mass"));
  }
}

void ParticleTracker::save_restart_header(std::ofstream& headerFile) {
  if (!usePT)
    return;

  std::string command_suffix = "_" + domainName + "\n";

  if (ParallelDescriptor::IOProcessor()) {
    headerFile << "#TESTPARTICLENUMBER" + command_suffix;
    for (auto& tp : parts) {
      headerFile << tp->init_particle_number() << "\n";
    }
    headerFile << "\n";
  }
}

void ParticleTracker::read_param(const std::string& command,
                                 ReadParam& readParam) {

  if (command == "#PARTICLETRACKER") {
    readParam.read_var("usePT", usePT);
  } else if (command == "#TPPARTICLES") {
    readParam.read_var("npcelx", nTPPerCell[ix_]);
    readParam.read_var("npcely", nTPPerCell[iy_]);
    readParam.read_var("npcelz", nTPPerCell[iz_]);
  } else if (command == "#TPCELLINTERVAL") {
    readParam.read_var("nIntervalX", nTPIntervalCell[ix_]);
    readParam.read_var("nIntervalY", nTPIntervalCell[iy_]);
    readParam.read_var("nIntervalZ", nTPIntervalCell[iz_]);
  } else if (command == "#TPREGION") {
    readParam.read_var("region", sPartRegion);
  } else if (command == "#TPSAVE") {
    readParam.read_var("IOUnit", sIOUnit);
    readParam.read_var("dnSave", dnSave);
  } else if (command == "#TPRELATIVISTIC") {
    readParam.read_var("isRelativistic", isRelativistic);
  } else if (command == "#TESTPARTICLENUMBER") {
    initPartNumber.clear();
    unsigned long int num;
    for (int iPart = 0; iPart < nSpecies; iPart++) {
      readParam.read_var("Number", num);
      initPartNumber.push_back(num);
    }
  }
}
