#include <algorithm>

#include "GridUtility.h"
#include "ParticleTracker.h"

using namespace amrex;

ParticleTracker::~ParticleTracker() {
  if (isGridEmpty || isNewGrid || !savectr)
    return;

  bool doSave = savectr->is_time_to(true);
  for (auto& tps : parts) {
    if (doSave) {
      tps->write_particles(tc->get_cycle());
    }
  }
}

void ParticleTracker::set_ic(Pic& pic) {
  if (isGridEmpty)
    return;

  complete_parameters();

  update_field(pic, ptRecordSize > 13);

  for (int i = 0; i < parts.size(); ++i) {
    auto& tps = parts[i];
    if (pInfo->doInitFromPIC) {
      tps->read_test_particle_list(pInfo->listFiles);
      tps->add_test_particles_from_pic(pic.get_particle_pointer(i));
    } else {
      tps->add_test_particles_from_fluid(pInfo->tpStates);
    }
    tps->update_initial_particle_number();

    Print() << printPrefix << " initial particle # is "
            << tps->init_particle_number() << " for species " << i << std::endl;

    // The initial state is special. Here, we do a fake update with dt=0, and
    // write the initial state to disk.
    bool doSave = true;
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      const MultiFab* jacPtr =
          (ptRecordSize > 13 && !nodeJacB.empty()) ? &nodeJacB[iLev] : nullptr;
      if (pic.useHybridPIC) {
        tps->move_and_save_particles_cell_centered(
            iLev, centerE[iLev], centerB[iLev], 0, 0, tc->get_time_si(), doSave,
            jacPtr);
      } else {
        tps->move_and_save_particles(iLev, nodeE[iLev], nodeB[iLev], 0, 0,
                                     tc->get_time_si(), doSave, jacPtr);
      }
    }
    if (tps->get_dt_save() > 0.0) {
      tps->advance_next_save();
    }
    tps->write_particles(tc->get_cycle());
  }
}

//==========================================================
void ParticleTracker::write_log(bool doForce, bool doCreateFile) {
  if (isGridEmpty)
    return;

  if (doCreateFile) {
    std::stringstream ss;
    ss << component << "/plots/log_pt_n" << std::setfill('0') << std::setw(8)
       << tc->get_cycle() << ".log";
    logFile = ss.str();
    if (ParallelDescriptor::IOProcessor()) {
      std::ofstream of(logFile.c_str());
      of << "time nStep";
      for (int i = 0; i < parts.size(); ++i)
        of << " mass_" << i << " moment_x_" << i << " moment_y_" << i
           << " moment_z_" << i << " energy_" << i;

      of << std::endl;
      of.close();
    }
  }

  if (tc->ptLog.is_time_to(doForce)) {

    Vector<std::array<Real, 5> > moments;
    moments.reserve(parts.size());
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

  bool needJacobian = false;
  for (int i = 0; i < parts.size(); ++i) {
    if (parts[i]->is_time_to_record(tc->get_time_si(), tc->get_cycle(),
                                    pInfo->dnSave[i])) {
      needJacobian = true;
      break;
    }
  }

  update_field(pic, needJacobian);

  bool doSave = savectr->is_time_to();
  for (int i = 0; i < parts.size(); ++i) {
    auto& tps = parts[i];
    bool doRecord = tps->is_time_to_record(tc->get_time_si(), tc->get_cycle(),
                                           pInfo->dnSave[i]);

    for (int iLev = 0; iLev < n_lev(); iLev++) {
      const MultiFab* jacPtr =
          (ptRecordSize > 13 && !nodeJacB.empty() && doRecord) ? &nodeJacB[iLev]
                                                               : nullptr;
      if (pic.useHybridPIC) {
        tps->move_and_save_particles_cell_centered(
            iLev, centerE[iLev], centerB[iLev], tc->get_dt(), tc->get_next_dt(),
            tc->get_time_si(), doRecord, jacPtr);
      } else {
        tps->move_and_save_particles(iLev, nodeE[iLev], nodeB[iLev],
                                     tc->get_dt(), tc->get_next_dt(),
                                     tc->get_time_si(), doRecord, jacPtr);
      }
    }

    if (doRecord && tps->get_dt_save() > 0.0) {
      tps->advance_next_save();
    }

    if (doSave) {
      auto nt = tps->TotalNumberOfParticles();
      auto n0 = tps->init_particle_number();
      Print() << printPrefix << "particle number of species " << i
              << ": initial = " << n0 << ". current = " << nt
              << ". ratio = " << (n0 > 0 ? (double)nt / n0 : 0.0) << std::endl;

      tps->write_particles(tc->get_cycle());

      // Refill test particles if necessary.
      if (pInfo->doInitFromPIC) {
        tps->add_test_particles_from_pic(pic.get_particle_pointer(i));
      } else if (nt < pInfo->launchThreshold[i] * n0) {
        tps->add_test_particles_from_fluid(pInfo->tpStates);
      }
    }
  }
}

void ParticleTracker::update_field(Pic& pic, bool needJacobian) {
  if (pic.useHybridPIC) {
    // Hybrid: gather from the live cell-centred fields.
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      MultiFab::Copy(centerE[iLev], pic.centerEhybrid[iLev], 0, 0,
                     centerE[iLev].nComp(), centerE[iLev].nGrow());
      MultiFab::Copy(centerB[iLev], pic.centerB[iLev], 0, 0,
                     centerB[iLev].nComp(), centerB[iLev].nGrow());
    }
  } else {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      MultiFab::Copy(nodeE[iLev], pic.nodeE[iLev], 0, 0, nodeE[iLev].nComp(),
                     nodeE[iLev].nGrow());
      MultiFab::Copy(nodeB[iLev], pic.nodeB[iLev], 0, 0, nodeB[iLev].nComp(),
                     nodeB[iLev].nGrow());
    }
  }

  // If magnetic field gradient is requested AND needed this step, compute
  // Jacobian from centerB
  if (needJacobian && ptRecordSize > 13) {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      jacobian_center_to_node(pic.centerB[iLev], nodeJacB[iLev],
                              Geom(iLev).InvCellSize());
      nodeJacB[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }
}

void ParticleTracker::post_process_param() {
  const int nSpecies = fi->get_nS();
  Real min_dtSave = -1.0;
  for (int i = 0; i < nSpecies; ++i) {
    if (pInfo->dtSave[i] > 0.0) {
      if (min_dtSave < 0.0 || pInfo->dtSave[i] < min_dtSave) {
        min_dtSave = pInfo->dtSave[i];
      }
    }
  }

  if (min_dtSave > 0.0) {
    savectr = std::make_unique<PlotCtr>(ParallelDescriptor::Communicator(), tc,
                                        gridID, nPTRecord * min_dtSave, -1);
  } else {
    int min_dnSave = *std::min_element(pInfo->dnSave.begin(),
                                       pInfo->dnSave.begin() + nSpecies);
    savectr = std::make_unique<PlotCtr>(ParallelDescriptor::Communicator(), tc,
                                        gridID, -1, nPTRecord * min_dnSave);
    savectr->set_multiple(min_dnSave);
  }
}

void ParticleTracker::pre_regrid() {
  for (auto& tp : parts) {
    tp->label_particles_outside_active_region();
  }
}

void ParticleTracker::post_regrid() {
  if (nodeB.empty()) {
    nodeB.resize(n_lev_max());
  }
  if (nodeE.empty()) {
    nodeE.resize(n_lev_max());
  }
  if (centerB.empty()) {
    centerB.resize(n_lev_max());
  }
  if (centerE.empty()) {
    centerE.resize(n_lev_max());
  }
  if (ptRecordSize > 13 && nodeJacB.empty()) {
    nodeJacB.resize(n_lev_max());
  }

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    distribute_FabArray(nodeE[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst, false);
    distribute_FabArray(nodeB[iLev], nGrids[iLev], DistributionMap(iLev), 3,
                        nGst, false);
    distribute_FabArray(centerE[iLev], cGrids[iLev], DistributionMap(iLev), 3,
                        nGst, false);
    distribute_FabArray(centerB[iLev], cGrids[iLev], DistributionMap(iLev), 3,
                        nGst, false);
    if (ptRecordSize > 13) {
      distribute_FabArray(nodeJacB[iLev], nGrids[iLev], DistributionMap(iLev),
                          9, nGst, false);
    }
  }

  distribute_grid_arrays();

  //--------------test particles-----------------------------------
  const int nSpecies = fi->get_nS();

  if (parts.empty()) {
    for (int i = 0; i < nSpecies; ++i) {
      auto ptr = std::make_unique<TestParticles>(
          this, fi, tc, i, fi->get_species_charge(i), fi->get_species_mass(i),
          gridID);
      ptr->set_ppc(pInfo->nTPPerCell);
      ptr->set_interval(pInfo->nTPIntervalCell);
      ptr->set_particle_region(pInfo->sRegion, tpShapes);
      ptr->set_relativistic(pInfo->isRelativistic);
      ptr->set_dt_save(pInfo->dtSave[i], tc->get_time_si());
      parts.push_back(std::move(ptr));
    }
    Print() << gridName << " pt: Number of test particle species: " << nSpecies
            << std::endl;
    Print() << gridName
            << " pt: TPSAVE parameters (nPTRecord, ptRecordSize): " << nPTRecord
            << ", " << ptRecordSize << std::endl;
    Print() << gridName << " pt: TPREGION: " << pInfo->sRegion << std::endl;
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
    if (iPart < (int)pInfo->initPartNumber.size())
      parts[iPart]->init_particle_number(pInfo->initPartNumber[iPart]);
    else
      parts[iPart]->init_particle_number(0);
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
  writer.set_plotString("3d fluid test_particle real4 " + pInfo->sIOUnit);
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
