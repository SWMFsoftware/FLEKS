#include "Domain.h"
#include "GridUtility.h"

using namespace amrex;

//========================================================
void Domain::init(double time, const int iDomain,
                  const std::string &paramString,
                  const amrex::Vector<int> &paramInt,
                  const amrex::Vector<double> &paramRegion,
                  const amrex::Vector<double> &paramComm) {
  if (AMREX_SPACEDIM != 3)
    Abort("Error: AMReX should be compiled with 3D configuration!!");

  tc->set_time_si(time);

  gridID = iDomain;
  gridName = std::string("FLEKS") + std::to_string(gridID);
  printPrefix = gridName + ": ";

  param = paramString;

  if (!paramInt.empty())
    if (paramInt[0] == 2)
      isFake2D = true;

  prepare_grid_info(paramRegion);

  if (initFromSWMF) {
    fi = std::make_shared<FluidInterface>(
        gm, amrInfo, nGst, gridID, "fi", paramInt,
        Vector<double>(paramRegion.begin() + 18, paramRegion.end()), paramComm);
  } else {
    fi = std::make_shared<FluidInterface>(gm, amrInfo, nGst, gridID, "fi");
  }

  pic = std::make_unique<Pic>(gm, amrInfo, nGst, fi, tc, gridID);

  pt = std::make_unique<ParticleTracker>(gm, amrInfo, nGst, fi, tc, gridID);

  read_param(false);

#ifdef _PT_COMPONENT_
  stateOH = std::make_shared<FluidInterface>(gm, amrInfo, nGst, gridID,
                                             "stateOH", fi.get());
#endif

  init_time_ctr();

  gridInfo.init(nCell[ix_], nCell[iy_], nCell[iz_], fi->get_nCellPerPatch());

  fi->print_info();

  if (stateOH)
    stateOH->print_info();

  pic->init_source(*fi);

  if (doRestart) {
    // Restoring the restart data before coupling with GM, because the PIC grid
    // may change again during coupling.
    read_restart();
  }
};

//========================================================
void Domain::update() {
  std::string funcName = "Domain::update";
  timing_func(funcName);

  bool doReport = tc->monitor.is_time_to();

  if (pic->is_grid_empty()) {
    if (tc->get_dt_si() <= 0) {
      tc->set_dt_si(tc->get_dummy_dt_si());
    }
  }

  const Real t0 = tc->get_time_si();
  // update time, step number.
  tc->update();

  if (doReport) {
    const Real t1 = tc->get_time_si();
    Print() << "\n==== " << printPrefix << " Cycle " << tc->get_cycle()
            << " from t = " << std::setprecision(6) << t0
            << " (s) to t = " << std::setprecision(6) << t1
            << " (s) with dt = " << std::setprecision(6) << tc->get_dt_si()
            << " (s) ====" << std::endl;
  }

  pic->update(doReport);

  write_plots();

  pic->write_log();

  pt->update(*pic);
  pt->write_log();
};

//========================================================
void Domain::update_param(const std::string &paramString) {
  param = paramString;
  read_param(false);
  init_time_ctr();
};

//========================================================
void Domain::prepare_grid_info(const amrex::Vector<double> &info) {

  read_param(true);
  param.roll_back();

  // If MHD is 2D, PIC has to be periodic in the z-direction.
  if (isFake2D)
    set_periodicity(iz_, true);

  if (!doRestart && initFromSWMF) {
    // If restart, the grid info will be read from restart.H

    Real si2noL = 1. / info[18];

    int n = 0;
    for (int i = 0; i < nDim; i++) {
      Real phyMin = info[n++] * si2noL; // Lmin
      Real phyMax = phyMin + info[n++] * si2noL;
      Real dx = info[n++] * si2noL; // dx
      nCell[i] = (int)((phyMax - phyMin) / dx + 0.5);

      if (isFake2D && i == iz_) {
        phyMin = 0;
        phyMax = dx;
      }

      domainRange.setLo(i, phyMin);
      domainRange.setHi(i, phyMax);
    }
  }

  for (int i = 0; i < nDim; i++) {
    centerBoxLo[i] = 0;
    centerBoxHi[i] = nCell[i] - 1;
  }

  centerBox.setSmall(centerBoxLo);
  centerBox.setBig(centerBoxHi);

  gm.define(centerBox, &domainRange, coord, periodicity);

  amrInfo.blocking_factor.clear();
  amrInfo.blocking_factor.push_back(IntVect(1, 1, 1));

  Print() << printPrefix << "Domain range = " << domainRange << std::endl;
}

//========================================================
void Domain::regrid() {

  std::string nameFunc = "Domain::regrid";

  // If the PIC grid does not change, then return.
  // If the PIC grid is empty at the beginning, gridInfo.is_grid_new() is false,
  // but it is still required to run the rest of the function to initialize
  // variables. That's why we need isGridInitialized here.
  if (!gridInfo.is_grid_new() && isGridInitialized)
    return;

  Print() << printPrefix << nameFunc << " is called" << std::endl;

  timing_func(nameFunc);

  BoxList bl;
  get_boxlist_from_region(bl, gridInfo, centerBoxLo, centerBoxHi);
  BoxArray activeRegionBA(bl);

  long nCellPic = 0;
  for (const auto &bx : bl) {
    nCellPic += bx.numPts();
  }

  BoxArray baPic(activeRegionBA);

  baPic.maxSize(maxBlockSize);
  Print() << "=====" << printPrefix << " Grid Information summary========="
          << "\n Number of Boxes to describe active region = "
          << activeRegionBA.size()
          << "\n Number of boxes                           = " << baPic.size()
          << "\n Number of active cells                    = " << nCellPic
          << "\n Number of domain cells                    = "
          << centerBox.numPts()
          << "\n Ratio: (active cell)/(Domain cell)        = "
          << nCellPic / centerBox.d_numPts()
          << "\n===================================================="
          << std::endl;

  DistributionMapping dmPic;
  if (!baPic.empty())
    dmPic.define(baPic);

  fi->regrid(baPic, dmPic);

  if (stateOH)
    stateOH->regrid(baPic, dmPic);

  pic->regrid(activeRegionBA, baPic, dmPic);

  pt->regrid(activeRegionBA, baPic, dmPic, *pic);

  iGrid++;
  iDecomp++;

  isGridInitialized = true;
}

//========================================================
void Domain::receive_grid_info(int *status) { gridInfo.set_status(status); }

//========================================================
void Domain::set_ic() {

  // If it is restart, the values should have been restored before coupling with
  // GM. See Domain::init().
  if (doRestart)
    return;

#ifdef _PT_COMPONENT_
  fi->set_node_fluid();
#endif

  pic->fill_new_cells();
  write_plots(true);
  pic->write_log(true, true);

  pt->set_ic(*pic);
  pt->write_log(true, true);
}

//========================================================
void Domain::set_state_var(double *data, int *index) {

  Print() << printPrefix << " GM -> PC coupling at t =" << tc->get_time_si()
          << " (s)" << std::endl;

  if (stateOH) {
    // PT mode
    stateOH->set_node_fluid(data, index);
  } else {
    fi->set_node_fluid(data, index);
    pic->update_cells_for_pt();
  }
}

//========================================================
int Domain::get_grid_nodes_number() { return pic->get_grid_nodes_number(); }

//========================================================
void Domain::get_grid(double *pos_DI) { pic->get_grid(pos_DI); }

//========================================================
void Domain::find_mpi_rank_for_points(const int nPoint,
                                      const double *const xyz_I,
                                      int *const rank_I) {
  pic->find_mpi_rank_for_points(nPoint, xyz_I, rank_I);
}

//========================================================
void Domain::get_fluid_state_for_points(const int nDim, const int nPoint,
                                        const double *const xyz_I,
                                        double *const data_I, const int nVar) {
  pic->get_fluid_state_for_points(nDim, nPoint, xyz_I, data_I, nVar);
}

//========================================================
void Domain::read_restart() {
  std::string restartDir = component + "/restartIN/";

  MultiFab tmp;
  VisMF::Read(tmp, restartDir + gridName + "_centerB");
  BoxArray baPic = tmp.boxArray();
  DistributionMapping dmPic = tmp.DistributionMap();

  pic->regrid(baPic, baPic, dmPic);
  fi->regrid(baPic, dmPic);

  // Assume dmPT == dmPIC so far.
  pt->regrid(baPic, baPic, dmPic, *pic);

  fi->read_restart();

  pic->read_restart();
  write_plots(true);
  pic->write_log(true, true);

  pt->read_restart();
  pt->write_log(true, true);
}

//========================================================
void Domain::save_restart() {
  save_restart_header();
  save_restart_data();
}

//========================================================
void Domain::save_restart_data() {
  fi->save_restart_data();
  pic->save_restart_data();
  pt->save_restart_data();
}

//========================================================
void Domain::save_restart_header() {

  if (ParallelDescriptor::IOProcessor()) {
    Print() << printPrefix
            << "Saving restart file at time = " << tc->get_time_si() << " (s)"
            << std::endl;

    VisMF::IO_Buffer ioBuffer(VisMF::IO_Buffer_Size);

    std::ofstream headerFile;

    headerFile.rdbuf()->pubsetbuf(ioBuffer.dataPtr(), ioBuffer.size());

    std::string headerFileName(component + "/restartOUT/" + gridName +
                               "_restart.H");

    headerFile.open(headerFileName.c_str(),
                    std::ofstream::out | std::ofstream::trunc);

    if (!headerFile.good()) {
      amrex::FileOpenFailed(headerFileName);
    }

    headerFile.precision(17);

    headerFile << std::scientific;

    /*
    The format of the header:

    #COMMAND
    value TAB TAB [TAB] name

    The third TAB is added if the length of value is less than 8 characters.
    Otherwise there are 2 TABS. The “name" should agree with the description in
    the PARAM.XML file.
    */

    headerFile << "Restart header\n\n";

    std::string command_suffix = "_" + gridName + "\n";

    headerFile << "#RESTART" + command_suffix;
    headerFile << (pic->is_grid_empty() ? "F" : "T") << "\t\t\tdoRestart\n";
    headerFile << "\n";

    headerFile << "#NSTEP" + command_suffix;
    headerFile << tc->get_cycle() << "\t\t\tnStep\n";
    headerFile << "\n";

    headerFile << "#TIMESIMULATION" + command_suffix;
    headerFile << tc->get_time_si() << "\t\ttimeSimulation\n";
    headerFile << "\n";

    headerFile << "#TIMESTEP" + command_suffix;
    bool useFixedDt = tc->get_cfl() <= 0;
    headerFile << (useFixedDt ? "T" : "F") << "\t\t\tuseFixedDt\n";
    if (useFixedDt) {
      headerFile << tc->get_dt_si() << "\t\tdt\n";
    } else {
      headerFile << tc->get_cfl() << "\t\tcfl\n";
    }
    headerFile << "\n";

    if (!useFixedDt) {
      headerFile << "#DT" + command_suffix;
      headerFile << tc->get_dt_si() << "\t\tdtSI\n";
      headerFile << tc->get_next_dt_si() << "\t\tdtNextSI\n";
    }

    headerFile << "\n";

    // Geometry
    headerFile << "#GEOMETRY" + command_suffix;
    headerFile << domainRange.lo(ix_) << "\t\txMin\n";
    headerFile << domainRange.hi(ix_) << "\t\txMax\n";
    headerFile << domainRange.lo(iy_) << "\t\tyMin\n";
    headerFile << domainRange.hi(iy_) << "\t\tyMax\n";
    headerFile << domainRange.lo(iz_) << "\t\tzMin\n";
    headerFile << domainRange.hi(iz_) << "\t\tzMax\n";
    headerFile << "\n";

    // Cell
    headerFile << "#NCELL" + command_suffix;
    headerFile << nCell[ix_] << "\t\t\tnCellX\n";
    headerFile << nCell[iy_] << "\t\t\tnCellY\n";
    headerFile << nCell[iz_] << "\t\t\tnCellZ\n";
    headerFile << "\n";

    pic->save_restart_header(headerFile);
    pt->save_restart_header(headerFile);

    headerFile << "\n";
  }
}

//========================================================
void Domain::init_time_ctr() {
  tc->set_si2no(fi->get_Si2NoT());

  { //----------Init plot data------------------------

    //------ Scalar parameters.----------
    std::vector<std::string> scalarName_I;
    std::vector<double> scalarVar_I;
    std::string ms = "mS", qs = "qS";
    const int nS = fi->get_nS();
    for (int i = 0; i < nS; ++i) {
      scalarName_I.push_back(ms + std::to_string(i));
      scalarName_I.push_back(qs + std::to_string(i));
      scalarVar_I.push_back(fi->get_species_mass(i));
      scalarVar_I.push_back(fi->get_species_charge(i) /
                            fi->get_scaling_factor());
    }
    scalarName_I.push_back("cLight");
    scalarVar_I.push_back(fi->get_cLight_SI());
    scalarName_I.push_back("rPlanet");
    scalarVar_I.push_back(fi->get_rPlanet_SI());
    //-------------------------------------

    for (auto &plot : tc->plots) {
      auto &writer = plot.writer;

      // Pass information to writers.
      writer.set_rank(ParallelDescriptor::MyProc());
      writer.set_nProcs(ParallelDescriptor::NProcs());
      writer.set_nDim(fi->get_fluid_dimension());
      writer.set_iRegion(gridID);
      writer.set_domainMin_D({ { domainRange.lo(ix_), domainRange.lo(iy_),
                                 domainRange.lo(iz_) } });

      writer.set_domainMax_D({ { domainRange.hi(ix_), domainRange.hi(iy_),
                                 domainRange.hi(iz_) } });

      const Real *dx = gm.CellSize();
      writer.set_dx_D({ { dx[ix_], dx[iy_], dx[iz_] } });
      writer.set_nSpecies(nS);
      writer.set_units(fi->get_No2SiL(), fi->get_No2SiV(), fi->get_No2SiB(),
                       fi->get_No2SiRho(), fi->get_No2SiP(), fi->get_No2SiJ(),
                       fi->get_rPlanet_SI());
      writer.set_No2NoL(fi->get_MhdNo2NoL());

      writer.set_scalarValue_I(scalarVar_I);
      writer.set_scalarName_I(scalarName_I);
      //--------------------------------------------------
      writer.init();
      // writer.print();
    }
  }
}

//========================================================
void Domain::read_param(const bool readGridInfo) {
  // The default values shoudl be set in the constructor.

  std::string command;

  param.set_verbose(false);
  param.set_command_suffix(gridName);
  param.set_component(component);
  while (param.get_next_command(command)) {

    bool isGridCommand = command == "#MAXBLOCKSIZE" ||
                         command == "#PERIODICITY" || command == "#GEOMETRY" ||
                         command == "#NCELL" || command == "#RESTART" ||
                         command == "#INITFROMSWMF";

    // Skip this command
    if (readGridInfo != isGridCommand)
      continue;

    param.set_verbose(ParallelDescriptor::IOProcessor());
    Print() << "\n"
            << component << ": " << command << " " << gridName << std::endl;

    if (command == "#DIVE" || command == "#EFIELDSOLVER" ||
        command == "#PARTICLES" || command == "#ELECTRON" ||
        command == "#DISCRETIZE" || command == "#DISCRETIZATION" ||
        command == "#RESAMPLING" || command == "#SMOOTHE" ||
        command == "#TESTCASE" || command == "#MERGEPARTICLE" ||
        command == "#SOURCE" || command == "#PIC" ||
        command == "#EXPLICITPIC") {
      pic->read_param(command, param);
    } else if (command == "#PARTICLETRACKER" ||
               command == "#TESTPARTICLENUMBER" || command == "#TPPARTICLES" ||
               command == "#TPCELLINTERVAL" || command == "#TPREGION" ||
               command == "#TPSAVE" || command == "#TPRELATIVISTIC" ||
               command == "#TPINITFROMPIC" || command == "#TPSTATESI") {
      pt->read_param(command, param);
    } else if (command == "#NORMALIZATION" || command == "#SCALINGFACTOR" ||
               command == "#BODYSIZE" || command == "#PLASMA" ||
               command == "#UNIFORMSTATE") {
      fi->read_param(command, param);
    } else if (command == "#RESTART") {
      param.read_var("doRestart", doRestart);
    } else if (command == "#INITFROMSWMF") {
      param.read_var("initFromSWMF", initFromSWMF);
    } else if (command == "#GEOMETRY") {
      for (int i = 0; i < nDim; ++i) {
        Real lo, hi;
        param.read_var("min", lo);
        param.read_var("max", hi);
        domainRange.setLo(i, lo);
        domainRange.setHi(i, hi);
      }
      if (!domainRange.ok())
        Abort("Error: invalid input!");
    } else if (command == "#NCELL") {
      for (int i = 0; i < nDim; ++i) {
        param.read_var("nCell", nCell[i]);
        if (nCell[i] <= 0)
          Abort("Error: invalid input!");
      }
      isFake2D = (nCell[iz_] == 1);
    } else if (command == "#NOUTFILE") {
      param.read_var("nFileField", nFileField);
      param.read_var("nFileParticle", nFileParticle);
    } else if (command == "#PARTICLESTAGGERING") {
      bool doStaggering;
      param.read_var("doStaggering", doStaggering);
      ParticleStaggering ps = doStaggering ? Staggered : NonStaggered;

      Particles<nPicPartReal, 0>::particlePosition = ps;
      Particles<nPTPartReal, nPTPartInt>::particlePosition = ps;

    } else if (command == "#MAXBLOCKSIZE") {
      // The block size in each direction can not larger than maxBlockSize.
      int tmp;
      for (int i = 0; i < nDim; i++) {
        param.read_var("maxBlockSize", tmp);
        maxBlockSize[i] = tmp;
      }
    } else if (command == "#TIMESTEP" || command == "#TIMESTEPPING") {
      bool useFixedDt;
      param.read_var("useFixedDt", useFixedDt);
      if (useFixedDt) {
        Real dtSI;
        param.read_var("dtSI", dtSI);
        tc->set_dt_si(dtSI);
        tc->set_next_dt_si(dtSI);
        tc->set_cfl(-1);
      } else {
        Real cfl;
        param.read_var("cfl", cfl);
        tc->set_cfl(cfl);
      }
    } else if (command == "#PERIODICITY") {
      for (int i = 0; i < nDim; i++) {
        bool isPeriodic;
        param.read_var("isPeriodic", isPeriodic);
        set_periodicity(i, isPeriodic);
      }

    } else if (command == "#SAVELOG") {
      int dn;
      param.read_var("dnSavePic", dn);
      tc->picLog.init(-1, dn);
      param.read_var("dnSavePT", dn);
      tc->ptLog.init(-1, dn);
    } else if (command == "#MONITOR") {
      int dn;
      param.read_var("dnReport", dn);
      tc->monitor.init(-1, dn);
    } else if (command == "#LOADBALANCE") {
      int dn;
      param.read_var("dn", dn);
      Real dt;
      param.read_var("dt", dt);
      tc->loadBalance.init(dt, dn);
    } else if (command == "#SAVEPLOT" || command == "#SAVEIDL") {

      /*
      Example:
      #SAVEPLOT
      6                                 nPlot
      z=0 var real4 planet              plotString
      -1                                dn
      20                                dt
      1                                 dx
      {fluid} numS0                     varName
      y=0 fluid real4 planet            plotString
      -100                              dn
      5                                 dt
      -1                                dx
      3d var amrex planet               plotString
      -1                                dn
      10                                dt
      1                                 dx
      {fluid} numS0                     varName
      3d fluid real4 planet compact     plotString
      -1                                dn
      10                                dt
      1                                 dx
      cut fluid real8 si                plotString
      -1                                dn
      100                               dt
      5                                 xMin
      10                                xMax
      -2                                yMin
      2                                 yMax
      -2                                zMin
      2                                 zMax
      1                                 dx
      3d particles0 amrex planet        plotString
      -1                                dn
      200                               dt
      1                                 dx
      */

      int nPlot;
      param.read_var("nPlotFile", nPlot);

      if (nPlot > 0) {
        tc->plots.clear();
      }

      for (int iPlot = 0; iPlot < nPlot; iPlot++) {

        std::string plotString;
        param.read_var("plotString", plotString);
        {
          std::string::size_type pos = plotString.find_first_not_of(' ');
          if (pos != std::string::npos)
            plotString.erase(0, pos);
        }

        int dnSave;
        param.read_var("dnSavePlot", dnSave);

        Real dtSave;
        param.read_var("dtSavePlot", dtSave);

        std::array<double, nDim> plotMin_D = { 1, 1, 1 },
                                 plotMax_D = { -1, 1 - 1 };
        if (plotString.find("cut") != std::string::npos) {
          // Output range is 'cut' type.
          for (int iDim = 0; iDim < nDim; iDim++) {
            param.read_var("plotMin", plotMin_D[iDim]);
            param.read_var("plotMax", plotMax_D[iDim]);
          }
        }

        int dxSave;
        param.read_var("dxSavePlot", dxSave);

        std::string plotVar;
        if (plotString.find("var") != std::string::npos) {
          param.read_var("plotVar", plotVar);
        }

        PlotCtr pcTmp(tc.get(), iPlot, dtSave, dnSave, plotString, dxSave,
                      plotVar, plotMin_D, plotMax_D);
        tc->plots.push_back(pcTmp);
      }
    } else if (command == "#OHMSLAW") {
      std::string sOhmU;
      Real eta;
      param.read_var("OhmU", sOhmU);
      param.read_var("resistivity", eta);

      fi->set_resistivity(eta);
      fi->set_ohm_u(sOhmU);
      //--------- The commands below exist in restart.H only --------
    } else if (command == "#NSTEP") {
      int nStep;
      param.read_var("nStep", nStep);
      tc->set_cycle(nStep);
    } else if (command == "#TIMESIMULATION") {
      Real time;
      param.read_var("time", time);
      tc->set_time_si(time);
    } else if (command == "#DT") {
      // NOTE: this command is useful only for CFL based time stepping.
      Real dtSI, dtNextSI;
      param.read_var("dtSI", dtSI);
      param.read_var("dtNextSI", dtNextSI);
      tc->set_dt_si(dtSI);
      tc->set_next_dt_si(dtNextSI);
    } else {
      Print() << "Error: command = " << command << std::endl;
      Abort("Can not find this command!");
    }
    //--------- The commands above exist in restart.H only --------
    param.set_verbose(false);
  } // While

  if (!readGridInfo) {
    pic->post_process_param();
    pt->post_process_param();
    fi->post_process_param();
  }

  VisMF::SetNOutFiles(nFileField);

  ParmParse pp("particles");
  pp.add("particles_nfiles", nFileParticle);
}

//========================================================
void Domain::write_plots(bool doForce) { pic->write_plots(doForce); }
