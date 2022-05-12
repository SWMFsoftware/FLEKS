#include "Domain.h"
#include "GridInfo.h"
#include "GridUtility.h"

using namespace amrex;

//========================================================
void Domain::update() {
  std::string funcName = "Domain::update";
  timing_func(funcName);

  bool doReport = tc->monitor.is_time_to();

  if (pic.is_grid_empty()) {
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

  pic.update(doReport);

  write_plots();

  pic.write_log();

  pt.update(pic);
};

//========================================================
void Domain::init(double time, const std::string &paramString, int *paramInt,
                  double *gridDim, double *paramReal, int iDomain) {
  if (AMREX_SPACEDIM != 3)
    Abort("Error: AMReX should be compiled with 3D configuration!!");

  tc->set_time_si(time);

  domainID = iDomain;
  {
    std::stringstream ss;
    ss << "FLEKS" << domainID;
    domainName = ss.str();
    printPrefix = domainName + ": ";
  }

  fluidInterface->init(domainID);
  fluidInterface->receive_info_from_gm(paramInt, gridDim, paramReal,
                                       paramString);

  pic.init(fluidInterface, tc, domainID);

  pt.init(fluidInterface, tc, domainID);

  read_param();

  fluidInterface->print_info();

  make_grid();

  init_time_ctr();

  pic.init_source(*fluidInterface);

  if (doRestart) {
    // Restoring the restart data before coupling with GM, because the PIC grid
    // may change again during coupling.
    read_restart();
  }
};

//========================================================
void Domain::make_grid() {
  nGst = 2;

  // If MHD is 2D, PIC has to be periodic in the z-direction.
  for (int iDim = fluidInterface->get_fluid_dimension(); iDim < nDim; iDim++)
    set_periodicity(iDim, true);

  if (!doRestart) {
    // If restart, these variables read from restart.H
    for (int iDir = ix_; iDir <= iz_; iDir++) {
      nCell[iDir] = fluidInterface->get_phy_cell_number(iDir);
    }

    for (int i = 0; i < nDim; i++) {
      domainRange.setLo(i, fluidInterface->get_phy_domain_min(i));
      domainRange.setHi(i, fluidInterface->get_phy_domain_max(i));
    }
  }

  const int nCellPerPatch = fluidInterface->get_nCellPerPatch();
  gridInfo.init(nCell[ix_], nCell[iy_], nCell[iz_], nCellPerPatch);

  for (int i = 0; i < nDim; i++) {
    centerBoxLo[i] = 0;
    centerBoxHi[i] = nCell[i] - 1;
  }

  centerBox.setSmall(centerBoxLo);
  centerBox.setBig(centerBoxHi);

  geom.define(centerBox, &domainRange, coord, periodicity);

  Print() << printPrefix << "Domain range = " << domainRange << std::endl;

  pic.set_geom(nGst, geom);
  fluidInterface->set_geom(nGst, geom);

  pt.set_geom(nGst, geom);
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
  BoxArray picRegionBA(bl);

  long nCellPic = 0;
  for (const auto &bx : bl) {
    nCellPic += bx.numPts();
  }

  BoxArray baPic(picRegionBA);

  baPic.maxSize(maxBlockSize);
  Print() << "=====" << printPrefix << " Grid Information summary========="
          << "\n Number of Boxes to describe active region = "
          << picRegionBA.size()
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

  fluidInterface->regrid(baPic, dmPic);

  pic.regrid(picRegionBA, baPic, dmPic);

  pt.regrid(picRegionBA, baPic, dmPic, pic);

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

  pic.fill_new_cells();
  write_plots(true);
  pic.write_log(true, true);

  pt.set_ic(pic);
}

//========================================================
void Domain::set_state_var(double *data, int *index) {
  pic.set_state_var(data, index);
}

//========================================================
int Domain::get_grid_nodes_number() { return pic.get_grid_nodes_number(); }

//========================================================
void Domain::get_grid(double *pos_DI) { pic.get_grid(pos_DI); }

//========================================================
void Domain::find_mpi_rank_for_points(const int nPoint,
                                      const double *const xyz_I,
                                      int *const rank_I) {
  pic.find_mpi_rank_for_points(nPoint, xyz_I, rank_I);
}

//========================================================
void Domain::get_fluid_state_for_points(const int nDim, const int nPoint,
                                        const double *const xyz_I,
                                        double *const data_I, const int nVar) {
  pic.get_fluid_state_for_points(nDim, nPoint, xyz_I, data_I, nVar);
}

//========================================================
void Domain::read_restart() {
  std::string restartDir = "PC/restartIN/";

  MultiFab tmp;
  VisMF::Read(tmp, restartDir + domainName + "_centerB");
  BoxArray baPic = tmp.boxArray();
  DistributionMapping dmPic = tmp.DistributionMap();

  pic.regrid(baPic, baPic, dmPic);
  fluidInterface->regrid(baPic, dmPic);

  // Assume dmPT == dmPIC so far.
  pt.regrid(baPic, baPic, dmPic, pic);

  fluidInterface->read_restart();
  pic.read_restart();
  pt.read_restart();

  write_plots(true);
  pic.write_log(true, true);
}

//========================================================
void Domain::save_restart() {
  save_restart_header();
  save_restart_data();
}

//========================================================
void Domain::save_restart_data() {
  VisMF::SetNOutFiles(64);
  fluidInterface->save_restart_data();
  pic.save_restart_data();
  pt.save_restart_data();
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

    std::string headerFileName("PC/restartOUT/" + domainName + "_restart.H");

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

    std::string command_suffix = "_" + domainName + "\n";

    headerFile << "#RESTART" + command_suffix;
    headerFile << (pic.is_grid_empty() ? "F" : "T") << "\t\t\tdoRestart\n";
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

    pic.save_restart_header(headerFile);
    pt.save_restart_header(headerFile);

    headerFile << "\n";
  }
}

//========================================================
void Domain::init_time_ctr() {
  tc->set_si2no(fluidInterface->get_Si2NoT());

  { //----------Init plot data------------------------

    //------ Scalar parameters.----------
    std::vector<std::string> scalarName_I;
    std::vector<double> scalarVar_I;
    std::string ms = "mS", qs = "qS";
    const int nS = fluidInterface->get_nS();
    for (int i = 0; i < nS; ++i) {
      scalarName_I.push_back(ms + std::to_string(i));
      scalarName_I.push_back(qs + std::to_string(i));
      scalarVar_I.push_back(fluidInterface->get_species_mass(i));
      scalarVar_I.push_back(fluidInterface->get_species_charge(i) /
                            fluidInterface->get_scaling_factor());
    }
    scalarName_I.push_back("cLight");
    scalarVar_I.push_back(fluidInterface->get_cLight_SI());
    scalarName_I.push_back("rPlanet");
    scalarVar_I.push_back(fluidInterface->get_rPlanet_SI());
    //-------------------------------------

    for (auto &plot : tc->plots) {
      auto &writer = plot.writer;

      // Pass information to writers.
      writer.set_rank(ParallelDescriptor::MyProc());
      writer.set_nProcs(ParallelDescriptor::NProcs());
      writer.set_nDim(fluidInterface->get_fluid_dimension());
      writer.set_iRegion(domainID);
      writer.set_domainMin_D({ { domainRange.lo(ix_), domainRange.lo(iy_),
                                 domainRange.lo(iz_) } });

      writer.set_domainMax_D({ { domainRange.hi(ix_), domainRange.hi(iy_),
                                 domainRange.hi(iz_) } });

      const Real *dx = geom.CellSize();
      writer.set_dx_D({ { dx[ix_], dx[iy_], dx[iz_] } });
      writer.set_nSpecies(nS);
      writer.set_units(
          fluidInterface->get_No2SiL(), fluidInterface->get_No2SiV(),
          fluidInterface->get_No2SiB(), fluidInterface->get_No2SiRho(),
          fluidInterface->get_No2SiP(), fluidInterface->get_No2SiJ(),
          fluidInterface->get_rPlanet_SI());
      writer.set_No2NoL(fluidInterface->get_MhdNo2NoL());

      writer.set_scalarValue_I(scalarVar_I);
      writer.set_scalarName_I(scalarName_I);
      //--------------------------------------------------
      writer.init();
      // writer.print();
    }
  }
}

//========================================================
void Domain::read_param() {
  // The default values shoudl be set in the constructor.

  std::string command;
  ReadParam &readParam = fluidInterface->readParam;
  readParam.set_verbose(ParallelDescriptor::MyProc() == 0);

  readParam.set_command_suffix(domainName);

  while (readParam.get_next_command(command)) {

    if (command == "#DIVE" || command == "#EFIELDSOLVER" ||
        command == "#PARTICLES" || command == "#ELECTRON" ||
        command == "#DISCRETIZE" || command == "#DISCRETIZATION" ||
        command == "#RESAMPLING" || command == "#SMOOTHE" ||
        command == "#TESTCASE" || command == "#MERGEPARTICLE" ||
        command == "#SOURCE" || command == "#PIC" ||
        command == "#EXPLICITPIC") {
      pic.read_param(command, readParam);
    } else if (command == "#PARTICLETRACKER" ||
               command == "#TESTPARTICLENUMBER" || command == "#TPPARTICLES" ||
               command == "#TPCELLINTERVAL" || command == "#TPREGION" ||
               command == "#TPSAVE" || command == "#TPRELATIVISTIC" ||
               command == "#TPINITFROMPIC" || command == "#TPSTATESI") {
      pt.read_param(command, readParam);
    } else if (command == "#OHMSLAW") {
      std::string sOhmU;
      Real eta;
      readParam.read_var("OhmU", sOhmU);
      readParam.read_var("resistivity", eta);

      fluidInterface->set_resistivity(eta);
      fluidInterface->set_ohm_u(sOhmU);
    } else if (command == "#RESTART") {
      readParam.read_var("doRestart", doRestart);
      pic.set_doRestart(doRestart);
    } else if (command == "#PARTICLESTAGGERING") {
      bool doStaggering;
      readParam.read_var("doStaggering", doStaggering);
      ParticleStaggering ps = doStaggering ? Staggered : NonStaggered;

      Particles<nPicPartReal, 0>::particlePosition = ps;
      Particles<nPTPartReal, nPTPartInt>::particlePosition = ps;

    } else if (command == "#MAXBLOCKSIZE") {
      // The block size in each direction can not larger than maxBlockSize.
      int tmp;
      for (int i = 0; i < nDim; i++) {
        readParam.read_var("maxBlockSize", tmp);
        maxBlockSize[i] = tmp;
      }
    } else if (command == "#TIMESTEP" || command == "#TIMESTEPPING") {
      bool useFixedDt;
      readParam.read_var("useFixedDt", useFixedDt);
      if (useFixedDt) {
        Real dtSI;
        readParam.read_var("dtSI", dtSI);
        tc->set_dt_si(dtSI);
        tc->set_next_dt_si(dtSI);
        tc->set_cfl(-1);
      } else {
        Real cfl;
        readParam.read_var("cfl", cfl);
        tc->set_cfl(cfl);
      }
    } else if (command == "#PERIODICITY") {
      for (int i = 0; i < nDim; i++) {
        bool isPeriodic;
        readParam.read_var("isPeriodic", isPeriodic);
        set_periodicity(i, isPeriodic);
      }

    } else if (command == "#SAVELOG") {
      int dn;
      readParam.read_var("dnSave", dn);
      tc->log.init(-1, dn);
    } else if (command == "#MONITOR") {
      int dn;
      readParam.read_var("dnReport", dn);
      tc->monitor.init(-1, dn);
    } else if (command == "#LOADBALANCE") {
      int dn;
      readParam.read_var("dn", dn);
      Real dt;
      readParam.read_var("dt", dt);
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
      readParam.read_var("nPlotFile", nPlot);

      for (int iPlot = 0; iPlot < nPlot; iPlot++) {

        std::string plotString;
        readParam.read_var("plotString", plotString);
        {
          std::string::size_type pos = plotString.find_first_not_of(' ');
          if (pos != std::string::npos)
            plotString.erase(0, pos);
        }

        int dnSave;
        readParam.read_var("dnSavePlot", dnSave);

        Real dtSave;
        readParam.read_var("dtSavePlot", dtSave);

        std::array<double, nDim> plotMin_D = { 1, 1, 1 },
                                 plotMax_D = { -1, 1 - 1 };
        if (plotString.find("cut") != std::string::npos) {
          // Output range is 'cut' type.
          for (int iDim = 0; iDim < nDim; iDim++) {
            readParam.read_var("plotMin", plotMin_D[iDim]);
            readParam.read_var("plotMax", plotMax_D[iDim]);
          }
        }

        int dxSave;
        readParam.read_var("dxSavePlot", dxSave);

        std::string plotVar;
        if (plotString.find("var") != std::string::npos) {
          readParam.read_var("plotVar", plotVar);
        }

        PlotCtr pcTmp(tc.get(), iPlot, dtSave, dnSave, plotString, dxSave,
                      plotVar, plotMin_D, plotMax_D);
        tc->plots.push_back(pcTmp);
      }
      //--------- The commands below exist in restart.H only --------
    } else if (command == "#NSTEP") {
      int nStep;
      readParam.read_var("nStep", nStep);
      tc->set_cycle(nStep);
    } else if (command == "#TIMESIMULATION") {
      Real time;
      readParam.read_var("time", time);
      tc->set_time_si(time);
    } else if (command == "#DT") {
      // NOTE: this command is useful only for CFL based time stepping.
      Real dtSI, dtNextSI;
      readParam.read_var("dtSI", dtSI);
      readParam.read_var("dtNextSI", dtNextSI);
      tc->set_dt_si(dtSI);
      tc->set_next_dt_si(dtNextSI);
    } else if (command == "#GEOMETRY") {
      for (int i = 0; i < nDim; ++i) {
        Real lo, hi;
        readParam.read_var("min", lo);
        readParam.read_var("max", hi);
        domainRange.setLo(i, lo);
        domainRange.setHi(i, hi);
      }
    } else if (command == "#NCELL") {
      for (int i = 0; i < nDim; ++i) {
        readParam.read_var("nCell", nCell[i]);
      }
    } else {
      Print() << "Error: command = " << command << std::endl;
      Abort("Can not find this command!");
    }
    //--------- The commands above exist in restart.H only --------
  }

  pic.post_process_param();
  pt.post_process_param();
}

//========================================================
void Domain::write_plots(bool doForce) { pic.write_plots(doForce); }
