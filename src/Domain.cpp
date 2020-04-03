#include "Domain.h"
#include "GridInfo.h"
#include "GridUtility.h"

using namespace amrex;

//========================================================
void Domain::update() {
  pic.update();

  write_plots();

  pic.write_log();
};

//========================================================
void Domain::init(amrex::Real timeIn, const std::string &paramString,
                  int *paramInt, double *gridDim, double *paramReal,
                  int iDomain) {

  tc->set_time_si(timeIn);

  fluidInterface->init();
  fluidInterface->receive_info_from_gm(paramInt, gridDim, paramReal,
                                       paramString);

  pic.init(timeIn, paramString, paramInt, gridDim, paramReal, fluidInterface,
           tc);

  read_param();

  fluidInterface->PrintFluidPicInterface();

  make_grid();

  init_time_ctr();
};

//========================================================
void Domain::make_grid() {
  set_nGst(2);

  // If MHD is 2D, PIC has to be periodic in the z-direction.
  for (int iDim = fluidInterface->getnDim(); iDim < nDimMax; iDim++)
    set_periodicity(iDim, true);

  if (!doRestart) {
    IntVect nCellTmp;
    RealBox boxRangeTmp;

    // If restart, these variables read from restart.H
    nCellTmp[ix_] = fluidInterface->getFluidNxc();
    nCellTmp[iy_] = fluidInterface->getFluidNyc();
    nCellTmp[iz_] = fluidInterface->getFluidNzc();

    for (int i = 0; i < nDim; i++) {
      boxRangeTmp.setLo(i, fluidInterface->getphyMin(i));
      boxRangeTmp.setHi(i, fluidInterface->getphyMax(i));
    }
    set_nCell(nCellTmp);
    set_boxRange(boxRangeTmp);
  }

  const int nCellPerPatch = fluidInterface->get_nCellPerPatch();
  gridInfo.init(nCell[ix_], nCell[iy_], nCell[iz_], nCellPerPatch);

  DomainGrid::init();

  pic.set_geom(nGst, geom);
  fluidInterface->set_geom(nGst, geom);
}

//========================================================
void Domain::regrid() {

  std::string nameFunc = "Domain::regrid";
  Print() << nameFunc << " is runing..." << std::endl;

  iGrid++;
  iDecomp++;

  BoxList bl;
  get_boxlist_from_region(bl, gridInfo, centerBoxLo, centerBoxHi);
  BoxArray picRegionBA(bl);

  BoxArray baPic(picRegionBA); //= resize_pic_ba(tc->get_cycle());
  baPic.maxSize(maxBlockSize);
  Print() << "Box # to describe PIC region = " << picRegionBA.size() << "\n"
          << "Total PIC box # = " << baPic.size() << std::endl;

  DistributionMapping dmPic(baPic);
  pic.regrid(picRegionBA, baPic, dmPic);
  fluidInterface->regrid(baPic, dmPic);

  if (doRestart && isInitializing) {
    // Restoring the restart data immediately after creating the data
    // structures, and injecting particles at boundaries. The field boundary
    // nodes will be overwritten later when GM pass information to PC.
    read_restart();
    isInitializing = false;
  }
}

//========================================================
void Domain::receive_grid_info(int *status) { gridInfo.set_status(status); }

//========================================================
void Domain::set_ic() {

  // If it is restart, the values should be read in immediately after creating
  // the data objects. See Domain::regrid().
  if (doRestart)
    return;

  pic.fill_new_cells();
  write_plots(true);
  pic.write_log(true, true);

  isInitializing = false;
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
  fluidInterface->read_restart();
  pic.read_restart();

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
}

//========================================================
void Domain::save_restart_header() {
  if (ParallelDescriptor::IOProcessor()) {
    Print() << "Saving restart file at time = " << tc->get_time_si() << " (s)"
            << std::endl;

    VisMF::IO_Buffer ioBuffer(VisMF::IO_Buffer_Size);

    std::ofstream headerFile;

    headerFile.rdbuf()->pubsetbuf(ioBuffer.dataPtr(), ioBuffer.size());

    std::string headerFileName("PC/restartOUT/restart.H");

    headerFile.open(headerFileName.c_str(),
                    std::ofstream::out | std::ofstream::trunc);

    if (!headerFile.good()) {
      amrex::FileOpenFailed(headerFileName);
    }

    headerFile.precision(17);

    headerFile << "Restart header \n\n";

    headerFile << "#RESTART\n";
    headerFile << "T"
               << "\t doRestart\n";
    headerFile << "\n";

    headerFile << "#NSTEP\n";
    headerFile << tc->get_cycle() << "\t nStep\n";
    headerFile << "\n";

    headerFile << "#TIMESIMULATION\n";
    headerFile << tc->get_time_si() << "\t TimeSimulation\n";
    headerFile << "\n";

    headerFile << "#TIMESTEP\n";
    headerFile << tc->get_dt_si() << "\t dt\n";
    headerFile << "\n";

    // Geometry
    headerFile << "#GEOMETRY\n";
    for (int i = 0; i < nDimMax; ++i) {
      headerFile << boxRange.lo(i) << "\t min\n";
      headerFile << boxRange.hi(i) << "\t max\n";
    }
    headerFile << "\n";

    // Cell
    headerFile << "#NCELL\n";
    for (int i = 0; i < nDimMax; ++i) {
      headerFile << nCell[i] << "\n";
    }
    headerFile << "\n";

    pic.save_restart_header(headerFile);

    headerFile << "\n";
  }
}

//========================================================
void Domain::init_time_ctr() {
  tc->set_si2no(fluidInterface->getSi2NoT());

  { //----------Init plot data------------------------

    //------ Scalar parameters.----------
    std::vector<std::string> scalarName_I;
    std::vector<double> scalarVar_I;
    std::string ms = "mS", qs = "qS";
    const int nS = fluidInterface->get_nS();
    for (int i = 0; i < nS; ++i) {
      scalarName_I.push_back(ms + std::to_string(i));
      scalarName_I.push_back(qs + std::to_string(i));
      scalarVar_I.push_back(fluidInterface->getMiSpecies(i));
      scalarVar_I.push_back(fluidInterface->getQiSpecies(i));
    }
    scalarName_I.push_back("cLight");
    scalarVar_I.push_back(fluidInterface->getcLightSI());
    scalarName_I.push_back("rPlanet");
    scalarVar_I.push_back(fluidInterface->getrPlanet());
    //-------------------------------------

    for (auto &plot : tc->plots) {
      auto &writer = plot.writer;

      // Pass information to writers.
      writer.set_rank(ParallelDescriptor::MyProc());
      writer.set_nProcs(ParallelDescriptor::NProcs());
      writer.set_nDim(fluidInterface->getnDim());
      writer.set_iRegion(0);
      writer.set_domainMin_D(
          { { boxRange.lo(ix_), boxRange.lo(iy_), boxRange.lo(iz_) } });

      writer.set_domainMax_D(
          { { boxRange.hi(ix_), boxRange.hi(iy_), boxRange.hi(iz_) } });

      const Real *dx = geom.CellSize();
      writer.set_dx_D({ { dx[ix_], dx[iy_], dx[iz_] } });

      writer.set_axisOrigin_D({ { 0, 0, 0 } });
      writer.set_nSpecies(nS);
      writer.set_units(fluidInterface->getNo2SiL(), fluidInterface->getNo2SiV(),
                       fluidInterface->getNo2SiB(),
                       fluidInterface->getNo2SiRho(),
                       fluidInterface->getNo2SiP(), fluidInterface->getNo2SiJ(),
                       fluidInterface->getrPlanet());
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
  while (readParam.get_next_command(command)) {

    if (command == "#DIVE" || command == "#EFIELDSOLVER" ||
        command == "#PARTICLES" || command == "#ELECTRON" ||
        command == "#DISCRETIZE" || command == "#DISCRETIZATION" ||
        command == "#RESAMPLING") {
      pic.read_param(command, readParam);
    } else if (command == "#RESTART") {
      readParam.read_var("doRestart", doRestart);
      pic.set_doRestart(doRestart);
    } else if (command == "#MAXBLOCKSIZE") {
      // The block size in each direction can not larger than maxBlockSize.
      int tmp;
      for (int i = 0; i < nDimMax; i++) {
        readParam.read_var("maxBlockSize", tmp);
        set_maxBlockSize(i, tmp);
      }
    } else if (command == "#TIMESTEP" || command == "#TIMESTEPPING") {
      Real dtSI;
      readParam.read_var("dtSI", dtSI);
      tc->set_dt_si(dtSI);
    } else if (command == "#PERIODICITY") {
      for (int i = 0; i < nDimMax; i++) {
        bool isPeriodic;
        readParam.read_var("isPeriodic", isPeriodic);
        set_periodicity(i, isPeriodic);
      }

    } else if (command == "#SAVELOG") {
      int dn;
      readParam.read_var("dnSave", dn);
      tc->log.init(-1, dn);
    } else if (command == "#LOADBALANCE") {
      int dn;
      readParam.read_var("dn", dn);
      Real dt;
      readParam.read_var("dt", dt);
      tc->loadBalance.init(dt, dn);
    } else if (command == "#SAVEPLOT" || command == "#SAVEIDL") {
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

        int dxSave;
        readParam.read_var("dxSavePlot", dxSave);

        std::string plotVar;
        if (plotString.find("var") != std::string::npos) {
          readParam.read_var("plotVar", plotVar);
        }

        PlotCtr pcTmp(tc.get(), iPlot, dtSave, dnSave, plotString, dxSave,
                      plotVar);
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
    } else if (command == "#GEOMETRY") {
      RealBox bxTmp;
      for (int i = 0; i < nDimMax; ++i) {
        Real lo, hi;
        readParam.read_var("min", lo);
        readParam.read_var("max", hi);
        bxTmp.setLo(i, lo);
        bxTmp.setHi(i, hi);
      }
      set_boxRange(bxTmp);
    } else if (command == "#NCELL") {
      IntVect tmp;
      for (int i = 0; i < nDimMax; ++i) {
        readParam.read_var("nCell", tmp[i]);
      }
      set_nCell(tmp);
    }
    //--------- The commands above exist in restart.H only --------
  }
}

//========================================================
void Domain::write_plots(bool doForce) { pic.write_plots(doForce); }