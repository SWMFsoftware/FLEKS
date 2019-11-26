#include "Domain.h"
#include "SWMFDomains.h"

using namespace amrex;

void Domain::set_state_var(double* data, int* index) {
  std::string nameFunc = "Domain::set_state_var";
  fluidInterface.set_couple_node_value(data, index);
  return;
}

int Domain::get_grid_nodes_number() {
  return fluidInterface.count_couple_node_number();
}

void Domain::get_grid(double* pos_DI) {
  std::string nameFunc = "Domain::get_grid";
  fluidInterface.get_couple_node_loc(pos_DI);
  return;
}

void Domain::find_mpi_rank_for_points(const int nPoint,
                                      const double* const xyz_I,
                                      int* const rank_I) {
  int nDimGM = fluidInterface.getnDim();
  amrex::Real si2nol = fluidInterface.getSi2NoL();
  for (int i = 0; i < nPoint; i++) {
    amrex::Real x = xyz_I[i * nDimGM + ix_] * si2nol;
    amrex::Real y = xyz_I[i * nDimGM + iy_] * si2nol;
    amrex::Real z = 0;
    if (nDimGM > 2)
      z = xyz_I[i * nDimGM + iz_] * si2nol;
    rank_I[i] = find_mpi_rank_from_coord(x, y, z);
  }
}

void Domain::get_fluid_state_for_points(const int nDim, const int nPoint,
                                        const double* const xyz_I,
                                        double* const data_I, const int nVar) {
  // (rho + 3*Moment + 6*p)*nSpecies+ 3*E + 3*B;
  const int nVarPerSpecies = 10;
  int nVarPIC = nSpecies * nVarPerSpecies + 6;
  double dataPIC_I[nVarPIC];

  const int iBx_ = nSpecies * nVarPerSpecies, iBy_ = iBx_ + 1, iBz_ = iBy_ + 1;
  const int iEx_ = iBz_ + 1, iEy_ = iEx_ + 1, iEz_ = iEy_ + 1;

  for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    double pic_D[3] = { 0 };
    for (int iDim = 0; iDim < nDim; iDim++) {
      pic_D[iDim] = xyz_I[iPoint * nDim + iDim] * fluidInterface.getSi2NoL();
    }

    const Real xp = pic_D[0];
    const Real yp = (nDim > 1) ? pic_D[1] : 0.0;
    const Real zp = (nDim > 2) ? pic_D[2] : 0.0;

    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      for (int iVar = iRho_; iVar <= iPyz_; iVar++) {
        const int iStart = iSpecies * nVarPerSpecies;
        dataPIC_I[iStart + iVar] =
            get_value_at_loc(nodePlasma[iSpecies], geom, xp, yp, zp, iVar);
      }

    for (int iDir = ix_; iDir <= iz_; iDir++) {
      dataPIC_I[iBx_ + iDir] = get_value_at_loc(nodeB, geom, xp, yp, zp, iDir);
    }

    for (int iDir = ix_; iDir <= iz_; iDir++) {
      dataPIC_I[iEx_ + iDir] = get_value_at_loc(nodeE, geom, xp, yp, zp, iDir);
    }

    // Combine PIC plasma data into MHD fluid data.
    fluidInterface.CalcFluidState(dataPIC_I, &data_I[iPoint * nVar]);
  }
}

void Domain::find_output_list(const PlotWriter& writerIn,
                              long int& nPointAllProc,
                              PlotWriter::VectorPointList& pointList_II,
                              std::array<double, nDimMax>& xMin_D,
                              std::array<double, nDimMax>& xMax_D) {
  const auto plo = geom.ProbLo();
  const auto plh = geom.ProbHi();

  Real xMinL_D[nDimMax] = { plh[ix_], plh[iy_], plh[iz_] };
  Real xMaxL_D[nDimMax] = { plo[ix_], plo[iy_], plo[iz_] };

  const auto dx = geom.CellSize();

  int iBlock = 0;
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.validbox();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    int iMax = hi.x, jMax = hi.y, kMax = hi.z;

    const bool isNode = true;
    if (isNode) {
      // Avoid double counting the shared edges.
      iMax--;
      jMax--;
      kMax--;

      if ((!geom.isPeriodic(ix_)) && box.bigEnd(ix_) == hi.x)
        iMax++;
      if ((!geom.isPeriodic(iy_)) && box.bigEnd(iy_) == hi.y)
        jMax++;
      if ((!geom.isPeriodic(iz_)) && box.bigEnd(iz_) == hi.z)
        kMax++;
    }

    for (int k = lo.z; k <= kMax; ++k) {
      const double zp = k * dx[iz_] + plo[iz_];
      for (int j = lo.y; j <= jMax; ++j) {
        const double yp = j * dx[iy_] + plo[iy_];
        for (int i = lo.x; i <= iMax; ++i) {
          const double xp = i * dx[ix_] + plo[ix_];

          if (writerIn.is_inside_plot_region(i, j, k, xp, yp, zp)) {
            pointList_II.push_back({ (double)i, (double)j, (double)k, xp, yp,
                                     zp, (double)iBlock });
            if (xp < xMinL_D[ix_])
              xMinL_D[ix_] = xp;
            if (yp < xMinL_D[iy_])
              xMinL_D[iy_] = yp;
            if (zp < xMinL_D[iz_])
              xMinL_D[iz_] = zp;

            if (xp > xMaxL_D[ix_])
              xMaxL_D[ix_] = xp;
            if (yp > xMaxL_D[iy_])
              xMaxL_D[iy_] = yp;
            if (zp > xMaxL_D[iz_])
              xMaxL_D[iz_] = zp;
          }
        }
      }
    }
    iBlock++;
  }

  nPointAllProc = pointList_II.size();
  ParallelDescriptor::ReduceLongSum(nPointAllProc);
  

  ParallelDescriptor::ReduceRealMin(xMinL_D, nDimMax);
  ParallelDescriptor::ReduceRealMax(xMaxL_D, nDimMax);

  for (int iDim = 0; iDim < nDimMax; ++iDim) {
    xMin_D[iDim] = xMinL_D[iDim];
    xMax_D[iDim] = xMaxL_D[iDim];
  }
}

void Domain::get_field_var(const VectorPointList& pointList_II,
                           const std::vector<std::string>& sVar_I,
                           MDArray<double>& var_II) {
  const int iBlk_ = 6;

  long nPoint = pointList_II.size();
  int nVar = sVar_I.size();

  int iBlockCount = 0;
  long iPoint = 0;
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.validbox();

    while (iPoint < nPoint) {
      const int ix = pointList_II[iPoint][ix_];
      const int iy = pointList_II[iPoint][iy_];
      const int iz = pointList_II[iPoint][iz_];
      const int iBlock = pointList_II[iPoint][iBlk_];

      if (iBlock == iBlockCount) {

        for (int iVar = 0; iVar < nVar; ++iVar) {
          var_II(iPoint, iVar) = get_var(sVar_I[iVar], ix, iy, iz, mfi);
        }
        iPoint++;
      } else {
        break;
      }
    }

    iBlockCount++;
  }
}

double Domain::get_var(std::string var, const int ix, const int iy,
                       const int iz, const MFIter& mfi) {

  auto get_is = [var]() {
    std::string::size_type pos;
    std::stringstream ss;
    int is;
    pos = var.find_first_of("0123456789");
    ss << var.substr(pos);
    ss >> is;
    return is;
  };

  double value;
  if (var.substr(0, 1) == "X") {
    const auto plo = geom.ProbLo();
    const auto dx = geom.CellSize();
    value = ix * dx[ix_] + plo[ix_];
  } else if (var.substr(0, 1) == "Y") {
    const auto plo = geom.ProbLo();
    const auto dx = geom.CellSize();
    value = iy * dx[iy_] + plo[iy_];
  } else if (var.substr(0, 1) == "Z") {
    const auto plo = geom.ProbLo();
    const auto dx = geom.CellSize();
    value = iz * dx[iz_] + plo[iz_];
  } else if (var.substr(0, 2) == "Ex") {
    const amrex::Array4<amrex::Real const>& arr = nodeE[mfi].array();
    value = arr(ix, iy, iz, ix_);
  } else if (var.substr(0, 2) == "Ey") {
    const amrex::Array4<amrex::Real const>& arr = nodeE[mfi].array();
    value = arr(ix, iy, iz, iy_);
  } else if (var.substr(0, 2) == "Ez") {
    const amrex::Array4<amrex::Real const>& arr = nodeE[mfi].array();
    value = arr(ix, iy, iz, iz_);
  } else if (var.substr(0, 2) == "Bx") {
    const amrex::Array4<amrex::Real const>& arr = nodeB[mfi].array();
    value = arr(ix, iy, iz, ix_);
  } else if (var.substr(0, 2) == "By") {
    const amrex::Array4<amrex::Real const>& arr = nodeB[mfi].array();
    value = arr(ix, iy, iz, iy_);
  } else if (var.substr(0, 2) == "Bz") {
    const amrex::Array4<amrex::Real const>& arr = nodeB[mfi].array();
    value = arr(ix, iy, iz, iz_);
  } else if (var.substr(0, 4) == "rhoS" || var.substr(0, 3) == "uxS" ||
             var.substr(0, 3) == "uyS" || var.substr(0, 3) == "uzS" ||
             var.substr(0, 4) == "pXXS" || var.substr(0, 4) == "pYYS" ||
             var.substr(0, 4) == "pZZS" || var.substr(0, 4) == "pXYS" ||
             var.substr(0, 4) == "pXZS" || var.substr(0, 4) == "pYZS" ||
             var.substr(0, 4) == "numS") {
    int iVar;

    if (var.substr(0, 4) == "rhoS")
      iVar = iRho_;
    if (var.substr(0, 3) == "uxS")
      iVar = iUx_;
    if (var.substr(0, 3) == "uyS")
      iVar = iUy_;
    if (var.substr(0, 3) == "uzS")
      iVar = iUz_;
    if (var.substr(0, 4) == "pXXS")
      iVar = iPxx_;
    if (var.substr(0, 4) == "pYYS")
      iVar = iPyy_;
    if (var.substr(0, 4) == "pZZS")
      iVar = iPzz_;
    if (var.substr(0, 4) == "pXYS")
      iVar = iPxy_;
    if (var.substr(0, 4) == "pXZS")
      iVar = iPxz_;
    if (var.substr(0, 4) == "pYZS")
      iVar = iPyz_;
    if (var.substr(0, 4) == "numS")
      iVar = iNum_;

    const amrex::Array4<amrex::Real const>& arr =
        nodePlasma[get_is()][mfi].array();
    value = arr(ix, iy, iz, iVar);
  } else if (var.substr(0, 2) == "pS") {
    const amrex::Array4<amrex::Real const>& arr =
        nodePlasma[get_is()][mfi].array();
    value = (arr(ix, iy, iz, iPxx_) + arr(ix, iy, iz, iPyy_) +
             arr(ix, iy, iz, iPzz_)) /
            3.0;
  } else {
    value = 0;
  }

  return value;
}

void Domain::save_restart() {
  save_restart_header();
  save_restart_data();
}

void Domain::save_restart_data() {
  VisMF::SetNOutFiles(64);

  std::string restartDir = "PC/restartOUT/";
  VisMF::Write(nodeE, restartDir + "nodeE");
  VisMF::Write(nodeB, restartDir + "nodeB");
  VisMF::Write(centerB, restartDir + "centerB");

  for (int iPart = 0; iPart < parts.size(); iPart++) {
    parts[iPart]->Checkpoint(restartDir, "particles" + std::to_string(iPart));
  }
}

void Domain::save_restart_header() {
  if (ParallelDescriptor::IOProcessor()) {
    Print() << "Saving restart file at time = " << tc.get_time_si() << " (s)"
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
    headerFile << tc.get_cycle() << "\t nStep\n";
    headerFile << "\n";

    headerFile << "#TIMESIMULATION\n";
    headerFile << tc.get_time_si() << "\t TimeSimulation\n";
    headerFile << "\n";

    headerFile << "#TIMESTEP\n";
    headerFile << tc.get_dt_si() << "\t dt\n";
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


    headerFile << "#ELECTRON\n";
    headerFile << qomEl << "\t qomEl\n";
    headerFile << "\n";

    headerFile << "#PARTICLES\n";
    for (int i = 0; i < nDimMax; ++i) {
      headerFile << nPartPerCell[i] << "\n";
    }
    headerFile << "\n";

  }
}

void Domain::read_restart() {

  std::string restartDir = "PC/restartIN/";
  VisMF::Read(nodeE, restartDir + "nodeE");
  VisMF::Read(nodeB, restartDir + "nodeB");
  VisMF::Read(centerB, restartDir + "centerB");

  nodeE.FillBoundary(geom.periodicity());
  nodeB.FillBoundary(geom.periodicity());
  centerB.FillBoundary(geom.periodicity());

  for (int iPart = 0; iPart < parts.size(); iPart++) {
    parts[iPart]->Restart(restartDir, "particles" + std::to_string(iPart));
  }
  sum_moments();
  sum_to_center(false);
}

void Domain::write_log(bool doForce, bool doCreateFile) {
  if (doCreateFile && ParallelDescriptor::IOProcessor()) {
    std::stringstream ss;
    int time = tc.get_time_si(); // double to int.
    ss << "PC/plots/log_n" << std::setfill('0') << std::setw(8)
       << tc.get_cycle() << ".log";
    logFile = ss.str();
    std::ofstream of(logFile.c_str());
    of << "time nStep Etot Ee Eb Epart ";
    for (int i = 0; i < nSpecies; i++)
      of << " Epart" << i << std::endl;
    of.close();
  }

  if (tc.log.is_time_to(doForce)) {
    Real eEnergy = calc_E_field_energy();
    Real bEnergy = calc_B_field_energy();
    if (ParallelDescriptor::IOProcessor()) {
      std::ofstream of(logFile.c_str(), std::fstream::app);
      of.precision(12);
      of << tc.get_time_si() << "\t" << tc.get_cycle() << "\t"
         << "\t" << (eEnergy + bEnergy + plasmaEnergy[iTot]) << "\t" << eEnergy
         << "\t" << bEnergy << "\t" << plasmaEnergy[iTot];
      for (int i = 0; i < nSpecies; i++)
        of << "\t" << plasmaEnergy[i];
      of << std::endl;
      of.close();
    }
  }
}

void find_output_list_caller(const PlotWriter& writerIn,
                             long int& nPointAllProc,
                             PlotWriter::VectorPointList& pointList_II,
                             std::array<double, nDimMax>& xMin_D,
                             std::array<double, nDimMax>& xMax_D) {
  MPICs->find_output_list(writerIn, nPointAllProc, pointList_II, xMin_D,
                          xMax_D);
}

void get_field_var_caller(const VectorPointList& pointList_II,
                          const std::vector<std::string>& sVar_I,
                          MDArray<double>& var_II) {
  MPICs->get_field_var(pointList_II, sVar_I, var_II);
}
