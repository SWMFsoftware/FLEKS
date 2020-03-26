#include <AMReX_PlotFileUtil.H>

#include "Pic.h"
#include "SWMFDomains.h"

using namespace amrex;

void Pic::set_state_var(double* data, int* index) {
  std::string nameFunc = "Pic::set_state_var";
  Print() << nameFunc << " begin" << std::endl;
  fluidInterface->set_couple_node_value(data, index);

  if (doNeedFillNewCell)
    fill_new_cells();

  Print() << nameFunc << " end" << std::endl;
  return;
}

int Pic::get_grid_nodes_number() {
  return fluidInterface->count_couple_node_number();
}

void Pic::get_grid(double* pos_DI) {
  std::string nameFunc = "Pic::get_grid";
  fluidInterface->get_couple_node_loc(pos_DI);
  return;
}

void Pic::find_mpi_rank_for_points(const int nPoint, const double* const xyz_I,
                                   int* const rank_I) {
  int nDimGM = fluidInterface->getnDim();
  amrex::Real si2nol = fluidInterface->getSi2NoL();
  for (int i = 0; i < nPoint; i++) {
    amrex::Real x = xyz_I[i * nDimGM + ix_] * si2nol;
    amrex::Real y = xyz_I[i * nDimGM + iy_] * si2nol;
    amrex::Real z = 0;
    if (nDimGM > 2)
      z = xyz_I[i * nDimGM + iz_] * si2nol;
    rank_I[i] = find_mpi_rank_from_coord(x, y, z);
  }
}

void Pic::get_fluid_state_for_points(const int nDim, const int nPoint,
                                     const double* const xyz_I,
                                     double* const data_I, const int nVar) {
  std::string nameFunc = "Pic::get_fluid_state_for_points";
  Print() << nameFunc << " begin" << std::endl;

  // (rho + 3*Moment + 6*p)*nSpecies+ 3*E + 3*B;
  const int nVarPerSpecies = 10;
  int nVarPIC = nSpecies * nVarPerSpecies + 6;
  double dataPIC_I[nVarPIC];

  const int iBx_ = nSpecies * nVarPerSpecies, iBy_ = iBx_ + 1, iBz_ = iBy_ + 1;
  const int iEx_ = iBz_ + 1, iEy_ = iEx_ + 1, iEz_ = iEy_ + 1;

  for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    double pic_D[3] = { 0 };
    for (int iDim = 0; iDim < nDim; iDim++) {
      pic_D[iDim] = xyz_I[iPoint * nDim + iDim] * fluidInterface->getSi2NoL();
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
    fluidInterface->CalcFluidState(dataPIC_I, &data_I[iPoint * nVar]);
  }
  Print() << nameFunc << " end" << std::endl;
}

void Pic::find_output_list(const PlotWriter& writerIn, long int& nPointAllProc,
                           PlotWriter::VectorPointList& pointList_II,
                           std::array<double, nDimMax>& xMin_D,
                           std::array<double, nDimMax>& xMax_D) {
  const auto plo = geom.ProbLo();
  const auto plh = geom.ProbHi();

  Real xMinL_D[nDimMax] = { plh[ix_], plh[iy_], plh[iz_] };
  Real xMaxL_D[nDimMax] = { plo[ix_], plo[iy_], plo[iz_] };

  const auto dx = geom.CellSize();

  int icount = 0;
  int iBlock = 0;
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.validbox();

    const auto& typeArr = nodeType[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k) {
      const double zp = k * dx[iz_] + plo[iz_];
      for (int j = lo.y; j <= hi.y; ++j) {
        const double yp = j * dx[iy_] + plo[iy_];
        for (int i = lo.x; i <= hi.x; ++i) {
          const double xp = i * dx[ix_] + plo[ix_];
          if (typeArr(i, j, k) == iHandle_ &&
              writerIn.is_inside_plot_region(i, j, k, xp, yp, zp)) {

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

  if (ParallelDescriptor::MyProc() == 0 && writerIn.get_plotDx() >= 0) {
    // Processor-0 output the inactive PIC nodes for structured output.
    Box gbx = convert(geom.Domain(), { 1, 1, 1 });

    if (writerIn.is_compact())
      gbx = convert(nodeBA.minimalBox(), { 1, 1, 1 });

    const auto lo = lbound(gbx);
    const auto hi = ubound(gbx);

    int kMax = hi.z;
    const bool is2D = geom.Domain().bigEnd(iz_) == geom.Domain().smallEnd(iz_);
    if (is2D)
      kMax = lo.z;

    for (int k = lo.z; k <= kMax; ++k) {
      const double zp = k * dx[iz_] + plo[iz_];
      for (int j = lo.y; j <= hi.y; ++j) {
        const double yp = j * dx[iy_] + plo[iy_];
        for (int i = lo.x; i <= hi.x; ++i) {
          const double xp = i * dx[ix_] + plo[ix_];
          if (writerIn.is_inside_plot_region(i, j, k, xp, yp, zp) &&
              !nodeBA.contains({ i, j, k })) {
            const int iBlock = -1;
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

void Pic::get_field_var(const VectorPointList& pointList_II,
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

      if (ParallelDescriptor::MyProc() == 0 && iBlock == -1) {
        // Processor-0 output the inactive PIC nodes for structured output.
        for (int iVar = 0; iVar < nVar; ++iVar) {
          var_II(iPoint, iVar) = get_var(sVar_I[iVar], ix, iy, iz, mfi, false);
        }
        iPoint++;
      } else if (iBlock == iBlockCount) {
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

double Pic::get_var(std::string var, const int ix, const int iy, const int iz,
                    const MFIter& mfi, bool isValidMFI) {

  auto get_is = [var]() {
    std::string::size_type pos;
    std::stringstream ss;
    int is;
    pos = var.find_first_of("0123456789");
    ss << var.substr(pos);
    ss >> is;
    return is;
  };

  double value = 0;
  if (isValidMFI || var.substr(0, 1) == "X" || var.substr(0, 1) == "Y" ||
      var.substr(0, 1) == "Z") {
    // If not isValidMFI, then it is not possible to output variables other than
    // 'X', 'Y', 'Z'
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

      if (var.substr(0, 1) == "u") {
        double rho = arr(ix, iy, iz, iRho_);
        if (rho != 0)
          value /= rho;
      }

    } else if (var.substr(0, 2) == "pS") {
      const amrex::Array4<amrex::Real const>& arr =
          nodePlasma[get_is()][mfi].array();
      value = (arr(ix, iy, iz, iPxx_) + arr(ix, iy, iz, iPyy_) +
               arr(ix, iy, iz, iPzz_)) /
              3.0;
    } else if (var.substr(0, 4) == "proc") {
      value = ParallelDescriptor::MyProc();
    } else if (var.substr(0, 5) == "block") {
      value = mfi.index();
    } else {
      value = 0;
    }
  }

  return value;
}

void Pic::save_restart_data() {
  VisMF::SetNOutFiles(64);

  std::string restartDir = "PC/restartOUT/";
  VisMF::Write(nodeE, restartDir + "nodeE");
  VisMF::Write(nodeB, restartDir + "nodeB");
  VisMF::Write(centerB, restartDir + "centerB");

  for (int iPart = 0; iPart < parts.size(); iPart++) {
    parts[iPart]->Checkpoint(restartDir, "particles" + std::to_string(iPart));
  }
}

void Pic::save_restart_header(std::ofstream& headerFile) {
  if (ParallelDescriptor::IOProcessor()) {
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

void Pic::read_restart() {

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

void Pic::write_log(bool doForce, bool doCreateFile) {
  if (doCreateFile && ParallelDescriptor::IOProcessor()) {
    std::stringstream ss;
    int time = tc->get_time_si(); // double to int.
    ss << "PC/plots/log_n" << std::setfill('0') << std::setw(8)
       << tc->get_cycle() << ".log";
    logFile = ss.str();
    std::ofstream of(logFile.c_str());
    of << "time nStep Etot Ee Eb Epart ";
    for (int i = 0; i < nSpecies; i++)
      of << " Epart" << i;
    of << std::endl;
    of.close();
  }

  if (tc->log.is_time_to(doForce)) {
    Real eEnergy = calc_E_field_energy();
    Real bEnergy = calc_B_field_energy();
    if (ParallelDescriptor::IOProcessor()) {
      std::ofstream of(logFile.c_str(), std::fstream::app);
      of.precision(15);
      of << std::scientific;
      of << tc->get_time_si() << "\t" << tc->get_cycle() << "\t"
         << "\t" << (eEnergy + bEnergy + plasmaEnergy[iTot]) << "\t" << eEnergy
         << "\t" << bEnergy << "\t" << plasmaEnergy[iTot];
      for (int i = 0; i < nSpecies; i++)
        of << "\t" << plasmaEnergy[i];
      of << std::endl;
      of.close();
    }
  }
}

void Pic::write_plots(bool doForce) {
  for (auto& plot : tc->plots) {
    if (plot.is_time_to(doForce)) {
      amrex::Print() << "Saving plot at time = " << tc->get_time_si()
                     << " (s) for " << plot.writer.get_plotString()
                     << std::endl;
      if (plot.writer.is_amrex_format()) {
        write_amrex(plot.writer, tc->get_time_si(), tc->get_cycle());
      } else {
        plot.writer.write(tc->get_time_si(), tc->get_cycle(),
                          find_output_list_caller, get_field_var_caller);
      }
    }
  }
}

void Pic::write_amrex(const PlotWriter& pw, double const timeNow,
                      int const iCycle) {
  Print() << "amrex::" << pw.get_amrex_filename(timeNow, iCycle) << std::endl;

  Geometry geomOut;
  { // Creating geomOut, which uses output length unit, for amrex format output.
    RealBox boxRangeOut;
    Real no2outL = pw.No2OutTable("X");
    for (int i = 0; i < nDimMax; i++) {
      boxRangeOut.setLo(i, geom.ProbLo(i) * no2outL);
      boxRangeOut.setHi(i, geom.ProbHi(i) * no2outL);
    }
    Array<int, nDim> periodicity;
    for (int i = 0; i < nDimMax; i++)
      periodicity[i] = geom.isPeriodic(i);
    geomOut.define(geom.Domain(), boxRangeOut, geom.Coord(), periodicity);
  }

  //-----lambda to write MultiFab in output units---------------
  auto write_single_level = [&](MultiFab& mf, std::string& filename,
                                Vector<std::string>& varNames) {
    // Save cell-centered, instead of the nodal, values, because the AMReX
    // document says some virtualiazaion tools assumes the AMReX format outputs
    // are cell-centered. 
    MultiFab centerMF;
    centerMF.define(centerBA, dm, mf.nComp(), 0);
    if (mf.is_nodal()) {
      average_node_to_cellcenter(centerMF, 0, mf, 0, mf.nComp(), 0);
    } else {
      MultiFab::Copy(centerMF, mf, 0, 0, centerMF.nComp(), 0);
    }

    // Convert to output unit.
    for (int i = 0; i < centerMF.nComp(); i++) {
      Real no2out = pw.No2OutTable(varNames[i]);
      centerMF.mult(no2out, i, 1);
    }

    WriteSingleLevelPlotfile(filename, centerMF, varNames, geomOut, timeNow,
                             iCycle);
  };
  //----------lambda end------------------------------------------

  { //------------nodeE-----------------
    Vector<std::string> varNames = { "Ex", "Ey", "Ez" };
    std::string filename = pw.get_amrex_filename(timeNow, iCycle) + "_E";
    write_single_level(nodeE, filename, varNames);
  }

  { //------------centerB------------------------
    Vector<std::string> varNames = { "Bx", "By", "Bz" };
    std::string filename = pw.get_amrex_filename(timeNow, iCycle) + "_B";
    write_single_level(centerB, filename, varNames);
  }

  { //-------------plasma---------------------

    // The order of the varname should be consistent with nodePlasma.
    Vector<std::string> varNames = { "rho", "ux",  "uy",  "uz",  "pxx",
                                     "pyy", "pzz", "pxy", "pxz", "pyz" };

    for (int i = 0; i < nSpecies; i++) {
      std::string filename = pw.get_amrex_filename(timeNow, iCycle) +
                             "_plasma" + std::to_string(i);

      MultiFab rho(nodePlasma[i], make_alias, iRho_, 1);
      if (rho.min(0) == 0)
        Abort("Error: density is zero!");

      // Get momentums
      MultiFab ux(nodePlasma[i], make_alias, iMx_, 1);
      MultiFab uy(nodePlasma[i], make_alias, iMy_, 1);
      MultiFab uz(nodePlasma[i], make_alias, iMz_, 1);

      // Convert momentum to velocity;
      MultiFab::Divide(ux, rho, 0, 0, 1, 0);
      MultiFab::Divide(uy, rho, 0, 0, 1, 0);
      MultiFab::Divide(uz, rho, 0, 0, 1, 0);

      MultiFab pl(nodePlasma[i], make_alias, iRho_, iPyz_ - iRho_ + 1);
      write_single_level(pl, filename, varNames);

      // Convert velocity to momentum
      MultiFab::Multiply(ux, rho, 0, 0, 1, 0);
      MultiFab::Multiply(uy, rho, 0, 0, 1, 0);
      MultiFab::Multiply(uz, rho, 0, 0, 1, 0);
    }
  }
}

void find_output_list_caller(const PlotWriter& writerIn,
                             long int& nPointAllProc,
                             PlotWriter::VectorPointList& pointList_II,
                             std::array<double, nDimMax>& xMin_D,
                             std::array<double, nDimMax>& xMax_D) {
  MPICs->pic.find_output_list(writerIn, nPointAllProc, pointList_II, xMin_D,
                              xMax_D);
}

void get_field_var_caller(const VectorPointList& pointList_II,
                          const std::vector<std::string>& sVar_I,
                          MDArray<double>& var_II) {
  MPICs->pic.get_field_var(pointList_II, sVar_I, var_II);
}
