#include <AMReX_PlotFileUtil.H>
#include <AMReX_RealVect.H>

#include "Pic.h"
#include "SWMFDomains.h"

using namespace amrex;

//==========================================================
void Pic::set_state_var(double* data, int* index) {
  std::string nameFunc = "Pic::set_state_var";

  Print() << nameFunc << std::endl;

  fluidInterface->set_couple_node_value(data, index);

  if (doNeedFillNewCell)
    fill_new_cells();

  return;
}

//==========================================================
int Pic::get_grid_nodes_number() {
  return fluidInterface->count_couple_node_number();
}

//==========================================================
void Pic::get_grid(double* pos_DI) {
  std::string nameFunc = "Pic::get_grid";
  fluidInterface->get_couple_node_loc(pos_DI);
  return;
}

//==========================================================
void Pic::find_mpi_rank_for_points(const int nPoint, const double* const xyz_I,
                                   int* const rank_I) {
  int nDimGM = fluidInterface->getnDim();
  amrex::Real si2nol = fluidInterface->getSi2NoL();
  const RealBox& range = geom.ProbDomain();
  for (int i = 0; i < nPoint; i++) {
    amrex::Real x = xyz_I[i * nDimGM + ix_] * si2nol;
    amrex::Real y = xyz_I[i * nDimGM + iy_] * si2nol;
    amrex::Real z = 0;
    if (nDimGM > 2)
      z = xyz_I[i * nDimGM + iz_] * si2nol;
    // Check if this point is inside this FLEKS domain.
    if (range.contains(RealVect(x, y, z), 1e-10))
      rank_I[i] = find_mpi_rank_from_coord(x, y, z);
  }
}

//==========================================================
void Pic::get_fluid_state_for_points(const int nDim, const int nPoint,
                                     const double* const xyz_I,
                                     double* const data_I, const int nVar) {
  std::string nameFunc = "Pic::get_fluid_state_for_points";
  Print() << nameFunc << std::endl;

  // (rho + 3*Moment + 6*p)*nSpecies+ 3*E + 3*B;
  const int nVarPerSpecies = 10;
  int nVarPIC = nSpecies * nVarPerSpecies + 6;
  double dataPIC_I[nVarPIC];

  const int iBx_ = nSpecies * nVarPerSpecies, iBy_ = iBx_ + 1, iBz_ = iBy_ + 1;
  const int iEx_ = iBz_ + 1, iEy_ = iEx_ + 1, iEz_ = iEy_ + 1;

  const RealBox& range = geom.ProbDomain();
  for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    double pic_D[3] = { 0 };
    for (int iDim = 0; iDim < nDim; iDim++) {
      pic_D[iDim] = xyz_I[iPoint * nDim + iDim] * fluidInterface->getSi2NoL();
    }

    const Real xp = pic_D[0];
    const Real yp = (nDim > 1) ? pic_D[1] : 0.0;
    const Real zp = (nDim > 2) ? pic_D[2] : 0.0;

    // Check if this point is inside this FLEKS domain.
    if (!range.contains(RealVect(xp, yp, zp), 1e-10))
      continue;

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
}

//==========================================================
void Pic::find_output_list(const PlotWriter& writerIn, long int& nPointAllProc,
                           PlotWriter::VectorPointList& pointList_II,
                           std::array<double, nDim>& xMin_D,
                           std::array<double, nDim>& xMax_D) {
  const auto plo = geom.ProbLo();
  const auto plh = geom.ProbHi();

  Real xMinL_D[nDim] = { plh[ix_], plh[iy_], plh[iz_] };
  Real xMaxL_D[nDim] = { plo[ix_], plo[iy_], plo[iz_] };

  const auto dx = geom.CellSize();

  int icount = 0;
  int iBlock = 0;
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    FArrayBox& fab = nodeE[mfi];
    const Box& box = mfi.validbox();

    const auto& typeArr = nodeAssignment[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k) {
      const double zp = k * dx[iz_] + plo[iz_];
      for (int j = lo.y; j <= hi.y; ++j) {
        const double yp = j * dx[iy_] + plo[iy_];
        for (int i = lo.x; i <= hi.x; ++i) {
          const double xp = i * dx[ix_] + plo[ix_];
          if (typeArr(i, j, k) == iAssign_ &&
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

  ParallelDescriptor::ReduceRealMin(xMinL_D, nDim);
  ParallelDescriptor::ReduceRealMax(xMaxL_D, nDim);

  for (int iDim = 0; iDim < nDim; ++iDim) {
    xMin_D[iDim] = xMinL_D[iDim];
    xMax_D[iDim] = xMaxL_D[iDim];
  }
}

//==========================================================
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

//==========================================================
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
    } else if (var.substr(0, 4) == "rank") {
      value = ParallelDescriptor::MyProc();
    } else if (var.substr(0, 5) == "block") {
      value = mfi.index();
    } else {
      value = 0;
    }
  }

  return value;
}

//==========================================================
void Pic::save_restart_data() {
  std::string restartDir = "PC/restartOUT/";
  VisMF::Write(nodeE, restartDir + domainName + "_nodeE");
  VisMF::Write(nodeB, restartDir + domainName + "_nodeB");
  VisMF::Write(centerB, restartDir + domainName + "_centerB");

  for (int iPart = 0; iPart < parts.size(); iPart++) {
    parts[iPart]->label_particles_outside_ba();
    parts[iPart]->Redistribute();
    parts[iPart]->Checkpoint(restartDir,
                             domainName + "_particles" + std::to_string(iPart));
  }
  inject_particles_for_boundary_cells();
}

//==========================================================
void Pic::save_restart_header(std::ofstream& headerFile) {
  std::string command_suffix = "_" + domainName + "\n";

  if (ParallelDescriptor::IOProcessor()) {
    headerFile << "#ELECTRON"+command_suffix;
    headerFile << qomEl << "\t qomEl\n";
    headerFile << "\n";

    headerFile << "#PARTICLES"+command_suffix;
    for (int i = 0; i < nDim; ++i) {
      headerFile << nPartPerCell[i] << "\n";
    }
    headerFile << "\n";
  }
}

//==========================================================
void Pic::read_restart() {
  Print() << "Pic::read_restart() start....." << std::endl;

  std::string restartDir = "PC/restartIN/";
  VisMF::Read(nodeE, restartDir + domainName + "_nodeE");
  VisMF::Read(nodeB, restartDir + domainName + "_nodeB");
  VisMF::Read(centerB, restartDir + domainName + "_centerB");

  for (int iPart = 0; iPart < parts.size(); iPart++) {
    parts[iPart]->Restart(restartDir,
                          domainName + "_particles" + std::to_string(iPart));
  }
  inject_particles_for_boundary_cells();

  sum_moments();
  sum_to_center(false);

  // The default of doNeedFillNewCell is true. The PIC cells have been filled in
  // here, so turn it of
  doNeedFillNewCell = false;
}

//==========================================================
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

//==========================================================
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

//==========================================================
void Pic::write_amrex(const PlotWriter& pw, double const timeNow,
                      int const iCycle) {

  if (pw.is_particle()) {
    write_amrex_particle(pw, timeNow, iCycle);
  } else {
    write_amrex_field(pw, timeNow, iCycle);
  }
}

//==========================================================
void Pic::write_amrex_particle(const PlotWriter& pw, double const timeNow,
                               int const iCycle) {
  const int nReal = 4;
  int iSpecies = pw.get_particleSpecies();

  std::string dirName = pw.get_amrex_filename(timeNow, iCycle);

  // Saving field/coordinates information. It is required by yt, but I do not
  // know whether it is necessary at this point. --Yuxi
  write_amrex_field(pw, timeNow, iCycle, std::string(), dirName);

  Vector<int> writeRealComp;
  for (int i = 0; i < nReal; ++i) {
    writeRealComp.push_back(1);
  }

  Vector<std::string> realCompNames;
  realCompNames.resize(nReal);
  realCompNames[Particles::iup_] = "velocity_x";
  realCompNames[Particles::ivp_] = "velocity_y";
  realCompNames[Particles::iwp_] = "velocity_z";
  realCompNames[Particles::iqp_] = "weight";

  Vector<int> writeIntComp;
  Vector<std::string> intCompNames;

  Geometry geomOut;
  set_IO_geom(geomOut, pw);
  Real no2outL = pw.No2OutTable("X");
  Real no2outV = pw.No2OutTable("u");
  Real no2outM = pw.No2OutTable("mass");

  RealBox outRange;

  if (pw.get_plotString().find("cut") != std::string::npos)
    for (int iDim = 0; iDim < nDim; iDim++) {
      outRange.setLo(iDim, pw.get_plotMin_D(iDim));
      outRange.setHi(iDim, pw.get_plotMax_D(iDim));
    }

  IOParticles particlesOut(*parts[iSpecies].get(), geomOut, no2outL, no2outV,
                           no2outM, outRange);

  particlesOut.WritePlotFile(dirName, "particle", writeRealComp, writeIntComp,
                             realCompNames, intCompNames);
}

void Pic::set_IO_geom(amrex::Geometry& geomIO, const PlotWriter& pw) {
  // Creating geomIO, which uses output length unit, for amrex format output.
  RealBox boxRangeOut;
  Real no2outL = pw.No2OutTable("X");
  for (int i = 0; i < nDim; i++) {
    boxRangeOut.setLo(i, geom.ProbLo(i) * no2outL);
    boxRangeOut.setHi(i, geom.ProbHi(i) * no2outL);
  }
  Array<int, nDim> periodicity;
  for (int i = 0; i < nDim; i++)
    periodicity[i] = geom.isPeriodic(i);
  geomIO.define(geom.Domain(), boxRangeOut, geom.Coord(), periodicity);
}

//==========================================================
void Pic::write_amrex_field(const PlotWriter& pw, double const timeNow,
                            int const iCycle, const std::string plotVars,
                            const std::string filenameIn) {
  Print() << "amrex::" << pw.get_amrex_filename(timeNow, iCycle) << std::endl;

  Geometry geomOut;
  set_IO_geom(geomOut, pw);

  int nVarOut = 0;

  if (plotVars.find("X") != std::string::npos)
    nVarOut += 3;

  if (plotVars.find("B") != std::string::npos)
    nVarOut += 3;

  if (plotVars.find("E") != std::string::npos)
    nVarOut += 3;

  if (plotVars.find("plasma") != std::string::npos)
    nVarOut += 10 * nSpecies;

  // Save cell-centered, instead of the nodal, values, because the AMReX
  // document says some virtualiazaion tools assumes the AMReX format outputs
  // are cell-centered.

  MultiFab centerMF;
  centerMF.define(centerBA, dm, nVarOut, 0);

  Vector<std::string> varNames;
  int iStart = 0;

  if (plotVars.find("X") != std::string::npos) {
    //--------------Coordinates-------------------------
    MultiFab xyz(centerMF, make_alias, 0, 3);

    for (MFIter mfi(xyz); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& cellArr = xyz[mfi].array();
      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k) {
        // Use the value returned from geom instead of geomOut, and it will be
        // converted to output unit just as other variables.
        const Real z0 = geom.CellCenter(k, iz_);
        for (int j = lo.y; j <= hi.y; ++j) {
          const Real y0 = geom.CellCenter(j, iy_);
          for (int i = lo.x; i <= hi.x; ++i) {
            const Real x0 = geom.CellCenter(i, ix_);
            cellArr(i, j, k, ix_) = x0;
            cellArr(i, j, k, iy_) = y0;
            cellArr(i, j, k, iz_) = z0;
          }
        }
      }
    }
    iStart += 3;
    varNames.push_back("X");
    varNames.push_back("Y");
    varNames.push_back("Z");
  }

  if (plotVars.find("B") != std::string::npos) {
    //------------------B---------------
    MultiFab::Copy(centerMF, centerB, 0, iStart, nodeB.nComp(), 0);
    iStart += nodeB.nComp();
    varNames.push_back("Bx");
    varNames.push_back("By");
    varNames.push_back("Bz");
  }

  if (plotVars.find("E") != std::string::npos) {
    //-----------------E-----------------------------
    average_node_to_cellcenter(centerMF, iStart, nodeE, 0, nodeE.nComp(), 0);
    iStart += nodeE.nComp();
    varNames.push_back("Ex");
    varNames.push_back("Ey");
    varNames.push_back("Ez");
  }

  if (plotVars.find("plasma") != std::string::npos) {
    //-------------plasma---------------------

    // The order of the varname should be consistent with nodePlasma.
    Vector<std::string> plasmaNames = { "rho", "ux",  "uy",  "uz",  "pxx",
                                        "pyy", "pzz", "pxy", "pxz", "pyz" };

    for (int i = 0; i < nSpecies; i++) {
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
      average_node_to_cellcenter(centerMF, iStart, pl, 0, pl.nComp(), 0);
      iStart += pl.nComp();
      for (auto& var : plasmaNames) {
        varNames.push_back(var + "s" + std::to_string(i));
      }

      // Convert velocity to momentum
      MultiFab::Multiply(ux, rho, 0, 0, 1, 0);
      MultiFab::Multiply(uy, rho, 0, 0, 1, 0);
      MultiFab::Multiply(uz, rho, 0, 0, 1, 0);
    }
  }

  for (int i = 0; i < centerMF.nComp(); i++) {
    Real no2out = pw.No2OutTable(varNames[i]);
    centerMF.mult(no2out, i, 1);
  }

  std::string filename = pw.get_amrex_filename(timeNow, iCycle);
  if (!filenameIn.empty()) {
    filename = filenameIn;
  }

  WriteSingleLevelPlotfile(filename, centerMF, varNames, geomOut, timeNow,
                           iCycle);

  if (ParallelDescriptor::IOProcessor()) {
    // Write FLEKS header
    const std::string headerName = filename + "/FLEKSHeader";
    std::ofstream headerFile;
    headerFile.open(headerName.c_str(), std::ios::out | std::ios::trunc);

    if (!headerFile.good())
      amrex::FileOpenFailed(headerName);

    headerFile << pw.get_plotString() << "\n";
    headerFile << fluidInterface->getrPlanet() << "\n";
  }
}

//==========================================================
void find_output_list_caller(const PlotWriter& writerIn,
                             long int& nPointAllProc,
                             PlotWriter::VectorPointList& pointList_II,
                             std::array<double, nDim>& xMin_D,
                             std::array<double, nDim>& xMax_D) {
  FLEKSs(FLEKSs.selected())
      .pic.find_output_list(writerIn, nPointAllProc, pointList_II, xMin_D,
                            xMax_D);
}

//==========================================================
void get_field_var_caller(const VectorPointList& pointList_II,
                          const std::vector<std::string>& sVar_I,
                          MDArray<double>& var_II) {
  FLEKSs(FLEKSs.selected()).pic.get_field_var(pointList_II, sVar_I, var_II);
}
