#include <AMReX_PlotFileUtil.H>
#include <AMReX_RealVect.H>

#include "GridUtility.h"
#include "Pic.h"
#include "SimDomains.h"

using namespace amrex;

//==========================================================
void Pic::update_cells_for_pt() {
  // If fill_new_cells is called when PIC component is off, it suggests the test
  // particle component is activated. The test particle component copies EM
  // field from PIC, so PIC EM field should be updated here.
  if (doNeedFillNewCell || !usePIC)
    fill_new_cells();

  return;
}

//==========================================================
void Pic::get_fluid_state_for_points(const int nDim, const int nPoint,
                                     const double* const xyz_I,
                                     double* const data_I, const int nVar) {
  std::string nameFunc = "Pic::get_fluid_state_for_points";

  Print() << printPrefix << component
          << " -> GM coupling at t =" << tc->get_time_si() << " (s)"
          << std::endl;

  if (!usePIC) {
    for (int i = 0; i < nVar * nPoint; i++) {
      data_I[i] = 0;
    }
    return;
  }

  // (rho + 3*Moment + 6*p)*nSpecies+ 3*E + 3*B;
  const int nVarPerSpecies = 10;
  int nVarPIC = nSpecies * nVarPerSpecies + 6;
  double dataPIC_I[nVarPIC];

  const int iBx_ = nSpecies * nVarPerSpecies, iBy_ = iBx_ + 1, iBz_ = iBy_ + 1;
  const int iEx_ = iBz_ + 1;

  const RealBox& range = Geom(0).ProbDomain();
  for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    double pic_D[3] = { 0 };
    for (int iDim = 0; iDim < nDim; iDim++) {
      pic_D[iDim] = xyz_I[iPoint * nDim + iDim] * fi->get_Si2NoL();
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
            get_value_at_loc(nodePlasma[iSpecies], Geom(0), xp, yp, zp, iVar);
      }

    for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
      for (int iDir = ix_; iDir <= iz_; iDir++) {
        dataPIC_I[iBx_ + iDir] =
            get_value_at_loc(nodeB[iLevTest], Geom(iLevTest), xp, yp, zp, iDir);
      }
      for (int iDir = ix_; iDir <= iz_; iDir++) {
        dataPIC_I[iEx_ + iDir] =
            get_value_at_loc(nodeE[iLevTest], Geom(iLevTest), xp, yp, zp, iDir);
      }
    }

    // Combine PIC plasma data into MHD fluid data.
    fi->calc_fluid_state(dataPIC_I, &data_I[iPoint * nVar]);
  }
}

//==========================================================
void Pic::find_output_list(const PlotWriter& writerIn, long int& nPointAllProc,
                           PlotWriter::VectorPointList& pointList_II,
                           std::array<double, nDim>& xMin_D,
                           std::array<double, nDim>& xMax_D) {
  if (isGridEmpty)
    return;
  // Loop not implemented correctly // Talha
  int iLevTest = 0;
  const auto plo = Geom(0).ProbLo();
  const auto plh = Geom(0).ProbHi();

  Real xMinL_D[nDim] = { plh[ix_], plh[iy_], plh[iz_] };
  Real xMaxL_D[nDim] = { plo[ix_], plo[iy_], plo[iz_] };

  const auto dx = Geom(0).CellSize();

  const Box& gbx = convert(Geom(0).Domain(), { 1, 1, 1 });

  const auto glo = lbound(gbx);
  const auto ghi = ubound(gbx);

  int iBlock = 0;
  for (MFIter mfi(nodeE[iLevTest]); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();

    const auto& typeArr = nodeShare[mfi].array();

    auto lo = box.loVect3d();
    auto hi = box.hiVect3d();

    // Do not output the rightmost nodes for periodic boundary.
    for (int iDim = 0; iDim < nDim; iDim++)
      if ((Geom(0).isPeriodic(iDim)) && gbx.bigEnd(iDim) == hi[iDim])
        hi[iDim]--;

    auto do_output_this_node = [&](int i, int j, int k) {
      int type = typeArr(i, j, k);

      if (type == iAssign_)
        return true;

      if (type != iIgnore_) {
        if (Geom(0).isPeriodic(ix_) && i == glo.x && !test_bit(type, ix_))
          return true;
        if (Geom(0).isPeriodic(iy_) && j == glo.y && !test_bit(type, iy_))
          return true;
        if (!isFake2D && Geom(0).isPeriodic(iz_) && k == glo.z &&
            !test_bit(type, iz_))
          return true;
      }

      return false;
    };

    for (int k = lo[iz_]; k <= hi[iz_]; ++k) {
      const double zp = k * dx[iz_] + plo[iz_];
      for (int j = lo[iy_]; j <= hi[iy_]; ++j) {
        const double yp = j * dx[iy_] + plo[iy_];
        for (int i = lo[ix_]; i <= hi[ix_]; ++i) {
          const double xp = i * dx[ix_] + plo[ix_];
          if (do_output_this_node(i, j, k) &&
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
    Box gbx = convert(Geom(0).Domain(), { 1, 1, 1 });

    if (writerIn.is_compact())
      gbx = convert(nGrids[0].minimalBox(), { 1, 1, 1 });

    const auto lo = lbound(gbx);
    const auto hi = ubound(gbx);

    int iMax = hi.x, jMax = hi.y, kMax = hi.z;
    if (Geom(0).isPeriodic(ix_))
      --iMax;
    if (Geom(0).isPeriodic(iy_))
      --jMax;
    if (Geom(0).isPeriodic(iz_))
      --kMax;

    if (isFake2D)
      kMax = lo.z;

    for (int k = lo.z; k <= kMax; ++k) {
      const double zp = k * dx[iz_] + plo[iz_];
      for (int j = lo.y; j <= jMax; ++j) {
        const double yp = j * dx[iy_] + plo[iy_];
        for (int i = lo.x; i <= iMax; ++i) {
          const double xp = i * dx[ix_] + plo[ix_];
          if (writerIn.is_inside_plot_region(i, j, k, xp, yp, zp) &&
              !nGrids[0].contains(IntVect{ i, j, k })) {
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

  // loop not implemented correctly // Talha
  int iLevTest = 0;
  const int iBlk_ = 6;

  long nPoint = pointList_II.size();
  int nVar = sVar_I.size();

  int iBlockCount = 0;
  long iPoint = 0;
  for (MFIter mfi(nodeE[iLevTest]); mfi.isValid(); ++mfi) {
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
      const auto plo = Geom(0).ProbLo();
      const auto dx = Geom(0).CellSize();
      value = ix * dx[ix_] + plo[ix_];
    } else if (var.substr(0, 1) == "Y") {
      const auto plo = Geom(0).ProbLo();
      const auto dx = Geom(0).CellSize();
      value = iy * dx[iy_] + plo[iy_];
    } else if (var.substr(0, 1) == "Z") {
      const auto plo = Geom(0).ProbLo();
      const auto dx = Geom(0).CellSize();
      value = iz * dx[iz_] + plo[iz_];
    } else if (var.substr(0, 2) == "Ex") {
      const amrex::Array4<amrex::Real const>& arr =
          nodeE[0][mfi].array(); // Talha- check // no loop
      value = arr(ix, iy, iz, ix_);
    } else if (var.substr(0, 2) == "Ey") {
      const amrex::Array4<amrex::Real const>& arr = nodeE[0][mfi].array();
      value = arr(ix, iy, iz, iy_);
    } else if (var.substr(0, 2) == "Ez") {
      const amrex::Array4<amrex::Real const>& arr = nodeE[0][mfi].array();
      value = arr(ix, iy, iz, iz_);
    } else if (var.substr(0, 2) == "Bx") {
      const amrex::Array4<amrex::Real const>& arr = nodeB[0][mfi].array();
      value = arr(ix, iy, iz, ix_);
    } else if (var.substr(0, 2) == "By") {
      const amrex::Array4<amrex::Real const>& arr = nodeB[0][mfi].array();
      value = arr(ix, iy, iz, iy_);
    } else if (var.substr(0, 2) == "Bz") {
      const amrex::Array4<amrex::Real const>& arr = nodeB[0][mfi].array();
      value = arr(ix, iy, iz, iz_);
    } else if (var.substr(0, 4) == "rhoS" || var.substr(0, 3) == "uxS" ||
               var.substr(0, 3) == "uyS" || var.substr(0, 3) == "uzS" ||
               var.substr(0, 4) == "pXXS" || var.substr(0, 4) == "pYYS" ||
               var.substr(0, 4) == "pZZS" || var.substr(0, 4) == "pXYS" ||
               var.substr(0, 4) == "pXZS" || var.substr(0, 4) == "pYZS" ||
               var.substr(0, 4) == "numS") {

      // The last element of nodePlasma is the sum of all species.
      if (get_is() >= nodePlasma.size() - 1) {
        value = 0;
      } else {
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
      }
    } else if (var.substr(0, 2) == "pS") {
      const amrex::Array4<amrex::Real const>& arr =
          nodePlasma[get_is()][mfi].array();
      value = (arr(ix, iy, iz, iPxx_) + arr(ix, iy, iz, iPyy_) +
               arr(ix, iy, iz, iPzz_)) /
              3.0;
    } else if (var.substr(0, 2) == "qc") {
      const amrex::Array4<amrex::Real const>& arr =
          centerNetChargeN[mfi].array();
      value = arr(ix, iy, iz);
    } else if (var.substr(0, 5) == "divEc") {
      const amrex::Array4<amrex::Real const>& arr = centerDivE[mfi].array();
      value = arr(ix, iy, iz);
    } else if (var.substr(0, 3) == "phi") {
      const amrex::Array4<amrex::Real const>& arr = centerPhi[mfi].array();
      value = arr(ix, iy, iz);
    } else if (var.substr(0, 7) == "smoothE") {
      const amrex::Array4<amrex::Real const>& arr = nodeSmoothCoef[mfi].array();
      value = arr(ix, iy, iz);
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
  if (isGridEmpty)
    return;

  std::string restartDir = component + "/restartOUT/";
  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    VisMF::Write(nodeE[iLevTest], restartDir + gridName + "_nodeE");
    VisMF::Write(nodeB[iLevTest], restartDir + gridName + "_nodeB");
    VisMF::Write(centerB[iLevTest], restartDir + gridName + "_centerB");
  }

  for (int iPart = 0; iPart < parts.size(); iPart++) {
    parts[iPart]->label_particles_outside_ba();
    parts[iPart]->Redistribute();
    parts[iPart]->Checkpoint(restartDir,
                             gridName + "_particles" + std::to_string(iPart));
  }
  inject_particles_for_boundary_cells();
}

//==========================================================
void Pic::save_restart_header(std::ofstream& headerFile) {
  std::string command_suffix = "_" + gridName + "\n";

  if (ParallelDescriptor::IOProcessor()) {
    headerFile << "#ELECTRON" + command_suffix;
    headerFile << qomEl << "\t\telectronChargePerMass\n";
    headerFile << "\n";

    headerFile << "#PARTICLES" + command_suffix;
    headerFile << nPartPerCell[ix_] << "\t\t\tnParticleX\n";
    headerFile << nPartPerCell[iy_] << "\t\t\tnParticleY\n";
    headerFile << nPartPerCell[iz_] << "\t\t\tnParticleZ\n";
    headerFile << "\n";
  }
}

//==========================================================
void Pic::read_restart() {
  Print() << "Pic::read_restart() start....." << std::endl;

  std::string restartDir = component + "/restartIN/";

  for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
    VisMF::Read(nodeE[iLevTest], restartDir + gridName + "_nodeE");
    VisMF::Read(nodeB[iLevTest], restartDir + gridName + "_nodeB");
    VisMF::Read(centerB[iLevTest], restartDir + gridName + "_centerB");
  }

  for (int iPart = 0; iPart < parts.size(); iPart++) {
    parts[iPart]->Restart(restartDir,
                          gridName + "_particles" + std::to_string(iPart));
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
  if (isGridEmpty || !usePIC)
    return;

  if (doCreateFile && ParallelDescriptor::IOProcessor()) {
    std::stringstream ss;
    ss << component << "/plots/log_pic_n" << std::setfill('0') << std::setw(8)
       << tc->get_cycle() << ".log";
    logFile = ss.str();
    std::ofstream of(logFile.c_str());
    of << "time nStep Etot Ee Eb Epart ";
    for (int i = 0; i < nSpecies; i++)
      of << " Epart" << i;
    of << std::endl;
    of.close();
  }

  if (tc->picLog.is_time_to(doForce)) {
    ParallelDescriptor::ReduceRealSum(plasmaEnergy.data(), plasmaEnergy.size(),
                                      ParallelDescriptor::IOProcessorNumber());

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
  if (isGridEmpty)
    return;
  for (auto& plot : tc->plots) {
    if (plot.is_time_to(doForce)) {
      amrex::Print() << printPrefix
                     << " Saving plot at time = " << tc->get_time_si()
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
  int iSpecies = pw.get_particleSpecies();

  std::string dirName = pw.get_amrex_filename(timeNow, iCycle);

  Vector<int> writeRealComp;
  for (int i = 0; i < nPicPartReal; ++i) {
    writeRealComp.push_back(1);
  }

  Vector<std::string> realCompNames;
  realCompNames.resize(nPicPartReal);
  realCompNames[Particles<>::iup_] = "velocity_x";
  realCompNames[Particles<>::ivp_] = "velocity_y";
  realCompNames[Particles<>::iwp_] = "velocity_z";
  realCompNames[Particles<>::iqp_] = "weight";

  Vector<int> writeIntComp;
  Vector<std::string> intCompNames;

  Real no2outL = pw.No2OutTable("X");
  Real no2outV = pw.No2OutTable("u");
  Real no2outM = pw.No2OutTable("mass");

  RealBox outRange;

  bool isCut = pw.get_plotString().find("cut") != std::string::npos;

  if (isCut)
    for (int iDim = 0; iDim < nDim; iDim++) {
      outRange.setLo(iDim, pw.get_plotMin_D(iDim) * no2outL);
      outRange.setHi(iDim, pw.get_plotMax_D(iDim) * no2outL);
    }

  BoxArray baIO;
  if (isCut) {
    // Find the box that contains the output region.
    const auto lo = outRange.lo();
    const auto hi = outRange.hi();

    const auto plo = Geom(0).ProbLo();
    const auto plh = Geom(0).ProbHi();
    const auto dx = Geom(0).CellSize();

    IntVect cellLo, cellHi;

    for (int i = 0; i < 3; i++) {
      cellLo[i] = fastfloor((lo[i] / no2outL - plo[i]) / dx[i]);
      cellHi[i] = fastfloor((hi[i] / no2outL - plo[i]) / dx[i]);
    }

    if (isFake2D) {
      cellLo[iz_] = 0;
      cellHi[iz_] = 0;
    }

    baIO.define(Box(cellLo, cellHi));
    baIO.maxSize(IntVect(8, 8, 8));
  } else {
    baIO = cGrid;
  }

  // Create a new grid for saving data in IO units
  Geometry geomOut;
  set_IO_geom(geomOut, pw);

  AmrInfo amrInfo;
  amrInfo.blocking_factor.clear();
  amrInfo.blocking_factor.push_back(IntVect(1, 1, 1));

  Grid gridIO(geomOut, amrInfo, 0, -gridID);
  gridIO.set_base_grid(baIO);
  gridIO.InitFromScratch(0.0);

  IOParticles particlesOut(*parts[iSpecies].get(), &gridIO, no2outL, no2outV,
                           no2outM, outRange);

  // Saving field/coordinates information. It is required by yt for unknown
  // reasons.
  write_amrex_field(pw, timeNow, iCycle, "E B", dirName, baIO);

  particlesOut.WritePlotFile(dirName, "particle", writeRealComp, writeIntComp,
                             realCompNames, intCompNames);
}

void Pic::set_IO_geom(amrex::Geometry& geomIO, const PlotWriter& pw) {
  // Creating geomIO, which uses output length unit, for amrex format output.
  RealBox boxRangeOut;
  Real no2outL = pw.No2OutTable("X");
  for (int i = 0; i < nDim; i++) {
    boxRangeOut.setLo(i, Geom(0).ProbLo(i) * no2outL);
    boxRangeOut.setHi(i, Geom(0).ProbHi(i) * no2outL);
  }
  Array<int, nDim> periodicity;
  for (int i = 0; i < nDim; i++)
    periodicity[i] = Geom(0).isPeriodic(i);
  geomIO.define(Geom(0).Domain(), boxRangeOut, Geom(0).Coord(), periodicity);
}

//==========================================================
void Pic::write_amrex_field(const PlotWriter& pw, double const timeNow,
                            int const iCycle, const std::string plotVars,
                            const std::string filenameIn,
                            const BoxArray baOut) {
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

  if (usePIC && plotVars.find("plasma") != std::string::npos)
    nVarOut += 10 * nSpecies;

  // Save cell-centered, instead of the nodal, values, because the AMReX
  // document says some virtualiazaion tools assumes the AMReX format outputs
  // are cell-centered.

  MultiFab centerMF;
  centerMF.define(cGrid, DistributionMap(0), nVarOut, 0);

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
        // Use the value returned from Geom(0) instead of geomOut, and it will
        // be converted to output unit just as other variables.
        const Real z0 = Geom(0).CellCenter(k, iz_);
        for (int j = lo.y; j <= hi.y; ++j) {
          const Real y0 = Geom(0).CellCenter(j, iy_);
          for (int i = lo.x; i <= hi.x; ++i) {
            const Real x0 = Geom(0).CellCenter(i, ix_);
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
    for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
      MultiFab::Copy(centerMF, centerB[iLevTest], 0, iStart,
                     nodeB[iLevTest].nComp(), 0);
      iStart += nodeB[iLevTest].nComp();
    }

    varNames.push_back("Bx");
    varNames.push_back("By");
    varNames.push_back("Bz");
  }

  if (plotVars.find("E") != std::string::npos) {
    //-----------------E-----------------------------
    for (int iLevTest = 0; iLevTest <= finest_level; iLevTest++) {
      average_node_to_cellcenter(centerMF, iStart, nodeE[iLevTest], 0, nodeE[iLevTest].nComp(), 0);
      iStart += nodeE[iLevTest].nComp();
    }

    varNames.push_back("Ex");
    varNames.push_back("Ey");
    varNames.push_back("Ez");
  }

  bool isDensityZero = false;
  int zeroI, zeroJ, zeroK;
  if (usePIC && plotVars.find("plasma") != std::string::npos) {
    //-------------plasma---------------------

    // The order of the varname should be consistent with nodePlasma.
    Vector<std::string> plasmaNames = { "rho", "ux",  "uy",  "uz",  "pxx",
                                        "pyy", "pzz", "pxy", "pxz", "pyz" };

    for (int i = 0; i < nSpecies; i++) {
      MultiFab rho(nodePlasma[i], make_alias, iRho_, 1);

      // Get momentums
      MultiFab ux(nodePlasma[i], make_alias, iMx_, 1);
      MultiFab uy(nodePlasma[i], make_alias, iMy_, 1);
      MultiFab uz(nodePlasma[i], make_alias, iMz_, 1);

      for (MFIter mfi(nodePlasma[i]); mfi.isValid(); ++mfi) {
        // Convert momentum to velocity;
        const Box& box = mfi.validbox();
        const Array4<Real>& plasmaArr = nodePlasma[i][mfi].array();
        const Array4<Real>& uxArr = ux[mfi].array();
        const Array4<Real>& uyArr = uy[mfi].array();
        const Array4<Real>& uzArr = uz[mfi].array();

        const auto lo = lbound(box);
        const auto hi = ubound(box);

        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              const Real rho = plasmaArr(i, j, k, iRho_);
              if (rho > 1e-99) {
                uxArr(i, j, k) = plasmaArr(i, j, k, iUx_) / rho;
                uyArr(i, j, k) = plasmaArr(i, j, k, iUy_) / rho;
                uzArr(i, j, k) = plasmaArr(i, j, k, iUz_) / rho;
              } else {
                isDensityZero = true;
                zeroI = i;
                zeroJ = j;
                zeroK = k;
              }
            }
      }

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

  if (!baOut.empty()) {
    distribute_FabArray(centerMF, baOut, DistributionMapping(baOut));
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
    headerFile << fi->get_rPlanet_SI() << "\n";
  }

  if (isDensityZero) {
    AllPrint()
        << "\n\n\n==========" << printPrefix
        << " Error ===========================\n"
        << "Density is zero at i = " << zeroI << " j = " << zeroJ
        << " k = " << zeroK << ".\n"
        << "Check the file " << filename << " to see what is going on. \n"
        << "Suggestions:\n"
        << "1) Use the #RESAMPLING command to control the particle number.\n"
        << "2) If 1) does not help, it is likely something is wrong at the PIC "
           "boundary."
        << "\n========================================================\n\n\n"
        << std::endl;
    Abort();
  }
}

//==========================================================
void find_output_list_caller(const PlotWriter& writerIn,
                             long int& nPointAllProc,
                             PlotWriter::VectorPointList& pointList_II,
                             std::array<double, nDim>& xMin_D,
                             std::array<double, nDim>& xMax_D) {
  fleksDomains(fleksDomains.selected())
      .pic->find_output_list(writerIn, nPointAllProc, pointList_II, xMin_D,
                             xMax_D);
}

//==========================================================
void get_field_var_caller(const VectorPointList& pointList_II,
                          const std::vector<std::string>& sVar_I,
                          MDArray<double>& var_II) {
  fleksDomains(fleksDomains.selected())
      .pic->get_field_var(pointList_II, sVar_I, var_II);
}
