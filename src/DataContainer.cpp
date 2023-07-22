
#include <fstream>
#include <iostream>
#include <vector>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "DataContainer.h"
#include "GridUtility.h"

using namespace amrex;

void AMReXDataContainer::read_header(std::string& headerName, int& nVar,
                                     int& nDim, amrex::Real& time,
                                     int& finest_level, amrex::RealBox& domain,
                                     amrex::Box& cellBox,
                                     amrex::Vector<std::string>& varNames) {
  std::ifstream HeaderFile;
  HeaderFile.open(headerName, std::ifstream::in);

  HeaderFile.precision(17);

  std::string versionName;

  HeaderFile >> versionName;
  HeaderFile >> nVar;

  varNames.clear();
  for (int ivar = 0; ivar < nVar; ++ivar) {
    std::string var;
    HeaderFile >> var;
    varNames.push_back(var);
  }

  HeaderFile >> nDim;
  HeaderFile >> time;
  HeaderFile >> finest_level;

  for (int i = 0; i < nDim; ++i) {
    amrex::Real lo;
    HeaderFile >> lo;
    domain.setLo(i, lo);
  }

  for (int i = 0; i < nDim; ++i) {
    amrex::Real hi;
    HeaderFile >> hi;
    domain.setHi(i, hi);
  }

  for (int i = 0; i < finest_level; ++i) {
    int refRatio;
    HeaderFile >> refRatio;
  }

  // Read the base grid and ignore refined leve grids.
  HeaderFile >> cellBox;
  HeaderFile.ignore(100000, '\n');

  HeaderFile.close();
}

void AMReXDataContainer::read_header() {

  std::string headerName = dirIn + "/Header";

  int finestLev;
  AMReXDataContainer::read_header(headerName, nVar, nDim, time, finestLev,
                                  domain, cellBox, varNames);

  SetFinestLevel(finestLev);
}

void AMReXDataContainer::read() {
  Print() << "Reading in " << dirIn << std::endl;

  read_header();

  Grid grid(Geom(0), get_amr_info(), nGst);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    VisMF::Read(mf[iLev], dirIn + "/Level_" + std::to_string(iLev) + "/Cell");

    grid.SetBoxArray(iLev, mf[iLev].boxArray());
    grid.SetDistributionMap(iLev, mf[iLev].DistributionMap());
  }
  grid.SetFinestLevel(n_lev() - 1);

  // Print() << "grids.ba0 = " << grid.boxArray(0) << std::endl;

  regrid(grid.boxArray(0), &grid);
}

void AMReXDataContainer::write() {

  int iLev = 0;

  if (!mf[iLev].is_cell_centered())
    Abort("Error: only support cell centered data!");

  nCell = loop_cell();

  nBrick = loop_brick();

  std::string outName = dirIn + ".dat";

  Print() << "Writing to " << outName << std::endl;

  outFile.open(outName.c_str(), std::ofstream::out | std::ofstream::trunc);

  //-----------Write header---------------
  outFile << "TITLE = " << '"' << outName << '"' << "\n";

  outFile << "VARIABLES = ";

  for (int i = 0; i < varNames.size(); ++i) {
    outFile << '"' << varNames[i] << '"';
    if (i != varNames.size() - 1) {
      outFile << ',' << " ";
    }
  }
  outFile << "\n";

  outFile << "ZONE "
          << " N=" << nCell << ", E=" << nBrick << ", F=FEPOINT, ET=BRICK"
          << "\n";
  //-----------------------------------------

  nCell = loop_cell(false);
  nBrick = loop_brick(false);

  if (outFile.is_open()) {
    outFile.close();
  }

  // Print() << "nCell = " << nCell << " nBrick = " << nBrick << std::endl;
}

int AMReXDataContainer::loop_cell(bool doCountOnly) {
  int iCount = 0;
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(iCell[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto& status = cell_status(iLev)[mfi].array();
      const Array4<Real>& cell = iCell[iLev][mfi].array();
      const Array4<Real>& data = mf[iLev][mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      for (int k = lo.z; k <= hi.z; ++k)
        for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i)
            if (!bit::is_refined(status(i, j, k))) {
              iCount++;
              cell(i, j, k) = iCount;
              if (!doCountOnly) {
                for (int iVar = 0; iVar < mf[iLev].nComp(); iVar++) {
                  outFile << data(i, j, k, iVar) << " ";
                }
                outFile << "\n";
              }
            }
    }

    iCell[iLev].FillBoundary();
  }

  for (int iLev = n_lev() - 2; iLev >= 0; iLev--) {
    fill_fine_lev_bny_cell_from_coarse(
        iCell[iLev], iCell[iLev + 1], 0, iCell[iLev].nComp(), ref_ratio[iLev],
        Geom(iLev), Geom(iLev + 1), cell_status(iLev + 1));
  }
  return iCount;
}

int AMReXDataContainer::loop_brick(bool doCountOnly) {
  int iBrick = 0;

  if (!doCountOnly)
    outFile.width(8);

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (MFIter mfi(iCell[iLev]); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const Array4<Real>& cell = iCell[iLev][mfi].array();
      const auto& status = cell_status(iLev)[mfi].array();

      const auto lo = lbound(box);
      const auto hi = ubound(box);

      // Loop over valid cells + low end ghost cells
      for (int k = lo.z - 1; k <= hi.z; ++k)
        for (int j = lo.y - 1; j <= hi.y; ++j)
          for (int i = lo.x - 1; i <= hi.x; ++i)
            if ((iLev > 0 && bit::is_lev_boundary(status(i, j, k))) ||
                (box.contains(i, j, k) && !bit::is_refined(status(i, j, k)))) {
              bool isBrick = true;

              for (int kk = k; kk <= k + 1; kk++)
                for (int jj = j; jj <= j + 1; jj++)
                  for (int ii = i; ii <= i + 1; ii++) {
                    if (cell(ii, jj, kk) == 0)
                      isBrick = false;
                  }

              if (isBrick) {
                iBrick++;

                if (!doCountOnly) {
                  outFile << (int)cell(i, j, k) << " ";
                  outFile << (int)cell(i + 1, j, k) << " ";
                  outFile << (int)cell(i + 1, j + 1, k) << " ";
                  outFile << (int)cell(i, j + 1, k) << " ";
                  outFile << (int)cell(i, j, k + 1) << " ";
                  outFile << (int)cell(i + 1, j, k + 1) << " ";
                  outFile << (int)cell(i + 1, j + 1, k + 1) << " ";
                  outFile << (int)cell(i, j + 1, k + 1) << "\n";
                }
              }
            }
    }
  }
  return iBrick;
}
