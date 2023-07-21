
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

  VisMF::Read(mf, dirIn + "/Level_0/Cell");
}

void AMReXDataContainer::write() {
  if (!mf.is_cell_centered())
    Abort("Error: only support cell centered data!");

  iCell.define(mf.boxArray(), mf.DistributionMap(), 1, 1);
  iCell.setVal(0);

  nCell = loop_cell();

  iCell.FillBoundary();

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
  for (MFIter mfi(iCell); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<int>& cell = iCell[mfi].array();
    const Array4<Real>& data = mf[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          iCount++;
          cell(i, j, k) = iCount;
          if (!doCountOnly) {
            for (int iVar = 0; iVar < mf.nComp(); iVar++) {
              outFile << data(i, j, k, iVar) << " ";
            }
            outFile << "\n";
          }
        }
  }
  return iCount;
}

int AMReXDataContainer::loop_brick(bool doCountOnly) {
  int iBrick = 0;

  if (!doCountOnly)
    outFile.width(8);

  for (MFIter mfi(iCell); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<int>& cell = iCell[mfi].array();

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
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
              outFile << cell(i, j, k) << " ";
              outFile << cell(i + 1, j, k) << " ";
              outFile << cell(i + 1, j + 1, k) << " ";
              outFile << cell(i, j + 1, k) << " ";
              outFile << cell(i, j, k + 1) << " ";
              outFile << cell(i + 1, j, k + 1) << " ";
              outFile << cell(i + 1, j + 1, k + 1) << " ";
              outFile << cell(i, j + 1, k + 1) << "\n";
            }
          }
        }
  }
  return iBrick;
}
