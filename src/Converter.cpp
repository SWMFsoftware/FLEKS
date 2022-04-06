#include <fstream>
#include <iostream>
#include <vector>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "Converter.h"
#include "GridUtility.h"

using namespace amrex;
int main(int argc, char* argv[]) {
  Initialize(MPI_COMM_WORLD);

  std::vector<std::string> commandLine;
  for (int i = 0; i < argc; ++i) {
    commandLine.push_back((std::string)(argv[i]));
  }

  std::array<std::string, 2> arg = { "-h", "-help" };
  if (argc > 1 && find(arg.begin(), arg.end(), commandLine[1]) != arg.end()) {
    std::cout << " \n"
              << " This exectuable converts a *_amrex file to a tecplot ascii "
                 "*.dat file. "
                 " \n\n Usage:\n "
                 " ./Converter.exe *_amrex\n\n";
    return 0;
  }

  for (std::vector<std::string>::size_type i = 1; i < commandLine.size(); i++) {
    std::cout << commandLine[i] << std::endl;
    Converter cv(commandLine[i]);
    cv.read();
    cv.write();
  }

  Finalize();

  return 1;
}

void Converter::read() {
  VisMF::Read(mf, dirIn + "/Level_0/Cell");
  Print() << "Reading in " << dirIn << std::endl;

  std::string headerName = dirIn + "/Header";
  std::ifstream HeaderFile;
  HeaderFile.open(headerName, std::ifstream::in);

  HeaderFile.precision(17);

  std::string versionName;
  HeaderFile >> versionName;

  int nVar;
  HeaderFile >> nVar;

  for (int ivar = 0; ivar < nVar; ++ivar) {
    std::string var;
    HeaderFile >> var;
    varNames.push_back(var);
  }

  int nDim;
  HeaderFile >> nDim;

  HeaderFile >> time;

  int nLevel;
  HeaderFile >> nLevel;
  for (int i = 0; i < nDim; ++i) {
    Real lo;
    HeaderFile >> lo;
    domainRange.setLo(i, lo);
  }

  for (int i = 0; i < nDim; ++i) {
    Real hi;
    HeaderFile >> hi;
    domainRange.setHi(i, hi);
  }

  HeaderFile.close();

  headerName = dirIn + "/FLEKSHeader";
  HeaderFile.open(headerName, std::ifstream::in);
  std::getline(HeaderFile, plot_string);
  HeaderFile >> rPlanet;
  HeaderFile.close();
}

void Converter::write() {
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

int Converter::loop_cell(bool doCountOnly) {
  int iCount = 0;
  for (MFIter mfi(iCell); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
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

int Converter::loop_brick(bool doCountOnly) {
  int iBrick = 0;

  if (!doCountOnly)
    outFile.width(8);

  for (MFIter mfi(iCell); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
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
