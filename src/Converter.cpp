#include <fstream>
#include <iostream>
#include <vector>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "Converter.h"
#include "GridUtility.h"

int main(int argc, char* argv[]) {
  amrex::Initialize(MPI_COMM_WORLD);

  std::vector<std::string> commandLine;
  for (int i = 0; i < argc; ++i) {
    commandLine.push_back((std::string)(argv[i]));
  }

  std::array<std::string, 2> arg = { "-h", "-help" };
  if (argc > 1 && find(arg.begin(), arg.end(), commandLine[1]) != arg.end()) {
    std::cout << " \n"
              << " This exectuable combines all the blocks in a *_amrex file "
                 "into one block. "
                 " \n\n Usage:\n "
                 " ./Converter.exe *_amrex\n\n";
    return 0;
  }

  for (std::vector<std::string>::size_type i = 1; i < commandLine.size(); i++) {
    std::cout << commandLine[i] << std::endl;
    Converter cv(commandLine[i]);
    cv.read();
    cv.convert();
    cv.write();
  }

  amrex::Finalize();

  return 1;
}

void Converter::read() {
  amrex::VisMF::Read(mf, dirIn + "/Level_0/Cell");
  amrex::Print() << "Reading in " << dirIn << std::endl;

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
    amrex::Real lo;
    HeaderFile >> lo;
    domainRange.setLo(i, lo);
  }

  for (int i = 0; i < nDim; ++i) {
    amrex::Real hi;
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

void Converter::convert() {
  amrex::BoxArray baIn = mf.boxArray();
  minBox = amrex::convert(baIn.minimalBox(), { 0, 0, 0 });

  amrex::BoxArray baOut;
  baOut.define(minBox);
  amrex::DistributionMapping dm;
  dm.define(baOut);

  distribute_FabArray(mf, baOut, dm, mf.nComp(), mf.nGrow(), true);

  int coord = 0;
  geom.define(minBox, &domainRange, coord);
}

void Converter::write() {
  std::string dirOut;
  std::size_t found = dirIn.find_last_of("_");
  if (found != std::string::npos) {
    dirOut = dirIn.substr(0, found) + "_single_block_amrex";
  } else {
    dirOut = "single_block_" + dirIn;
  }
  amrex::Print() << "Saving to " << dirOut << std::endl;

  amrex::WriteSingleLevelPlotfile(dirOut, mf, varNames, geom, time, 0);

  const std::string headerName = dirOut + "/FLEKSHeader";
  std::ofstream headerFile;
  headerFile.open(headerName.c_str(), std::ios::out | std::ios::trunc);
  if (!headerFile.good())
    amrex::FileOpenFailed(headerName);

  headerFile << plot_string << "\n";
  headerFile << rPlanet << "\n";
  headerFile.close();
}
