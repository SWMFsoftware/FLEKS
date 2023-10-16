#include <fstream>
#include <iostream>
#include <vector>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "Converter.h"
#include "GridUtility.h"

using namespace amrex;
int main(int argc, char* argv[]) {
  std::vector<std::string> cdl;
  for (int i = 0; i < argc; ++i) {
    cdl.push_back((std::string)(argv[i]));
  }

  std::vector<std::string> fileNames;
  FileType sType = FileType::UNSET;
  FileType dType = FileType::UNSET;

  int nSmooth = 0;

  size_t i = 1;
  while (i < cdl.size()) {
    if (cdl[i] == "-f") {
      i++;
      while (i < cdl.size() && cdl[i][0] != '-') {
        fileNames.push_back(cdl[i++]);
      }

      if (fileNames.empty()) {
        std::cout << "Error: -f option requires an argument.\n";
        return 0;
      }
    } else if (cdl[i] == "-h") {
      i++;

      printf("Convert FLEKS/BATSRUS data to other formats.\n\n");
      printf(" Usage:\n ");
      printf(
          " ./Converter.exe -f filename -d dest_format -s source_format\n\n");

      printf(" Options:\n");
      printf("  -h        : Print help message.\n");
      printf("  -f        : Specify the file name to convert. Multiple files "
             "can be converted at a time.\n");
      printf("  -d        : Specify the destination file format.\n");
      printf("               Options : VTK, TEC\n");
      printf("  -s [optional]: Specify the source file format.\n");
      printf("               Options: AMReX, IDL\n");
      printf("  -smooth n : Smooth the data n times\n");

      printf("\n");

      printf(" Examples:\n");
      printf("  ./Converter.exe -f f1_amrex f2_amrex -d VTK\n");
      printf("  ./Converter.exe -f 3d*_amrex -d TEC\n");
      printf("  ./Converter.exe -f 3d*_amrex -d TEC -smooth 3\n");

      return 0;
    } else if (cdl[i] == "-s") {
      i++;
      if (i >= cdl.size()) {
        std::cout << "Error: -s option requires an argument.\n";
        return 0;
      } else {
        sType = stringToFileType.at(cdl[i++]);
      }
    } else if (cdl[i] == "-d") {
      i++;
      if (i >= cdl.size()) {
        std::cout << "Error: -d option requires an argument.\n";
        return 0;
      } else {
        dType = stringToFileType.at(cdl[i++]);
        std::cout << "dtype = " << cdl[i - 1] << " dtype = " << (int)dType
                  << "\n";
      }
    } else if (cdl[i] == "-smooth") {
      i++;
      if (i >= cdl.size()) {
        std::cout << "Error: -smooth option requires an argument.\n";
        return 0;
      } else {
        nSmooth = std::stoi(cdl[i++]);
      }
    }
  }

  if (dType == FileType::UNSET) {
    std::cout << "Error: destination file format is required! Set the format "
                 "with -d option.\n";
    return 0;
  }

  Initialize(MPI_COMM_WORLD);

  for (size_t i = 0; i < fileNames.size(); i++) {
    Converter cv(fileNames[i], sType, dType);
    cv.read();
    if (nSmooth > 0) {
      cv.smooth(nSmooth);
    }
    cv.write();
  }

  Finalize();

  return 1;
}
