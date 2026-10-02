#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "Converter.h"

using namespace amrex;

extern "C" {
void timing_start_c(size_t*, char*) {}
void timing_stop_c(size_t*, char*) {}
}

int main(int argc, char* argv[]) {
  const std::vector<std::string> cdl(argv, argv + argc);
  std::vector<std::string> fileNames;
  FileType sType = FileType::UNSET;
  FileType dType = FileType::UNSET;

  bool deleteSource = false;

  int nSmooth = 0;

  size_t i = 1;
  while (i < cdl.size()) {
    if (cdl[i] == "-f") {
      ++i;
      const size_t first = i;
      while (i < cdl.size() && !cdl[i].empty() && cdl[i][0] != '-') {
        fileNames.push_back(cdl[i++]);
      }

      if (i == first) {
        std::cout << "Error: -f option requires an argument.\n";
        return EXIT_FAILURE;
      }
    } else if (cdl[i] == "-h") {
      ++i;

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
      printf("               Options: AMReX, IDL, TEC (ASCII .dat)\n");
      printf("  -D        : Delete the source files\n");
      printf("  -smooth n : Smooth the data n times\n");

      printf("\n");

      printf(" Examples:\n");
      printf("  ./Converter.exe -f f1_amrex f2_amrex -d VTK\n");
      printf("  ./Converter.exe -f 3d*_amrex -d TEC\n");
      printf("  ./Converter.exe -f 3d*_amrex -d TEC -smooth 3\n");
      printf("  ./Converter.exe -f 3d.dat -s TEC -d VTK\n");

      return 0;
    } else if (cdl[i] == "-s") {
      ++i;
      if (i >= cdl.size()) {
        std::cout << "Error: -s option requires an argument.\n";
        return EXIT_FAILURE;
      } else {
        const auto type = stringToFileType.find(cdl[i++]);
        if (type == stringToFileType.end() ||
            (type->second != FileType::AMREX && type->second != FileType::IDL &&
             type->second != FileType::TEC)) {
          std::cerr << "Error: source format must be AMReX, IDL, or TEC.\n";
          return EXIT_FAILURE;
        }
        sType = type->second;
      }
    } else if (cdl[i] == "-d") {
      ++i;
      if (i >= cdl.size()) {
        std::cout << "Error: -d option requires an argument.\n";
        return EXIT_FAILURE;
      } else {
        const auto type = stringToFileType.find(cdl[i++]);
        if (type == stringToFileType.end() ||
            (type->second != FileType::VTK && type->second != FileType::TEC)) {
          std::cerr << "Error: destination format must be VTK or TEC.\n";
          return EXIT_FAILURE;
        }
        dType = type->second;
      }
    } else if (cdl[i] == "-D") {
      ++i;
      deleteSource = true;
    } else if (cdl[i] == "-smooth") {
      ++i;
      if (i >= cdl.size()) {
        std::cout << "Error: -smooth option requires an argument.\n";
        return EXIT_FAILURE;
      } else {
        try {
          size_t parsed = 0;
          nSmooth = std::stoi(cdl[i], &parsed);
          if (parsed != cdl[i].size() || nSmooth < 0)
            throw std::invalid_argument("invalid smoothing count");
          ++i;
        } catch (const std::exception&) {
          std::cerr << "Error: -smooth requires a nonnegative integer.\n";
          return EXIT_FAILURE;
        }
      }
    } else {
      std::cerr << "Error: unknown option: " << cdl[i] << "\n";
      return EXIT_FAILURE;
    }
  }

  if (dType == FileType::UNSET) {
    std::cout << "Error: destination file format is required! Set the format "
                 "with -d option.\n";
    return EXIT_FAILURE;
  }
  if (fileNames.empty()) {
    std::cerr << "Error: specify input files with -f.\n";
    return EXIT_FAILURE;
  }

  Initialize(MPI_COMM_WORLD);

  int status = EXIT_SUCCESS;
  for (const auto& filename : fileNames) {
    const bool isTec = sType == FileType::TEC ||
                       (sType == FileType::UNSET &&
                        std::filesystem::path(filename).extension() == ".dat");
    if (isTec && nSmooth > 0) {
      std::cerr << "Error: -smooth is not supported for TEC input: " << filename
                << "\n";
      status = EXIT_FAILURE;
      continue;
    }
    try {
      Converter cv(filename, sType, dType);
      if (cv.read() == iFail) {
        std::cerr << "Error: reading file failed: " << filename << "\n";
        status = EXIT_FAILURE;
        continue;
      }
      if (nSmooth > 0)
        cv.smooth(nSmooth);
      if (cv.write() == iFail) {
        std::cerr << "Error: writing file failed: " << filename << "\n";
        status = EXIT_FAILURE;
        continue;
      }
      if (deleteSource) {
        std::cout << "Deleting source file: " << filename << "\n";
        std::error_code error;
        std::filesystem::remove_all(filename, error);
        if (error) {
          std::cerr << "Error: deleting source failed: " << error.message()
                    << "\n";
          status = EXIT_FAILURE;
        }
      }
    } catch (const std::exception& error) {
      std::cerr << "Error: converting " << filename << ": " << error.what()
                << "\n";
      status = EXIT_FAILURE;
    }
  }

  Finalize();

  return status;
}
