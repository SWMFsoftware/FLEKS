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

void print_help() {
  std::cout << "Convert FLEKS and BATSRUS simulation data between formats.\n\n"
            << "Usage:\n"
            << "  ./bin/converter.exe -f <file...> -d <dest_format> [-s "
               "<source_format>] [options]\n\n"
            << "Required options:\n"
            << "  -f <file...>       Specify input file(s) or directory(s) to "
               "convert.\n"
            << "                     Multiple files can be specified.\n"
            << "  -d <format>        Destination file format: VTK, TEC, VTM\n\n"
            << "Optional options:\n"
            << "  -s <format>        Source file format: AMReX, IDL, TEC "
               "(ASCII .dat)\n"
            << "                     Inferred automatically from extension if "
               "omitted:\n"
            << "                       *_amrex -> AMReX\n"
            << "                       *.out   -> IDL\n"
            << "                       *.dat   -> TEC\n"
            << "  -D                 Delete the source file(s) after "
               "successful conversion.\n"
            << "  -smooth <n>        Smooth the data n times (nonnegative "
               "integer; AMReX only).\n"
            << "  -h, --help         Print this help message.\n\n"
            << "Notes:\n"
            << "  - Format names are case-sensitive.\n"
            << "  - Output paths append an extension to the input path: "
               "<file>.vtk, <file>.dat, or <file>.vtm.\n"
            << "  - Output will not overwrite or alias the source file.\n\n"
            << "Examples:\n"
            << "  ./bin/converter.exe -f 3d.dat -d VTK\n"
            << "  ./bin/converter.exe -f 3d.dat -d VTM\n"
            << "  ./bin/converter.exe -f 3d.dat -s TEC -d VTK\n"
            << "  ./bin/converter.exe -f first.dat second.dat -d VTK\n"
            << "  ./bin/converter.exe -f 3d.dat -d TEC\n"
            << "  ./bin/converter.exe -f f1_amrex f2_amrex -d VTK\n"
            << "  ./bin/converter.exe -f 3d*_amrex -d TEC -smooth 3\n"
            << "  ./bin/converter.exe -f 3d.out -d VTK\n";
}

int main(int argc, char* argv[]) {
  if (argc <= 1) {
    print_help();
    return 0;
  }

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
    } else if (cdl[i] == "-h" || cdl[i] == "--help") {
      print_help();
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
            (type->second != FileType::VTK && type->second != FileType::TEC &&
             type->second != FileType::VTM)) {
          std::cerr << "Error: destination format must be VTK, TEC, or VTM.\n";
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
