#ifndef _DATAWRITER_H_
#define _DATAWRITER_H_

#include <cassert>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>

#include "DataContainer.h"

class DataWriter {
public:
  DataWriter(DataContainer* dcIn, const std::string& filenameIn) {
    dc = dcIn;
    filename = filenameIn;
    sourceFilename = filenameIn;
  };

  virtual ~DataWriter() {};

  virtual int write() = 0;

  std::string type_string() { return fileTypeString.at(fType); }

  void print() {
    std::cout << "========DataWriter========\n";
    std::cout << "Write data to: " << filename << "\n";
    std::cout << "File type: " << type_string() << "\n";
    std::cout << "================================" << std::endl;
  }

protected:
  bool can_write() const {
    std::error_code error;
    // weakly_canonical resolves symlinks and normalises without requiring
    // that the output file already exists (unlike std::filesystem::equivalent).
    const auto src = std::filesystem::weakly_canonical(sourceFilename, error);
    const auto dst = std::filesystem::weakly_canonical(filename, error);
    if (src == dst) {
      std::cerr << "Error: output aliases the input file: " << filename << "\n";
      return false;
    }
    return true;
  }

  DataContainer* dc;
  std::string filename;
  std::string sourceFilename;
  FileType fType;
  std::ofstream outFile;

  ZoneType zoneType;
};

class TECWriter : public DataWriter {
public:
  TECWriter(DataContainer* dcIn, const std::string& filenameIn)
      : DataWriter(dcIn, filenameIn) {
    fType = FileType::TEC;
    filename = filenameIn + ".dat";
  };

  ~TECWriter() {};

  int write() override {

    if (!can_write())
      return iFail;

    zoneType.set_type(dc->zone_type());
    if (zoneType.n_vertex() <= 0) {
      std::cerr << "Error: unsupported zone topology.\n";
      return iFail;
    }

    size_t nCell = dc->count_cell();
    size_t nBrick = dc->count_zone();

    amrex::Vector<float> vars;
    vars.resize(nCell * dc->n_var());
    dc->get_cell(vars);

    amrex::Vector<size_t> zones;

    zones.resize(nBrick * zoneType.n_vertex());

    dc->get_zones(zones);

    outFile.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
    if (!outFile) {
      std::cerr << "Error: cannot open output: " << filename << "\n";
      return iFail;
    }
    outFile.precision(std::numeric_limits<float>::max_digits10);

    //-----------Write header---------------
    outFile << "TITLE = " << std::quoted(filename) << "\n";

    outFile << "VARIABLES = ";

    auto varNames = dc->var_names();
    for (int i = 0; i < varNames.size(); ++i) {
      outFile << std::quoted(varNames[i]);
      if (i != varNames.size() - 1) {
        outFile << ',' << " ";
      }
    }
    outFile << "\n";

    outFile << "ZONE "
            << " N=" << nCell << ", E=" << nBrick
            << ", F=FEPOINT, ET=" << zoneType.tec_string() << "\n";
    //-----------------------------------------

    // Write cell data
    for (size_t i = 0; i < nCell; ++i) {
      for (int j = 0; j < dc->n_var(); ++j) {
        outFile << vars[i * dc->n_var() + j] << " ";
      }
      outFile << "\n";
    }

    // Write zone data
    for (size_t i = 0; i < nBrick; ++i) {
      for (int j = 0; j < zoneType.n_vertex(); ++j) {
        outFile << zones[i * zoneType.n_vertex() + j] << " ";
      }
      outFile << "\n";
    }

    if (outFile.is_open()) {
      outFile.close();
    }

    if (!outFile) {
      std::cerr << "Error: writing output failed: " << filename << "\n";
      return iFail;
    }
    return iSuccess;
  }
};

class VTKWriter : public DataWriter {
public:
  VTKWriter(DataContainer* dcIn, const std::string& filenameIn)
      : DataWriter(dcIn, filenameIn) {
    fType = FileType::VTK;
    filename = filenameIn + ".vtk";

    saveBinary = true;
  };

  ~VTKWriter() {};

  int write() override {
    if (!can_write())
      return iFail;
    zoneType.set_type(dc->zone_type());

    size_t nCell = dc->count_cell();
    size_t nBrick = dc->count_zone();
    const int nVertex = zoneType.n_vertex();
    const size_t maxInt = std::numeric_limits<int>::max();
    if (nVertex <= 0 || nCell > maxInt / 3 ||
        nBrick > maxInt / static_cast<size_t>(nVertex + 1)) {
      std::cerr << "Error: mesh exceeds legacy VTK integer limits.\n";
      return iFail;
    }

    amrex::Vector<float> vars;
    vars.resize(nCell * dc->n_var());
    dc->get_cell(vars);

    amrex::Vector<float> xyz;
    // If nDim == 2, set the coordinates of the third dimension to 0.
    xyz.resize(3 * nCell);
    dc->get_loc(xyz);
    assert(xyz.size() == 3 * nCell);

    //=== Brick data ===
    amrex::Vector<size_t> zones;
    amrex::Vector<int> bricksInt;
    zones.resize(nBrick * zoneType.n_vertex());
    dc->get_zones(zones);
    bricksInt.reserve(zones.size());
    for (const auto index : zones) {
      if (index == 0 || index > nCell) {
        std::cerr << "Error: invalid mesh connectivity.\n";
        return iFail;
      }
      bricksInt.push_back(static_cast<int>(index - 1));
    }
    int* brickData = bricksInt.data();
    //===================

    amrex::Vector<int> brickType(nBrick, zoneType.vtk_index());

    // All variables are scalars, so vardim is 1.
    amrex::Vector<int> vardim(dc->n_var(), 1);

    // Cell-based: 0  (This is connectivity cell, NOT simulation cell)
    // Point-based: 1
    amrex::Vector<int> centering(dc->n_var(), 1);

    amrex::Vector<amrex::Vector<char> > varnameStorage(dc->n_var());
    amrex::Vector<char*> varnames(dc->n_var());
    const auto varNames = dc->var_names();
    std::set<std::string> usedNames;
    for (int i = 0; i < dc->n_var(); ++i) {
      std::string base = varNames[i];
      for (auto& ch : base)
        if (std::isspace(static_cast<unsigned char>(ch)))
          ch = '_';
      if (base.empty())
        base = "variable";
      std::string name = base;
      for (size_t suffix = 2; !usedNames.insert(name).second; ++suffix)
        name = base + "_" + std::to_string(suffix);
      if (name != varNames[i])
        std::cout << "VTK variable: \"" << varNames[i] << "\" -> \"" << name
                  << "\"\n";
      varnameStorage[i].assign(name.begin(), name.end());
      varnameStorage[i].push_back('\0');
      varnames[i] = varnameStorage[i].data();
    }

    amrex::Vector<amrex::Vector<float> > valueStorage(
        dc->n_var(), amrex::Vector<float>(nCell));
    amrex::Vector<float*> v(dc->n_var());
    for (int i = 0; i < dc->n_var(); ++i) {
      v[i] = valueStorage[i].data();
      for (size_t j = 0; j < nCell; ++j) {
        v[i][j] = vars[j * dc->n_var() + i];
      }
    }

    if (!write_unstructured_mesh(
            filename.c_str(), saveBinary, static_cast<int>(nCell), xyz.data(),
            static_cast<int>(nBrick), brickType.data(), brickData, dc->n_var(),
            vardim.data(), centering.data(), varnames.data(), v.data())) {
      std::cerr << "Error: writing output failed: " << filename << "\n";
      return iFail;
    }

    return iSuccess;
  }

private:
  bool saveBinary;
};

#endif
