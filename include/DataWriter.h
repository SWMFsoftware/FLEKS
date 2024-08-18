#ifndef _DATAWRITER_H_
#define _DATAWRITER_H_

#include "DataContainer.h"

class DataWriter {
public:
  DataWriter(DataContainer* dcIn, const std::string& filenameIn) {
    dc = dcIn;
    filename = filenameIn;
  };

  virtual ~DataWriter() {};

  virtual void write() = 0;

  std::string type_string() { return fileTypeString.at(fType); }

  void print() {
    std::cout << "========DataWriter========\n";
    std::cout << "Write data to: " << filename << "\n";
    std::cout << "File type: " << type_string() << "\n";
    std::cout << "================================" << std::endl;
  }

protected:
  DataContainer* dc;
  std::string filename;
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

  void write() override {

    zoneType.set_type(dc->zone_type());

    size_t nCell = dc->count_cell();
    size_t nBrick = dc->count_zone();

    amrex::Vector<float> vars;
    vars.resize(nCell * dc->n_var());
    dc->get_cell(vars);

    amrex::Vector<size_t> zones;

    zones.resize(nBrick * zoneType.n_vertex());

    dc->get_zones(zones);

    outFile.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);

    //-----------Write header---------------
    outFile << "TITLE = " << '"' << filename << '"' << "\n";

    outFile << "VARIABLES = ";

    auto varNames = dc->var_names();
    for (int i = 0; i < varNames.size(); ++i) {
      outFile << '"' << varNames[i] << '"';
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

  void write() override {
    zoneType.set_type(dc->zone_type());

    size_t nCell = dc->count_cell();
    size_t nBrick = dc->count_zone();

    amrex::Vector<float> vars;
    vars.resize(nCell * dc->n_var());
    dc->get_cell(vars);

    amrex::Vector<float> xyz;
    // If nDim == 2, set the coordinates of the third dimension to 0.
    xyz.resize(3 * dc->n_var());
    dc->get_loc(xyz);

    //=== Brick data ===
    amrex::Vector<size_t> zones;
    amrex::Vector<int> bricksInt;
    zones.resize(nBrick * zoneType.n_vertex());
    dc->get_zones(zones);
    for (int i = 0; i < zones.size(); ++i) {
      bricksInt.push_back(zones[i] - 1);
    }
    int* brickData = bricksInt.data();
    //===================

    amrex::Vector<int> brickType(nBrick, zoneType.vtk_index());

    // All variables are scalars, so vardim is 1.
    amrex::Vector<int> vardim(dc->n_var(), 1);

    // Cell-based: 0  (This is connectivity cell, NOT simulation cell)
    // Point-based: 1
    amrex::Vector<int> centering(dc->n_var(), 1);

    char** varnames;
    varnames = new char*[dc->n_var()];
    for (int i = 0; i < dc->n_var(); ++i) {
      varnames[i] = new char[dc->var_names()[i].size() + 1];
      strcpy(varnames[i], dc->var_names()[i].c_str());
    }

    float** v;
    v = new float*[dc->n_var()];
    for (int i = 0; i < dc->n_var(); ++i) {
      v[i] = new float[nCell];
      for (size_t j = 0; j < nCell; ++j) {
        v[i][j] = vars[j * dc->n_var() + i];
      }
    }

    write_unstructured_mesh(filename.c_str(), saveBinary, nCell, xyz.data(),
                            nBrick, brickType.data(), brickData, dc->n_var(),
                            vardim.data(), centering.data(), varnames, v);

    for (int i = 0; i < dc->n_var(); ++i) {
      delete[] varnames[i];
      delete[] v[i];
    }
    delete[] varnames;
    delete[] v;
  }

private:
  bool saveBinary;
};

#endif