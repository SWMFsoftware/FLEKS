#ifndef _DATACONTAINER_H_
#define _DATACONTAINER_H_

#include <cassert>
#include <string>

#include <AMReX_Box.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include "Grid.h"
#include "VisitWriter.h"

enum class FileType {
  AMREX = 0,
  IDL,
  TECPLOT,
  VTK,
  HDF5,
  ADIOS2,
  UNSET,
  UNKNOWN
};

static const std::map<FileType, std::string> fileTypeString = {
  { FileType::AMREX, "AMReX" },     { FileType::IDL, "IDL" },
  { FileType::TECPLOT, "TECPLOT" }, { FileType::VTK, "VTK" },
  { FileType::HDF5, "HDF5" },       { FileType::ADIOS2, "ADIOS2" },
  { FileType::UNSET, "UNSET" },     { FileType::UNKNOWN, "UNKNOWN" }
};

static const std::map<std::string, FileType> stringToFileType = {
  { "AMReX", FileType::AMREX },     { "IDL", FileType::IDL },
  { "TECPLOT", FileType::TECPLOT }, { "VTK", FileType::VTK },
  { "HDF5", FileType::HDF5 },       { "ADIOS2", FileType::ADIOS2 },
  { "UNSET", FileType::UNSET },     { "UNKNOWN", FileType::UNKNOWN }
};

// Virtual base
class DataContainer {
public:
  virtual ~DataContainer() = default;

  virtual void read() = 0;

  virtual size_t count_cell() = 0;

  virtual size_t count_zone() = 0;

  virtual void get_cell(amrex::Vector<float>& vars) = 0;

  virtual void get_zones(amrex::Vector<size_t>& zones) = 0;

  virtual void get_loc(amrex::Vector<float>& vars) = 0;

  std::string type_string() { return fileTypeString.at(dataType); }

  int n_var() { return nVar; }

  int n_dim() { return nDim; }

  amrex::Vector<std::string> var_names() { return varNames; }

  void print() {
    std::cout << "========DataContainer========\n";
    std::cout << "Read data from: " << dirIn << "\n";
    std::cout << "Data type: " << type_string() << "\n";
    std::cout << "================================" << std::endl;
  }

protected:
  std::string dirIn;
  FileType dataType;

  size_t nCell;
  size_t nBrick;
  int nVar;
  int nDim;
  int iter;
  amrex::Real time;

  amrex::Vector<std::string> varNames;

  amrex::Real rPlanet;
};

class IDLDataContainer : public DataContainer {

  enum class IDLFileType { ASCII = 0, REAL4, REAL8 };

private:
  static const int i4 = 4;
  static const int i8 = 8;
  IDLFileType idlType;
  int nReal;

  std::string unit;

  int nParam;

  int nSize[3] = { 1, 1, 1 };
  std::vector<double> param_I;
  std::vector<std::string> paramName_I;
  double** data_II = nullptr;

  amrex::BaseFab<float> fab;

  amrex::BaseFab<int> iCell;

public:
  IDLDataContainer(const std::string& in) {
    dirIn = in;
    dataType = FileType::IDL;
    idlType = get_file_type();
  }

  ~IDLDataContainer() {
    if (data_II != nullptr) {
      for (int i = 0; i < nCell; ++i) {
        delete[] data_II[i];
      }
      delete[] data_II;
    }
  };

  IDLFileType get_file_type() {
    std::ifstream inFile;
    inFile.open(dirIn.c_str(), std::ifstream::in | std::ifstream::binary);

    int lenHead;
    static_assert(i4 == sizeof(int), "Error: the size of integer is not 4!");
    static_assert(i4 == sizeof(float), "Error: the size of float is not 4!");
    static_assert(i8 == sizeof(double), "Error: the size of double is not 8!");

    IDLFileType t;
    inFile.read(reinterpret_cast<char*>(&lenHead), i4);

    if (lenHead != 79 && lenHead != 500) {
      t = IDLFileType::ASCII;
    } else {
      // Read the length of the second line;
      char* tmp;
      tmp = new char[lenHead + 4];
      // Skip the rest of the first line;
      inFile.read(tmp, lenHead + 4);
      delete[] tmp;

      int len;
      inFile.read(reinterpret_cast<char*>(&len), i4);

      switch (len) {
        case 20:
          t = IDLFileType::REAL4;
          break;
        case 24:
          t = IDLFileType::REAL8;
          break;
        default:
          abort();
      }
    }

    if (inFile.is_open()) {
      inFile.close();
    }

    return t;
  }

  void read() override {
    if (idlType == IDLFileType::ASCII) {
      read_ascii();
    } else if (idlType == IDLFileType::REAL4) {
      nReal = 4;
      read_binary<float>();
    } else if (idlType == IDLFileType::REAL8) {
      nReal = 8;
      read_binary<double>();
    }
  };

  size_t count_cell() override { return nCell; }

  size_t count_zone() override {
    int n = 1;
    for (int i = 0; i < nDim; ++i) {
      n *= nSize[i] - 1;
    }
    return n;
  }

  void get_cell(amrex::Vector<float>& vars) override { loop_cell(vars, false); }

  void get_loc(amrex::Vector<float>& vars) override { loop_cell(vars, true); }

  size_t loop_cell(amrex::Vector<float>& vars, bool doStoreLoc) {
    std::string funcName = "IDLDataContainer::loop_cell()";
    BL_PROFILE(funcName);

    vars.clear();

    const amrex::Box& box = fab.box();
    const amrex::Array4<float>& data = fab.array();

    const amrex::Array4<int>& cell = iCell.array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int iCount = 0;
    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          iCount++;
          cell(i, j, k) = iCount;
          if (doStoreLoc) {
            for (int iDim = 0; iDim < 3; iDim++) {
              if (iDim < nDim) {
                vars.push_back(data(i, j, k, iDim));
              } else {
                vars.push_back(0.0);
              }
            }
          } else {
            for (int iVar = 0; iVar < fab.nComp(); iVar++)
              vars.push_back(data(i, j, k, iVar));
          }
        }
  }

  void get_zones(amrex::Vector<size_t>& zones) override {
    std::string funcName = "IDLDataContainer::get_zone()";
    BL_PROFILE(funcName);

    const amrex::Box& box = iCell.box();
    const amrex::Array4<int>& cell = iCell.array();

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    zones.clear();
    if (nDim == 2) {
      for (int j = lo.y; j <= hi.y - 1; ++j)
        for (int i = lo.x; i <= hi.x - 1; ++i) {
          int k = lo.z;
          zones.push_back(cell(i, j, k));
          zones.push_back(cell(i + 1, j, k));
          zones.push_back(cell(i + 1, j + 1, k));
          zones.push_back(cell(i, j + 1, k));
        }
    } else if (nDim == 3) {
      for (int k = lo.z; k <= hi.z - 1; ++k)
        for (int j = lo.y; j <= hi.y - 1; ++j)
          for (int i = lo.x; i <= hi.x - 1; ++i) {
            zones.push_back(cell(i, j, k));
            zones.push_back(cell(i + 1, j, k));
            zones.push_back(cell(i + 1, j + 1, k));
            zones.push_back(cell(i, j + 1, k));
            zones.push_back(cell(i, j, k + 1));
            zones.push_back(cell(i + 1, j, k + 1));
            zones.push_back(cell(i + 1, j + 1, k + 1));
            zones.push_back(cell(i, j + 1, k + 1));
          }
    }
  }

  //=================================================

  template <typename real> void read_binary() {
    std::ifstream inFile;
    inFile.open(dirIn.c_str(), std::ifstream::in | std::ifstream::binary);

    // Lambda reads 4 bytes integer;
    auto read_int = [&]() {
      int nRec;
      inFile.read(reinterpret_cast<char*>(&nRec), i4);
      return nRec;
    };

    // Lambda expression reads real4 or real8
    auto read_float = [&]() {
      real f;
      inFile.read(reinterpret_cast<char*>(&f), nReal);
      return f;
    };

    int nRec;
    { // Read unit.

      nRec = read_int();
      auto* unitTmp = new char[nRec];
      inFile.read(unitTmp, nRec);

      std::string sTmp(unitTmp);
      std::stringstream ss;
      ss << sTmp;
      ss >> unit;
      delete[] unitTmp;

      nRec = read_int();
    }

    {
      nRec = read_int();
      iter = read_int();
      time = read_float();
      nDim = read_int();
      nParam = read_int();
      nVar = read_int();
      nRec = read_int();

      nVar += nDim;
    }

    if (nDim <= 0) {
      std::cout << "Error: unstructured grid is not supported yet!"
                << std::endl;
      exit(1);
    }

    {
      nCell = 0;
      nRec = read_int();
      for (int iDim = 0; iDim < nDim; ++iDim) {
        int size = read_int();
        nSize[iDim] = size;
        nCell = nCell > 0 ? nCell * size : size;
      }
      nRec = read_int();

      amrex::Box bx(amrex::IntVect(0),
                    amrex::IntVect(nSize[0] - 1, nSize[1] - 1, nSize[2] - 1));
      fab.clear();
      fab.resize(bx, nVar);

      iCell.clear();
      iCell.resize(bx, 1);
    }

    {
      nRec = read_int();
      for (int i = 0; i < nParam; ++i) {
        param_I.push_back(read_float());
      }
      nRec = read_int();
    }

    {
      nRec = read_int();

      char* tmp;
      tmp = new char[nRec];
      inFile.read(tmp, nRec);
      std::stringstream ss;
      ss << tmp;
      delete[] tmp;

      for (int i = 0; i < nVar; ++i) {
        std::string sVar;
        ss >> sVar;
        varNames.push_back(sVar);
      }

      for (int i = 0; i < nParam; ++i) {
        std::string sVar;
        ss >> sVar;
        paramName_I.push_back(sVar);
      }

      nRec = read_int();
    }

    if (data_II == nullptr) {
      data_II = new double*[nCell];
      for (int i = 0; i < nCell; ++i) {
        data_II[i] = new double[nVar];
      }
    }

    const amrex::Box box = fab.box();

    const auto& data = fab.array();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    {
      // Read coordinates.
      nRec = read_int();
      assert(nRec / nDim / nReal == nCell);

      real* x;
      x = new real[nCell];
      for (int iDim = 0; iDim < nDim; ++iDim) {
        inFile.read(reinterpret_cast<char*>(x), nRec / nDim);
        int iPoint = 0;
        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              data(i, j, k, iDim) = x[iPoint];
              iPoint++;
            }
      }

      delete[] x;
      nRec = read_int();
    }

    {
      real* w;
      w = new real[nCell];

      for (int iVar = nDim; iVar < nVar; ++iVar) {
        nRec = read_int();
        inFile.read(reinterpret_cast<char*>(w), nRec);

        int iPoint = 0;
        for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
              data(i, j, k, iVar) = w[iPoint];
              iPoint++;
            }
        nRec = read_int();
      }
      delete[] w;
    }

    if (inFile.is_open()) {
      inFile.close();
    }
  }

  void read_ascii() {
    std::ifstream inFile;
    inFile.open(dirIn.c_str(), std::ifstream::in);
    inFile >> unit >> iter >> time >> nDim >> nParam >> nVar;

    nVar += nDim;

    printf("read_ascii");

    if (nDim <= 0) {
      std::cout << "Error: unstructured grid is not supported yet!"
                << std::endl;
      exit(1);
    }

    nCell = 0;
    for (int i = 0; i < nDim; ++i) {
      int size;
      inFile >> size;
      nSize[i] = size;
      nCell = nCell > 0 ? nCell * size : size;
    }

    amrex::Box bx(amrex::IntVect(0),
                  amrex::IntVect(nSize[0] - 1, nSize[1] - 1, nSize[2] - 1));
    fab.clear();
    fab.resize(bx, nVar);

    iCell.clear();
    iCell.resize(bx, 1);

    // Read parameters.
    for (int i = 0; i < nParam; ++i) {
      double param;
      inFile >> param;
      param_I.push_back(param);
    }

    // Read variable name.
    for (int i = 0; i < nVar; ++i) {
      std::string s;
      inFile >> s;
      varNames.push_back(s);
    }

    // Read parameter name.
    for (int i = 0; i < nParam; ++i) {
      std::string s;
      inFile >> s;
      paramName_I.push_back(s);
    }

    if (data_II == nullptr) {
      data_II = new double*[nCell];
      for (int i = 0; i < nCell; ++i) {
        data_II[i] = new double[nVar];
      }
    }

    const amrex::Box box = fab.box();

    const auto& data = fab.array();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          for (int iVar = 0; iVar < nVar; ++iVar) {
            float f;
            inFile >> f;
            data(i, j, k, iVar) = f;
            // printf("data(%d, %d, %d, %d) = %f\n", i, j, k, iVar, f);
          }
        }

    if (inFile.is_open()) {
      inFile.close();
    }
  }
  //===============================================================
};

class AMReXDataContainer : public Grid, public DataContainer {
private:
  amrex::Vector<amrex::MultiFab> mf;
  amrex::Vector<amrex::MultiFab> iCell;

  amrex::RealBox domain;
  amrex::Box cellBox;

public:
  AMReXDataContainer(const std::string& in, const amrex::Geometry& gm,
                     const amrex::AmrInfo& amrInfo)
      : Grid(gm, amrInfo, 2), DataContainer() {
    dirIn = in;
    dataType = FileType::AMREX;

    mf.resize(n_lev_max());
    iCell.resize(n_lev_max());
  }
  ~AMReXDataContainer(){};

  static void read_header(std::string& headerName, int& nVar, int& nDim,
                          amrex::Real& time, int& finest_level,
                          amrex::RealBox& domain, amrex::Box& cellBox,
                          amrex::Vector<std::string>& varNames);
  void read_header();

  void post_regrid() override {
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      distribute_FabArray(iCell[iLev], cGrids[iLev], DistributionMap(iLev), 1,
                          nGst, false);
    }
    distribute_grid_arrays();
  }

  void read() override;

  size_t count_cell() override {
    amrex::Vector<float> vars;
    return loop_cell(false, vars);
  }

  size_t count_zone() override {
    amrex::Vector<size_t> zones;
    return loop_zone(false, zones);
  }

  void get_cell(amrex::Vector<float>& vars) override { loop_cell(true, vars); }

  void get_loc(amrex::Vector<float>& vars) override {
    loop_cell(true, vars, true);
  }

  void get_zones(amrex::Vector<size_t>& zones) override {
    loop_zone(true, zones);
  }

  size_t loop_cell(bool doStore, amrex::Vector<float>& vars,
                   bool doStoreLoc = false);
  size_t loop_zone(bool doStore, amrex::Vector<size_t>& zones);
};

#endif