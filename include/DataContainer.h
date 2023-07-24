#ifndef _DATACONTAINER_H_
#define _DATACONTAINER_H_

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
  virtual void write() = 0;

  virtual size_t count_cell() = 0;

  virtual size_t count_brick() = 0;

  virtual void get_cell(amrex::Vector<float>& vars) = 0;

  virtual void get_bricks(amrex::Vector<size_t>& bricks) = 0;

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

  amrex::Vector<amrex::MultiFab> mf;
  amrex::Vector<amrex::MultiFab> iCell;

  size_t nCell;
  size_t nBrick;
  int nVar;
  int nDim;

  amrex::RealBox domain;
  amrex::Box cellBox;

  std::string plot_string;
  amrex::Real rPlanet;

  amrex::Real time;

  amrex::Vector<std::string> varNames;

  std::ofstream outFile;
};

class AMReXDataContainer : public Grid, public DataContainer {
private:
public:
  AMReXDataContainer(const std::string& in, const amrex::Geometry& gm,
                     const amrex::AmrInfo& amrInfo)
      : Grid(gm, amrInfo, 2), DataContainer() {
    dirIn = in;

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
  void write() override;

  size_t count_cell() override {
    amrex::Vector<float> vars;
    return loop_cell(false, false, vars);
  }

  size_t count_brick() override {
    amrex::Vector<size_t> bricks;
    return loop_brick(false, false, bricks);
  }

  void get_cell(amrex::Vector<float>& vars) override {
    loop_cell(false, true, vars);
  }

  void get_loc(amrex::Vector<float>& vars) override {
    loop_cell(false, true, vars, true);
  }

  void get_bricks(amrex::Vector<size_t>& bricks) override {
    loop_brick(false, true, bricks);
  }

  void write_cell() {
    amrex::Vector<float> vars;
    loop_cell(true, false, vars);
  }

  void write_brick() {
    amrex::Vector<size_t> bricks;
    loop_brick(true, false, bricks);
  }

  size_t loop_cell(bool doWrite, bool doStore, amrex::Vector<float>& vars,
                   bool doStoreLoc = false);
  size_t loop_brick(bool doWrite, bool doStore, amrex::Vector<size_t>& bricks);
};

#endif