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

#include "grid.h"

enum class FileType { AMREX = 0, IDL, TECPLOT, VTK, HDF5, ADIOS2 };

static const std::map<FileType, std::string> fileTypeString = {
  { FileType::AMREX, "AMReX" },     { FileType::IDL, "IDL" },
  { FileType::TECPLOT, "TECPLOT" }, { FileType::VTK, "VTK" },
  { FileType::HDF5, "HDF5" },       { FileType::ADIOS2, "ADIOS2" }
};

// Virtual base
class DataContainer {
public:
  virtual ~DataContainer() = default;

  virtual void read(){};
  virtual void write(){};

  std::string type_string() { return fileTypeString.at(dataType); }

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

  int nCell;
  int nBrick;
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
  int loop_cell(bool doCountOnly = false);
  int loop_brick(bool doCountOnly = false);
};

#endif