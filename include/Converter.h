#ifndef _CONVERTER_H_
#define _CONVERTER_H_

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

class Converter {
private:
  std::string dirIn;
  amrex::MultiFab mf;
  amrex::iMultiFab iCell;
  int nCell;
  int nBrick;
  amrex::RealBox domainRange;
  std::string plot_string;
  amrex::Real rPlanet;

  amrex::Real time;

  amrex::Vector<std::string> varNames;

  std::ofstream outFile;

public:
  Converter(const std::string& in) { dirIn = in; }
  ~Converter() = default;

  void read();
  void write();
  int loop_cell(bool doCountOnly = false);
  int loop_brick(bool doCountOnly = false);
};

#endif
