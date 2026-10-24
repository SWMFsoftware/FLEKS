#ifndef DEMO_GRID_H_
#define DEMO_GRID_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>

#include "DemoConfig.h"

class DemoGrid {
public:
  explicit DemoGrid(const DemoConfig& config);

  amrex::Geometry geom;
  amrex::BoxArray ba;
  amrex::BoxArray nodeBa;
  amrex::DistributionMapping dm;
};

#endif
