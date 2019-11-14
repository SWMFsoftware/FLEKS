#ifndef _FLUIDINTERFACE_H_
#define _FLUIDINTERFACE_H_

#include "FluidPicInterface.h"

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

class FluidInterface : public FluidPicInterface {

private:
  // ------Grid info----------
  amrex::DistributionMapping dm;
  amrex::Geometry geom;
  amrex::BoxArray centerBA;
  amrex::BoxArray nodeBA;
  //------------------------

  amrex::MultiFab nodeFluid;

  int nGst;

public:
  void receive_info_from_gm(int *paramInt, double *gridDim,
                            double *paramDouble);

  void make_grid(const amrex::DistributionMapping &dmIn,
                 const amrex::Geometry geomIn, const amrex::BoxArray centerBAIn,
                 const amrex::BoxArray nodeBAIn, const int nGstIn);

  void set_node_value(double *data, int *index);
};

#endif