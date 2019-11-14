#include "FluidInterface.h"

using namespace amrex;

void FluidInterface::receive_info_from_gm(int *paramInt, double *gridDim,
                                          double *paramDouble) {
  std::stringstream *ss = nullptr;
  ReadFromGMinit(paramInt, gridDim, paramDouble, ss);
}

void FluidInterface::make_grid(const amrex::DistributionMapping &dmIn,
                               const amrex::Geometry geomIn,
                               const amrex::BoxArray centerBAIn,
                               const amrex::BoxArray nodeBAIn,
                               const int nGstIn) {
  dm = dmIn;
  geom = geomIn;
  centerBA = centerBAIn;
  nodeBA = nodeBAIn;
  nGst = nGstIn;
  nodeFluid.define(nodeBA, dm, nVarCoupling, nGst);
  nodeFluid.setVal(0);  
}



void FluidInterface::set_node_value(double *data, int *index){


}