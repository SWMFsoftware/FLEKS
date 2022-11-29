#ifndef _FLUIDSOURCE_H_
#define _FLUIDSOURCE_H_

#include "FluidInterface.h"

extern "C" {
void get_source_wrapper(double xyzSI[3], double sourceSI[6]);
}

class FluidSource : public FluidInterface {
public:
  FluidSource(const FluidInterface& other, int id, std::string tag,
              FluidType typeIn = SourceFluid)
      : FluidInterface(other, id, tag, typeIn) {
    initFromSWMF = false;
  }

  void set_node_fluid(const FluidInterface& other) override {
    std::string nameFunc = "FS:set_node_fluid";
    amrex::Print() << nameFunc << " is called.";
    
    FluidInterface::set_node_fluid(other);

  }
};

#endif