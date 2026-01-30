#ifndef _SOURCEINTERFACE_H_
#define _SOURCEINTERFACE_H_

#include "FluidInterface.h"

// An abstract class for source implementations.
class SourceInterface : public FluidInterface {
public:
  SourceInterface(const FluidInterface& other, int id, std::string tag,
                  FluidType typeIn = SourceFluid)
      : FluidInterface(other, id, tag, typeIn) {
    initFromSWMF = false;
  }

  virtual ~SourceInterface() = default;

  virtual void sum_to_single_source() {
    amrex::Print()
        << "Warning: SourceInterface::sum_to_single_source is called but not "
           "implemented.";
  };

  virtual void set_source(const FluidInterface& other) {
    amrex::Print() << "Warning: SourceInterface::set_source is called but not "
                      "implemented.";
  };
};

#endif