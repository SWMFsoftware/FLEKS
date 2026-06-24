#ifndef _SOURCEINTERFACE_H_
#define _SOURCEINTERFACE_H_

#include "FluidInterface.h"

// An abstract class for source implementations.
class SourceInterface : public FluidInterface {
protected:
  std::string info = "SourceInterface class";
  bool useFluidSource = false;

public:
  SourceInterface(const FluidInterface& other, int id, std::string tag,
                  FluidType typeIn = SourceFluid)
      : FluidInterface(other, id, tag, typeIn) {
    initFromSWMF = false;
  }

  virtual ~SourceInterface() = default;

  virtual std::string get_info() const { return info; }

  virtual void sum_to_single_source() {
    amrex::Print()
        << "Warning: SourceInterface::sum_to_single_source is called but not "
           "implemented."
        << std::endl;
  };

  virtual void set_source(const FluidInterface& other) {
    if (!useFluidSource)
      return;

    amrex::Print() << "Warning: SourceInterface::set_source is called but not "
                      "implemented."
                   << std::endl;
  };

  /// Get total neutral exosphere density at radial distance r (SI units).
  /// Default returns 0.0 — override in UserSource for exosphere profiles.
  virtual double get_exosphere_density(double r) const { return 0.0; }

  /// Get single-component neutral exosphere density at radial distance r.
  /// Default returns 0.0 — override in UserSource for exosphere profiles.
  virtual double get_exosphere_component_density(double r, int iC) const {
    return 0.0;
  }
};

#endif
