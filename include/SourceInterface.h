#ifndef _SOURCEINTERFACE_H_
#define _SOURCEINTERFACE_H_

#include "FluidInterface.h"

// An abstract class for source implementations.
// Ionization parameters are stored here (not in FluidInterface) so that
// FluidInterface is only touched for neutral-profile data that the MHD
// coupling layer needs.
class SourceInterface : public FluidInterface {
protected:
  std::string info = "SourceInterface class";
  bool useFluidSource = false;

  // ---- Ionization process flags (set by #PHOTOIONIZATION, etc.) ----
  bool usePhotoIonization = false;
  bool useElectronImpact = false;
  bool useChargeExchange = false;

  // ---- Photoionization (#PHOTOIONIZATION command) ----
  int nPhotoComponent = 0;
  amrex::Vector<double> photoNu0; // ionization rate at planet surface [s^-1]

  // ---- Electron impact ionization (#ELECTRONIMPACT command) ----
  int nImpactComponent = 0;
  amrex::Vector<double> impactEIon; // ionization energy [eV]
  amrex::Vector<double> impactA;    // Voronov A coefficient [cm^3/s]
  amrex::Vector<double> impactK;    // Voronov K coefficient
  amrex::Vector<double> impactX;    // Voronov X coefficient

  // ---- Charge exchange (#CHARGEEXCHANGE command) ----
  int nCXComponent = 0;
  amrex::Vector<double> cxSigma; // cross section [cm^2]

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

  /// Read ionization-related parameter commands.
  /// Override in UserSource to handle #PHOTOIONIZATION, #ELECTRONIMPACT,
  /// #CHARGEEXCHANGE.
  virtual void read_param(const std::string& command, ReadParam& param) {
    amrex::ignore_unused(command, param);
  }

  /// Validate consistency across ionization commands after all parameters
  /// have been read.  Called after read_param() for all commands.
  virtual void post_process_param() {}

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
