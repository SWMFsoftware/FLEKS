#ifndef _SOURCEINTERFACE_H_
#define _SOURCEINTERFACE_H_

#include <string>
#include <vector>

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
  amrex::Vector<amrex::Real> photoNu0; // ionization rate at planet surface
                                       // [s^-1]

  // ---- Shadow cylinder (#SHADOWCYLINDER command) ----
  bool useShadowCylinder = false;
  amrex::Real solarDir[3] = { 0.0, 0.0, 0.0 }; // unit vector toward the Sun
  amrex::Real shadowCylinderRadius = 0.0;      // shadow cylinder radius [m]
  amrex::Real shadowCylinderHalfHeight = 0.0;  // half-height anti-solar [m]

  // ---- Electron impact ionization (#ELECTRONIMPACT command) ----
  amrex::Vector<amrex::Real> impactEIon; // ionization energy [eV]
  amrex::Vector<amrex::Real> impactA;    // Voronov A coefficient [cm^3/s]
  amrex::Vector<amrex::Real> impactK;    // Voronov K coefficient
  amrex::Vector<amrex::Real> impactX;    // Voronov X coefficient

  // ---- Charge exchange (#CHARGEEXCHANGE command) ----
  int nCXIonSpecies = 0; // number of ion species that exchange charge
  // Cross-section matrix [cm^2], flattened as [iC * nCXIonSpecies + iIon].
  amrex::Vector<amrex::Real> cxSigma;

  // ---- Recombination (#RECOMBINATION command) ----
  bool useRecombination = false;
  amrex::Vector<int> recombIonIndex;        // ion species index (iSp)
  amrex::Vector<amrex::Real> recombRate0;   // base rate coefficient k0 [cm^3/s]
  amrex::Vector<amrex::Real> recombTempExp; // temperature exponent alpha
  amrex::Vector<amrex::Real> recombRefTemp; // reference temperature T_ref [K]

  // ---- General chemistry (#CHEMISTRY command) ----
  bool useChemistry = false;
  struct ChemistryReaction {
    int reactantIon;      // 0 = none, 1+ = ion species index
    int productIon;       // 0 = none, 1+ = ion species index
    int neutralComp;      // -1 = none, 0+ = exosphere component
    int rateType;         // 0 = thermal k(T), 1 = photoionization (1/r^2)
    amrex::Real rateCoef; // k0 [cm^3/s] for thermal, nu0 [s^-1] for photo
    amrex::Real tempExp;  // alpha: k = k0 * (Tref/Te)^alpha
    amrex::Real refTemp;  // T_ref [K]
  };
  amrex::Vector<ChemistryReaction> chemReactions;

  // ---- Loss term storage ----
  // Loss frequency [code units of 1/time] of every species at the nodes.
  amrex::Vector<amrex::MultiFab> nodeLossFluid;

  // Time step [s] the next set_source() is applied over. The stiff
  // reactions use it to keep the products consistent with the loss.
  amrex::Real dtSI = 0.0;

public:
  SourceInterface(Grid& gridIn, const FluidInterface& other, int id,
                  std::string tag, FluidType typeIn)
      : FluidInterface(gridIn, other, id, tag, typeIn) {
    initFromSWMF = false;
  }

  virtual ~SourceInterface() = default;

  virtual std::string get_info() const { return info; }

  // Register source-specific PARAM.in commands with the domain dispatcher.
  // Implementations own the command list and may extend it without editing
  // Domain.cpp.
  virtual void register_parameter_commands(
      std::vector<std::string>& commands) const {
    (void)commands;
  }

  virtual void sum_to_single_source() {
    amrex::Print()
        << "Warning: SourceInterface::sum_to_single_source is called but not "
           "implemented."
        << std::endl;
  };

  void set_dt_si(amrex::Real dt) { dtSI = dt; }
  amrex::Real get_dt_si() const { return dtSI; }

  /// Radius [m] of the inner body (#BODY / #INNERBOUNDARY): no neutral reacts
  /// inside it.
  virtual void set_inner_radius_si(amrex::Real r) { amrex::ignore_unused(r); }

  /// Radius [m] of the photochemical equilibrium boundary (#INNERBOUNDARY
  /// photochemical).
  virtual void set_photochem_radius_si(amrex::Real r) {
    amrex::ignore_unused(r);
  }

  virtual void set_source(const FluidInterface& other) {
    if (!useFluidSource)
      return;

    amrex::Print() << "Warning: SourceInterface::set_source is called but not "
                      "implemented."
                   << std::endl;
  };

  /// Read ionization-related parameter commands.
  virtual void read_param(const std::string& command, ReadParam& param) {
    amrex::ignore_unused(command, param);
  }

  /// Validate consistency across ionization commands after all parameters
  /// have been read.  Called after read_param() for all commands.
  virtual void post_process_param() {}

  /// Get total neutral exosphere density at radial distance r (SI units).
  virtual amrex::Real get_exosphere_density(amrex::Real r) const { return 0.0; }

  /// Get single-component neutral exosphere density at radial distance r.
  virtual amrex::Real get_exosphere_component_density(amrex::Real r,
                                                      int iC) const {
    return 0.0;
  }

  // ---- Loss term management ----
  // Requires species ordering: species 0 = electron, 1..nS-1 = ions.

  void post_regrid() override {
    FluidInterface::post_regrid(); // distributes nodeFluid
    distribute_loss_arrays();
  }

  void distribute_loss_arrays() {
    if (!use_loss_source())
      return;
    if (nodeLossFluid.empty())
      nodeLossFluid.resize(n_lev_max());
    if (nS == 0)
      return;
    const bool doCopy = true;
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      distribute_FabArray(nodeLossFluid[iLev], grid.node_box_array(iLev),
                          DistributionMap(iLev), nS, grid.get_n_ghost(),
                          doCopy);
    }
  }

  void set_node_loss_fluid_to_zero() {
    if (!use_loss_source())
      return;
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      if (!nodeLossFluid[iLev].empty())
        nodeLossFluid[iLev].setVal(0.0);
    }
  }

  void fill_loss_boundary() {
    if (!use_loss_source())
      return;
    // Use FillBoundary (copy) instead of SumBoundary (sum) for loss rates.
    // Each node's loss rate is computed independently and should NOT be
    // summed across shared/periodic nodes.  FillBoundary fills ghost cells
    // with the correct values from the interior, handling periodic BCs.
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      if (!nodeLossFluid[iLev].empty())
        nodeLossFluid[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }

  void sum_loss_boundary() { fill_loss_boundary(); }

  /// Loss frequency [code units of 1/time] of species iSp at the position
  /// xyz, interpolated from the nodes. A particle of this species keeps the
  /// fraction exp(-nu*dt) of its weight over a step dt.
  amrex::Real get_loss_frequency(const amrex::MFIter& mfi,
                                 const amrex::RealVect& xyz, const int iSp,
                                 const int iLev = 0) const {
    return get_value_at_loc(nodeLossFluid[iLev], mfi, Geom(iLev), xyz, iSp);
  }

  /// Survival factor of a particle of species iSp over timestep dt at position
  /// xyz, multilinearly interpolated from the nodal decay factors exp(-nu*dt).
  /// This conservative discretization ensures exact balance with the stiff
  /// product source (UserSource::set_source) even when reaction rates vary
  /// spatially.
  amrex::Real get_loss_decay(const amrex::MFIter& mfi,
                             const amrex::RealVect& xyz, const int iSp,
                             const amrex::Real dt, const int iLev = 0) const {
    return get_decay_at_loc(nodeLossFluid[iLev], mfi, Geom(iLev), xyz, iSp, dt);
  }

  amrex::MultiFab& get_node_loss_fluid(int iLev) { return nodeLossFluid[iLev]; }

  bool use_loss_source() const {
    return useRecombination || useChemistry || useChargeExchange;
  }

  /// Check whether nodeLossFluid is allocated and non-empty at iLev.
  bool has_loss_array(int iLev) const {
    return use_loss_source() && iLev < static_cast<int>(nodeLossFluid.size()) &&
           !nodeLossFluid[iLev].empty();
  }
};

#endif
