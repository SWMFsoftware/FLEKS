#ifndef _DOMAINPARAMETERS_H_
#define _DOMAINPARAMETERS_H_

#include <string>

#include "FleksDistributionMap.h"

enum class PicSolverMode {
  Implicit, // Semi-implicit full PIC (default)
  Explicit, // Explicit full PIC
  Hybrid,   // Hybrid PIC (fluid electrons + kinetic ions)
  None      // No field solve (particle tracer / test particles)
};

// Configuration for a FLEKS Domain, owned by Domain and shared
// with its child components.
struct DomainParameters {
  // Restart control.
  bool doRestart = false;
  bool doRestartPT = false;
  bool doRestartFIOnly = false;
  bool doRegridOnRestart = false;

  // Coupling / initialization mode.
  bool isStandalone = false;
  bool initFromSWMF = true;
  bool receiveICOnly = false;
  bool useHybridPIC = false;
  PicSolverMode fieldSolverMode = PicSolverMode::Implicit;

  void set_field_solver_mode(const std::string &modeStr) {
    std::string s = modeStr;
    for (char &c : s)
      c = std::tolower(c);
    if (s == "implicit" || s == "semi-implicit") {
      fieldSolverMode = PicSolverMode::Implicit;
      useHybridPIC = false;
    } else if (s == "explicit" || s == "explicit-pic") {
      fieldSolverMode = PicSolverMode::Explicit;
      useHybridPIC = false;
    } else if (s == "hybrid" || s == "hybrid-pic") {
      fieldSolverMode = PicSolverMode::Hybrid;
      useHybridPIC = true;
    } else if (s == "none" || s == "off" || s == "false") {
      fieldSolverMode = PicSolverMode::None;
      useHybridPIC = false;
    } else {
      amrex::Abort("Invalid #FIELDSOLVER mode '" + modeStr +
                   "'. Expected 'implicit', 'explicit', 'hybrid', or 'none'.");
    }
  }

  // Component toggles; gate child construction in Domain::init().
  bool usePT = false;
  bool useSource = false;

  // Number of files per AMReX output.
  int nFileField = 64;
  int nFileParticle = 256;

  // Load-balancing strategy.
  BalanceStrategy balanceStrategy = BalanceStrategy::Cell;
  int cellWeight = 10;

  // Reporting options.
  int dnReportLB = 100;
  bool doCompact = true;

  // Return the first invalid domain-level combination, or an empty string.
  // Parsing and application remain separate from semantic validation.
  std::string validation_error() const {
    if (nFileField < 1 || nFileParticle < 1)
      return "nFileField and nFileParticle must be positive.";
    if (cellWeight < 1)
      return "cellWeight must be positive.";
    if (doRestart && receiveICOnly)
      return "#RESTART and #RECEIVEICONLY cannot both be enabled.";
    return {};
  }
};

#endif // _DOMAINPARAMETERS_H_
