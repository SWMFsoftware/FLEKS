#ifndef _DOMAINPARAMETERS_H_
#define _DOMAINPARAMETERS_H_

#include <string>

#include "FleksDistributionMap.h"

// Configuration for a FLEKS Domain, owned by Domain and shared
// with its child components.
struct DomainParameters {
  // Restart control.
  bool doRestart = false;
  bool doRestartPT = false;
  bool doRestartFIOnly = false;

  // Coupling / initialization mode.
  bool isStandalone = false;
  bool initFromSWMF = true;
  bool receiveICOnly = false;

  // Component toggles; gate child construction in Domain::init().
  bool usePT = false;
  bool useSource = false;

  // Number of files per AMReX output.
  int nFileField = 64;
  int nFileParticle = 256;

  // Load-balancing strategy.
  BalanceStrategy balanceStrategy = BalanceStrategy::Cell;
  int cellWeight = 10;

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
