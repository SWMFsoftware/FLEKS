#ifndef _DOMAINPARAMETERS_H_
#define _DOMAINPARAMETERS_H_

#include "FleksDistributionMap.h"

// Shared, rarely-changing configuration for a FLEKS Domain.
//
// Owned privately by Domain (see Domain.h); only Domain mutates it, via
// Domain::read_param(). Children (pic/pt/fi/source) never own or copy it.
// Future step: inject a narrow const reference into the consumers that need a
// slice, instead of reaching through Domain.
struct DomainParameters {
  // Restart control.
  bool doRestart = false;
  bool doRestartPT = false;
  bool doRestartFIOnly = false;

  // Coupling / initialization mode.
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
};

#endif // _DOMAINPARAMETERS_H_
