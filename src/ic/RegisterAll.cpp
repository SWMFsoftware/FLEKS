#include "InitialCondition.h"

#include "BeamIC.h"
#include "TopHatIC.h"
#include "WaveIC.h"

void register_all_initial_conditions() {
  // Phase 5: the four wave tests are now ONE parameterized WaveIC, registered
  // under each legacy name with its preset profile (so existing PARAM.in files
  // keep working unchanged), plus a generic "waveic" for brand-new wave tests
  // configured entirely through a #WAVEIC block (zero C++ to add a wave test).
  ICRegistry::instance().register_ic("lightwave", []() {
    return std::unique_ptr<InitialCondition>(new WaveIC(WaveIC::LightWave));
  });
  ICRegistry::instance().register_ic("hybridwave", []() {
    return std::unique_ptr<InitialCondition>(new WaveIC(WaveIC::HybridWave));
  });
  ICRegistry::instance().register_ic("convectionwave", []() {
    return std::unique_ptr<InitialCondition>(new WaveIC(WaveIC::ConvectionWave));
  });
  ICRegistry::instance().register_ic("ionacousticwave", []() {
    return std::unique_ptr<InitialCondition>(new WaveIC(WaveIC::IonAcousticWave));
  });
  ICRegistry::instance().register_ic("waveic", []() {
    return std::unique_ptr<InitialCondition>(new WaveIC(WaveIC::Generic));
  });

  // Non-wave tests keep their dedicated plugins.
  ICRegistry::instance().register_ic("beam", []() {
    return std::unique_ptr<InitialCondition>(new BeamIC());
  });
  ICRegistry::instance().register_ic("tophat", []() {
    return std::unique_ptr<InitialCondition>(new TopHatIC());
  });
}
