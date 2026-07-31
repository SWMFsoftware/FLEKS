#include "InitialCondition.h"

#include "BeamIC.h"
#include "HybridWaveIC.h"
#include "IonAcousticWaveIC.h"
#include "LightWaveIC.h"
#include "TopHatIC.h"

void register_all_initial_conditions() {
  // Phase 1
  ICRegistry::instance().register_ic("lightwave", []() {
    return std::unique_ptr<InitialCondition>(new LightWaveIC());
  });

  // Phase 2
  ICRegistry::instance().register_ic("beam", []() {
    return std::unique_ptr<InitialCondition>(new BeamIC());
  });
  ICRegistry::instance().register_ic("ionacousticwave", []() {
    return std::unique_ptr<InitialCondition>(new IonAcousticWaveIC());
  });
  ICRegistry::instance().register_ic("hybridwave", []() {
    return std::unique_ptr<InitialCondition>(new HybridWaveIC(true, 0.02));
  });
  ICRegistry::instance().register_ic("convectionwave", []() {
    return std::unique_ptr<InitialCondition>(new HybridWaveIC(false, 0.20));
  });
  ICRegistry::instance().register_ic("tophat", []() {
    return std::unique_ptr<InitialCondition>(new TopHatIC());
  });
}
