#include "BeamIC.h"
#include "FadeevIC.h"
#include "InitialCondition.h"
#include "TopHatIC.h"
#include "WaveIC.h"

void register_all_initial_conditions() {
  // The four wave tests share one parameterized WaveIC (names + generic config)
  ICRegistry::instance().register_ic("lightwave", []() {
    return std::make_unique<WaveIC>(WaveIC::LightWave);
  });
  ICRegistry::instance().register_ic("hybridwave", []() {
    return std::make_unique<WaveIC>(WaveIC::HybridWave);
  });
  ICRegistry::instance().register_ic("convectionwave", []() {
    return std::make_unique<WaveIC>(WaveIC::ConvectionWave);
  });
  ICRegistry::instance().register_ic("ionacousticwave", []() {
    return std::make_unique<WaveIC>(WaveIC::IonAcousticWave);
  });
  ICRegistry::instance().register_ic(
      "waveic", []() { return std::make_unique<WaveIC>(WaveIC::Generic); });

  ICRegistry::instance().register_ic(
      "beam", []() { return std::make_unique<BeamIC>(); });
  ICRegistry::instance().register_ic(
      "tophat", []() { return std::make_unique<TopHatIC>(); });
  ICRegistry::instance().register_ic(
      "fadeev", []() { return std::make_unique<FadeevIC>(); });
}
