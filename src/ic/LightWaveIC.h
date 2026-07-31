#ifndef _LIGHTWAVE_IC_H_
#define _LIGHTWAVE_IC_H_

#include "InitialCondition.h"

// Left/right-hand circularly polarized transverse EM wave seeded on top of the
// uniform (zero) field, with a fixed wavelength of 48 cells.
class LightWaveIC : public InitialCondition {
public:
  void set_fields(PicICFields& fields) const override;
  std::string name() const override { return "LightWave"; }
};

#endif
