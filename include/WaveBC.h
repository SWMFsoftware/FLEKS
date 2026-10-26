#ifndef FLEKS_WAVEBC_H
#define FLEKS_WAVEBC_H

#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "GpuBackend.h"
#include "ReadParam.h"

class FluidInterface;

// One spectral component of a wave injected at a boundary face.
struct WaveComponent {
  enum {
    kMono = 0,    // A sin(k.x - w t + phi)
    kPacket = 1,  // mono * Gaussian envelope
    kCustom = 2,  // std::function hook
    kSpectral = 3 // superposition of sub-components
  };

  int iField = 0;    // 0 = B, 1 = E, 2 = velocity, 3 = density weight
  int iSpecies = -1; // species for velocity/density; -1 = EM field

  amrex::Real amplitude = 0.0;
  amrex::Real frequency = 0.0;  // omega
  amrex::Real waveLength = 0.0; // 0 => pump (uniform along dir)
  amrex::Real phase = 0.0;
  amrex::Real dir[3] = { 1, 0, 0 };

  amrex::Real tStart = 0.0;
  amrex::Real tEnd = -1.0;    // -1 => forever
  amrex::Real rampTime = 0.0; // cos^2 ramp-in
  amrex::Real tCenter = 0.0;  // packet center (kPacket)
  amrex::Real tWidth = 1.0;   // packet width  (kPacket)

  int profile = kMono;

  std::vector<WaveComponent> comps; // sub-components (kSpectral)
  int seed = 1;                     // random-phase seed (kSpectral)
  amrex::Real spectralIndex = 0.0;  // power-law index (unused in v1)

  std::function<amrex::Real(const WaveComponent&, amrex::Real,
                            const amrex::Real*)>
      custom = nullptr; // kCustom

  amrex::Real pol[3] = { 0, 1, 0 };   // polarisation vector
  amrex::Real k_vec[3] = { 0, 0, 0 }; // precomputed (2*pi/waveLength) * dir
};

// Per-face emitter stack.
struct WaveFace {
  int direction = 0; // 0,1,2 = x,y,z
  int side = 0;      // 0 = lo, 1 = hi

  std::vector<WaveComponent> comps;
};

class WaveProfile {
public:
  virtual ~WaveProfile() = default;
  virtual amrex::Real value(const WaveComponent& c, amrex::Real t,
                            const amrex::Real* pos) const = 0;
};

class MonoWave : public WaveProfile {
public:
  amrex::Real value(const WaveComponent& c, amrex::Real t,
                    const amrex::Real* pos) const override;
};

class WavePacket : public MonoWave {
public:
  amrex::Real value(const WaveComponent& c, amrex::Real t,
                    const amrex::Real* pos) const override;
};

class SpectralWave : public WaveProfile {
public:
  amrex::Real value(const WaveComponent& c, amrex::Real t,
                    const amrex::Real* pos) const override;
};

class WaveRegistry {
public:
  static std::unique_ptr<WaveProfile> create(const WaveComponent& c);
};

// Runtime-editable wave-boundary manager held by Pic.
class WaveBoundaryManager {
public:
  bool active = false;
  std::vector<WaveFace> faces; // indexed by (direction, side)

  void read_param(ReadParam& param, const FluidInterface* fi);
  amrex::Real value(const WaveComponent& c, amrex::Real t,
                    const amrex::Real* pos) const;

  WaveFace& face(int direction, int side);
  void clear();
  void add_component(int direction, int side, const WaveComponent& c);

  amrex::Real envelope(const WaveComponent& c, amrex::Real t) const;

  bool siInput = true;
  amrex::Real maxAmplitude = 0.0; // 0 => amplitude guard disabled
};

// The velocity-kick components (iField == 2) of the wave boundary at one
// time, in a form a kernel can capture by value. The time-only factors (the
// envelope and the packet shape) are evaluated when the view is made; kick()
// evaluates the rest with the same arithmetic as WaveBoundaryManager::value.
struct WaveKickComp {
  int profile = WaveComponent::kMono;
  amrex::Real amplitude = 0.0;
  amrex::Real k_vec[3] = { 0, 0, 0 };
  amrex::Real frequency = 0.0;
  amrex::Real phase = 0.0;
  amrex::Real packet = 1.0; // packet shape exp(-tau^2) (kPacket)
  amrex::Real env = 0.0;    // envelope
  amrex::Real pol[3] = { 0, 0, 0 };
  int sub0 = 0, nSub = 0;              // sub-components (kSpectral)
  const WaveComponent* orig = nullptr; // for a custom (host) profile
};

struct WaveKickView {
  KernelArray<WaveKickComp> comps;
  KernelArray<WaveKickComp> subs; // mono sub-components of kSpectral comps

  AMREX_GPU_HOST_DEVICE static amrex::Real mono(const WaveKickComp& c,
                                                amrex::Real t,
                                                const amrex::Real* pos) {
    if (c.amplitude == 0.0)
      return 0.0;
    const amrex::Real kdotx =
        c.k_vec[0] * pos[0] + c.k_vec[1] * pos[1] + c.k_vec[2] * pos[2];
    const amrex::Real arg = kdotx - c.frequency * t + c.phase;
    return c.amplitude * std::sin(arg);
  }

  // Add the velocity kick at position pos and time t (the time of the view).
  AMREX_GPU_HOST_DEVICE void kick(const amrex::Real* pos, amrex::Real t,
                                  amrex::Real& dvx, amrex::Real& dvy,
                                  amrex::Real& dvz) const {
    dvx = dvy = dvz = 0.0;
    for (int i = 0; i < comps.size(); ++i) {
      const WaveKickComp& c = comps[i];
      if (c.env == 0.0)
        continue;
      amrex::Real val = 0.0;
      if (c.profile == WaveComponent::kCustom) {
        // A custom profile is a host function (rejected in a GPU run).
        AMREX_IF_ON_HOST((if (c.orig && c.orig->custom) {
          val = c.orig->custom(*c.orig, t, pos) * c.env;
        }))
      } else if (c.profile == WaveComponent::kSpectral) {
        amrex::Real sum = 0.0;
        for (int j = c.sub0; j < c.sub0 + c.nSub; ++j)
          sum += mono(subs[j], t, pos);
        val = sum * c.env;
      } else if (c.profile == WaveComponent::kPacket) {
        const amrex::Real carrier = mono(c, t, pos);
        val = (carrier == 0.0 ? 0.0 : carrier * c.packet) * c.env;
      } else {
        val = mono(c, t, pos) * c.env;
      }
      if (val == 0.0)
        continue;
      dvx += val * c.pol[0];
      dvy += val * c.pol[1];
      dvz += val * c.pol[2];
    }
  }
};

// The velocity kick of `wave` at time t (empty if the wave BC is inactive).
WaveKickView make_wave_kick_view(const WaveBoundaryManager& wave,
                                 amrex::Real t);

#endif // FLEKS_WAVEBC_H
