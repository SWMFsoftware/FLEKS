# Wave-Injection Boundary-Condition Test

## Purpose

Verify the new **wave-injection boundary condition** (`BC::wave` + `#WAVEBC`).
A monochromatic transverse magnetic wave is injected at the `x = xLo` boundary
face into a pure-EM (lightwave-style) domain and propagates into the interior.

## Physics Setup

- **Domain**: pure-EM vacuum (no particles, `#TESTCASE lightwave` sets
  `nPartPerCell = 0`), so the injected EM wave propagates cleanly without
  plasma dispersion.
- **Boundaries**: `#BFIELDBOXBOUNDARY` sets the `x-lo` face to `wave`, the rest
  periodic.
- **`#WAVEBC`**: one component on the `x-lo` face:
  - `iField = 0` → B field
  - `profile = 0` → monochromatic
  - `amplitude = 1e-4 T` (SI)
  - `frequency = 0.1 Hz` (SI, converted to code units)
  - `waveLength = 0` → pump / uniform along the propagation direction
  - `dir = (1,0,0)` → propagates along +x
  - `pol = (0,0,1)` → Bz polarization
  - continuous in time (`tStart=0`, `tEnd=-1`, no ramp)

## Validation

`validate.py`:
1. **Energy log**: `Etot` finite, positive, and approximately conserved
   (no blow-up or collapse).
2. **Plot output**: the injected `Bz` perturbation is non-zero in the interior
   (middle third of the domain), confirming the wave propagated in from the
   injection face rather than staying confined to the boundary ghost cells.

## Running

```bash
make EXE
python3 tests/validate_tests.py --test=wave_inject
```

## Implementation

- `include/WaveBC.h` / `src/WaveBC.cpp`: `WaveComponent`, `WaveEnvelope`,
  `WaveFace`, `WaveProfile` (mono/packet/spectral), `WaveBoundaryManager`
  (`#WAVEBC` parsing, SI→code conversion, amplitude guard).
- `include/BC.h`: `BC::wave = 6`.
- `src/Pic.cpp`: `Pic::apply_wave_field` (hard source), `Pic::wave_velocity_kick`;
  `#WAVEBC` parsed in `Pic::read_param`; wave applied at the field-fill and
  E-update sites.
- `include/Particles.h` / `src/Particles.cpp`: `waveVelocityKick` hook so
  boundary-injected particles carry the matching bulk-velocity perturbation.
- `include/FluidInterface.h`: `get_No2SiE()` accessor for SI→code E conversion.
