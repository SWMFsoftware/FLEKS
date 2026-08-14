# Absorbing Field Boundary Test

Verify the **`absorb`** field boundary condition: an outgoing
electromagnetic pulse leaves the domain through an absorbing face with
**minimal reflection**.

The `tophat` testcase seeds a discontinuous `Ey = 1, Bz = 1` EM step in the
central half of the domain with **zero macroparticles** (pure EM, field solver
ON).  The step immediately splits into a right-going and a left-going pulse that
travel toward the `+x` and `-x` faces.  Both x-faces use the **absorbing**
(outgoing-characteristic) field BC, so the pulses are carried out of the domain
and the **total EM energy decays**.

## Physics Setup

- **Domain**: `Lx = 20` (code units), fake-2D (`1×1` cells in y/z, periodic).
  The short domain means a light-speed pulse (`c = 1`) reaches a wall in
  `t ≈ 5`, so over `TimeMax = 40` it transits the box many times and is fully
  drained.
- **Initial condition**: `#TESTCASE tophat` seeds `Ey = Bz = 1` in the central
  half of the x-domain, zero everywhere else, no macroparticles.
- **Solver**: pure EM (`solveEM = T`), implicit `theta = 0.5`.
- **Boundaries**:
  - `-x` / `+x`: **`absorb`** (outgoing-characteristic matched-impedance).
  - `y` / `z`: periodic.
- **Absorber tuning**: `#ABSORB charSpeed` defaults to `0` → auto (light speed
  `c = 1` in FLEKS code units).

## Validation

`validate.py` checks:

1. **Energy log (primary, decisive)**: total EM energy (`Ee + Eb`) is finite and
   > 0 initially, and **decays** over the run (`Etot_final < 25%` of initial).
   A reflecting wall or a no-op absorber would NOT decay like this.
2. **Plot output (secondary)**: the interior EM field amplitude (`|Bz|`, `|Ey|`)
   at late time is clearly reduced below the seeded amplitude (`1.0`) — a
   shrinking wavefront rather than a standing wave.

## Running

```bash
python3 tests/validate_tests.py --test=absorb_fields
```
