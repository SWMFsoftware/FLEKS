# Absorbing Particle Boundary Test

Verify the **`absorb`** particle boundary condition: a macroparticle that
crosses an absorbing face is **removed** (and tallied) instead of being
reflected or refluxed back into the domain.

The electromagnetic field solver is **OFF** (`solveEM = F`, `useHybridPIC = F`),
so the fields stay frozen and the run exercises only the particle pusher and the
absorbing boundary.  A single kinetic ion species (`q = m = 1`) is initialised
with a bulk flow toward `+x` and **zero magnetic field**, so the ions move
ballistically.

Both x faces are **`absorb`**.  Using `absorb` on *both* x faces is essential:
a reflecting face would re-inject particles and hide the absorption (see the
`reflect` re-injection path in `inject_particles_for_boundary_cells`).  With no
face re-injecting, the ions drain out through the x walls and the **particle
kinetic energy `Epart` decays toward zero** — the opposite of the
`reflecting_wall` test, where specular reflection conserves `Epart`.

## Physics Setup

- **Domain**: `Lx = 6.4` (code units), fake-2D/1D (`1×1` cells in y/z,
  periodic).
- **Solver**: field OFF (`solveEM = F`, `useHybridPIC = F`).
- **Particles**: one ion species (`q = m = 1`), bulk `ux = 40 km/s` toward
  `+x`, `T = 10^4 K`.
- **Boundaries**: both x faces **`absorb`**; y/z periodic.
- **Run**: `TimeMax = 20`, `dt = 0.02` (fixed).  With `ux_code ≈ 0.8` the ions
  traverse the box many times and are fully drained through the x faces.

## Validation

`validate.py` checks:

1. **Energy log (primary)**: total energy is finite and positive; particle
   kinetic energy **`Epart` decays** (`Epart_final < 50%` of initial).  A
   reflecting wall (or a no-op absorber) would NOT decay `Epart` — it would
   conserve it (see `tests/reflecting_wall`).
2. **Plot output (secondary)**: if the deposited ion density `rhoS0` is present
   in the plot, it must be finite and non-negative (no NaN).

## Running

```bash
python3 tests/validate_tests.py --test=absorb_particles
```
