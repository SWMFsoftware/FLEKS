# Reflecting-Wall Boundary Test (Pure Particle Push)

## Purpose

Verify the new **`reflect`** particle boundary condition in isolation.  The
electromagnetic field solver is **turned off** (`solveEM = F`, `useHybridPIC =
F`), so the fields stay frozen at the uniform initial state and the run reduces
to the particle pusher + boundary reflection.

> **Why no field solver?** A full-PIC run with a flowing plasma crossing a
> conducting wall builds an inconsistent motional (`u × B`) electric field: the
> bulk flow through a tangential guide field generates a tangential E that a
> conducting wall (requiring `E_t = 0`) cannot reconcile in the bulk.  By
> disabling the field solve, we isolate and robustly validate the reflecting
> particle boundary, which is the novel physics.

## Physics Setup

- **Domain**: a closed box in x (`8×1×1` cells, `Lx = 6.4 d_i`), fake-2D
  periodic in y/z.  Both x-faces are **reflecting**, so there is no particle
  escape path in x.
- **Plasma**: a single kinetic ion species (`q = m = 1`) streaming toward the
  `+x` wall (`ux > 0`) with a small oblique `uy != 0`.
- **Guide field `Bz`** (frozen): the Boris pusher conserves kinetic energy
  (magnetic fields do no work), so Epart conservation cleanly tests reflection.
- **Boundaries**:
  - `-x` / `+x`: **reflect** (specular: normal velocity reversed, position
    mirrored `x → 2·x_w − x`).
  - `y` / `z`: periodic.

## Test Cases Covered

| Case | How it is exercised |
|------|---------------------|
| **Normal incidence** | Ions with `ux > 0`, `uy = 0` hit the `+x` wall head-on; `vx → -vx`. |
| **Oblique incidence** | The small `uy != 0` gives a glancing component; only `vx` reverses, `vy`/`vz` are preserved (specular). |
| **Edge case — fast crossing** | A large `ux` (with `dt = 0.02`) can carry a particle past the wall in one step; the reflection must still bring it back inside. |
| **Edge case — wall cell** | A particle landing exactly in the wall ghost cell is mirrored correctly by `2·face − x`. |

## Validation

The plot uses the **`var`** form (`z=0 var ascii pic`) with an explicit
`plotVar` list (`rhoS0 Bx By Bz Ex Ey Ez uxS0 uyS0 uzS0`) so only species-0
variables are output — the `{fluid}` macro's hard-coded empty `S1` columns are
avoided.

`validate.py` checks:

1. **Energy log**: `Etot` is finite; `Epart` does **not** decay toward zero
   (specular reflection conserves kinetic energy, so a large Epart loss would
   mean particles are leaking out instead of being reflected).
2. **Plot output**: the deposited ion density `rhoS0` stays positive and finite
   across the domain (particles remain confined, no vacuum blow-out at the
   walls).  In the no-solver mode the node moment deposit is coupled to the
   E-solver and may be zero, in which case this check is skipped.

## Hybrid-PIC variant

A second variant is auto-discovered (`BC_REFLECTING (HYBRID)`).  It is the
hybrid-PIC counterpart: the same single ion species streams toward the `+x`
wall, but the field advance uses the Ohm's-law + Faraday hybrid solver
(`useHybridPIC = T`).  The **`reflect`** particle BC is exercised in the hybrid
cell-centred mover via the *same* shared hook
(`Particles::reflect_or_delete_particle()`), so `Epart` must again be retained
(specular reflection).  The Ohm's law is kept minimal (ideal frozen-in: no Hall,
no resistivity, no electron-pressure-gradient) so the motional field is small
and reflection is isolated from hybrid field-evolution artifacts.

The validator is solver-agnostic: the energy-log `Epart` conservation check and
the `rhoS0` confinement check apply unchanged to both variants.

## Running

Both variants run together under the test name:

```bash
python3 tests/validate_tests.py --test=bc_reflecting          # both variants
```

To run a **single variant**, append its token to the test name:

```bash
python3 tests/validate_tests.py --test=bc_reflecting.full      # full PIC only (PARAM.in)
python3 tests/validate_tests.py --test=bc_reflecting.hybrid    # hybrid only (PARAM.in.hybrid)
```
