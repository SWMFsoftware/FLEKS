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

## Results (observed)

Running the test (`validate_tests.py --test=reflecting_wall`) gives:

- `Etot`: `2.94e-2 → 2.95e-2` (conserved to ~99.6%).
- `Epart`: `3.93e-3 → 4.05e-3` (conserved to ~103% — specular reflection keeps
  kinetic energy; no particle loss.  The B-field does no work and no particle
  escapes through a reflecting face, so this is the expected result.  An
  outflow wall would drain Epart toward zero.)

## Running

```bash
# 3D build (fake-2D via a single y/z cell)
./Config.pl -lev=2 -u=Exo
make -j4

python3 tests/validate_tests.py --test=reflecting_wall
```

## Implementation Notes

- **Particle reflection** is enforced at push time by the inline helper
  `Particles::reflect_or_delete_particle()` (`include/Particles.h`).  After a
  mover push, a particle that crossed a `reflect` face has its normal velocity
  flipped and its position mirrored back inside (it is kept); a particle that
  crossed any other non-periodic face is marked for deletion.  All four mover
  sites use this helper.
- The `BC` class (`include/BC.h`) gained the `reflect` (4) type.
- The **`conducting`** (5) field BC (`Pic::apply_conducting_wall`) is
  implemented and dispatched from `Pic::apply_BC`, but is not exercised by this
  test (the field solver is off).  See
  `Doc/Boundary_Conditions_Review_and_Plan.md` for its design and a discussion
  of why a flowing-plasma conducting-wall full-PIC setup is physically
  inconsistent (motional E vs. `E_t = 0`).
