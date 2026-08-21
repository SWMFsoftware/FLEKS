# tests/shock

1D oblique magnetized shock for the **HYBRID** (Hall) solver, exercising the
`fixed` (reflect-field) left wall together with a clean `inflow` right boundary.

## Setup

* Domain: 1D in x (`nCellX = 64`, `maxBlockSizeX = 16`), y/z periodic and
  single-cell.  `xMin = 0`, `xMax = 64` in code units.
* Flow: bulk velocity `ux = -713 km/s` (toward `-x`), density
  `n = 5 amu/cc`, ion temperature `T = 6.28e5 K`.
* Upstream magnetic field inclined **45 deg** to the shock normal (the x axis):
  `Bx = Bz = B0/sqrt(2) = 7.382e-9 T`, `B0 = 1.044e-8 T`.
* Hall term on (`useHallTerm T`), averaged B (`useAvgFieldB T`), electron
  temperature `Te = 141 eV`, `electronDensity0 = 5.0`.
* Boundaries:
  * **-x (left) wall = `fixed`** (reflect-field BC, the same EM-field machinery
    as `inflow`): the ghost B and E are pinned to the upstream state, so
    particles reflect off a wall that holds the upstream field while the
    tangential motional E `E = -u x B` is preserved.
  * **+x (right) = `inflow`**: the upstream state is supplied via the `#INFLOW`
    command (SI units) and re-seeded each step, replenishing the flow so a steady
    shock can form at the left wall.
  * y/z = `periodic`.

The upstream motional electric field is `E = -u x B = (0, -ux*Bz, 0)`:
`ux = -713e3 m/s`, `Bz = 7.382e-9 T` gives `Ey = +5.263e-3 V/m`, which is set
explicitly in `#UNIFORMSTATE` and matches the `#INFLOW` upstream state.

## What it validates

A shock forms at the left wall and is continuously fed from the right (inflow)
face.  With the Hall term on, this is a non-trivial boundary-condition
interaction: the left wall must hold the upstream field while the right face
keeps injecting the upstream plasma.

**Regression guarded:** the original `CONDUCTING`/`PEC` left wall forced
`E_t = 0`, which is incompatible with a flowing magnetized plasma (the motional
field `E = -u x B` demands a non-zero tangential E).  With that wall the
**field** energy `Eb` diverged *first* (e.g. 18.5 -> 6.1e4 within ~960 steps, and
~5e18 by the end) while the particle energy `Epart` stayed matched to its initial
value until the very end — i.e. the blow-up was field-driven, not
particle-driven.  The fix pins the left-wall EM field to the upstream state via
the `fixed` BC (same EM machinery as `inflow`) instead of forcing `E_t = 0`.

## Validation

The run is validated in two ways, both implemented in `validate.py`
(see that file's module docstring for the regression context).  The run must
also finish with exit code 0 (enforced by the test harness), so a clean
completion through `TimeMax` is part of the pass criteria.

### `validate_log` — PIC energy log

* `Etot`, `Ee`, `Eb`, `Epart` finite every frame (no NaN/Inf).
* `Eb` and `Epart` growth ratios stay within `[0.2, 1e3]` — generous for a
  legitimate shock, but violated by the Inf/NaN divergence the fix removed.
* `Ee` stays below `1e2 * Eb` (no spurious electric-field build-up).

### `validate_plot` — final `#SAVEPLOT` snapshot

With `ux < 0` the inflow is at `+x`, so the rightmost 25 % is the **upstream**
region and the leftmost 50 % is the **downstream** (wall) region.  Checks are
self-consistent (the upstream target is read from the snapshot's own steady
inflow region — no brittle SI→code conversions) and verify the two things this
test exists to confirm:

* **Inflow upstream stable (#2):** in the upstream region `Bx`, `rhoS0`,
  `uxS0`, and `Ey` (and `|B|`) are each steady and uniform — relative std < 0.15
  (the boundary field is not being eroded and the injected state is not
  drifting).
* **Shock formed, left = downstream, right = upstream (#1):** the downstream
  mean `rhoS0` (> 1.2×) and mean `|B|` (> 1.1×) exceed the upstream means, and
  the peak density reaches at least `1.5×` the upstream mean (a clear jump, not
  uniform flow). Any NaN/Inf in the wall region also fails.

## Running

```
python3 tests/validate_tests.py --test=shock
```

The final snapshot is written to `run_test/PC/plots/y=0_fluid_region1_0_t*
_n*.out`; if `#SAVEPLOT` is disabled the plot check is skipped gracefully and
validation falls back to the energy-log checks above.
