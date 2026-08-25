# tests/shock

1D magnetized shock test for FLEKS, exercising a clean `inflow` boundary together
with a reflecting wall.  Two solver variants are provided:

* **`PARAM.in`** — the **FULL PIC** solver.
* **`PARAM.in.hybrid`** — the **HYBRID** solver, which holds
  the electrons as a massless fluid.

Both variants share the same physics; only the field solver (`solveEM` /
`useHybridPIC`) differs.

## Setup

* Domain: 1D in x (`nCellX = 64`, `maxBlockSizeX = 16`), y/z periodic and
  single-cell.  `xMin = 0`, `xMax = 64` in code units.
* Flow: bulk velocity `ux = -713 km/s` (toward `-x`), density `n = 5 amu/cc`.
* Upstream magnetic field inclined **45 deg** to the shock normal (x axis):
  `Bx = Bz = B0/sqrt(2) = 7.382e-9 T`, `B0 = 1.044e-8 T`.
* Boundaries:
  * **-x (left) wall = `fixed`**: ghost B and E are pinned to the upstream state,
    so particles reflect off a wall holding the upstream field while the
    tangential motional E `E = -u x B` is preserved.
  * **+x (right) = `inflow`**: the upstream state is supplied via `#INFLOW` (SI
    units) and re-seeded each step, feeding the flow so a steady shock forms at
    the left wall.
  * y/z = `periodic`.

The upstream motional electric field `E = -u x B = (0, -ux*Bz, 0)` for
`ux = -713e3 m/s`, `Bz = 7.382e-9 T` gives `Ey = -5.263e-3 V/m`, set in
`#UNIFORMSTATE` and matching the `#INFLOW` state.

## Validation

Both variants are validated identically by `validate.py` (see its module
docstring for the regression context), and the harness requires a clean exit
code 0 through `TimeMax`.  Two checks:

### `validate_log` — PIC energy log

* `Etot`, `Ee`, `Eb`, `Epart` finite every frame (no NaN/Inf).
* `Eb` and `Epart` growth ratios within `[0.2, 1e3]` (bounded; no divergence).
* `Ee` below `1e2 * Eb` (no spurious electric-field build-up).

This guards the regression the boundary fix removed: the original
`CONDUCTING`/`PEC` left wall forced `E_t = 0`, incompatible with a flowing
magnetized plasma (the motional E `E = -u x B` needs non-zero tangential E).
With that wall the **field** energy `Eb` diverged first (18.5 -> 6.1e4 within
~960 steps) while `Epart` stayed matched to its initial value until the end —
a field-driven, not particle-driven, blow-up.  The fix pins the left-wall EM
field to the upstream state via the `fixed` BC instead of forcing `E_t = 0`.

### `validate_plot` — final `#SAVEPLOT` snapshot

With `ux < 0` the inflow is at `+x`, so the rightmost 25 % is the **upstream**
region and the leftmost 50 % is the **downstream** (wall) region.  Checks are
self-consistent (the upstream target is read from the snapshot's own steady
inflow region):

* **Inflow upstream stable:** `Bx`, `rhoS0`, `uxS0`, `Ey`, `|B|` steady and
  uniform upstream (relative std < 0.15) — the boundary field is not eroded and
  the injected state does not drift.
* **Shock formed:** downstream mean `rhoS0` (> 1.2×) and `|B|` (> 1.1×) exceed
  the upstream means, and peak density reaches at least `1.5×` the upstream mean
  (a clear jump, not uniform flow).  Any NaN/Inf in the wall region also fails.

## Running

```
python3 tests/validate_tests.py --test=shock          # both (PARAM.in + hybrid)
python3 tests/validate_tests.py --test=shock.full     # full PIC only (PARAM.in)
python3 tests/validate_tests.py --test=shock.hybrid   # hybrid only (PARAM.in.hybrid)
```
