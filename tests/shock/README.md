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
  * **-x (left) wall = `fixed`**: the ghost EM field is a zero-gradient copy of
    the adjacent edge cell (same filler as `inflow`, see `apply_inflow_wall`), so
    the tangential motional E `E = -u x B` is preserved and particles reflect off
    a wall holding a continuous field.
  * **+x (right) = `inflow`**: the upstream plasma state is supplied via
    `#INFLOW` (SI units: density, velocity, temperature).  The ghost EM field
    is again a zero-gradient copy of the edge cell; the upstream flow is fed in
    by the flux-weighted particle injection each step, so a steady shock forms
    at the left wall.
  * y/z = `periodic`.

The upstream magnetic field `B = (7.382, 0, 7.382)e-9 T` and the upstream
motional electric field `E = -u x B = (0, -ux*Bz, 0)` for `ux = -713e3 m/s`,
`Bz = 7.382e-9 T` give `Ey = -5.263e-3 V/m`, all set in `#UNIFORMSTATE`
(`#INFLOW` does not carry a field state).

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
a field-driven, not particle-driven, blow-up.  The fix replaces the left-wall
`CONDUCTING`/PEC BC (which forces `E_t = 0`) with the `fixed` BC, whose
zero-gradient ghost copy keeps the tangential motional E `-u x B` intact.

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

## 2D variant: not stable yet (deliberately *not* in the suite)

Extending the deck to a 2D plane (256×1×64, x–z, z periodic, same physics and the
same `reflect` + `inflow` pair) was tried on 2026-08-31 to see whether the
particle/field boundary-condition separation of `feat/boundary-conditions`
stabilizes the upstream inflow. **It does not, in either solver**:

| | full PIC (t = 8 s) | hybrid (t = 30 s requested, died at 29.6 s) |
| --- | --- | --- |
| γ(By k8 at the inflow face) | 1.46 s⁻¹ (baseline 1.41) | 0.108 s⁻¹ (baseline 0.109) |
| when the upstream is lost | saturates t ≈ 4.75 s, destroyed by t = 8 s | usable to t ≈ 18 s, dead by t ≈ 29 s |
| upstream mean drift | ux −9 %, Bz 1.81 → 28.8 | ux −10 %, Bz −13 %, ρ +5 % |
| end state | run completes, fully turbulent upstream | hangs (CFL 24, max(uth) 4853); Eb → 3.3e8 |

Both re-runs used the historical decks byte-for-byte (only the field-BC command
rename), so the comparison against the pre-refactor baselines is exact: the
refactor is numerically neutral for this instability (γ differs by 1–4 %).

The 1D decks here are unaffected — the mode needs a transverse direction.

One caveat before calling 1D "fine": a 1D control run (hybrid oblique, 256 cells,
t = 30 s, `run_hyb_1d_t30s/`) stays **numerically healthy** to 30 s (CFL 0.063,
clean exit) but drifts in the same way the 2D runs do (ρ +4.6 %, ux −7.4 %,
versus ρ +5.3 %, ux −10.3 % in 2D at the same time).  So the slow upstream
degradation is a k = 0 effect present in 1D as well — probably the
quasi-parallel foreshock (the archived 1D *perpendicular* hybrid is clean to
20 s) — and only the exponentially growing transverse mode plus the terminal
blow-up are genuinely 2D.

The evidence lives in `run_shock_2d/ANALYSIS_inflow_2d.md` §10 and the
`run_shock_2d_bcsep/` run directory (both are development artifacts, outside
version control).  The next thing to try is the higher-rate / lower-weight
boundary injection (or transverse smoothing of the injected flux) in
`inject_flux_at_inflow_faces`; `npcelz = 4` is the PARAM-only stop-gap that
halves γ at 4× cost.  Until one of those lands, no 2D full-PIC shock variant is
added to `tests/`.
