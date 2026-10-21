# tests/shock

Magnetized collisionless shock test for FLEKS: a super-magnetosonic plasma
streams into a reflecting wall, exercising a reflecting wall and an upstream
inflow boundary in the same deck.  Four variants are provided:

* **`PARAM.in.implicit`** — the **SEMI-IMPLICIT FULL PIC** solver (kinetic ions + kinetic electrons,
  reduced mass ratio `m_i/m_e = 25`, `#FIELDSOLVER implicit`).
* **`PARAM.in.hybrid`** — the **HYBRID** solver (kinetic ions, massless
  electron fluid at `Te ≈ 140 eV`, `#FIELDSOLVER hybrid`).
* **`PARAM.in.hybrid_2d`**, **`PARAM.in.hybrid_2d_hyper`** — the hybrid
  solver in 2D at twice the resolution and four times the run length, where
  grid-scale noise must be damped (see [2D variants](#2d-variants)).

All variants share the same upstream physics.

## Setup

The upstream is a moderately strong super-critical shock:

| quantity | value |
| --- | --- |
| shock-normal angle θ_Bn | **45°** (oblique) |
| Alfvén Mach number `M_A` | **≈ 7** (flow speed ≈ 7 `v_A`) |
| fast-magnetosonic Mach number `M_f` | **≈ 5** |
| ion plasma beta `β_i` | **≈ 1** (thermal ≈ magnetic pressure) |
| ion inertial length `d_i` | ≈ 1 code unit |
| domain / cell size | ≈ 63 `d_i` long, `dx ≈ 1 d_i` (ion-scale resolution) |
| run length | 5 s ≈ **0.8 ion gyroperiods** (`Ω_ci ≈ 1` in code units) |

At `M_f ≈ 5` the shock is well beyond the second critical Mach number, so ion
reflection drives a genuine super-critical shock: the downstream is compressed
by a factor of a few and the wall region carries the reflected-ion foot.

* Boundaries:
  * **-x (left) wall = `outflow` (field) + `reflect` (particles)**: the ghost
    EM field is a zero-gradient copy of the edge cell, so the wall imposes no
    constraint -- the tangential motional E `E = -u x B` stays continuous and B
    is free to compress as the shock forms.  This is the physical pairing with a
    reflecting wall: `conducting`/`symmetry` would force `E_t = 0` or `B_n = 0`
    (wrong for an oblique flowing plasma), and `fixed` would pin a compressed
    wall to the upstream state.
  * **+x (right) = `fixed` field + `inflow` particles**: the upstream plasma
    state is supplied by `#INFLOW` through the flux-weighted particle injection
    each step, while the ghost EM field is Dirichlet-pinned to the uniform state.
  * y/z = `periodic` (single-cell).

## Validation

Both variants are validated identically by `validate.py` (see its module
docstring for the regression context), and the harness requires a clean exit
code 0 through `TimeMax`.  Two checks:

### PIC energy log

* `Etot`, `Ee`, `Eb`, `Epart` finite every frame (no NaN/Inf).
* `Eb` and `Epart` growth ratios within `[0.2, 1e3]` (bounded; no divergence).
* `Ee` below `1e2 * Eb` (no spurious electric-field build-up).

This guards the regression the boundary fix removed: the original
`conducting`/PEC left wall forced `E_t = 0`, incompatible with a flowing
magnetized plasma (the motional E `E = -u x B` needs non-zero tangential E).
With that wall the **field** energy `Eb` diverged first (18.5 -> 6.1e4 within
~960 steps) while `Epart` stayed matched to its initial value until the end --
a field-driven, not particle-driven, blow-up.  The wall now uses `outflow`
(zero-gradient), which leaves the tangential motional E `-u x B` intact.

### Final snapshot

With `ux < 0` the inflow is at `+x`, so the rightmost 25 % is the **upstream**
region and the leftmost 50 % is the **downstream** (wall) region.  Checks are
self-consistent (the upstream target is read from the snapshot's own steady
inflow region):

* **Inflow upstream stable:** `Bx`, `rhoS0`, `uxS0`, `Ey`, `|B|` steady and
  uniform upstream (relative std < 0.15) -- the boundary field is not eroded and
  the injected state does not drift.
* **Shock formed:** downstream mean `rhoS0` (> 1.2×) and `|B|` (> 1.1×) exceed
  the upstream means, and peak density reaches at least `1.5×` the upstream mean
  (a clear jump, not uniform flow).  Any NaN/Inf in the wall region also fails.

## 2D variants

`PARAM.in.hybrid_2d_*` run the same upstream on a 64 × 8 `d_i` domain with
128 × 16 cells (`dx = 0.5 d_i`, periodic in y), 100 ions per cell and
4 field sub-cycles, to `t = 20` (about three ion gyroperiods).  They differ only
in how the grid-scale noise of the centered Ohm's law is damped:

* **`hybrid_2d`**: the default upwind dissipation of `#HYBRIDUPWIND` (MC limiter
  with the grid-scale sensor, flow plus fast-wave speed);
* **`hybrid_2d_hyper`**: grid-mode `#HYPERRESISTIVITY` with `C_h = 0.01`, upwind
  off.

The table below predates the sensor (plain MC upwind); with the sensor the 2D
deck gives overshoot 1.63 and high-k fraction 0.059 on 2 ranks.

Without either term the deck does not survive: the shock grows a `|B|`
overshoot of about 3, grid-scale noise rises to 18 % of the y-varying field
power, the field energy grows 11× (8× when damped), and the noise erodes the
upstream.  At `t = 20` (overshoot = max of the y-averaged `|B|` over its
downstream mean; high-k = fraction of the y-varying B power above half the
Nyquist wavenumber):

| damping | shock x | overshoot | high-k | `Eb` growth |
| --- | --- | --- | --- | --- |
| none | 49.0 | 2.87 | 0.175 | 11.1× |
| hyper-resistivity `C_h = 0.003` | 46.0 | 1.63 | 0.025 | 8.0× |
| hyper-resistivity `C_h = 0.01` | 43.0 | 1.49 | 0.014 | 7.9× |
| upwind, `upwindCoefWhistler = 0` | 44.5 | 1.66 | 0.033 | 7.2× |
| upwind, `upwindCoefWhistler = 0.25` | 43.5 | 1.53 | 0.025 | 7.6× |
| upwind, `upwindCoefWhistler = 1` | 41.5 | 1.53 | 0.028 | 7.2× |

So the upwind term damps the shock as well as hyper-resistivity, without the
deck-dependent `C_h`.  In the 1D deck (`PARAM.in.hybrid`, `t = 5`), both
terms lower the field-energy growth from 5.0× to 3.0–3.3× and pass;
hyper-resistivity there gives the smaller overshoot (1.14 vs 1.30–1.36).

The 2D variants run the same checks as the 1D ones, with two changes:

* The upstream window is the rightmost **15 %** (x > 54), because the shock
  foot reaches x ≈ 50 by `t = 20`.
* **Grid noise damped:** overshoot ≤ **2.2** and high-k fraction ≤ **0.08**.

## Running

```
python3 tests/validate_tests.py --test=shock                    # all variants
python3 tests/validate_tests.py --test=shock.implicit           # full PIC only
python3 tests/validate_tests.py --test=shock.hybrid             # 1D hybrid only
python3 tests/validate_tests.py --test=shock.hybrid_2d          # 2D, default upwind
python3 tests/validate_tests.py --test=shock.hybrid_2d_hyper    # 2D, hyper-resistivity
```
