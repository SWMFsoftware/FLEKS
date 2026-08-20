# 1D Oblique Shock Test

## Physics

A uniform magnetized hybrid plasma (kinetic ions + massless fluid electrons)
streams along `-x` (upstream `ux = -713 km/s`) through a 1D domain (64 cells in
x, `xMin = 0`, `xMax = 64`).  The shock reflects at the left wall (`x = lo`):

* **left wall (`-x`):** particle `reflect` (specular mirror) + **field `fixed`**
  (`BC::fixed`) — the ghost EM field is pinned to the upstream `#INFLOW` state
  (B = upstream B, E = -u_in x B_in).  This is the same EM machinery as the
  `inflow` field BC.
* **right wall (`+x`):** `inflow` field + `inflow` particles (open, upstream-fed).

The upstream magnetic field is inclined **45 deg to the shock normal (x)**:

```
B0      = 1.044e-8 T
Bx = Bz = B0 * cos(45) = 7.382e-9 T     (By = 0)
```

so this is an *oblique* shock with a non-zero normal field component `Bx != 0`.
That distinguishes it from `run_shock_1d` (perpendicular, `Bx = 0`): the 2D
oblique shock showed that a PEC/conducting wall at the reflection face is
incompatible with a flowing magnetized plasma, because forcing `E_t = 0`
contradicts the bulk motional E = -u x B and drives a field-energy runaway.
Hence the left field uses `fixed` (pinned to upstream), not `conducting`.

## Initial / upstream motional E

With `uy = uz = 0`, `By = 0`, only the `y` component of `E = -u x B` is
non-zero:

```
Ey = -ux * Bz = -(-713e3 m/s) * (7.382e-9 T) = +5.263e-3 V/m
```

This is set in `#UNIFORMSTATE` (interior) and enforced at the `fixed`/`inflow`
faces by `Pic::apply_inflow_wall`.

## The `#INFLOW` command

```
#INFLOW
7.382e-9    bx [T]
0.0         by [T]
7.382e-9    bz [T]
5.0         rho [amu/cc]
-713.0      ux [km/s]
0.0         uy [km/s]
0.0         uz [km/s]
6.28e5      T [K]
```

Because the left face is `fixed` (not `inflow`) for particles, outgoing particles
that cross `-x` outward are removed/tallied by `reflect_or_delete_particle`
(BC::reflect mirror), while the *field* on that face is pinned to upstream via the
`inflow`-EM path.  A `fixed` face is a no-op for the field without an `#INFLOW`
block (it then falls back to the open-boundary copy).

## Running

```bash
python3 tests/validate_tests.py --test=shock_1d_oblique
```

## Validation

The oblique shock is *not* a uniform steady state, so unlike `bc_inflow` the
checker does not require a uniform field everywhere.  It guards against exactly
the failure the `fixed` BC was introduced to fix:

* **No blow-up:** ion kinetic energy `Epart` and magnetic energy `Eb` stay finite
  and bounded.  A `conducting` left wall drives a field-energy runaway (Eb ->
  ~1e18), so this is the decisive regression check.
* **Reflect-field boundary holds the upstream state:** near the left (`-x`) face
  the field components `Bx`, `Bz` stay at the upstream `#INFLOW` value (`7.382e-9`)
  with no spurious erosion, and the motional `Ey` is maintained.  (A PEC wall would
  zero `Ey` and disagree with `-u x B`.)
* **Inflow (+x) face maintains upstream:** the upstream density `rhoS0` and bulk
  velocity `uxS0` at the rightmost cells stay within tolerance of the prescribed
  `#INFLOW` state (no draining toward vacuum, no spurious source).
* **No spurious electric field:** mean `Ex`, `Ez` stay near zero; `Ey` tracks the
  upstream motional value.

If the run blows up (non-finite energies) the test FAILS loud; otherwise it PASS.
