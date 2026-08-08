# HYBRID_WHISTLER

Parallel whistler–Alfvén wave test in FLEKS. This directory now runs **two
variants of the same eigenmode**:

| Variant | Config file | Field solver | Species |
| --- | --- | --- | --- |
| **Full PIC** | `PARAM.in` | Maxwell/GMRES (`solveEM=T`) | kinetic ions **and** kinetic electrons |
| **Hybrid** | `PARAM.in.hybrid` | generalized Ohm's law (`useHybridPIC=T`) | kinetic ions + fluid electrons |

Both variants are seeded with the same x-aligned, transverse, circularly
polarized whistler wave (`#TESTCASE HybridWave`, `frac=0.02`) and the same
Alfvén velocity kick, so they exercise the identical whistler eigenmode with
two different field/B solvers. This lets the full-PIC dispersive result (with
kinetic electrons) be compared directly against the Hall-only hybrid result.

## Physics / normalization

All quantities are in normalized CGS/code units chosen so ion scales are O(1):

* Box length `Lx = 6.4 d_i` (`nCellX = 64`, `maxBlockSizeX = 64`).
* `lNormSI = d_i ≈ 1.0e5 m`, `uNormSI = v_A ≈ 5.0e4 m/s`, so `d_i ~ 1`,
  `v_A ~ 1`, and (with `q/(m c) = 1` for the proton) `B = Ω_i ~ 1`.
* Timestep `dt = 0.02` → `dt·Ω_i = 0.02 ≪ 1`, numerically stable for the
  whistler mode (bounded by `Ω_i`, since
  `ω/Ω_i = (k d_i)² / (1 + (k d_i)²)`).
* Physical inputs are solar-wind-like: `n = 5 amu/cc`, `B = 5.22 nT`,
  `v_A ≈ 51 km/s`, `d_i ≈ 102 km`.
* `TimeMax = 13.0` ≈ one whistler period for the `n = 1` mode.

In the full-PIC variant the kinetic electron species has
`q = -1`, `m = 1/1836` (proton mass ratio), so the PIC Maxwell solve carries
the full electron+ion whistler dispersion. The guide field and electron density
match the ions (`n_e = n_i`).
This real electron mass does **not** affect the timestepping in the test:
the timestep here is fixed (`useFixedDt = T`, `dt = 0.02`) and is set by the
ion-scale whistler frequency (`ω/Ω_i ≲ 1`, i.e. `dt·Ω_i ≪ 1`).

## Running

Both variants are driven automatically by the test runner
(`tests/validate_tests.py`): `PARAM.in` runs as the **FULL PIC** variant and
`PARAM.in.hybrid` as the **HYBRID** variant in the same directory.

Manually:

```bash
# Full PIC
make -C <build> run PARAMPATH=tests/whistler/PARAM.in
# Hybrid
make -C <build> run PARAMPATH=tests/whistler/PARAM.in.hybrid
```

## Validation

`validate.py` runs the whistler field checks for **both** variants: it measures
the whistler frequency from the transverse B perturbation and compares it with
the analytic dispersion `ω/Ω_i = (k d_i)² / (1 + (k d_i)²)`, and checks the
Alfvén-speed propagation. For the **hybrid** variant it additionally runs the
shared hybrid energy-log checks (see the [hybrid README](../_shared/hybrid.md)).
