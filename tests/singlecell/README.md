# Single-Cell Hybrid Hall Test (`tests/singlecell/`)

## Purpose

Regression test for the hybrid field solver on a degenerate single-cell domain.
The grid is exactly **one cell** with periodic boundaries in all three
directions.  On a single cell every finite-difference gradient is identically
zero, so `curl B == 0`.  Because the Hall term of the generalized Ohm's law is

```
E_Hall = - (J x B) / (n_e e),    J = curl B / (4 pi),
```

the Hall term is **exactly zero**.  With a uniform plasma (`U_i = 0`) the
convective term `U_i x B` also vanishes, so the total electric field `E = 0` and
Faraday's law gives `dB/dt = -curl E = 0`: the guide field is frozen.

This guards against a class of bugs where the Hall or convective term is
evaluated even when it must vanish by construction (e.g. a spurious single-cell
contribution that makes `Eb` drift or blow up).

## Setup

* 1x1x1 periodic Cartesian grid (`#NCELL 1 1 1`, `#MAXBLOCKSIZE 1 1 1`).
* Hybrid PIC (`#SOLVEEM F`, `#HYBRIDPIC T`), Hall term **ON**
  (`#HALLTERM T`) — the test proves the term contributes exactly zero.
* Uniform ion loading (4 particles/cell), uniform guide field `Bx0 ~ 1`
  (code units, `B_code ~ 1` from the `lNormSI=1e5` / `uNormSI=5e4` device scale).
* Fixed `dt = 0.02`, runs to `TimeMax = 2.0`.

## Validation

`validate_singlecell` (in `tests/validate_tests.py`) reads the pic-log history
and requires:

1. **`Eb` conserved** to round-off (`0.9999 < ratio < 1.0001`) — no spurious
   Hall/convective evolution.
2. **`Ee` (electric field energy) ~ 0** throughout (`Ee < 1e-2 * Eb`) — confirms
   the field is truly frozen, not merely energy-conserving (a non-dispersive
   propagating wave would also conserve `Eb` but have `Ee ~ Eb`).  The 1e-2
   relative floor accommodates single-cell round-off (residual Ee ~ 1e-3 * Eb).

A passing run shows `Eb ratio: 1.000000` and `Ee (electric, max): ~0`.

## Run

```
python3 tests/validate_tests.py --test=singlecell
```
