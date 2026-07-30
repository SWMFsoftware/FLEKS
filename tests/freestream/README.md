# Free-stream — uniform-equilibrium field-solver test

**Test directory:** `tests/freestream/`

A 1D periodic **free-stream**: the initial condition is a *perfectly uniform*
plasma (uniform number density, zero bulk velocity, uniform guide field `Bx0`
along x) on a coarse 32-cell periodic grid. No `#TESTCASE` is set, so there is
no wave seed — the only structure is PIC shot-noise. A correct field solver
must keep this trivial state **exactly steady**: any growth in the magnetic or
particle energy is purely numerical error or shot-noise.

This is the cheapest regression test and the baseline for the seeded-wave
tests (`hybrid_whistler`, `hybrid_ohm`, `hybrid_convection_wave`): those
*seed* a perturbation on top of this same equilibrium, while the free-stream
checks the equilibrium itself does not move.

## Both solvers, one test directory

The directory holds **two parameter files** that differ only in the
field-solver block:

| File | Field solver | Ohm's-law terms |
| --- | --- | --- |
| `PARAM.in` | Full PIC implicit Maxwell/GMRES (`#SOLVEEM T`, `#HYBRIDPIC F`) | n/a (full E–B solve) |
| `PARAM.in.hybrid` | Hybrid (Ohm + Faraday) (`#SOLVEEM F`, `#HYBRIDPIC T`, `#HALLTERM F`) | convection only (no sub-cycling) |

When the test runner (`tests/validate_tests.py`) reaches this test it runs the
directory **twice**, once with `PARAM.in` and once with `PARAM.in.hybrid`.
Everything outside the solver block — grid, plasma state, normalization,
timestepping, particles, output cadence — is identical between the two files,
so the two runs form a direct cross-comparison of the two solvers on the same
trivial equilibrium.

## Validation

Both variants are checked by `validate_hybrid`: the magnetic energy `Eb` and ion
energy `Epart` must be steady (within a small growth/tolerance factor) from the
first to the last logged step. Expected: `Eb` and `Epart` constant to ~1e-6
(full PIC) / ~1e-4 (hybrid, residual PIC shot-noise) relative.

## Physics scales

- Normalization chosen so `d_i ~ v_A ~ 1`: `lNormSI = 1e5 m` (= `d_i`),
  `uNormSI = 5e4 m/s` (= `v_A`).
- Single kinetic ion species (`q = m = 1`), `n = 5 /cc`, uniform guide field
  `Bx0 = 5.2195e-9 T` (~52.2 nT), `T = 10000 K`.
- Grid spans `Lx = 6.4 d_i` with 32 cells (periodic, 1 cell in y/z).
- Fixed `dt = 0.02` for 500 steps (`tmax = 10 s`).

## Relationship to other tests

- `hybrid_whistler/`, `hybrid_ohm/`, `hybrid_convection_wave/` — seeded wave
  tests built on the same equilibrium.
- `whistler_solveem/` — the full-PIC SOLVEEM *wave* test (seeded wave, not a
  free-stream).
