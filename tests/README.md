# FLEKS Unit Tests

This directory is the future home for standalone unit testing of the FLEKS codebase (independent of the overarching SWMF macro regression tests).

## Framework Guidance
The recommended framework for internal kinetic math and grid testing is **GoogleTest (gtest)** or **Catch2**. 

When adding tests here in the future:
1. Initialize the framework of your choice (e.g., via a standalone CMake build system targeting this folder).
2. Tests should execute entirely isolated from SWMF Fortran wrappers (focus purely on `src/` classes like `Particles.cpp`, `LinearSolver.cpp`, and mathematical `GridUtility` routines).
3. If leveraging AMReX data arrays, ensure your test runner initializes AMReX via `amrex::Initialize()` before the test suite and calls `amrex::Finalize()` during test tear-down.

| Test              | Dir                 | Description                                              |
|-------------------|---------------------|----------------------------------------------------------|
| Beam instability  | `beam/`             | 1D ion beam EM instability: cyclotron wave growth,       |
|                   |                     | transverse B-field amplification, energy conservation    |
| Exosphere profile | `exosphere/`        | Chamberlain neutral profile with all 3 ionization        |
|                   |                     | processes (photo, impact, CX) for H+O atmosphere         |
| Electron impact   | `electron_impact/`  | Voronov 1997 e-impact rate with hot electrons            |
|                   |                     | (T ~ 100,000 K), only electron impact enabled            |
| Charge exchange   | `charge_exchange/`  | Constant CX cross-section with flowing solar wind ions,  |
|                   |                     | only charge exchange enabled                             |
| Performance       | `performance/`      | Beam-based scaling benchmark (excluded from CI suite)    |

### Ionization Parameter Commands

Each ionization process is enabled via a dedicated command in PARAM.in:

- **`#PHOTOIONIZATION`**: per-component rates at planet surface [s^-1],
  diluted geometrically as `(rPlanet / r)^2`
- **`#ELECTRONIMPACT`**: Voronov 1997 formula: `sigmav(T) = A*(T/EI)^K /
  [X+(T/EI)] * exp(-EI/T)` [cm^3/s], parameters per component
- **`#CHARGEEXCHANGE`**: constant cross-section: `sigmav(u) = sigmaCX * |u_i|`
  [cm^3/s], per-component sigmaCX [cm^2]

All three can be combined (as in `tests/exosphere/`) or tested individually
(as in `tests/electron_impact/` and `tests/charge_exchange/`).

## Architecture

Ionization parameters are stored in `SourceInterface` (not `FluidInterface`)
and read by `UserSource::read_param()` in `userfiles/ExoSource.h`. The Domain
routes `#PHOTOIONIZATION`, `#ELECTRONIMPACT`, and `#CHARGEEXCHANGE` commands
to the source object rather than to `FluidInterface`. This keeps the MHD
coupling layer uncluttered by ionization-specific data.

## Running the Tests

A unified Python runner is provided to dynamically discover, run, and validate
the standalone tests:

```bash
python3 tests/validate_tests.py
```

The script:
1. Scans subdirectories of `tests/` for a `PARAM.in` file,
   excluding `performance/`.
2. Creates the execution directory `run_test/` with necessary symlinks to the
   FLEKS executable and post-processing tools.
3. Copies `PARAM.in` from the test folder, executes `./FLEKS.exe`, and runs
   `./PostProc.pl`.
4. Parses diagnostic outputs and validates the physics checks for each test.
5. Writes a summary table to `tests/summary.md`.
