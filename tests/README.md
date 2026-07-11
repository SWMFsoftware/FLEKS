# FLEKS Tests

This directory contains standalone tests for the FLEKS (Flexible Exascale Kinetic
Simulator) particle-in-cell (PIC) solver, independent of SWMF coupling.

## Directory Structure

Each test case is contained within its own dedicated subdirectory containing a
`PARAM.in` configuration file and a README describing the expected behavior:

| Test              | Dir                | Description                                              |
|-------------------|--------------------|----------------------------------------------------------|
| Beam instability  | `beam/`            | 1D ion beam EM instability: cyclotron wave growth,       |
|                   |                    | transverse B-field amplification, energy conservation    |
| Photoionization   | `photoionization/` | Chamberlain neutral profile with photoionization of      |
|                   |                    | H+O atmosphere                                           |
| Electron impact   | `electronimpact/`  | Voronov 1997 e-impact rate with hot electrons            |
|                   |                    | (T ~ 100,000 K), only electron impact enabled            |
| Charge exchange   | `chargeexchange/`  | Constant CX cross-section with flowing solar wind ions,  |
|                   |                    | only charge exchange enabled                             |
| Hybrid wave       | `hybrid_wave/`     | Hybrid PIC (kinetic ions / fluid electrons) whistler–     |
|                   |                    | Alfven wave; validates the Hall term & CGS units          |
| Performance       | `performance/`     | Beam-based scaling benchmark (excluded from CI suite)    |

### Ionization Parameter Commands

Each ionization process is enabled via a dedicated command in PARAM.in:

- **`#PHOTOIONIZATION`**: per-component rates at planet surface [s^-1],
  diluted geometrically as `(rPlanet / r)^2`
- **`#ELECTRONIMPACT`**: Voronov 1997 formula: `sigmav(T) = A*(T/EI)^K / [X+(T/EI)] * exp(-EI/T)` [cm^3/s], parameters per component
- **`#CHARGEEXCHANGE`**: constant cross-section: `sigmav(u) = sigmaCX * |u_i|` [cm^3/s], cross-section matrix `sigmaCX(neutral, ion)` [cm^2]; each neutral component exchanges with all ion species and the frequency is summed over ions

All three can be combined (as in `tests/photoionization/`) or tested individually (as in `tests/electronimpact/` and `tests/chargeexchange/`).

## Architecture

Ionization parameters are stored in `SourceInterface` (not `FluidInterface`) and read by `UserSource::read_param()` in `userfiles/ExoSource.h`. The Domain routes `#PHOTOIONIZATION`, `#ELECTRONIMPACT`, and `#CHARGEEXCHANGE` commands to the source object rather than to `FluidInterface`. This keeps the MHD coupling layer uncluttered by ionization-specific data.

## Running the Tests

### Standard Test Suite

A unified Python runner is provided to dynamically discover, run, and validate
the standalone tests:

```bash
# Run all tests in serial mode (no MPI):
python3 tests/validate_tests.py

# Run a single test by name (e.g. beam):
python3 tests/validate_tests.py --test=beam

# Equivalently, with a space:
python3 tests/validate_tests.py --test beam

# Run with N MPI processes:
python3 tests/validate_tests.py -n 4

# Or equivalently:
python3 tests/validate_tests.py --nprocs 4
```

When `-n 1` (or the flag is omitted), the executable is invoked directly as
`./FLEKS.exe` without `mpirun`. When `-n N` with `N > 1`, it uses
`mpirun -n N ./FLEKS.exe`.

The `--test NAME` (or `--test=NAME`) option selects a single test to run from
the available test subdirectories (`beam`, `photoionization`, `electronimpact`,
`chargeexchange`). If the given name does not match any test, the script exits
with an error listing the available tests. When the flag is omitted, all tests
are run (the default behavior). The flag may be combined with `-n`/`--nprocs`.

The script:
1. Scans subdirectories of `tests/` for a `PARAM.in` file, excluding `performance/`.
2. Creates the execution directory `run_test/` with necessary symlinks to the
   FLEKS executable and post-processing tools.
3. Copies `PARAM.in` from the test folder, executes `./FLEKS.exe` (or via
   `mpirun`), and runs `./PostProc.pl`.
4. Parses diagnostic outputs and validates the physics checks for each test.
5. Writes a summary table to `tests/summary.md`.

### Performance Benchmark

```bash
python3 tests/validate_performance.py
```

The script:
1. Runs the beam test with 1 and 2 MPI processes (3 runs each for statistical
   robustness).
2. Parses AMReX TinyProfiler output to extract particle mover and field solver
   timings.
3. Computes particle-step rates (μs/part-step) and parallel speedup.
4. Validates against baseline targets and writes a report to
   `tests/performance_summary.md`.

## Test Details

Each test subdirectory contains its own README with detailed physics setup,
expected results, and instructions for manual execution:

| Subdirectory       | Readme |
|--------------------|--------|
| `beam/`            | [README](beam/README.md)             |
| `photoionization/` | [README](photoionization/README.md)  |
| `electronimpact/`  | [README](electronimpact/README.md)   |
| `chargeexchange/`  | [README](chargeexchange/README.md)   |
| `performance/`     | — (see `validate_performance.py`)    |
