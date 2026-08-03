# FLEKS Tests

This directory contains standalone tests for the FLEKS (Flexible Exascale Kinetic
Simulator) particle-in-cell (PIC) solver, independent of SWMF coupling.

## Directory Structure

Each test case is contained within its own dedicated subdirectory containing a
`PARAM.in` configuration file and a README describing the expected behavior:

| Test              | Dir                | Description                                              | Readme |
|-------------------|--------------------|----------------------------------------------------------|--------|
| Beam instability  | `beam/`            | 1D ion beam EM instability: cyclotron wave growth,       | [README](beam/README.md)            |
|                   |                    | transverse B-field amplification, energy conservation    |        |
| Photoionization   | `photoionization/` | Chamberlain neutral profile with photoionization of      | [README](photoionization/README.md)  |
|                   |                    | H+O atmosphere                                           |        |
| Electron impact   | `electronimpact/`  | Voronov 1997 e-impact rate with hot electrons            | [README](electronimpact/README.md)   |
|                   |                    | (T ~ 100,000 K), only electron impact enabled            |        |
| Charge exchange   | `chargeexchange/`  | Constant CX cross-section with flowing solar wind ions,  | [README](chargeexchange/README.md)   |
|                   |                    | only charge exchange enabled                             |        |
| Hybrid whistler   | `hybrid_whistler/` | Hybrid PIC (kinetic ions / fluid electrons) whistler–     | [README](hybrid_whistler/README.md)  |
|                   |                    | Alfven wave; Hall term only (resistivity & e-pressure off)|       |
| Hybrid Ohm's law  | `hybrid_ohm/`      | Hybrid PIC full generalized Ohm's law: convection + Hall  | [README](hybrid_ohm/README.md)       |
|                   |                    | + resistive + electron-pressure-gradient, all active       |        |
| Free-stream       | `freestream/`      | 1D uniform free-stream, run with both solvers (full PIC and hybrid Hall-off) | [README](freestream/README.md) |
| Convection wave   | `hybrid_convection_wave/` | Hybrid PIC convection-only (`#HALLTERM F`) transverse wave advected by uniform flow; checks phase advance and amplitude conservation | [README](hybrid_convection_wave/README.md) |
| Light wave        | `lightwave/`       | 3D vacuum transverse EM (light) wave on a periodic AMR  | [README](lightwave/README.md)        |
|                   |                    | grid; energy-conservation check (needs `nLevMax >= 2`)   |        |
| Performance       | `performance/`     | Beam-based scaling benchmark (excluded from CI suite)    | — (see `validate_performance.py`)    |

### Ionization Parameter Commands

Each ionization process is enabled via a dedicated command in PARAM.in:

- **`#PHOTOIONIZATION`**: per-component rates at planet surface [s^-1],
  diluted geometrically as `(rPlanet / r)^2`
- **`#ELECTRONIMPACT`**: Voronov 1997 formula: `sigmav(T) = A*(T/EI)^K / [X+(T/EI)] * exp(-EI/T)` [cm^3/s], parameters per component
- **`#CHARGEEXCHANGE`**: constant cross-section: `sigmav(u) = sigmaCX * |u_i|` [cm^3/s], cross-section matrix `sigmaCX(neutral, ion)` [cm^2]; each neutral component exchanges with all ion species and the frequency is summed over ions

All three can be combined (as in `tests/photoionization/`) or tested individually (as in `tests/electronimpact/` and `tests/chargeexchange/`).

## Architecture

Ionization parameters are stored in `SourceInterface` (not `FluidInterface`) and read by `UserSource::read_param()` in `userfiles/ExoSource.h`. The Domain routes `#PHOTOIONIZATION`, `#ELECTRONIMPACT`, and `#CHARGEEXCHANGE` commands to the source object rather than to `FluidInterface`. This keeps the MHD coupling layer uncluttered by ionization-specific data.

## Building the Test Executable

All standalone tests run with a **single FLEKS executable** built with the
Exosphere user source (`-u=Exo`, needed by the ionization tests) and two grid
levels (`-lev=2`, needed by the `lightwave` AMR test):

```bash
cd PC/FLEKS            # (or the FLEKS root)
./Config.pl -lev=2 -u=Exo
make -j4
```

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
python3 tests/validate_tests.py -n 2

# Or equivalently:
python3 tests/validate_tests.py --nprocs 2
```

When `-n 1` (or the flag is omitted), the executable is invoked directly as
`./FLEKS.exe` without `mpirun`. When `-n N` with `N > 1`, it uses
`mpirun -n N ./FLEKS.exe`.

The `--test NAME` (or `--test=NAME`) option selects a single test to run from
the available test subdirectories (`beam`, `photoionization`, `electronimpact`,
`chargeexchange`). If the given name does not match any test, the script exits
with an error listing the available tests. When the flag is omitted, all tests
are run (the default behavior). The flag may be combined with `-n`/`--nprocs`.

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
