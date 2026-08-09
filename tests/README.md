# FLEKS Tests

This directory contains standalone tests for the FLEKS (Flexible Exascale Kinetic
Simulator) particle-in-cell (PIC) solver, independent of SWMF coupling.

## Directory Structure

Each test case is contained within its own dedicated subdirectory containing a
`PARAM.in` configuration file and a README file:

| Test              | Dir                | Description                                              | Readme |
|-------------------|--------------------|----------------------------------------------------------|--------|
| Beam instability  | `beam/`            | 1D ion beam EM instability: cyclotron wave growth,       | [README](beam/README.md)            |
|                   |                    | transverse B-field amplification, energy conservation.   |        |
|                   |                    | Run with both solvers (full PIC and hybrid PIC)          |        |
| Photoionization   | `photoionization/` | Chamberlain neutral profile with photoionization of      | [README](photoionization/README.md)  |
|                   |                    | H+O atmosphere                                           |        |
| Electron impact   | `electronimpact/`  | Voronov 1997 e-impact rate with hot electrons            | [README](electronimpact/README.md)   |
|                   |                    | (T ~ 100,000 K), only electron impact enabled            |        |
| Charge exchange   | `chargeexchange/`  | Constant CX cross-section with flowing solar wind ions,  | [README](chargeexchange/README.md)   |
|                   |                    | only charge exchange enabled                             |        |
| Whistler wave     | `whistler/`        | Whistler–Alfven wave; full PIC (ions+electrons) + hybrid  | [README](whistler/README.md)         |
|                   |                    | Alfven wave; Hall term only (resistivity & e-pressure off)|       |
| Hybrid Ohm's law  | `ohm/`             | Hybrid PIC full generalized Ohm's law: convection + Hall  | [README](ohm/README.md)              |
|                   |                    | + resistive + electron-pressure-gradient, all active       |        |
| Free-stream       | `freestream/`      | 1D uniform free-stream, run with both solvers (full PIC and hybrid Hall-off) | [README](freestream/README.md) |
| Light wave        | `lightwave/`       | 3D vacuum transverse EM (light) wave on a periodic AMR  | [README](lightwave/README.md)        |
|                   |                    | grid; energy-conservation check (needs `nLevMax >= 2`)   |        |
| PCAI              | `pcai/`            | Hybrid PIC proton-cyclotron anisotropy instability:     | [README](pcai/README.md)             |
|                   |                    | bi-Maxwellian `T_perp/T_par = 3`, `beta_par = 1`;       |        |
|                   |                    | transverse-wave growth rate `gamma/Omega_ci = 0.162`    |        |
| Performance       | `performance/`     | Beam-based scaling benchmark (excluded from CI suite)    | — (see `validate_performance.py`)    |

### Ionization Parameter Commands

Each ionization process is enabled via a dedicated command in PARAM.in:

- **`#PHOTOIONIZATION`**: per-component rates at planet surface [s^-1],
  diluted geometrically as `(rPlanet / r)^2`
- **`#ELECTRONIMPACT`**: Voronov 1997 formula: `sigmav(T) = A*(T/EI)^K / [X+(T/EI)] * exp(-EI/T)` [cm^3/s], parameters per component
- **`#CHARGEEXCHANGE`**: constant cross-section: `sigmav(u) = sigmaCX * |u_i|` [cm^3/s], cross-section matrix `sigmaCX(neutral, ion)` [cm^2]; each neutral component exchanges with all ion species and the frequency is summed over ions

## Architecture

Ionization parameters are stored in `SourceInterface` and read by `UserSource::read_param()` in `userfiles/ExoSource.h`. The Domain routes specific commands to the source object rather than to `FluidInterface`. This keeps the MHD coupling layer uncluttered by ionization-specific data.

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
python3 tests/validate_tests.py --nprocs 2

# Show detailed per-check diagnostics (energy numbers, ratios, etc.):
python3 tests/validate_tests.py --verbose
python3 tests/validate_tests.py -v
```

When `-n 1` (or the flag is omitted), the executable is invoked directly as
`./FLEKS.exe` without `mpirun`. When `-n N` with `N > 1`, it uses
`mpirun -n N ./FLEKS.exe`.

The `--test NAME` (or `--test=NAME`) option selects a single test to run from
the available test subdirectories (`beam`, `photoionization`, `electronimpact`,
`chargeexchange`, ...). If the given name does not match any test, the script
exits with an error listing the available tests. When the flag is omitted, all
tests are run (the default behavior). The flag may be combined with `-n`/`--nprocs`.

When a test directory contains both `PARAM.in` and `PARAM.in.hybrid` (currently
`beam/` and `freestream/`), the runner executes the test once per field solver,
listing both variants in the summary table (e.g. `BEAM` and `BEAM (HYBRID)`).

### Performance Benchmark

```bash
python3 tests/validate_performance.py
```

The script benchmarks both the full-PIC beam test (`performance/PARAM.in`) and
the hybrid-PIC whistler test (`performance/PARAM.in.hybrid`), and writes the
results to `tests/performance_summary.md`.
