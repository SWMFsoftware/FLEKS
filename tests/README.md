# FLEKS Tests

This directory contains standalone tests for the FLEKS (Flexible Exascale Kinetic
Simulator) particle-in-cell (PIC) solver, independent of SWMF coupling.

## Directory Structure

Each test case is contained within its own dedicated subdirectory containing a
`PARAM.in` configuration file and a README file:

| Test                  | Dir                   | Description                                                  | Readme |
|-----------------------|-----------------------|--------------------------------------------------------------|--------|
| Beam instability      | `beam/`               | 1D ion-beam cyclotron instability; energy conservation        | [README](beam/README.md) |
| Photoionization       | `photoionization/`    | Chamberlain neutral atmosphere with photoionization           | [README](photoionization/README.md) |
| Electron impact       | `electronimpact/`     | Voronov 1997 electron-impact ionization of hot electrons      | [README](electronimpact/README.md) |
| Charge exchange       | `chargeexchange/`     | Constant-cross-section charge exchange with solar wind ions   | [README](chargeexchange/README.md) |
| Whistler wave         | `whistler/`           | Whistler–Alfven wave (Hall term)                              | [README](whistler/README.md) |
| Ohm's law             | `ohm/`                | Full generalized Ohm's law check | [README](ohm/README.md) |
| Free-stream           | `freestream/`         | 1D uniform free-stream                                        | [README](freestream/README.md) |
| Light wave            | `lightwave/`          | 3D vacuum light wave on a periodic AMR grid | [README](lightwave/README.md) |
| PCAI                  | `pcai/`               | Proton-cyclotron anisotropy instability (`T_perp/T_par=3`, `gamma/Omega_ci=0.162`) | [README](pcai/README.md) |
| Reconnection          | `reconnection/`       | Current-sheet reconnection: Fadeev, classic GEM challenge, and asymmetric reconnection (<30s serial CI) | [README](reconnection/README.md) |
| AMR reconnection      | `reconnection_amr/`   | Fadeev current-sheet reconnection on a two-level AMR grid | [README](reconnection_amr/README.md) |
| Reflecting & PEC BC   | `bc_reflecting/`      | Specular reflecting particle walls + conducting (PEC) field walls (4 variants: full/hybrid fields/particles) | [README](bc_reflecting/README.md) |
| Absorbing BC          | `bc_absorb/`          | Absorbing field + particle boundaries (4 variants: full/hybrid fields/particles) | [README](bc_absorb/README.md) |
| Wave injection        | `bc_wave/`            | Grouped wave-injection tests: mono Bz wave + shear Alfvén wave via `#WAVEBC` (one `PARAM.in.<suffix>` per variant) | [README](bc_wave/README.md) |
| Oblique shock         | `shock/`              | 1D oblique magnetized shock | [README](shock/README.md) |
| Inner body            | `body/`               | Absorbing spherical inner boundary (`#BODY`): particle absorption, empty interior, wake | [README](body/README.md) |
| Performance           | `performance/`        | Beam-based scaling benchmark         | — (see `validate_performance.py`) |

### Ionization Parameter Commands

Each ionization process is enabled via a dedicated command in PARAM.in:

- **`#PHOTOIONIZATION`**: per-component unattenuated rates [s^-1]. The
  production rate is `n_i(r) * nu0_i * A(r)`, where the attenuation `A(r)` is
  the `#SHADOWCYLINDER` mask (0 inside the cylinder, 1 outside) by default, or
  `exp(-tau/mu)` when `#OPTICALDEPTH` is used. The `nu0` values follow the
  BATSRUS `ModUserMars` `Rate_I` convention.
- **`#ELECTRONIMPACT`**: Voronov 1997 formula: `sigmav(T) = A*(T/EI)^K / [X+(T/EI)] * exp(-EI/T)` [cm^3/s], parameters per component
- **`#CHARGEEXCHANGE`**: constant cross-section: `sigmav(u) = sigmaCX * |u_i|` [cm^3/s], cross-section matrix `sigmaCX(neutral, ion)` [cm^2]; each neutral component exchanges with all ion species and the frequency is summed over ions
- **`#SHADOWCYLINDER`** / **`#OPTICALDEPTH`**: two mutually exclusive EUV
  attenuation models (FLEKS aborts if both are present). `#OPTICALDEPTH`
  implements the BATSRUS optical-depth model `exp(-tau/mu)` with
  `tau = sum_c n_c * sigma_c * H_c`, optionally with the curved-atmosphere
  Chapman function

### Units in the ionization decks

The two length conventions coexist and must not be mixed:

- **grid** — `#GEOMETRY` (and any other grid-oriented command such as
  `#REGION`) is in *code units*, i.e. multiples of `#NORMALIZATION lNormSI`
  metres. Standalone, the numbers are used verbatim as AMReX coordinates.
- **physics input** — `#BODYSIZE`, `#EXOSPHERE` (`n0`, `H0`, `T0`),
  `#PHOTOIONIZATION`, `#ELECTRONIMPACT`, `#CHARGEEXCHANGE`,
  `#SHADOWCYLINDER` and `#OPTICALDEPTH` are all **SI**.
  `UserSource` converts the code-unit cell position to metres (`No2SiL`) before
  it evaluates the neutral profile or the shadow/optical-depth geometry.

This is also what a GM-coupled run looks like: there `No2SiL` is the length
of one GM length unit (one planetary radius), so the deck of a coupled run and
the deck of a standalone run describe the same setup as long as each keeps to
its own convention. A validator must therefore use the plot unit of the
`PLANETARY` output, namely one `#BODYSIZE` radius, to turn plot coordinates
back into metres.

## Architecture

Ionization parameters are stored in `SourceInterface` and read by `UserSource::read_param()` in `userfiles/ExoSource.h`. The Domain routes specific commands to the source object rather than to `FluidInterface`. This keeps the MHD coupling layer uncluttered by ionization-specific data.

## Building the Test Executable

### 3D AMReX build (default)

Most standalone tests run with a **single FLEKS executable** built with the
Exosphere user source (`-u=Exo`, needed by the ionization tests) and two grid
levels (`-lev=2`, needed by the AMR tests).  This is the default build and links
the **3D** AMReX library:

```bash
cd PC/FLEKS            # (or the FLEKS root)
./Config.pl -lev=2 -u=Exo
make -j4
```

### 2D AMReX build

The suite can also be built against the **true-2D AMReX library**
(`AMReX_SPACEDIM = 2`):

```bash
cd PC/FLEKS            # (or the FLEKS root)
./Config.pl -amrex2d -lev=2 -u=Exo
make -j4
```

The runner auto-detects the built AMReX dimension and skips any test that needs
the other one (e.g. the true-3D `lightwave` and the source-term /
hyper-resistivity tests under a 2D build).  Use `-amrex3d` to switch back to the
default 3D build.

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
`chargeexchange`, ...), or a specific variant such as `--test=reconnection.forcefree`.
If the given name does not match any test, the script exits with an error listing
the available tests. When the flag is omitted, all standard tests are run (the default behavior).
The flag may be combined with `-n`/`--nprocs`.

When a test directory contains multiple `PARAM.in.<variant>` files, the runner can execute
each variant (e.g. `PARAM.in.hybrid`, `PARAM.in.forcefree`). Certain long-running or research
tests (such as `reconnection.forcefree`, `beam.instability`, and `iaw.landau`) are designated as
expensive and excluded from the default quick test run. They can be executed specifically via:
```bash
# Run a specific variant directly:
python3 tests/validate_tests.py --test=reconnection.forcefree
python3 tests/validate_tests.py --test=beam.instability
python3 tests/validate_tests.py --test=iaw.landau

# Or include all expensive tests in the full suite:
python3 tests/validate_tests.py --include-expensive
# (or python3 tests/validate_tests.py --all)
```

### Performance Benchmark

```bash
python3 tests/validate_performance.py
```

The script benchmarks the full-PIC beam test (`performance/PARAM.in`),
the hybrid-PIC whistler test (`performance/PARAM.in.hybrid`), and the
particle tracker test (`performance/PARAM.in.pt`), and writes the results to
`tests/performance_summary.md`.

### Memory Regression (master vs PR)

Every standalone run already reports memory, so no new instrumentation is
needed: the AMReX TinyProfiler report gives the allocation count and peak bytes
per `BL_PROFILE` region, and the FLEKS load-balance report gives process RSS.
Two scripts capture and compare them — timings are not compared:

```bash
python3 tests/capture_memory.py --out mine.json                # step 1: your tree
python3 tests/capture_memory.py --out base.json --ref master   # step 1: reference
python3 tests/compare_memory.py base.json mine.json            # step 2: compare
```

`compare_memory.py` prints a table and exits non-zero on a regression. Run
either script with `--help` for the options; the ones worth knowing are
`capture_memory.py --list` (which decks are captured — beam,
performance.hybrid, reconnection.fadeev_pic, shock and 2-rank beam, ~25 s in
total), `--timeout` (per-test wall-clock limit, 900 s by default: a deck that
diverges into an endless loop is killed and recorded as an error instead of
hanging the job), and `--verify`, which reports how reproducible memory is on
your machine and is what `--rss-tol` should be tuned from. The report
identifies the exact baseline and candidate commits, includes absolute and
relative changes, and warns if both profiles come from the same commit.

Arena counts and peak bytes are exact integers across runs, so any increase
fails. RSS drifts by up to 0.4 MB between runs, so it gets a tolerance
(`--rss-tol`, default 2 MB) and only the per-rank maximum is gated.

CI runs both steps on the same runner and comments on the PR:
`.github/workflows/memory_test.yml`.

Two notes:

* Do **not** enable `#MEMORY` to obtain the RSS series — it is not a reporting
  switch. Every `dnMemory` cycles it also calls `Pic::free_memory()`, which runs
  `CArena::freeUnused()`, `ShrinkToFit()` on every particle container and
  `malloc_trim(0)`, perturbing exactly what is being measured. The report is
  printed anyway whenever `doReport` is set.
* The arena profiler only sees `MultiFab`/`FArrayBox`/particle-tile traffic. RSS
  covers the rest (`std::vector`, `new`) but cannot attribute it to a function.

To inspect a single run report outside the regression workflow:

```bash
python3 tests/profiler.py prof.txt --top 10          # timing + arena memory
python3 tests/profiler.py run.log --load-balance     # RSS series
```
