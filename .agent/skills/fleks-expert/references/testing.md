# FLEKS Testing Reference

> Canonical: `tests/README.md` (catalogue and runner options) and
> `doc/DEVELOPING.md` §6. This page is the agent-oriented summary, including
> the build configurations each test needs.

FLEKS has two independent test systems: **standalone tests** (no SWMF, under
`tests/`) and **SWMF coupled tests** (GM-PC / MHD-AEPIC, run from the SWMF
root). Use standalone tests for solver-level changes (fields, particles, BCs,
sources); use coupled tests for interface changes (`FleksInterface.cpp`,
`PC_wrapper.f90`, `GM_couple_pc.f90`).

## Standalone tests (`tests/`)

One directory per scenario, each with a `PARAM.in`, a `README.md` and (for most
tests) its own `validate.py` holding that test's validators. The unified runner
`tests/validate_tests.py` keeps only the shared infrastructure and dynamically
loads the `validate.py` of the test it is about to run, so one test never
imports another's code. Code shared by several tests lives in `tests/_shared/`.

| Test | Dir | Physics |
|---|---|---|
| Beam instability | `beam/` | Ion beam cyclotron instability; full PIC + hybrid |
| Whistler wave | `whistler/` | Whistler–Alfvén wave (Hall term) |
| Ion acoustic wave | `iaw/` | Electron-pressure branch + ion Landau damping |
| Hybrid Ohm's law | `ohm/` | Convection + Hall + resistive + electron-pressure terms |
| Free-stream | `freestream/` | 1D uniform equilibrium; full PIC + hybrid (Hall off) |
| Light wave | `lightwave/` | 3D vacuum transverse wave on a periodic AMR grid |
| PCAI | `pcai/` | Proton-cyclotron anisotropy instability (hybrid benchmark) |
| Reconnection | `reconnection/` | Fadeev current sheet, uniform grid |
| AMR reconnection | `reconnection_amr/` | Fadeev current sheet, two-level AMR grid |
| Oblique shock | `shock/` | Magnetized shock; wall/inflow boundary study |
| Single-cell Hall | `singlecell/` | 1-cell periodic; curl B = 0 ⇒ Hall term exactly 0 |
| Zero-current wave | `zerocurrent/` | No particles (J = 0); Hall term vanishes |
| Hyper-resistivity | `hyper_resistivity/` | 4th-order `-eta_h ∇²J` Ohm term (smoke test) |
| Photoionization | `photoionization/` | Chamberlain profile, photoionization |
| Electron impact | `electronimpact/` | Voronov-1997 impact rate |
| Charge exchange | `chargeexchange/` | Constant cross-section CX with solar-wind ions |
| Chemistry | `chemistry/` | Chemical loss / source network |
| Recombination | `recombination/` | Recombination loss |
| Reflecting BC | `bc_reflecting/` | Specular reflecting particle wall |
| Absorbing BC | `bc_absorb/` | Grouped: EM-pulse absorber + particle absorber |
| Inflow BC | `bc_inflow/` | `#INFLOW`-driven inflow face |
| Wave injection BC | `bc_wave/` | Grouped: mono Bz wave + shear Alfvén wave via `#WAVEBC` |
| Performance | `performance/` | Beam-based scaling benchmark (see `validate_performance.py`) |

- Directories may ship a `PARAM.in.hybrid` (or `PARAM.in.<suffix>`) variant; the
  runner executes the test once per variant and lists both in the summary
  (e.g. `BEAM` and `BEAM (HYBRID)`).
- Results are written to `tests/summary.md`.

### Build the test executable

```bash
./Config.pl -lev=2 -u=Exo     # Exo source (ionization tests) + 2 AMR levels
make -j4
./Config.pl -amrex2d -lev=2 -u=Exo   # true-2D AMReX build (optional)
make -j4
```

The runner auto-detects the AMReX dimension and skips tests that need the other
one (e.g. the true-3D `lightwave` under a 2D build). Use `-amrex3d` to switch
back.

### Run

```bash
python3 tests/validate_tests.py                 # all tests, serial
python3 tests/validate_tests.py --test=beam     # one test
python3 tests/validate_tests.py -n 2            # 2 MPI ranks (mpirun -n 2)
python3 tests/validate_tests.py --verbose       # per-check diagnostics
python3 tests/validate_performance.py           # scaling benchmark
```

With `-n 1` (or omitted) the executable is invoked directly as `./FLEKS.exe`;
with `-n N > 1` via `mpirun -n N ./FLEKS.exe`.

### Ionization commands

- **`#PHOTOIONIZATION`** — per-component rates at the planet surface [s⁻¹],
  diluted as `(rPlanet / r)²`.
- **`#ELECTRONIMPACT`** — Voronov 1997:
  `σv(T) = A (T/EI)^K / [X + (T/EI)] · exp(-EI/T)` [cm³/s].
- **`#CHARGEEXCHANGE`** — constant cross-section `σv(u) = σ_CX |u_i|`; the
  matrix `σ_CX(neutral, ion)` is summed over ion species.

Architecture note: these parameters are stored in `SourceInterface` and read by
`UserSource::read_param()` (`userfiles/ExoSource.h`); `Domain` routes the
commands to the source object, keeping `FluidInterface` free of
ionization-specific data.

## Adding a test case / initial condition

Initial conditions are plug-ins resolved by name from `#TESTCASE` through the
`ICRegistry` (`include/InitialCondition.h`, `src/ic/RegisterAll.cpp`). An
unknown `#TESTCASE` name aborts loudly listing the registered names, so a typo
can never silently fall back to a uniform plasma.

Each test directory ships a `PARAM.in`, a `README.md` and (for most tests) its
own `validate.py`; the runner loads only that module for the test it runs.

- **Wave / sinusoidal tests need no C++.** All transverse and sinusoidal seeds
  are one parameterized plug-in (`WaveIC`), registered under the aliases
  `lightwave`, `hybridwave`, `convectionwave`, `ionacousticwave` and a generic
  `waveic`. Add a `#TESTCASE waveic` plus a `#WAVEIC` block in the test's
  `PARAM.in`; sub-parameters (`seedE`, `seedB`, `seedWeight`, `oblique`, `dir`,
  `waveLength`, `guideField`, `velKick`, `frac`, `pert`, `waveMode`) are all
  optional via `read_optional`.
- **Non-wave tests** (beam, tophat, fadeev) keep dedicated plug-ins: subclass
  `InitialCondition` (override `read_param`, `set_fields`, and the per-particle
  `modify_particle_weight` / `modify_particle_velocity` hooks, plus `name()`),
  register it in `src/ic/RegisterAll.cpp`, and add the `.cpp` to `SRCS` in
  `src/Makefile`.

## SWMF coupled tests (GM-PC / MHD-AEPIC)

Run from the **SWMF root**:

```bash
make test16_3d_compile
make test16_3d          # or: make test16_3d_rundir / _run / _check
```

- `test16_2d` and `test16_3d` are the Cartesian GM-PC regression tests.
- Default MPI is `mpiexec -n 2`; override with `NP=4`. Bless new reference
  results with `make test16_3d_check BLESS=YES`.
- Empty/near-empty `.diff` = PASS; results land in `run_test/RESULTS/3d/PC/`,
  references in `output/test16/`.
- Planetary coupled setups (e.g. Mars) use CON-level commands in `PARAM.in`
  (`#PLANET`, `#ROTATION`, `#STOP`, …); standalone-only commands
  (`#TIMEACCURATE`, `#SAVERESTART`, `#IDEALAXES`, …) are rejected inside
  `BEGIN_COMP GM`.
