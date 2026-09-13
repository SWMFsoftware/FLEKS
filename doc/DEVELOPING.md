# Developing FLEKS

This is the developer-oriented companion to the [README](../README.md): where
things live, how the build works, and how to extend the code. For submission
rules (formatting, commit style, PR checklist) see
[CONTRIBUTING.md](../CONTRIBUTING.md) and
[Coding_standards.md](Coding_standards.md). Mathematical derivations live in
[Algorithm.tex](Algorithm.tex).

FLEKS is a Particle-In-Cell and test-particle code on top of AMReX. It serves
as the **PC** (PIC) and **PT** (particle tracker) components inside SWMF
(MHD-AEPIC, coupled to BATS-R-US/GM) and also builds as a standalone
executable. Two field solvers are selected per run:

- **Full PIC** — kinetic ions and electrons, semi-implicit θ-scheme
  (default θ = 0.51), implicit E solved by GMRES.
- **Hybrid PIC** — kinetic ions + massless fluid electrons, generalized Ohm's
  law, explicit RK4/SSPRK3 Faraday advance (`#HYBRIDPIC`).

---

## 1. Repository Layout

```
FLEKS/
├── include/           Public headers (.h)
├── src/               C++ sources (.cpp), src/Makefile (SRCS list), main.cpp
│   └── ic/            Initial-condition plug-ins (private headers live here)
├── srcInterface/      SWMF coupling layer: PC_wrapper.f90, PT_wrapper.f90,
│                      FleksInterface.cpp
├── doc/               Algorithm.tex, Coding_standards.md, DEVELOPING.md, Tex/
├── tools/             Post-processing, conversion and formatting scripts
├── tests/             Standalone test suite (one directory per scenario)
├── userfiles/         Selectable user-source templates (*Source.h)
├── .agent/            Agent skills, workflows and knowledge base
├── Config.pl          Perl configuration (AMReX, AMR levels, user source, …)
├── Makefile           Top-level makefile
├── Makefile.def.FLEKS Default FLEKS makefile definitions
├── PARAM.XML          Parameter command reference (source of the user manual)
└── .clang-format      Mozilla-based style, 80 columns, 2-space indent
```

## 2. Build System

### Configuration

```bash
./Config.pl                  # show current settings
./Config.pl -install         # generate Makefile.conf if missing
./Config.pl -lev=2           # max AMR levels
./Config.pl -u=Exo           # select userfiles/ExoSource.h -> include/UserSource.h
./Config.pl -u               # list available user sources
./Config.pl -tp=PBE          # test-particle output (P, PB, PBE, PBEG)
./Config.pl -amrex2d         # true-2D AMReX build (-amrex3d to switch back)
```

`Config.pl` locates `share`, `util`, `lib` and AMReX. Inside an SWMF tree
(`SWMF/PC/FLEKS`) a bare `./Config.pl -install` reuses the parent SWMF paths; in
a pure FLEKS checkout it can clone or use local dependencies.

### Targets

| Target | Result |
|---|---|
| `make EXE -j8` (alias `make FLEKS`) | standalone `bin/FLEKS.exe` |
| `make LIB -j8` | SWMF component library `src/libFLEKS.a` + wrappers |
| `make CONVERTER` | `bin/converter.exe` |
| `make compile_commands` | regenerate `compile_commands.json` (IDE index) |
| `make clean` / `make distclean` | object files / full reset |
| `make PDF` | `doc/USERMANUAL.pdf` (needs `pdflatex`, `makeindex`, `fvextra`) |

`make LIB` is only valid inside a built SWMF tree (needs `libSHARE.a` and
`con_comp_param.mod`); use `make EXE` for standalone work. `EXE`, `FLEKS` and
`LIB` all invoke `compile_commands` and `CONVERTER` automatically.

### Dependencies

| Dependency | Required | Notes |
|---|---|---|
| AMReX | yes | `../../util/AMREX/InstallDir/` in SWMF, `util/AMREX/InstallDir/` standalone |
| MPI | yes | `mpicxx`/`mpif90` on PATH |
| SWMF share | yes | `share/Scripts`, `share/Library`, `libSHARE.a` |
| HDF5 | optional | parallel HDF5 output |
| Perl | yes | `Config.pl` and `share/Scripts/` |

### src/Makefile

Every `.cpp` must be listed in the `SRCS` variable — nothing is
auto-discovered. Useful variables:

| Variable | Purpose |
|---|---|
| `SRCS` | sources compiled into `libFLEKS.a` |
| `SEARCH_C` | include search paths (e.g. `-I../include`) |
| `FLAGC_EXTRA` | `-D_${COMPONENT}_COMPONENT_` |
| `LIBFLEKS` | output library name |

Preprocessor flags: `_PC_COMPONENT_` (PIC build), `_PT_COMPONENT_` (particle
tracker), `_USE_HDF5_` (HDF5 output).

## 3. Code Layout and Entry Points

Do not rely on exhaustive file catalogues — they go stale. Use `ls` plus
semantic navigation; the table below lists the stable entry points.

| Area | Entry point | Neighbours |
|---|---|---|
| Time loop, regridding, orchestration | `src/Pic.cpp` (`Pic::update`) | `src/Domain.cpp` |
| Parameter parsing | `src/PicParam.cpp` (`Pic::read_param`) | `src/Domain.cpp` |
| Field boundary conditions | `src/PicBC.cpp` | `include/BC.h`, `src/BC.cpp` |
| Full-PIC field solve | `src/PicFieldSolver.cpp` | `src/PicDivE.cpp` (div E cleaning) |
| Hybrid solver | `src/PicHybrid.cpp` | `include/Pic.h` |
| Particles | `src/Particles.cpp` | `src/ParticlesInit.cpp`, `src/ParticlesBC.cpp`, `src/ParticlesMoments.cpp`, `src/ParticlesMassMatrix.cpp`, `src/ParticlesMover.cpp`, `src/ParticlesResample.cpp`, `src/ParticlesReactions.cpp` |
| Fluid / coupling state | `src/FluidInterface.cpp` | `include/FluidInterface.h` |
| Output | `src/PlotWriter.cpp` | `src/PicIO.cpp`, `src/DataContainer.cpp` |
| Linear algebra | `src/LinearSolver.cpp` | `include/LinearSolver.h` |
| Grid, AMR, load balancing | `src/Grid.cpp` | `src/FleksDistributionMap.cpp` |
| Standalone driver | `src/main.cpp` | — |
| SWMF entry points | `srcInterface/FleksInterface.cpp` | `PC_wrapper.f90`, `PT_wrapper.f90` |

### Class hierarchy

```
Domain                        — top-level simulation manager
 ├── DomainGrid               — grid information container
 ├── Pic : Grid : AmrCore     — PIC solver (fields + particle push on AMR grid)
 ├── ParticleTracker : Grid   — test particle tracker
 ├── FluidInterface : Grid    — MHD/fluid state on grid (coupling data)
 ├── SourceInterface          — source terms
 │    └── UserSource          — selected user source (userfiles/*Source.h)
 ├── OHInterface              — outer-heliosphere coupling data
 └── TimeCtr                  — time stepping, event control, plot scheduling
      ├── EventCtr            — periodic event trigger (dn or dt based)
      └── PlotCtr             — plot scheduling (EventCtr + PlotWriter)
```

Ionization parameters live in **SourceInterface** (physics in `UserSource`),
never in `FluidInterface`, which stays reserved for MHD coupling data.

### Non-obvious layout rules

1. **`src/ic/` keeps its own headers.** `WaveIC.h`, `BeamIC.h`, `TopHatIC.h`
   are private plug-in headers: the public surface is
   `include/InitialCondition.h` (abstract base class + `ICRegistry`). If an IC
   header is ever needed outside `src/ic/`, promote it to `include/` and
   regenerate dependencies (`make DEPEND`).
2. **Generated files — never edit the output.**
   `include/Constants.h` ← `Constants.h.orig` (via `Config.pl -tp=`, `-lev=`);
   `include/UserSource.h` ← `userfiles/<Name>Source.h` (via `Config.pl -u=`);
   `include/show_git_info.h` is created at build time.
3. **Split `Particles` translation units.** `Particles<NStructReal,NStructInt>`
   is implemented in independently compiled `src/Particles*.cpp` files and
   instantiated only for `PicParticles` and `PTParticles`. Every new
   out-of-line method needs explicit member-function instantiations for both
   aliases **in the file that defines it**; whole-class instantiation belongs
   only in `src/Particles.cpp` (repeating it elsewhere gives duplicate symbols
   with the Mach-O linker). Keep the two aliases synchronized when a signature
   changes.
4. **Fake 2D** = a single cell in z; FLEKS is always 3D internally unless AMReX
   itself is built 2D (`./Config.pl -amrex2d`).
5. **Component macros.** Code guarded by `_PC_COMPONENT_` / `_PT_COMPONENT_`
   must stay valid for both builds.

### Coding conventions

Full text in [Coding_standards.md](Coding_standards.md); essentials:
`PascalCase` files/classes, `camelCase` variables/members, `snake_case`
functions, `#ifndef _FILENAME_H_` guards, header order std → AMReX → project,
no `using namespace` in headers, `nullptr`, `const` everywhere, smart pointers
for ownership, 80-column limit.

Formatting is CI-enforced:

```bash
pip install clang-format==20.1.8 findent
python3 tools/format_all.py
```

## 4. Extending FLEKS

### Add a source file

1. Create `include/NewFeature.h` (include guards, no `using namespace`).
2. Create `src/NewFeature.cpp`.
3. Add `NewFeature.cpp` to `SRCS` in `src/Makefile`.
4. `make LIB -j8` (and `make EXE -j8` if the standalone path uses it).

### Add a parameter command

1. Document it in `PARAM.XML` (`<command>` block; use `multiple="T"` and
   `COMMANDNAME_FLEKS{0,1,2}` aliases for per-domain overrides).
2. Add the member variable to the owning class header
   (`Pic.h`, `FluidInterface.h`, `Domain.h`, `TimeCtr.h`, …).
3. Parse it in the matching `read_param()` (`src/PicParam.cpp`,
   `src/Domain.cpp`, …).
4. Run `make test16_3d` from the SWMF root to check nothing breaks.

Detail: `.agent/workflows/add-param.md`.

### Add a test case / initial condition

Initial conditions are plug-ins resolved by name from `#TESTCASE` through the
`ICRegistry` (`include/InitialCondition.h`, `src/ic/RegisterAll.cpp`). An
unknown `#TESTCASE` name aborts loudly listing the registered names, so a typo
can never silently fall back to a uniform plasma.

- **Wave / sinusoidal tests need no C++.** All transverse and sinusoidal seeds
  are one parameterized plug-in (`WaveIC`), registered under the aliases
  `lightwave`, `hybridwave`, `convectionwave`, `ionacousticwave` and a generic
  `waveic`. Add a `#TESTCASE waveic` plus a `#WAVEIC` block in the test's
  `PARAM.in`; sub-parameters (`seedE`, `seedB`, `seedWeight`, `oblique`, `dir`,
  `waveLength`, `guideField`, `velKick`, `frac`, `pert`, `waveMode`) are all
  optional via `read_optional`.
- **Non-wave tests** (beam, tophat) keep dedicated plug-ins: subclass
  `InitialCondition` (override `read_param`, `set_fields`, and the per-particle
  `modify_particle_weight` / `modify_particle_velocity` hooks, plus `name()`),
  register it in `src/ic/RegisterAll.cpp`, and add the `.cpp` to `SRCS`.

### Add a user source

Add `userfiles/NewSource.h` (name must end in `Source.h`; the selection name is
the prefix), override `set_source()` when it depends on fluid fields, then
select it with `./Config.pl -u=New`, which copies it to
`include/UserSource.h`. Enable with `#SOURCE` in `PARAM.in`.

### Modify the SWMF interface

1. Add the C++ function in `srcInterface/FleksInterface.cpp` with
   `extern "C"` linkage (trailing underscore for Fortran interop).
2. Declare it in `include/FleksInterface.h`.
3. Add the interface block and wrapper subroutine in `PC_wrapper.f90` /
   `PT_wrapper.f90` using `iso_c_binding`.

Coupling layer details, entry-point tables and Fortran/C++ pitfalls:
`.agent/skills/fleks-expert/references/coupling.md`. Adding an exchanged GM↔PC
variable: `.agent/workflows/add-coupling-var.md`.

### Modify standalone behaviour

Update `src/main.cpp` for driver-level changes, keep `Domain` logic shared with
SWMF whenever possible, rebuild with `make EXE -j8`, and run from a directory
containing `PARAM.in`. Standalone runs use domain name `FLEKS1` (plots and
restarts under `FLEKS1/`) and need `#INITFROMSWMF F` plus `#NORMALIZATION`,
`#PLASMA` and `#UNIFORMSTATE`.

## 5. Parameters

All input commands are documented in `PARAM.XML`; `make PDF` renders
`doc/USERMANUAL.pdf`. Major groups:

| Group | Commands |
|---|---|
| Output | `#SAVEPLOT`, `#MONITOR`, `#SAVELOG`, `#NOUTFILE` |
| Scheme | `#PIC`, `#TIMESTEPPING`, `#DISCRETIZATION`, `#EFIELDSOLVER`, `#DIVE`, `#DIVB` |
| Hybrid PIC | `#HYBRIDPIC`, `#RESISTIVITY`, `#HALLTERM`, `#ELECTRONTEMPERATURE`, `#HYPERRESISTIVITY`, `#BSUBCYCLE`, `#MINIMUMDENSITY`, `#FIELDINTEGRATOR`, `#AVGFIELDB`, `#SMOOTHMOMENTS` |
| Particles | `#PARTICLES`, `#RESAMPLING`, `#FASTMERGE`, `#VACUUM`, `#PARTICLETRACKER` |
| Initial / boundary | `#GEOMETRY`, `#NCELL`, `#REGION`, `#BC` |
| Coupling | `#OHION`, `#CHARGEEXCHANGE`, `#MAXCHARGEEXCHANGERATE` |

Do not use `#` for referencing a command inside a `PARAM.in` comment — it would
be parsed as a command.

## 6. Testing

Two independent systems:

- **Standalone** (`tests/`, no SWMF) — the catalogue and runner options are in
  [tests/README.md](../tests/README.md); run with
  `python3 tests/validate_tests.py [--test=NAME] [-n N] [--verbose]`. Build the
  executable first with `./Config.pl -lev=2 -u=Exo && make -j4`.
- **SWMF coupled** (GM-PC / MHD-AEPIC) — `make test16_3d` from the SWMF root.

Use standalone tests for solver-level changes (fields, particles, BCs,
sources); use coupled tests for interface changes (`FleksInterface.cpp`,
`PC_wrapper.f90`, `GM_couple_pc.f90`).

## 7. Tools and Post-Processing

| Script | Purpose |
|---|---|
| `tools/amrex2tec.py` | AMReX plot → Tecplot |
| `tools/amrex2vtk.sh` | AMReX plot → VTK (ParaView) |
| `tools/tec2vtk.sh` | Tecplot → VTK |
| `tools/converter.py` | general data conversion |
| `tools/clean_dat.py` | clean up `.dat` output |
| `tools/format_all.py` | bulk reformat (C++ + Fortran), required before a PR |
| `tools/generate_compile_commands.py` | regenerate `compile_commands.json` |
| `tools/install_skill.sh` | install `.agent/skills/*` into CodeBuddy |
| `tools/check_docs.py` | verify the agent documentation tree (single `AGENT.md`, skill frontmatter, resolvable paths, canonical docs) — runs in CI |

Output is AMReX block-structured (plus IDL/`.h` for SWMF post-processing); do
not use generic NetCDF loaders. Use `flekspy` (`pip install flekspy`) in Python
or `Batsrus.jl` in Julia. ParaView can read HDF5 directly but handles AMReX
block boundaries better with the BATSRUSReader plugin from `flekspy`.

## 8. Documentation

| File | Build |
|---|---|
| `doc/Algorithm.tex` | `cd doc && pdflatex Algorithm.tex` (twice for references) |
| `PARAM.XML` | `make PDF` → `doc/USERMANUAL.pdf` |
| `doc/DEVELOPING.md` | this file |

Keep `doc/Coding_standards.md` and this file in sync with structural changes:
new directories, new generated files, new extension points.
