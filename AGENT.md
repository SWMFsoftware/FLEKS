# FLEKS — FLexible Exascale Kinetic Simulator

## Overview

FLEKS is a Particle-In-Cell (PIC) and particle tracker code built on top of the
[AMReX](https://amrex-codes.github.io/) framework. It serves as the **PC**
(Particle-in-Cell) and **PT** (Particle Tracker) components within the
[SWMF](https://github.com/SWMFsoftware) (Space Weather Modeling Framework),
coupling with BATS-R-US (the MHD solver) via the MHD-AEPIC algorithm. It can
also run standalone for kinetic plasma simulations.

**Primary use cases:**
- Implicit Particle-In-Cell simulations (semi-implicit θ-scheme with GMRES
  field solver)
- Test particle tracking in electromagnetic fields
- Coupled MHD-PIC simulations within SWMF (MHD-AEPIC)
- Outer-heliosphere pickup-ion simulations (OH-PT coupling)
- Solar Energetic Particle (SEP) transport

**Language:** C++ (C++14/17), with Fortran 90 interface wrappers for SWMF
coupling.

**License:** Apache 2.0 (University of Michigan)

---

## Repository Layout

```
FLEKS/
├── AGENT.md              ← You are here
├── include/              ← All header files (.h)
├── src/                  ← All implementation files (.cpp) and src Makefile
├── srcInterface/         ← SWMF coupling layer (Fortran wrappers + C++ interface)
├── documents/            ← Algorithm docs (LaTeX) and coding standards
├── tools/                ← Post-processing and conversion scripts (Python/bash)
├── bin/                  ← Built executables (created by build)
├── .agent/skills/        ← AI agent skills (build, add-source, etc.)
├── Config.pl             ← Perl configuration script (AMReX, test particles, etc.)
├── configure.py          ← Standalone configure (writes Makefile.def)
├── Makefile              ← Top-level Makefile
├── Makefile.def.FLEKS    ← Default Makefile definitions for standalone
├── PARAM.XML             ← Parameter command documentation (XML, ~1500 lines)
├── .clang-format         ← Clang-format configuration (Mozilla-based, 80-col)
└── .gitignore
```

### Key directories

| Directory        | Purpose                                                     |
|------------------|-------------------------------------------------------------|
| `include/`       | All `.h` header files — class declarations, inline methods  |
| `src/`           | All `.cpp` files — implementations, `main.cpp`, src Makefile|
| `srcInterface/`  | SWMF coupling: `PC_wrapper.f90`, `PT_wrapper.f90`, `FleksInterface.cpp` |
| `documents/`     | `Algorithm.tex` (math/physics), `Coding_standards.md`       |
| `tools/`         | `amrex2tec.py`, `amrex2vtk.sh`, `converter.py`, `generate_compile_commands.py` |

---

## Architecture

### Class Hierarchy

```
Domain                        — Top-level simulation manager
 ├── DomainGrid               — Grid information container
 ├── Pic : Grid : AmrCore     — PIC solver (fields + particle push on AMR grid)
 ├── ParticleTracker : Grid   — Test particle tracker
 ├── FluidInterface : Grid    — MHD/fluid state on grid (coupling data)
 ├── SourceInterface          — Source terms for coupling
 ├── OHInterface              — Outer-heliosphere coupling data
 └── TimeCtr                  — Time stepping, event control, plot scheduling
      ├── EventCtr            — Periodic event trigger (dn or dt based)
      └── PlotCtr             — Plot scheduling (combines EventCtr + PlotWriter)
```

### Core Classes

| Class              | Header              | Purpose                                                   |
|--------------------|----------------------|-----------------------------------------------------------|
| `Domain`           | `Domain.h`           | Top-level orchestrator: owns Pic, FluidInterface, TimeCtr |
| `Pic`              | `Pic.h`              | PIC solver: field solve, particle push, moments, I/O      |
| `Grid`             | `Grid.h`             | AMR grid management (inherits `AmrCore`)                  |
| `Particles`        | `Particles.h`        | Templated particle container (PIC or PT particles)        |
| `TestParticles`    | `TestParticles.h`    | Test particle class (trajectory recording, I/O)           |
| `FluidInterface`   | `FluidInterface.h`   | Fluid/MHD state variables on PIC grid                     |
| `LinearSolver`     | `LinearSolver.h`     | GMRES Krylov solver for implicit E-field                  |
| `TimeCtr`          | `TimeCtr.h`          | Time step management, CFL, event scheduling               |
| `PlotWriter`       | `PlotWriter.h`       | Output formatting (IDL, AMReX, HDF5, VTK, Tecplot)        |
| `DataContainer`    | `DataContainer.h`    | Data reading (IDL format, binary/ascii)                   |
| `GridUtility`      | `GridUtility.h`      | Discrete operators (curl, div, grad, averaging)           |
| `BC`               | `BC.h`               | Boundary conditions (periodic, coupled, outflow, vacuum)  |
| `Constants`        | `Constants.h`        | Physical/numerical constants, particle record sizes       |

### SWMF Interface Layer (`srcInterface/`)

| File                  | Purpose                                                    |
|-----------------------|------------------------------------------------------------|
| `PC_wrapper.f90`      | Fortran wrapper for PC component (called by SWMF coupler)  |
| `PT_wrapper.f90`      | Fortran wrapper for PT component                           |
| `FleksInterface.cpp`  | C++ entry points called by Fortran wrappers via `bind(C)`  |

The interface layer follows a pattern:
1. SWMF calls `PC_wrapper.f90` subroutines (e.g., `PC_run`, `PC_put_from_gm`)
2. Fortran wrapper calls C functions in `FleksInterface.cpp` (e.g., `fleks_run_`)
3. C++ functions operate on the global `fleksDomains` object

### Algorithm Summary

- **PIC Method:** Semi-implicit θ-scheme with Boris particle mover
- **Field Solver:** GMRES iterative solver for the implicit electric field
- **Time Stepping:** CFL-based or fixed Δt; configurable θ parameter (default 0.51)
- **Divergence Cleaning:** Accurate div(E) correction via particle position adjustment; optional hyperbolic div(B) cleaning
- **Particle Management:** Splitting, merging, fast merge algorithm with Lagrange multipliers
- **AMR:** Block-structured AMR via AMReX with configurable refinement levels and regions

---

## Build System

### Dependencies

| Dependency | Required | Notes                                    |
|------------|----------|------------------------------------------|
| AMReX      | Yes      | At `../../util/AMREX/InstallDir/`        |
| MPI        | Yes      | `mpicxx` must be in PATH                 |
| HDF5       | Optional | Parallel HDF5 for HDF5 output support    |
| Perl       | Yes      | For `Config.pl` and `share/Scripts/`      |

### Build Targets

```bash
# Standalone executable
make FLEKS -j8

# Library for SWMF integration
make LIB -j8

# Converter tool
make CONVERTER

# Regenerate compile_commands.json (for IDE support)
make compile_commands

# Clean
make clean       # Object files only
make distclean   # Full reset
```

### Configuration

```bash
# Standalone setup (clones share/ and util/ if missing)
./Config.pl -s                    # Show current settings
./Config.pl -tp=PBE              # Set test particle output (P, PB, PBE, PBEG)
./Config.pl -lev=2               # Set max AMR levels
```

### Build Artifacts

| File                    | Description                         |
|-------------------------|-------------------------------------|
| `src/libFLEKS.a`       | Static library                      |
| `bin/FLEKS.exe`         | Standalone executable               |
| `bin/converter.exe`     | Data format converter               |
| `compile_commands.json` | Compilation database for IDE        |

---

## Coding Conventions

See `documents/Coding_standards.md` for the full list. Key points:

### Naming

| Element        | Convention       | Example                    |
|----------------|------------------|----------------------------|
| Files          | PascalCase       | `GridUtility.cpp`          |
| Classes        | PascalCase       | `FluidInterface`           |
| Functions      | snake_case       | `apply_float_boundary()`   |
| Variables      | camelCase         | `nCellPerPatch`            |
| Constants      | camelCase/UPPER  | `cProtonMassSI`, `nDim3`   |
| Private members| camelCase         | `doRestart`                |

### Style Rules

1. **Smart pointers** — Use `unique_ptr`/`shared_ptr` for ownership; raw
   pointers only when there is no ownership
2. **`using namespace amrex`** — Allowed only in `.cpp` files, never in headers
3. **Header order** — std → AMReX → project headers
4. **`nullptr`** over `NULL`
5. **`const`** wherever possible
6. **80-column limit** — Enforced by `.clang-format` (Mozilla-based style)
7. **Lambdas** — OK for short/local use; prefer named functions for reusable logic
8. **Commits** — Follow [conventional commits](https://www.conventionalcommits.org/)
9. **Include guards** — `#ifndef _FILENAME_H_` / `#define _FILENAME_H_` pattern

### Formatting

**C++ Formatting:**
The project uses Clang-Format with a Mozilla-based configuration. Key settings:
- `IndentWidth: 2`
- `ColumnLimit: 80`
- `BreakBeforeBraces: Attach`
- `UseTab: Never`

**Fortran Formatting:**
For `.F90`/`.f90` files, use `findent` configured to match Emacs' `f90-mode` indentation. This is important to remember when modifying Fortran source files in `srcInterface/` or other components.

---

## Parameter System

Simulation parameters are controlled via `PARAM.in` files using the SWMF
parameter reading system. All supported commands are documented in `PARAM.XML`.

### Key Command Groups

| Group                        | Example Commands                                   |
|------------------------------|----------------------------------------------------|
| **Output**                   | `#SAVEPLOT`, `#MONITOR`, `#SAVELOG`, `#NOUTFILE`   |
| **Scheme**                   | `#PIC`, `#TIMESTEPPING`, `#DISCRETIZATION`, `#EFIELDSOLVER`, `#DIVE`, `#DIVB` |
| **Particles**                | `#PARTICLES`, `#RESAMPLING`, `#FASTMERGE`, `#VACUUM`, `#PARTICLETRACKER` |
| **Initial/Boundary Cond.**   | `#GEOMETRY`, `#NCELL`, `#REGION`, `#PARTICLES`, `#BC` |
| **Coupling**                 | `#OHION`, `#CHARGEEXCHANGE`, `#MAXCHARGEEXCHANGERATE` |

---

## Testing & Running

When coupled with SWMF, FLEKS is initialized and driven by the SWMF framework.
For standalone use, the `main.cpp` sets up AMReX and a `Domain` object.

### Run Directory Structure

```
run/PC/
├── restartIN/
├── restartOUT/
└── plots/
```

### Tools & Post-Processing

| Script                      | Purpose                                |
|-----------------------------|----------------------------------------|
| `tools/amrex2tec.py`       | Convert AMReX output to Tecplot format |
| `tools/amrex2vtk.sh`       | Convert AMReX output to VTK format     |
| `tools/converter.py`       | General data conversion                |
| `tools/generate_compile_commands.py` | Generate `compile_commands.json` |

**Visualization:** While ParaView can read HDF5 data natively, it may struggle with continuous rendering of AMReX block boundaries unless you use the `BATSRUSReader` ParaView plugin provided by the `flekspy` Python toolkit. Alternatively, converting output to `.vtm` or Tecplot formats using the tools above is often preferred. For Python scripting, the official tool is the **`flekspy`** package (`pip install flekspy`), which natively parses FLEKS AMReX data and integrates with Matplotlib and YT.

---

## AI Agent Skills

To assist with common tasks, specialized instructions (skills) are defined in `.agent/skills/`. You are encouraged to view these Markdown files when performing related tasks:

| Skill             | Description                                          | Path                                           |
|-------------------|------------------------------------------------------|------------------------------------------------|
| **Build FLEKS**   | Compile project and regenerate compile_commands.json | `.agent/skills/build-fleks/SKILL.md`           |
| **Add New Source**| Create new `.cpp`/`.h` files following conventions   | `.agent/skills/add-new-source/SKILL.md`        |
| **Code Cleanup**  | Formatting, unused variables, standard checks        | `.agent/skills/code-cleanup/SKILL.md`          |
| **Debug Session** | Setup and launch a debug session (lldb / VS Code)    | `.agent/skills/debug-session/SKILL.md`         |
| **Generate Docs** | Build HTML/PDF documentation                         | `.agent/skills/generate-docs/SKILL.md`         |

---

## Common Workflows

### Adding a New Source File

1. Create `include/NewFeature.h` with include guards
2. Create `src/NewFeature.cpp`
3. Add `NewFeature.cpp` to `SRCS` in `src/Makefile`
4. Build: `make LIB -j8`

### Adding a New Parameter Command

1. Document in `PARAM.XML`
2. Parse in the relevant `read_param()` method (usually in `Domain.cpp` or `Pic.cpp`)
3. Add the corresponding member variable

### Modifying the SWMF Interface

1. Add C++ function in `srcInterface/FleksInterface.cpp`
2. Add Fortran wrapper in `PC_wrapper.f90` or `PT_wrapper.f90`
3. Ensure C-compatible calling convention (`extern "C"` or `bind(C)`)

---

## Important Notes

- FLEKS is always 3D internally. 2D is achieved with a single cell in the
  z-direction ("fake 2D").
- The `Constants.h` file is generated from `Constants.h.orig` and should not be
  edited directly — use `Config.pl` to change settings like `ptRecordSize` or
  `nLevMax`.
- `UserSource.h` is similarly copied from `UserSource.h.orig` and is intended
  for user-customizable source terms.
- The `show_git_info.h` file is auto-generated at build time with git revision
  information.
