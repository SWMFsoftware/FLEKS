---
name: fleks-expert
description: Architecture, build, test, coupling and coding guidance for FLEKS (AMReX-based PIC / particle tracker; the PC and PT components of SWMF). Use for any FLEKS task — implementing parameters, sources or boundary conditions, building (make LIB / EXE), running standalone tests or GM-PC coupled tests, editing the SWMF interface, or debugging a PIC run.
---

# FLEKS Expert

FLEKS is a Particle-In-Cell and test-particle code built on AMReX. It runs as
the **PC** (PIC) and **PT** (particle tracker) components inside SWMF
(MHD-AEPIC, coupled to BATS-R-US/GM) and also standalone via `bin/FLEKS.exe`.

Two field solvers, selected per run (`PARAM.in` vs a `PARAM.in.<suffix>`
variant):

- **Full PIC** — kinetic ions *and* electrons, semi-implicit θ-scheme
  (default θ = 0.51), implicit E solved by GMRES.
- **Hybrid PIC** — kinetic ions + massless fluid electrons, generalized Ohm's
  law, explicit RK4/SSPRK3 Faraday advance (`#HYBRIDPIC`).

## Single source of truth

| Content | Location |
|---|---|
| Human overview / quick start | `README.md` |
| Build, layout, parameters, tools, extension recipes | `.agent/skills/fleks-expert/references/` |
| Algorithms and math | `doc/Algorithm.tex` |
| Coding standards | `doc/Coding_standards.md` |
| Parameter reference | `PARAM.XML` (`make PDF` → `doc/USERMANUAL.pdf`) |
| Test catalogue & runner options | `tests/README.md` |
| Task recipes | `.agent/workflows/`, `.agent/skills/*/SKILL.md` |

The root `AGENT.md` is only a **router** to these files; there are no
per-directory `AGENT.md` files. Extend the files listed above instead of adding
new ones.

## Commands

```bash
./Config.pl -lev=2 -u=Exo       # configure a standalone build
make EXE -j8                    # standalone bin/FLEKS.exe (+ converter, compile_commands)
make LIB -j8                    # SWMF component library src/libFLEKS.a
python3 tests/validate_tests.py [--test=beam] [-n 2] [--verbose]
make test16_3d                  # GM-PC regression — run from the SWMF root
python3 tools/format_all.py     # CI-enforced formatting (C++ + Fortran)
```

## Hard constraints

Violating any of these gives silently wrong results or link errors:

1. `Particles<NStructReal,NStructInt>` is split across `src/Particles*.cpp`,
   each compiled independently. Every new out-of-line method needs **explicit
   instantiations for both `PicParticles` and `PTParticles`** in the file that
   defines it; whole-class instantiation belongs only in `src/Particles.cpp`
   (repeating it elsewhere gives duplicate symbols).
2. `include/Constants.h` and `include/UserSource.h` are **generated** — from
   `Constants.h.orig` (via `Config.pl -tp=`, `-lev=`) and from
   `userfiles/*Source.h` (via `Config.pl -u=...`). Edit the sources, never the
   generated files.
3. FLEKS is always 3D internally; 2D is one cell in z ("fake 2D"). A true-2D
   AMReX build requires `./Config.pl -amrex2d`.
4. Builds define `_PC_COMPONENT_` / `_PT_COMPONENT_`; guarded code must stay
   valid for both.
5. Standalone runs use domain name `FLEKS1` (plots/restarts under `FLEKS1/`)
   and require `#INITFROMSWMF F` plus `#NORMALIZATION` / `#PLASMA` /
   `#UNIFORMSTATE`.
6. New `.cpp` files must be listed in `SRCS` in `src/Makefile` — nothing is
   auto-discovered.
7. Formatting is CI-enforced: `python3 tools/format_all.py`
   (`clang-format==20.1.8` + `findent`).
8. Ionization parameters live in `SourceInterface` / `UserSource`, not in
   `FluidInterface` (which is reserved for MHD coupling data).

## Task routing

| Task | Go to |
|---|---|
| Multi-file change or long build/test loop | spawn the `fleks` agent (`.agent/agents/fleks.md`) |
| Build / compile / link errors | `.agent/skills/build-fleks/SKILL.md` |
| Build configuration, targets, Makefile | `references/build.md` |
| Parameter groups / command lookup | `references/parameters.md` |
| Post-processing / conversion tools | `references/tools.md` |
| New `.cpp`/`.h` or user source | `.agent/skills/add-new-source/SKILL.md` |
| Formatting, cleanup, conventions | `.agent/skills/code-cleanup/SKILL.md` |
| Debugger session (gdb/lldb) | `.agent/skills/debug-session/SKILL.md` |
| LaTeX / Doxygen / PARAM.XML docs | `.agent/skills/generate-docs/SKILL.md` |
| Add a new `#COMMAND` | `.agent/workflows/add-param.md` |
| Add a GM↔PC coupling variable | `.agent/workflows/add-coupling-var.md` |
| Run the GM-PC regression | `.agent/workflows/run-test.md` |
| Standalone test catalogue | `references/testing.md`, `tests/README.md` |
| Which class / file to edit | `references/file-layout.md` |
| Solver & physics details | `references/architecture.md` |
| SWMF Fortran↔C++ interface | `references/coupling.md` |
| Naming / style checklist | `references/standards.md` |

## References

Load only what the current task needs. Each page names the canonical source it
summarizes (`PARAM.XML`, `doc/Coding_standards.md`, `tests/README.md`) — open
that source when the task needs the full detail:

- `references/architecture.md` — class hierarchy, full vs hybrid solver, time
  stepping, divergence cleaning.
- `references/file-layout.md` — repository layout and where to look for what.
- `references/build.md` — `Config.pl` options, make targets, dependencies,
  `src/Makefile`, standalone runs.
- `references/parameters.md` — command groups and `PARAM.XML` conventions.
- `references/tools.md` — scripts and reading FLEKS output.
- `references/standards.md` — naming, style, SWMF interface patterns.
- `references/testing.md` — standalone and coupled test suites, adding a test
  case or initial condition.
- `references/coupling.md` — `srcInterface/` layer, entry-point tables,
  Fortran/C++ interoperability pitfalls.

## Output data

FLEKS writes AMReX block-structured output (plus IDL/`.h` for SWMF
post-processing). Do not use generic NetCDF loaders — use `flekspy`
(`pip install flekspy`) in Python or `Batsrus.jl` in Julia.
