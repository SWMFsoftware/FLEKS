# FLEKS — FLexible Exascale Kinetic Simulator

FLEKS is a Particle-In-Cell and test-particle code built on AMReX. It runs as
the **PC** (PIC) and **PT** (particle tracker) components inside SWMF
(MHD-AEPIC, coupled to BATS-R-US/GM) and standalone via `bin/FLEKS.exe`.

Two field solvers, selected per run:

- **Full PIC** — kinetic ions and electrons, semi-implicit θ-scheme, GMRES
  E-field solve.
- **Hybrid PIC** — kinetic ions + massless fluid electrons, generalized Ohm's
  law, explicit Faraday advance (`#HYBRIDPIC`).

## Single source of truth

| Content | Location |
|---|---|
| Human overview / quick start | `README.md` |
| Build, layout, parameters, tools, extension recipes | `.agents/skills/fleks-expert/references/` |
| Coding standards | `doc/Coding_standards.md` |
| Algorithms and math | `doc/Algorithm.tex` |
| Parameter reference | `PARAM.XML` (`make PDF` → `doc/USERMANUAL.pdf`) |
| Standalone test catalogue & runner | `tests/README.md` |
| Agent knowledge base | `.agents/skills/fleks-expert/` (+ `references/`) |
| Task recipes | `.agents/skills/*/SKILL.md` |
| Multi-step guides | `.agents/workflows/` |
| Submission rules (formatting, commits) | `CONTRIBUTING.md` |

## Commands

```bash
./Config.pl -lev=2 -u=Exo       # configure a standalone build
make EXE -j8                    # standalone bin/FLEKS.exe (+ Converter/)
make LIB -j8                    # SWMF component library src/libFLEKS.a
make CONVERTER                  # only bin/converter.exe (make -C Converter)
python3 tests/validate_tests.py [--test=beam] [-n 2] [--verbose]
make test16_3d                  # GM-PC regression — from the SWMF root
python3 tools/format_all.py     # CI-enforced formatting (C++ + Fortran)
```

## Hard constraints

Violating any of these gives silently wrong results or link errors:

1. `Particles` is split across independently compiled `src/Particles*.cpp`
   files. Every new out-of-line method needs **explicit instantiations for both
   `PicParticles` and `PTParticles`** in the defining file; whole-class
   instantiation belongs only in `src/Particles.cpp`.
2. `include/Constants.h` and `include/UserSource.h` are **generated** (from
   `Constants.h.orig` and `userfiles/*Source.h`). Edit the sources.
   Sync `include/Constants.h` via `Config.pl` or `make install`. Sync `include/UserSource.h`
   via `./Config.pl -u=<Name>`. Custom source commands must be registered via
   `register_parameter_commands()` and documented in `PARAM.XML`.
3. FLEKS is always 3D internally; 2D is one cell in z ("fake 2D"). A true-2D
   AMReX build needs `./Config.pl -amrex2d`.
4. Builds define `_PC_COMPONENT_` / `_PT_COMPONENT_`; guarded code must stay
   valid for both.
5. Standalone runs use domain name `FLEKS1` and require `#INITFROMSWMF F` plus
   `#NORMALIZATION` / `#PLASMA` / `#UNIFORMSTATE`.
6. `src/Makefile` picks up every `.cpp` in `src/` and `src/ic/` with a
   wildcard; only `src/main.cpp` is excluded (it has its own `main()`). A new
   FLEKS source needs no Makefile edit — but it will be compiled, so do not
   leave scratch `.cpp` files in `src/`.
   The format converter is **not** part of the library: it lives in
   `Converter/` with its own `Makefile` and its own `main()`.
7. Formatting is CI-enforced (`python3 tools/format_all.py`).
8. Ionization parameters belong in `SourceInterface` / `UserSource`, not in
   `FluidInterface`.

## Agent skills and workflows

`.agents/` is versioned with the code. Antigravity discovers `skills/` there
natively; other tools need a symlink into their own skills directory — see
`.agents/README.md`.

| Entry point | Use it for |
|---|---|
| `agents/fleks.md` | Dedicated agent for context-heavy work (multi-file changes, long build/test loops) — it reads the docs for you and reports a short summary |
| `skills/fleks-expert/` | Knowledge base — architecture, file layout, standards, testing, coupling (references loaded on demand) |
| `skills/build-fleks/` | Compiling and `compile_commands.json` |
| `skills/add-new-source/` | New `.cpp`/`.h` or user-source templates |
| `skills/code-cleanup/` | Formatting, unused variables, conventions |
| `skills/debug-session/` | gdb / lldb sessions |
| `skills/generate-docs/` | LaTeX, Doxygen, PARAM.XML docs |
| `workflows/add-param.md` | Add a new `#COMMAND` end-to-end |
| `workflows/add-coupling-var.md` | Add a GM↔PC exchange variable |
| `workflows/run-test.md` | Run the `test16_3d` GM-PC regression |

## Output data

FLEKS writes AMReX block-structured output (plus IDL/`.h` for SWMF
post-processing). Do not use generic NetCDF loaders — use `flekspy`
(`pip install flekspy`) in Python or `Batsrus.jl` in Julia.
