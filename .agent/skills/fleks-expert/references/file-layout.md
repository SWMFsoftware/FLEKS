# FLEKS File Layout

## Directories

| Directory | Contents |
|---|---|
| `include/` | All public headers (`.h`). |
| `src/` | All C++ sources (`.cpp`) and `src/Makefile` (`SRCS` list). |
| `src/ic/` | Initial-condition plug-ins; **private** IC headers live here. |
| `srcInterface/` | SWMF coupling layer: `PC_wrapper.f90`, `PT_wrapper.f90`, `FleksInterface.cpp`. |
| `Doc/` | `Algorithm.tex`, `Coding_standards.md`, user manual build. |
| `tools/` | Post-processing / conversion / formatting scripts. |
| `tests/` | Standalone test suite (one directory per scenario). |
| `userfiles/` | Selectable user-source templates (`*Source.h`). |
| `.agent/` | Agent skills, workflows and this knowledge base. |
| `Config.pl` | Perl configuration (AMReX, AMR levels, user source, test particles). |
| `PARAM.XML` | Parameter command reference (source of `Doc/USERMANUAL.pdf`). |

## Where to look for what

Do not rely on exhaustive file catalogues — they go stale. Use `ls` plus
semantic lookup (`documentSymbol` / `workspaceSymbol`); the entries below are
the stable entry points.

| Area | Entry point | Neighbours |
|---|---|---|
| Time loop, regridding, orchestration | `src/Pic.cpp` (`Pic::update`) | `src/Domain.cpp` |
| Parameter parsing | `src/PicParam.cpp` (`Pic::read_param`) | `src/Domain.cpp` (`Domain::read_param`) |
| Field boundary conditions | `src/PicBC.cpp` | `include/BC.h`, `src/BC.cpp` |
| Hybrid solver | `src/PicHybrid.cpp` | `src/PicFieldSolver.cpp` (full PIC) |
| div(E) cleaning | `src/PicDivE.cpp` | — |
| Particles | `src/Particles.cpp` (lifecycle + whole-class instantiation) | `src/ParticlesInit.cpp`, `ParticlesBC.cpp`, `ParticlesMoments.cpp`, `ParticlesMassMatrix.cpp`, `ParticlesMover.cpp`, `ParticlesResample.cpp`, `ParticlesReactions.cpp` |
| Fluid / coupling state | `src/FluidInterface.cpp` | `include/FluidInterface.h` |
| Output | `src/PlotWriter.cpp` | `src/DataContainer.cpp`, `src/PicIO.cpp` |
| Linear algebra | `src/LinearSolver.cpp` | `include/LinearSolver.h` |
| Grid / AMR | `src/Grid.cpp` | `src/FleksDistributionMap.cpp` |
| Standalone driver | `src/main.cpp` | — |
| SWMF entry points | `srcInterface/FleksInterface.cpp` | `PC_wrapper.f90`, `PT_wrapper.f90` |

## Non-obvious layout rules

1. **`src/ic/` keeps its own headers.** `WaveIC.h`, `BeamIC.h`, `TopHatIC.h`
   are private plug-in headers, intentionally outside `include/`. The public
   surface is `include/InitialCondition.h` (base class + `ICRegistry`). If an
   IC header is ever needed outside `src/ic/`, promote it to `include/`.
2. **Generated headers.** `include/Constants.h` ← `Constants.h.orig`;
   `include/UserSource.h` ← `userfiles/<Name>Source.h` via `Config.pl -u=<Name>`;
   `include/show_git_info.h` is created at build time. Never edit the outputs.
3. **Split `Particles` translation units** each need their own explicit
   two-alias member instantiations — see the hard constraints in `SKILL.md`.
4. **Fake 2D** = a single cell in z; `nDim` stays 3 unless AMReX itself is
   built 2D (`./Config.pl -amrex2d`).
5. **Component macros:** `_PC_COMPONENT_` for PIC builds, `_PT_COMPONENT_` for
   the particle tracker.
