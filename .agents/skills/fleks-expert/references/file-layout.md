# FLEKS File Layout

> Navigation aid for the source tree. Pair it with `references/build.md` (how
> the tree is built), `references/architecture.md` (what the classes do) and
> `references/tools.md` (post-processing).

## Repository layout

```
FLEKS/
├── include/           Public headers (.h)
├── src/               C++ sources (.cpp), src/Makefile (SRCS list), main.cpp
│   └── ic/            Initial-condition plug-ins (private headers live here)
├── srcInterface/      SWMF coupling layer: PC_wrapper.f90, PT_wrapper.f90,
│                      FleksInterface.cpp
├── doc/               Algorithm.tex, Coding_standards.md, Tex/
├── tools/             Post-processing, conversion and formatting scripts
├── tests/             Standalone test suite (one directory per scenario)
├── userfiles/         Selectable user-source templates (*Source.h)
├── .agents/            Agent skills, workflows and knowledge base
├── Config.pl          Perl configuration (AMReX, AMR levels, user source, …)
├── Makefile           Top-level makefile
├── Makefile.def.FLEKS Default FLEKS makefile definitions
├── PARAM.XML          Parameter command reference (source of the user manual)
└── .clang-format      Mozilla-based style, 80 columns, 2-space indent
```

## Where to look

| Need | Open |
|---|---|
| Time loop / regridding / orchestration | `src/Pic.cpp` (`Pic::update`), `src/Domain.cpp` |
| Parameter parsing | `src/PicParam.cpp`, `src/Domain.cpp` (`read_param`) |
| Field boundary conditions | `src/PicBC.cpp`, `include/BC.h` |
| Full-PIC field solve / div(E) | `src/PicFieldSolver.cpp`, `src/PicDivE.cpp` |
| Hybrid solver | `src/PicHybrid.cpp` |
| Particles | `src/Particles.cpp` plus the `src/Particles*.cpp` split files |
| Fluid / coupling state | `src/FluidInterface.cpp` |
| Output | `src/PlotWriter.cpp`, `src/PicIO.cpp` |
| AMR grid / load balancing | `src/Grid.cpp`, `src/FleksDistributionMap.cpp` |
| Standalone driver | `src/main.cpp` |
| SWMF entry points | `srcInterface/FleksInterface.cpp` |

## Rules agents most often get wrong

1. **`src/ic/` headers are private.** IC plug-in headers stay next to their
   `.cpp`; the public surface is `include/InitialCondition.h`. Promote a header
   to `include/` only if something outside `src/ic/` needs it.
2. **Generated files are not editable.** `include/Constants.h` ←
   `Constants.h.orig`; `include/UserSource.h` ← `userfiles/*Source.h`;
   `include/show_git_info.h` is created at build time.
3. **New `.cpp` files must be added to `SRCS`** in `src/Makefile` — nothing is
   auto-discovered, and the omission only surfaces at link time.
4. **Fake 2D** is one cell in z; a true-2D AMReX build needs
   `./Config.pl -amrex2d`.
5. **Component macros** `_PC_COMPONENT_` / `_PT_COMPONENT_` — guarded code must
   compile both ways.
