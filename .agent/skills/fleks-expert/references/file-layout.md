# FLEKS File Layout

> Canonical: `Doc/DEVELOPING.md` §1 (Repository Layout) and §3 (Code Layout and
> Entry Points). Open the doc for the full tables — this reference only holds
> the shortcuts and the rules that are easy to get wrong.

## Where to look

| Need | Open |
|---|---|
| Time loop / regridding / orchestration | `src/Pic.cpp` (`Pic::update`), `src/Domain.cpp` |
| Parameter parsing | `src/PicParam.cpp`, `src/Domain.cpp` (`read_param`) |
| Field boundary conditions | `src/PicBC.cpp`, `include/BC.h` |
| Full-PIC field solve / div(E) | `src/PicFieldSolver.cpp`, `src/PicDivE.cpp` |
| Hybrid solver | `src/PicHybrid.cpp` |
| Particles | `src/Particles.cpp` + `src/Particles{Init,BC,Moments,MassMatrix,Mover,Resample,Reactions}.cpp` |
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
