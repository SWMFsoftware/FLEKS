# FLEKS Coding Standards Reference

> Canonical: `Doc/Coding_standards.md` plus `Doc/DEVELOPING.md` §3
> (Coding conventions). This page adds the FLEKS-specific rules that are not in
> the generic style guide.

Full text for humans: `Doc/Coding_standards.md`. Summary below.

## Naming

| Element | Convention | Example |
|---|---|---|
| Files | PascalCase | `GridUtility.cpp` |
| Classes | PascalCase | `FluidInterface` |
| Functions | snake_case | `apply_float_boundary()` |
| Variables | camelCase | `nCellPerPatch` |
| Constants | camelCase / UPPER | `cProtonMassSI`, `nDim3` |
| Private members | camelCase | `doRestart` |

## Style rules

1. **Smart pointers** — `unique_ptr`/`shared_ptr` for ownership; raw pointers
   only when there is no ownership. Never raw `new`.
2. **`using namespace amrex`** — allowed only in `.cpp` files, never in headers.
3. **Header order** — std → AMReX → project headers.
4. **`nullptr`**, never `NULL`.
5. **`const`** wherever possible.
6. **Include guards** — `#ifndef _FILENAME_H_` / `#define _FILENAME_H_`.
7. **80-column limit**, 2-space indent, Mozilla-based `.clang-format`; Fortran
   indented with `findent` to match Emacs `f90-mode`.
8. **Commits** — [Conventional Commits](https://www.conventionalcommits.org/).

## FLEKS-specific rules

- **Split `Particles` files:** every new out-of-line `Particles` method needs
  explicit instantiations for both `PicParticles` and `PTParticles` in the
  file that defines it. Only `src/Particles.cpp` may hold whole-class
  instantiations.
- **Generated files:** edit `Constants.h.orig` / `userfiles/*Source.h`, never
  `include/Constants.h` / `include/UserSource.h`.
- **Source list:** add every new `.cpp` to `SRCS` in `src/Makefile`.
- **Ionization data** belongs in `SourceInterface`/`UserSource`, not
  `FluidInterface`.
- **No `#` in comments** inside `PARAM.in` when the comment references a
  command (it would be parsed as a command).

## Formatting commands

```bash
python3 tools/format_all.py                    # C++ (clang-format) + Fortran (findent)
clang-format --dry-run -Werror src/Pic.cpp     # check only
pip install clang-format==20.1.8 findent       # match CI versions
```

## SWMF interface patterns

1. SWMF calls `PC_wrapper.f90` / `PT_wrapper.f90` subroutines.
2. The Fortran wrapper calls C functions in `FleksInterface.cpp`
   (e.g. `fleks_run_`).
3. C++ functions must use C-compatible conventions: `extern "C"` in C++,
   `bind(C)` in Fortran, and they operate on the global `fleksDomains` object.

More detail in `references/coupling.md`.
