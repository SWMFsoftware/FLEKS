# FLEKS Build Reference

> The actionable recipe (compile, diagnose, verify) is
> `.agents/skills/build-fleks/SKILL.md`; this page holds the configuration
> surface: what `Config.pl` can set, which target produces what, and how
> `src/Makefile` is wired.

## Configuration

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

## Targets

| Target | Result |
|---|---|
| `make EXE -j8` (alias `make FLEKS`) | standalone `bin/FLEKS.exe` |
| `make LIB -j8` | SWMF component library `src/libFLEKS.a` + wrappers |
| `make CONVERTER` | `bin/converter.exe` (delegates to `make -C Converter`) |
| `make compile_commands` | regenerate `compile_commands.json` (IDE index) |
| `make clean` / `make distclean` | object files / full reset |
| `make PDF` | `doc/USERMANUAL.pdf` (needs `pdflatex`, `makeindex`, `fvextra`) |

`make LIB` is only valid inside a built SWMF tree (needs `libSHARE.a` and
`con_comp_param.mod`); use `make EXE` for standalone work. `EXE`, `FLEKS` and
`LIB` invoke `compile_commands` automatically. `CONVERTER` is built separately.

## Dependencies

| Dependency | Required | Notes |
|---|---|---|
| AMReX | yes | `../../util/AMREX/InstallDir/` in SWMF, `util/AMREX/InstallDir/` standalone |
| MPI | yes | `mpicxx`/`mpif90` on PATH |
| SWMF share | yes | `share/Scripts`, `share/Library`, `libSHARE.a` |
| HDF5 | optional | parallel HDF5 output |
| Perl | yes | `Config.pl` and `share/Scripts/` |

## src/Makefile

`SRCS` is a wildcard over `src/*.cpp` and `src/ic/*.cpp`, minus `main.cpp`
(which has its own `main()`), wrapped in `$(sort)` so the archive order is
reproducible. A new FLEKS source needs no Makefile edit — and a stray `.cpp`
in `src/` will be compiled. Useful variables:

| Variable | Purpose |
|---|---|
| `SRCS` | sources compiled into `libFLEKS.a` (auto-discovered) |
| `SEARCH_C` | include search paths (e.g. `-I../include`) |
| `FLAGC_EXTRA` | `-D_${COMPONENT}_COMPONENT_` |
| `LIBFLEKS` | output library name |

## Converter/Makefile

The format converter is a separate program with its own `main()` and its own
Makefile at `Converter/Makefile`; it is **not** part of `libFLEKS.a`.

```bash
make CONVERTER              # from the root (delegates to -C Converter)
make -C Converter CONVERTER # directly, e.g. while iterating on the tool
```

| Variable | Purpose |
|---|---|
| `SRCS` | converter-only sources: `Converter.cpp`, `DataContainer.cpp`, `VisitWriter.cpp` |
| `SHARED_SRCS` | `../src/Grid.cpp`, `../src/FleksDistributionMap.cpp` |

`SHARED_SRCS` exist because `AMReXDataContainer` derives from `Grid`. They are
compiled into `Converter/` rather than taken from `libFLEKS.a`: the full
archive also holds the PT/OH coupling objects, which reference OH symbols that
only exist inside SWMF. Their objects must stay local — do not let make write
`../src/Grid.o`, it would clobber the library's own object.

Gotcha: `default:` must be the first target in `Converter/Makefile`.
`../Makefile.conf` defines a `cleanfiles:` target and `include` expands in
place, so anything after the includes is no longer the default goal.

Preprocessor flags: `_PC_COMPONENT_` (PIC build), `_PT_COMPONENT_` (particle
tracker), `_USE_HDF5_` (HDF5 output).

## Standalone runs

`src/main.cpp` is the driver. Update it for driver-level behaviour, keep
`Domain` logic shared with SWMF whenever possible, rebuild with `make EXE -j8`,
and run from a directory containing `PARAM.in`. Standalone runs use domain name
`FLEKS1` (plots and restarts under `FLEKS1/`) and need `#INITFROMSWMF F` plus
`#NORMALIZATION`, `#PLASMA` and `#UNIFORMSTATE`.
