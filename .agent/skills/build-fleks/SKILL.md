---
name: Build FLEKS
description: Compile the FLEKS project with proper configuration and regenerate compile_commands.json
---

# Build FLEKS

This skill compiles the FLEKS project and ensures IntelliSense stays up to date.

## Prerequisites

- `mpicxx` must be available in PATH
- AMReX must be installed; the path is set via `AMREXDIR` in
  `Makefile.def` (typically `../../util/AMREX/InstallDir/`)
- The SWMF `Makefile.conf` must exist (referenced by `./Makefile.conf`)

## Steps

### 1. Working Directory

All commands should be run from the FLEKS project root (i.e., the
directory containing this `.agent/` folder).

### 2. Clean Build (Optional)

If a clean build is requested, run:
```bash
make clean
```

For a complete reset (removes all generated files):
```bash
make distclean
```

### 3. Build the Project

For standalone executable:
```bash
make FLEKS -j8
```

For library (used with SWMF):
```bash
make LIB -j8
```

Note: The `compile_commands` target is automatically invoked by both
`FLEKS` and `LIB` targets.

### 4. Check for Errors

- Review the build output for any compilation errors
- Common issues:
  - Missing headers → Check include paths in `src/Makefile` (`SEARCH_C`)
    and `.vscode/c_cpp_properties.json`
  - Undefined symbols → Check if all required `.cpp` files are listed in
    the `SRCS` variable in `src/Makefile`
  - AMReX errors → Verify AMReX installation path in `Makefile.def`
    (`AMREXDIR` / `AMREXSEARCH`)

### 5. Verify compile_commands.json

After successful build, verify that `compile_commands.json` was regenerated:
```bash
ls -la compile_commands.json
```

This file is used by VS Code's C/C++ IntelliSense.

## Targets Reference

| Target | Description |
|--------|-------------|
| `make FLEKS` | Build standalone executable (`bin/FLEKS.exe`) |
| `make LIB` | Build library for SWMF integration (`src/libFLEKS.a`) |
| `make CONVERTER` | Build converter tool (`bin/converter.exe`) |
| `make clean` | Remove object files |
| `make distclean` | Remove all generated files |
| `make compile_commands` | Regenerate `compile_commands.json` only |

## Key Compilation Flags

The C++ compilation flags are set in `Makefile.conf` (via SWMF):

| Variable | Typical Value | Purpose |
|----------|---------------|---------|
| `DEBUGC` | `-g -Wall -Wextra -Wno-unused-parameter` | Debug symbols and warnings |
| `OPT3` | `-O0` (debug) or `-O3` (release) | Optimization level |
| `FLAGC_EXTRA` | `-D_PC_COMPONENT_` | Component preprocessor flag |

To enable a debug build, ensure `OPT3` is set to `-O0` in the SWMF
`Makefile.conf` (it is `-O0` by default during development).

## Success Criteria

- Build completes without errors
- `src/libFLEKS.a` is updated (check timestamp)
- `compile_commands.json` is present and current
