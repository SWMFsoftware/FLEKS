---
name: Build FLEKS
description: Compile the FLEKS project with proper configuration and regenerate compile_commands.json
---

# Build FLEKS

This skill compiles the FLEKS project and ensures IntelliSense stays up to date.

## Prerequisites

- `mpicxx` must be available in PATH
- AMReX must be installed at `../../util/AMREX/InstallDir/`

## Steps

### 1. Navigate to FLEKS Directory

All commands should be run from `/Users/yuxichen/shock/SWMF/PC/FLEKS`

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

Note: The `compile_commands` target is automatically invoked by both `FLEKS` and `LIB` targets.

### 4. Check for Errors

- Review the build output for any compilation errors
- Common issues:
  - Missing headers → Check include paths in `c_cpp_properties.json`
  - Undefined symbols → Check if all required `.cpp` files are compiled
  - AMReX errors → Verify AMReX installation at `../../util/AMREX/InstallDir/`

### 5. Verify compile_commands.json

After successful build, verify that `compile_commands.json` was regenerated:
```bash
ls -la compile_commands.json
```

This file is used by VS Code's C/C++ IntelliSense.

## Targets Reference

| Target | Description |
|--------|-------------|
| `make FLEKS` | Build standalone executable |
| `make LIB` | Build library for SWMF integration |
| `make CONVERTER` | Build converter tool |
| `make clean` | Remove object files |
| `make distclean` | Remove all generated files |
| `make compile_commands` | Regenerate compile_commands.json only |

## Success Criteria

- Build completes without errors
- `src/libFLEKS.a` is updated (check timestamp)
- `compile_commands.json` is present and current
