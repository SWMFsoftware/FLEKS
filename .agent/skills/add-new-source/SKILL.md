---
name: Add New Source
description: Create new .cpp/.h files following FLEKS project conventions and coding standards
---

# Add New Source

This skill creates new source files following FLEKS conventions.

## File Locations

| File Type | Directory |
|-----------|-----------|
| Headers (`.h`) | `include/` |
| Implementation (`.cpp`) | `src/` |
| Interface code | `srcInterface/` |

## Naming Conventions

- **File names**: `PascalCase` (e.g., `FluidInterface.cpp`, `GridUtility.h`)
- **Class names**: `PascalCase` (e.g., `class FluidInterface`)
- **Function names**: `snake_case` (e.g., `void apply_float_boundary()`)
- **Variable names**: `camelCase` (e.g., `int nCellPerPatch`)
- **Private members**: `camelCase` (e.g., `bool doRestart`)

## Header File Template

Create `include/NewFeature.h`:

```cpp
#ifndef _NEWFEATURE_H_
#define _NEWFEATURE_H_

// Standard headers
#include <memory>
#include <vector>

// AMReX headers
#include <AMReX.H>

// User headers
#include "Utility.h"

class NewFeature {
public:
  NewFeature() = default;
  ~NewFeature() = default;

  // Public methods
  void initialize();
  void do_something();

  // Getters (use const)
  int get_value() const { return value_; }

private:
  // Member variables (trailing underscore for truly private state)
  int value_ = 0;
};

#endif // _NEWFEATURE_H_
```

## Implementation File Template

Create `src/NewFeature.cpp`:

```cpp
#include "NewFeature.h"

// Standard headers
#include <iostream>

// AMReX headers (can use 'using namespace amrex' in .cpp only)
using namespace amrex;

void NewFeature::initialize() {
  // Implementation
}

void NewFeature::do_something() {
  // Implementation
}
```

## Steps to Add a New File

### 1. Create Header File

Create `include/NewFeature.h` with the template above. Ensure proper
include guards (`#ifndef _NEWFEATURE_H_`).

### 2. Create Implementation File

Create `src/NewFeature.cpp` with the template above.

### 3. Register in Makefile

**Important:** Add the new `.cpp` file to the `SRCS` variable in
`src/Makefile`. The source list starts at line 5:

```makefile
SRCS := \
	Domain.cpp \
	...
	FleksDistributionMap.cpp \
	NewFeature.cpp           # <-- add here
```

This is required — files in `src/` are NOT auto-discovered.

### 4. Rebuild and Verify

```bash
make LIB -j8
```

### 5. Regenerate compile_commands.json

This happens automatically with the build, but you can force it:
```bash
make compile_commands
```

## Header Order Standard

Always order includes as:

```cpp
// 1. Standard library headers
#include <algorithm>
#include <memory>
#include <vector>

// 2. AMReX headers
#include <AMReX.H>
#include <AMReX_MultiFab.H>

// 3. Project headers
#include "Grid.h"
#include "Utility.h"
```

## Key Reminders

1. **Use `nullptr`** instead of `NULL`
2. **Use smart pointers** (`unique_ptr`, `shared_ptr`) for ownership
3. **Use `const`** wherever possible
4. **No `using namespace`** in header files
5. **80-column limit** (enforced by `.clang-format`)

## Integrating with Existing Classes

If your new class needs to work with existing code:

1. Check if there's a base class or interface to inherit from
2. Look at similar classes for patterns (e.g., `Grid.h`, `Domain.h`)
3. Common base classes:
   - `Grid` — for classes that need AMR grid access (inherits `AmrCore`)
   - `FluidInterface` — for fluid coupling data
   - `SourceInterface` — for source terms

## Verification Checklist

- [ ] Header has include guards (`#ifndef _FILENAME_H_`)
- [ ] File names follow `PascalCase`
- [ ] Class name matches file name
- [ ] Header order is correct (std → AMReX → user)
- [ ] No `using namespace` in header
- [ ] `const` used where applicable
- [ ] Smart pointers for ownership
- [ ] **File added to `SRCS` in `src/Makefile`**
- [ ] Code compiles without warnings
- [ ] `compile_commands.json` updated
