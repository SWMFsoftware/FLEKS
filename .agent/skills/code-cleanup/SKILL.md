---
name: Code Cleanup
description: Find unused variables, check formatting, apply fixes following FLEKS coding standards
---

# Code Cleanup

This skill performs code cleanup tasks following FLEKS coding standards.

## Coding Standards Reference

From `documents/Coding_standards.md`:

1. **Memory Management**: Use `shared_ptr` or `unique_ptr`, avoid raw `new`
2. **Naming Conventions**:
   - Files: `FileName.cpp`, `HeaderName.h`
   - Classes: `ClassName`
   - Variables: `variableName`
   - Functions: `this_is_a_function_name`
3. **Namespace**: `using namespace amrex` allowed only in `.cpp` files
4. **Header Order**: std headers → AMReX headers → user headers
5. **Pointers**: Use `nullptr`, not `NULL`
6. **Const**: Always use `const` when possible
7. **Lambdas**: Prefer regular functions for universal or long functions
8. **Commits**: Follow [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/)

## Clang-Format

The project uses `.clang-format` with Mozilla-based style (2-space indent, 80-column limit).

### Format Single File

```bash
clang-format -i src/FileName.cpp
```

### Format All Source Files

```bash
find src include -name "*.cpp" -o -name "*.h" | xargs clang-format -i
```

### Check Formatting Without Changing

```bash
clang-format --dry-run -Werror src/FileName.cpp
```

## Finding Issues

### 1. Unused Variables

Use compiler warnings:
```bash
# Add to compile command
-Wall -Wextra -Wunused-variable -Wunused-parameter
```

Or use clang-tidy:
```bash
clang-tidy src/FileName.cpp -- -I../include -I../../util/AMREX/InstallDir/include
```

### 2. Find Raw `new` Usage

```bash
grep -rn "\bnew\b" src/ include/ --include="*.cpp" --include="*.h"
```

### 3. Find `NULL` Usage (should be `nullptr`)

```bash
grep -rn "\bNULL\b" src/ include/ --include="*.cpp" --include="*.h"
```

### 4. Find `using namespace` in Headers

```bash
grep -rn "using namespace" include/ --include="*.h"
```

### 5. Check Naming Conventions

Look for:
- Functions not using `snake_case`
- Variables not using `camelCase`
- Classes not using `PascalCase`

## Automated Cleanup Tasks

### Task 1: Fix NULL to nullptr

```bash
# Preview changes
grep -rn "\bNULL\b" src/ include/

# Apply fix (carefully review first!)
find src include -name "*.cpp" -o -name "*.h" | xargs sed -i '' 's/\bNULL\b/nullptr/g'
```

### Task 2: Remove Trailing Whitespace

```bash
find src include -name "*.cpp" -o -name "*.h" | xargs sed -i '' 's/[[:space:]]*$//'
```

### Task 3: Ensure Newline at End of File

```bash
for f in src/*.cpp include/*.h; do
  [ -n "$(tail -c1 "$f")" ] && echo >> "$f"
done
```

## Review Checklist

When reviewing cleanup changes:

- [ ] No functional changes introduced
- [ ] All files still compile
- [ ] Naming conventions followed
- [ ] No new warnings introduced
- [ ] Header order is correct (std → AMReX → user)
- [ ] `const` used where applicable
- [ ] No raw pointers managing ownership

## Commit Message Format

For cleanup commits, use conventional commit format:
```
refactor: remove unused variables in Pic.cpp

- Removed unused `tempVar` variable
- Fixed NULL → nullptr
- Applied clang-format
```
