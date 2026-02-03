---
name: Generate Docs
description: Build documentation from code comments and LaTeX sources
---

# Generate Docs

This skill handles documentation generation for FLEKS.

## Documentation Locations

| Type | Location |
|------|----------|
| Algorithm documentation | `documents/Algorithm.tex` |
| Coding standards | `documents/Coding_standards.md` |
| Parameter documentation | `PARAM.XML` |

## LaTeX Documentation

### Build Algorithm PDF

```bash
cd documents
pdflatex Algorithm.tex
```

For full build with references:
```bash
cd documents
pdflatex Algorithm.tex
pdflatex Algorithm.tex  # Run twice for references
```

### View Generated PDF

```bash
open documents/Algorithm.pdf  # macOS
```

## Code Documentation with Doxygen

If you want to add Doxygen-style documentation:

### 1. Install Doxygen

```bash
brew install doxygen  # macOS
```

### 2. Create Doxyfile

```bash
cd /Users/yuxichen/shock/SWMF/PC/FLEKS
doxygen -g
```

### 3. Configure Doxyfile

Key settings to modify:
```
PROJECT_NAME = "FLEKS"
INPUT = include src
FILE_PATTERNS = *.h *.cpp
EXTRACT_ALL = YES
GENERATE_LATEX = NO
```

### 4. Generate Documentation

```bash
doxygen Doxyfile
```

Output will be in `html/index.html`.

## Doxygen Comment Style

### For Classes

```cpp
/**
 * @brief A brief description of the class.
 * 
 * A more detailed description of what this class does,
 * its purpose, and usage patterns.
 */
class MyClass {
```

### For Functions

```cpp
/**
 * @brief Brief description of the function.
 * 
 * @param param1 Description of first parameter
 * @param param2 Description of second parameter
 * @return Description of return value
 * 
 * @note Any important notes about usage
 * @warning Any warnings about side effects
 */
int my_function(int param1, double param2);
```

### For Member Variables

```cpp
class MyClass {
private:
  int count_;  ///< Number of items processed
  double tolerance_;  ///< Convergence tolerance for solver
};
```

## PARAM.XML Documentation

The `PARAM.XML` file contains input parameter documentation. Format:

```xml
<command name="PARAMETER_NAME">
  <parameter name="Value" type="real" default="1.0"/>
  <description>
    Description of what this parameter does.
  </description>
</command>
```

## README Updates

### Adding Usage Examples

```markdown
## Quick Start

1. Configure the simulation:
   ```bash
   ./Config.pl -install
   ```

2. Build FLEKS:
   ```bash
   make FLEKS -j8
   ```

3. Run:
   ```bash
   cd ../run
   ./SWMF.exe
   ```
```

## Documentation Checklist

- [ ] All public classes have brief descriptions
- [ ] All public functions have parameter documentation
- [ ] Complex algorithms have detailed comments
- [ ] `PARAM.XML` is up-to-date with new parameters
- [ ] `Algorithm.tex` reflects current implementation
- [ ] README has updated build/run instructions
