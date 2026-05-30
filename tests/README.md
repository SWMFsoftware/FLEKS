# FLEKS Testing Framework

This directory contains standalone validation and physics tests for **FLEKS (Flexible Exascale Kinetic Simulator)**. These tests verify the core Particle-in-Cell (PIC) solver's mathematical models, physical transport, boundary conditions, and collision/ionization processes completely independent of SWMF coupling.

---

## 🚀 Running the Tests

A unified Python runner is provided to dynamically discover, execute, and validate the test suite:

```bash
python3 tests/test_standalone/validate_tests.py
```

### What the test runner does:
1. Scans the `tests/test_standalone/` directory for any subdirectories containing a `PARAM.in` configuration file.
2. Initializes the execution directory `run_test/` with necessary symlinks to the FLEKS executable and post-processing tools.
3. Copies the respective `PARAM.in` file, runs `./FLEKS.exe`, and invokes the post-processing pipeline `./PostProc.pl`.
4. Parses diagnostic outputs and validates physical correctness against analytical models or exit-code verifications.

---

## 🛠️ Prerequisites & Compilation

Before running the tests, both the main **FLEKS executable** and the **IDL post-processor (`PostIDL.exe`)** must be compiled.

### 1. Compile the FLEKS Executable
Ensure the standalone FLEKS executable is compiled first:
```bash
make EXE
```

### 2. Compile the Post-Processor (`PIDL`)
The test runner executes post-processing scripts which require the `PostIDL.exe` executable. If you encounter the following error:
> `ERROR in pIDL: .../PC/PostIDL.exe is not available, please make PIDL`

You can build the post-processor from the FLEKS root directory by running:
```bash
make -C share/Library/src PIDL
```
*(This compiles `PostIDL.exe` and places the binary under the `bin/` directory, resolving the symlink dependency).*

---

## 📂 Test Catalog

The test suite in `tests/test_standalone/` covers a wide range of physical and numerical configurations:

*   **`box/`**: 0D equilibrium test validating spatial particle injection, population growth, and uniform loading.
*   **`chamber/`**: 1D inflow and free outflow boundary test verifying convergence to steady-state particle count.
*   **`photoionization/`**: Validates pair-injection of ions and electrons to ensure local and global charge conservation.
*   **`electron_impact/`**: Verifies collision probabilities and Lotz cross-sections using a Monte Carlo Collision (MCC) algorithm.
*   **`charge_exchange/`**: Validates ion drift deceleration and thermal cooling via MCC velocity replacements.
*   **`exosphere/`**: Combined test verifying exospheric photoionization, electron impact, and charge exchange running concurrently.
*   **`beam/`**: 1D ion beam electromagnetic instability test validating cyclotron wave growth and energy conservation.
*   **`tophat/`**: 1D electromagnetic wave propagation test validating $B_z$ and $E_y$ field dynamics under the EM solver with zero particles.
*   **Wave & Plasma Modes**: Includes tests for `alfven_wave`, `fast_wave`, `slow_wave`, `sound_wave`, and `langmuir` oscillations to verify plasma physics modes.
