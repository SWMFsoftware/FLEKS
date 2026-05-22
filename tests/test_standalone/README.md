# Standalone FLEKS Tests

This directory contains standalone tests for FLEKS (Flexible Exascale Kinetic Simulator) to verify exosphere injection, physical transport, boundary conditions, and pickup ion acceleration independent of SWMF coupling.

These tests are executed as part of the verification and CI workflow to ensure the core particle-in-cell (PIC) solver remains functional and accurate.

## Directory Structure

Each test case is contained within its own dedicated subdirectory containing a `PARAM.in` configuration file and a dedicated `README.md` describing the physics, analytical model, and expected results:

*   **[box/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/box/)**: 0D equilibrium test validating spatial injection rates, particle population growth, and grid cell uniform loading.
*   **[chamber/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/chamber/)**: 1D chamber flow test checking inflow boundary, free outflow boundary, and convergence to a steady-state particle count.
*   **[pickup/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/pickup/)**: Multi-species pickup ion acceleration test verifying leapfrog mover physics ($E \times B$ drift) for H+ and O+ alongside charge-neutralizing electrons under constant perpendicular background fields.

## Running the Tests

A unified Python runner is provided to dynamically discover, run, and validate all tests:

```bash
python3 tests/test_standalone/validate_tests.py
```

The script:
1. Scans all subdirectories of `tests/test_standalone/` dynamically for any directory containing a `PARAM.in` file.
2. Creates the execution directory `run_test/` with necessary symlinks to the FLEKS executable and post-processing tools.
3. Copies `PARAM.in` from the test folder, executes `./FLEKS.exe`, and runs `./PostProc.pl`.
4. Parses diagnostic outputs and validates results using custom analytical models (e.g. discrete time-integration for pickup ions, steady-state checks for chamber flow) or general exit-code verification.

