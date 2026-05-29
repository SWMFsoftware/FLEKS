# Standalone FLEKS Tests

This directory contains standalone tests for FLEKS (Flexible Exascale Kinetic Simulator) to verify exosphere injection, physical transport, boundary conditions, and pickup ion acceleration independent of SWMF coupling.

These tests are executed as part of the verification and CI workflow to ensure the core particle-in-cell (PIC) solver remains functional and accurate.

## Directory Structure

Each test case is contained within its own dedicated subdirectory containing a `PARAM.in` configuration file and a dedicated `README.md` describing the physics, analytical model, and expected results:

*   **[box/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/box/)**: 0D equilibrium test validating spatial injection rates, particle population growth, and grid cell uniform loading.
*   **[chamber/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/chamber/)**: 1D chamber flow test checking inflow boundary, free outflow boundary, and convergence to a steady-state particle count.
*   **[photoionization/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/photoionization/)**: Exospheric photoionization test validating automatic pair-injection of ions and neutralizing electrons to ensure strict local and global charge conservation.
*   **[electron_impact/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/electron_impact/)**: Exospheric electron impact ionization test verifying collision probabilities, Lotz cross-sections, and Opal-Beaty energy partitioning using a Monte Carlo Collision (MCC) algorithm.
*   **[charge_exchange/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/charge_exchange/)**: Exospheric charge exchange test validating ion drift deceleration and thermal cooling through in-place velocity replacements using a Monte Carlo Collision (MCC) algorithm.
*   **[exosphere/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/exosphere/)**: Combined exosphere processes test validating exospheric photoionization pair-injection, electron impact ionization MCC, and exospheric charge exchange MCC running concurrently.
*   **[beam/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/beam/)**: 1D ion beam electromagnetic instability test validating cyclotron wave growth, transverse magnetic field amplification, and energy conservation.
*   **[tophat/](file:///home/hyzhou/simulation/SWMF/PC/FLEKS_myversion/tests/test_standalone/tophat/)**: 1D electromagnetic top-hat wave propagation test validating $B_z$ and $E_y$ field dynamics under the EM solver with zero particles.


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

