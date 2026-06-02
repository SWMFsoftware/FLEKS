# Standalone FLEKS Tests

This directory contains standalone tests for FLEKS (Flexible Exascale Kinetic
Simulator). The active standalone test suite focuses on the core
particle-in-cell (PIC) solver independent of SWMF coupling.

## Directory Structure

Each test case is contained within its own dedicated subdirectory containing a
`PARAM.in` configuration file and a README describing the expected behavior:

* **beam/**: 1D ion beam electromagnetic instability test validating cyclotron
  wave growth, transverse magnetic field amplification, and energy conservation.
* **performance/**: beam-based performance benchmark used by the dedicated
  performance workflow.

## Running the Tests

A unified Python runner is provided to dynamically discover, run, and validate
the standalone tests:

```bash
python3 tests/test_standalone/validate_tests.py
```

The script:
1. Scans subdirectories of `tests/test_standalone/` for a `PARAM.in` file,
   excluding `performance/`.
2. Creates the execution directory `run_test/` with necessary symlinks to the
   FLEKS executable and post-processing tools.
3. Copies `PARAM.in` from the test folder, executes `./FLEKS.exe`, and runs
   `./PostProc.pl`.
4. Parses diagnostic outputs and validates the beam physics checks.
