# FLEKS Unit Tests

This directory is the future home for standalone unit testing of the FLEKS codebase (independent of the overarching SWMF macro regression tests).

## Framework Guidance
The recommended framework for internal kinetic math and grid testing is **GoogleTest (gtest)** or **Catch2**. 

When adding tests here in the future:
1. Initialize the framework of your choice (e.g., via a standalone CMake build system targeting this folder).
2. Tests should execute entirely isolated from SWMF Fortran wrappers (focus purely on `src/` classes like `Particles.cpp`, `LinearSolver.cpp`, and mathematical `GridUtility` routines).
3. If leveraging AMReX data arrays, ensure your test runner initializes AMReX via `amrex::Initialize()` before the test suite and calls `amrex::Finalize()` during test tear-down.

## Macro Regression Tests
If you simply want to test whether FLEKS correctly couples with BATS-R-US or global SWMF targets, please navigate to the root SWMF directory and run `make test16_3d` or other pre-compiled XML routines in the SWMF test catalog.
