# Recent AMR Comments and test23 Verification Plan

**Goal:** Explain the recent dynamic AMR implementation at its non-obvious decision points and verify both coupled fast-wave AMR tests on the current code.

**Architecture:** Keep comments close to the code they explain, without changing behavior. Run the existing 2D and 3D SWMF targets in separate run directories, then check that session two regrids and that the refined mesh changes in addition to checking the standard numerical comparisons.

**Tech Stack:** C++, AMReX, SWMF `make test23_2d` / `test23_3d`.

---

### Task 1: Comment recent code

1. Review FLEKS commits from the prior 48 hours and the affected functions.
2. Add concise why/order comments around persistent shape definitions, selector changes, fake-2D ratios, regrid propagation, fine-cell initialization, and session handling.
3. Run `clang-format --dry-run --Werror` on changed C++ files and `git diff --check`.

### Task 2: Coupled 2D verification

1. Run `make test23_2d TESTDIR=run_codex_test23_2d` from the SWMF root.
2. Inspect the run log for both sessions and their refinement grids; inspect the numerical diff and run status.
3. If it fails, trace the first error before changing behavior.

### Task 3: Coupled 3D verification

1. Run `make test23_3d TESTDIR=run_codex_test23_3d` from the SWMF root.
2. Apply the same session, mesh, diff, and status checks.
3. Preserve existing unrelated working-tree changes.
