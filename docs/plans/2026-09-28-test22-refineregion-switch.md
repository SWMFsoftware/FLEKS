# Test22 Refinement Selector Switch Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make the OH–PT/FLEKS test22 case move refinement by changing `#REFINEREGION` between sessions.

**Architecture:** Define immutable named boxes in the first PT session. Keep the existing first-session selectors and switch level 1 from `+box2` to `+box3` in session 2. Run the complete coupled test and inspect both the numerical comparison and mesh logs.

**Tech Stack:** SWMF parameter files, Makefile.test, FLEKS AMReX logging.

---

### Task 1: Add the selector transition

**Files:** Modify `../../Param/PARAM.in.test.OHPT.FLEKS.outerhelio.swhpui`.

1. Confirm the current local parameter file lacks a session-2 selector change (`rg -n '#RUN|#REFINEREGION' ...`); this is the failing baseline.
2. Add `box2` and `box3` definitions with distinct names in the initial PT section; preserve the previous box geometry and use the shifted x extent for `box3`.
3. Set initial level 1 to `+box2`, then use `#REFINEREGION` level 1 `+box3` after `#RUN`.
4. Preserve a one-year first endpoint and two-year final endpoint.

### Task 2: Verify the coupled regression

**Files:** No source changes unless the run identifies a FLEKS defect.

1. Run `make test22 TESTDIR=run_codex_test22 COMPILE.mpicxx=/opt/local/libexec/mpich-gcc15/mpicxx LINK.f90=/opt/local/libexec/mpich-gcc15/mpif90` from the SWMF root.
2. Check `SWMF.SUCCESS`, `test22_ohpt.diff`, and the two PT fine-grid summaries in `RESULTS/swhpui/runlog`.
3. If reference output differs because the test now exercises a different mesh, inspect the diff and update only the relevant reference after confirming the run is valid.
4. Run `git diff --check` and review repository status before committing the relevant changes.
