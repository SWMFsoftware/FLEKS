# AMR Region Selection Implementation Plan

**Goal:** Change refinement only when `#REFINEREGION` selects a different region, and reject any second definition of a named shape.

**Architecture:** Keep named shapes in a persistent registry across parameter sessions. A `#REGION` declaration inserts a new name or fails if that name already exists. `#REFINEREGION` selectors remain level-specific and are compared with the previous selector; only selector or grid-efficiency changes mark refinement modified. Existing active-grid updates remain independent.

**Tech Stack:** C++, AMReX, standalone FLEKS regression tests.

---

### Task 1: Duplicate shape regression

**Files:** Create `tests/dynamic_amr/test_region_names.py`.

1. Add a standalone test with two parameter sessions that define the same `#REGION` name.
2. Run it against the current executable and confirm it fails because the redefinition is accepted.
3. Keep the test for both identical and changed definitions.

### Task 2: Immutable shapes and selector-driven regrid

**Files:** Modify `include/DomainGrid.h` and `src/Domain.cpp`.

1. Replace shape signature/upsert handling with a persistent set of shape names and an insert helper that aborts on duplicates.
2. Remove changed-shape tracking from `read_param`; retain selector and grid-efficiency change tracking.
3. Build the standalone executable and rerun the duplicate-name regression.

### Task 3: Positive dynamic AMR regression and documentation

**Files:** Modify `tests/dynamic_amr/PARAM.in`, `tests/dynamic_amr/README.md`, and `PARAM.XML`.

1. Define `left` and `right` once, then change the level-0 selector from `+left` to `+right` and finally to `none`.
2. Update documentation to state that a shape name cannot be redefined.
3. Run the focused standalone dynamic AMR test and inspect regrid diagnostics.
4. Check formatting, diff, and working-tree status; preserve the pre-existing `srcInterface/FleksInterface.cpp` change.
