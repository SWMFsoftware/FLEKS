# Dynamic AMR Multi-Session Regression Test

This test verifies the dynamic AMR lifecycle across multiple parameter sessions:

1. **Session 1 (Cycles 1 to 5):**
   - The grid starts with a two-level hierarchy where the left half of the domain (`x in [-1000, 0]`) is refined to level 1.
   - Initial conditions seed fluid and macroparticles.
   - The EM field solver and particle pusher advance for 5 cycles.

2. **Session 2 (Cycles 6 to 10):**
   - A subsequent session defines a new `right` shape and changes `#REFINEREGION` from `+left` to `+right`, moving refinement to `x in [0, 1000]`.
   - `Domain::update_param()` detects the changed selector and marks `refineRegions` modified.
   - `Domain::update()` executes `Domain::regrid()`.
   - `FluidInterface::fill_new_cells()` interpolates coarse fluid state to newly refined nodes, ensuring valid non-zero fluid state.
   - Newly refined cells receive injected particles.
   - The field solver and particle pusher continue advancing for 5 more cycles without instability.

3. **Session 3 (Cycles 11 to 15):**
   - Refinement is cleared with `#REFINEREGION 0 none`.
   - The mesh is derefined back to a single level.
   - The field solver and particle pusher advance on the base mesh for 5 cycles.

Named `#REGION` shapes are immutable. A later session cannot redefine `left`,
even with the same coordinates; use a new name and update `#REFINEREGION`.

Run the focused checks with:

```bash
python3 -m unittest tests.dynamic_amr.test_region_names -v
python3 tests/validate_tests.py --test=dynamic_amr
```

The first check also verifies that defining a new shape alone does not regrid.
