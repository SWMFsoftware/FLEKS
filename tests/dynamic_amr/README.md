# Dynamic AMR Multi-Session Regression Test

This test verifies the dynamic AMR lifecycle across multiple parameter sessions:

1. **Session 1 (Cycles 1 to 5):**
   - The grid starts with a two-level hierarchy where the left half of the domain (`x in [-1000, 0]`) is refined to level 1.
   - Initial conditions seed fluid and macroparticles.
   - The EM field solver and particle pusher advance for 5 cycles.

2. **Session 2 (Cycles 6 to 10):**
   - A subsequent session moves the refinement region to the right half of the domain (`x in [0, 1000]`).
   - `Domain::update_param()` updates the shapes, detects that refinement criteria changed, and marks `refineRegions` modified.
   - `Domain::update()` executes `Domain::regrid()`.
   - `FluidInterface::fill_new_cells()` interpolates coarse fluid state to newly refined nodes, ensuring valid non-zero fluid state.
   - Newly refined cells receive injected particles.
   - The field solver and particle pusher continue advancing for 5 more cycles without instability.

3. **Session 3 (Cycles 11 to 15):**
   - Refinement is cleared with `#REFINEREGION 0 none`.
   - The mesh is derefined back to a single level.
   - The field solver and particle pusher advance on the base mesh for 5 cycles.
