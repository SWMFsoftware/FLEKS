# Semi-Implicit AMR Inflow + Div(E) + Co-Moving Boundary Tests (`tests/amr_dive_inflow`)

Regression suite for semi-implicit PIC AMR boundary, `div(E)`, `#BODY`, and
co-moving frame synchronization when `#PRESPLITTING`, `#RESAMPLING`, `#DIVE`,
and `#COMOVING` are active together:

1. **Inflow face `target_PPC` injection scaling (`src/ParticlesBC.cpp`)**:
   Injects inflow particles using per-cell `target_PPC(0)` (64 PPC inside the
   2-cell pre-splitting collar when base PPC is 8) with matching macroparticle
   weights so the inflow boundary does not drain or inject overweight particles.
2. **Unified `div(E)` boundary, body, & coarse-fine stencil guard (`include/GridUtility.h`, `src/ParticlesMover.cpp`, `src/ParticlesMoments.cpp`, `src/PicDivE.cpp`, `src/AmrexLinearSolver.cpp`)**:
   - Keeps `L0` physical `is_lev_edge` cells active in `div(E)` cleaning on
     non-coupled boundaries while preserving `is_lev_edge` skipping on
     MHD-coupled (`FieldBC::coupled` / `ParticleBC::coupled`) boundaries.
   - Uses `dive_particle_cleaning_weight()` and
     `skip_particle_for_dive_cleaning()` across `sum_to_center_amr()` and
     `divE_correct_position()`: `0.25 * dx` coarse (`is_refined_neighbour`)
     margin excluded from both `sum_to_center_amr()` and
     `divE_correct_position()` (retaining the Dirichlet anchor across refined
     inner patches), and `0.5 * dx` fine (`iLev > 0` non-domain `is_lev_edge`)
     half-cell margin skipped in `divE_correct_position()` while contributing a
     `0.35` Robin/Dirichlet anchor weight to `centerMM` in
     `sum_to_center_amr()` so active fine particles never differentiate across
     unsolved `centerPhi == 0` ghost cells yet interfaces remain locally
     grounded.
   - Pre-shift dual-cell guarding skips only `bit::is_domain_boundary` cells
     (when `iLev == 0 || !bc.has(ParticleBC::coupled)`), keeping interior
     coarse-fine edges and `#BODY`-adjacent cells active with a post-shift
     revert whenever a shift would cross into domain/level ghost cells or
     inside `#BODY`.
   - Unconditionally refreshes `AmrexDivESolver::ensure_mask()` contents so
     dynamic regrids on `iLev + 1` update `bit::is_refined` masking on `iLev`.
3. **Startup / regrid clearance warning (`src/Pic.cpp`)**:
   Emits `WARNING [Pic::SetTargetPPC]` when `L1` is placed fewer than 4 `L0`
   cells from a `ParticleBC::inflow` domain face (exercised by `PARAM.in` with
   2 `L0` cells, and clean in `PARAM.in.comoving_b0` with 4 `L0` cells).
4. **Co-moving frame `uBg`/`eBg` AMR sync & intrinsic `B0` (`src/Pic.cpp`)**:
   Pre-fills covered coarse `uBg` before smoothing, refreshes coarse-fine ghosts
   between `smooth_multifab` passes, restricts smoothed `uBg`/`eBg`
   fine-to-coarse onto covered coarse nodes (`bit::is_refined`), includes
   intrinsic `nodeB0` in `eBg = -uBg x (B1 + B0)` consistently with the
   particle mass matrix, and syncs coarse-fine ghosts without overwriting
   fine-interface midpoint `eBg`.

## Variants

- **`PARAM.in` (`amr_dive_inflow.full`)**:
  - Domain `x in [-2.0, 1.5]`, `y, z in [-2.0, 2.0]` (`14 x 16 x 16`, `dx_0 = 0.25`,
    `y, z` periodic), `L1` `[-1.0, 1.0]^3` (only `2` `L0` cells upstream of `L1`),
    uniform field-aligned `B = (-2, 0, 0) nT`, `80` steps (`TimeMax = 8.0`).
- **`PARAM.in.comoving_b0` (`amr_dive_inflow.comoving_b0`)**:
  - Non-periodic 3D domain `[-2.0, 2.0]^3` (`16^3`, `dx_0 = 0.25`), `L1`
    `[-1.0, 1.0]^3` (`4` `L0` cells upstream of `L1`), `ambipolar`/`linetied`
    spherical `#BODY` (`r = 0.45`), `#DIPOLE` `30 nT`, solar wind `n = 2 /cc`,
    `u = (-300, 0, 0) km/s`, `By = 1.5 nT` (`Ez = +450 uV/m`), `50` steps
    (`TimeMax = 5.0`).

## Running

```bash
python3 tests/validate_tests.py --test=amr_dive_inflow -n 4 --omp-threads 2 -v
```
