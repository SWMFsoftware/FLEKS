# Semi-Implicit AMR Inflow + Div(E) Boundary Test (`tests/amr_dive_inflow`)

Regression test for three coupled inflow/AMR defects in the semi-implicit PIC
solver when `#PRESPLITTING`, `#RESAMPLING`, and `#DIVE` are active together:

1. **Inflow face `target_PPC` injection scaling (`src/ParticlesBC.cpp`)**:
   Injects inflow particles using per-cell `target_PPC(0)` (64 PPC inside the
   2-cell pre-splitting collar when base PPC is 8) with matching macroparticle
   weights so the inflow boundary does not drain or inject overweight particles.
2. **Unified `div(E)` boundary & coarse-fine stencil guard (`include/GridUtility.h`, `src/ParticlesMover.cpp`, `src/ParticlesMoments.cpp`, `src/PicDivE.cpp`, `src/AmrexLinearSolver.cpp`)**:
   Keeps `L0` physical `is_lev_edge` cells active in `div(E)` cleaning on
   non-coupled boundaries while preserving the `is_lev_edge` skip on
   MHD-coupled (`FieldBC::coupled` / `ParticleBC::coupled`) boundaries, and
   applies the `2^D` `is_lev_boundary` / body particle-cleaning guard across all
   levels `iLev >= 0` for non-coupled and `#BODY` runs so neither the outer
   inflow face nor coarse-fine interfaces accumulate charge-separation spikes.
3. **Startup / regrid clearance warning (`src/Pic.cpp`)**:
   Emits `WARNING [Pic::SetTargetPPC]` when `L1` is placed fewer than 4 `L0`
   cells from a `ParticleBC::inflow` domain face (here 2 `L0` cells).

## Physical & Grid Setup

- **Domain**: `x in [-2.0, 1.5]`, `y, z in [-2.0, 2.0]` (`14 x 16 x 16` base mesh,
  `dx_0 = 0.25`, `maxBlockSize = 8`).
- **Boundaries**: `x-` `outflow`/`open`, `x+` `inflow`/`inflow`, `y, z` `periodic`.
- **Refinement**: Stationary `L1` box `[-1.0, 1.0]^3` (`dx_1 = 0.125`), leaving
  only `2` `L0` cells `[1.0, 1.5]` upstream of `L1`.
- **Plasma**: Uniform magnetized electron (`m = 0.0625`) + proton (`m = 1.0`)
  flow at `n = 1.0 /cc`, `u = (-300, 0, 0) km/s`, `T = 60,000 K`,
  field-aligned `B = (-2.0, 0, 0) nT` (`E = -u x B = 0`).
- **Duration**: `dt = 0.1`, `TimeMax = 8.0` (`80` steps).

## Running

```bash
python3 tests/validate_tests.py --test=amr_dive_inflow -n 4 --omp-threads 2 -v
```
