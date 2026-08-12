# AMR Magnetic Reconnection Test

This test demonstrates **adaptive mesh refinement (AMR)** applied to the Fadeev
current-sheet reconnection problem.  It is the AMR counterpart of the uniform
full-PIC `reconnection` test (`tests/reconnection/PARAM.in`), using the same
`FadeevIC` equilibrium, the same full-PIC Maxwell/GMRES field solver, and the
same physical parameters — but with a **two-level AMR grid** that refines the
thin layer around the central reconnecting current sheet via `#REFINEREGION`.

The physics goal is to **resolve the reconnection zone (the central current
sheet and the growing islands) at the same resolution as the validated
uniform-grid test while leaving the surrounding region coarse**, so the
simulation captures the dynamics of magnetic reconnection accurately at a
fraction of the uniform-grid cost.

## Grid hierarchy

The test is built against a **true-2D AMReX library** (`AMReX_SPACEDIM = 2`,
`./Config.pl -amrex2d -lev=2`), so the grid has only `x` and `y` (the
reconnection plane).  `x` is the periodic drive direction, `y` the
current-sheet normal (the Fadeev profile coordinate).  The AMR grid has two
levels:

| Level | Cells `(nx, ny)` | Resolution `(dx, dy)` | Coverage |
|-------|------------------|-----------------------|----------|
| 0 (base) | `32 x 16` | `1.9625 d_i` (code units) | whole box |
| 1 (fine) | `64 x 16` | `0.9813 d_i` (code units) | refined layer only |

- Box: `Lx = 2*pi*L*num_islands ~ 62.8 d_i`, `Ly = Lx/2 ~ 31.4 d_i`.
- **Base grid** (`#NCELL 32 16`): `dx = dy = 62.8/32 = 1.9625 d_i` — the same
  grid as the uniform-grid full-PIC reconnection test (`32 x 16`).
- **Refined level-1 grid** (`#REFINEMENTRATIO 2`): `dx = dy = 0.9813 d_i`
  inside the layer, so the reconnection zone is resolved **twice as finely**
  as the uniform test — enough to resolve the current-sheet scale `L = 5 d_i`
  and the reconnecting islands accurately.
- **`#MAXBLOCKSIZE 16 2`**: the base grid is decomposed into `2 x 8` blocks of
  `16 x 2`.  `maxBlockSizeY = 2` is chosen so the refined layer aligns with
  whole base blocks (see *Refinement criteria* below).

## Refinement criteria

- **`#REFINEREGION`**: the box region `layer`
  (`x in [-31.4, 31.4]`, `y in [-7.85, 7.85]`) is tagged for refinement at
  level 0 (`iLev = 0`), i.e. base cells whose **centers** fall inside the layer
  are tagged.  The layer spans the full `x` extent (the whole current sheet is
  the reconnecting region in `x`) and is **thin in `y`** (the sheet normal) so
  it is a genuine "layer near the central reconnecting current sheet".
  `iLev = 0` (rather than 1) is required: with `nLevMax = 2` the `#REFINEREGION`
  index must be `< max_level`, and tagging level-0 cells creates the single
  refined level-1 grid.
- **`#REFINEMENTRATIO 2`**: each tagged coarse cell is replaced by `2 x 2`
  fine cells (ratio 2), giving `64 x 16` fine resolution inside the layer.
- **`#GRIDEFFICIENCY 1.0`**: a base block is refined only if **all** of its
  cells are tagged (AMReX's `gridEfficiency` threshold is the fraction of a
  block that must be tagged for the whole block to be refined).  With
  `maxBlockSizeY = 2`, the base `y`-blocks are `[0,1]`, `[2,3]`, `[4,5]`,
  `[6,7]`, `[8,9]`, `[10,11]`, `[12,13]`, `[14,15]`; the layer covers cells
  `j = 4..11` (`|y| < 7.85 d_i`), so blocks `[4,5]`..`[10,11]` (and both
  `x`-blocks) refine, and the boundary blocks stay coarse.  This yields a
  **clean, block-aligned** refined band `y in [-7.85, 7.85] d_i` — no
  partially-refined blocks.  (A lower `gridEfficiency`, e.g. the AMReX default
  0.7, would also work but can produce ragged fine/coarse boundaries.)
- **`nLevMax = 2`** (from the build flag `-lev=2`): this test, like `lightwave`,
  requires AMR support.  With `nLevMax = 1` the `#REFINEREGION` command aborts
  at startup (`iLev` must be smaller than the max level index).

## Boundary conditions

- **`#PERIODICITY`**: periodic in `x`; **not** periodic in `y` (the
  cross-sheet / current-sheet-normal direction).
- **`#BFIELDBOXBOUNDARY`**: `periodic` on the `x` (Lo/Hi) faces and
  **`outflow`** (float B) on the `y` (Lo/Hi) faces.  This is the same BC as the
  uniform reconnection test.  The current sheet is far from the `y`-walls
  (`Ly/2 = 15.7 d_i` vs `L = 5 d_i`), so the walls do not affect the interior
  reconnection, and the AMR refinement region does not touch the `y`-faces
  (the refined band `|y| < 7.85 d_i` is well inside the box).

## Initial condition

`FadeevIC` (`src/ic/FadeevIC.h/.cpp`) seeds the Fadeev force-free equilibrium
and the m=1 perturbation that drives reconnection, exactly as in the uniform
reconnection test.  With `#TESTCASE fadeev`, `L = 5`, `eps = 0.4`,
`num_islands = 2`, `bg = 0`, `perturb = -0.05`, `nb = 0.2`:

```
Bx = b0*sinh(y/L)/(cosh(y/L) + eps*cos(x/L))   + m=1 perturbation
By = b0*eps*sin(x/L)/(cosh(y/L) + eps*cos(x/L)) + m=1 perturbation
Bz = 0
n(x,y) = nb + (1-nb)*profile(x,y)/profile_max
```

The field and density profiles are set on both levels by `FadeevIC::set_fields`
and the particle weights are boosted to the sheet profile on every level.

## Time stepping

- Full-PIC solver, `m_i/m_e = 25`, fixed `dt = 0.005` (resolves the electron
  timescale), `TimeMax = 3` (600 steps), matching the uniform full-PIC
  reconnection test.
- `#PARTICLES 5 5` = 25 macroparticles per cell per species.  The refined layer
  carries finer cells (2x in each of `x` and `y`), so more macroparticles are
  resolved inside the reconnection zone while the outer region stays coarse;
  the total particle count keeps the runtime within the standalone serial
  budget.

## AMR control: standalone PC vs coupled GM-PC

This standalone test refines the grid with PC's native AMR commands, which are
the direct counterpart of the static-geometric refinement used in the coupled
GM-PC reconnection runs (e.g. `Param/PARAM.in.test.GMPC.FLEKS.AMR.FASTWAVE.2D`).

| GM (BATSRUS) AMR control | PC (FLEKS) equivalent | Used here |
|--------------------------|-----------------------|-----------|
| `#GRID` / `#PICGRID` (base grid + cell size) | `#NCELL` + `#MAXBLOCKSIZE` | `32 x 16`, blocks `16 x 2` |
| `#AMRREGION` (box/sphere/shell/paraboloid) | `#REGION` (same shapes) | `layer` = box `\|y\| < 7.85 d_i` |
| `#AMRCRITERIALEVEL` `RefineTo`/`CoursenTo` | `#REFINEREGION` (`refineLev` + region) | `iLev 0` + `+layer` |
| `#GRIDRESOLUTION` per area | `#REFINEMENTRATIO` | `2` |
| block-efficiency threshold | `#GRIDEFFICIENCY` | `1.0` (clean block-aligned band) |
| — | `#CONSTANTPPV` | `F` (per-cell particle count) |
| — | `#RESAMPLING` | particle splitting/merging at coarse-fine interface |

GM additionally supports *dynamic* AMR (`#DOAMR` + `#AMRCRITERIA` curlB/gradlogP/
J2, and the adaptive-PIC `#PICCRITERIA` `j/bperp` + `#PICADAPT`) that grows or
shrinks the refined/PIC region as the physics evolves.  PC's `#REGION` /
`#REFINEREGION` is static-geometric only; it cannot adapt to a moving criterion.
For this reconnection problem the central current sheet is fixed at `y = 0`, so
a static refined layer near `y = 0` is the correct and sufficient control —
exactly the choice made by the coupled GM-PC fast-wave AMR test on the PC side.

## Running

This test is built and run against a **true-2D AMReX library**
(`AMReX_SPACEDIM = 2`), which avoids the coarse/fine z-inconsistency that makes
fake-2D (nCellZ = 1 backed by 3D AMReX) multi-level runs produce NaN in the
field energy.  The binary needs AMR support (`nLevMax >= 2`) as well:

```bash
./Config.pl -amrex2d -lev=2 && make -j4
```

Then run the automated check:

```bash
python3 tests/validate_tests.py --test=reconnection_amr
```

## Validation

The automated check (`validate.py`) verifies:

1. **Energy-log sanity** (from `log_pic`): finite `Eb` / `Epart`, bounded `Eb`
   (no blow-up).  Because magnetic reconnection converts stored magnetic energy
   into ion kinetic energy, `Epart` is expected to grow while `Eb` stays
   bounded — the signature of reconnection rather than a numerical instability.
2. **AMR refinement present** (from the `z=0` x-y plane `.out` files): the plot
   output contains **both** the coarse (`dx ~ 1.9625 d_i`) and the fine
   (`dx ~ 0.98 d_i`) cell spacings, with the fine cells forming the central
   current-sheet layer `|y| < 7.85 d_i` — confirming the layer is actually
   refined.
3. **Active reconnection** (from the `.out` files, along the refined midplane
   `y ~ 0`): the equilibrium carries the X-point (By null near `x ~ 0`), and the
   seeded m=1 perturbation grows (the midplane `By` profile changes) as the
   instability drives reconnection inside the refined layer.

The `#SAVEPLOT` writes the two-level grid per-cell (`dx = -1`, non-uniform), so
the validator reads the `.out` files line-by-line (using the column names from
the header) instead of reshaping them into one rectangular grid.
