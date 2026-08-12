# AMR Magnetic Reconnection Test

This test is the AMR counterpart of the uniform reconnection test
(`tests/reconnection`): the same `FadeevIC` equilibrium and Maxwell/GMRES
solver, but on a **two-level AMR grid** that refines the thin layer around the
central reconnecting current sheet, so the reconnection zone is resolved more
finely than the surrounding region at a fraction of the uniform-grid cost.

## Grid hierarchy

The grid is built against a **2D AMReX library** (`AMReX_SPACEDIM = 2`,
`./Config.pl -amrex2d -lev=2`); only `x` and `y` exist (`x` = periodic drive
direction, `y` = current-sheet normal).  Two levels:

| Level | Cells `(nx, ny)` | Resolution `(dx, dy)` |
|-------|------------------|-----------------------|
| 0 (base) | `32 x 16` | `1.9625 d_i` |
| 1 (fine) | `64 x 16` | `0.9813 d_i` |

- Box: `Lx = 2*pi*L*num_islands ~ 62.8 d_i`, `Ly = Lx/2 ~ 31.4 d_i`.
- The base grid matches the uniform reconnection test; level 1 refines the
  current-sheet layer to **2x resolution** (`#REFINEMENTRATIO 2`).
- `#MAXBLOCKSIZE 16 2` makes the refined layer align with whole base blocks so
  it refines cleanly.

## Refinement criteria

- **`#REGION`** `layer` = box `x in [-31.4, 31.4]`, `y in [-7.85, 7.85]`: the
  central current-sheet layer (full `x`, thin `y`).
- **`#REFINEREGION`** `iLev 0 +layer`: tags that layer for refinement at level 0
  (`iLev` must be `< max_level`, so with `nLevMax = 2` it tags level 0).
- **`#REFINEMENTRATIO 2`**: each tagged coarse cell becomes `2 x 2` fine cells.
- **`#GRIDEFFICIENCY 1.0`**: only base blocks lying **entirely** inside the
  layer refine.  With `maxBlockSizeY = 2` this yields a clean, block-aligned
  refined band `|y| < 7.85 d_i` (no partially-refined blocks).

## Boundary conditions

- **`#PERIODICITY`**: periodic in `x`, not in `y`.
- **`#BFIELDBOXBOUNDARY`**: `periodic` on the `x` faces, **`outflow`** (float B)
  on the `y` faces — the same BC as the uniform reconnection test.  The sheet
  (`L = 5 d_i`) is far from the `y`-walls (`Ly/2 = 15.7 d_i`), so the walls and
  the refined band (`|y| < 7.85 d_i`) do not interact.

## Initial condition

`FadeevIC` seeds the same force-free Fadeev equilibrium and m=1 perturbation as
the uniform reconnection test (`L = 5`, `eps = 0.4`, `num_islands = 2`,
`bg = 0`, `perturb = -0.05`, `nb = 0.2`), on both levels.

## Time stepping

Full-PIC, `m_i/m_e = 25`, fixed `dt = 0.005`, `TimeMax = 3` (600 steps), and
`#PARTICLES 5 5` (25 macroparticles/cell/species), matching the uniform test.

## Running

Build with the 2D AMReX library and AMR support, then run the check:

```bash
./Config.pl -amrex2d -lev=2 && make -j4
python3 tests/validate_tests.py --test=reconnection_amr
```

## Validation

The automated check (`validate.py`) verifies:

1. **Energy-log sanity**: finite, bounded `Eb` (no blow-up), with `Epart`
   growing as magnetic energy is converted into kinetic energy — the signature
   of reconnection.
2. **AMR refinement present**: the `z=0` `.out` plots contain both the coarse
   and the fine cell spacings, with the fine cells forming the central
   current-sheet layer — confirming the layer is actually refined.
3. **Active reconnection**: on the refined midplane `y ~ 0`, the equilibrium
   X-point (By null near `x ~ 0`) is present and the seeded m=1 perturbation
   grows (the midplane `By` profile changes).
