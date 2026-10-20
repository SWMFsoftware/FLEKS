# AMR Magnetic Reconnection Test

This is the AMR counterpart of the uniform reconnection test
(`tests/reconnection`): the same `FadeevIC` equilibrium, but on a **two-level
AMR grid** that refines the thin layer around the central reconnecting current
sheet, so the reconnection zone is resolved more finely than the surrounding
region at a fraction of the uniform-grid cost.

Two solver variants share this grid and setup — see `#FIELDSOLVER`:

- **Full PIC** (`PARAM.in.implicit`): Maxwell + implicit GMRES solve (`#FIELDSOLVER implicit`).
- **Hybrid PIC** (`PARAM.in.hybrid`): kinetic ions + fluid electrons (`#FIELDSOLVER hybrid`).

The runner auto-detects both variants, listing them as `RECONNECTION_AMR (IMPLICIT)` and `RECONNECTION_AMR (HYBRID)` in `tests/summary.md`.

Three short hybrid decks check the coarse-fine interface itself rather than
reconnection:

- **`PARAM.in.forcefree_hybrid`**: the force-free sheet of
  `../reconnection/PARAM.in.forcefree_hybrid` (`lambda = 1 d_i`, `Bg = 0.3`,
  uniform density) on a `10 x 5 d_i` box, base `32 x 16`
  (`dx = 0.3125 d_i`) with the patch `|x| < 2.5`, `|y| < 1.25 d_i` refined 2x,
  for 500 steps (half an ion gyroperiod). `|B|` and `n` are uniform in this
  equilibrium, so whatever the interface imprints on the field shows up.
- **`PARAM.in.gem_hybrid`**: a `1.6 x 3.2 d_i` Harris sheet with the GEM
  parameters (`lambda0 = 0.5 d_i`, `nb = 0.2`, conducting `y` walls), base
  `16 x 32` (`dx = 0.1 d_i`) with the patch `|x| < 0.4`, `|y| < 1.2 d_i`
  refined to `0.05 d_i`, for 30 steps. Its `y` interfaces lie in the
  low-density lobes, where injecting the fine E into the covered coarse nodes
  drove a `2 dx_coarse` mode along the interface (`max |B|` reached 3x B0
  within 26 steps here, and the full GEM AMR deck blew up).
- **`PARAM.in.lobe_hybrid`**: the GEM lobe plasma alone (`n = 0.2`,
  `B = B0 x`, `Ti/Te = 5`, no current sheet) in a periodic `3.2 x 3.2 d_i`
  box, base `32 x 32` (`dx = 0.1 d_i`) with the patch `|x|, |y| < 0.8 d_i`
  refined to `0.05 d_i`, for 50 steps. Its `x` interfaces are normal to `B`,
  so parallel whistlers cross them. With the coarse E on the interface nodes,
  the covered-cell relaxation closed a feedback loop that grew a mode ~6 fine
  cells long at ~40-70 `Omega_ci` (`max |B|` reached 2.3 B0 here, and the
  mode filled the refined patch of the full GEM AMR deck). The transverse
  field must stay at the particle noise level (`max |B_perp| < 0.1 B0`).

All three use a single block per level patch: with many small blocks the
coarse-fine fills dominate the run time.

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
- The refined layer (`|y| < 7.85 d_i`) covers the current sheet (`L = 5 d_i`)
  with margin; the ion inertial length `d_i ~ 1` is resolved at `0.98 d_i` on
  the fine level.

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
- **`#FIELDBOXBOUNDARY`**: `periodic` on the `x` faces, **`outflow`** (float B)
  on the `y` faces — the same BC as the uniform reconnection test.  The sheet
  (`L = 5 d_i`) is far from the `y`-walls (`Ly/2 = 15.7 d_i`), so the walls and
  the refined band (`|y| < 7.85 d_i`) do not interact.

## Initial condition

`FadeevIC` seeds the same force-free Fadeev equilibrium and m=1 perturbation as
the uniform reconnection test (`L = 5`, `eps = 0.4`, `num_islands = 2`,
`bg = 0`, `perturb = -0.05`, `nb = 0.2`), on both levels.

## Solver variants

| Parameter | Full-PIC | Hybrid |
|-----------|------------------------|----------------------------|
| `#SOLVEEM` | `T` | `F` |
| `#HYBRIDPIC` | `F` | `T` |
| `#PLASMA` | 2 species (ions + electrons, `m_i/m_e=25`) | 1 species (ions only; fluid electrons) |
| `#RESISTIVITY` | — | `4.0e5 m^2/s` (etaCode ~ 0.001) |
| `#ELECTRONTEMPERATURE` | — | `5 eV`, isothermal |
| `#HYPERRESISTIVITY` | — | grid mode, `C_h = 0.001` |
| `#BSUBCYCLE` | — | `4` |
| `#FIELDINTEGRATOR` | — | `rk4` |
| `#PARTICLES` | `4 4` (16 macroparticles/cell/species) | `5 5` (25 macroparticles/cell) |
| `#TIMESTEPPING dt` | `0.005` (semi-implicit; resolves electron timescale for `m_i/m_e=25`) | `0.02` (explicit; ion-scale step) |
| `#STOP TimeMax` | `3` | `20` (ion-scale onset is slower) |

## Validation

The automated check (`validate.py`) verifies the same criteria for both
variants, with solver-aware thresholds:

1. **Energy-log sanity**: finite, bounded `Eb` (no blow-up), finite `Epart`
   growing as magnetic energy is converted into kinetic energy.
2. **AMR refinement present**: coarse + fine cell spacings in the `z=0` `.out`
   plots, with the fine cells forming the central current-sheet layer.
3. **Active reconnection**: on the refined midplane `y ~ 0`, the equilibrium
   X-point (By null near `x ~ 0`) is present at `t=0`, and the seeded m=1
   perturbation grows (midplane `By` profile changes by `> 0.05`).
4. **Perturbation growth**: `|delta By|` grows past a solver-aware threshold —
   `0.1` for full-PIC, `0.05` for hybrid (the ion-scale onset is slower but the
   run is longer).

For every hybrid variant `validate_log` also reads the `divB-AMR` lines
(`#DIVB alwaysComputeDivB`) from the saved stdout. They must be there, each
with the level-0 and level-1 records and the `dGhost=...|B|=...` field (only
the full-PIC variant, which prints no such report, may skip the check): the fine-level interface
div(B) must stay within 10x the interior maximum, and the drift of the first
fine ghost layer from the coarse interpolation (`dGhost/|B|`) below 0.2 at the
first cycle and 0.5 afterwards. The ghost limit catches coarse-fine ghost
cells that are overwritten by the boundary fill (`dGhost/|B| ~ 1` from the
first cycle). For the force-free and GEM-like decks `max |B|` must in
addition stay within 1.5x its initial value, and the plot check is replaced by:

- **force-free**: `|B|^2` uniform over both levels (`std/mean < 0.05` at
  `t = 0`, `< 0.10` at the end) and the mean `|B|^2` of the refined patch within
  3% of the coarse level;
- **GEM-like**: the Harris sheet at `t = 0` (`Bx = +-1` at the walls, density
  peak ~1.2 and lobes ~0.2) and a finite field with `max |B| < 1.5 B0` at the end.

The hybrid test is **not** expected to reproduce the full-PIC reconnection rate
quantitatively (different electron physics), but it must show a stable, bounded
run with clear reconnection signatures.

## Running

This test needs the **2D AMReX** build with AMR support (`-amrex2d -lev=2`; see
`tests/README.md` for the build steps).  Then run (both variants together):

```bash
python3 tests/validate_tests.py --test=reconnection_amr
```

To run a **single variant**, append its token to the test name:

```bash
python3 tests/validate_tests.py --test=reconnection_amr.implicit          # full PIC only
python3 tests/validate_tests.py --test=reconnection_amr.hybrid            # hybrid only
python3 tests/validate_tests.py --test=reconnection_amr.forcefree_hybrid  # interface: force-free
python3 tests/validate_tests.py --test=reconnection_amr.gem_hybrid        # interface: GEM-like
```

## Long-horizon div(B) / ghost-drift measurements

A 100-time-unit run of the hybrid variant with the `#DIVB` diagnostic enabled
settled two questions about the coarse-fine treatment:

- every AMR div(B) bucket (level-0 interface / covered / interior and the
  fine-level interior) is **bit-for-bit flat** over 5000 cycles, so the
  interface adds no div(B) after the initial condition;
- the `dGhost` diagnostic (drift of the evolved first fine ghost layer from the
  coarse interpolation) stays bounded. (Measured during #160 it saturated near
  4x the peak |B|; with the current code it stays below 0.15 |B| over the
  20-unit hybrid run.)

It also shows that the level-1 *physical* domain-edge bucket grows steadily,
which is the outflow boundary condition rather than anything AMR does.
