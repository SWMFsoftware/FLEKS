# Plan: Direct Cell-Centered Output for the Hybrid-PIC Solver

**Goal**: eliminate the cell → node conversion entirely on the hybrid-PIC
(`useHybridPIC`) output path, so plot / restart / diagnostics read the live
cell-centered fields (`centerB`, `centerEhybrid`, `centerPlasma`,
`centerPlasmaSum`, `centerJ`) directly, and the node mirrors (`nodeE`,
`nodeB`, `nodePlasma`) are no longer materialized for output.

**Non-goal**: the full-PIC path is untouched and must stay bit-identical
(guarded by the `tests/beam` energy diagnostics).

## Status

- **Phase 0**: DONE. `get_var` exposes the hybrid total current
  (`centerJ`) for `jxS/jyS/jzS` (was silently 0); `ParticleTracker` now
  syncs `nodeE` for hybrid (was stale).
- **Phase 1**: DONE. Structured (IDL/ascii) output for hybrid reads the
  live cell-centered fields directly; the per-step `nodeB` projection and
  the `nodePlasma`/`nodeE` output mirrors are no longer materialized for
  hybrid structured plots. Load balancing and Mach number are computed
  from `centerPlasmaSum` directly. All 18 standalone tests pass
  (including the `beam` full-PIC bit-identity guard).

  Implementation notes / discoveries:
  - `cellStatus` never has the owner bit set (only `nodeStatus` does, in
    `update_node_status`); cells are uniquely owned by construction, so the
    hybrid branch checks only refinement, not ownership.
  - `PlotWriter::nDim` (from `get_fluid_dimension()`) can be **less than
    `Pic::nDim`** for effectively-2D problems (e.g. `z=0` cut with nCellZ=1).
    `PlotWriter::init` only sets plot bounds for the writer's `nDim`; the
    extra dimension keeps sentinel `plotMin/Max` (`1`/`-1`). `get_var`'s
    cell-center coordinates and `is_inside_cell_plot_region` must therefore
    use the writer's `nDim` (added `PlotWriter::get_nDim()`), not `Pic::nDim`,
    or every point is rejected.
  - `mMach` is allocated on the **cell** grid for hybrid (was node grid) so
    `get_var` (which iterates `centerB`) can read it; `calc_mach_number`
    reads `centerPlasmaSum[nSpecies]` for hybrid.

---

## 1. Current Status

The hybrid solver core is already purely cell-centered:

- The RK4 / SSPRK3 sub-cycle loop in `Pic::update_B_hybrid`
  (`src/Pic.cpp:2666-2831`) reads/writes only `centerB`, `centerEhybrid`,
  `centerJ`, `centerPlasma*` — no center↔node round-trips inside the loop.
- The particle mover gathers from cell-centered fields
  (`charged_particle_mover_cell_centered`, `src/Particles.cpp:1791+`).
- Moments are deposited cell-centered into `centerPlasma[iSpecies]`
  (`Pic::sum_moments`, `src/Pic.cpp:1171-1199`).

Paths that are **already** cell-centered for hybrid (no conversion):

| Path | Location | Note |
|---|---|---|
| amrex / hdf5 plotfile | `Pic::write_amrex_field`, `src/PicIO.cpp:819-1101` | Direct `MultiFab::Copy` from `centerB` / `centerEhybrid` / `centerPlasma` (`saveNode=false` default) |
| Restart write / read | `src/PicIO.cpp:454-469, 512-529` | Writes `centerEhybrid` + `centerB`; node fields intentionally not saved |
| GM coupling | `Pic::get_fluid_state_for_points`, `src/PicIO.cpp:61-91` | Reads cell-centered fields via `get_value_at_loc` |

Residual cell → node conversions (the subject of this refactor):

1. **Structured IDL/ascii plots** — `Pic::write_plots` materializes the node
   mirrors before every ascii/IDL write (`src/PicIO.cpp:658-669`):
   - `sync_node_plasma_output(needMach)` — `average_center_to_node` for
     `centerPlasma → nodePlasma` over `nSpecies+1` MultiFabs × `nMoments`(11)
     components, each with `FillBoundary` MPI exchanges
     (`src/Pic.cpp:1230-1249`).
   - `sync_node_E_output()` — `average_center_to_node` for
     `centerEhybrid → nodeE` (`src/Pic.cpp:1252-1265`).
   - The whole structured path (`find_output_list`, `get_var`,
     `src/PicIO.cpp:99-444`) iterates the node grid (`nodeE` MFIter,
     `nodeStatus`, `LoEdge` coordinates).
2. **Per-step `centerB → nodeB` projection** at the end of
   `Pic::update_B_hybrid` (`src/Pic.cpp:2839-2868`). Not needed by the solver;
   consumed only by `dBdt` diagnostics, the test-particle tracker, and IDL
   `Bx/By/Bz` output.
3. **Test-particle tracking** — `ParticleTracker` copies `pic.nodeE` /
   `pic.nodeB` (`src/ParticleTracker.cpp:33, 105, 132-135`).
4. **Load balancing** — `Pic::calc_cost_per_cell` calls
   `sync_node_plasma_output(false)` then
   `average_node_to_cellcenter(nodePlasma[nSpecies] → cellCost)`
   (`src/Pic.cpp:1329-1346`): a cell→node→cell double conversion.

---

## 2. Bottleneck

Cost per IDL/ascii plot event (hybrid):

- `(nSpecies+1) × 11 + 3` components of `average_center_to_node` averaging
- Multiple `FillBoundary` rounds (MPI halo exchanges) per MultiFab
- Plus, on every step regardless of output: one 3-component
  `centerB → nodeB` projection + `FillBoundary` + boundary-condition fill
  (`src/Pic.cpp:2852-2868`).

The impact scales with output cadence: minor for sparse production output,
significant for output-dense hybrid regression tests and debugging runs.

Secondary (correctness) issues with the current mirror-based output:

- **`{all}` variables silently zero**: hybrid `{all}` expansion adds
  `qS kXXS ... jxS jyS jzS` (`src/PlotWriter.cpp:289-294`), but `get_var`
  has no branch for them → they are written as 0 (`src/PicIO.cpp:438-439`).
- **Cell quantities indexed as nodes**: `qc`, `divEc`, `phi` in `get_var`
  read cell-centered MultiFabs with node indices (`src/PicIO.cpp:416-427`),
  a half-cell misalignment.
- **Stale nodeE for test particles**: `ParticleTracker` never calls
  `sync_node_E_output()`, so hybrid test-particle runs track a stale E
  mirror (B is fine — projected every step).
- The node averaging is a smoothing operation: mirrored output has lower
  effective resolution than the live cell data.

---

## 3. Proposed Changes

### Phase 0 — Independent bug fixes (no output-format change)

1. **`get_var` hybrid branches for `jxS/jyS/jzS`, `qS`, `kXXS` …** — expose
   `centerJ` (and document that `kXXS` tensor / net-charge output is
   node-only for full-PIC) so `{all}` hybrid output is not silently zero.
2. **`ParticleTracker` E sync** — call `pic.sync_node_E_output()` (or read
   `centerEhybrid` directly) before tracking so hybrid test particles see
   the live field.

### Phase 1 — Cell-centered structured (IDL/ascii) output for hybrid

Gated entirely on `useHybridPIC`; full-PIC code paths unchanged.

1. **`Pic::find_output_list`** — hybrid branch iterating `cGrids[iLev]` with
   `cellStatus`, coordinates from `Geom(iLev).CellCenter()` instead of
   `LoEdge`. The periodic "drop the rightmost node" special case disappears
   (cell centers have no duplicated periodic point).
2. **`Pic::get_var`** — hybrid branches reading `centerEhybrid` (Ex/Ey/Ez),
   `centerB` (Bx/By/Bz), `centerPlasma` (rhoS/uxS/.../ppcS), `centerJ`
   (jxS/jyS/jzS). Node-only diagnostics (`jHat*`, `nMM*`, `dB*dt`, `E0*`,
   `u0*`) either stay node-based (full-PIC only) or are documented as
   unavailable for hybrid.
3. **`Pic::write_plots`** — skip `sync_node_plasma_output` /
   `sync_node_E_output` on the hybrid path.
4. **`Pic::calc_mach_number`** — hybrid branch reading
   `centerPlasmaSum[nSpecies]` directly (removes the nodePlasma dependency
   for the `mach` plot variable).
5. **`Pic::calc_cost_per_cell`** — hybrid branch computing `cellCost` from
   `centerPlasmaSum[nSpecies]` `iNum_` directly; delete the
   cell→node→cell double conversion.
6. **Cut planes (`x=`/`y=`/`z=`)** — snap the requested cut coordinate to the
   nearest cell-center plane. Today `PlotWriter::init` sets
   `plotMax = plotMin + 1e-10` and `is_inside_plot_region` uses a 0.01·dx
   tolerance, which matches node planes exactly but would select **zero**
   cell rows once coordinates shift by dx/2. Implement snapping in the
   hybrid branch of `find_output_list` (round the cut coordinate to the
   nearest `CellCenter` index before the region test), and widen the
   tolerance to 0.5·dx for hybrid.
7. **Header metadata** — `#PLOTRANGE` correction (±0.4·dx in
   `PlotWriter::write_idl`) happens to yield the correct cell count for
   cell centers as well; verify the reported min/max and `#CELLSIZE` remain
   consistent for PostIDL.

### Phase 2 — Lazy per-step `nodeB` projection

- Move the `centerB → nodeB` block at `src/Pic.cpp:2839-2868` behind a
  `nodeBStale` flag (mirroring the existing `nodePlasmaStale` / `nodeEStale`
  pattern) and a `sync_node_B_output()` helper.
- Materialize only when a consumer needs it: IDL plot requesting
  `B*` / `dB*dt`, test-particle tracker active, or `useHyperbolicCleaning` /
  `useUpwindB` paths that read `divB`/`nodeB`.

### Phase 3 — Restart / reference data / tests

- Restart is already cell-native for hybrid; no change.
- Regenerate hybrid test reference outputs (`.out` files): coordinates shift
  by dx/2 per axis. Note the point count is unchanged under periodic BCs
  (the current node output drops the duplicated periodic node, so it already
  emits `nCell` points), and `tests/_shared/hybrid.py` validators are
  index-based (spatial DFT, phase-fit dispersion) and insensitive to the
  dx/2 shift — logic changes should not be needed, but all reference data
  must be regenerated.
- Run the full hybrid suite (`beam`, `pcai`, `ohm`, `whistler`, `iaw`,
  `freestream`, `singlecell`, `zerocurrent`) plus the full-PIC `beam`
  bit-identity check.

### Deferred / out of scope

- Test-particle tracker on cell-centered fields natively (Phase 2 keeps the
  nodeB/nodeE mirror only when the tracker is active; a cell-centered
  tracker gather is a separate change).
- Multi-level structured output (`plotDx >= 0` already aborts for
  `n_lev() > 1`; cell-centered output inherits this limitation).

---

## 4. Risks and Mitigations

| Risk | Mitigation |
|---|---|
| Cut-plane semantics change (empty cuts) | Explicit snap-to-nearest-cell-center with 0.5·dx tolerance; document in `PARAM.XML` |
| Reference-data churn for all hybrid tests | Single regeneration commit; validators are shift-insensitive |
| Full-PIC regression | All changes gated on `useHybridPIC`; `tests/beam` bit-identity guard |
| Consumers expecting node grids (`nodeStatus`-based ownership) | Hybrid branch uses `cellStatus` consistently; verify `bit::is_owner` semantics on cell status |
| PostIDL compatibility of header block | Verify `#PLOTRANGE`/`#NCELL`/`#CELLSIZE` against a 1-cell and multi-cell periodic case |

---

## 5. Expected Outcome

- Zero cell→node conversions on the hybrid plot/restart/diagnostic path
  (except when the test-particle tracker is active, Phase 2).
- Higher-fidelity output (no averaging-induced smoothing of the mirrors).
- Fixes the silent-zero `{all}` current output and the half-cell
  misalignment of `qc`/`divEc`/`phi` as side effects.
- Lower per-plot cost: removes `(nSpecies+1)×11 + 3` components of
  averaging plus halo exchanges per IDL write, and one 3-component
  projection per step when no output/tracker needs `nodeB`.
