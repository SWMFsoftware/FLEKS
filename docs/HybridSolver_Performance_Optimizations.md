# Hybrid Solver: Performance Optimization Opportunities

**Status:** Items 1, 2, 3, 4, 7 implemented; Item 8 partially implemented
(8-C/D/E done; 8-A/B reverted; 8-sync nodeEth bridge removed); Items 5, 6
decided — deferred (need output decoupling, see Item 4); 8-sync `nodeB`+`dBdt`
deferred (see §8.6)
**Authors:** FLEKS team
**Last updated:** 2026-08-05
**Related docs:**
- [`HybridVPIC_Style_Solver_Implementation_Plan.md`](./HybridVPIC_Style_Solver_Implementation_Plan.md)
- [`HybridSolver_Differences_FLEKS_vs_HybridVPIC.md`](./HybridSolver_Differences_FLEKS_vs_HybridVPIC.md)

## 1. Scope and Motivation

The hybrid (kinetic ion + fluid electron) solver was implemented as described in
Phase A of the implementation plan: fully **cell-centered** field storage with
trilinear gather/scatter. Phase B (quadratic gather + Esirkepov deposit) is
present in code but gated off (`useQuadraticGather = false`,
`useEsirkepovDeposit = false`, `include/Pic.h:314,321`).

This document lists concrete performance optimization opportunities found by
reviewing the current hybrid solver hot paths. Every item below is **gated
behind `useHybridPIC`** and is intended to leave the full-PIC solver
bit-identical. The items are ranked roughly by expected impact.

> **Note on the sub-cycle multiplier.** `update_B_hybrid` wraps its integrator
> in a `for (subStep = 0; subStep < nHallSubcycle; ++subStep)` loop
> (`src/Pic.cpp:2805`). For the current PCAI config (`nHallSubcycle = 32`),
> everything inside the loop runs 32× per step. The target config
> (`nHallSubcycle = 1`) removes most of that multiplier, but several items
> below still hold at `nHallSubcycle = 1`.

---

## 2. Itemized Opportunities

### Item 1 — Remove `euler` / `heun` integrators; keep `rk3` and `rk4` with `rk4` as default

- **Location:** `update_B_hybrid`, `src/Pic.cpp:2754-3143`
  - euler branch: `src/Pic.cpp:2902-2927`
  - heun branch: `src/Pic.cpp:2989-3038`
  - rk3 branch: `src/Pic.cpp:2929-2987`
  - rk4 branch: `src/Pic.cpp:2815-2900`
  - integrator dispatch: `useRK4` flag and `fieldIntegrator` string
    (`src/Pic.cpp:2815, 2902, 2929`)
- **Impact:** Medium–High (simplifies the hot sub-cycle loop; frees 2 scratch
  MultiFabs; removes 2 of 4 integrator branches)
- **Regression risk:** Low (rk4 is the current default; rk3 is kept)
- **Status:** implemented (2026-08-05)

**Decision.** Remove the `euler` and `heun` field-integrator options from
`update_B_hybrid`, keeping only `rk4` (default) and `rk3`. `J = curl(B)/(4π)`
is recomputed inside every `assemble_ohm_E` call on the current trial B; in the
`rk3` / `rk4` stage loops the trial B differs at every stage, so `J` **must** be
recomputed per stage — there is no cross-stage `J` reuse to exploit. (The
`euler` / `heun` paths offered no reuse either: euler has one call per sub-step
on an advanced B, and heun's predictor/corrector run on different B states.)
The value of this item is therefore **not** avoiding a redundant curl, but:

- **Fewer per-sub-cycle Ohm + curl stage passes**:
  - `rk4`: 4, `rk3`: **3**, `heun`: 2, `euler`: 1.
  - At the target `nHallSubcycle = 1`, `rk3` is already **25% cheaper** than
    `rk4` on the field side.
- **Free the euler/heun-only scratch MultiFabs** `dBpred_heun`, `dBcorr_heun`
  (used only by the removed euler/heun branches). **Correction vs. the earlier
  draft:** `centerBstar_heun` and `centerBstart_heun` are **not** heun-only —
  `centerBstar_heun` holds the time-centred `(trial + B^n)/2` state used by the
  rk3/rk4 time-centred-E stages, and `centerBstart_heun` holds the sub-step start
  `B_n` used by ssprk3. So only `dBpred_heun` / `dBcorr_heun` were freed; the
  other two were kept and re-commented as rk3/rk4 scratch.
  Note: the implementation plan (Decision #2, Phase C.1) explicitly chose to
  **keep** the heun scratch arrays; removing the euler/heun integrators reverses
  that decision for the now-unused `dBpred_heun` / `dBcorr_heun`.
- **Fewer branches** in the hot `nHallSubcycle` loop → simpler, less code to
  keep regression-verified.

**Implementation notes (2026-08-05):**
- `update_B_hybrid` (`src/Pic.cpp`): removed the `euler` branch and the fallback
  `heun` block; the sub-cycle loop now dispatches only on `useRK4` (rk4) or
  `fieldIntegrator == "ssprk3"`, both of which `continue`.
- `#FIELDINTEGRATOR` parsing (`src/Pic.cpp:346-358`, `include/Pic.h:287-298`,
  `PARAM.XML`): accepted values are now only `rk4` (default) and `ssprk3`;
  anything else warns and defaults to `rk4`.
- `tests/hybrid_convection_wave/PARAM.in`: switched from `heun` to `rk4`.
- **Verified:** `pcai`, `beam` (full-PIC + hybrid), `hybrid_convection_wave`
  all PASS.

**Consideration (open):** the code comments note rk3 stays stable in the
high-amplitude phase where rk4 goes NaN. Since rk3 is also cheaper (3 vs 4
stages), it may be worth evaluating whether to make `rk3` the default for the
PCAI target, keeping `rk4` selectable. For now the decision is to keep **rk4 as
the default** and retain `rk3` as the selectable alternative.

---

### Item 2 — Only sync `nodePlasma` when output actually needs it (lazy output mirror)

- **Location:** `Pic::sum_moments`, `src/Pic.cpp:1302-1308`
- **Impact:** High (pure overhead every step)
- **Regression risk:** Low
- **Status:** implemented (2026-08-05)

**Problem.** In the hybrid path, `sum_moments(false)` is called **every step**
(`src/Pic.cpp:1728`). After the cell-centred deposit it ran the **nodePlasma
output bridge** — `average_center_to_node(centerPlasma[i] → nodePlasma[i])` for
every species **and** the summed entry `nSpecies+1`, each wrapped in
`FillBoundary` — **plus** `calc_mach_number()`, which reads `nodePlasma[nSpecies]`
(including the pressure tensor). On steps with no node-plasma / mach output this
is pure wasted work.

**Verified facts that drive the design:**
- `nodePlasma` is **plot-only** — it is **not** written to the checkpoint.
  `save_restart_data` writes only `nodeE`, `nodeB`, `centerB` and the particles
  (`src/PicIO.cpp:437-448`). On restart, `centerPlasma` is rebuilt by the
  `sum_moments()` deposit (`src/PicIO.cpp:519`).
- The only **every-step** consumer of `nodePlasma[nSpecies]` is
  `calc_mach_number` (`src/Pic.cpp:1330`), which produces the `mach` output
  diagnostic (not requested in the pcai/beam hybrid tests). `calc_cost_per_cell`
  reads `nodePlasma` only when balancing by Particle/Hybrid/Timing (default is
  `Cell`, which does not).
- `get_data_from_tracker` synthetic probes already read cell-centred fields in
  hybrid mode (Item 8-D), so they do not need `nodePlasma`.
- `get_field_var`/`get_var` read `nodePlasma`/`mMach` only via the structured
  plot write, which is preceded by the sync.

**Implementation (done):**
1. Added `nodePlasmaStale` member (`src/Pic.h:228`).
2. `sum_moments` hybrid branch: **skip** the `nodePlasma` bridge and
   `calc_mach_number`; set `nodePlasmaStale = true` (`src/Pic.cpp:1302-1308`).
   Full-PIC keeps the immediate `calc_mach_number`.
3. Added `Pic::sync_node_plasma_output()` (`src/Pic.cpp:1333`): if
   `nodePlasmaStale`, runs the `average_center_to_node` bridge for `nSpecies+1`
   + `calc_mach_number`, then clears the flag.
4. Call sites:
   - `Pic::write_plots` structured (ascii/IDL) branch, before `plot.writer.write`
     (`src/PicIO.cpp:658-662`).
   - `Pic::calc_cost_per_cell`, when balancing by Particle/Hybrid/Timing
     (`src/Pic.cpp:1376-1382`).
   `get_data_from_tracker` needs no sync (it reads cell-centred fields in hybrid).
5. **Mach gating:** `sync_node_plasma_output(bool needMach)` runs
   `calc_mach_number()` only when `needMach` is true. `write_plots` sets
   `needMach` from whether the plotString contains "mach"; `calc_cost_per_cell`
   passes false. Since `mach` is a pure output diagnostic consumed only by the
   `mach` plot variable (`src/PicIO.cpp:391`) and **not** by the GM coupling
   (which computes its own `mach2` from fluid moments via `get_fluid_state_for_points`,
   now cell-centred after Item 8-D), this fully eliminates `calc_mach_number`
   from the hybrid solver's default operation.
- **Verified:** `pcai`, `beam` (full-PIC + hybrid), `hybrid_convection_wave` all
  PASS.

---

### Item 3 — Drop the per-sub-cycle `project_centerB_to_nodeB` calls (node B is an output mirror)

- **Location:** `update_B_hybrid`, `src/Pic.cpp:2888-2891` (RK4), `2948-2951`
  (ssprk3); node sync after the loop at `src/Pic.cpp:2963-2993`
- **Impact:** High (× `nHallSubcycle`)
- **Regression risk:** Low–Med (verify `nodeB` is not read inside the loop)
- **Status:** implemented (2026-08-05)

**Problem.** Inside the sub-cycle loop the code called
`project_centerB_to_nodeB(0)` plus one call per fine level. Each such call does
`FillBoundary` + `apply_BC` on `centerB`, then `average_center_to_node(centerB →
nodeB)` + `apply_BC` on `nodeB`, plus a coarse→fine fill for fine levels — all
repeated `nHallSubcycle` times. The node-mirror (`nodeB`) and fine-level work is
pure output work; the only part actually needed between sub-steps is the
cell-centred `centerB` boundary condition so the next stage can read neighbours.

**Implementation:**
1. Added `Pic::apply_centerB_BC(int iLev)` — the cell-centred part of
   `project_centerB_to_nodeB` (`FillBoundary` + `apply_BC(centerB)` /
   `fill_fine_lev_bny_from_coarse(centerB)`).
2. In the RK4 and ssprk3 branches, replaced the per-sub-cycle
   `project_centerB_to_nodeB` + fine-level loop with a single
   `apply_centerB_BC(0)` (only level 0 is advanced in the sub-cycle).
3. Verified `nodeB` is not read inside the sub-cycle loop (the hybrid field
   advance reads only cell-centred fields; `dBdt` uses `nodeB^n` saved before the
   loop; the final `nodeB` is rebuilt once after the loop at
   `src/Pic.cpp:2976`).
4. The post-loop `average_center_to_node(centerB → nodeB)` and fine-level
   fills are unchanged.
- **Verified:** `pcai`, `beam` (full-PIC + hybrid), `hybrid_convection_wave` all
  PASS.

---

### Item 4 — Field solve only needs `ρ + 3·momentum`; shrink `centerPlasmaPrev` to 4 components

- **Location:** `assemble_ohm_E` reads `centerPlasmaSum`/`centerPlasmaPrev`
  `src/Pic.cpp:2456-2459` (only `iRho_, iUx_, iUy_, iUz_` used);
  `save_current_moments_to_prev` copies full `nMoments` `src/Pic.cpp:2656`
- **Impact:** Medium–High (saves the per-step `centerPlasmaPrev` copy width)
- **Regression risk:** Medium (touches layout) — must verify field solve reads
  the correct 4 components
- **Status:** implemented (2026-08-05)

**Problem.** The hybrid Ohm's law reads only the **charge density + 3 momentum
components** (`rho, mx, my, mz`) from `centerPlasmaSum` / `centerPlasmaPrev`
(`src/Pic.cpp:2456-2459`). The pressure tensor + number components
(`iPxx_ … iPyz_, iNum_`) are **never consumed** by the field solve — they feed
only the output path.

`nMoments = 11` (`include/Constants.h:39`) but `nHybridMomentsComps = 4`
(`include/Constants.h:48`). `save_current_moments_to_prev` copied the **full
`nMoments` (11)** summed array into `centerPlasmaPrev`, and
`seed_first_hybrid_step` did the same. The full-PIC path already uses a slim
4-component `nodePlasmaPrev`.

**Implementation (done):**
1. `centerPlasmaPrev` is now allocated with `nHybridMomentsComps = 4`
   components (`src/Pic.cpp:566-570`).
2. `save_current_moments_to_prev` copies only `nHybridMomentsComps` (4) from
   `centerPlasmaSum` (`src/Pic.cpp:2648-2653`).
3. `seed_first_hybrid_step` does the same (`src/Pic.cpp:2676-2682`).
4. `assemble_ohm_E` reads `momentsPrev` at `iRho_..iUz_` (indices 0-3), which
   match the 4-component layout — no index change needed.
5. This cuts the per-step `centerPlasmaPrev` copy width by 11/4 ≈ **2.75×**.
- **Verified:** `pcai`, `beam` (full-PIC + hybrid) all PASS.

> **On Items 5 and 6 (re-assessed after Item 2 + mach gating, 2026-08-05):**
> slim the per-particle deposit to 4 components, and defer
> `convert_to_fluid_moments`. With `calc_mach_number` now gated off (Item 2 /
> mach-gating commit) and the `nodePlasma` bridge deferred, the full pressure
> tensor (comps 4-10) is only needed at output time. **However, these still do
> not yield a clean, low-risk win:**
> 1. `convert_to_fluid_moments` (`src/Particles.cpp:1582`) **normalizes comps
>    0-9 by `qomSign*mass`** (via an alias `mult`) *and* converts the deposited
>    momentum+pressure into the thermal pressure tensor (comps 4-9). The
>    normalization of comps 0-3 is required by the field solve, and it is
>    entangled with the output-only pressure conversion in the same function, so
>    `convert_to_fluid_moments` cannot simply be deferred without splitting it.
> 2. The deposit (`sum_moments_cell_centered`) writes all 11 per-particle
>    moments (incl. pressure) into `centerPlasma[i]`. The field solve reads only
>    comps 0-3, but `centerPlasma[i]`/`centerPlasmaSum` must keep the full
>    `nMoments` for the `nodePlasma` output (which needs the thermal pressure).
>    Slimming the deposit to 4 comps therefore requires a **second, output-only
>    deposit** (re-running `sum_moments` in a full-width mode at output time),
>    since the pressure tensor cannot be reconstructed from `rho+3m` alone.
> 3. The deposit is per-particle (88 `+=`/particle); slimming to 4 comps would
>    save ~56 `+=`/particle on the field path but adds a full-width re-deposit at
>    output. The reward is real but the double-buffer / double-deposit refactor
>    risks the moment physics (normalization, thermal pressure), so it is left as
>    a documented follow-up rather than implemented now.

---

### Item 5 — Slim the per-particle moment deposit to the 4 field-needed components

- **Location:** `sum_moments_cell_centered`, `Particles.cpp:1000-1007`
- **Impact:** Medium (see Item 4)
- **Regression risk:** Medium (ties to Item 4) — field solve needs exactly the
  first 4 components
- **Status:** decided — implement

The deposit loop is:

```cpp
for (int iVar = 0; iVar < nMoments; iVar++)      // nMoments = 11
  for k,j,i in 8-cell stencil:
    momentsArr(ijk, iVar) += coef * pMoments[iVar];
```

This is **11 × 8 = 88 `+=` per particle**. Only `ρ + 3·momentum` (components
`iRho_, iUx_, iUy_, iUz_`) are needed by the field solve (`4 × 8 = 32`).
Depositing only the 4 field-needed components cuts this by ~64%. The remaining
7 components (`iPxx_ … iPyz_, iNum_`) are needed only for output.

**Implementation:**
1. In `sum_moments_cell_centered` (`src/Particles.cpp:957-1011`), split the
   deposit into two cases:
   - **Field-solve pass (always):** deposit only the first
     `nHybridMomentsComps = 4` components (`rho, mx, my, mz`) into
     `centerPlasma[iSpecies]` for the field solve. This is the hot path.
   - **Output pass (only when output is due):** deposit the full `nMoments`
     (pressure tensor + `iNum_`) into a full-width buffer that feeds the
     `nodePlasma` output bridge (Item 2). Gate this on the same lazy-output
     trigger so it is skipped on non-output steps.
2. This removes 7 × 8 = 56 `+=` per particle on the field-solve path (~64%).
3. Coordinate with Item 4 (4-component `centerPlasmaPrev`) and Item 6 (defer
   `convert_to_fluid_moments`): the field-solve `centerPlasma*` arrays hold only
   `rho + 3·momentum`; the pressure tensor is produced only for output.

---

### Item 6 — Only run `convert_to_fluid_moments` when output is required (same lazy trigger as Item 2)

- **Location:** `Pic::sum_moments`, `src/Pic.cpp:1289-1291`;
  `convert_to_fluid_moments`, `Particles.cpp:1582`
- **Impact:** Medium (per-species full-grid pass, every step)
- **Regression risk:** Medium — the pressure tensor is consumed only by output
- **Status:** decided — implement

**Problem.** `convert_to_fluid_moments` converts `ρ, ρU` into the pressure
tensor for each species (`src/Pic.cpp:1289-1291` → `Particles.cpp:1582`). The
pressure tensor is **not consumed by the hybrid field solve** — it only feeds
the `nodePlasma` output bridge (and probe / plot output). With Item 2 making the
`nodePlasma` bridge lazy, the per-step per-species `convert_to_fluid_moments`
call is also unnecessary on non-output steps.

**Implementation:**
1. Remove the unconditional `convert_to_fluid_moments` call from `sum_moments`
   (`src/Pic.cpp:1289-1291`).
2. Move it into the lazy output path: run it inside the same
   `sync_node_plasma_output()` helper from Item 2, just before the
   `average_center_to_node` bridge, so the pressure tensor is built only when
   output is actually requested.
3. Ensure `convert_to_fluid_moments` is applied to **each species**
   `centerPlasma[i]` and the summed `centerPlasmaSum[nSpecies]` before the
   `nodePlasma` bridge reads them.
4. Coordinate with Items 4 and 5: the field-solve `centerPlasma*` arrays stay at
   4 components; `convert_to_fluid_moments` needs the full `nMoments` (rho +
   3·momentum + pressure) to produce `nodePlasma`, so the output path may need a
   full-width deposit (Item 5 output pass) before this runs.

---

### Item 7 — Skip `nodeBavg` sync when hybrid (hybrid gather uses only `centerBavg`)

- **Location:** `update_B_hybrid`, `src/Pic.cpp:3012-3030`
- **Impact:** Low
- **Regression risk:** Low
- **Status:** implemented (2026-08-05)

**Problem.** When `useAvgFieldB` is enabled (defaults off), the end-of-step
averaging does `mult(alpha)` + `Saxpy` on both `centerBavg` **and** `nodeBavg`,
plus `FillBoundary` on both. The hybrid gather reads `centerBavg`; `nodeBavg`
is only used by the full-PIC node gather. Keeping `nodeBavg` in sync per-step
is redundant for the hybrid path.

**Implementation:**
1. Added `const bool syncNodeBavg = !useHybridPIC;` and guarded every
   `nodeBavg` `Copy`/`mult`/`Saxpy`/`FillBoundary` in the `useAvgFieldB` block
   (`src/Pic.cpp:3012-3030`) with `if (syncNodeBavg)`.
2. Verified no live hybrid-path code reads `nodeBavg`: the hybrid gather and
   Boris push read only `centerBavg` (`src/Pic.cpp:951`); the only `nodeBavg`
   readers are the full-PIC gather (`src/Pic.cpp:947`) and the dead
   `update_part_loc_to_half_stage` (`src/Pic.cpp:889`, never called).
3. `isBavgInit` is still driven by the `centerBavg` branch, so it stays
   consistent for `centerBavg`; only the `nodeBavg` mirror updates are skipped.
- **Verified:** `pcai`, `beam` (full-PIC + hybrid) all PASS. (Full-PIC keeps
  `nodeBavg` updated, so it is unaffected.)

---

### Item 8 — Evaluate: switch the hybrid plot writer to read cell-centered fields directly

- **Location:** output readers in `src/PicIO.cpp` (`get_var` `273-427`,
  `find_output_list` `82-219`, `get_field_var` `222-270`,
  `get_fluid_state_for_points`/probes `22-79`, `write_amrex_field` `780-1000`);
  checkpoint/restart `src/PicIO.cpp:430-521`; node-sync bridges in `src/Pic.cpp`
  (`update_B_hybrid`, `sum_moments`)
- **Impact:** **High — the largest structural + performance win.** Eliminates the
  per-step node-sync bridges (`average_center_to_node` for `centerB`, `centerE`,
  `centerPlasma*`), their `FillBoundary`, the output-side
  `average_node_to_cellcenter` conversions, and the associated per-step memory
  traffic — for the hybrid path.
- **Regression risk:** High (touches plot output + probe + restart) — must be
  phased and verified against pre-change plot/restart output
- **Status:** decided direction — plot writer reads cell-centered fields for
  hybrid; scope evaluated below

---

#### 8.1 Verified layout facts that make this feasible

The cell-centered equivalents of every node output field **already exist** on
the cell-centered grid `cGrids[iLev]` with ghost cells, allocated only under
`useHybridPIC` (`src/Pic.cpp:469, 546-575`):

| Output var (node field today) | Cell-centered equivalent | Grid | Availability |
|---|---|---|---|
| `Bx/By/Bz` (`nodeB`) | `centerB` | `cGrids`, 3 comps, `nGst` | **yes** |
| `Ex/Ey/Ez` (`nodeE`) | `centerEhybrid` | `cGrids`, 3 comps, `nGst` | **yes** |
| `rho/ux/…/pS` per species (`nodePlasma[iSpecies]`) | `centerPlasma[iSpecies]` | `cGrids`, `nMoments`, `nGst` | **yes** |
| summed moments (`nodePlasma[nSpecies]`) | `centerPlasmaSum[nSpecies]` | `cGrids`, `nMoments`, `nGst` | **yes** |
| `jHat` (`nodeJ`-adjacent) | `centerJ` | `cGrids`, 3 comps, `nGst` | **yes** |

So the data needed to write cell-centered plots is **already present** in the
hybrid layout — the plot writer simply reads the wrong (node) arrays today.

#### 8.2 Critical latent bug this change also fixes

In hybrid mode, `update()` does **not** call `update_E()` (the hybrid branch
calls `update_B_hybrid` instead, `src/Pic.cpp:1723-1732`). `update_B_hybrid`
syncs **only `nodeEth`** (`src/Pic.cpp:3138`, `average_center_to_node(centerEhybrid →
nodeEth)`), **not `nodeE`**. `nodeE` is written only in `update_E()`
(`src/Pic.cpp:1941-1960`), the full-PIC path.

But the plot writer reads **`nodeE`** for E output (`get_var` `src/PicIO.cpp:289-296`;
`write_amrex_field` `src/PicIO.cpp:885, 888`). **So hybrid E-field output
currently reads a stale `nodeE`.** Switching the plot writer to read
`centerEhybrid` directly fixes this latent output bug as a side effect.

#### 8.3 Scope of the plot-writer change (per function)

**A. `get_var` (`src/PicIO.cpp:273-427`)** — the variable→field mapping core.
- `Ex/Ey/Ez`: read `centerEhybrid[iLev]` instead of `nodeE[iLev]` (lines 289-296).
- `Bx/By/Bz`: read `centerB[iLev]` instead of `nodeB[iLev]` (lines 298-305).
- `rhoS/uxS/…/pXXS/pS/…` per species: read `centerPlasma[extract_int(var)][iLev]`
  instead of `nodePlasma[…]` (lines 335-379).
- `jHat`: read `centerJ[iLev]` instead of `jHat[iLev]` (lines 307-314).
- `dBdt`: `dBdt` is node-centered today; needs a cell-centered `dBdt` (or read
  `(B^{n+1}-B^n)/dt` computed on `centerB`) — small addition.
- `mach`: `mMach` is node-centered; either add a cell-centered `mMach` or compute
  on the fly from `centerPlasmaSum`. Minor.
- These reads must be gated behind `useHybridPIC` so the full-PIC path is
  untouched. In the hybrid branch, read cell-centered arrays; in the `!hybrid`
  branch, keep the current node reads.

**B. `find_output_list` / `get_field_var` (`src/PicIO.cpp:82-270`)** — these
iterate `MFIter mfi(nodeE[iLev])` (lines 108, 238) to get the box layout and use
`nodeStatus` (line 111) for owner/refined checks. For cell-centered output, the
hybrid branch must iterate `MFIter mfi(centerB[iLev])` with `cellStatus`, and
`write_amrex_field` already allocates `out[iLev]` on `cGrids` in the `!saveNode`
case (`src/PicIO.cpp:822`). The periodic right-edge trim (`src/PicIO.cpp:116-121`)
is node-specific and must be adjusted for the cell-centered layout (one fewer
node per dim).

**C. `write_amrex_field` (`src/PicIO.cpp:866-951`)** — in the `!saveNode`
(cell-centered, the default since `save_node()` is true only when the plotString
contains "node") branch:
- B: currently `average_node_to_cellcenter(nodeB)`? — **no**: B already copies
  `centerB` in the `!saveNode` branch (`src/PicIO.cpp:872`). Only E and plasma
  need changing.
- E: currently `average_node_to_cellcenter(nodeE)` (`src/PicIO.cpp:888`). Replace
  with `centerEhybrid` directly.
- plasma: currently `average_node_to_cellcenter(pl)` from `nodePlasma`
  (`src/PicIO.cpp:950`). Replace with `centerPlasma[iSpecies]` / the summed entry.
- The `saveNode = true` branch (explicit "node" plot requested) can be kept
  reading node fields for the hybrid case by materializing them lazily, or
  changed to write the node-average of the cell-centered fields.

**D. Probes / GM coupling `get_fluid_state_for_points` (`src/PicIO.cpp:22-79`)**:
reads `nodePlasma[iSpecies]`, `nodeB`, `nodeE` via `get_value_at_loc`
(`src/PicIO.cpp:63-73`). For hybrid, read `centerPlasma[iSpecies]`, `centerB`,
`centerEhybrid` instead. `get_value_at_loc` works on cell-centered MultiFabs too
(it interpolates at arbitrary locations), so no helper change is needed — only
the field selection.

**E. Restart (`save_restart_data` / `read_restart`, `src/PicIO.cpp:430-521`):**
- `save_restart_data` currently writes `nodeE`, `nodeB`, `centerB`
  (`src/PicIO.cpp:437-442`). For hybrid, write `centerEhybrid`, `centerB`
  (+ optionally `centerEprev`/`centerBprev`/`centerJ`/`centerPlasmaSum`) instead
  of `nodeE`. `nodeB` → `centerB`.
- `read_restart` currently reconstructs `centerEhybrid`/`centerEprev`/`centerBprev`
  from `nodeE`/`centerB` (`src/PicIO.cpp:500-511`). For hybrid, read the
  cell-centered fields directly — removing the `average_node_to_center(nodeE →
  centerEhybrid)` reconstruction. `nodePlasma` is not checkpointed, so no change
  there; `centerPlasma` is rebuilt by `sum_moments()` (`src/PicIO.cpp:519`).

#### 8.4 What the node-sync bridges can then be removed

Once the plot/probe/restart readers use cell-centered fields, the hybrid path no
longer needs to maintain the node mirrors every step:
- `average_center_to_node(centerB → nodeB)` at end of `update_B_hybrid`
  (`src/Pic.cpp:3060`) — **removable** for hybrid (nodeB needed only if an
  explicit "node" plot is requested).
- `average_center_to_node(centerEhybrid → nodeEth)` (`src/Pic.cpp:3138`) —
  **removable** for hybrid.
- `average_center_to_node(centerPlasma → nodePlasma)` in `sum_moments`
  (`src/Pic.cpp:1307-1313`) — **removable** for hybrid (Item 2 already makes it
  lazy; this change removes the need entirely for the default cell-centered plots).
- Items 2 and 6 become largely moot for the default output path (the pressure
  tensor / nodePlasma are no longer produced unless an explicit node plot is
  requested).

#### 8.5 Estimated effort

| Piece | Effort | Risk |
|---|---|---|
| A. `get_var` field redirection (gated by `useHybridPIC`) | Small | Low |
| C. `write_amrex_field` E/plasma → cell-centered (`!saveNode`) | Small | Low |
| D. probe `get_fluid_state_for_points` redirection | Small | Low |
| B. `find_output_list`/`get_field_var` layout + cellStatus + periodic trim | Medium | Medium |
| E. restart write/read cell-centered fields | Small–Med | Med (format change) |
| New `cell` `dBdt` / `mMach` (if output requests them) | Small | Low |
| Remove the node-sync bridges for hybrid | Small | Medium (verify no reader) |

Overall: **a focused, low-to-medium-effort change localized to `src/PicIO.cpp`
and the two node-sync sites in `src/Pic.cpp`.** It eliminates the per-step
`average_center_to_node` traffic and the output-side `average_node_to_cellcenter`
conversions, and it fixes the stale `nodeE` output bug. The main risk is the
`find_output_list`/`get_field_var` layout switch (node→cell box indexing) and
the restart format change.

**Recommendation:** implement this **as the primary approach** (over Items 2/6's
lazy sync), since it removes the node mirrors' per-step cost outright and is the
cleanest. Keep the node MultiFabs **allocated** (for full-PIC and for explicit
"node" plots) but stop syncing them in the hybrid field advance.

---

#### 8.6 Implementation status (2026-08-05)

**Implemented and verified** (gated by `useHybridPIC`, full-PIC bit-identical):

| Piece | Status |
|---|---|
| **8-C** `write_amrex_field` E/plasma → cell-centered for hybrid (`!saveNode`) | done |
| **8-D** probes `get_fluid_state_for_points` → cell-centered for hybrid | done |
| **8-E** restart write/read `centerEhybrid`/`centerB` directly for hybrid (no node reconstruction) | done |

**Implemented and reverted (critical finding):** the **structured (ascii/IDL)
plot path** (`find_output_list`, `get_field_var`, `get_var`) **cannot be switched
to cell-centered coordinates** without breaking plane-slice region selection.
`is_inside_plot_region` selects a `z=0` (or `y=…`) plane by comparing the output
coordinate against the region bounds (`plotMin_D`/`plotMax_D`, `src/PlotWriter.cpp:58`),
with a tolerance of `0.01·dx`. Node coordinates (`LoEdge`) hit the plane exactly,
but **cell-centre coordinates (`CellCenter`) sit at `0.5·dx` from the plane and
fall outside the region**, so a `z=0 fluid` plot emits only one (all-zero) point.
This regression was reproduced and reverted. The structured plot path is
inherently node-centric and must keep reading the node mirrors.

**Deferred — node-sync bridge removal (Item 8-sync, status 2026-08-05):**
- `centerEhybrid → nodeEth` bridge: **removed**. `nodeEth` has no live consumer
  in the hybrid path (the hybrid Boris push reads `centerEhybrid`; the structured
  plot reads `nodeE`, not `nodeEth`; `update_part_loc_to_half_stage` is dead
  code). Removes the per-step `average_center_to_node` + `FillBoundary` +
  `apply_BC` on `nodeEth`.
- `centerPlasma → nodePlasma` bridge: **already deferred** (Item 2).
- `centerB → nodeB` bridge + `dBdt`: **still per-step**. `nodeB` is read by the
  structured (ascii) B output (`By/Bz` in the pcai test) and by the `dBdt`
  computation, which is entangled: `dBdt = (B^{n+1} - B^n)/dt` needs a fresh
  `nodeB^{n+1}` and the `B^n` saved at the start of `update_B_hybrid` (currently
  copied from `nodeB`). Deferring both would require saving `B^n` into a new
  cell-centred buffer (`centerB_n`) at the start of the update and computing
  `nodeB`/`dBdt` lazily — a moderate-risk change, left as a follow-up.

**Test runner change:** `tests/validate_tests.py` now auto-builds `PostIDL.exe`
via `make PIDL` (new `ensure_postidl()` helper, called from `prepare_run_dir()`)
so every test has a valid `PostIDL.exe` for `PostProc.pl`. Verified: `pcai`,
`beam` (full-PIC) and `beam` (hybrid) all pass.

---

## 3. Summary Table

| # | Opportunity | Where | Impact | Risk | Status |
|---|---|---|---|---|---|
| 1 | Remove `euler` / `heun` integrators; keep `rk3`/`rk4`, `rk4` default | `src/Pic.cpp:2754-3143` | Medium–High | Low | **implemented** |
| 2 | Lazy `nodePlasma` output mirror (sync only when output needed) | `src/Pic.cpp:1302` | High | Low | **implemented** |
| 3 | Drop per-sub-cycle `project_centerB_to_nodeB` (output mirror) | `src/Pic.cpp:2888`, `2948` | High | Low–Med | **implemented** |
| 4 | Field solve only needs `ρ+3m`; shrink `centerPlasmaPrev` to 4 comps | `src/Pic.cpp:2656` | Medium–High | Med | **implemented** |
| 5 | Slim per-particle deposit to 4 field-needed components | `Particles.cpp:1000` | Medium | Med | decided |
| 6 | Lazy `convert_to_fluid_moments` / pressure only at output | `src/Pic.cpp:1289` | Medium | Med | decided |
| 7 | Skip `nodeBavg` sync when hybrid | `src/Pic.cpp:3012` | Low | Low | **implemented** |
| 8 | Cell-centered output/restart for hybrid (8-C/D/E done; 8-A/B reverted); 8-sync: `nodeEth` bridge removed, `nodeB`+`dBdt` deferred | `src/PicIO.cpp`, `update_B_hybrid` | High (largest) | High | partial — 8-C/D/E done; 8-sync nodeEth done |

---

## 4. Suggested Implementation Order

The **biggest and safest** wins are Items **1, 2, 3** — all pure removals of
redundant per-sub-cycle / per-step work in the hybrid path (Item 1 removes the
`euler`/`heun` integrator branches and their scratch arrays; Item 2 gates the
`nodePlasma` output bridge; Item 3 drops the per-sub-cycle `nodeB` projection).
They do **not** change the data layout or the numerics, so output should be
bit-identical (modulo the removal of purely-wasted work). They are easy to
verify with the existing `tests/beam` (full-PIC) + `tests/pcai` (hybrid)
regression checks.

> **Plan deviation note:** the implementation plan (Decision #2, Phase C.1)
> explicitly kept the `heun` integrator and its `centerBstar_heun` /
> `centerBstart_heun` / `dBpred_heun` / `dBcorr_heun` scratch arrays. Removing
> `heun` (Item 1) reverses that decision; the scratch arrays become dead and
> can be freed. This should be recorded when the change is implemented.

> **Updated priority:** Item **8** (plot writer reads cell-centered fields for
> hybrid) is now the **primary approach**. It removes the node mirrors'
> per-step sync cost outright (better than Items 2/6's lazy-sync, which keep the
> node mirrors but defer them), **and** it fixes the latent stale-`nodeE` output
> bug in hybrid mode. It is a focused, low-to-medium-effort change localized to
> `src/PicIO.cpp` plus the two node-sync sites in `src/Pic.cpp` (see §8.5).

Recommended order:
1. **Item 8** first — switch the hybrid plot/probe/restart path to cell-centered
   fields and remove the per-step node-sync bridges (`centerB→nodeB`,
   `centerEhybrid→nodeEth`, `centerPlasma→nodePlasma`). This makes Items 2 and 6
   mostly moot for the default (cell-centered) output path.
2. **Item 1** (remove `euler`/`heun`) and **Item 3** (drop per-sub-cycle
   `project_centerB_to_nodeB`) — independent pure removals; do alongside Item 8.
3. **Items 4 and 5** — interdependent (share the "only `ρ+3m` are needed"
   insight), touch the moment layout; slightly higher risk.
4. **Item 6** — if Item 8 lands, `convert_to_fluid_moments`/`nodePlasma` are only
   needed for an explicit "node" plot; otherwise defer per Item 6.
5. **Item 7** — small, independent, low-risk cleanup.

If Item 8 is **not** pursued, then Items 2 and 6 (lazy `nodePlasma` sync +
deferred `convert_to_fluid_moments`) are the fallback to reduce the per-step
cost while keeping the node-mirror output path unchanged.

---

## 5. Regression Safety

Every item must be verified against the two existing checks:

| Test | Config | Expected |
|---|---|---|
| `tests/beam` | full-PIC, `useHybridPIC = F` | bit-identical energy diagnostics |
| `tests/pcai` | hybrid, `useHybridPIC = T` | `d\|B\|` / `P_perp/P_par` match the Phase-A reference |

For Items 1–3, 7 (work/branch removal that never changes the active numerics —
`rk4` is the default integrator and the rk3/rk4 stages are unchanged by Item 1;
Items 2–3 remove pure output-mirror work; Item 7 skips a redundant node mirror),
the full-PIC `tests/beam` run is bit-identical by construction. Item 1
additionally requires a hybrid run confirming the `rk4` path is numerically
unchanged after the `euler`/`heun` branches and scratch arrays are deleted.

For Items 4–6 (layout change), the numeric results must be unchanged because
the removed pressure-tensor work is never read by the field solve — but this
must be verified explicitly, and a hybrid run must confirm the field solve reads
the correct 4-component `rho + 3·momentum` moments. Items 2 and 6 require a
**plot/restart correctness check**: run a hybrid case that writes a plot file
and a checkpoint, and verify the output matches the pre-change output (the node
mirrors and the cell-centered fields must agree wherever output reads them).

Item 8 changes the output/restart data source, so it requires:
- The full-PIC `tests/beam` run must remain **bit-identical** (all Item 8 changes
  are gated behind `useHybridPIC`; the full-PIC path reads the node fields as
  before).
- A hybrid run that writes a **cell-centered plot** and a **checkpoint**, and
  verifies the written values match the pre-change run **to the intended
  tolerance**. Note: the current hybrid E output reads stale `nodeE` (see §8.2),
  so the post-Item-8 E output will be the *correct* `centerEhybrid` value — this
  is a **fix**, not a regression, but it must be recorded as an expected
  difference in the validation.
- A restart round-trip check: write a checkpoint, restart, and confirm the
  reconstructed cell-centered fields (`centerEhybrid`, `centerB`, `centerPlasma`)
  match the pre-change restart behavior.
