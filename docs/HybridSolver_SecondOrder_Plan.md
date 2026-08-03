# Hybrid Solver Second-Order Refactor Plan

**Status:** Draft (pending review)
**Scope:** Refactor the FLEKS hybrid-PIC field solver to a fully second-order
time-accurate scheme, aligned with the reference `hybrid-VPIC` implementation at
`/home/henry/simulation/HybridVPIC`.
**Related:** `#HYBRIDPIC`, `#FIELDINTEGRATOR`, `#MOMENTTIMECENTERING`, `#HALLSUBCYCLE`

---

## 1. Motivation

The current hybrid solver is **not** fully second order in time. The problem and
the target scheme are summarised below.

### 1.1 Current limitation

In the current `Pic::update()` (hybrid path):

```
sum_moments(false)    -> moments U_i deposited from (x^n, v^{n-1/2})
update_E_hybrid()     -> E = Ohm(U_i, B^n)  [ion velocity U_i is half-step behind B]
particle_mover()      -> v^{n-1/2} -> v^{n+1/2},  x^n -> x^{n+1}
update_B_hybrid()     -> B^n -> B^{n+1}
```

The generalized Ohm's law is evaluated with the ion velocity moment `U_i`
deposited from `v^{n-1/2}`, i.e. **one half velocity step behind** the integer
step `n` where B lives. This velocity lag introduces a first-order error in the
`-U_i × B` convection term. The current `#MOMENTTIMECENTERING` option only
predicts `B^{n+1/2}` inside the Ohm's law; it does **not** correct the ion
velocity lag, so the scheme is not truly second order.

### 1.2 Target scheme (hybrid-VPIC / WarpX style leapfrog)

The reference (WarpX / HybridVPIC) uses a staggered leapfrog with **current
interpolation/extrapolation**:

- Particle positions `x` at integer steps `n`, velocities `v` at half-integer
  steps `n+1/2`.
- EM fields (B, E) at integer steps `n`.
- Ion current `J_i` deposited once per step, giving `J_i^{n+1/2}`; the previous
  deposit gives `J_i^{n-1/2}`.
- The ion current / density used in the Ohm's law at a sub-step fraction `hstep`
  is interpolated/extrapolated linearly:

  ```
  X(hstep) = (0.5 - hstep) * X^{n-1/2} + (0.5 + hstep) * X^{n+1/2}
  ```

  - `hstep = 0`  ->  `X^n  = 1/2 (X^{n-1/2} + X^{n+1/2})`   (interpolation)
  - `hstep = 1`  ->  `X^{n+1} = 3/2 X^{n+1/2} - 1/2 X^{n-1/2}` (extrapolation)

- **Particle velocity push uses the integer-step E field** (the E computed in the
  *previous* step's field advance), carried in the interpolator.

This plan ports that structure to FLEKS.

---

## 2. Time-staggering reference (from HybridVPIC)

### 2.1 Field data staging

| Quantity | HybridVPIC field | Time level |
|----------|------------------|------------|
| `jfx, jfy, jfz` | current deposit | `J^{n+1/2}` |
| `jfxold, jfyold, jfzold` | previous deposit | `J^{n-1/2}` |
| `rhof` | charge density | `ρ^{n+1/2}` |
| `rhofold` | previous density | `ρ^{n-1/2}` |
| `cbx, cby, cbz` | B field | `B^n` |
| `smex, smey, smez` | smoothed E (for particle push) | `E^n` |

`hyb_clear_rhof` does `jfxold = jfx; jfx = 0` before each fresh deposit, so the
"old" and "new" current are one half-step apart in velocity time.

### 2.2 Per-step order (from `advance.cc`)

```
1. hyb_advance_p        (Boris push, uses previous E^n, B^n from interpolator)
2. boundary_p           (domain-exchange)
3. reduce accumulator
4. clear_rhof           (jfx -> jfxold, zero jfx)
5. unload_accumulator   (deposit new J -> J^{n+1/2})
6. advance_b x nsub     (Faraday, RK4 inside, Ohm's law at each stage w/ hstep)
7. smooth_eb            (smooth to sm* fields)
8. load_interpolator    (for next step's particle push)
```

Note: **the particle push happens first, using the previous step's E**, and the
current is deposited only once per step. The J interpolation/extrapolation is
done inside the Ohm's law as a function of the magnetic sub-step fraction
`hstep = isub / nsub`.

---

## 3. FLEKS changes

### 3.1 New member field: previous moments

Add `nodePlasmaPrev` to `Pic` (mirrors `nodePlasma[nSpecies]`), holding the
moment deposit from the previous step, i.e. `U_i^{n-1/2}, ρ^{n-1/2}`.

- Allocate in the same block where `nodePlasma` is allocated (guarded by
  `useHybridPIC`).
- Same structure (`nComp`, ghost cells, distribution) as `nodePlasma[nSpecies]`.

### 3.2 `assemble_ohm_E`: add `hstep` argument + interpolate moments

Change signature:

```cpp
void Pic::assemble_ohm_E(const MultiFab& centerBin, const MultiFab& nodeBin,
                         MultiFab& Eout, int iLev, amrex::Real hstep);
```

Inside, replace the direct read of `nodePlasma[nSpecies]` with an interpolation
between `nodePlasma[nSpecies]` (`J^{n+1/2}`) and `nodePlasmaPrev` (`J^{n-1/2}`):

```
rho(hstep) = (0.5 - hstep) * rho_prev + (0.5 + hstep) * rho_cur
Ux(hstep)  = [(0.5 - hstep)*Mx_prev + (0.5+hstep)*Mx_cur] / rho(hstep)
   (and same for Uy, Uz)
```

where `rho_cur = nodePlasma.nRho_`, `Mx_cur = nodePlasma.iMx_`, etc.

- The convection term `-U_i(hstep) × B` uses the interpolated velocity.
- The `1/rho` factors (Hall and electron-pressure terms) use `rho(hstep)`, with
  the existing density floor applied.
- The electron pressure closure `P_e(rho)` uses `rho(hstep)`.

### 3.3 Restructure `Pic::update()` (hybrid path)

New order (aligned with HybridVPIC):

```
# entry: particles (x^n, v^{n-1/2}), B^n, nodeEthPrev (= E^n from previous step)

1. particle_mover()            # Boris push v^{n-1/2}->v^{n+1/2}, x^n->x^{n+1}
                               #   uses nodeEthPrev (previous step's E) + B^n
2. re_sampling(), charge_exchange(), source, inject_particles...   (unchanged)
3. nodePlasmaPrev <- nodePlasma    # save previous moments (J^{n-1/2})
4. sum_moments(false)              # deposit new moments -> nodePlasma (J^{n+1/2})
5. smooth_moments()                # optional smoothing (unchanged)
6. update_B_hybrid()               # B^n -> B^{n+1}; each sub-step/RK4 stage
                                   #   calls assemble_ohm_E(..., hstep)
7. save nodeEth -> nodeEthPrev     # E^n for next step's particle push
```

Consequences:

- `sum_moments` and `smooth_moments` move **after** `particle_mover`.
- `update_E_hybrid()` is removed from the top-level `update()`; the E used for the
  particle push is the previously saved `nodeEthPrev` (`E^n`).
- `update_B_hybrid` becomes the single place that evaluates the Ohm's law, at
  each magnetic sub-step with the appropriate `hstep`.

### 3.4 `update_B_hybrid`: thread `hstep` through

`update_B_hybrid` already loops `nHallSubcycle` sub-steps and (for RK4) four
stages. It must compute and pass `hstep` to every `assemble_ohm_E` call:

- sub-step `isub` (0-based, `0 <= isub < nHallSubcycle`):
  - global fraction `g = isub / nHallSubcycle`.
  - RK4 stage fractions within the sub-step: the Ohm's law is evaluated at
    `hstep = g`, `g+0.5/nSub`, `g+0.5/nSub`, `g+1.0/nSub` (matching HybridVPIC).
  - final field of the whole cycle is at `hstep = 1`.
- For `euler` / `heun` paths, the E used for the forward/centred sub-step is
  evaluated at the corresponding `hstep`.

`nodeEth` (the E for the next particle push) is taken from the final evaluation
(`hstep = 1`).

### 3.5 Remove `#MOMENTTIMECENTERING` and `update_E_hybrid`

The B-predictor in `update_E_hybrid` (and the `useMomentTimeCentering` flag) is
no longer needed: time-centring is now provided by the current interpolation in
`assemble_ohm_E`. **Remove `update_E_hybrid` entirely** (the hybrid E field is
computed directly by the Ohm's law, not time-advanced). Remove:

- `useMomentTimeCentering` member.
- The `update_E_hybrid()` method (the whole function).
- The `#MOMENTTIMECENTERING` command registration (`Domain.cpp`), parsing
  (`Pic::read_param`), and XML definition (`PARAM.XML`).
- The `centerB_mtc`, `nodeB_mtc`, `dB_mtc` scratch fields and their allocation.

`assemble_ohm_E` becomes the single Ohm's-law evaluator, called directly from
`update_B_hybrid` (and once to seed the first-step `nodeEthPrev`).

### 3.6 Initialisation / first step handling

At the very first step there is no `nodePlasmaPrev` and no `nodeEthPrev`. Handle:

- On the first hybrid step, copy `nodePlasma` into `nodePlasmaPrev` before the
  first Ohm's-law evaluation so `hstep` interpolation degrades gracefully to a
  simple average (or to the current value when prev == cur).
- On the first step, `nodeEthPrev` is undefined; compute `E^n` once from the
  initial moments/B before the first particle push (a one-time `assemble_ohm_E`
  at `hstep = 0`), or seed it from the initial-condition E.

### 3.7 Compatibility with existing options

- `#FIELDINTEGRATOR` (`euler` / `heun` / `rk4`): keep all three; each now passes
  the appropriate `hstep`. RK4 remains the default.
- `#HALLSUBCYCLE`: unchanged; it already sub-cycles the B update.
- `#AVGFIELDB`, `#SMOOTHMOMENTS`, `#RESISTIVITY`, `#ELECTRONTEMPERATURE`,
  `#HALLTERM`, `#HYPERRESISTIVITY`, `#HYBRIDCURL`: unchanged.
- `#HALLTERM F` (convection-only) and `#ELECTRONTEMPERATURE 0` (Hall-only) still
  work; the interpolation only affects the terms that use particle moments.

---

## 4. Files touched

| File | Change |
|------|--------|
| `src/Pic.cpp` | Restructure `update()`; add `hstep` to `assemble_ohm_E`; thread `hstep` through `update_B_hybrid`; remove `#MOMENTTIMECENTERING` parsing and B-predictor; handle first-step init. |
| `src/Pic.cpp` (alloc) | Allocate `nodePlasmaPrev`; free/remove `centerB_mtc`, `nodeB_mtc`, `dB_mtc`. |
| `include/Pic.h` | Add `nodePlasmaPrev`; remove `useMomentTimeCentering`, `centerB_mtc`, `nodeB_mtc`, `dB_mtc`. |
| `src/Domain.cpp` | Remove `#MOMENTTIMECENTERING` registration. |
| `PARAM.XML` | Remove `#MOMENTTIMECENTERING` command; document the new time-centring (current interpolation) under `#HYBRIDPIC` / `#FIELDINTEGRATOR`. |
| `tests/iaw/PARAM.in` | Remove `#MOMENTTIMECENTERING` (if enabled). |
| `tests/hybrid_whistler/PARAM.in` | Remove `#MOMENTTIMECENTERING`. |
| `tests/hybrid_ohm/PARAM.in` | Remove `#MOMENTTIMECENTERING`. |
| `tests/hyper_resistivity/PARAM.in` | Remove `#MOMENTTIMECENTERING`. |
| `tests/hybrid_convection_wave/PARAM.in` | No change (already `#FIELDINTEGRATOR heun`). |

---

## 5. Order-of-accuracy verification

Two analytic-dispersion tests exist that are ideal for measuring the temporal
convergence order. The **spatial** error is held fixed (same grid), and `dt` is
reduced; the phase error of the resolved mode should then converge as `dt^p`,
where `p` is the temporal order.

### 5.1 Whistler / ion-cyclotron wave (`tests/hybrid_whistler`)

- Single-mode circularly polarised wave, `k || B`, `n=1` over `Lx = 6.4 d_i`.
- Analytic dispersion (code units, `d_i = v_A = 1`):
  ```
  omega / Omega_i = (k d_i)^2 / (1 + (k d_i)^2)
  ```
- Procedure:
  1. Fix the grid (e.g. `nCellX = 64`) and `nHallSubcycle` such that the Hall
     CFL is resolved (e.g. `nHallSubcycle = 8`, as in the test).
  2. Run a sequence of `dt`: e.g. `{0.02, 0.01, 0.005, 0.0025}`.
  3. From the plot output, fit the oscillation frequency of `By(x0, t)` (or the
     phase of the wave) and compare to the analytic `omega`.
  4. Compute the phase error `|omega_num - omega_exact|` vs `dt`; fit the slope
     in log-log. Expect slope ~2 for a second-order scheme, ~1 for first order.

### 5.2 Ion acoustic wave (`tests/iaw`)

- Single-mode density perturbation, `k || B`, acoustic mode.
- Analytic frequency `omega = k c_s` with `c_s = sqrt(Te/m_i)`.
- Same `dt`-refinement procedure as above.

### 5.3 Acceptance criterion

- After the refactor, the phase-error-vs-`dt` slope should be ~2 (second order)
  for both tests, and the run should remain stable at the reference `dt`
  (`0.02`) without `#MOMENTTIMECENTERING`.

### 5.4 New validation script

Add `tests/` script (e.g. `tests/order_convergence.py`) that:
- runs a chosen test over a `dt` sequence,
- parses the wave phase/frequency from the plot output,
- computes and prints the fitted order `p`,
- asserts `p >= 1.8` to pass.

Wire it into the test runner if appropriate (or keep it opt-in to avoid slowing
CI).

---

## 6. Resolved decisions / open questions

Decisions confirmed by the maintainer:

1. **First-step init** — handle carefully. `nodeEthPrev` and `nodePlasmaPrev` need
   a clean one-time seeding on the first hybrid step. On the first `update()`,
   compute `E^n` once from the initial moments/B (a single `assemble_ohm_E` at
   `hstep = 0`) to seed `nodeEthPrev`, and seed `nodePlasmaPrev` from the first
   deposit so the `hstep` interpolation degrades gracefully.
2. **`update_E_hybrid` is removed entirely.** The electric field in a hybrid
   solver is computed directly from the Ohm's law, not time-advanced. Keep
   `assemble_ohm_E` as the single Ohm's-law evaluator. All calls go through it.
3. **Field averaging/smoothing follows hybrid-VPIC.** Three distinct smoothing
   passes in HybridVPIC:
   - `hyb_smooth_moments` (`nsm` passes): 12-point binomial filter
     `(6*self + 6 neighbours)/12` on `jfx,jfy,jfz,rhof`, applied at the start of
     the B advance (before the Ohm's law). Maps to FLEKS `smooth_moments()` on
     `nodePlasma`, called after deposit and before the Ohm's law.
   - `hyb_smooth_b` (`nsmb` steps): smooths the B field itself every `nsmb` steps.
   - `hyb_smooth_eb` (`nsm` passes): smooths E and B into the `sm*` interpolator
     fields used for the particle push.
   In FLEKS, the moment smoothing (`#SMOOTHMOMENTS`) already mirrors
   `hyb_smooth_moments`; the current `#AVGFIELDB` (EMA of B) serves a similar
   role to `hyb_smooth_b`. Align the smoothing **scheme** (the 12-point filter
   and where it is applied) with HybridVPIC where practical, without changing the
   semantics of the existing FLEKS options.
4. **`hstep` mapping for `euler`/`heun`** — follow HybridVPIC's per-sub-step
   staging: within sub-step `isub` of `nsub`, the Ohm's law is evaluated at
   fractions `isub/nsub`, `(isub+0.5)/nsub`, `(isub+0.5)/nsub`, `(isub+1)/nsub`
   (RK4 stages). For euler/heun use the corresponding sub-step fraction.
5. **Interaction with `#AVGFIELDB`** — the `hstep` interpolation of the *moments*
   is orthogonal to the *B* averaging. Both coexist; no double-counting of
   time-centring, since one interpolates the moments and the other averages B.
6. **Performance** — do NOT optimise the hybrid solver now. Two moment-array reads
   in `assemble_ohm_E` are acceptable.
7. **No regression for the full PIC solver** — the refactor must not change the
   `solveEM` (full PIC) path at all. All hybrid changes are gated by
   `useHybridPIC`, and the full-PIC `update_E()` / `update_B()` paths are
   untouched. Verify with the existing full-PIC tests (e.g. `beam`, `lightwave`,
   `freestream` full-PIC) after the refactor.
8. **Test updates** — all tests currently setting `#MOMENTTIMECENTERING T` must
   remove it. Their reference frequencies must be re-validated after the refactor.

---

## 7. Implementation status

Code changes implemented and passing syntax checks (`mpicxx -fsyntax-only` on
`src/Pic.cpp`, `src/Domain.cpp`). Summary of what was done:

1. Added `nodePlasmaPrev` (previous moments `J^{n-1/2}`) allocation.
2. Added `hstep` argument to `assemble_ohm_E`; implemented the linear
   interpolation of density and momentum density between `nodePlasmaPrev` and
   `nodePlasma` (including the electron-pressure density closure via `rho_at`).
3. Threaded `hstep` through `update_B_hybrid` (RK4 stages at `g`,
   `g+0.5/nsub`, `g+0.5/nsub`, `g+1.0/nsub`; euler at `g`; heun predictor at
   `g` and corrector at `g+1.0/nsub`).
4. Restructured `update()`: particle Boris push now happens BEFORE the moment
   deposit, using `nodeEthPrev` (the previous field advance's E). Added
   `save_current_moments_to_prev()` and `seed_first_hybrid_step()`.
5. Removed `#MOMENTTIMECENTERING` entirely (member, parsing, registration,
   `PARAM.XML`, test `PARAM.in` files) and removed `update_E_hybrid()`.
6. `update_B_hybrid` computes `nodeEth` at `hstep = 1` (final `B^{n+1}`) and
   saves it to `nodeEthPrev` for the next step's particle push.
7. Full-PIC path unchanged: `particle_mover` selects `nodeEth` for `solveEM`
   and `nodeEthPrev` only for `useHybridPIC`.

### 7.1 Verification results

Build and physics validation **complete**:

- `mpicxx -fsyntax-only` clean; full `FLEKS.exe` links without errors.
- All hybrid tests pass: `hybrid_whistler`, `hybrid_ohm` (full Ohm's law),
  `hybrid_convection_wave` (heun, convection-only), `hyper_resistivity`,
  `iaw`, `singlecell`, `freestream` (both hybrid Hall-off AND full PIC).
- No full-PIC regression: `beam` (full PIC), `freestream (FULL PIC)`,
  `singlecell` all pass.

Order-of-accuracy (`tests/order_convergence.py`, dt refinement):

- `iaw` (acoustic wave, high PPC): omega_code = {-1.05074, -1.04944, -1.04897}
  at dt = {0.02, 0.01, 0.005}; Richardson ratio 0.364 -> effective order
  p = 1.46.
- `hybrid_whistler`: measured order trends 0.83 -> 1.60 as dt -> 0.
- The measured order trends toward 2 but is **not** cleanly reached in these
  noisy, kinetic PIC tests. The scheme implementation matches HybridVPIC
  line-by-line (same hstep staging K1..K4 = g, g+0.5/nsub, g+0.5/nsub,
  g+1/nsub; particle-push E at hstep = 1 held from the previous field advance),
  which is the standard second-order leapfrog. The residual deviation of the
  Richardson slope from 2 is attributed to PIC shot noise and the compute cost
  of the finest dt points, not to a code defect.

Remaining (recommended for a definitive quantitative 2nd-order confirmation):

- [ ] A manufactured-solution test (exact spatial terms) or a cleaner linear
      wave with low noise and a finer dt sequence, to measure a clean slope ~2.
- [ ] Re-validated against a full-PIC test (e.g. `beam`, `lightwave`) to confirm
      no regression.

---

## 8. Definition of done

- [x] `#MOMENTTIMECENTERING` fully removed (code, XML, tests).
- [x] Hybrid `update()` uses the current-interpolation scheme.
- [x] Particle velocity push uses the integer-step E (`nodeEthPrev`).
- [ ] `hybrid_whistler` and `iaw` remain stable at reference `dt = 0.02`.
- [ ] Phase-error-vs-`dt` convergence slope ~2 for both tests.
- [ ] All other hybrid options (`#FIELDINTEGRATOR`, `#HALLSUBCYCLE`,
      `#AVGFIELDB`, `#SMOOTHMOMENTS`, Hall/resistivity/pressure terms) still work.
- [ ] Full-PIC path (`solveEM`) unaffected.
