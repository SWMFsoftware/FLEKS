# Hybrid PIC Whistler–Alfven Wave Standalone Test (Hall-only)

This test validates the **Hybrid PIC (kinetic ions, massless fluid electrons)**
solver in standalone FLEKS, with emphasis on the **Hall term** of the
generalized Ohm's law and its consistency with FLEKS's normalized CGS units
(see `docs/Algorithm.tex`).

It is the **Hall-only** member of the hybrid-PIC wave test family: the
resistive (`#RESISTIVITY`) and electron-pressure-gradient
(`#ELECTRONTEMPERATURE`) terms are disabled (`0.0`), so only the convection
(`-U_i×B`) and Hall (`(J×B)/(en_e)`) terms are active. The companion
[`hybrid_ohm`](../hybrid_ohm/README.md) test enables *all* four terms.

## Physics & Configuration

A 1D domain (64 cells, periodic in x) carries a uniform guide field `Bx0`
along x and a small circularly-polarized transverse perturbation
`By = b1 cos(k0 x)`, `Bz = b1 sin(k0 x)` (seeded by the `HybridWave`
initializer). A matching **Alfvén-wave ion velocity perturbation**
`δU_⊥ = -v_A · δB_⊥ / B0` is applied to the particles (in code units
`v_A = B0 = 1`, so `δU_y = -b1 cos(k0 x)`, `δU_z = -b1 sin(k0 x)`) via
`Pic::perturb_hybrid_wave_velocities`, which runs after `fill_particles` and
before `sum_moments`. This self-consistent eigenmode IC avoids the large
transient that a B-only seed would produce.

### Active solvers (require the hybrid port on `feature/hybrid-pic`)
- `#SOLVEEM = F` — standard Maxwell/GMRES solver disabled.
- `#HYBRIDPIC = T` — Ohm's-law + Faraday hybrid field advance.
- `#RESISTIVITY = 0.0` [m²/s] — Hall-only (isolates dispersion).
- `#ELECTRONTEMPERATURE = 0.0` [eV] — no electron-pressure-gradient term.
- `#HALLSUBCYCLE = 32` — sub-cycles the B update so `dt_sub = dt/32`
  satisfies the Hall (whistler) CFL; see *Implementation notes* below.

`#PLASMA` declares a single kinetic ion species (`q = m = 1`); the electron is
an implicit neutralizing fluid, so no electron particle species is present.

## Unit conventions (verified against `docs/Algorithm.tex` and master code)

- `E` and `B` share the same code unit (`E* = B*`, confirmed by
  `fill_lightwaves` with `|E| = |B|` and by `update_B`'s `∂B/∂t = -∇×E`).
- The plasma current density obeys the CGS Ampère law
  `∇×B = 4π J` (see `Pic::update_E_rhs`: `∂E/∂t = ∇×B − 4π·jHat`, with
  `jHat = Σ q v` from `PicParticles::calc_jhat`). Therefore the current used in
  Ohm's law **must** be `J = ∇×B / (4π)`, **not** raw `∇×B`.
- The Hall/pressure denominator `ρ` is `nodePlasma[...].iRho_`, which
  `PicParticles::sum_moments` stores as the **charge density** `ρ_q = e n_e`
  (`iRho_ = qp`). This matches the CGS Ohm's-law denominators
  `(J×B)/(e n_e)` and `∇p_e/(e n_e)` exactly — so no extra `q/m` factor is
  needed there.

**Consequence for the draft on `origin/hybrid`:** its `update_E_hybrid`
computes `J = ∇×B` (raw curl) and divides the Hall and resistive terms by it.
That makes both terms **4π too large**. The fix (roadmap Step 4) is to define
`J = ∇×B / fourPI` inside `update_E_hybrid`. The convection (`-U_i×B`, identical
to `Pic::update_U0_E0`) and pressure-gradient (`∇Pe/ρ_q`) terms are already
unit-consistent and need no change.

## Expected Results & Success Criteria

1. **Stability** — FLEKS exits cleanly; no NaN/Inf; magnetic energy `Eb` and
   ion energy `Epart` remain finite and bounded (the `4π` bug or insufficient
   `nHallSubcycle` makes the whistler blow up).
2. **Particle-number (not energy) conservation** — with periodic BCs and no
   source/loss the ion *number* is conserved exactly; the ion *kinetic energy*
   is **not** (it continuously exchanges with the electromagnetic field, and
   the first-order `E^n` push documented in the roadmap adds a small numerical
   drift).  The automated check therefore only guards against *gross*
   non-conservation (ion-energy ratio outside `[0.2, 10]`); the decisive
   stability guards are the magnetic-energy blow-up check (`Eb` < 5× initial)
   and the bounded transverse-wave check (Section 4).
3. **Wave propagation / Alfvén speed** — the transverse pulse propagates at the
   Alfvén speed `v_A`; for parallel propagation the low-frequency branch matches
   the MHD Alfvén speed.
4. **Hall (whistler) dispersion** — for parallel propagation the measured
   frequency obeys the normalized whistler relation
   ```
   ω / Ω_i = (k d_i)^2 / (1 + (k d_i)^2)
   ```
   where `Ω_i = q B / m` is the ion gyrofrequency and `d_i = v_A / Ω_i` the ion
   inertial length. This dimensionless relation is **unit-invariant** and is the
   decisive test of the Hall term: with the `4π` bug present, the effective Hall
   coefficient is `4π` too large, so the measured `ω` is `4π` too high.
   *Verification:* run with time-resolved `Bx,By,Bz` output (see `#SAVEPLOT`),
   Fourier-transform `By/Bz` in space and time, and fit `ω(k)` to the relation
   above. The automated `validate_tests.py` check performs a two-stage
   verification: (a) **early time** — the spatial DFT of `By` must peak at mode
   `n=1` (one wavelength in the box, confirming the seed), and (b) **late time**
   — the transverse amplitude must stay below 10×B₀ (catching a missing-4π
   blow-up or CFL violation). At moderate PPC the Hall term amplifies grid-scale
   particle noise over long times (a well-known hybrid-PIC limitation); the
   early-time `n=1` check validates the solver physics, while the late-time
   bound catches genuine instabilities. The `ω(k)` fit remains the reference
   discriminator for manual validation.

## Implementation notes

- `#TESTCASE HybridWave` must be added to the `TestCase` enum and an initializer
  (circularly-polarized `By,Bz` perturbation on `Bx0`) implemented in
  `Pic::fill_new_cells`, mirroring `Pic::fill_lightwaves`. The `#TESTCASE`
  string in `PARAM.in` must match the code comparison **exactly**
  (`"HybridWave"`, case-sensitive) — a mismatch silently falls through to
  `RegularSimulation` and the wave is never seeded.
- **Particle count (PPC).** The test uses 2000 particles per cell in X. At
  lower PPC (e.g. 100) the Hall term amplifies grid-scale density noise into a
  noisy instability that can dominate the seeded `n=1` wave at late times.
  2000 PPC keeps the growth factor to ~16× (vs ~190× at 100 PPC) and ion-energy
  conservation to ~0.4%.
- The `HYBRIDPIC` / `RESISTIVITY` / `ELECTRONTEMPERATURE` / `HALLSUBCYCLE`
  commands must be present in `PARAM.XML` (added during the hybrid port).
- **Stability / Hall CFL (the binding constraint).** The Hall term in Faraday,
  `∂B/∂t = -∇×(J×B)/ρ`, is a *dispersive* term whose explicit forward-Euler
  advance is limited by the **discrete** operator, not the continuous
  `ω ≤ Ω_i` bound.  The grid-scale (short-wavelength) modes behave like a
  magnetic diffusion with coefficient `D_H = v_A·d_i`, which gives the 3-D
  stability limit
  ```
  dt_sub ≲ dx² / (6 · v_A · d_i)        with  dt_sub = dt / nHallSubcycle
  ```
  With `dx = 6.4/64 = 0.1` and `v_A = d_i = 1` this requires
  `dt_sub ≲ 1.7e-3`; for `dt = 0.02` that means `nHallSubcycle ≥ 12`.  The test
  uses `nHallSubcycle = 32` (`dt_sub = 6.25e-4`, comfortable margin).  Note this
  is **independent of the physics frequency** `ω ≤ Ω_i`: even though
  `dt·Ω_i ≪ 1` is satisfied, too few sub-cycles still makes the whistler blow
  up.  (The per-sub-step `update_E_hybrid()` only re-curls `B` to refresh `J`;
  the ion moments `U_i, ρ_q, P_e` are held fixed during the sub-cycle, so extra
  sub-cycles add negligible cost.)
- **Normalization (critical).** The standalone test MUST use *ion-scale*
  normalization so that in CODE (normalized) units the ion scales are `O(1)`:
  set `#NORMALIZATION lNormSI ≈ d_i` and `uNormSI ≈ v_A`.  Then
  `B_code ≈ Ω_i ≈ d_i ≈ v_A ≈ 1`, `dt·Ω_i ≈ 0.02 ≪ 1`, and the Hall CFL above is
  achievable.  If instead a naive `lNormSI = uNormSI = 1` (SI) normalization is
  used, `Si2NoB ≈ 9.6e7` makes `B_code ≈ 1e8` and `Ω_i_code ≈ 1e8`, so
  `dt·Ω_i ≈ 1e6` and the explicit scheme diverges in a **single** step — no
  amount of sub-cycling can cure it.  The GM-PC `beam` test avoids this by using
  realistic SI inputs (`B = 5 nT`, `n = 5 /cc`, `lNormSI = 1000 m`,
  `uNormSI = 1e6 m/s`) so its code units are also ion-scale `O(1)`.
- **Plot format.** The `#SAVEPLOT` plot string must use the **ascii** format
  (e.g. `z=0 fluid ascii pic`), **not** `real4`/`real8`: `validate_tests.py`
  reads the text `.out` header to extract `By`/`Bz` for the bounded-wave check,
  and the manual `ω(k)` fit also needs text `.out`.  Binary plot formats write
  Fortran-record `.h` files that the harness cannot parse.
