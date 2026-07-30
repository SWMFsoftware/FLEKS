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

The solver-selection commands (`#SOLVEEM`, `#HYBRIDPIC`, `#RESISTIVITY`,
`#ELECTRONTEMPERATURE`, `#HALLSUBCYCLE`, `#PLASMA`) are documented as comments
directly in [`PARAM.in`](PARAM.in).

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
