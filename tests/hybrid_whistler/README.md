# Hybrid PIC Whistler–Alfven Wave Standalone Test (Hall-only)

This test validates the **Hall term** of the generalized Ohm's law in the
Hybrid PIC solver (kinetic ions, massless fluid electrons) against the
whistler dispersion relation. It is the Hall-only member of the hybrid-PIC
wave test family: resistivity and electron-pressure-gradient are disabled, so
only the convection (`-U_i×B`) and Hall (`(J×B)/(en_e)`) terms are active. The
companion [`hybrid_ohm`](../hybrid_ohm/README.md) test enables all four terms.

## Physics & Configuration

A 1D domain (64 cells, periodic in x) carries a uniform guide field `Bx0` along
x and a small **circularly-polarized transverse** perturbation
`δB = (0, b1 cos(k0 x), b1 sin(k0 x))`, seeded by the `HybridWave` initializer
(`src/ic/WaveIC.cpp`). A matching Alfvén ion velocity kick
`δU = -v_A δB_⊥ / B0` is applied to the particles, forming a self-consistent
eigenmode. Parameters are chosen so that in code units `d_i ~ v_A ~ Ω_i ~ 1`
with `dt = 0.02`; physical inputs (`n = 5 /cc`, `B = 5.22 nT`) are
solar-wind-like.

Solver-selection commands are documented as comments in
[`PARAM.in`](PARAM.in).

## Unit conventions

- `E` and `B` share the same code unit; the plasma current obeys CGS Ampère's
  law `∇×B = 4π J`, so the current used in Ohm's law is `J = ∇×B/(4π)`.
- The Hall/pressure denominator is the charge density `ρ_q = e n_e`
  (`nodePlasma[...].iRho_` from `PicParticles::sum_moments`), matching the CGS
  denominators exactly — no extra `q/m` factor is needed.

## How to Run

```bash
python3 tests/validate_tests.py --test=hybrid_whistler
```

## Expected Results & Success Criteria

1. **Stability** — clean exit, no NaN/Inf, bounded `Eb`/`Epart`. The `4π` bug
   (Hall/resistive terms using raw `∇×B` instead of `∇×B/(4π)`) or a bad
   `nHallSubcycle` makes the whistler blow up.
2. **Hall (whistler) dispersion** — for parallel propagation in a hybrid
   (kinetic-ion) model the correct branch is the **finite-frequency
   (ion-cyclotron-bounded) whistler**
   ```
   ω / Ω_i = (k d_i)^2 / (1 + (k d_i)^2)
   ```
   with `Ω_i = q B / m` and `d_i = v_A / Ω_i`. This reduces to the cold
   Hall-MHD form `ω/Ω_i = (k d_i)^2` in the `ω ≪ Ω_i` (`k d_i ≪ 1`) limit and
   is bounded by `Ω_i` as `k d_i → ∞`. It is the physically correct branch for
   the seeded `n=1` mode, where `k d_i ≈ 1` (`ω/Ω_i ≈ 0.5`, well into the
   near-cyclotron regime): the cold form overpredicts `ω` by ~40% there.
   The automated `validate_tests.py` check tracks the seeded mode's
   circularly-polarized phase over the time-resolved `.out` frames and compares
   the fitted `|ω|` to `(k d_i)^2/(1+(k d_i)^2)` with a ~25% relative tolerance —
   tight enough to be a real check of the Hall term (a missing/factor-4π Hall
   current shifts `ω` by the same factor, far outside the window) while still
   absorbing the residual kinetic-ion / finite-Larmor corrections at
   `k d_i ≈ 1`. For the `n=1` box mode the measured `|ω|/Ω_i ≈ 0.54` matches
   `(k d_i)^2/(1+(k d_i)^2) = 0.49` to ~9% (phase fit `r² ≈ 0.90`); the cold
   Hall-MHD value would have been `(k d_i)^2 = 0.96`, a ~40% overshoot.

## Sub-cycling (`#HALLSUBCYCLE`)

`nHallSubcycle` sub-cycles the Faraday/Hall B update by `subDt = dt/nHallSubcycle`.
At this test's `dt = 0.02` the advance is already RK4-stable without
sub-cycling (the compact curl operators annihilate the Nyquist mode, so the
grid-scale Hall mode is far less stiff than the naive continuum estimate). The
value `8` is a conservative margin that still exercises the sub-cycling
machinery.
