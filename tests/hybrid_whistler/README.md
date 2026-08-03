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
2. **Hall (whistler) dispersion** — the solver is Hall-MHD (no electron
   inertia), so for parallel propagation
   ```
   ω / Ω_i = (k d_i)^2
   ```
   with `Ω_i = q B / m` and `d_i = v_A / Ω_i`. The automated
   `validate_tests.py` check tracks the seeded mode's circularly-polarized
   phase over the time-resolved `.out` frames and compares the fitted `|ω|`
   to `(k d_i)^2` (50% tolerance, which catches the factor-4π error). At the
   `n=1` box mode (`k d_i ≈ 1`) the measured `|ω|/Ω_i ≈ 0.59` vs `0.96` for
   cold Hall-MHD — a kinetic-ion / near-ion-cyclotron shift, not the 4π error.
   The cleanest quantitative check is the `n=2` mode, where the measured
   `|ω|/Ω_i = 4.49` matches `(k d_i)^2 = 3.86` to ~16% (phase fit `r² = 1.00`).

## Sub-cycling (`#HALLSUBCYCLE`)

`nHallSubcycle` sub-cycles the Faraday/Hall B update by `subDt = dt/nHallSubcycle`.
At this test's `dt = 0.02` the advance is already RK4-stable without
sub-cycling (the compact curl operators annihilate the Nyquist mode, so the
grid-scale Hall mode is far less stiff than the naive continuum estimate). The
value `8` is a conservative margin that still exercises the sub-cycling
machinery.
