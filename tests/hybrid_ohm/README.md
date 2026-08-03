# Hybrid PIC Full Generalized Ohm's Law Standalone Test

This test exercises the **complete** Hybrid PIC (kinetic ions, massless fluid
electrons) field solver in standalone FLEKS: all four terms of the generalized
Ohm's law are active in `Pic::assemble_ohm_E` (see `docs/Algorithm.tex`):

1. **Convection** — `E = -U_i × B` (always active).
2. **Hall** — `(J × B) / (e n_e)`, sub-cycled (`#HALLSUBCYCLE`).
3. **Resistive** — `η J`, enabled by `#RESISTIVITY > 0`.
4. **Electron pressure gradient** — `∇P_e / (e n_e)`, enabled by
   `#ELECTRONTEMPERATURE > 0`.

This is the companion to [`hybrid_whistler`](../hybrid_whistler/README.md),
which disables terms 3 and 4. Both share the same geometry, ion-scale
normalization (`lNormSI ≈ d_i`, `uNormSI ≈ v_A` so `d_i ≈ v_A ≈ 1`), and the
circularly-polarized `HybridWave` seed. See the `hybrid_whistler` README for
the unit-convention discussion; the Hall-CFL sub-cycling rationale is
documented in the `PARAM.in` comments.

The solver-selection commands (`#SOLVEEM`, `#HYBRIDPIC`, `#RESISTIVITY`,
`#ELECTRONTEMPERATURE`, `#HALLSUBCYCLE`) and the conversion of the SI/eV
inputs to code units are documented as comments directly in
[`PARAM.in`](PARAM.in).

## Purpose & Success Criteria

This test is a **stability / no-blow-up regression check** for the *full*
solver. The `HybridWave` seed only perturbs `B` and the ion velocity
(eigenmode IC), so the pressure-gradient term has no large coherent source;
the test instead verifies that enabling the resistive and pressure-gradient
terms does not destabilize the coupled field advance. The automated check
reuses `validate_hybrid`:

1. **Stability** — FLEKS exits cleanly; `Eb` and `Epart` finite and bounded
   (no NaN/Inf, no runaway from an active `η J` or `∇P_e` term).
2. **Seeded wavelength** — early-time DFT of `By` still peaks at mode `n=1`
   (the wave is correctly seeded, same as `hybrid_whistler`).
3. **Bounded amplitude** — late-time `|B_⊥|` stays below 10×`B₀`
   (resistivity damps; a broken `∇P_e` or `η` term would blow up instead).
