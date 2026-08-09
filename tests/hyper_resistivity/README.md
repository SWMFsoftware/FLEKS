# Hybrid PIC Hyper-Resistivity (fourth-order) Test

Exercises the fourth-order hyper-resistive term `eta_h * nabla^2 J` (equivalently
`-(eta_h/4*pi) * nabla x (nabla^2 B)`) in the Hybrid PIC generalized Ohm's law, on
top of the full solver (convection + Hall + `eta J` + `grad P_e`). See
`../ohm/README.md` for the shared geometry / normalization and the
success criteria inherited from that test.

## Validation

1. **Stability** — FLEKS exits cleanly; `Eb`/`Epart` finite and bounded with the
   fourth-order term active (no NaN/Inf, no checkerboard blow-up).
2. **Damping of grid-scale noise** — the hyper term damps the highest resolved
   modes (`k -> k_Nyquist`), so late-time `|B_⊥|` and the magnetic-energy drift
   stay below the no-hyper-resistivity baseline.

## Parameters

- `#HYPERRESISTIVITY`: `etaHyperSI` (SI `[m^4/s]`), `etaHyperMode` (`"si"` or
  `"grid"`), `etaHyperCh` (the dimensionless `C_h` for grid mode).
- `"grid"` (default here) is robust to resolution/time-step changes:
  `eta_h = 4*pi * C_h * dx_min^4 / dt_sub`. `"si"` uses a direct physical value.
- `#MINIMUMDENSITY`: minimum `rho` for the `1/rho` factors; `<= 0` = auto
  (`1e-6 * electronDensity0`).

## Dispersion validation

For a clean quantitative check use `etaHyperMode = "si"` with a known `eta_h` and
measure the decay rate of an eigenmode:
- Resistive: `gamma = eta * k^2` (second order).
- Hyper-resistive: `gamma = eta_h * k^4` (fourth order) — verify the `k^4` scaling
  by repeating at two wavenumbers.
- Checkerboard test: seed a `(-1)^i` Nyquist-mode B perturbation; with the compact
  current it must decay, while a collocated-curl implementation would leave it
  untouched — this is the direct demonstration of why Phase 1 (compact current) is
  a prerequisite for the hyper-resistive term.
