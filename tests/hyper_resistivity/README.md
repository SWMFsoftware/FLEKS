# Hybrid PIC Hyper-Resistivity (fourth-order) Test

Exercises the fourth-order hyper-resistive term `eta_h * nabla^2 J` (equivalently
`-(eta_h/4*pi) * nabla x (nabla^2 B)`) in the Hybrid PIC generalized Ohm's law, on
top of the full solver (convection + Hall + `eta J` + `grad P_e`). See
`../ohm/README.md` for the shared geometry / normalization.

**Smoke test only.** There is no `validate.py`, so the runner falls back to a
generic no-op check: it only confirms FLEKS runs and exits cleanly with the term
active. Nothing is checked quantitatively.

## Parameters

- `#HYPERRESISTIVITY`: `etaHyperSI` (SI `[m^4/s]`), `etaHyperMode` (`"si"` or
  `"grid"`), `etaHyperCh` (the dimensionless `C_h` for grid mode).
- `"grid"` (default here) is robust to resolution/time-step changes:
  `eta_h = 4*pi * C_h * dx_min^4 / dt_sub`. `"si"` uses a direct physical value.
- `#MINIMUMDENSITY`: minimum `rho` for the `1/rho` factors; `<= 0` = auto
  (`1e-6 * electronDensity0`).

## Implementation note

The term is built in `Pic::assemble_ohm_E` as a compact Laplacian of `B` followed
by the same `2*dx` curl that Faraday's law uses, mirroring Hybrid-VPIC
(`hyb_hypereta.cc`). That curl vanishes at the Nyquist wavenumber, so the term is
not expected to damp the pure `(-1)^i` checkerboard mode.

## TODO

1. Add a `validate.py`: check `Eb`/`Epart` stay finite and bounded, and compare
   against a no-hyper-resistivity baseline run (requires adding the paired case).
2. Add a quantitative check of the damping rate using `etaHyperMode = "si"` with
   `#RESISTIVITY` disabled. This likely needs its own case rather than reusing
   this `PARAM.in`, which is configured as a full-solver integration run.
3. Verify the high-`k` behaviour noted above, and if grid-scale damping is needed,
   consider an opt-in `nabla^4 B` stencil (keep it opt-in to preserve
   comparability with Hybrid-VPIC).
