# Hybrid PIC Full Generalized Ohm's Law Standalone Test

This test exercises the **complete** Hybrid PIC (kinetic ions, massless fluid
electrons) field solver in standalone FLEKS: all four terms of the generalized
Ohm's law are active in `Pic::assemble_ohm_E` (see `Doc/Algorithm.tex`):

1. **Convection** — `E = -U_i × B` (always active).
2. **Hall** — `(J × B) / (e n_e)`, sub-cycled (`#BSUBCYCLE`).
3. **Resistive** — `η J`, enabled by `#RESISTIVITY > 0`.
4. **Electron pressure gradient** — `∇P_e / (e n_e)`, enabled by
   `#ELECTRONTEMPERATURE > 0`.

This is the companion to [`whistler`](../whistler/README.md),
which disables terms 3 and 4. Both share the same geometry, ion-scale
normalization (`lNormSI ≈ d_i`, `uNormSI ≈ v_A` so `d_i ≈ v_A ≈ 1`), and the
circularly-polarized `HybridWave` seed. See the `whistler` README for
the unit-convention discussion; the Hall-CFL sub-cycling rationale is
documented in the `PARAM.in` comments.

The solver-selection commands (`#SOLVEEM`, `#HYBRIDPIC`, `#RESISTIVITY`,
`#ELECTRONTEMPERATURE`, `#BSUBCYCLE`) and the conversion of the SI/eV
inputs to code units are documented as comments directly in
[`PARAM.in`](PARAM.in).

## Validation

Stability / no-blow-up regression check for the *full* solver: enabling the
resistive and electron-pressure-gradient terms must not destabilize the coupled
field advance (the `whistler` test leaves both disabled). The automated check
reuses `validate_hybrid` (`tests/_shared/hybrid.py`): stable exit, finite
`Eb`/`Epart`, seeded mode `n=1` at early time, and bounded late-time amplitude.

On top of that, `validate.py` adds a **resistive-damping bound**: the resistivity
is large enough (`#RESISTIVITY` = 2.0e9 m²/s, `etaCode` ~ 2.51) to hold back the
late-time transverse amplitude, so the late/early `max|B_perp|` growth factor
must stay below 1.9 (measured 1.65 with the term active, 2.16 with a negligible
`eta`). This is the check that catches a silently disabled resistive term — e.g.
a broken SI→code unit conversion that multiplies `eta` by zero.
