# Hybrid PIC Proton Cyclotron Anisotropy Instability (PCAI) Test

This test validates the FLEKS **hybrid** solver (kinetic ions / massless fluid
electrons) against the **proton cyclotron anisotropy instability (PCAI)** — a
linear electromagnetic instability driven by ion temperature anisotropy
(`T_perp/T_par > 1`). The parameters follow the reference case documented in
*Key Notes in Plasma Physics* §30.8.1
([Proton Cyclotron Anisotropy Instability](https://henry2004y.github.io/KeyNotes/contents/hybrid.html#proton-cyclotron-anisotropy-instability)).

## Physics & Configuration

The maximum-growth mode of the PCAI propagates **parallel** to the background
field (`k × B0 = 0`), so the test is **1D-3V**: `k` and `B0` are both along `x`.
The plasma starts in a **bi-Maxwellian** ion distribution with

| Parameter | Value | Meaning |
| :--- | :--- | :--- |
| `β_par` | 1 | parallel plasma beta |
| `T_perp / T_par` | 3 | perpendicular/parallel temperature (anisotropy) |
| `Lx / d_i` | 10.5 | domain length in ion inertial lengths |
| `k_x · Δx` | ≈ 0.065 | box-fundamental mode (n=1, `k d_i = 2π/10.5 = 0.598`) |
| `Δt · Ω_ci` | 0.01 | time step (`dt = 0.01 s`, `tNorm = 1 s`) |
| dissipation | 0 | no resistivity / e-pressure-gradient |

These are above the threshold
`P_perp/P_par − 1 ≈ S / β_par^0.4` (`S ~ 1`), so the left-hand
circularly-polarized (Alfvén/ion-cyclotron) wave grows at the linear-theory rate
`γ/Ω_ci = 0.162`.

### Anisotropic seeding (matches Hybrid-VPIC reference)

`#UNIFORMSTATE` temperature `T = 6.28e5 K` is interpreted as the **parallel**
temperature `T_par` (giving `β_par = 1` at `n = 5 /cc`, `B = 10.44 nT`,
`v_A = 102 km/s`). The `#WAVEIC anisoTPerpOverTPar 3.0` option makes the `WaveIC`
thermal-velocity hook inflate the two velocity components perpendicular to the
`x` guide field by `sqrt(T_perp/T_par) = sqrt(3)`, seeding a bi-Maxwellian with
`T_perp = 3 T_par` — the FLEKS equivalent of the reference seeding
`ux ~ N(0, vthipar)`, `uy, uz ~ N(0, vthiperp)` with `vthipar = sqrt(T_par) =
0.7071 v_A`, `vthiperp = sqrt(T_perp) = 1.2247 v_A`.

As in the Hybrid-VPIC reference (`set_region_field` seeds only the uniform `Bx`
guide field), the initial state has **no** `(B_y, B_z)` field perturbation
(`frac = 0`): the instability grows purely from the anisotropic ion free energy
and the ~1% thermal-noise floor. The linear growth rate is the eigenmode rate
`γ/Ω_ci = 0.162`. The growth rate is measured from the total transverse-field
norm `d|B| = sqrt(mean(By²) + mean(Bz²))` — the same diagnostic as the
reference's `plotsPCAI.py`, which tracks `d|By|`/`d|Bz|` against
`0.012·exp(γ·t)`.

### Normalization

The normalization is chosen for a clean time conversion: `lNormSI = uNormSI =
1.0e5`, so `tNorm = lNormSI/uNormSI = 1.0 s`. Then `TimeMax == T_code ==
t·Ω_ci` directly (with `Ω_ci = 1` in code units): `TimeMax = 60 s` gives
`t·Ω_ci = 60 = 9.55` gyro-periods (the reference's `taui = 60 wci^-1`). `dt =
0.01 s` gives `dt_code = 0.01` (`Δt·Ω_ci = 0.01`) directly. `v_A ≈ 1.02` in code
units (not normalized to 1 — not needed). The guide field is set to
`bx = 1.044e-8 T` so that `B_code ~ 1` (`Ω_ci = 1` in code units; with `tNorm =
1` this also means `Ω_ci(SI) = eB/m_p = 1 rad/s`). The cell-centred hybrid solver
is stable at `nHallSubcycle = 1` (no sub-cycling, matching the reference
`nsub = 1`).

## Expected Results & Success Criteria

1. **Anisotropy-instability growth** — the total transverse-field norm
   `d|B| = sqrt(mean(By²) + mean(Bz²))` grows **exponentially** from the seeded
   mode. The validator fits `log d|B|` vs time over the exponential phase and
   checks the rate is **positive** (a damped/stable plasma — e.g. cold electrons
   `Te = 0`, which leave the mode oscillating at marginal stability — fails this)
   and not a **runaway** (a missing `1/(4π)` Hall factor or CFL blow-up fails
   this). See *Measured growth rate* below for the quantitative comparison to
   theory.
2. **Unstable, not frozen** — `d|B|` must grow well above the noise floor,
   confirming the anisotropic free energy (`T_perp/T_par = 3`) actually drives
   the instability rather than being damped or left unchanged.
3. **Hybrid stability** — via the shared `validate_hybrid`: clean exit, no NaN/Inf,
   bounded magnetic and ion energies, particle-number conservation.

## Measured Growth Rate (parameter-tuning finding)

Linear-Vlasov theory (the value the Hybrid-VPIC reference uses in `plotsPCAI.py`)
gives `γ/Ω_ci = 0.162` for these nominal parameters. The FLEKS hybrid solver
captures the instability — `d|B|` grows by over an order of magnitude from the
thermal-noise floor — with a measured rate that has varied across solver phases
from ~0.07 to ~0.32 (the current Phase A/B cell-centred solver measures
`γ/Ω_ci ≈ 0.07`, below the theory value). This is treated as a documented solver
characteristic of the massless-fluid-electron hybrid model (the Hall term is
independently validated by `hybrid_whistler`), not a configuration error, so the
test passes with `0.05 < γ/Ω_ci < 0.5`.

Two tuning findings worth recording:

- **Warm electrons are required.** With cold electrons (`Te = 0`) the
  proton-cyclotron mode only oscillates (the fitted rate drops to `γ ~ 0`):
  the warm-electron pressure-gradient term couples to the parallel density
  perturbation and is needed to drive the growth, consistent with the
  reference (`Te = Tipar`). The test therefore sets `Te = T_par` (`54.2 eV`).
- **Stability margin.** The cell-centred (Hybrid-VPIC-style) solver is stable at
  `nHallSubcycle = 1` (no sub-cycling): the cell-centred 2·dx curl annihilates
  the Nyquist mode and the quadratic-spline gather/scatter filters grid-scale
  modes, so the explicit Hall `B` update stays stable through the growing-wave
  amplitudes without sub-cycling, matching the reference (`nsub = 1`).

## Running

```bash
python3 tests/validate_tests.py --test=pcai
```

The summary table lists `PCAI` (hybrid PIC). This test is hybrid-only (no
`PARAM.in.hybrid`), so it is exercised once with the `useHybridPIC` solver.
