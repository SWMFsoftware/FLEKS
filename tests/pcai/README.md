# Hybrid PIC Proton Cyclotron Anisotropy Instability (PCAI) Test

Validates the FLEKS **hybrid** solver (kinetic ions / massless fluid electrons)
against the **proton cyclotron anisotropy instability (PCAI)**, a linear
instability driven by ion temperature anisotropy `T_perp/T_par > 1`. Parameters
follow the reference case in *Key Notes in Plasma Physics* §30.8.1
([Proton Cyclotron Anisotropy Instability](https://henry2004y.github.io/KeyNotes/contents/hybrid.html#proton-cyclotron-anisotropy-instability)).

## Configuration

The maximum-growth mode propagates **parallel** to the background field, so the
test is **1D-3V** (`k ∥ B0 ∥ x`) with a bi-Maxwellian ion distribution:

| Parameter | Value | Meaning |
| :--- | :--- | :--- |
| `β_par` | 1 | parallel plasma beta |
| `T_perp / T_par` | 3 | temperature anisotropy |
| `Lx / d_i` | 10.5 | box-fundamental mode `k d_i = 2π/10.5 ≈ 0.6` |
| `Δt · Ω_ci` | 0.01 | time step |
| `nParticle` | 600 | macro-particles per cell (kept low so the test runs ~0.5 min on 1 core) |
| `TimeMax` / `dn` | 40 / 500 | run `t·Ω_ci = 40`; ~8 output frames |
| `nHallSubcycle` | 1 | no Hall sub-cycling (stable with the cell-centred solver) |

The anisotropy is seeded as a bi-Maxwellian (`T_perp = 3 T_par` via the
`anisoTPerpOverTPar` option); as in the Hybrid-VPIC reference, there is **no
initial `B_y`/`B_z` perturbation** (`frac = 0`), so the left-hand
circularly-polarized wave grows from the thermal-noise floor at the linear-theory
rate `γ/Ω_ci = 0.162`. Normalization is `lNormSI = uNormSI = 1e5` (`tNorm = 1 s`,
so `TimeMax` directly equals `t·Ω_ci`), with `bx = 1.044e-8 T` for `Ω_ci = 1` in
code units.

## Validation

The growth rate is measured from the total transverse-field norm
`d|B| = sqrt(mean(By²) + mean(Bz²))` vs time. The validator:

- Fits `log(d|B|)` vs `t` over the exponential window `[0.2, 0.85]×amp_max` and
  requires `0.06 ≤ γ/Ω_ci ≤ 0.35` (must clearly grow — not damped/stable — and
  not run away).
- Checks clean exit, no NaN/Inf, bounded magnetic and ion energies, and strict
  ion-energy conservation (`Epart0` ratio within `[0.5, 2.0]`).

The measured rate is near the linear-theory `γ/Ω_ci = 0.162` (the FLEKS hybrid
model under-predicts somewhat). The Hall term is independently validated by the
`hybrid_whistler` test.

## Running

```bash
python3 tests/validate_tests.py --test=pcai
```
