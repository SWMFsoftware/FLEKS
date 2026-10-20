# Explicit PIC Thermal-Plasma Energy Test

A 1D uniform, unmagnetized, thermal electron-ion plasma (m_i/m_e = 100,
256 cells, 64 particles per cell, periodic) advanced 500 steps with the
explicit full-PIC solver (`#FIELDSOLVER explicit`). Nothing drives the plasma,
so the total energy must stay constant. This is the reproducer of FLEKS issue
170, where the explicit solver heated it by a factor of ~560.

| Variant | File | Integrator | Pass condition |
|---------|------|------------|----------------|
| `thermal_plasma_leapfrog` | `PARAM.in.leapfrog` | leapfrog | max \|Etot/Etot0 - 1\| < 1% |
| `thermal_plasma_leapfrog_mag` | `PARAM.in.leapfrog_mag` | leapfrog, uniform B at 45 deg (electron w_ce = 0.3 w_pe) | max \|Etot/Etot0 - 1\| < 1% |
| `thermal_plasma_pc` | `PARAM.in.pc` | predictor-corrector | max \|Etot/Etot0 - 1\| < 2% |

Parameters in code units (c = 1, n = 1): omega_pe = 35.4, v_te = 0.354,
lambda_De = 0.010, dx = 2.8 lambda_De, c dt/dx = 0.36, omega_pe dt = 0.35.
The leapfrog variant keeps Etot within 1e-4 of its initial value; the
predictor-corrector variant loses about 0.6% in 500 steps.

```bash
python3 tests/validate_tests.py --test=thermal_plasma
```
