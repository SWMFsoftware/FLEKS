# Langmuir Test (full-PIC particle-field coupling)

A 1D uniform, thermal electron-ion plasma (m_i/m_e = 100, 256 cells,
64 particles per cell, periodic) with nothing driving it, advanced 1600 steps.
The thermal noise excites Langmuir oscillations at every wavenumber, which
makes this a direct check of how the particles' current drives E:

1. **Energy**: the total (field + particle) energy must stay constant.
2. **Langmuir branch**: the omega-k spectrum of E_x(x, t) must peak at the
   Bohm-Gross frequency omega^2 = omega_pe^2 (1 + 3 k^2 lambda_De^2) for the
   12 lowest modes (k lambda_De <= 0.11).

A current with the wrong sign or size shifts the branch; one centered at the
wrong time makes it grow. FLEKS issue 170 (the explicit solver used the mass
flux, centered at t^n) heated this plasma by a factor of ~560 in 500 steps.

| Variant | File | Solver | Energy bound | Langmuir branch |
|---------|------|--------|--------------|-----------------|
| `langmuir_explicit` | `PARAM.in.explicit` | explicit, leapfrog | max \|Etot/Etot0 - 1\| < 1% | each mode < 4%, mean < 1% |
| `langmuir_explicit_sub` | `PARAM.in.explicit_sub` | explicit, leapfrog, 2 field subcycles | < 1% | each < 4%, mean < 1% |
| `langmuir_explicit_pc` | `PARAM.in.explicit_pc` | explicit, predictor-corrector | < 3% | each < 6%, mean < 3% |
| `langmuir_explicit_mag` | `PARAM.in.explicit_mag` | explicit, leapfrog, uniform B at 45 deg (electron w_ce = 0.3 w_pe) | < 1% | not checked (the oblique field splits the branch) |
| `langmuir_implicit` | `PARAM.in.implicit` | semi-implicit (reference) | Etot/Etot0 in [0.80, 1.01] | each < 6%, mean < 3% |

Each spectral peak must also exceed 100 times the median of its spectrum. A
mode's frequency is the power-weighted centroid of its spectrum within 8% of
omega_pe around the peak. With the 16 s record one noise realization
scatters single modes by about 1%; an 8 s record, or the bare peak bin,
scatters them by up to 4%. The predictor-corrector and implicit schemes run
1-2% low, more at high k (numerical dispersion), and both lose energy
steadily: -1.5% and -11% over the run (the implicit scheme uses
theta = 0.51).

Parameters in code units (c = 1, n = 1): omega_pe = 35.4, v_te = 0.354,
lambda_De = 0.010, dx = 2.8 lambda_De, c dt/dx = 0.36, omega_pe dt = 0.35.
Each variant takes about 25 s on one core.

Reference results (2D build; frequency mean / max deviation from
Bohm-Gross, energy at the end):

| Variant | Langmuir branch | Etot/Etot0 |
|---------|-----------------|------------|
| explicit | +0.36% / 1.25% | 1.0000 (within 4e-4 throughout) |
| explicit_sub | +0.52% / 1.44% | 0.9998 |
| explicit_pc | -1.12% / 2.85% | 0.985 |
| explicit_mag | - | 1.0001 |
| implicit | -1.31% / 2.04% | 0.887 |

```bash
python3 tests/validate_tests.py --test=langmuir
python3 tests/validate_tests.py --test=langmuir.explicit   # one variant
```
