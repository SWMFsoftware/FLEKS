# Langmuir Test (full-PIC particle-field coupling)

A 1D uniform, thermal electron-ion plasma (m_i/m_e = 100, 256 cells,
64 particles per cell, periodic) with nothing driving it, advanced 800 steps.
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
| `langmuir_explicit_pc` | `PARAM.in.explicit_pc` | explicit, predictor-corrector | < 2% | each < 6%, mean < 3% |
| `langmuir_explicit_mag` | `PARAM.in.explicit_mag` | explicit, leapfrog, uniform B at 45 deg (electron w_ce = 0.3 w_pe) | < 1% | not checked (the oblique field splits the branch) |
| `langmuir_implicit` | `PARAM.in.implicit` | semi-implicit (reference) | Etot/Etot0 in [0.90, 1.01] | each < 6%, mean < 3% |

Each spectral peak must also exceed 100 times the median of its spectrum.
The leapfrog scatter of about 2% is the frequency resolution of the 8 s
record; the predictor-corrector and implicit schemes run 1-2% low, more at
high k (numerical dispersion), and the implicit scheme (theta = 0.51) loses
about 6% of the energy over the run.

Parameters in code units (c = 1, n = 1): omega_pe = 35.4, v_te = 0.354,
lambda_De = 0.010, dx = 2.8 lambda_De, c dt/dx = 0.36, omega_pe dt = 0.35.
Each variant takes about 15 s on one core.

Reference results (2D build): explicit +0.15% mean / 2.6% max, energy within
3e-4; explicit_sub -0.07% / 1.3%; explicit_pc -1.8% / 4.7%, energy -0.9%;
implicit -1.3% / 3.9%, energy -6.3%.

```bash
python3 tests/validate_tests.py --test=langmuir
python3 tests/validate_tests.py --test=langmuir.explicit   # one variant
```
