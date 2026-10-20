# Ion Beam-Beam Instability Test

## Description
This is a standalone test designed to verify the simulation of the ion beam-beam instability in FLEKS. The test initializes a uniform background plasma moving in the $-x$ direction and a beam of ions of a smaller density moving in the $+x$ direction, solving for electromagnetic field dynamics. The test directory holds two parameter files that exercise the same physics with the two field solvers:

- **`PARAM.in.implicit`** — semi-implicit full PIC: kinetic electrons and ions (`#FIELDSOLVER implicit`).
- **`PARAM.in.explicit`** — explicit full PIC: kinetic electrons and ions with explicit Maxwell advance (`#FIELDSOLVER explicit`).
- **`PARAM.in.hybrid`** — hybrid PIC: kinetic ions / massless fluid electrons (`#FIELDSOLVER hybrid`).

All variants share the identical grid, guide field, beam ratio/speeds, density, normalization and output, so the same FFT-based validation applies.

## Physics & Solver Setup
- **Geometry & Boundaries**: 1D periodic grid (64 cells in $X$, 1 in $Y$ and $Z$).
- **Physical Model**: A single charged species ($m = 1.0$, $q = 1.0$) with uniform background density $5.0\,\text{amu/cc}$ and background velocity $-400\,\text{km/s}$ ($u_x = -0.4$ in normalized units).
- **Beam setup**: $1\%$ of the particles have a beam velocity of $+400\,\text{km/s}$ ($u_x = +0.4$ in normalized units).
- **Electromagnetic Fields**:
  - Semi-implicit full PIC: `#FIELDSOLVER implicit` with a guiding magnetic field $B_x = 5\,\text{nT}$ ($5.0 \times 10^{-9}$ T).
  - Explicit full PIC: `#FIELDSOLVER explicit`, Boris-Maxwell leapfrog advance.
  - Hybrid PIC: `#FIELDSOLVER hybrid`, Hall term ON, resistivity and electron-pressure-gradient OFF so the instability grows undamped; the guide field is dynamically weak ($v_A \ll u_{\mathrm{beam}}$) so the explicit hybrid advance is not stiff at this resolution.
- **Purpose**: Verifies that the particle count and mean velocity are conserved, and tracks the growth of cyclotron waves in the transverse magnetic fields $B_y$ and $B_z$ from seed noise.

## Validation

The validator checks that the total particle count is conserved, the mean
velocity `MeanVx` is conserved within 1%, and that transverse `By`/`Bz`
waves grow from noise (power concentrated in low-order box-compatible modes). The hybrid variant additionally runs the shared hybrid energy checks.

## Running

All variants run together under the test name:

```bash
python3 tests/validate_tests.py --test=beam          # all non-expensive variants
```

To run a **single variant**, append its token to the test name:

```bash
python3 tests/validate_tests.py --test=beam.implicit  # semi-implicit full PIC
python3 tests/validate_tests.py --test=beam.explicit  # explicit full PIC
python3 tests/validate_tests.py --test=beam.hybrid    # hybrid PIC
```

The output table lists `BEAM (IMPLICIT)`, `BEAM (EXPLICIT)`, and `BEAM (HYBRID)` results.
