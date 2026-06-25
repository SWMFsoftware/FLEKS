# Photoionization Source Test (H + O, Chamberlain Neutral Profile)

## Description

This standalone test verifies the `#PHOTOIONIZATION` command for photoionization
of analytical neutral exosphere density profiles with **two species** (H and O)
using the **Chamberlain model**. The photoionization source injects heavy ions
(species 1) into the PIC simulation domain.

## Physics & Solver Setup

- **Geometry & Boundaries**: Quasi-1D periodic grid (32 cells in X, 1 in Y, Z),
  domain spans [-1000, 1000] m.
- **Exosphere Model**: Chamberlain profile with 2 neutral components:
  | Component | n0 [m^-3] | H0 [m] | Description                  |
  |-----------|-----------|--------|------------------------------|
  | H         | 1.0×10^10 | 1000   | Light species, large scale height |
  | O         | 3.0×10^10 | 100    | Heavy species, small scale height |

  Density: \( n(r) = n_0 \exp\left[-H_0\left(\frac{1}{R_p} - \frac{1}{r}\right)\right] \) for \( r \geq R_p \).
  Planet radius \( R_p = 500 \) m.

- **Ionization**: Photoionization only, enabled via `#PHOTOIONIZATION` command,
  rate \( \nu_{\text{photo}} = 1.0 \times 10^{-6} \) s^-1 per component.
  Source mass rate per component:
  \( S_{\rho,i} = n_i(r) \cdot \nu_{\text{photo},i} \cdot (R_p/r)^2 \).
- **Other Processes**: `#ELECTRONIMPACT` and `#CHARGEEXCHANGE` commands are not
  present, so those source processes are skipped.

- **Plasma Species**:
  | Species | Mass [amu] | Charge [e] | Role                       |
  |---------|------------|------------|----------------------------|
  | 0 (H+)  | 1.0        | +1         | Background light ions      |
  | 1 (O+)  | 16.0       | +1         | Heavy ions (gets exosphere source) |

- **Electromagnetic Fields**: Enabled (`solveEM = T`) with guiding field
  \( B_x = 5.0 \times 10^{-9} \) T.

## Expected Results

1. **Particle Count Growth**: Species 1 (O+) particle count increases monotonically
   over time due to continuous photoionization of the exosphere neutrals.

2. **Total Energy Growth**: Kinetic energy of species 1 increases as more
   particles are injected.

3. **No Crashes**: The simulation completes all 10 iterations without errors.

## Photoionization Rate Profile

At planet surface (r = 500 m):
- \( n_H = 1.0 \times 10^{10} \) m^-3, source rate: \( 1.0 \times 10^4 \) m^-3 s^-1
- \( n_O = 3.0 \times 10^{10} \) m^-3, source rate: \( 3.0 \times 10^4 \) m^-3 s^-1
- \( S_{\rho,\text{total}} = 4.0 \times 10^4 \) m^-3 s^-1

At domain boundary (r ≈ 1414 m):
- \( n_H \approx 2.75 \times 10^9 \) m^-3, source: \( \approx 2.75 \times 10^3 \) m^-3 s^-1
- \( n_O \approx 2.64 \times 10^{10} \) m^-3, source: \( \approx 2.64 \times 10^4 \) m^-3 s^-1

## Running

From the FLEKS root directory (requires compiled `bin/FLEKS.exe`):

```bash
python3 tests/test_standalone/validate_tests.py
```

Or manually:

```bash
mkdir -p run_test/PC/plots run_test/PC/restartOUT
ln -sf ../bin/FLEKS.exe run_test/FLEKS.exe
cp tests/test_standalone/photoionization/PARAM.in run_test/PARAM.in
cd run_test && ./FLEKS.exe
```
