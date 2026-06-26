# Charge Exchange Source Test (H + O, Chamberlain Neutral Profile)

## Description

This standalone test verifies the `#CHARGEEXCHANGE` command for charge exchange
ionization of analytical neutral exosphere density profiles with **two species**
(H and O) using the **Chamberlain model**. The charge exchange source injects
heavy ions (species 1) into the PIC simulation domain.

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

- **Ionization**: Charge exchange, enabled via `#CHARGEEXCHANGE` command.
  Uses a constant cross-section model:
  \( \langle\sigma v\rangle = \sigma_{\text{CX}} \cdot u_{\text{rel}} \),
  where \( u_{\text{rel}} \) is the ion-neutral relative speed.
  | Component | sigmaCX [cm^2] |
  |-----------|---------------|
  | H         | 2.0×10^-15    |
  | O         | 1.0×10^-15    |
  Source mass rate per component:
  \( S_{\rho,i} = n_i(r) \cdot n_{\text{ion}} \cdot \sigma_{\text{CX},i} \cdot u_{\text{rel}} \).
- **Other Processes**: `#PHOTOIONIZATION` and `#ELECTRONIMPACT` commands are not
  present, so those source processes are skipped.

- **Plasma Species**:
  | Species | Mass [amu] | Charge [e] | Role                       |
  |---------|------------|------------|----------------------------|
  | 0 (H+)  | 1.0        | +1         | Background ions            |
  | 1 (O+)  | 16.0       | +1         | Ion species for CX (ux=100 km/s), gets exosphere source |

- **Electromagnetic Fields**: Enabled (`solveEM = T`) with guiding field
  \( B_x = 5.0 \times 10^{-9} \) T.

## Expected Results

1. **Particle Count Growth**: Species 1 (O+) particle count increases monotonically
   over time due to continuous charge exchange of the exosphere neutrals.

2. **Total Energy Growth**: Kinetic energy of species 1 increases as more
   particles are injected.

3. **No Crashes**: The simulation completes all 10 iterations without errors.

## Running

From the FLEKS root directory (requires compiled `bin/FLEKS.exe`):

```bash
python3 tests/validate_tests.py --test=chargeexchange
```
