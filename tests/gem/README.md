# GEM Challenge & Asymmetric Magnetic Reconnection Test

This directory provides 2D collisionless magnetic reconnection test cases for
standalone FLEKS using the semi-implicit Full-PIC solver (Maxwell/GMRES) with
kinetic ions and electrons ($m_i/m_e = 25$).

## Available Configurations

- **`PARAM.in`** (Default):
  Classic GEM challenge benchmark (Birn et al. 2001) with a single Harris current
  sheet in a conducting box ($L_x = 25.6\,d_i, L_y = 12.8\,d_i$) and modal magnetic
  perturbation ($B_0 = 1.0, \lambda_0 = 0.5, \text{apert} = 0.1$).

- **`PARAM.in.asym`**:
  2D asymmetric magnetic reconnection benchmark with double current sheets in
  a periodic domain ($L_x = 64.0\,d_i, L_y = 28.0\,d_i$) and Gaussian perturbation
  seeds:
  - Asymmetric fields and temperatures: $B_1 = 1.0, B_2 = 2.0$, $T_1 = 1.33, T_2 = 3.33$.
  - Current sheet half-thickness: $\lambda_0 = 1.0\,d_i$.
  - Periodic boundary conditions in $x, y, z$.

## Running the Test

Run the automated test runner:

```bash
# Run both GEM and GEM (ASYM) variants:
python3 tests/validate_tests.py --test=gem

# Run only the asymmetric variant:
python3 tests/validate_tests.py --test=gem.asym

# Run with 2 MPI processes:
python3 tests/validate_tests.py --test=gem -n 2
```
