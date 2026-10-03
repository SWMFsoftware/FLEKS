# AMR Plasma Thermal Equilibrium Test

This test verifies the stability and accuracy of a uniform plasma thermal equilibrium across a stationary coarse-fine AMR interface using both the Full PIC solver (`PARAM.in`) and the Hybrid PIC solver (`PARAM.in.hybrid`).

## Physical Setup

- **Hybrid PIC (`PARAM.in.hybrid`)**:
  - Kinetic ions ($q = 1, m = 1$, $n_0 = 1.0\ \text{amu/cc}$, $T_i = 10\ \text{eV}$) and neutralizing massless fluid electrons ($T_e = 10\ \text{eV}$, isothermal $\gamma = 1.0$).
  - Perpendicular uniform guide field: $B_{z0} = 1.0\ \text{nT}$.
  - Tests the generalized Ohm's law and ambipolar field $\mathbf{E}_{\text{ambi}} = -\nabla P_e / (e n_e)$ across AMR interfaces.
- **Full PIC (`PARAM.in`)**:
  - Kinetic ions ($m_i = 1.0, q_i = 1$) and kinetic electrons ($m_e = 0.04, q_e = -1$, $m_i/m_e = 25$).
  - Charge neutral ($n_e = n_i = 1.0$) and equal initial temperature ($T_i = T_e = 10\ \text{eV}$).
  - Semi-implicit field solver with GMRES.

## Mesh Refinement Geometry

- **Domain**: 2D periodic box with $x \in [-16, 16]$, $y \in [-8, 8]$, $z \in [-1, 1]$.
- **Base Level (Level 0)**: $32 \times 16 \times 1$ cells ($\Delta x_0 = \Delta y_0 = 1.0$).
- **Refinement Region (Level 1)**: Central box (`center_slab`) covering $x \in [-8, 8]$, $y \in [-4, 4]$ with refinement ratio 2:
  - Fine grid cell size: $\Delta x_1 = \Delta y_1 = 0.5$.
  - Internal coarse-fine interfaces at $x = \pm 8.0$ and $y = \pm 4.0$.
  - Periodic domain boundaries at $x = \pm 16.0, y = \pm 8.0$ remain purely coarse-to-coarse, cleanly decoupling the AMR interfaces from periodic boundary wrapping.

## Validation Criteria

The test validator (`validate.py`) enforces:
1. **Energy Conservation & Stability**:
   - Total energy drift $|\Delta E_{\text{tot}}| / E_0 < 0.1\%$ across the run.
   - Magnetic energy $E_b$ remains bounded without unphysical growth.
2. **AMR Hierarchy Detection**:
   - Both coarse ($\Delta x \approx 1.0$) and fine ($\Delta x \approx 0.5$) spacings are identified.
   - Fine cells are strictly confined to the central slab $|x| \le 8.0$.
3. **Equilibrium Preservation**:
   - Mean density profile $\bar{\rho}(x)$ across the interfaces remains flat within PIC shot noise ($\max |\bar{\rho} - 1.0| < 0.10$).
   - Midplane reflection symmetry: $|\bar{\rho}(x) - \bar{\rho}(-x)| < 0.08$.

## Running the Tests

```bash
# Run both variants:
python3 tests/validate_tests.py --test=amr_equilibrium -v

# Run only Full PIC:
python3 tests/validate_tests.py --test=amr_equilibrium.full -v

# Run only Hybrid PIC:
python3 tests/validate_tests.py --test=amr_equilibrium.hybrid -v
```
