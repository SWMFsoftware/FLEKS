# Hybrid PIC AMR Thermal Equilibrium Test

This test verifies the stability and accuracy of a uniform plasma thermal equilibrium across a stationary coarse-fine AMR interface using the Hybrid PIC solver (`#HYBRIDPIC T`, `#SOLVEEM F`).

## Physical Setup

- **Plasma**: Single kinetic ion species ($q = 1, m = 1$) in a uniform Maxwellian thermal distribution with neutralizing massless fluid electrons ($n_e = n_i$).
- **Equilibrium**:
  - Uniform ion density: $n_0 = 1.0\ \text{amu/cc}$.
  - Ion temperature: $T_i = 10\ \text{eV}$ ($116045\ \text{K}$).
  - Electron temperature: $T_e = 10\ \text{eV}$ (isothermal $\gamma = 1.0$).
  - Perpendicular uniform guide field: $B_{z0} = 1.0\ \text{nT}$ ($B_0 \sim 1.0$ in normalized code units).
  - Bulk velocity: $\mathbf{u}_i = 0$, current $\mathbf{J} = 0$, ambipolar field $\mathbf{E}_{\text{ambi}} = -\nabla P_e / (e n_e) = 0$.
- **Analytic state**: Exact steady state ($\partial/\partial t = 0$). Macroparticle thermal motion tests kinetic fluxes and discrete noise across the AMR interface.

## Mesh Refinement Geometry

- **Domain**: 2D periodic box with $x \in [-16, 16]$, $y \in [-8, 8]$, $z \in [-1, 1]$.
- **Base Level (Level 0)**: $32 \times 16 \times 1$ cells ($\Delta x_0 = \Delta y_0 = 1.0$).
- **Refinement Region (Level 1)**: Central slab (`center_slab`) covering $x \in [-8, 8]$, $y \in [-8, 8]$ with refinement ratio 2:
  - Fine grid cell size: $\Delta x_1 = \Delta y_1 = 0.5$.
  - Two symmetric, internal coarse-fine interfaces at $x = -8.0$ and $x = +8.0$.
  - Periodic domain boundaries at $x = \pm 16.0$ remain purely coarse-to-coarse, cleanly decoupling the AMR interface from periodic boundary wrapping.

## Validation Criteria

The test validator (`validate.py`) enforces:
1. **Energy Conservation & Stability**:
   - Total energy drift $|\Delta E_{\text{tot}}| / E_0 < 2\%$ across the run.
   - Magnetic energy $E_b$ remains bounded without unphysical Hall/whistler growth.
2. **AMR Hierarchy Detection**:
   - Both coarse ($\Delta x \approx 1.0$) and fine ($\Delta x \approx 0.5$) spacings are identified.
   - Fine cells are strictly confined to the central slab $|x| \le 8.0$.
3. **Equilibrium Preservation**:
   - Mean density profile $\bar{\rho}(x)$ across the interfaces remains flat within PIC shot noise ($\max |\bar{\rho} - 1.0| < 0.10$).
   - Midplane reflection symmetry: $|\bar{\rho}(x) - \bar{\rho}(-x)| < 0.08$.

## Running the Test

```bash
python3 tests/validate_tests.py --test=amr_equilibrium -v
```
