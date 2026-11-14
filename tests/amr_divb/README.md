# Semi-Implicit AMR Div(B) Preservation Test (`tests/amr_divb`)

Tests the divergence-free upwind dissipation of the semi-implicit (full PIC)
field solver across a 3D coarse-fine AMR interface under periodic boundary
conditions.

## Physical Setup

A drifting, magnetized electron-ion plasma traverses a central stationary
refined cube:
- **Plasma**: Drift velocity $\mathbf{u} = (-300, 0, 0)\ \text{km/s}$,
  density $n = 1.0\ \text{cm}^{-3}$, temperatures $T_e = T_i = 60{,}000\ \text{K}$,
  mass ratio $m_i / m_e = 16$.
- **Magnetic field**: $\mathbf{B} = (1.06, 1.06, 0.0)\ \text{nT}$ ($45^\circ$ to the flow).
- **Motional electric field**: $\mathbf{E} = -\mathbf{u} \times \mathbf{B} = (0, 0, 3.18 \times 10^{-4})\ \text{V/m}$.
- **Domain**: 3D periodic box $x, y, z \in [-2.0, 2.0]\ R_M$. Base mesh $16^3$ cells ($\Delta x_0 = 0.25\ R_M$).
- **AMR hierarchy**: Central cube $[-1.0, 1.0]^3\ R_M$ refined by $2\times$ (Level 1, $\Delta x_1 = 0.125\ R_M$).

## Variants

| File | Upwind Scheme | Cleaning | Expected $\max \|\nabla \cdot \mathbf{B}\|$ | Purpose |
| :--- | :--- | :--- | :--- | :--- |
| `PARAM.in` | Curl upwind (`typeUpwindB = curl`) | None | $< 10^{-10}$ (machine zero, $\sim 10^{-15}$) | Standard CI regression test for AMR solenoidal preservation |
| `PARAM.in.inflow` | Curl upwind (`typeUpwindB = curl`) | None | $< 10^{-10}$ (interior/interface) | Non-periodic inflow/outflow AMR regression test with `#PRESPLITTING` and `#RESAMPLING` verifying suppression of refinement-boundary $E$-field artifacts |
| `PARAM.in.scalar` | Legacy scalar diffusion (`typeUpwindB = scalar`) | Hyperbolic (`hypDecay = 0.1`) | $\sim 4 \times 10^{-2}$ | Negative control demonstrating coarse-fine interface artifact |
| `PARAM.in.uniform` | Curl upwind (`typeUpwindB = curl`) | None | $< 10^{-10}$ ($\sim 10^{-15}$) | Verifies divergence-free evolution on uniform mesh |

## Numerical Background

### Legacy Scalar Upwind
The legacy scheme applies cell-centered diffusive flux $\Delta \mathbf{B}_c = \sum_d (\mathcal{F}_{d+1/2} - \mathcal{F}_{d-1/2})$.
Because scalar diffusion does not commute with divergence, non-solenoidal fluxes accumulate at coarse-fine
boundaries where grid spacing changes by $2\times$. Even with hyperbolic cleaning, $\nabla \cdot \mathbf{B}$ spikes
at the refinement boundary to $\sim 4 \times 10^{-2}\ \text{nT}/R_M$.

### Curl Upwind
The curl upwind scheme computes dissipation as an effective electric field on nodes:
$$
\mathbf{E}_{\rm upwind} = \sum_d \frac{a_d}{2} \hat{\mathbf{e}}_d \times (\mathbf{B}^R_d - \mathbf{B}^L_d)
$$
Faraday's law advances $\mathbf{B}$ exclusively through the discrete curl $\mathbf{B}^{n+1} = \mathbf{B}^n - \Delta t\, \nabla \times (\mathbf{E} + \mathbf{E}_{\rm upwind})$.
Since $\mathrm{div}_{\rm node \to center} \circ \mathrm{curl}_{\rm node \to center} \equiv 0$ on each level, the update preserves
$\nabla \cdot \mathbf{B} = 0$ to machine precision without requiring any hyperbolic cleaning.

## Running

```bash
# Run standard AMR curl test:
python3 tests/validate_tests.py --test=amr_divb -v

# Run side-by-side comparison variants:
python3 tests/validate_tests.py --test=amr_divb.inflow -v
python3 tests/validate_tests.py --test=amr_divb.scalar -v
python3 tests/validate_tests.py --test=amr_divb.uniform -v
```
