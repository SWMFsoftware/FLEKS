# Planar Solar Wind Discontinuity Boundary Condition Test (`bc_discontinuity`)

Verifies the spatio-temporal `#DISCONTINUITY` boundary condition (`include/DiscontinuityBC.h`, `src/DiscontinuityBC.cpp`) across `inflow` field and particle boundaries.

## Variants

1. **`PARAM.in` (`BC_DISCONTINUITY` — 2D Tilted Tangential Discontinuity)**:
   - A $45^\circ$ tilted tangential discontinuity ($\hat{\mathbf{n}} = (1, 1, 0)/\sqrt{2}$) enters a 2D domain $[0, 8] \times [0, 8]$ across both the `x-lo` and `y-lo` `inflow` boundaries.
   - State 1 (`#UNIFORMSTATE`, with `#INFLOW` omitted so State 1 moments are automatically populated from `#UNIFORMSTATE`): $n_1 = 5\text{ cm}^{-3}$, $\mathbf{u}_1 = (100, 100, 0)\text{ km/s}$, $T_1 = 10^5\text{ K}$, $\mathbf{B}_1 = (0, 0, +5)\text{ nT}$.
   - State 2 (`#DISCONTINUITY` `td`): $n_2 = 10\text{ cm}^{-3}$, $\mathbf{u}_2 = (100, 100, 0)\text{ km/s}$, $T_2 = 96110.182\text{ K}$, $\mathbf{B}_2 = (0, 0, -3)\text{ nT}$, satisfying exact total thermal + magnetic pressure balance $P_{\text{tot},1} = P_{\text{tot},2}$ and internal finite-width temperature compensation.

2. **`PARAM.in.rd` (`BC_DISCONTINUITY (RD)` — Rotational Discontinuity)**:
   - A rotational discontinuity ($\hat{\mathbf{n}} = (1, 0, 0)$) with normal field $B_x = 3\text{ nT}$ and a $180^\circ$ rotation of the tangential field $B_y: +4\text{ nT} \to -4\text{ nT}$ at constant $|\mathbf{B}| = 5\text{ nT}$ and exact Walén velocity jump $\Delta u_y = \Delta B_y / \sqrt{\mu_0 \rho_{m,1}}$.

## Running

```bash
python3 tests/validate_tests.py --test=bc_discontinuity
```

