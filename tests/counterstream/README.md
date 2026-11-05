# Counter-Streaming Two-Ion Hybrid Solver Benchmark (Phobos Foreshock Setup)

This test adapts the simulation setup described in:
> *Backstreaming-Ion Effects on the Phobos Wake in an Idealized Martian Foreshock: Wake Refilling and Dual-Wake Formation*

## Physical Context
The simulation investigates the interaction of a spherical unmagnetized absorbing obstacle (Phobos, $R = 11\text{ km}$) with two kinetic $H^+$ populations in an idealized Martian foreshock:
1. **Solar Wind Protons (Species 0)**:
   - Bulk density $n_{\text{sw}} = 3.0\text{ cm}^{-3}$
   - Bulk velocity $V_{\text{sw}} = (-400, 0, 0)\text{ km/s}$ (anti-parallel to $B_0$)
   - Thermal temperature $T_{\text{sw}} = 7.45\text{ eV} \approx 86,440\text{ K}$
   - Injected dynamically across the right boundary ($\text{Hi-}x$, $+x$)

2. **Counter-Streaming Beam Protons (Species 1)**:
   - Density ratio $\alpha = n_b / n_{\text{sw}} = 0.125 \implies n_b = 0.375\text{ cm}^{-3}$
   - Bulk velocity $V_b = (+400, 0, 0)\text{ km/s}$ (parallel to $B_0$)
   - Thermal temperature $T_b = 29.8\text{ eV} \approx 345,760\text{ K}$
   - Injected dynamically across the left boundary ($\text{Lo-}x$, $-x$)

3. **Background Magnetic Field**:
   - Uniform guide field aligned with $+x$: $B_0 = (3.0, 0, 0)\text{ nT}$.

4. **Electron Fluid**:
   - Massless neutralizing fluid ($n_e = n_0 + n_1$) with isothermal closure ($\gamma_e = 1.0, T_e = 7.45\text{ eV}$).
   - Low-density Ohm's law floor $\rho_{\text{min}} = 0.15\text{ cm}^{-3}$.

## 2D Automated Test Setup
- **Domain**: $x \in [-16, 16]$ code units ($32 R$), $y \in [-7, 7]$ code units ($14 R$), $z \in [-0.5, 0.5]$ (1 cell).
- **Resolution**: $64 \times 28 \times 1$ cells ($\Delta x = \Delta y = 0.5$ code units).
- **Obstacle**: Center at $(0, 0, 0)$, radius $R = 1.0$ code unit (diameter 4 cells across in $y$).
- **Boundaries**:
  - $x$: Field $\text{Lo-}x$ is outflow (float), $\text{Hi-}x$ is fixed inflow ($B_0 = 3\text{ nT}$).
  - $x$ Particles: Per-species `#PARTICLEBOXBOUNDARY` blocks:
    - Species 0: `outflow` at $\text{Lo-}x$, `inflow` at $\text{Hi-}x$.
    - Species 1: `inflow` at $\text{Lo-}x$, `outflow` at $\text{Hi-}x$.
  - $y, z$: Periodic for both fields and particles.

## Full 3D Production Scaling
To scale from this 2D regression test to the full 3D run from the paper:
- `nCellX = 128`, `nCellY = 56`, `nCellZ = 56`
- `zMin = -7.0`, `zMax = 7.0`
- `dx = dy = dz = 0.25` code units ($2.75\text{ km}$, 8 cells across Phobos diameter $2 R = 22\text{ km}$).
- `#PERIODICITY`: $y$ and $z$ periodic.

