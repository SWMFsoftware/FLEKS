# 2D Cylindrical Expansion-into-Vacuum (Rarefaction) Benchmark

## Description

This benchmark evaluates inner body boundary conditions in the sub-grid sheath regime ($\Delta x \gg \lambda_{De}$).
An absorbing cylinder of radius $R_B = 1.5\,d_i$ sits at the origin of a $12\,d_i \times 12\,d_i$ 2D Cartesian domain ($[-6, 6]^2$).
The surrounding plasma is uniform, at rest ($u = 0$), quasi-neutral ($n_e = n_i = n_0$), with $T_e = T_i$ and mass ratio $m_i / m_e = 25$ ($v_{th,e} / v_{th,i} = 5$).
The interior of the cylinder enforces $\mathbf{E} = 0$ and $\partial \mathbf{B} / \partial t = 0$.

Because $\Delta x / \lambda_{De} = 10 \gg 1$, the macroscopic plasma outside $r > R_B$ should behave as a **quasi-neutral ambipolar rarefaction wave** propagating outward at the ion-acoustic speed $c_s \approx 0.016\,c$:

$$n(r, t) \approx n_0 \exp\left( -1 - \frac{r - R_B}{c_s t} \right) \quad \text{for } R_B < r \le R_B + c_s t$$

$$u_r(R_B^+) = -c_s, \quad n(R_B^+) \approx n_0 e^{-1} \approx 0.37\,n_0$$

$$E_r^{\text{amb}}(r, t) \approx -\frac{k_B T_e}{e c_s t}$$

With naive particle absorption (`#BODYBOUNDARY absorb`), the $5\times$ thermal velocity disparity causes rapid electron drainage, unphysical positive charging around the cylinder, and spurious electric fields that push ions outward ($u_{r,i} > 0$) rather than pulling them inward at $-c_s$.

## Benchmark Configurations

1. **`PARAM.in.unmagnetized` (Config 1: $\mathbf{B}_0 = 0$)**:
   - Pure radial cylindrical symmetry in 2D Cartesian coordinates.
   - Tests radial ambipolar balance, cut-cell staircasing, and convergence to the analytical rarefaction profile.

2. **`PARAM.in.magnetized` (Config 2: $\mathbf{B}_0 = B_0 \hat{\mathbf{x}}$)**:
   - Guide field $B_0 = 0.02$ code units ($v_A = 0.02\,c$, $\Omega_{ci} = 0.02\,\omega_{pi0}$, $\Omega_{ce} = 0.5\,\omega_{pi0}$).
   - Electron gyroradius $\rho_{ce} = 0.1\,\Delta x$ (magnetized, sub-grid).
   - Ion gyroradius $\rho_{ci} = 5\,\Delta x$ (unmagnetized / large gyroradius).
   - Exposes the anisotropic charging disparity:
     - **Poles** ($x \approx \pm R_B, y = 0$): Field-aligned electron loss along $\mathbf{B}_0 \parallel \hat{\mathbf{r}}$.
     - **Equator** ($x = 0, y \approx \pm R_B$): Electrons are cross-field trapped ($\rho_{ce} < \Delta x$), while unmagnetized ions impact the surface.

## Running

```bash
# Run with the standalone test validator
python3 tests/validate_tests.py --test=body_rarefaction
```
