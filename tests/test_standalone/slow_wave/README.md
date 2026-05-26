# 1D Adiabatic MHD Slow Wave Test

## Description
This is a standalone test designed to verify the propagation of a 1D linear adiabatic MHD slow magnetosonic wave in FLEKS. The test initializes self-consistent perturbations in both the kinetic ion particles and electromagnetic fields on a periodic grid, simulating wave propagation over exactly one wave period.

## Physics & Solver Setup
- **Geometry & Boundaries**: 1D periodic grid (64 cells in $X$, 1 in $Y$ and $Z$).
- **Background State (Normalized)**: 
  - Ion density $\rho_0 = 1.0$, pressure $p_0 = 0.6$ ($\gamma = 5/3$), bulk velocity $\vec{u}_0 = (0.0, 0.0, 0.0)$.
  - Guiding magnetic field $\vec{B}_0 = (1.0, \sqrt{2.0}, 0.5)$.
- **Perturbations (Eigenfunctions)**:
  - Perturbation amplitude $\text{amp} = 10^{-4}$, wavelength $\lambda = 1.0$.
  - $\delta \rho = 0.8944272\,\text{amp} \sin(kx)$
  - $\delta u_x = -0.4472136\,\text{amp} \sin(kx)$
  - $\delta u_y = -0.8432740\,\text{amp} \sin(kx)$
  - $\delta u_z = -0.2981424\,\text{amp} \sin(kx)$
  - $\delta p = 0.8944272\,\text{amp} \sin(kx)$
  - $\delta B_y = -0.4216370\,\text{amp} \sin(kx)$
  - $\delta B_z = -0.1490712\,\text{amp} \sin(kx)$
- **Solver Setup**: Treated in hybrid PIC mode with kinetic ion particles ($nSpecies = 1$) and background electron fluid. Field solver is enabled (`solveEM = T`).
- **Wave Propagation**: The slow magnetosonic speed is $v_{\text{slow}} = 0.5$. For wavelength $\lambda = 1.0$, the wave period is exactly $T_{\text{max}} = 2.0$.

## Expected Results
- **Particle Count Conservation**: The total number of particles must remain exactly constant (deviation $< 10^{-5}$).
- **Kinetic Energy Stability**: The total kinetic energy must remain stable over the entire wave period (variation $< 10\%$), as the linear wave propagates without excessive numerical dissipation or grid heating.
