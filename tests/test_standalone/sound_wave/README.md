# 1D Adiabatic Hydrodynamic (Sound) Wave Test

## Description
This is a standalone test designed to verify the propagation of a 1D linear adiabatic hydrodynamic sound wave in FLEKS. The test initializes a spatial sine perturbation on the grid, which naturally loads self-consistent perturbations in the ion particles and electromagnetic fields, and simulates wave propagation over exactly one wave period.

## Physics & Solver Setup
- **Geometry & Boundaries**: 1D periodic grid (64 cells in $X$, 1 in $Y$ and $Z$).
- **Background State (Normalized)**: 
  - Ion density $\rho_0 = 1.0$, pressure $p_0 = 0.6$ ($\gamma = 5/3$), and zero bulk velocity.
  - Guiding magnetic field $\vec{B}_0 = (0.0, 0.0, 0.0)$.
- **Perturbations (Eigenfunctions)**:
  - Perturbation amplitude $\text{amp} = 10^{-4}$, wavelength $\lambda = 1.0$.
  - $\delta \rho = \text{amp} \sin(kx)$
  - $\delta u_x = -\text{amp} \sin(kx)$
  - $\delta p = \text{amp} \sin(kx)$
- **Solver Setup**: Treated in hybrid PIC mode with kinetic ion particles ($nSpecies = 1$) and background electron fluid. Field solver is enabled (`solveEM = T`).
- **Wave Propagation**: The phase speed is $v_{\text{phase}} = \sqrt{\gamma p_0 / \rho_0} = 1.0$, meaning the wave travels exactly one period at $T_{\text{max}} = 1.0$.

## Expected Results
- **Particle Count Conservation**: The total number of particles must remain exactly constant (deviation $< 10^{-5}$).
- **Kinetic Energy Stability**: The total kinetic energy must remain stable over the entire wave period (variation $< 10\%$), as the linear wave propagates without excessive numerical dissipation or grid heating.
