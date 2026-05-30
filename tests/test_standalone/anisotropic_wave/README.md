# 1D Anisotropic MHD Wave Test

## Description
This is a standalone test designed to verify the propagation of a 1D linear anisotropic MHD wave perpendicular to the magnetic field (the fast wave in anisotropic MHD, also known as the CGL/anisotropic fast magnetosonic wave) as described in Daldorff et al. (2014) and Key Notes in Plasma Physics.

## Physics & Solver Setup
- **Geometry & Boundaries**: 1D periodic grid ($64$ cells in $X$, $1$ in $Y$ and $Z$, with domain length $L_x = 32.0$).
- **Background State (Normalized)**: 
  - Ion density $\rho_0 = 1.0$, pressure $p_0 = 4.5 \times 10^{-4}$ ($\gamma = 5/3$), bulk velocity $\vec{u}_0 = (0.0, 0.0, 0.0)$.
  - Guiding magnetic field $\vec{B}_0 = (0.0, 0.04, 0.0)$.
- **Perturbations (Eigenfunctions)**:
  - Perturbation amplitude $\delta = 0.1$ (moderately non-linear to be clearly visible above the PIC noise), wavelength $\lambda = 32.0$.
  - $\rho = \rho_0 [1 + \delta \sin(kx)]$
  - $u_x = V_f \delta \sin(kx)$
  - $p = p_0 [1 + \gamma \delta \sin(kx)]$
  - $p_\parallel = p_0 [1 + \delta \sin(kx)]$
  - $B_y = B_0 [1 + \delta \sin(kx)]$
  - Where the CGL perpendicular fast magnetosonic speed is:
    $$V_f = \sqrt{\frac{B_0^2 + 2 p_0}{\rho_0}} = 0.05$$
- **Solver Setup**: Full PIC mode with kinetic electrons (Species 0) and kinetic ions (Species 1). Field solver is enabled (`solveEM = T`).

## Expected Results
- **Particle Count Conservation**: The total number of particles must remain exactly constant (deviation $< 10^{-5}$).
- **Kinetic Energy Stability**: The total kinetic energy must remain stable (variation $< 20\%$), as the linear wave propagates without excessive numerical dissipation or grid heating.
