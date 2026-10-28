# Open Particle Boundary Test

## Description

A uniform plasma of H+ and O+ at T = 1e5 K (thermal speeds 28.7 and 7.2 km/s,
a ratio of 4:1 matching electron-ion disparity at m_i/m_e = 16) fills a 2D box
(24x24 cells, 6 s), hybrid PIC with B = 0 and Te = 0, so E = 0.

* `rest`: The plasma at rest, periodic in x and open on y- and y+ faces.
* `drift`: The plasma drifting along +x at 2 sigma(H+) = 8 sigma(O+), entering
  through an inflow face (#INFLOW) at x- and leaving through x+ with flanks at
  y- and y+, with the predictor-corrector step (`rk4-pc`).

With naive `outflow` (pure particle deletion / vacuum sink), the lighter species
H+ rapidly drains out of the box due to thermal motion, losing ~23.5% of its
energy in 6 s (0.4% per step), and heavily cavitating near the boundary faces.

With an `open` particle boundary, inward thermal half-flux is maintained at the
open faces so that mass and temperature are conserved.
