# Open Particle Boundary Test (`outflowthermal`)

## Description

A uniform plasma of H+ and O+ at T = 1e5 K (thermal speeds 28.7 and 7.2 km/s,
a ratio of 4 like electrons and ions of equal temperature with m_i/m_e = 16)
fills a 3D box (24^3 cells), hybrid PIC with B = 0 and Te = 0, so E = 0.

* `rest`: the plasma at rest, `outflowthermal` on all six faces.
* `drift`: the plasma drifting along +x at 2 sigma(H+) = 8 sigma(O+), entering
  through an `inflow` face (#INFLOW) at x- and leaving through
  `outflowthermal` faces elsewhere (outflow at x+, parallel flow at the
  flanks).

`outflowthermal` re-injects at each face cell the inward half-flux of the
species' local Maxwellian, measured one node layer inside the face, with only
the outward part of the normal drift (the re-entered particles bias the
measured drift inward, which otherwise feeds back and runs away).

## Validation

The mass of each species stays within 3% of its initial value over 10 s
(measured: rest 1.1%, drift 1.7%). With plain `outflow`, H+ loses 0.87% per
step, about half of its mass over the run.

Needs a 3D AMReX build.
