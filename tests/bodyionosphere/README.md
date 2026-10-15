# Ionospheric Body Boundary Test (`#BODYIONOSPHERE`)

## Description

A spherical body (`#BODY`, radius 1.5, absorbing) carries an ionospheric
boundary: before every push its one-cell-thick surface layer is refilled with a
Maxwellian (n0 = 1 /cc, T = 1e5 K, H+, uniform in solar zenith angle). The run
is hybrid PIC with B = 0 and Te = 0, so E = 0 and the ions fly ballistically
into vacuum (3D, 24^3 cells, outflow walls, 30 s).

Collisionless effusion from a sphere gives the escape rate
`Gamma = n0 vth/sqrt(2 pi) 4 pi R^2` and the density
`n(r) = n0/2 (1 - sqrt(1 - R^2/r^2))`.

## Validation

From the PIC log (`mBodyInject`, `mBodyAbsorb`, `Mpart0`):

* escape rate / box mass vs analytic, within 5% (independent of units);
* with n0 taken from one refill of the ghost layer, the escape rate and the
  box mass are 0-20% / 0-25% above the sphere values (the staircase body
  emits slightly more than a smooth sphere; measured 1.08 and 1.12).

Needs a 3D AMReX build.
