# Floating Body Test (`#FLOATINGBODY`)

## Description

An absorbing sphere (`#BODY`, radius 1.5) in a closed box (reflecting walls)
of electron-proton plasma at rest, T = 1e5 K, mi/me = 16, full PIC with the
fields off (`#SOLVEEM F`), so the particles fly ballistically (3D, 24^3
cells, 2x2x2 particles per cell; 8 s, absorb 5 s). The electrons reach the sphere 4 times faster than the ions,
and in the closed box the net charge of the plasma is the opposite of the net
charge the body keeps.

* `absorb`: plain absorbing body (reference).
* `reemit`: the net absorbed electron charge of each step is re-emitted as
  thermal electrons over a sphere just above the body.

## Validation

* absorb: the plasma charges up (+10% of N(H+) in 5 s, accepted > 5%), equal
  to the absorbed electron minus ion charge.
* reemit: re-emitted = absorbed electron minus ion charge, plasma neutral
  within 1e-3.

Needs a 3D AMReX build.
