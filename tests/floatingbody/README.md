# Floating Body Test (`#FLOATINGBODY`)

## Description

An absorbing sphere (`#BODY`, radius 1.5) in a closed box (reflecting walls)
of electron-proton plasma at rest, T = 1e5 K, mi/me = 16, full PIC with the
fields off (`#SOLVEEM F`), so the particles fly ballistically (3D, 24^3
cells, 20 s). The electrons reach the sphere 4 times faster than the ions,
and in the closed box the net charge of the plasma is the opposite of the net
charge the body keeps.

* `absorb`: plain absorbing body (reference).
* `sheath`: subgrid sheath: an electron entering the sphere is absorbed only
  if its normal energy exceeds the body potential, else reflected; the
  potential relaxes (1 s) toward zero net absorbed current.
* `reemit`: the net absorbed electron charge of each step is re-emitted as
  thermal electrons over a sphere just above the body.

## Validation

* absorb: the plasma charges up (+34% of N(H+) in 20 s), equal to the absorbed
  electron minus ion charge.
* sheath: after the potential has relaxed (5 s), zero net absorbed current
  (within 10%) and a flat plasma charge (0.5% of N(H+)); the potential
  (measured 7.8 - 9 V) is of order kT/e ln 4 = 12 V, lower because the box
  depletes (accepted 0.7 - 1.6 kT/e).
* reemit: re-emitted = absorbed electron minus ion charge, plasma neutral
  within 1e-3.

Needs a 3D AMReX build.
