# Hybrid PIC Hyper-Resistivity (fourth-order) Test

Quantitative test of the fourth-order hyper-resistive term `-eta_h * nabla^2 J`
(equivalently `-(eta_h/4*pi) * curl(nabla^2 B)`) in the Hybrid PIC generalized
Ohm's law.

## Setup

Frozen plasma: **no macroparticles** are loaded (`#PARTICLES 0 0 0`), so the ion
charge/mass density `rho` is zero everywhere and every other term of the Ohm's
law (convection, Hall, `eta*J`, `grad Pe`) is short-circuited by the `rho > 0`
guard.  The hyper-resistive term is then the only thing driving `B`.  A single
circularly polarized transverse mode

    B = (Bx0, B1 cos kx, B1 sin kx),   k = 2*pi*waveMode/Lx

is seeded on the uniform guide field `Bx0`, whose energy is constant, so the
magnetic energy obeys

    Eb(t) - Eb_guide = A * exp(-2*gamma*t)

with `gamma` the discrete rate below.  This is what makes a quantitative check
possible: in a particle run the late-time amplitude is noise driven (see
`../ohm/README.md`), so it is not a usable observable.

## Verified discrete damping rate

The term is built as `centerLapB = nabla^2 B`, then `centerHyperE =
curl(centerLapB)`, then `E -= (eta_h/4*pi) * centerHyperE`.  The collocated
`2*dx` curl has symbol `i*sin(theta)/dx` and the 3-point Laplacian has symbol
`-4 sin^2(theta/2)/dx^2`, so a transverse Fourier mode with `theta = k*dx`
decays at

```
gamma = (eta_h/4*pi) * 4 sin^2(theta) sin^2(theta/2) / dx^4
```

which peaks at `theta = acos(-1/3) = 1.911` (a wavelength of **3.29 cells**,
where `gamma = 2.37 * (eta_h/4*pi) / dx^4`) and **vanishes at the Nyquist
mode** (`sin(pi) = 0`), i.e. a two-cell checkerboard is not damped at all.

`validate.py` derives this rate from the PARAM.in values (`etaHyperSI`,
`lNormSI`/`uNormSI`, `dx`, `waveMode`, `dt`, with the RK4 amplification
`-ln(R(-gamma*dt))/dt` folded in) and asserts:

1. the term is dissipative (`Eb` decreases, stays finite);
2. the decay is a clean single exponential -- residual `< 1e-6` of the
   amplitude.  A broken stencil or stale RK-stage ghosts show up here long
   before they show up as a wrong rate (the ghost bug gave `2e-2`);
3. the fitted rate matches the analytic symbol within `0.5%`
   (measured agreement: `1e-10`, residual `1.6e-12`).

Together these catch a disabled or mis-converted `eta_h` (which is exactly how
the SI→code conversion bug showed up), a wrong sign, and any regression in the
stencils or in the RK trial-state ghost handling.

## Parameters

- `#HYPERRESISTIVITY`: `etaHyperSI` (SI `[m^4/s]`), `etaHyperMode` (`"si"` or
  `"grid"`), `etaHyperCh` (the dimensionless `C_h` for grid mode).
- `"si"` mode (used here): `eta_h_code = 4*pi * etaHyperSI * Si2NoV * Si2NoL^3`,
  converted once by `Pic::convert_resistivity()` after the normalization is
  final.  The value `4.0e17` gives `eta_h_code = 5.03e-2` and `gamma = 5.0` for
  the seeded mode 8 of 32 (`theta = pi/2`, `dx = 0.2`), i.e. `gamma*dt = 0.01`.
- `"grid"` mode: `eta_h = 4*pi * C_h * dx_min^4 / dt`.  Not exercised here (see
  below) but covered by `../reconnection/PARAM.in.hybrid`.

## Known limitations

* **Grid mode is very weak.**  With `C_h = 1e-3` the fastest-damped mode loses
  only `2.4e-3` of its amplitude per step, so over a typical run the term does
  not visibly clean up PIC shot noise.  Use `"si"` mode with a physical value
  if the damping has to matter.
* **Grid-scale noise is out of reach.**  The `2*dx` curl annihilates the
  two-cell checkerboard, so the term cannot damp the shortest wavelength the
  grid supports -- the very thing it is usually added for.  An opt-in
  `nabla^4 B` stencil (damping `~ -eta_h k^4`, no Nyquist blind spot) would fix
  that; it should stay opt-in to preserve comparability with Hybrid-VPIC, and
  the expected rate in `validate.py` would have to be updated.
* The Nyquist null space is documented above rather than tested: it needs its
  own run (one seeded mode per run), and the rate check already fails if the
  stencil changes.
* Full-solver coverage -- hyper-resistivity together with convection, Hall,
  `eta*J` and `grad Pe`, with kinetic ions -- comes from
  `../reconnection/PARAM.in.hybrid` (grid mode, completion check only).
