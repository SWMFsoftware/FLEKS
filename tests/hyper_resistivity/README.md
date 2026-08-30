# Hybrid PIC Hyper-Resistivity Test

Quantitative test of the fourth-order hyper-resistive term `-eta_h * nabla^2 J`
(= `-(eta_h/4*pi) * curl(nabla^2 B)`) in the Hybrid PIC generalized Ohm's law.

## Setup

Frozen plasma: no macroparticles (`#PARTICLES 0 0 0`), so `rho = 0` and the
`rho > 0` guard short-circuits every other Ohm term (convection, Hall, `eta*J`,
`grad Pe`) — the hyper-resistive term is then the only thing driving `B`.  A
single circularly polarized transverse mode

    B = (Bx0, B1 cos kx, B1 sin kx),   k = 2*pi*waveMode/Lx

is seeded on the uniform guide field `Bx0` (mode 8 of 32 cells, so
`theta = k*dx = pi/2` with `dx = 0.2`; `frac = 0.5`).  The guide-field energy
is constant, hence

    Eb(t) - Eb_guide = A * exp(-2*gamma*t)

with the discrete rate of the term

    gamma = (eta_h/4*pi) * 4 sin^2(theta) * sin^2(theta/2) / dx^4

which follows from the two operators it is built from: the collocated `2*dx`
curl (symbol `i*sin(theta)/dx`) and the 3-point Laplacian (symbol
`-4 sin^2(theta/2)/dx^2`).  `etaHyperSI = 4.0e17` m^4/s converts to
`eta_h_code = 4*pi*etaHyperSI*Si2NoV*Si2NoL^3 = 5.03e-2`, giving `gamma = 5.0`
and `gamma*dt = 0.01` at `dt = 0.002`.

The rate peaks at `theta = acos(-1/3) = 1.911` (wavelength 3.29 cells,
`gamma = 2.37*(eta_h/4*pi)/dx^4`) and vanishes at the Nyquist mode
(`sin(pi) = 0`) — see Limitations.

## Validation

`validate.py` derives the expected rate from the PARAM.in values
(`etaHyperSI`, `lNormSI`/`uNormSI`, `dx`, `waveMode`, `dt`, with the RK4
amplification `-ln(R(-gamma*dt))/dt` folded in), fits `Eb(t)` (a single mode
decays exactly geometrically, so `rho = exp(-2*gamma*dt)` follows from three
consecutive log samples and the constant part from two) and asserts:

1. **dissipative** — `Eb` decreases and stays finite;
2. **single exponential** — residual `< 1e-6` of the amplitude.  A broken
   stencil or stale RK-stage ghosts show up here first (the ghost bug gave
   `2e-2`, a correct run gives `1.6e-12`);
3. **rate** — within `0.5%` of the analytic symbol (measured agreement `1e-10`).

Together these catch a disabled or mis-converted `eta_h` (exactly how the
SI→code conversion bug presented), a wrong sign, and any regression in the
stencils or in the RK trial-state ghost handling.  Verified to fail when
`etaHyperSI = 0`.

## Run

```bash
python3 tests/validate_tests.py --test=hyper_resistivity
```

Excluded from the 2D AMReX suite.

## Limitation

* **Grid mode is very weak.**  `eta_h = 4*pi*C_h*dx_min^4/dt` with `C_h = 1e-3`
  damps the fastest mode by only `2.4e-3` per step, so it does not visibly
  clean up PIC shot noise.  Use `"si"` mode with a physical value when the
  damping has to matter.  Grid mode is exercised (completion check only) by
  `../reconnection/PARAM.in.hybrid`.
* **The two-cell checkerboard is out of reach.**  The `2*dx` curl annihilates
  the Nyquist mode, so the term cannot damp the shortest wavelength on the grid
  — the usual reason for adding it.  An opt-in `nabla^4 B` stencil
  (`~ -eta_h k^4`, no Nyquist blind spot) would fix that; keep it opt-in for
  Hybrid-VPIC comparability and update the expected rate in `validate.py` if it
  is added.
* The Nyquist null space is **documented, not tested**: one seeded mode per run
  means it needs its own case, and the rate check already fails if the stencil
  changes.
* Full-solver coverage (hyper-resistivity together with convection, Hall,
  `eta*J`, `grad Pe` and kinetic ions) comes from
  `../reconnection/PARAM.in.hybrid`, where the late-time amplitude is
  noise driven and therefore only a completion check is possible.
