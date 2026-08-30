# Hybrid PIC Hyper-Resistivity (fourth-order) Test

Exercises the fourth-order hyper-resistive term `-eta_h * nabla^2 J`
(equivalently `-(eta_h/4*pi) * curl(nabla^2 B)`) in the Hybrid PIC generalized
Ohm's law.

Three variants, one per PARAM file (the runner executes each and reports it
separately):

| PARAM file | base_name | What it is | What it checks |
|---|---|---|---|
| `PARAM.in` | `hyper_resistivity` | Full solver: convection + Hall + `eta*J` + `grad P_e` + grid-mode hyper-resistivity, with particles | Completion only (integration smoke test) |
| `PARAM.in.damping` | `hyper_resistivity_damping` | Frozen plasma (no macroparticles) + `"si"` mode hyper-resistivity + one seeded transverse mode | **Quantitative**: the measured exponential decay rate against the analytic discrete symbol (2% tolerance) |
| `PARAM.in.damping_nyquist` | `hyper_resistivity_damping_nyquist` | Same, seeded at the Nyquist mode | The field must **not** decay (stencil null space) |

## Discrete damping rate (verified)

The term is built as `centerLapB = nabla^2 B` followed by `centerHyperE =
curl(centerLapB)` and `E -= (eta_h/4*pi) * centerHyperE`. With the collocated
`2*dx` curl (symbol `i*sin(theta)/dx`) and the 3-point Laplacian (symbol
`-4 sin^2(theta/2)/dx^2`), a transverse Fourier mode with `theta = k*dx` decays
at

```
gamma = (eta_h/4*pi) * 4 sin^2(theta) sin^2(theta/2) / dx^4
```

so `Eb(t) - Eb_guide = A * exp(-2*gamma*t)`. `PARAM.in.damping` measures that
rate and compares it with the expression above (the RK4 amplification
`-ln(R(-gamma*dt))/dt` is folded in). Measured agreement: 4.99983 vs an
analytic 5.0.

Consequences of this symbol, both pinned down by the tests:

* it peaks at `theta = acos(-1/3) = 1.911`, i.e. a wavelength of **3.29 cells**,
  where `gamma = 2.37 * (eta_h/4*pi) / dx^4`;
* it vanishes identically at the Nyquist mode (`sin(pi) = 0`), so a
  two-cell checkerboard is **not damped at all**.

## Parameters

- `#HYPERRESISTIVITY`: `etaHyperSI` (SI `[m^4/s]`), `etaHyperMode` (`"si"` or
  `"grid"`), `etaHyperCh` (the dimensionless `C_h` for grid mode).
- `"si"` mode: `eta_h_code = 4*pi * etaHyperSI * Si2NoV * Si2NoL^3`, converted
  by `Pic::convert_resistivity()` once the normalization is final.
- `"grid"` mode (used by the base `PARAM.in`): `eta_h = 4*pi * C_h * dx_min^4 /
  dt`, robust to resolution / time-step changes.
- `#MINIMUMDENSITY`: minimum `rho` for the `1/rho` factors; `<= 0` = auto
  (`1e-6 * electronDensity0`). Only relevant for the particle-carrying variant.

## Known limitations

* **Grid mode is very weak.** With `C_h = 1e-3` the fastest-damped mode loses
  only `2.4e-3` of its amplitude per step, so over a typical run the term does
  not visibly clean up PIC shot noise. Use `"si"` mode with a physical value if
  the damping has to matter.
* **Grid-scale noise is out of reach.** Because the `2*dx` curl annihilates the
  two-cell checkerboard (see above), the term cannot damp the shortest
  wavelength the grid supports -- the very thing it is usually added for. If
  that is needed, an opt-in `nabla^4 B` stencil (damping `~ -eta_h k^4` with no
  Nyquist blind spot) would be the fix; it should stay opt-in to preserve
  comparability with Hybrid-VPIC, and the two damping variants above would then
  need their expected rates updated.
* The base `PARAM.in` variant remains a smoke test: with kinetic ions the
  late-time state is noise driven, so its amplitude is not a reliable
  quantitative observable (the same reason `tests/ohm` needs a relative bound
  rather than an absolute decay rate).
