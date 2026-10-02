# Evolved electron pressure (`#ELECTRONPRESSURE`)

Tests the hybrid-PIC scalar electron pressure PDE that `#ELECTRONPRESSURE T`
switches on in place of the algebraic polytropic closure
`Pe = P0 (rho/rho0)^gamma`:

```
dPe/dt + div(u_e Pe) + (gamma_e - 1) Pe div(u_e)
    = (gamma_e - 1) [ div(kappa_hat . grad(Te)) + H_ei ]
```

with `u_e = U_i - J/(e n_e)`, `Te = Pe/n_e` and the Spitzer conductivity
`kappa = kappa0 Te^2.5`. `H_ei` is a documented hook and is zero for now.

All three decks are the `tests/iaw` ion acoustic wave (48x1x1, periodic,
2000 particles/cell, TimeMax 5) with a 20% single-wavelength density
perturbation along x, so `Pe` starts from the polytropic closure and
`Pe/rho^gamma` is spatially uniform.

## Decks

| deck          | kappa0 [W/(m K^3.5)] | B           | expected closure |
|---------------|----------------------|-------------|------------------|
| `adiabatic`   | 0                    | 0           | `Pe ~ rho^gamma` |
| `conduction`  | 1.0e-19 (reduced)    | 0           | `Pe ~ rho` (Te flattened) |
| `crossfield`  | 1.0e-19 (reduced)    | `by` = 1e-10 T | `Pe ~ rho^gamma` (conduction suppressed) |

The conduction coefficient is far below the physical Spitzer value
(9.2e-12) on purpose: with this normalization the physical value is so
conductive that Te would flatten within a single step. The chosen value still
violates the explicit parabolic limit `dt < dx^2/(2 chi)` by orders of
magnitude, which is what the point-implicit update is there to survive.

`fieldAlignedBMinSI` (1e-15 T) is the field strength below which the field
direction is treated as meaningless. It matters here: an unmagnetized run
still carries round-off level `B` (~1e-18 in code units), and normalising it
gives an arbitrary direction that would silently steer the heat flux away
from `grad(Te)` and switch the conduction off. Real fields (the `crossfield`
deck has `|B|` ~ 1e-2 in code units) are many orders of magnitude above the
cutoff.

## Validation

Both scale-free ratios

- `A = Pe / rho^gamma`  (flat for adiabatic electrons)
- `I = Pe / rho`        (flat for isothermal electrons)

are evaluated pointwise over the plot frames, and `validate.py` checks that
the expected-flat one has a much smaller relative spread (std/mean) than the
other. Because only the *relative* spread is compared, the check does not
depend on the output units of `Pe` or `rho`.

The wave is standing and Landau damped, and the conduction needs ~2 s to
flatten the initial adiabatic transient (the point-implicit update is stable
for any step but converges at one Jacobi sweep per step). The frame used is
therefore the one where the two closure hypotheses differ most,
`max |spreadA - spreadI|`, which skips both the transient and the frames
where the wave has damped away.

Run:

```
python3 tests/validate_tests.py --test=electron_pressure
```
