# Optical-depth attenuation of exospheric photoionization

Validates the `#OPTICALDEPTH` EUV attenuation model of BATSRUS `ModUserMars`:

```
tau(r) = sum_c n_c(r) * sigma_c * H_c(r)      vertical optical depth
mu     = max(x . solarDir / r, cosSzaFloor)   cosine of the solar zenith angle
A(r)   = exp(-tau/mu)                         plain model
A(r)   = exp(-tau_* chap(r/H_*, cosSZA))      chapmanFunction = T
S_i(r) = n_i(r) * nu0_i * A(r)                production of neutral i
```

`nu0` is the **unattenuated** rate, i.e. exactly what BATSRUS calls `Rate_I`
for the selected `#SOLARCONDITION`; the solarmax values `7.3e-7` (CO2),
`2.734e-7` (O) and `8.59e-8 s^-1` (H) are identical to the `RateDim_I` table
in `ModUserMars.f90`.

## Deck

`PARAM.in` (plain `exp(-tau/mu)`) and `PARAM.in.chapman` (Chapman function).

* 2D, 64x64 cells over `[-3, 3] rPlanet`, periodic, 1 s at `dt = 0.1 s`.
* Two exosphere components: H (`sigma = 0`, `H0 = 3000 km`) and O
  (`sigma = 8e-21 m^2`, `H0 = 1200 km`), so `tau(rPlanet) = 0.96` and the
  optical depth falls below 0.06 at `r = 2 rPlanet` — the whole transition is
  inside the box.
* Because H does not absorb but is still attenuated by the O column, the H+
  channel verifies that `tau` sums over **all** neutral components.
* The plain deck disables every floor (`minProduction = 0`, `tauFloor = 0`,
  `tauCutoff = 1e30`) so the attenuation is exactly `exp(-tau/mu)` on the
  dayside and exactly zero on the nightside.

## What the validator checks

1. the O+ source profile along the subsolar line follows `exp(-tau/mu)`
   (relative agreement better than 25% at 1.2, 1.5, 1.8, 2.2, 2.6 rPlanet),
2. the H+ profile does too,
3. the nightside production is suppressed below 10% of the dayside,
4. the frame contains only finite values.

The source signal is measured as the difference between the final and the
initial frame, which cancels the uniform background plasma and makes the check
independent of the plot density unit.

Typical measured/model ratios are 0.99-1.22, and the nightside-to-dayside
ratio is `< 1e-6` for the plain deck.

## Running

```sh
python3 tests/validate_tests.py --test=opticaldepth -v
```

## Note on the Chapman function

The Smith & Smith (1972) fit used by BATSRUS becomes **negative** deep on the
nightside (sin SZA -> 0), where the slant column is in fact optically thick.
`exp(-tau*chap)` then exceeds unity, i.e. it would amplify the EUV flux.
BATSRUS never evaluates that regime because its own optical depth exceeds the
`13.8` guard first; with a thin column the guard does not trigger, so FLEKS
treats a non-positive Chapman value as opaque (`minProduction` floor).  The
validator mirrors this convention.
