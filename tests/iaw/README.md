# Ion Acoustic Wave (IAW) test

A hybrid-PIC regression test for an ion acoustic wave seeded as a single
sinusoidal ion density perturbation in a uniform, Maxwellian, zero-bulk-velocity
proton plasma.

## Physics

The initial condition is

```
  n(x) = n0 * (1 + pert * sin(kx * x)),    kx = 2*pi*waveMode / Lx
  u_i  = 0   (Maxwellian about zero)
  E    = 0
  B    = B0  (constant guide field)
```

with `pert = 0.1` and one wavelength across the domain (`waveMode = 1`).  The
restoring force is the electron pressure-gradient term in the generalized
Ohm's law

```
  E = -(1/rho) ∇P_e  =  -(T_e / rho) ∇n_e ,       (isothermal, gamma = 1)
```

so the acoustic speed is `c_s = sqrt(T_e / m_i)`.  With the chosen
normalization (`u_norm = 5e4 m/s`) the parameters give `c_s ≈ v_A ≈ ω_ci ≈ 1`
in code units and `c_s / v_thi ≈ 4`, i.e. a weakly Landau-damped wave that
oscillates and decays over the run.

The guide field is aligned with the propagation direction (`k || B`, `B` along
`x`), which makes the mode a clean longitudinal ion-acoustic wave (the Hall and
convective `u × B` terms vanish for a parallel perturbation).  This is a minor,
physics-motivated deviation from setups that place `B` along `z`; the seeded
density perturbation still lives along `x`.

## Implementation

The density perturbation is seeded by the C++ `IonAcousticWave` test case (in
`Pic` / `Particles`): each loaded particle's weight is multiplied by
`(1 + pert * sin(kx * x))`, so the deposited ion density (charge/mass/pressure
moments) carries the sinusoidal profile while the bulk velocity stays uniform.

Activate it in `PARAM.in` with

```
#TESTCASE
testCase  IonAcousticWave
pert      0.1
waveMode  1
```

`pert` is the relative density-perturbation amplitude (must be `< 1`) and
`waveMode` is the number of wavelengths across the x-domain.

## Run

```
python3 tests/validate_tests.py --test=iaw
```

This builds and runs the test, then validates that the run is stable (no NaN;
mass conserved) and that the seeded density profile is a clean single-mode
sinusoid (`rhoS0` in the plot output).  Plot-file checks require `PostIDL.exe`
(`make PIDL`) and `PostProc.pl`; if those are unavailable the profile check is
skipped and only the stability checks are enforced.

## Reference

Modeled on the hybrid-VPIC-style wave test (`tests/hybrid_whistler`), with the
transverse-wave seeding (`HybridWave` test case) replaced by a longitudinal ion
density perturbation (`IonAcousticWave` test case) and the electron
pressure-gradient (`#ELECTRONTEMPERATURE`) providing the restoring force.
