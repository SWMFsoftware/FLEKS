# Hybrid PIC Convection-Term Wave Test (advection of a transverse wave)

This test validates the **convection term** `E = -U_i × B` of the generalized
Ohm's law in the FLEKS hybrid solver (kinetic ions, massless fluid electrons),
in **isolation** from every other term. It is the wave-propagation companion of
the static [`hybrid_freestream`](../hybrid_freestream/README.md)
baseline: instead of merely checking that a uniform state stays uniform, it
checks that the convection term produces the **correct dynamics** — rigid
advection of a transverse magnetic wave at the bulk-flow speed.

## Physics & Configuration

With the Hall term OFF (`#HALLTERM F`), resistivity `0.0`, and electron
pressure `0.0`, the only E-field source is `E = -U_i × B`, so Faraday's law
reduces to the ideal-advection equation for a uniform bulk flow `U`:

```
dB/dt = ∇×(U×B) = -U·∇B
```

A **true 1D** periodic domain (32×1×1 cells, `Lx = 6.4` code units) carries:

- a **circularly polarized transverse wave** `By = B1 cos(kx x)`,
  `Bz = B1 sin(kx x)` with `kx = 2π/Lx` and `B1 = 0.2·Bx0`, seeded by
  `#TESTCASE ConvectionWave` (`Pic::fill_convection_wave`, which reuses
  `fill_hybrid_wave(0.2)`) with **no velocity perturbation**;
- a **uniform bulk flow** `ux = 10 km/s` (`U_code = 0.2`) from `#UNIFORMSTATE`.

### Why the guide field is deliberately weak

`Bx0 = 5.2195e-11 T` (`B_code ≈ 0.01`), 100× below the usual `v_A ≈ 1`
scaling of the other hybrid tests. This makes the field **dynamically
passive**: `v_A_code ≈ 0.01 ≪ U_code = 0.2` and `Ω_i·T_code ≈ 0.05 ≪ 1`, so
the kinetic ions cannot react on the run time scale. A δB-only seed on a
*strong* guide field is **not** rigidly advected — with no matching δu it
splits into two counter-propagating Alfvén/ion-cyclotron branches
(`δu = ∓δB/√(μ0ρ)` eigenmodes) and the pure-advection phase prediction fails.
Weakening `B` freezes the ion response and isolates the convection term
exactly, with residual Alfvén contamination only
`O(1 − cos(kx·v_A·T_code)) ≈ 0.1%` in amplitude.

### Time units (important)

`#TIMESTEPPING dt`, `#SAVEPLOT dt` and `#STOP TimeMax` are in **SI seconds**,
while wave phase speeds are per **code time unit**; the conversion is
`tNorm = lNormSI/uNormSI = 1.0e5/5.0e4 = 2.0 s`. The run lasts **one
advection period**: the wavelength `λ = Lx·lNormSI = 6.4·1.0e5 m = 640 km`, so
`T_period = λ/u_x = 640 km / 10 km/s = 64 s`; `TimeMax = 64 s` therefore
corresponds to `T_code = 32.0`, i.e. a full `-2π` phase turn (the pattern
returns to its initial shape). Plots are written every `16 s` (`#SAVEPLOT dt`),
i.e. at **quarter periods** `t = 0, 16, 32, 48, 64 s`, so the `k = 1` phase
advances by exactly `-π/2` between consecutive frames.
(Get the SI↔code time conversion wrong and the wave *appears* to move at
exactly half speed — the factor is precisely `tNorm = 2`.)

## Expected Results & Success Criteria

At a fixed position the complex transverse field `C = By + i·Bz` advances in
phase by

```
Δφ = -kx · U_code · T_code = -(2π/6.4) · 0.2 · 32.0 = -2π rad
```

so over the full run the transverse pattern makes one complete revolution and
returns to its initial shape, with **no amplitude change** (advection neither
damps nor grows the wave).

The automated check (`validate_tests.py :: _check_convection_advection`)
extracts the 1D `By, Bz` x-profile from each `z=0` plot frame and takes the
spatial DFT `k = 1` mode of `C`, then applies **two complementary checks**:

1. **Rate check** — over the two *consecutive* frames with the shortest time
   gap (one quarter period, `Δφ_expected = -π/2`), the measured phase shift
   must match `-kx·U_code·Δt_code` within `max(0.15 rad, 10%)`, converting the
   frame separation from SI seconds to code time via `tNorm`. The
   quarter-period baseline keeps the shift well-resolved and far from any `2π`
   wrap.
2. **Return-to-start check** — between the *first and last* frames the wave
   has advected exactly one wavelength, so the wrapped phase shift must be
   `≈ 0` (expected total `-2π`): the endpoint pattern must coincide with the
   initial one. This closes the loop the rate check alone leaves open — any
   per-quarter rate error accumulates 4× here.

Both checks require the transverse amplitude ratio in `[0.7, 1.3]` (advection
neither damps nor grows the wave).

The standard hybrid stability gates (clean exit, no NaN/Inf, bounded energies,
particle-number conservation) also apply via `validate_hybrid`.

## Running

```bash
python3 tests/validate_tests.py --test=hybrid_convection_wave
```
