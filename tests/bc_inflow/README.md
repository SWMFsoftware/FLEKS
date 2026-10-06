# Inflow / Outflow Open-Boundary Hybrid Test

## Physics

A uniform magnetized hybrid plasma (kinetic ions + massless fluid electrons)
streams along `+x` through a 1D domain. The `-x` face is an **inflow**
boundary (`BC::inflow`) and the `+x` face is an **outflow** boundary
(`BC::outflow`).

The plasma has **two kinetic ion species**, so the per-species form of the
`#INFLOW` command is exercised:

| Species | m [m_p] | q [e] | n [/cc] | T [K] | `#UNIFORMSTATE` rho | `#INFLOW` rho |
|---|---|---|---|---|---|---|
| 0: solar-wind H+ | 1 | +1 | 5 | 10000 | 5.0 | 5.0 |
| 1: heavy minor ion | 16 | +1 | 1 | 40000 | 1.0 | 1.0 |

`rho` is a **number density [1/cc]** in both commands, so the heavy species
reads 1.0 in both places (its mass density is m·n = 16 amu/cc). Only the first
`#INFLOW` block is mandatory; a species without a block of its own reuses the
last declared one, so dropping the second block injects the heavy ion at the
*proton* density and temperature.

A uniform two-species streaming plasma in which every species shares the same
bulk velocity is the exact steady-state solution, so this test validates that
the inflow/outflow pair keeps a uniform state uniform **species by species**:

* Upstream densities `rhoS0` / `rhoS1` and bulk velocities `uxS0` / `uxS1` are
  preserved at the inflow face.
* Each species keeps its own temperature, i.e. the injected thermal speed is
  `sqrt(T/m)` with the *species* mass: `T1/T0` must stay at the prescribed
  factor 4 (a proton-mass `vth` for the heavy ion would give 64).
* The guide field `Bx` stays uniform and unchanged.
* No spurious electric field develops.
* Ion kinetic energy `Epart` and magnetic energy `Eb` stay finite and bounded.

## Running

```bash
python3 tests/validate_tests.py --test=bc_inflow
```

## Validation

`validate.py` checks:

1. **Energy log** (`validate_log`): `Etot`/`Ee`/`Eb`/`Epart` finite; `Epart`
   bounded (open BCs inject/remove particles every step, so the ratio
   tolerance is loose); `Eb` within ~20% of its initial value (uniform `Bx`); `Ee`
   negligible (no spurious E build-up).
2. **Plot output** (`validate_plot`), for every species in the output:
   * `<uxS{i}>` conserved across the domain,
   * inflow-side (first third of the domain) `rhoS{i}` preserved to ~10%,
   * inflow-side temperature `T_i = pS{i} * m_i / rhoS{i}` preserved to ~15%,
   * inflow-side pressure isotropic (`Pxx` ~ `Pyy` ~ `Pzz`, `Pzz` > 0),
   * guide field `Bx` uniform and unchanged,
   * no spurious mean `Ex`/`Ey`/`Ez`.
3. **Cross-species ratios**, read back from the deck (`#PLASMA` /
   `#UNIFORMSTATE`) so the checks follow the deck:
   * `rhoS1/rhoS0` equals the prescribed mass-density ratio
     `n1*m1 / (n0*m0)` = 3.2 (a single broadcast `#INFLOW` block gives 16),
   * `T1/T0` equals the prescribed temperature ratio (4.0).

Both ratio checks are the discriminating ones: they fail by a factor of order
`m1/m0` if the upstream state is not per-species or if the thermal speed
ignores the species mass.
