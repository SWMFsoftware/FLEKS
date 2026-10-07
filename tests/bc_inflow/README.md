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
last declared one.

A uniform two-species plasma in which every species shares the same bulk
velocity is the exact steady-state solution, so the inflow/outflow pair has to
keep it uniform **species by species** — see Validation.

## Running

```bash
python3 tests/validate_tests.py --test=bc_inflow
```

## Validation

`validate.py` reads the expected per-species state back from the deck
(`#PLASMA`, `#UNIFORMSTATE`), so the checks follow the deck instead of
hard-coded numbers.

* **Energy log** — `Etot`/`Ee`/`Eb`/`Epart` finite; `Epart` bounded (open BCs
  inject and remove particles every step, so that tolerance is loose); `Eb`
  within 20% of its initial value (the `Bx` guide field is uniform); `Ee`
  negligible, i.e. no spurious `E` build-up at the faces.
* **Per species, from the `.out` frames** — `<uxS{i}>` conserved across the
  domain; inflow-side (first third) `rhoS{i}` within 10% and temperature
  `T_i = pS{i} * m_i / rhoS{i}` within 15% of the first frame; inflow-side
  pressure isotropic with `Pzz > 0`.
* **Fields** — `Bx` uniform and unchanged; no spurious mean `Ex`/`Ey`/`Ez`.
* **Cross-species ratios** — `rhoS1/rhoS0 = n1*m1 / (n0*m0)` (3.2 here) and
  `T1/T0` (4.0 here).

The last item is the discriminating one: both ratios break by `O(m1/m0)` if the
upstream state is not per-species, or if the injected thermal speed ignores the
species mass. Dropping the second `#INFLOW` block, for instance, gives
`rhoS1/rhoS0 = 14.2` and `T1/T0 = 1.11` and fails the test.
