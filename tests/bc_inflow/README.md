# Inflow / Outflow Open-Boundary Hybrid Test

## Physics

A uniform magnetized hybrid plasma (kinetic ions + massless fluid electrons)
streams along `+x` through a 1D domain. The `-x` face is an **inflow**
boundary (`BC::inflow`) and the `+x` face is an **outflow** boundary
(`BC::outflow`).

A uniform streaming plasma is the exact steady-state solution, so this test
validates that the inflow/outflow pair keeps a uniform state uniform:

* Upstream density `rhoS0` and bulk velocity `uxS0` are preserved at the inflow face.
* The guide field `Bx` stays uniform and unchanged.
* No spurious electric field develops.
* Ion kinetic energy `Epart` and magnetic energy `Eb` stay finite and bounded.

The upstream state is prescribed by the `#INFLOW` command (see `PARAM.XML` for
command syntax and parameter descriptions).

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
2. **Plot output** (`validate_plot`):
   * `<uxS0>` conserved across the domain,
   * guide field `Bx` uniform and unchanged,
   * inflow-side (first third of the domain) `rhoS0` preserved to ~10%
     (the inflow maintains the upstream state rather than draining),
   * no spurious mean `Ex`/`Ey`/`Ez`.
