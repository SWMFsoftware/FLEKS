# Inflow / Outflow Open-Boundary Hybrid Test

## Physics

A uniform magnetized hybrid plasma (kinetic ions + massless fluid electrons)
streams along `+x` through a 1D domain.  The `-x` face is an **inflow**
boundary (`BC::inflow`) and the `+x` face is an **outflow** boundary
(`BC::outflow`).

A perfectly uniform streaming plasma is the exact steady solution, so this test
validates that the inflow/outflow pair keeps a uniform state uniform:

* the upstream density `rhoS0` and bulk velocity `uxS0` are preserved at the
  inflow face,
* the guide field `Bx` stays uniform and unchanged,
* no spurious electric field develops,
* ion kinetic energy `Epart` and magnetic energy `Eb` stay finite and bounded.

## The `#INFLOW` command

This test uses the `#INFLOW` command to prescribe the upstream state at the
inflow face (in SI units, mirroring `#UNIFORMSTATE`):

```
#INFLOW
1.0439e-8    bx [T]
0.0          by [T]
0.0          bz [T]
5.0          rho [amu/cc]
20.0         ux [km/s]
0.0          uy [km/s]
0.0          uz [km/s]
10000.0      T [K]
```

When an `inflow` face is present and this block is supplied:

1. **Field BC** — `Pic::apply_inflow_wall` (in `src/Pic.cpp`) pins the ghost
   B to the prescribed `(bx, by, bz)` and the ghost tangential E to the upstream
   motional E = -u × B; the normal E is left free (copied from the mirrored
   valid cell).  This is a TRUE inflow field BC: the user controls the boundary
   field rather than the ghost inheriting the interior value.

2. **Particle BC** — `Particles::inject_particles_at_boundary` (in
   `src/Particles.cpp`) samples the re-seeded inflow particles from a drifting
   Maxwellian with the prescribed upstream density, bulk velocity, and
   temperature, routed through the existing `add_particles_cell` Maxwellian
   injector via a per-species `Vel` override (the `InflowVel` struct on
   `FluidInterface`).

Outgoing particles that cross the inflow face outward are removed and tallied
(in `Particles::reflect_or_delete_particle`), matching the `absorb` behaviour.

Without an `#INFLOW` block, an `inflow` face falls back to the zero-gradient
ghost copy (the original open-boundary behaviour): `Pic::use_float` copies the
adjacent physical cell's field, and the particles are sampled from the live
fluid-interface state.

## Running

```bash
# run just this test
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
