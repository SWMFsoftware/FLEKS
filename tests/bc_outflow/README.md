# Outflow Field-Boundary Tests

Tests of the `outflow` field faces (`#FIELDBOXBOUNDARY outflow`), whose ghost
cells and nodes copy the boundary value (zero normal gradient).

## Variants

| Variant | File | Physics | Validation |
|---------|------|---------|------------|
| `bc_outflow_largedt_fields` | `PARAM.in.largedt.fields` | Vacuum GEM field between outflow y-faces, semi-implicit solver at `c*dt/dx = 7`, `theta = 0.5`, 60 steps | Field energy bounded: last step within 5% of the maximum of the first 30 steps, at most 100x the initial energy |

## Physics Setup

- **Domain**: the GEM box (25.6 x 12.8 d_i) on a coarse 16 x 8 grid, periodic in
  x, outflow field faces in y, no particles (`#PARTICLES 0`).
- **Field**: the Harris-sheet magnetic field with the GEM perturbation; without
  the sheet current it launches electromagnetic waves toward the y-faces.
- **Time step**: `dt = 3.16`, i.e. `c*dt/dx = 7`, with `theta = 0.5` and
  `ratioDivC2C = 0`.

## Validation

The zero-gradient outflow condition extrapolates both E and B, so it is not an
energy-conserving or a strictly absorbing boundary: the field energy rises by
about 13x in the first steps and then stays flat. The test checks that it does
not keep growing. This guards the boundary treatment inside the GMRES operator
of the implicit E solve: the correction must see the zero-gradient ghost nodes,
the same as the full boundary fill applied to E afterwards. With zero ghost
nodes in the operator (the earlier behavior) the energy grew exponentially and
reached 1e22 times its initial value in 60 steps.

## Running

```bash
python3 tests/validate_tests.py --test=bc_outflow
```
