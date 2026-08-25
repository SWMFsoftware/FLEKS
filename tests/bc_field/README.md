# bc_field

Grouped test for the hybrid-PIC **field** boundary conditions.
Unlike the full-PIC `bc_absorb_fields` (a pure-EM, zero-particle `tophat` pulse),
these variants run **hybrid-PIC** (`useHybridPIC=T`) so the generalized Ohm's-law
field advance is exercised — which requires kinetic ions to be present, so the
`tophat` testcase cannot be reused here.

## Variants

Both variants use the identical geometry, seed and solver settings; they differ
only in the `#BFIELDBOXBOUNDARY` at the two x faces:

| Variant | Field BC (x-lo/x-hi) | Particle BC (x) |
|---------|----------------------|-----------------|
| `bc_field (HYBRID.CONDUCTING)` | `conducting` (PEC: `E_t=0`, `B_n=0`) | `reflect` |
| `bc_field (HYBRID.ABSORB)`     | `absorb` (outgoing characteristic)    | `reflect` |

- A transverse shear-Alfvén wave is seeded with the `HybridWave` testcase (guide
  field `Bx`, `δBy`/`δBz` perturbation + velocity kick) in a finite-x cavity.
- The x **field** faces are the object under test; the x **particle** faces are
  `reflect` so the ions (and hence the Ohm's-law moments) stay inside.
- y and z are thin (fake-2D/1D) and periodic.

## Purpose

Add hybrid-PIC regression coverage for the `conducting` and `absorb` **field**
boundary conditions, which currently have no hybrid test:
- `conducting` exercises `Pic::apply_conducting_wall` on the cell-centred
  hybrid `E`/`B` (`centerEhybrid`, `centerB`) via `assemble_ohm_E` and
  `apply_centerB_BC`.
- `absorb` exercises `Pic::apply_absorbing_wall` on the same cell-centred
  fields.

The validator confirms the runs are **finite and stable** (no NaN, bounded
energy), catching any field-wall regression that would destabilise the hybrid
advance.  A clean physics discriminator (e.g. standing-wave energy decay for
`absorb`) is not asserted because a domain-filling `HybridWave` standing wave is
a weak absorber test; the primary value here is stability/regression coverage
plus the PEC/absorbing-wall wiring check.

## Validation

`validate.py`:
1. **Energy log**: `Etot` and `Epart` finite; `Etot` positive and bounded
   (`Etot` does not blow up by `>1e3x` from its initial value).
2. **Plot**: `By/Bz/Ey/Ez` finite and bounded (`peak < 1e6`) in the last two
   frames — no NaN / no blow-up near the walls.

## Running

All variants run together under the test name:

```bash
python3 tests/validate_tests.py --test=bc_field          # both variants
```

To run a **single variant**, append its token to the test name:

```bash
python3 tests/validate_tests.py --test=bc_field.hybrid.absorb       # hybrid absorb only
python3 tests/validate_tests.py --test=bc_field.hybrid.conducting   # hybrid conducting only
```
