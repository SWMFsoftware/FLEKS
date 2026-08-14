# Absorbing Boundary-Condition Tests

This directory groups the absorbing-boundary tests.  Each variant is a
`PARAM.in.<suffix>` file; the runner discovers all of them and runs each once.
A single `validate.py` branches on the variant name.

## Variants

| Variant | File | Physics | Validation |
|---------|------|---------|------------|
| `bc_absorb_fields` | `PARAM.in.fields` | Tophat EM pulse (Ey/Bz) absorbed at the x-faces (`#BFIELDBOXBOUNDARY absorb`, full-PIC) | EM energy decays; late interior field reduced |
| `bc_absorb_particles` | `PARAM.in.particles` | Ions with bulk +x flow drain out through absorbing x-walls (`#PARTICLEBOXBOUNDARY absorb`, pure particle push) | Epart decays toward zero |

## Purpose

Verify the **absorbing** boundary condition for both fields and particles.

- **`bc_absorb_fields`**: a `tophat` EM pulse seeded in the interior splits
  into left/right-going pulses that leave through the absorbing (outgoing-
  characteristic) x-faces.  The pulses carry their EM energy out, so the total
  EM energy **decays** — the decisive absorber signature (conducting/periodic
  walls would reflect the pulses and keep the energy).
- **`bc_absorb_particles`**: a single ion species with a bulk `+x` flow and
  zero magnetic field streams ballistically toward the absorbing x-walls; both
  x faces are `absorb` (particles removed + tallied), so the ions drain out and
  **Epart decays** — the opposite of the reflecting-wall test.

## Validation

`validate.py` branches on the variant name:
1. **Energy log**:
   - `bc_absorb_fields`: `Etot` (Ee+Eb) decays below 25% of its initial value.
   - `bc_absorb_particles`: `Epart` decays below 50% of its initial value.
2. **Plot**:
   - `bc_absorb_fields`: late `|Bz|,|Ey|` are much reduced from the seeded 1.0
     (no reflected standing wave).
   - `bc_absorb_particles`: `rhoS0` is finite/positive (skipped if not
     deposited in no-solver mode).

## Running

```bash
# run all variants in this directory
python3 tests/validate_tests.py --test=bc_absorb
```
