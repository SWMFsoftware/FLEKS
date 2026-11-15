# Absorbing Boundary-Condition Tests

This directory groups the absorbing-boundary tests.  Each variant is a
`PARAM.in.<suffix>` file; the runner discovers all of them and runs each once.

## Variants

| Variant | File | Physics | Validation |
|---------|------|---------|------------|
| `bc_absorb_fields` | `PARAM.in.fields` | Tophat EM pulse (Ey/Bz) absorbed at the x-faces (`#FIELDBOXBOUNDARY absorb`, full-PIC) | EM energy decays; late interior field reduced |
| `bc_absorb_hybrid_fields` | `PARAM.in.hybrid.fields` | Localized Gaussian shear-Alfvén pulse (`AlfvenPulse`) absorbed at x-faces (`#FIELDBOXBOUNDARY absorb`, hybrid-PIC) | Net energy decreases (`Etot_final < Etot_initial`); late `max|By|` < 20% of initial peak |
| `bc_absorb_particles` | `PARAM.in.particles` | Ions with bulk +x flow drain out through absorbing x-walls (`#PARTICLEBOXBOUNDARY absorb`, pure particle push) | Epart decays toward zero |
| `bc_absorb_hybrid_particles` | `PARAM.in.hybrid.particles` | Same ion drain through absorbing x-walls, but in the hybrid mover (`useHybridPIC = T`) | Epart decays toward zero |

## Validation

All four variants verify the **absorbing** boundary condition for fields and
particles. `validate.py` branches on the variant name:
1. **Energy log**:
   - `bc_absorb_fields`: a `tophat` EM pulse seeded in the interior leaves through
     the absorbing (outgoing-characteristic) x-faces; the pulse carries EM energy
     out, so `Etot` (Ee+Eb) decays below 25% of its initial value — the decisive
     absorber signature (conducting/periodic walls would reflect the pulses and
     keep the energy).
   - `bc_absorb_hybrid_fields`: verifies that a localized `AlfvenPulse` exits
     and is absorbed at the `absorb` x-faces without reflection or instability
     (`Etot` decreases cleanly as wave energy leaves).
   - `bc_absorb_particles` / `bc_absorb_hybrid_particles`: a single ion species
     with a bulk `+x` flow streams ballistically toward the absorbing x-walls;
     particles are removed on exit, so `Epart` decays below 1% of its initial
     value (full evacuation).
2. **Plot**:
   - `bc_absorb_fields`: late `|Bz|,|Ey|` are much reduced from the seeded 1.0
     (max fields < 0.7; no reflected standing wave).
   - `bc_absorb_hybrid_fields`: fields (`Bx`, `By`, `Bz`, `Ex`, `Ey`, `Ez`) remain
     finite across all frames and late `max|By|` drops below 20% of its initial
     seeded peak (confirming >96% perturbation magnetic energy absorption).
   - `bc_absorb_particles` / `bc_absorb_hybrid_particles`: late `rhoS0` peak is
     evacuated to zero (peak < 0.1), verifying that particles were absorbed rather
     than trapped or reflected.

## Running

```bash
python3 tests/validate_tests.py --test=bc_absorb                   # all 4 variants
python3 tests/validate_tests.py --test=bc_absorb.fields            # full-PIC fields only
python3 tests/validate_tests.py --test=bc_absorb.hybrid.fields     # hybrid-PIC fields only
python3 tests/validate_tests.py --test=bc_absorb.particles         # full-PIC particles only
python3 tests/validate_tests.py --test=bc_absorb.hybrid.particles  # hybrid+particles only
```
