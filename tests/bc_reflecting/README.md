# Reflecting & Conducting Wall Boundary-Condition Tests

This directory groups tests for reflecting particle boundaries (`#PARTICLEBOXBOUNDARY reflect`)
and perfectly conducting field walls (PEC, `#FIELDBOXBOUNDARY conducting`).

## Variants

| Variant | File | Physics | Validation |
|---------|------|---------|------------|
| `bc_reflecting_fields` | `PARAM.in.fields` | Tophat EM pulse (Ey/Bz) reflecting off conducting (PEC) x-faces (`#SOLVEEM T`) | EM energy conserved ($0.8 \le E_1/E_0 \le 1.25$); fields bounded |
| `bc_reflecting_largedt_fields` | `PARAM.in.largedt.fields` | The same pulse in vacuum with the semi-implicit solver at `c*dt/dx = 3.5` and `theta = 0.5` | EM energy conserved at every step ($\max|E/E_0 - 1| \le 10^{-4}$) |
| `bc_reflecting_hybrid_fields` | `PARAM.in.hybrid.fields` | Transverse shear-Alfvén wave (`HybridWave`) in a PEC cavity (`#HYBRIDPIC T`) | Run finite and stable; fields bounded near walls |
| `bc_reflecting_particles` | `PARAM.in.particles` | Kinetic ions streaming toward +x wall specularly reflect (pure particle push, `solveEM = F`) | Particle kinetic energy conserved ($E_{part} > 50\%$); particles confined |
| `bc_reflecting_hybrid_particles` | `PARAM.in.hybrid.particles` | Kinetic ions streaming toward +x wall specularly reflect in hybrid mover (`useHybridPIC = T`) | Particle kinetic energy conserved ($E_{part} > 50\%$); particles confined |

## Physics Setup

### Particle Reflection (`PARAM.in.particles` / `PARAM.in.hybrid.particles`)
- **Domain**: a closed box in x (`8×1×1` cells, `Lx = 6.4 d_i`), fake-2D periodic in y/z. Both x-faces are **reflecting**.
- **Plasma**: a single kinetic ion species (`q = m = 1`) streaming toward the `+x` wall (`ux > 0`) with an oblique `uy != 0`.
- **Specular Reflection**: $v_n \to -v_n$, position mirrored $x \to 2\cdot x_w - x$.
- **Validation**:
  1. **Energy log**: $E_{part}$ does not drain toward zero (stays $> 50\%$ of initial).
  2. **Plot output**: ion density $\rho_{S0}$ stays positive and confined across domain and boundary cells; normal velocity $ux_{S0}$ is finite.

### Conducting Field Walls (`PARAM.in.fields` / `PARAM.in.hybrid.fields`)
- **`PARAM.in.fields` (Full-PIC)**:
  - An interior tophat EM pulse (Ey/Bz step) propagates outward and reflects back and forth between PEC boundaries ($E_t = 0, B_n = 0$).
  - **Validation**: unlike absorbing boundaries (where energy decays $< 25\%$), conducting walls reflect the pulses and strictly conserve total electromagnetic energy ($E_{tot}$ conserved within discretization bounds, $0.8 \le E_1/E_0 \le 1.25$).
- **`PARAM.in.largedt.fields` (Full-PIC, large time step)**:
  - The pulse of `PARAM.in.fields` in vacuum (`#PARTICLES 0`), with `dt = 1.75` (`c*dt/dx = 3.5`), `theta = 0.5` and `ratioDivC2C = 0`, for which the semi-implicit scheme conserves the EM energy exactly; 80 steps cover 7 wall-to-wall transits.
  - **Validation**: the energy stays within $10^{-4}$ of its initial value at every step (about $3\times10^{-6}$ is the GMRES-tolerance drift). This guards the conducting wall inside the GMRES operator: the correction must see $E_t = 0$ on the wall and mirrored ghost nodes, with identity rows for the wall unknowns. Without them the energy grew at each reflection ($5\times10^{-3}$ after 10 steps, $10^{19}$ after 80).
- **`PARAM.in.hybrid.fields` (Hybrid-PIC)**:
  - A transverse shear-Alfvén wave (`HybridWave`) in a finite-x domain with conducting field walls and reflecting particle walls.
  - **Validation**: verifies that the hybrid field advance (`assemble_ohm_E` and `apply_centerB_BC`) stays finite, positive, and bounded without numerical instability.

## Running

All five variants run together under the test name:

```bash
python3 tests/validate_tests.py --test=bc_reflecting                      # all 5 variants
```

To run a **single variant**, append its token to the test name:

```bash
python3 tests/validate_tests.py --test=bc_reflecting.fields               # full-PIC conducting field wall
python3 tests/validate_tests.py --test=bc_reflecting.largedt.fields       # same, semi-implicit at c*dt/dx = 3.5
python3 tests/validate_tests.py --test=bc_reflecting.hybrid.fields        # hybrid-PIC conducting field wall
python3 tests/validate_tests.py --test=bc_reflecting.particles            # full-PIC particle reflection
python3 tests/validate_tests.py --test=bc_reflecting.hybrid.particles     # hybrid-PIC particle reflection
```

Legacy tokens `--test=bc_reflecting.full` and `--test=bc_reflecting.hybrid` are mapped to `particles` and `hybrid.particles` respectively for backward compatibility.
