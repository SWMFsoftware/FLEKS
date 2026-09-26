# Inner Body Test (`#BODY`: Absorbing Sphere)

## Description

This standalone test verifies the `#BODY` command, which adds a **real inner
boundary** to a Cartesian PIC domain: a uniform plasma streams in +x past an
absorbing sphere at the origin. Particles that are pushed into a body cell are
removed and tallied, no particle is created inside the body, the electric field
is pinned to zero there, and the particle moments are written as zero, so the
body appears as an empty hole with a wake behind it.

Note that `#PLANETRADIUS` (formerly `#BODYSIZE`) does **not** create such a
boundary — it is only the exosphere reference radius and the `planet` output
unit. `#BODY` is the command that defines the PIC inner boundary.

## Physics & Solver Setup

- **Geometry & Boundaries**: 2D periodic grid (64 × 64 × 1), domain spans
  [-3.2, 3.2] code units in x and y, so `dx = 0.1` and the sphere radius is
  12 cells. One cell in z (fake 2D).

- **Inner Body**: declared with `#BODY` (radius 1.2, center at the origin,
  **code units** like `#REGION`). The mask is built from the cell centers
  (staircase boundary, resolution `dx/2`), and:
  | Quantity | Treatment |
  |----------|-----------|
  | Particles | removed (and tallied) when pushed into a body cell |
  | Injection | initial fill, source and boundary injection skip body cells |
  | E | body nodes excluded from the implicit solve (Dirichlet) ⇒ `E ≡ 0` |
  | B | not forced to zero; frozen at its initial value in the interior |
  | jHat / nodeMM | not masked — their body rows are dropped from the solve, so no charge is silently discarded |
  | div(E) | residual zeroed in body cells |
  | Output | `body` mask variable; particle moments reported as zero inside |

  The absorption test is the **cell-based** boundary, not the exact radius:
  the deposition is node-centred CIC, so an exact-radius test would let
  surface particles deposit into masked nodes every step (a charge sink with
  no bookkeeping) while the field boundary is still the staircase.

- **Plasma Species**:
  | Species | Mass [amu] | Charge [e] | n [amu/cc] | ux [km/s] | T [K] | Role |
  |---------|------------|------------|------------|-----------|-------|------|
  | 0 | 1.0 | +1 | 5.0 | 100 | 314000 | Ions |
  | 1 | 0.04 | −1 | 0.2 | 100 | 314000 | Electrons (n_e = n_i) |

- **Electromagnetic Fields**: enabled (`solveEM = T`) with `B_z = 3.0e-9` T
  and the matching motional field `E = -u × B` = (0, 3.0e-4, 0) V/m, so the
  body has to pin a non-zero ambient E to zero. Solver: implicit GMRES
  (`theta = 0.5`, tol 1e-8, 30 iterations), comoving frame
  (`solveFieldInCoMov = T`, 5 smoothing passes), Lax-Friedrichs upwind
  viscosity on both the E and the B equation (limiter theta = 1, which also
  enables hyperbolic div-B cleaning) and digital-filter smoothing of J
  (1 pass, coefficient 0.5).

- **Time Stepping**: `dt = 0.02` fixed, `TimeMax = 2.0` (100 steps),
  4 × 4 × 1 particles per cell (≈ 131k macroparticles).

## Validation

From the pic-log history:

- `nBodyAbsorb` (cumulative number of absorbed macroparticles, appended after
  the fixed columns) is present, non-decreasing and non-zero at the end;
- all energies stay finite and `Etot` does not grow by more than 10×.

From the last plot frame:

- on every point with `body == 1`, `rhoS0`, `rhoS1`, `Ex`, `Ey`, `Ez` are
  **exactly zero**;
- the region just downstream of the body is depleted to below 50% of the
  upstream density (a wake forms).

Reference numbers of the 4-rank run: 437 of 4096 points inside the body
(geometric expectation πR²/(Lx·Ly) = 452), wake/upstream density = 0.22,
`nBodyAbsorb` → 5.85e4 by t = 2.0 (the count is dominated by the electron
thermal flux, which is much faster than the bulk flow).

## Running

From the FLEKS root directory (requires compiled `bin/FLEKS.exe`):

```bash
python3 tests/validate_tests.py --test=body
```

> Only the full-PIC solver is supported: `#BODY` with `#HYBRIDPIC` aborts, and
> AMR (with refinement regions) is not verified yet.
