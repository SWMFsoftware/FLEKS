# 3D AMR Alfvén Pulse Test (`tests/pulse_amr_3d`)

Propagates a localized Gaussian Alfvén pulse across a 3D coarse-fine AMR interface
with the hybrid-PIC field solver (`#HYBRIDPIC` with kinetic ions and massless
fluid electrons).

## Physics & Numerical Setup

- **Domain**: 3D periodic box $x \in [-4.0, 4.0]$, $y \in [-2.0, 2.0]$, $z \in [-1.0, 1.0]$.
- **Grids**: Level 0 grid with $32 \times 16 \times 8$ cells; Level 1 ($2\times$ refinement)
  covers the central 3D box $|x| \le 2.0$, $|y| \le 1.0$, $|z| \le 0.5$.
- **Initial condition**: Localized transverse Gaussian pulse (`#TESTCASE AlfvenPulse`)
  $B_y(x) = B_1 \exp(-(x/\sigma)^2)$ on a background guide field $B_{x0} = 10^{-9}$ T.
  Analytically divergence-free everywhere ($\nabla \cdot \mathbf{B} = 0$).
- **AMR div(B) treatment**:
  - `syncEmfAmr = T`: interface nodes share synchronized nodal electric fields.
  - `ctRestrictB = T`: covered coarse cells are advanced with injected EMF and relaxed
    towards fine-cell averages using the true 3D nodal vector potential $\mathbf{A}$
    ($\mathbf{B} \mathrel{+}= \nabla \times \mathbf{A}$).
  - `evolveGhostB = T`: first fine ghost layer is advanced by Faraday's law rather
    than re-interpolated from coarse level.
- **Divergence monitoring**:
  - `#DIVB` with `alwaysComputeDivB = T` outputs `divB-AMR` diagnostics at every
    `#MONITOR dnReport` step.

## Validation

The test validates:
1. **Energy conservation**: Total energy drift is bounded ($< 5\%$) and magnetic
   energy remains stable without whistler/Hall instability.
2. **Divergence preservation**: Coarse-fine interface and interior $\nabla \cdot \mathbf{B}$
   remain at machine round-off ($< 10^{-10}$, typically $\sim 10^{-15}$).
3. **Spatial structure**: PostIDL .out frames detect both coarse and fine grid levels
   and confirm the propagation of the $B_y$ pulse.

## Running

```bash
python3 tests/validate_tests.py --test=pulse_amr_3d -v
```
