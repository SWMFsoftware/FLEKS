# 3D AMR Traveling Wave Test (`tests/wave_amr_3d`)

Propagates a unidirectional transverse circularly polarized wave across a 3D coarse-fine AMR interface
under periodic boundary conditions. Two variants are provided:

- **Full PIC** (`PARAM.in`): Kinetic ions and kinetic electrons ($m_i/m_e = 25$) with the electromagnetic field solver (`#SOLVEEM T`).
- **Hybrid PIC** (`PARAM.in.hybrid`): Kinetic ions and massless fluid electrons with the generalized Ohm's law solver (`#HYBRIDPIC T`).

## Physics & Numerical Setup

- **Domain**: 3D periodic box $x \in [-4.0, 4.0]$, $y \in [-2.0, 2.0]$, $z \in [-1.0, 1.0]$.
- **Grids**: Level 0 grid with $32 \times 16 \times 8$ cells; Level 1 ($2\times$ refinement)
  covers the central 3D box $|x| \le 2.0$, $|y| \le 1.0$, $|z| \le 0.5$.
- **Initial condition**: Transverse circularly polarized eigenmode (`#TESTCASE HybridWave`)
  with $k_x = 2\pi / L_x$ ($\lambda = 8.0$), $B_y = B_1 \cos(k_x x)$, $B_z = -B_1 \sin(k_x x)$ on a guide field $B_{x0} = 10^{-9}$ T.
  Matching Walén velocity kick and motional $\mathbf{E} = -\mathbf{u} \times \mathbf{B}$ produce a unidirectional wave propagating along $+x$.
- **AMR div(B) treatment (Hybrid)**:
  - Synchronized nodal electric fields across coarse-fine interface nodes.
  - Covered coarse cells advanced with the coarse EMF and relaxed towards fine-cell
    averages using the true 3D nodal vector potential $\mathbf{A}$
    ($\mathbf{B} \mathrel{+}= \nabla \times \mathbf{A}$).
  - First fine ghost layer advanced by Faraday's law rather than re-interpolated
    from coarse level.
- **Divergence monitoring**:
  - `#DIVB` with `alwaysComputeDivB = T` outputs `divB-AMR` diagnostics at every
    `#MONITOR dnReport` step in the hybrid run.

## Validation

The test validates:
1. **Energy conservation**: Total energy drift is bounded ($< 5\%$) and magnetic
   energy remains stable without whistler/Hall instability.
2. **Divergence preservation**: In the hybrid run, coarse-fine interface and interior $\nabla \cdot \mathbf{B}$
   remain at machine round-off ($< 10^{-10}$, typically $\sim 10^{-15}$).
3. **Spatial structure**: PostIDL .out frames detect both coarse and fine grid levels
   and confirm the propagation of the wave.

## Running

```bash
# Run both full-PIC and hybrid variants
python3 tests/validate_tests.py --test=wave_amr_3d -v

# Or run a specific variant:
python3 tests/validate_tests.py --test=wave_amr_3d -v --variant=hybrid
```
