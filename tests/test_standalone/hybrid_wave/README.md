# Hybrid Wave Propagation Standalone Test

This test verifies the physical correctness and numerical stability of the **Hybrid PIC (kinetic ions, fluid electrons)** solver in standalone FLEKS.

## Physics & Configuration

The test initializes a 1D grid of 64 cells with:
- 1 kinetic ion species (charge $q = 1$, mass $m = 1$).
- A neutralizing massless electron fluid.
- A uniform background magnetic field $B_x$.
- Small perturbations to test wave propagation.

### Active Solvers:
- `#HYBRIDPIC` = T (treats electrons as a fluid and advances ions kinetically).
- `#SOLVEEM` = F (disables the standard Maxwell solver).
- `#RESISTIVITY` = 0.01 ($\eta$ resistivity term active).
- `#ELECTRONTEMPERATURE` = 1.5 ($T_{e0} = 1.5$, isothermal $\gamma_e = 1.0$, reference density $n_{e0} = 5.0$).
- `#HALLSUBCYCLE` = 2 (subcycling active for the magnetic field solver).

## Expected Results & Success Criteria
1. The solver must complete the execution successfully without encountering segmentation faults or NaN values.
2. The conservation of total ion particle counts must hold with zero particle loss (due to periodic boundaries).
3. The convective electric field $\mathbf{E} = -\mathbf{U}_i \times \mathbf{B}$ must be correctly evolved.
