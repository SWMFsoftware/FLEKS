# 2D Alfven Wave Standalone Test - Hybrid PIC

This test case implements a standalone 2D Alfven Wave propagation simulation using the **Hybrid PIC** solver in FLEKS.

## Configuration
- **Active Solver**: Hybrid PIC (`#SOLVEEM = F`, `#HYBRIDPIC = T`).
- **Test Case**: `AlfvenWave` (which sets up sinusoidal density, velocity, and magnetic field perturbations corresponding to SWMF test17).
- **Grid**: $32 \times 24 \times 1$ cells on a periodic 2D domain.
- **Species**: 1 kinetic species (kinetic ions only). Electrons are treated as a fluid.
