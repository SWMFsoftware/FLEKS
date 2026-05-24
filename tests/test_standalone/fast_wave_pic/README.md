# 2D Alfven Wave Standalone Test - Full PIC

This test case implements a standalone 2D Alfven Wave propagation simulation using the **Full PIC** solver in FLEKS.

## Configuration
- **Active Solver**: Full PIC (`#SOLVEEM = T`, `#HYBRIDPIC = F`).
- **Test Case**: `AlfvenWave` (which sets up sinusoidal density, velocity, and magnetic field perturbations corresponding to SWMF test17).
- **Grid**: $32 \times 24 \times 1$ cells on a periodic 2D domain.
- **Species**: 2 kinetic species (kinetic ions and kinetic electrons).
