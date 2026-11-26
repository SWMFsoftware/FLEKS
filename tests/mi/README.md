# 3D Hybrid-PIC Magnetosphere-Ionosphere (MI) Coupling Test

Fast standalone validation test for FLEKS's self-contained `-mi` ionosphere
potential solver, FAC extraction, and surface-E mapping (`./Config.pl -amrex3d -mi`).

- **Grid**: `16 x 16 x 16` cells with `8 x 8 x 8` blocks (8 boxes) on `[-18, 14] x [-14, 14]^2 R_E`.
- **Physics**: Hybrid PIC (`#FIELDSOLVER hybrid`), dipole field (`-2700 nT` at `R_E`), absorbing inner sphere (`R = 3 R_E`), `#BODYBOUNDARY ionosphere`, auroral conductance on a `46 x 91` hemispherical grid, 20 steps.
- **Checks**: Finite PIC energy log and `ionosphere_n0000020.dat` validation (CPCP range and North/South symmetry, non-zero FACs, and auroral Pedersen enhancement above the 3.0 S floor).
