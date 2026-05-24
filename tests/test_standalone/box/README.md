# 0D Box Equilibrium Test

## Description
This is a zero-dimensional (0D) equilibrium test designed to verify the chemical and kinetic balance of the exosphere model in FLEKS.

## Physics & Solver Setup
- **Geometry & Boundaries**: Fully periodic 3D domain (from $-1.0$ to $1.0$ in all directions).
- **Physical Model**: A single charged species ($m = 1.0$, $q = 1.0$) with exponential exosphere density profile ($H_0 = 0.03$, $T_0 = 200\,\text{K}$).
- **Electromagnetic Fields**: Field solver is disabled (`solveEM = F`).
- **Purpose**: Verifies that the macroparticle injection matches the total physical production rate over time under stationary equilibrium.

## Expected Results
- **Particle Count Verification**: The number of physical particles $N_{phys}$ should follow:
  $$N_{phys}(t) = \text{TotalProductionRate} \times (t + dt)$$
- The validation script checks that the relative difference between actual and expected particle counts is within 2%.
