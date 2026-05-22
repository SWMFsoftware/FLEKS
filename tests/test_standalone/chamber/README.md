# 1D Chamber Asymmetric Inflow/Outflow Test

## Description
This is a fake 1D chamber test designed to simulate asymmetric inflow and outflow boundary conditions.

## Physics & Solver Setup
- **Geometry & Boundaries**: Non-periodic in the X direction (from $-1.0$ to $1.0$), periodic in Y and Z directions.
- **Physical Model**: A single charged species ($m = 1.0$, $q = 1.0$) with exosphere density source.
- **Purpose**: Verifies that the system achieves a stable steady-state balancing the particle injection and particles escaping through the non-periodic X boundaries.

## Expected Results
- **Steady State Check**: Over the last 1.5 seconds of the 5.0-second run (from $t=3.5$ to $t=5.0$), the particle number should stabilize.
- The validation script checks that the relative change in the total physical particle number during this steady-state phase is less than 5%.
