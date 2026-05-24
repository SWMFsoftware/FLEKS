# Ion Beam-Beam Instability Test

## Description
This is a standalone test designed to verify the simulation of the ion beam-beam instability in FLEKS. The test initializes a uniform background plasma moving in the $-x$ direction and a beam of ions of a smaller density moving in the $+x$ direction, solving for electromagnetic field dynamics.

## Physics & Solver Setup
- **Geometry & Boundaries**: 1D periodic grid (64 cells in $X$, 1 in $Y$ and $Z$).
- **Physical Model**: A single charged species ($m = 1.0$, $q = 1.0$) with uniform background density $5.0\,\text{amu/cc}$ and background velocity $-400\,\text{km/s}$ ($u_x = -0.4$ in normalized units).
- **Beam setup**: $1\%$ of the particles have a beam velocity of $+400\,\text{km/s}$ ($u_x = +0.4$ in normalized units).
- **Electromagnetic Fields**: Field solver is enabled (`solveEM = T`) with a guiding magnetic field $B_x = 5\,\text{nT}$ ($5.0 \times 10^{-9}$ T).
- **Purpose**: Verifies that the particle count and mean velocity are conserved, and tracks the growth of cyclotron waves in the transverse magnetic fields $B_y$ and $B_z$ from seed noise.

## Expected Results
- **Particle Count Conservation**: The total number of particles remains exactly constant.
- **Mean Velocity**: The expected species mean velocity `MeanVx` is:
  $$\text{MeanVx} = 0.01 \times (0.4) + 0.99 \times (-0.4) = -0.392$$
  which must be conserved within $1\%$ tolerance.
- **Wave Growth**: Transverse components $B_y$ and $B_z$ grow over time due to the instability seeding from thermal noise, indicating the successful propagation and growth of cyclotron waves.
