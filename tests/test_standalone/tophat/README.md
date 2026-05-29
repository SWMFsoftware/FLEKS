# 1D Electromagnetic Top-Hat Wave Propagation Test

## Description
This is a standalone test designed to verify the physical propagation of a 1D electromagnetic square (top-hat) pulse in a periodic domain, benchmarked against low-dispersion/low-diffusion properties of the field solver in FLEKS.

## Physics & Solver Setup
- **Geometry & Boundaries**: 1D periodic grid ($x \in [-1.0, 1.0]$, 64 cells in $X$, 1 in $Y$ and $Z$).
- **Physical Model**: Pure electromagnetic field propagation. No kinetic particles are initialized (`nPartPerCell = IntVect::Zero`), meaning the plasma density is effectively zero and there are no particle effects.
- **Initial Fields**: A square top-hat pulse is initialized in the central 50% region ($[-0.5, 0.5]$):
  - $E_y = 1.0$
  - $B_z = 1.0$
  - Outside this interval, all fields are set to $0.0$.
- **Electromagnetic Fields**: Field solver is enabled (`solveEM = T`) with standard explicit/implicit discretization parameters.
- **Propagation Dynamics**: In normalized units with speed of light $c = 1$, the initial state represents a pure right-going transverse wave pulse. Over time, the top-hat pulse propagates to the right at speed $v = 1.0$ while keeping its shape.

## Expected Results
- **Amplitude Retention**: The peak magnetic field $B_z$ must remain close to $1.0$ (subject to minor numerical dispersion/diffusion of the square shape on a discrete grid, typically in the range $[0.8, 1.2]$).
- **Decoupled Fields**: The transverse magnetic field component $B_y$ must remain negligible (near $0.0$) throughout the run, confirming there is no unphysical cross-coupling of wave modes.
