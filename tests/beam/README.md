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
- **Transverse-Wave Resonance Check**: At the final output frame ($t \approx 0.1$), the validation performs a DFT of the spatial $B_y$/$B_z$ profile and compares the dominant mode to the theoretical cyclotron-resonant wavenumber $k_{\mathrm{res}} = \Omega_i / \Delta v$, where $\Omega_i = q_p B_x / m_p$ and $\Delta v = 800\,\text{km/s}$. For this test the resonant wavelength ($\sim 10^4\,\text{km}$) exceeds the $2\,\text{km}$ periodic box, so the resonant mode cannot fit; the check instead verifies that the wave has grown above the noise floor and that its power is concentrated in low-order spatial modes consistent with the box-limited instability.

## Running

From the FLEKS root directory (requires compiled `bin/FLEKS.exe`):

```bash
python3 tests/validate_tests.py --test=beam
```
