# Ion Beam-Beam Instability Test

## Description
This is a standalone test designed to verify the simulation of the ion beam-beam instability in FLEKS. The test initializes a uniform background plasma moving in the $-x$ direction and a beam of ions of a smaller density moving in the $+x$ direction, solving for electromagnetic field dynamics. The test directory holds two parameter files that exercise the same physics with the two field solvers:

- **`PARAM.in`** — full-PIC: kinetic electrons and ions (Maxwell/GMRES `solveEM` solver).
- **`PARAM.in.hybrid`** — hybrid-PIC: kinetic ions / massless fluid electrons (Ohm's-law + Faraday `useHybridPIC` solver).

Both variants share the identical grid, guide field, beam ratio/speeds, density, normalization and output, so the same FFT-based validation applies to both.

## Physics & Solver Setup
- **Geometry & Boundaries**: 1D periodic grid (64 cells in $X$, 1 in $Y$ and $Z$).
- **Physical Model**: A single charged species ($m = 1.0$, $q = 1.0$) with uniform background density $5.0\,\text{amu/cc}$ and background velocity $-400\,\text{km/s}$ ($u_x = -0.4$ in normalized units).
- **Beam setup**: $1\%$ of the particles have a beam velocity of $+400\,\text{km/s}$ ($u_x = +0.4$ in normalized units).
- **Electromagnetic Fields**:
  - Full PIC: field solver enabled (`solveEM = T`) with a guiding magnetic field $B_x = 5\,\text{nT}$ ($5.0 \times 10^{-9}$ T).
  - Hybrid PIC: `useHybridPIC = T`, `solveEM = F`, Hall term ON (single sub-step), resistivity and electron-pressure-gradient OFF so the instability grows undamped; the guide field is dynamically weak ($v_A \ll u_{\mathrm{beam}}$) so the explicit hybrid advance is not stiff at this resolution.
- **Purpose**: Verifies that the particle count and mean velocity are conserved, and tracks the growth of cyclotron waves in the transverse magnetic fields $B_y$ and $B_z$ from seed noise.

## Expected Results
- **Particle Count Conservation**: The total number of particles remains exactly constant.
- **Mean Velocity**: The expected species mean velocity `MeanVx` is:
  $$\text{MeanVx} = 0.01 \times (0.4) + 0.99 \times (-0.4) = -0.392$$
  which must be conserved within $1\%$ tolerance.
- **Wave Growth**: Transverse components $B_y$ and $B_z$ grow over time due to the instability seeding from thermal noise, indicating the successful propagation and growth of cyclotron waves.
- **Transverse-Wave Resonance Check**: At the final output frame ($t \approx 0.1$), the validation performs a DFT of the spatial $B_y$/$B_z$ profile and compares the dominant mode to the theoretical cyclotron-resonant wavenumber $k_{\mathrm{res}} = \Omega_i / \Delta v$, where $\Omega_i = q_p B_x / m_p$ and $\Delta v = 800\,\text{km/s}$. For this test the resonant wavelength ($\sim 10^4\,\text{km}$) exceeds the $2\,\text{km}$ periodic box, so the resonant mode cannot fit; the check instead verifies that the wave has grown above the noise floor and that its power is concentrated in low-order spatial modes consistent with the box-limited instability.
- **Hybrid stability (energy check)**: the hybrid variant additionally runs the shared hybrid-family energy checks — magnetic and ion energies stay finite and the magnetic energy does not blow up (catches a NaN / missing-Hall-factor runaway).

## Running

```bash
python3 tests/validate_tests.py --test=beam
```

The output table lists both `BEAM` (full PIC) and `BEAM (HYBRID)` (hybrid PIC) results.
