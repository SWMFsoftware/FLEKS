# Photoionization (Automatic Pair-Injection) Test

This test verifies the exospheric photoionization process in FLEKS, specifically validating the **automatic pair-injection** of ions and neutralizing electrons to ensure strict local charge conservation.

## What is being tested?
*   **Automatic Pair-Injection**: For every ion macroparticle spawned from the exospheric neutral density profile, a neutralizing electron macroparticle is simultaneously generated at the *exact same spatial coordinate* with the *exact same macroparticle weight*. 
*   **Charge Conservation**: Strict local and global charge neutrality is preserved throughout the injection process, preventing grid charging and artificial electrostatic noise.
*   **Leapfrog Mover Physics**: Newly spawned pickup ions (H+ and O+) are accelerated by perpendicular electric ($E_y=1.0$) and magnetic ($B_z=1.0$) fields, executing classical cycloidal $E \times B$ drift motion.

## How is it verified?
The test is automatically executed and validated by the test runner `validate_tests.py` using:

1.  **Production Rate Verification**: Asserts that the physical particle population $N(t)$ for each species (0 = H+, 1 = O+, 2 = e-) at each diagnostic output time step matches the configured exospheric production rates ($8.0 \times 10^{23}$ for H+, $2.0 \times 10^{23}$ for O+, and $1.0 \times 10^{24}$ for electrons) under a 2% relative tolerance:
    $$N_{expected}(t) = S_{production} \cdot (t + \Delta t)$$
2.  **Velocity History Validation**: Compares the mean drift velocity of H+ ($\omega_c = 1.0$) and O+ ($\omega_c = 0.0625$) against the exact discrete analytical leapfrog mover model:
    $$v_{y, expected} = \frac{1}{N_{steps}}\sum_{k=0}^{N_{steps}} \sin(\omega_c k \Delta t)$$
    $$v_{x, expected} = \frac{1}{N_{steps}}\sum_{k=0}^{N_{steps}} (1 - \cos(\omega_c k \Delta t))$$
    Asserts that the computed mean velocities match the simulation results within a strict tolerance of $0.02$.
