# Electron Impact Ionization (MCC) Test

This test verifies the exospheric electron impact ionization module in FLEKS, which is implemented using a **Monte Carlo Collision (MCC)** algorithm.

## What is being tested?
*   **Probability Evaluation (Lotz Cross-Sections)**: Evaluates the energy-dependent collision probability $P_j$ for each electron moving through cells with background neutral exosphere densities $n_n$:
    $$P = 1 - \exp\left(-n_n \sigma(E_k) v_e \Delta t\right)$$
    where $\sigma(E_k)$ is computed dynamically using the Lotz formulation based on electron energy.
*   **Kinematic Energy Partitioning (Opal-Beaty Formulation)**: When a collision is triggered, the incident energy is partitioned into primary and secondary electrons using the Opal-Beaty semi-empirical distribution:
    $$E_{sec} = B \tan\left(R \arctan\left(\frac{E_k - E_{iz}}{B}\right)\right)$$
    $$E_{pri} = E_k - E_{iz} - E_{sec}$$
*   **Neutralizing Pair Spawning**: Generates the new secondary electron and ion at the *exact same spatial coordinate* and with matching macroparticle weights, ensuring strict charge conservation.
*   **Blowout Mitigation**: Triggers conservative moment-preserving particle merging using the `#ELECTRONIMPACT` configured threshold (`blowoutLimitRatio = 2.0`) to avoid numerical particle explosion.

## How is it verified?
The test is automatically verified by `validate_tests.py` using:

1.  **Ionization Yield Verification**: Asserts that the actual physical particle count for all species at $t = 0.2$ is strictly higher than the nominal exosphere photoionization rates:
    $$\text{Yield}_j = N_{actual, j} - N_{nominal\_exosphere, j} > 0$$
    Due to the high neutral density ($n_0 = 1.0\times 10^{15}$ m$^{-3}$) and high thermal energy ($T_0 = 500,000$ K), electrons undergo significant impact ionization, resulting in a large, positive physical yield that validates all MCC and pair-spawning code paths.
