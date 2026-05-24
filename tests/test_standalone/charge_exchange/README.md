# Exospheric Charge Exchange (MCC) Test

This test verifies the exospheric charge-exchange process in FLEKS, which is implemented using an efficient in-place **Monte Carlo Collision (MCC)** algorithm.

## What is being tested?
*   **Probability Evaluation (Lindsay-Stebbings Formula)**: Computes the energy-dependent charge exchange cross-section $\sigma_{cx}$ and collision probability $P_j$ for each hot drifting ion moving through the cold neutral exosphere:
    $$\sigma_{cx} = (4.15 - 0.531 \ln(E_{rel}))^2 \times 10^{-20} \text{ m}^2$$
    $$P_j = n_{n,j} \sigma_{cx} v_{rel\_SI} \Delta t_{SI}$$
*   **In-Place Cold Velocity Sampling**: When a collision is triggered, the hot ion is converted into a cold exospheric pickup ion by **updating its velocity vector in-place** with a cold thermal velocity sampled from the exosphere temperature $T_0$ (using Box-Muller transform). This avoids spawning or deleting particles, ensuring perfect charge and particle number conservation.
*   **Blowout Mitigation**: Triggers conservative moment-preserving merging based on the `#CHARGEEXCHANGE` configured threshold (`cxBlowoutLimitRatio = 2.0`).

## How is it verified?
The test is automatically verified by `validate_tests.py` using:

1.  **Ion Beam Cooling / Deceleration Verification**: The simulation initializes a high-velocity drifting plasma beam ($u_x = 400.0$ km/s) traversing a cold exospheric neutral background ($T_0 = 200.0$ K). The validator asserts that both H+ (Species 0) and O+ (Species 1) experience significant deceleration and cooling (their mean $v_x$ is reduced below the initial $400.0$ km/s):
    $$v_{x, final} < 400.0 \text{ km/s}$$
    This cooling verifies that the LS cross-sections, Monte Carlo selection, and in-place cold Maxwellian velocity replacements are functioning correctly.
