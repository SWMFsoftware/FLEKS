# Combined Exosphere Processes Test

This test verifies all three exospheric processes running concurrently in a single simulation: exospheric photoionization (automatic pair-injection), electron impact ionization MCC, and exospheric charge exchange MCC.

## What is being tested?
1.  **Exospheric Photoionization**: Converts exospheric neutral density profiles into neutralizing electron-ion pairs, validating strict local charge conservation.
2.  **Electron Impact Ionization MCC**: High-temperature electrons undergo collisions with exospheric neutrals (Lotz cross-section and Opal-Beaty kinematically partitioned energy), spawning secondary electrons/ions and resulting in physical plasma yields.
3.  **Exospheric Charge Exchange MCC**: Hot drifting plasma ions ($u_x = 400.0$ km/s) undergo collisions with the cold neutral background, replacing their velocities in-place with cold exosphere thermal velocities (Lindsay-Stebbings cross-section and thermal Box-Muller Maxwellian sampling), causing significant deceleration/cooling of the plasma species.

## How is it verified?
The test is automatically verified by `validate_tests.py` using:

1.  **Multi-Process Yield Verification**: Asserts that the actual physical particle count for all species (H+, O+, and e-) at the end of the simulation exceeds the nominal photoionization rate due to positive electron impact ionization yields.
2.  **Multi-Process Deceleration/Cooling Verification**: Asserts that the mean velocity $v_x$ of both H+ and O+ is significantly reduced (cooled) below the initial uniform drift velocity ($400.0$ km/s) due to charge exchange events.
