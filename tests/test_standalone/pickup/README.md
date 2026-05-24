# Uniform Injection and Acceleration (Pickup) Test with Charge Neutrality

## Description
This is a multi-species particle pickup test designed to verify the physical correctness of the charged particle mover under uniform electric and magnetic fields, while maintaining physical charge neutrality.

## Physics & Solver Setup
- **Geometry & Boundaries**: Periodic 3D domain.
- **Physical Model**: Three species representing:
  1. **H+** ($m = 1.0$, $q = 1.0$)
  2. **O+** ($m = 16.0$, $q = 1.0$)
  3. **e-** ($m = 0.01$, $q = -1.0$)
- **Electromagnetic Fields**: Uniform static background fields: $E_y = 1.0$, $B_z = 1.0$.
- **Charge Neutrality**: Maintained dynamically by matching the total production rate of positive ions to the negative electron production rate:
  $$R_{e^-} = R_{H^+} + R_{O^+}$$
  With $R_{H^+} = 8.0 \times 10^{23}$, $R_{O^+} = 2.0 \times 10^{23}$, and $R_{e^-} = 1.0 \times 10^{24}$, the net injected charge is exactly zero.

## Analytical Validation
For continuously injected particles under uniform electric and magnetic fields, the discrete time-integration of the Boris-like mover yields precise velocity updates for each age cohort.
- At $t=0.2$ (with $dt=0.1$):
  - **H+ (Species 1, $\omega_c = 1.0$):**
    - Cohort 2 (0 steps): $v_x = 0.0, v_y = 0.0$
    - Cohort 1 (1 step): $v_x \approx 0.0049875$
    - Cohort 0 (2 steps): $v_x \approx 0.0199004$
    - Mean Diagnostics: $V_x \approx 0.0083$, $V_y \approx 0.0994$
  - **O+ (Species 2, $\omega_c = 0.0625$):**
    - Cohort 2 (0 steps): $v_x = 0.0, v_y = 0.0$
    - Cohort 1 (1 step): $v_x \approx 0.0000195$
    - Cohort 0 (2 steps): $v_x \approx 0.0000781$
    - Mean Diagnostics: $V_x \approx 0.0000326$, $V_y \approx 0.00625$
- The validation script checks that both species match their respective physical production rates and analytical discrete velocity averages within tight tolerances, and that electron counts correctly ensure charge balance.

