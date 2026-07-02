# Mars Chemistry Test (10 BATSRUS ModUserMars Reactions, 4 Ion Species)

## Description

This standalone test verifies the `#CHEMISTRY` command with all 10 chemical
reactions from the BATSRUS `ModUserMars.f90` implementation. It uses 4 ion
species (H+, O+, O2+, CO2+) and 3 neutral exosphere components (H, O, CO2).

This test exercises **cross-species charge exchange** — where a reactant ion
is consumed and a different product ion is produced (e.g., H+ + O -> O+ + H).
The product ion inherits the reactant ion's velocity and temperature, which
is the key extension over the existing `#CHARGEEXCHANGE` command.

## Chemistry Reactions

All 10 reactions from BATSRUS ModUserMars, with rate coefficients inflated
1e5x to produce detectable signals in the short 20-step simulation:

| # | Reaction | Type | Rate (BATSRUS) | Rate (test) |
|---|----------|------|----------------|-------------|
| 1 | CO2 + hv -> CO2+ + e- | Photoionization | 2.47e-7 s^-1 | 2.47e-2 s^-1 |
| 2 | O + hv -> O+ + e- | Photoionization | 8.89e-8 s^-1 | 8.89e-3 s^-1 |
| 3 | CO2+ + O -> O2+ + CO | Charge exchange | 1.64e-10 cm^3/s | 1.64e-5 cm^3/s |
| 4 | O+ + CO2 -> O2+ + CO | Charge exchange | 1.1e-9 cm^3/s | 1.1e-4 cm^3/s |
| 5 | CO2+ + O -> O+ + CO2 | Charge exchange | 9.60e-11 cm^3/s | 9.60e-6 cm^3/s |
| 6 | O2+ + e- -> O + O | Recombination | 7.38e-8 cm^3/s | 7.38e-3 cm^3/s |
| 7 | CO2+ + e- -> CO + O | Recombination | 3.1e-7 cm^3/s | 3.1e-2 cm^3/s |
| 8 | H+ + O -> O+ + H | Charge exchange | 5.084e-10 cm^3/s | 5.084e-5 cm^3/s |
| 9 | O+ + H -> H+ + O | Charge exchange | 6.4e-10 cm^3/s | 6.4e-5 cm^3/s |
| 10 | H + hv -> H+ + e- | Photoionization | 5.58e-8 s^-1 | 5.58e-3 s^-1 |

## Physics & Solver Setup

- **Plasma Species** (5 total):
  | Species | Mass [amu] | Charge [e] | Role |
  |---------|------------|------------|------|
  | 0 (e-)  | 0.01       | -1         | Electron (quasi-neutral) |
  | 1 (H+)  | 1.0        | +1         | Ion (reactant/product in CX) |
  | 2 (O+)  | 16.0       | +1         | Ion (reactant/product in CX) |
  | 3 (O2+) | 32.0       | +1         | Ion (product of CX, lost by recombination) |
  | 4 (CO2+) | 44.0      | +1         | Ion (product of photoionization, lost by CX/recombination) |

- **Exosphere**: Chamberlain profile with 3 components (H, O, CO2).

- **Temperature exponent = 0**: All reactions use constant rate coefficients
  (no Te dependence) to avoid the #UNIFORMSTATE pressure normalization issue
  in standalone tests.

## Expected Results

With the cross-species CX source term working correctly, the expected energy
changes over the 20-step simulation are:

1. **O2+ (Epart3) energy INCREASES significantly**: O2+ is produced by
   cross-species CX (reactions 3, 4) and consumed by recombination
   (reaction 6).  The CX source rate (~6.3 s^-1, driven by the large
   exosphere neutral density ~5e10 m^-3) vastly exceeds the recombination
   loss rate (~4e-17 s^-1, limited by the small plasma electron density
   in SI units).  Therefore O2+ energy must increase.  This is the critical
   validation for the cross-species CX source term.

2. **CO2+ (Epart4) energy changes**: CO2+ is produced by photoionization
   (reaction 1) and consumed by CX (reactions 3, 5) and recombination
   (reaction 7).  Both source and loss are active.

3. **H+ (Epart1) energy changes**: H+ is produced by photoionization
   (reaction 10) and CX (reaction 9), and consumed by CX (reaction 8).

4. **O+ (Epart2) energy changes**: O+ is produced by photoionization
   (reaction 2) and CX (reactions 5, 8), and consumed by CX (reactions 4, 9).

The validation checks that ALL 4 ion species show significant energy changes
(> 0.1%), and specifically requires O2+ to increase — proving the
cross-species CX source term is functional.

## Running

From the FLEKS root directory (requires compiled `bin/FLEKS.exe`):

```bash
python3 tests/validate_tests.py --test=chemistry
```
