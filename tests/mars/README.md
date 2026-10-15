# Mars Smoke Test (coarse 3D, semi-implicit full PIC)

## Description

The standalone Martian magnetosphere of `run_mars/PARAM.in.fullpic`
(5 species e-, H+, O2+, O+, CO2+; H/O/CO2 exosphere; 10-reaction
photochemistry with EUV optical depth; absorbing, conducting body at 100 km
altitude; Parker-spiral IMF) on a coarse 24 x 32 x 32 grid
(dx = 0.25 R_M), run for 20 s.

It is not a physics benchmark (the ionosphere is unresolved at this
resolution); it guards the setup against the bugs found while building it:

* `#INFLOW` blocks must start with a line named `rho`; a relabeled line made
  the heavy ions inherit the solar-wind block and flow in from upstream;
* the reactions must use the particle plasma state, so that O2+ (made only by
  ion-neutral chemistry) becomes the dominant ionospheric ion;
* the run must stay quasi-neutral and the field energies bounded.

Cost: ~4 min on 8 ranks. Marked expensive (`--test=mars` or `--all`).
