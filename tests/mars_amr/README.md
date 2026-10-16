# Mars Test with AMR near the Body (semi-implicit full PIC)

The setup of `tests/mars` (coarse base grid dx = 0.25 R_M, absorbing and
conducting body at 100 km altitude, 10-reaction photochemistry, dual-cell source
averages, `#PHOTOCHEMICALBOUNDARY` at 200 km) with two refinement levels: 0.125
R_M inside 2 R_M and 0.0625 R_M (212 km) in the shell 0.95-1.3 R_M around the
body. `#CONSTANTPPV T` with `minPPCPerDim 2` and `#RESAMPLING` keep 64/8/8
particles per cell. Runs 20 s; the checks are those of `tests/mars`.

Measured (8 ranks, ~1.3 h): net charge < 0.3% of the solar-wind H+, Eb/Eb0 = 1.43
(a 15-19 nT pileup layer next to the body, 4.5 nT on the coarse grid), ~16M
particles, 7.8 GB. A third level (106 km) does not fit in 15 GB: the node mass
matrix of the many small fine boxes alone needs ~9 GB.

Needs a 3D build with nLevMax >= 3 (`./Config.pl -amrex3d -lev=4 -u=Exo`).
Marked expensive (`--test=mars_amr` or `--all`).
