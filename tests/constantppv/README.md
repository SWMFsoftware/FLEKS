# Constant Particles per Volume on AMR

A uniform thermal plasma (full PIC, periodic, 16^3 cells) with two nested
refined boxes (3 levels), `#CONSTANTPPV T` with `minPPCPerDim 2` (64, 8 and 8
particles per cell) and `#RESAMPLING`. The splitting and merging goals must be
the targets of the levels: the particle number stays within 15% and the energy
within 1%. Needs a 3D AMReX build with nLevMax >= 3 (`./Config.pl -lev=4`).
