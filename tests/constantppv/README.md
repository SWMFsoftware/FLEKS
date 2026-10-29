# Constant Particles per Volume on AMR

A uniform thermal plasma (full PIC, periodic, 16^3 cells) with a central
refined box (2 levels), `#CONSTANTPPV T` with `#MINPARTICLESPERCELL 2` (64 and 8
particles per cell) and `#RESAMPLING`. The splitting and merging goals must be
the targets of the levels: the particle number stays within 15% and the energy
within 1%.
