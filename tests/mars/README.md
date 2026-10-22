# Mars Test (coarse 3D, semi-implicit full PIC)

## Description

The standalone Martian magnetosphere on a coarse 24 x 32 x 32 grid
(dx = 0.25 R_M) with the full boundary setup:

* absorbing, conducting body at 3590 km (200 km altitude), made floating by
  `#FLOATINGBODY` (the net electron charge it absorbs is re-emitted);
* rhybrid-style emission shell at 3790 km (`#IONOSPHEREEMISSION`: O2+
  2e25/s, O+ 1.4e25/s, 20000 K), photochemical boundary at the shell;
* `Mixed` exosphere with the rhybrid hot H and O coronae (pickup ions);
* inflow at +x, `outflowthermal` faces elsewhere.

Checks (log): quasi-neutral within 0.5% of N(H+) (measured 0.13%); the
re-emitted charge equals the sum of the per-step electron excesses of the
absorbed charge; the H+ budget from the face columns closes on the corona
production (+1.2 code mass/s); O2+ accumulates (emitted minus what falls
back); bounded field and particle energies. The bow shock and the dayside
pileup (30 - 40 nT) need 30 - 40 s and are not checked here.

`PARAM.in.fullpic_emission` is `run_mars/final_coarse` (stable to 200 s,
bow shock at 1.75 R_M from 40 s on) shortened to 10 s, with 2 x 2 x 2
particles per cell. The ionospheric layer (scale heights of 6.5-13 km) is
unresolved on 849 km cells: below the emission shell it is represented by the
boundary condition. The speed of light of a stand-alone run is uNormSI: it is
3000 km/s here, 10x the solar wind.

Cost: ~3 min on 8 ranks. Marked expensive (`--test=mars.fullpic_emission`).
