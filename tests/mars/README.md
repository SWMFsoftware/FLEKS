# Mars Smoke Test (coarse 3D, semi-implicit full PIC)

## Description

The standalone Martian magnetosphere of `run_mars/PARAM.in.fullpic`
(5 species e-, H+, O2+, O+, CO2+; H/O/CO2 exosphere; 10-reaction
photochemistry with EUV optical depth; absorbing, conducting body at 100 km
altitude; Parker-spiral IMF) on a coarse 24 x 32 x 32 grid
(dx = 0.25 R_M), run for 20 s.

It is not a physics benchmark: the ionospheric layer (scale heights of
6.5-13 km) is unresolved on 849 km cells. The sources are averaged over the
dual cell of every node, so the column production is right, and
`#PHOTOCHEMICALBOUNDARY` (200 km altitude) drops the ionization below the
photochemical boundary, where the ions recombine locally; injected into a
cell-diluted layer they would barely recombine and pile up. The test guards
the setup against the bugs found while building it:

* `#INFLOW` blocks must start with a line named `rho`; a relabeled line made
  the heavy ions inherit the solar-wind block and flow in from upstream;
* the reactions must use the particle plasma state, so that O2+ (made only by
  ion-neutral chemistry) becomes the dominant ionospheric ion;
* the run must stay quasi-neutral and the field energies bounded.

The speed of light of a stand-alone run is uNormSI: it is 3000 km/s here, 10x the
solar wind (the decks used 100 km/s at first, i.e. a flow 3x faster than light).

Cost: ~4 min on 8 ranks. Marked expensive (`--test=mars` or `--all`).

## fullpic_emission: the full boundary setup (10 s)

`PARAM.in.fullpic_emission` is `run_mars/final_coarse` (stable to 200 s,
bow shock at 1.75 R_M from 40 s on) shortened to 10 s, with 2 x 2 x 2
particles per cell:

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

Cost: ~3 min on 8 ranks. Marked expensive (`--test=mars.fullpic_emission`).
