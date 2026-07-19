# 3D Light Wave in Periodic PIC with Adaptive Mesh Refinement

## Description

This test is a standalone migration of the coupled SWMF test
`Param/PARAM.in.test.FLEKS.AMR.LightWave.3D` (SWMF `Makefile.test` **test25**).

In the coupled run the GM/MHD component only supplied a uniform background
state and the PIC region; the PC component then initialised a transverse
electromagnetic (light) wave via `testCase = lightwave`.  This standalone
version reproduces the same physics directly: a vacuum transverse EM wave
(`nPartPerCell = 0`, so no macroparticles are loaded) propagating on a
**periodic AMR grid**.

## Physics & Solver Setup

- **Geometry & Boundaries**: 3D periodic box
  $x \in [-40, 40]$, $y \in [-30, 30]$, $z \in [-4, 4]$ (all in metres, with
  the normalized PIC units $l_{\rm norm}=1\,{\rm m}$, $u_{\rm norm}=1\,{\rm m/s}$,
  i.e. $c=1$).  Resolution $\Delta x = \Delta y = \Delta z = 1.0$.
- **Initial condition**: `fill_lightwaves(48.0)` fills the node E and B fields
  with an analytic transverse plane wave of wavelength $\lambda = 48$
  (normalized length units), propagating diagonally in the $x$–$y$ plane
  (wave-vector $\propto 0.6\hat{x} + 0.8\hat{y}$).  Both E and B have
  amplitudes up to $0.8$.
- **Time stepping**: fixed $\Delta t = 0.4$ (CFL $\approx 0.4$).
- **Discretization**: $\theta = 0.5$ (Crank–Nicolson-like), no explicit
  diffusion (${\rm coefDiff}=0$), divergence correction of E disabled
  (`doCorrectDivE = F`).
- **AMR**: a refinement region `region1` ($x\in[-20,20]$, $y\in[-15,15]$,
  $z\in[-2,2]$) is tagged for refinement at level 1 with
  `gridEfficiency = 1.0`.  Base blocks are $20\times20\times2$ so that several
  base blocks lie entirely inside `region1` and are refined.

## Requirements

AMR requires FLEKS to be built with `nLevMax >= 2` (the default is `nLevMax = 1`,
which disables refinement).  Enable it and rebuild before running:

```bash
cd PC/FLEKS && ./Config.pl -lev=2 && make -j4
```

With `nLevMax = 1` the `#REFINEREGION` command aborts at startup
(`iLev` must be smaller than the max level index), so the test cannot run
without AMR support.

## Expected Results

- The transverse EM wave propagates across the periodic domain and is
  refined inside `region1`.
- The total electromagnetic energy $E_{\rm tot} = E_e + E_b$ (written to
  `log_pic_n*.log`) is non-zero at $t=0$ and remains approximately conserved
  over the run (a vacuum Maxwell wave should not decay to zero or blow up).
- The plot output (`z=0` slice of `Ex Ey Ez Bx By Bz`) shows a non-zero field
  amplitude, i.e. the wave is present at the final time.

## Running

From the FLEKS root directory (requires `bin/FLEKS.exe` built with AMR):

```bash
python3 tests/validate_tests.py --test=lightwave
```

The validation checks (1) energy conservation from the `log_pic` energy log
($0.3 \le E_{\rm tot}^{\rm final}/E_{\rm tot}^{\rm initial} \le 3.0$, both
non-zero and finite) and (2) that the final plot frame contains a non-zero
magnetic-field amplitude (the wave is present).
