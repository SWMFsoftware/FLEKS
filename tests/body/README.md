# Inner body (`#BODY`)

A uniform plasma streams in +x through a periodic 2D box past an absorbing
sphere of radius 1.2 (code units) at the origin, declared with `#BODY`.

What is exercised:

* **Particle absorption** — a particle is removed (and tallied) as soon as it
  is pushed into a body cell. The boundary is the cell-based staircase, the
  same object as the grid mask, so nothing is silently discarded at the
  surface (see `PARAM.XML` for the trade-offs of this choice).
* **No creation inside the body** — the initial fill, the source particles and
  the boundary injection all skip body cells.
* **E = 0 inside the body** — the body nodes are excluded from the implicit
  E solve (Dirichlet), and `nodeE`/`nodeEth` are zeroed there after the last
  operation that can write into the body.
* **B is not forced to zero** — with E = 0 and no particles inside, the
  interior field simply stays at its initial value.
* **div(E) cleaning** — the residual of the div(E) correction is zeroed in the
  body cells so the correction does not move particles to compensate for the
  charge that the CIC tails of the surrounding plasma deposit inside.
* **Output** — `body` reports the mask, and the particle moments are reported
  as zero inside the body.

## Verification (`validate.py`)

1. `log_pic`: `nBodyAbsorb` (cumulative number of absorbed macroparticles,
   appended after the fixed columns) is present, strictly non-decreasing and
   non-zero at the end; all energies stay finite and `Etot` does not blow up.
2. `log_pic`: the GMRES solve is unaffected — the relative error stays at the
   configured tolerance.
3. plot: on every point with `body == 1`, `rhoS0`, `rhoS1`, `Ex`, `Ey`, `Ez`
   are exactly zero.
4. plot: the region just downstream of the body is depleted to below 50% of
   the upstream density, i.e. a wake forms.

## Notes

* The absorbed count is dominated by the electron thermal flux: the electrons
  of this deck are much faster than the bulk flow, so they hit the body at
  their thermal rate while the ions arrive with the bulk flow.
* Run with `python3 tests/validate_tests.py --test=body -n 4`.
