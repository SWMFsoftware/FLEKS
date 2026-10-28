#!/usr/bin/env python3
"""Validator for the unresolved photoproduction test (tests/photoproduction).

An exponential atmosphere n0 exp(-(r - R)/H) with H = dx/33 surrounds a body
of radius R. The photoproduction per time is
  P = nu n0 4 pi (R^2 H + 2 R H^2 + 2 H^3).
The source sampled at the nodes would alias with the grid (a node either
hits the thin layer or not); with the dual-cell averages the injected mass
per step must match P within 10% for both grid offsets.

The mass unit is the uniform seed: Mpart1(0) = m n_seed V_seed, with V_seed
the volume of the cells outside the body.

The corona variant: a transparent Chamberlain corona n0 exp(-H0 (1/R - 1/r))
with the optical depth on. Only the geometric shadow of the planet (x < 0,
y^2 + z^2 < R^2) is dark, so P = nu * integral of n over the box outside the
body and the shadow (5%); a dark nightside would give about half of it.
"""
import logging
import math

import numpy as np

logger = logging.getLogger(__name__)

R, H = 1.5e5, 1.0e3        # [m]
N0, NU = 1.0e12, 1.0e-3    # [m^-3], [1/s]
NSEED = 1.0e3              # [m^-3]
DT = 0.1                   # [s]
NCELL, LNORM = 24, 1.0e5


def _seed_volume(lo):
    dx = 8.0 / NCELL
    c = lo + (np.arange(NCELL) + 0.5) * dx
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    outside = x * x + y * y + z * z >= (R / LNORM)**2
    return outside.sum() * (dx * LNORM)**3


H0_CORONA = 1.5e6           # [m]


def _corona_production(lo, n=320):
    """nu * integral of the corona over the box, outside body and shadow."""
    edges = np.linspace(lo, lo + 8.0, n + 1) * LNORM
    c = 0.5 * (edges[1:] + edges[:-1])
    dv = (edges[1] - edges[0])**3
    y, z = np.meshgrid(c, c, indexing="ij")
    perp2 = y * y + z * z
    total = 0.0
    for x in c:
        r = np.sqrt(x * x + perp2)
        lit = (r > R) & ~((x < 0) & (perp2 < R * R))
        nCor = N0 * np.exp(-H0_CORONA * (1 / R - 1 / np.maximum(r, R)))
        total += (nCor * lit).sum() * dv
    return NU * total


def _non_finite(diags):
    """Reason string if a log value is NaN/inf (comparisons with NaN are false,
    so a tolerance check alone would pass a broken run), else None."""
    keys = sorted({k for e in (diags or []) for k, v in e.items()
                   if isinstance(v, float) and not math.isfinite(v)})
    if keys:
        return f"non-finite log values in {', '.join(keys[:6])}"
    return None


def _validate_source_ppc(d):
    """sourceppc: with #CONSTANTPPV the cells carry a target particle count
    (2 x 2 x 2 here); #SOURCEPARTICLES 1 x 1 x 1 must still add one particle
    per cell and step (the target overrode it: 8 per cell). The source covers
    the 24^3 box; the cold seed barely moves, so the macroparticle count grows
    by about one per cell per step (5%)."""
    if "nParticles" not in d[0] or len(d) < 3:
        return False, "missing nParticles log column"
    perStep = [(b["nParticles"] - a["nParticles"]) for a, b in zip(d, d[1:])]
    ratio = sum(perStep) / len(perStep) / NCELL**3
    logger.debug("    source particles per cell and step: %.3f", ratio)
    if not math.isfinite(ratio) or abs(ratio - 1) > 0.05:
        return False, (f"source particles per cell and step = {ratio:.3f} "
                       "(expected 1, tol 5%)")
    return True, "Passed"


def validate_log(pic_diags=None, test_name=None):
    bad = _non_finite(pic_diags)
    if bad:
        return False, bad
    d = pic_diags
    if not d or len(d) < 3 or "Mpart1" not in d[0]:
        return False, "missing log or Mpart1 column"
    if test_name and test_name.endswith("sourceppc"):
        return _validate_source_ppc(d)
    lo = -3.9 if test_name and test_name.endswith("shifted") else -4.0

    seed = NSEED * _seed_volume(lo)          # seed ions
    unit = d[0]["Mpart1"] / seed             # code mass per ion
    perStep = [(b["Mpart1"] - a["Mpart1"]) / unit for a, b in zip(d, d[1:])]
    produced = np.mean(perStep) / DT
    corona = test_name and test_name.endswith("corona")
    if corona:
        P = _corona_production(lo)
    else:
        P = NU * N0 * 4 * math.pi * (R**2 * H + 2 * R * H**2 + 2 * H**3)
    tol = 0.05 if corona else 0.10
    ratio = produced / P
    logger.debug("    production %.4e /s vs analytic %.4e /s: ratio %.3f",
                 produced, P, ratio)
    if not math.isfinite(ratio) or abs(ratio - 1) > tol:
        return False, f"production = {ratio:.3f} x analytic (tol {tol:.0%})"
    return True, "Passed"
