#!/usr/bin/env python3
"""Validator for the ionospheric body boundary (tests/bodyionosphere).

A sphere of radius R emits an isotropic Maxwellian (#BODYIONOSPHERE) into
vacuum, without fields. Collisionless effusion gives the escape rate
Gamma = n0 vth/sqrt(2 pi) 4 pi R^2 and the steady density
n(r) = n0/2 (1 - sqrt(1 - R^2/r^2)). The log provides the escape rate
d(mBodyInject - mBodyAbsorb)/dt and the mass in the box, Mpart0.

* escape rate / box mass (no normalization involved) vs analytic: 5%;
* n0 follows from one refill of the ghost layer (mBodyInject of step 1);
  the staircase body emits a little more than the sphere, so the escape
  rate and the box mass must exceed their analytic values by < 20 / 25%.
"""
import logging
import math

import numpy as np

logger = logging.getLogger(__name__)

R = 1.5        # body radius [code units]
HALF = 4.0     # half width of the box
VTH = math.sqrt(1.380649e-23 * 1.0e5 / 1.67262192e-27) / 1.0e5  # T = 1e5 K


def _box_integral():
    """Integral of (1/2)(1 - sqrt(1 - R^2/r^2)) over the box outside R."""
    g = np.linspace(-HALF, HALF, 241)
    c = 0.5 * (g[1:] + g[:-1])
    dv = (g[1] - g[0])**3
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    r = np.sqrt(x * x + y * y + z * z)
    f = np.where(r > R, 0.5 * (1 - np.sqrt(np.clip(1 - (R / r)**2, 0, 1))), 0)
    return f.sum() * dv


def _ghost_volume(nCell=24):
    dx = 2 * HALF / nCell
    c = -HALF + (np.arange(nCell) + 0.5) * dx
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    body = x * x + y * y + z * z < R * R
    pad = np.pad(body, 1, constant_values=False)
    inner = (pad[2:, 1:-1, 1:-1] & pad[:-2, 1:-1, 1:-1] & pad[1:-1, 2:, 1:-1] &
             pad[1:-1, :-2, 1:-1] & pad[1:-1, 1:-1, 2:] & pad[1:-1, 1:-1, :-2])
    return (body & ~inner).sum() * dx**3


def _non_finite(diags):
    """Reason string if a log value is NaN/inf (comparisons with NaN are false,
    so a tolerance check alone would pass a broken run), else None."""
    keys = sorted({k for e in (diags or []) for k, v in e.items()
                   if isinstance(v, float) and not math.isfinite(v)})
    if keys:
        return f"non-finite log values in {', '.join(keys[:6])}"
    return None


def validate_log(pic_diags=None, test_name=None):
    bad = _non_finite(pic_diags)
    if bad:
        return False, bad
    d = pic_diags
    if not d or len(d) < 100 or "mBodyInject" not in d[0]:
        return False, "missing log or mBodyInject column"

    reasons = []

    def check(cond, msg):
        logger.debug("    %s: %s", "OK  " if cond else "FAIL", msg)
        if not cond:
            reasons.append(msg)

    # Steady state: the last third of the run.
    a, b = d[2 * len(d) // 3], d[-1]
    escape = ((b["mBodyInject"] - b["mBodyAbsorb"]) -
              (a["mBodyInject"] - a["mBodyAbsorb"])) / (b["time"] - a["time"])
    mass = np.mean([e["Mpart0"] for e in d[2 * len(d) // 3:]])
    integral = _box_integral()
    sphere = VTH / math.sqrt(2 * math.pi) * 4 * math.pi * R**2

    # n0 in code units from one refill of the ghost layer (one per step),
    # whose volume follows from the 24^3 grid: body cells (center inside R)
    # that touch a cell outside the body.
    n0 = (d[2]["mBodyInject"] - d[1]["mBodyInject"]) / _ghost_volume()
    # The log is written with the ghost layer filled; its content is part of
    # Mpart0 but not of the plasma outside the body.
    mass -= n0 * _ghost_volume()

    ratio = (escape / mass) / (sphere / integral)
    check(abs(ratio - 1) < 0.05,
          f"escape rate / box mass = {ratio:.3f} x analytic (tol 5%)")

    escRatio = escape / (n0 * sphere)
    check(1.0 < escRatio < 1.2,
          f"escape rate = {escRatio:.3f} x n0 vth/sqrt(2pi) 4pi R^2 "
          f"(staircase surface: expected in (1, 1.2))")
    massRatio = mass / (n0 * integral)
    check(1.0 < massRatio < 1.25,
          f"box mass = {massRatio:.3f} x analytic effusion profile "
          f"(expected in (1, 1.25))")

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"
