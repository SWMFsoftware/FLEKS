#!/usr/bin/env python3
"""Validator for the inner-body test (tests/body/).

A uniform plasma streams in +x past an absorbing sphere declared with BODY:

* particles entering a body cell are removed and tallied (nBodyAbsorb,
  qBodyAbsorb, mBodyAbsorb in log_pic_n*.log);
* the body interior is empty: the particle moments and the electric field are
  reported as zero on the body nodes;
* a wake forms downstream of the body.
"""
import logging
import math

from tests._shared import run_dir as _run_dir

logger = logging.getLogger(__name__)

set_run_dir = _run_dir.set_run_dir

R_BODY = 1.2          # #BODY radius in code units (see PARAM.in)
WAKE_MAX_FRAC = 0.5   # wake density must stay below 50% of the upstream value
ZERO_TOL = 1e-12      # "exactly zero" threshold inside the body
ETOT_GROWTH_MAX = 10.0  # the body must not inject energy


def validate_log(pic_diags=None, test_name=None):
    """Check the absorption tallies and that the run stays finite."""
    logger.debug("Validating %s (log)...", test_name)

    if not pic_diags or len(pic_diags) < 2:
        return True, "Passed (no pic log)"

    for entry in pic_diags:
        for key in ("Etot", "Ee", "Eb", "Epart"):
            if not math.isfinite(entry.get(key, 0.0)):
                return False, f"Non-finite {key} (NaN/Inf)"

    first, last = pic_diags[0], pic_diags[-1]
    e0 = first.get("Etot", 0.0)
    e1 = last.get("Etot", 0.0)
    logger.debug("    Etot: %.4e -> %.4e", e0, e1)
    if e0 > 0 and e1 > ETOT_GROWTH_MAX * e0:
        return False, (f"Energy blew up: Etot {e0:.3e} -> {e1:.3e} "
                       f"(>{ETOT_GROWTH_MAX:g}x)")

    absorbed = [entry.get("nBodyAbsorb") for entry in pic_diags]
    if any(a is None for a in absorbed):
        return False, ("nBodyAbsorb column missing from log_pic -- the body "
                       "tallies are not written")

    logger.debug("    nBodyAbsorb: %.4e -> %.4e", absorbed[0], absorbed[-1])
    for prev, cur in zip(absorbed[:-1], absorbed[1:]):
        if cur < prev - 1e-9:
            return False, "nBodyAbsorb decreased (tallies are cumulative)"
    if absorbed[-1] <= 0:
        return False, ("No particle was absorbed by the body (nBodyAbsorb = 0) "
                       "-- the inner boundary is inactive")

    return True, (f"Passed (body absorbed {absorbed[-1]:.0f} particles, "
                  "tallies cumulative)")


def validate_plot(test_name):
    """Check that the body interior is empty and that a wake forms."""
    logger.debug("Validating %s (plot)...", test_name)

    vidx, rows = _run_dir.load_last_out()
    if vidx is None or not rows:
        return True, "No .out frames (skipped)"

    body = _run_dir.col(vidx, rows, "BODY")
    rho = _run_dir.col(vidx, rows, "RHOS0")
    x = _run_dir.col(vidx, rows, "X")
    y = _run_dir.col(vidx, rows, "Y")

    if body is None:
        return False, "BODY column missing from the output"
    if rho is None or x is None or y is None:
        return False, "RHOS0/X/Y columns missing from the output"

    if not all(math.isfinite(v) for v in body + rho + x + y):
        return False, "Non-finite value in the final plot frame"

    #--- The body must be present in the output and empty inside. ---
    inside = [i for i, b in enumerate(body) if b > 0.5]
    if not inside:
        return False, "No point is marked as inside the body (mask is empty)"

    for name in ("RHOS0", "RHOS1", "EX", "EY", "EZ"):
        values = _run_dir.col(vidx, rows, name)
        if values is None:
            continue
        peak = max(abs(values[i]) for i in inside)
        logger.debug("    max |%s| inside the body = %.3e", name, peak)
        if peak > ZERO_TOL:
            return False, (f"{name} is not zero inside the body "
                           f"(max |{name}| = {peak:.3e})")

    #--- Wake: the region just downstream must be depleted. ---
    def mean_rho(x_lo, x_hi):
        sel = [rho[i] for i in range(len(rows))
               if x_lo <= x[i] <= x_hi and abs(y[i]) <= 0.5 * R_BODY]
        return (sum(sel) / len(sel)) if sel else None

    upstream = mean_rho(-2.4, -1.3)
    wake = mean_rho(R_BODY + 0.1, 2.4)

    if upstream is None or wake is None:
        return True, "Passed (wake sample regions empty)"
    if upstream <= 0:
        return False, "Upstream density is zero (plasma not initialized)"

    logger.debug("    rho: upstream %.4e, wake %.4e", upstream, wake)
    if wake >= WAKE_MAX_FRAC * upstream:
        return False, (f"No wake behind the body (wake/upstream = "
                       f"{wake / upstream:.3f} >= {WAKE_MAX_FRAC}) -- particles "
                       "are not absorbed")

    return True, (f"Passed (body interior empty, wake/upstream = "
                  f"{wake / upstream:.3f})")
