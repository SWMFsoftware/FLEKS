#!/usr/bin/env python3
"""Validator for constant particles per volume on AMR (tests/constantppv).

A uniform thermal plasma on a 3-level 3D AMR grid with #CONSTANTPPV and
#RESAMPLING. The split/merge goals must be the targets of every level, so the
number of macroparticles stays put (the base-grid goal split the refined
levels back up to 64 particles per cell, doubling them every step), and the
energy is conserved like in amr_equilibrium.
"""
import logging
import math

logger = logging.getLogger(__name__)


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
    if not d or len(d) < 10 or "nParticles" not in d[0]:
        return False, "missing log or nParticles column"
    n0 = d[0]["nParticles"]
    worst = max(abs(e["nParticles"] / n0 - 1) for e in d)
    logger.debug("    particle number: %d -> %d, worst change %.3f",
                 n0, d[-1]["nParticles"], worst)
    if worst > 0.15:
        return False, f"particle number changed by {worst * 100:.1f}% (tol 15%)"
    e0 = d[0]["Etot"]
    drift = max(abs(e["Etot"] / e0 - 1) for e in d)
    if not all(math.isfinite(e["Etot"]) for e in d) or drift > 0.01:
        return False, f"energy drift {drift * 100:.2f}% (tol 1%)"
    return True, "Passed"
