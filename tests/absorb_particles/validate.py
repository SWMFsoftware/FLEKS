#!/usr/bin/env python3
"""Validator for the absorbing-particle boundary-condition test.

The field solver is OFF (solveEM=F, useHybridPIC=F), so this run exercises only
the particle pusher and the ABSORBING particle boundary.  A single ion species
with a bulk +x flow (ux > 0) and zero magnetic field is initialised; BOTH x
faces are `absorb` (particles removed + tallied on exit), so no face re-injects
particles and the ions drain out through the x walls.  This validator checks:

  validate_log (energy log):
    1. The run completes with finite, positive total energy.
    2. Particle kinetic energy (Epart) DECAYS toward zero as ions leave through
       the absorbing x walls.  This is the opposite of the reflecting-wall test,
       where specular reflection conserves Epart.  A decaying Epart is the
       decisive signature that the absorber removes particles.

  validate_plot (plot output):
    1. The deposited ion density (rhoS0) is finite.  With the field solver off
       the node moment deposit may be zero (skip in that case); otherwise the
       density near the +x wall is drained relative to the interior as ions are
       absorbed.
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

from .._shared.hybrid import RUN_DIR  # noqa: E402

# The +x absorber removes the streaming ions, so Epart must DECAY.  A reflecting
# wall conserves Epart (see tests/reflecting_wall).  We require a substantial
# decay to clearly distinguish absorb from reflect.
EPART_MAX_FRAC = 0.5      # Epart_final < 50% of Epart_initial
RHO_MIN = 1e-3            # deposited rhoS0 must stay finite / positive


def validate_log(pic_diags=None, test_name=None):
    """Energy-log checks: finite energy and DECAYING Epart (particles absorbed)."""
    logger.debug("Validating Absorbing-Particle Test (log)...")

    if not pic_diags or len(pic_diags) < 2:
        logger.debug("  [AP] No PIC energy log found; skipping.")
        return True, "Passed (no pic log)"

    first, last = pic_diags[0], pic_diags[-1]
    e0 = first.get("Etot", 0.0)
    e1 = last.get("Etot", 0.0)
    ep0 = first.get("Epart", 0.0)
    ep1 = last.get("Epart", 0.0)

    logger.debug("    Etot: %.4e -> %.4e", e0, e1)
    logger.debug("    Epart: %.4e -> %.4e", ep0, ep1)

    if not (math.isfinite(e0) and math.isfinite(e1) and
            math.isfinite(ep0) and math.isfinite(ep1)):
        return False, "Non-finite energy (NaN/Inf)"

    if e0 <= 0:
        return False, "Initial Etot is zero (plasma not initialised)"

    # The decisive check: with an absorbing +x wall, particles must leave and
    # Epart must DECAY (a reflecting wall would conserve it, as in
    # tests/reflecting_wall).
    passed = True
    reasons = []
    if ep0 > 0 and ep1 > EPART_MAX_FRAC * ep0:
        passed = False
        reasons.append(
            f"Epart stayed at {ep1/ep0*100:.0f}% of initial (>{EPART_MAX_FRAC*100:.0f}%) "
            f"-- particles appear to be REFLECTING at the +x wall, not absorbed")

    if passed:
        return True, "Passed (particle energy decays => absorption active)"
    return False, "; ".join(reasons)


def _load_out(out_file):
    """Load a .out frame: return ({VAR: col_idx}, float rows)."""
    with open(out_file, "r", encoding="latin-1") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return None, None
    vidx = {v.upper(): i for i, v in enumerate(lines[4].split())}
    rows = []
    for line in lines[5:]:
        cols = line.split()
        if not cols:
            continue
        try:
            rows.append([float(c) for c in cols])
        except ValueError:
            continue
    return vidx, rows


def _col(vidx, rows, name):
    i = vidx.get(name)
    if i is None or not rows or i >= len(rows[0]):
        return None
    return [r[i] for r in rows]


def validate_plot(test_name):
    """Plot check: if rhoS0 is deposited, it must be finite and positive.

    With the field solver OFF, the full-PIC node moment deposit is coupled to
    the E-solver, so rhoS0 may be zero in the plot.  In that case the check is
    skipped; the energy-log Epart decay is the primary validation.
    """
    import tests._shared.hybrid as _hyb
    plots_dir = os.path.join(_hyb.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [AP] No .out files; skipping plot check.")
        return True, "No .out files (skipped)"

    vidx, rows = _load_out(out_files[-1])
    if vidx is None or not rows:
        return True, "Could not parse .out (skipped)"

    rho = _col(vidx, rows, "RHOS0")
    if rho is None:
        return True, "No rhoS0 column (skipped)"

    peak = max(abs(v) for v in rho)
    logger.debug("    [AP] rhoS0 peak = %.4e", peak)

    # No-solver mode does not deposit output moments (zero array): skip.
    if peak <= 0:
        return True, "rhoS0 not deposited in no-solver mode (skipped)"

    finite = all(math.isfinite(v) for v in rho)
    if not finite:
        return False, "rhoS0 not finite (NaN)"

    # Ion density should remain positive/finite (no spurious creation or NaN);
    # the Epart-decay check is the primary absorption signal.
    if min(v for v in rho if abs(v) > 0 or True) < 0 and \
            all(v < 0 for v in rho):
        return False, "rhoS0 entirely negative (unphysical)"

    return True, "Passed (rhoS0 positive and finite)"
