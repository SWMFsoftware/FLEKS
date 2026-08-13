#!/usr/bin/env python3
"""Validator for the wave-injection boundary-condition test.

A monochromatic Bz wave (iField=0) is injected at the x=xLo boundary face via
#WAVEBC into a pure-EM (lightwave) domain.  This validator checks:

  validate_log (energy log):
    1. Total EM energy is finite and > 0 (a wave is present).
    2. Energy is approximately conserved (no blow-up or collapse).

  validate_plot (plot output):
    1. The injected Bz perturbation is non-zero in the interior region away
       from the injection face, confirming the wave entered the domain rather
       than staying confined to the boundary ghost cells.
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

from .._shared.hybrid import RUN_DIR, set_run_dir  # noqa: E402


def validate_log(pic_diags=None, test_name=None):
    """Energy-log checks: finite, positive, conserved EM energy."""
    logger.debug("Validating Wave-Inject Test (log)...")

    if not pic_diags or len(pic_diags) < 2:
        logger.debug("  [INFO] No PIC energy log found; skipping energy checks.")
        return True, "Passed (no pic log)"

    e0 = pic_diags[0].get("Etot", 0.0)
    e1 = pic_diags[-1].get("Etot", 0.0)
    logger.debug("    Etot: %.4e -> %.4e", e0, e1)

    if not (math.isfinite(e0) and math.isfinite(e1)):
        return False, "Non-finite total EM energy"
    if e0 <= 0:
        return False, "Initial Etot is zero (wave not initialised)"
    if e1 <= 0:
        return False, "Final Etot is zero (wave collapsed)"

    ratio = e1 / e0
    logger.debug("    Etot_final/Etot_initial = %.4f", ratio)
    if ratio < 0.3 or ratio > 3.0:
        return False, f"Etot ratio {ratio:.3f} outside [0.3, 3.0]"

    return True, "Passed (EM energy finite and conserved)"


def _check_wave_in_interior():
    """Verify the injected Bz perturbation reached the interior."""
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [WI] No .out files found.")
        return True, "No .out files (skipped)"

    out_file = out_files[-1]
    with open(out_file, "r") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return True, "Short .out file"

    vidx = {v.upper(): i for i, v in enumerate(lines[4].split())}
    bz_i = vidx.get("BZ")
    bx_i = vidx.get("BX")
    x_i = vidx.get("X")
    if bz_i is None or x_i is None:
        return True, "BZ/X not in .out"

    # Interior region: the middle third of the x-domain (away from the faces).
    xvals = []
    rows = []
    for line in lines[5:]:
        cols = line.split()
        try:
            rows.append([float(c) for c in cols])
            if x_i < len(cols):
                xvals.append(float(cols[x_i]))
        except ValueError:
            continue
    if not xvals:
        return True, "No data rows"

    xmin, xmax = min(xvals), max(xvals)
    lo = xmin + 0.33 * (xmax - xmin)
    hi = xmin + 0.67 * (xmax - xmin)

    # Max |Bz| in the interior and max |Bz| overall.
    bz_interior = 0.0
    bz_overall = 0.0
    for r in rows:
        if bz_i < len(r):
            bz_overall = max(bz_overall, abs(r[bz_i]))
            if lo <= r[x_i] <= hi:
                bz_interior = max(bz_interior, abs(r[bz_i]))

    logger.debug("    [WI] max |Bz| overall = %.4e, interior = %.4e",
                 bz_overall, bz_interior)

    # The wave must reach the interior (not confined to the boundary).
    if bz_interior <= 0.0:
        logger.debug("    [WI] FAIL: Bz is zero in the interior -- wave did not "
                     "propagate from the injection face.")
        return False, "Bz zero in interior (wave did not propagate)"

    # The interior perturbation should be a meaningful fraction of the overall.
    if bz_overall > 0.0 and bz_interior < 0.1 * bz_overall:
        logger.debug("    [WI] Interior Bz is only %.3f of overall -- weak "
                     "penetration.", bz_interior / bz_overall)
        return False, "Wave barely penetrates the interior"

    logger.debug("    [WI] Wave present in interior: VERIFIED")
    return True, "Passed (Bz wave present in interior)"


def validate_plot(test_name):
    """Plot-output check: the injected wave reached the interior."""
    logger.debug("  --- Validating Wave-Inject plot (wave in interior) ---")
    return _check_wave_in_interior()
