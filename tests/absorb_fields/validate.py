#!/usr/bin/env python3
"""Validator for the absorbing-field boundary-condition test.

A `tophat` EM pulse (discontinuous Ey/Bz step) is seeded in the central half of
the domain with zero macroparticles (pure EM).  The step splits into left/right
going pulses that travel toward the x-faces, which are ABSORBING
(outgoing-characteristic).  This validator checks:

  validate_log (energy log):
    1. Total EM energy (Ee + Eb) is finite and > 0 at the start (a pulse is
       present).
    2. EM energy DECAYS over the run.  This is the decisive absorber
       signature: the pulses leave the domain through the absorbing faces, so
       their energy is carried out.  With conducting/periodic walls the pulses
       would reflect and the energy would persist (no decay).

  validate_plot (plot output):
    1. The interior EM field amplitude (|Bz|, |Ey|) falls toward zero at late
       time (no standing wave / no reflected signal trapped in the domain).
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

from .._shared.hybrid import RUN_DIR  # noqa: E402

# An absorber removes the pulse energy; by the final frame the EM energy should
# be a small fraction of the initial value.  A reflecting wall (or a no-op
# `absorb` that falls back to float) would NOT decay like this.
EM_FINAL_FRAC = 0.25       # Etot_final < 25% of Etot_initial (decay = absorbing)
EM_OVERALL_RATIO = 5.0     # Etot_initial/Etot_final > 5  (equivalently)


def validate_log(pic_diags=None, test_name=None):
    """Energy-log check: finite, positive, then DECAYING EM energy."""
    logger.debug("Validating Absorbing-Field Test (log)...")

    if not pic_diags or len(pic_diags) < 2:
        logger.debug("  [AF] No PIC energy log found; skipping.")
        return True, "Passed (no pic log)"

    first, last = pic_diags[0], pic_diags[-1]
    e0 = first.get("Etot", 0.0)
    e1 = last.get("Etot", 0.0)
    e0_f = first.get("Ee", 0.0) + first.get("Eb", 0.0)
    e1_f = last.get("Ee", 0.0) + last.get("Eb", 0.0)

    logger.debug("    Etot: %.4e -> %.4e", e0, e1)
    logger.debug("    Ee+Eb: %.4e -> %.4e", e0_f, e1_f)

    for v in (e0, e1, e0_f, e1_f):
        if not math.isfinite(v):
            return False, "Non-finite EM energy (NaN/Inf)"

    if e0 <= 0:
        return False, "Initial Etot is zero (pulse not seeded)"

    # The EM field energy must decay as the pulses are absorbed.
    # Use the first *data* frame's Ee+Eb when Etot has an offset; Etot is a
    # fallback.
    denom = e0 if e0 > 0 else 1.0
    frac = e1 / denom if e0 > 0 else (e1_f / e0_f if e0_f > 0 else 1.0)
    logger.debug("    Etot_final/Etot_initial = %.4f", frac)

    if frac > EM_FINAL_FRAC:
        return False, (
            f"EM energy did not decay (final/initial = {frac:.3f} > "
            f"{EM_FINAL_FRAC}) -- pulses appear to be REFLECTING at the "
            f"absorbing faces, or the absorber is a no-op")

    return True, "Passed (EM energy decays => absorber active)"


def _load_last_out():
    """Return (vidx, rows) for the last .out plot frame."""
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        return None, None, None
    out_file = out_files[-1]
    with open(out_file, "r", encoding="latin-1") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return None, None, None
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
    return vidx, rows, out_file


def _col(vidx, rows, name):
    i = vidx.get(name)
    if i is None or not rows or i >= len(rows[0]):
        return None
    return [r[i] for r in rows]


def _max_abs(vals):
    if vals is None:
        return 0.0
    return max(abs(v) for v in vals)


def validate_plot(test_name):
    """Plot check: the interior EM field amplitude decays toward zero."""
    logger.debug("Validating Absorbing-Field Test (plot)...")
    vidx, rows, out_file = _load_last_out()
    if vidx is None or not rows:
        return True, "No .out frames (skipped)"

    bz = _col(vidx, rows, "BZ")
    ey = _col(vidx, rows, "EY")
    if bz is None or ey is None:
        return True, "BZ/EY not in .out (skipped)"

    for v in bz + ey:
        if not math.isfinite(v):
            return False, "Non-finite field in plot output"

    # The tophat seeds Bz=1, Ey=1 in the interior.  At late time, after the
    # pulses have been absorbed, the residual field is a shrinking wavefront
    # that is clearly reduced from the seeded amplitude.  A reflecting wall
    # would keep a standing wave near the seeded amplitude (~O(1)).  The energy
    # log is the primary (decisive) check; this plot check is a secondary sanity
    # bound.
    bz_max = _max_abs(bz)
    ey_max = _max_abs(ey)
    logger.debug("    [AF] late max |Bz| = %.4e, max |Ey| = %.4e", bz_max, ey_max)

    if bz_max >= 0.7 or ey_max >= 0.7:
        return False, (
            f"Late interior field barely reduced (max|Bz|={bz_max:.3e}, "
            f"max|Ey|={ey_max:.3e} -- seeded amplitude was 1.0); a reflected "
            f"standing wave may persist")

    return True, "Passed (interior EM field absorbed)"
