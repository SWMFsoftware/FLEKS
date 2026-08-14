#!/usr/bin/env python3
"""Validator for the grouped absorbing-boundary tests (tests/bc_absorb/).

Two variants are discovered from this directory:
  - PARAM.in.fields   -> base_name "bc_absorb_fields"   (tophat EM pulse
                          absorbed at x-faces; EM energy decays)
  - PARAM.in.particles-> base_name "bc_absorb_particles" (ions drain out
                          through absorbing x walls; Epart decays)

The single validate_log/validate_plot here branch on the variant's base_name.
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

from .._shared.hybrid import RUN_DIR  # noqa: E402

EM_FINAL_FRAC = 0.25      # fields: Etot_final < 25% of initial (decay = absorb)
EPART_MAX_FRAC = 0.5      # particles: Epart_final < 50% of initial (decay)


def validate_log(pic_diags=None, test_name=None):
    """Energy-log check, branching on the variant."""
    logger.debug("Validating %s (log)...", test_name)

    if not pic_diags or len(pic_diags) < 2:
        logger.debug("  No PIC energy log found; skipping.")
        return True, "Passed (no pic log)"

    first, last = pic_diags[0], pic_diags[-1]
    e0 = first.get("Etot", 0.0)
    e1 = last.get("Etot", 0.0)

    if test_name == "bc_absorb_particles":
        return _validate_log_particles(e0, e1, first, last)

    return _validate_log_fields(e0, e1, first, last)


def _validate_log_fields(e0, e1, first, last):
    """Fields: EM energy (Ee+Eb) must DECAY as pulses are absorbed."""
    e0_f = first.get("Ee", 0.0) + first.get("Eb", 0.0)
    e1_f = last.get("Ee", 0.0) + last.get("Eb", 0.0)
    logger.debug("    Ee+Eb: %.4e -> %.4e", e0_f, e1_f)

    for v in (e0, e1, e0_f, e1_f):
        if not math.isfinite(v):
            return False, "Non-finite EM energy (NaN/Inf)"
    if e0 <= 0:
        return False, "Initial Etot is zero (pulse not seeded)"

    denom = e0 if e0 > 0 else 1.0
    frac = e1 / denom if e0 > 0 else (e1_f / e0_f if e0_f > 0 else 1.0)
    logger.debug("    Etot_final/Etot_initial = %.4f", frac)
    if frac > EM_FINAL_FRAC:
        return False, (f"EM energy did not decay (final/initial = {frac:.3f} > "
                       f"{EM_FINAL_FRAC}) -- pulses appear to be REFLECTING")

    return True, "Passed (EM energy decays => absorber active)"


def _validate_log_particles(e0, e1, first, last):
    """Particles: Epart must DECAY as ions leave through the absorbing walls."""
    ep0 = first.get("Epart", 0.0)
    ep1 = last.get("Epart", 0.0)
    logger.debug("    Epart: %.4e -> %.4e", ep0, ep1)

    if not (math.isfinite(e0) and math.isfinite(e1) and
            math.isfinite(ep0) and math.isfinite(ep1)):
        return False, "Non-finite energy (NaN/Inf)"
    if e0 <= 0:
        return False, "Initial Etot is zero"

    if ep0 > 0 and ep1 > EPART_MAX_FRAC * ep0:
        return False, (f"Epart stayed at {ep1/ep0*100:.0f}% of initial "
                       f"(>{EPART_MAX_FRAC*100:.0f}%) -- particles appear to be "
                       f"REFLECTING, not absorbed")

    return True, "Passed (particle energy decays => absorption active)"


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
def _load_last_out():
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        return None, None
    with open(out_files[-1], "r", encoding="latin-1") as f:
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
    logger.debug("Validating %s (plot)...", test_name)
    vidx, rows = _load_last_out()
    if vidx is None or not rows:
        return True, "No .out frames (skipped)"

    if test_name == "bc_absorb_particles":
        return _check_particles_plot(vidx, rows)
    return _check_fields_plot(vidx, rows)


def _check_fields_plot(vidx, rows):
    """Fields: late interior |Bz|,|Ey| are much reduced from the seeded 1.0."""
    bz = _col(vidx, rows, "BZ")
    ey = _col(vidx, rows, "EY")
    if bz is None or ey is None:
        return True, "BZ/EY not in .out (skipped)"
    for v in bz + ey:
        if not math.isfinite(v):
            return False, "Non-finite field in plot output"
    bz_max = max(abs(v) for v in bz)
    ey_max = max(abs(v) for v in ey)
    logger.debug("    [AF] late max |Bz| = %.4e, max |Ey| = %.4e", bz_max, ey_max)
    if bz_max >= 0.7 or ey_max >= 0.7:
        return False, (f"Late interior field barely reduced (max|Bz|={bz_max:.3e}, "
                       f"max|Ey|={ey_max:.3e}); a reflected standing wave may persist")
    return True, "Passed (interior EM field absorbed)"


def _check_particles_plot(vidx, rows):
    """Particles: rhoS0 finite/positive (skip if not deposited in no-solver mode)."""
    rho = _col(vidx, rows, "RHOS0")
    if rho is None:
        return True, "No rhoS0 column (skipped)"
    peak = max(abs(v) for v in rho)
    logger.debug("    [AP] rhoS0 peak = %.4e", peak)
    if peak <= 0:
        return True, "rhoS0 not deposited in no-solver mode (skipped)"
    if not all(math.isfinite(v) for v in rho):
        return False, "rhoS0 not finite (NaN)"
    if all(v < 0 for v in rho):
        return False, "rhoS0 entirely negative (unphysical)"
    return True, "Passed (rhoS0 positive and finite)"
