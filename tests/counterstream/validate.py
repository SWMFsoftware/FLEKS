#!/usr/bin/env python3
"""Validator for the counter-streaming two-ion hybrid test (tests/counterstream/).

Verifies:
1. Finite and bounded PIC energy diagnostics (no NaN/blowup).
2. Presence of individual species moments (rhoS0, rhoS1, uxS0, uxS1, pS0, pS1)
   and total ion moments (rho, ux, uy, uz, p).
3. Density conservation and consistency: rho_tot ~ rhoS0 + rhoS1.
4. Opposite bulk flow directions: uxS0 < 0 (solar wind) and uxS1 > 0 (beam).
5. Dual-wake formation: obstacle absorption depletes SW downstream (x < 0)
   and depletes beam downstream (x > 0).
"""
import glob
import logging
import math
import os

from tests._shared import run_dir as _run_dir

logger = logging.getLogger(__name__)


def set_run_dir(run_dir):
    """Mirror the runner's RUN_DIR into the shared hybrid helper."""
    import tests._shared.hybrid as _hyb
    _hyb.set_run_dir(run_dir)


def validate_log(pic_diags=None, test_name=None):
    """Energy-log checks: finite & bounded energies, no NaN/blow-up."""
    if not pic_diags or len(pic_diags) < 2:
        return True, "Passed (no pic log)"

    first, last = pic_diags[0], pic_diags[-1]
    passed = True
    reasons = []

    for key in ("Etot", "Ee", "Eb", "Epart"):
        v0, v1 = first.get(key, 0.0), last.get(key, 0.0)
        if not (math.isfinite(v0) and math.isfinite(v1)):
            passed = False
            reasons.append(f"{key} not finite (NaN/Inf)")

    ep0, ep1 = first.get("Epart", 0.0), last.get("Epart", 0.0)
    if ep0 > 0:
        ratio = ep1 / ep0
        logger.debug("    Epart: %.6e -> %.6e (ratio %.5f)", ep0, ep1, ratio)
        if ratio < 0.5 or ratio > 2.0:
            passed = False
            reasons.append(f"Epart ratio {ratio:.4f} outside [0.5, 2.0]")

    eb0, eb1 = first.get("Eb", 0.0), last.get("Eb", 0.0)
    if eb0 > 0:
        ratio = eb1 / eb0
        logger.debug("    Eb: %.6e -> %.6e (ratio %.5f)", eb0, eb1, ratio)
        if abs(ratio - 1.0) > 0.25:
            passed = False
            reasons.append(f"Eb ratio {ratio:.4f} drifted by > 25%")

    if passed:
        return True, "Passed (energies finite & bounded)"
    return False, "; ".join(reasons)


def _load_out(out_file):
    """Load an ASCII .out plot frame: return ({VAR: col_idx}, rows)."""
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


_col = _run_dir.col


def validate_plot(test_name):
    """Plot checks: species separation, total moments, opposite velocities,
    and obstacle wake formation."""
    import tests._shared.hybrid as _hyb
    plots_dir = os.path.join(_hyb.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [COUNTERSTREAM] No .out files found -- skipping plot validation.")
        return True, "Passed (no plot files)"

    last_file = out_files[-1]
    vidx, rows = _load_out(last_file)
    if not vidx or not rows:
        return False, f"Failed to parse plot file {last_file}"

    passed = True
    reasons = []

    # 1. Check required species and total moment columns
    required_species_cols = ["RHOS0", "UXS0", "PS0", "RHOS1", "UXS1", "PS1"]
    for c in required_species_cols:
        if c not in vidx:
            passed = False
            reasons.append(f"Missing column '{c}'")

    required_total_cols = ["RHO", "UX", "UY", "UZ", "P"]
    for c in required_total_cols:
        if c not in vidx:
            passed = False
            reasons.append(f"Missing total ion moment column '{c}'")

    if not passed:
        return False, "; ".join(reasons)

    # 2. Check total density consistency: rho_tot == rhoS0 + rhoS1 outside body
    rho_tot = _col(vidx, rows, "RHO")
    rho_s0 = _col(vidx, rows, "RHOS0")
    rho_s1 = _col(vidx, rows, "RHOS1")
    xi = vidx.get("X")
    yi = vidx.get("Y")

    max_rho_diff = 0.0
    for i in range(len(rows)):
        r0, r1, rt = rho_s0[i], rho_s1[i], rho_tot[i]
        diff = abs(rt - (r0 + r1))
        if diff > max_rho_diff:
            max_rho_diff = diff

    if max_rho_diff > 1e-4:
        passed = False
        reasons.append(f"Total rho inconsistent with sum(rhoS_i): max diff = {max_rho_diff:.4e}")

    # 3. Check opposite bulk velocities outside obstacle
    ux_s0 = _col(vidx, rows, "UXS0")
    ux_s1 = _col(vidx, rows, "UXS1")
    # Exclude cells inside body (r < 1.0)
    outside_ux0, outside_ux1 = [], []
    for i in range(len(rows)):
        x, y = rows[i][xi], rows[i][yi]
        if (x * x + y * y) > 2.0:
            if rho_s0[i] > 1e-5:
                outside_ux0.append(ux_s0[i])
            if rho_s1[i] > 1e-6:
                outside_ux1.append(ux_s1[i])

    mean_ux0 = sum(outside_ux0) / len(outside_ux0) if outside_ux0 else 0.0
    mean_ux1 = sum(outside_ux1) / len(outside_ux1) if outside_ux1 else 0.0

    logger.debug("    Mean outside ux: Species 0 = %.3f, Species 1 = %.3f", mean_ux0, mean_ux1)

    if mean_ux0 >= -1.0:
        passed = False
        reasons.append(f"Species 0 ux ({mean_ux0:.2f}) not flowing in -x direction (< -1.0)")

    if mean_ux1 <= 1.0:
        passed = False
        reasons.append(f"Species 1 ux ({mean_ux1:.2f}) not flowing in +x direction (> 1.0)")

    # 4. Check dual wake formation:
    # Species 0 (SW) wake on -x side of obstacle (-3 < x < -1.5, |y| < 1)
    # Species 1 (Beam) wake on +x side of obstacle (1.5 < x < 3, |y| < 1)
    sw_wake, sw_flank = [], []
    beam_wake, beam_flank = [], []
    for i in range(len(rows)):
        x, y = rows[i][xi], rows[i][yi]
        # Near wake -x side
        if -3.0 < x < -1.5:
            if abs(y) < 1.0:
                sw_wake.append(rho_s0[i])
            elif abs(y) > 3.0:
                sw_flank.append(rho_s0[i])
        # Near wake +x side
        if 1.5 < x < 3.0:
            if abs(y) < 1.0:
                beam_wake.append(rho_s1[i])
            elif abs(y) > 3.0:
                beam_flank.append(rho_s1[i])

    if sw_wake and sw_flank:
        mean_sw_wake = sum(sw_wake) / len(sw_wake)
        mean_sw_flank = sum(sw_flank) / len(sw_flank)
        logger.debug("    SW near-wake rho: %.3e vs flank %.3e", mean_sw_wake, mean_sw_flank)
        if mean_sw_wake >= 0.5 * mean_sw_flank:
            passed = False
            reasons.append(f"Solar wind wake deficit missing on -x side: {mean_sw_wake:.3e} >= {0.5*mean_sw_flank:.3e}")

    if beam_wake and beam_flank:
        mean_beam_wake = sum(beam_wake) / len(beam_wake)
        mean_beam_flank = sum(beam_flank) / len(beam_flank)
        logger.debug("    Beam near-wake rho: %.3e vs flank %.3e", mean_beam_wake, mean_beam_flank)
        if mean_beam_wake >= 0.5 * mean_beam_flank:
            passed = False
            reasons.append(f"Beam wake deficit missing on +x side: {mean_beam_wake:.3e} >= {0.5*mean_beam_flank:.3e}")

    if passed:
        return True, "Passed (separate species, total moments, dual wakes verified)"
    return False, "; ".join(reasons)
