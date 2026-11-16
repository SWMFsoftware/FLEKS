#!/usr/bin/env python3
"""Validator for tests/amr_dive_inflow (AMR + div(E) + tight inflow BC)."""
import glob
import logging
import math
import os

import numpy as np

RUN_DIR = "run_test"
NHEADER = 5
ME = 0.0625
MI = 1.0
MAX_ETOT_DRIFT = 0.15
WARNING_NEEDLE = (
    "WARNING [Pic::SetTargetPPC]: L1 refinement is only 2 L0 cell(s) "
    "from the inflow face (x+)"
)

logger = logging.getLogger(__name__)


def set_run_dir(run_dir):
    """Point the validator at the current run directory."""
    global RUN_DIR
    RUN_DIR = run_dir


def _load_frame(path):
    """Parse one PostIDL ASCII .out frame with upper-cased column keys."""
    with open(path, "r", encoding="latin-1") as f:
        lines = f.read().splitlines()
    if len(lines) <= NHEADER:
        return None
    names = [tok.upper() for tok in lines[NHEADER - 1].split()]
    rows = []
    for ln in lines[NHEADER:]:
        cols = ln.split()
        if len(cols) >= 5:
            rows.append(cols)
    if len(rows) < 4:
        return None
    # PostIDL appends scalar footer parameter names after column names on line 5.
    ncol = min(len(rows[0]), len(names))
    arr = np.array([r[:ncol] for r in rows], dtype=float)
    return {names[i]: arr[:, i] for i in range(ncol)}


def _per_x_means(x, ne, ni, mask):
    """Return per-x-plane transverse-averaged (ne, ni) arrays on *mask*."""
    xs = np.round(x[mask], 6)
    unique_x = np.unique(xs)
    ne_m = np.array([float(np.mean(ne[mask][xs == xv])) for xv in unique_x])
    ni_m = np.array([float(np.mean(ni[mask][xs == xv])) for xv in unique_x])
    return unique_x, ne_m, ni_m


def validate_log(pic_diags=None, test_name=None):
    """Verify energy boundedness and SetTargetPPC inflow clearance warning."""
    if not pic_diags or len(pic_diags) < 2:
        n_got = len(pic_diags) if pic_diags else 0
        return False, f"Insufficient log entries (expected >= 2, got {n_got})"

    e0 = pic_diags[0].get("Etot", 0.0)
    if not math.isfinite(e0) or e0 <= 0.0:
        return False, f"Invalid initial Etot={e0}"

    for r in pic_diags:
        cyc = r["cycle"]
        etot = r.get("Etot", 0.0)
        ee = r.get("Ee", 0.0)
        eb = r.get("Eb", 0.0)
        ep = r.get("Epart", 0.0)
        for name, val in [("Etot", etot), ("Ee", ee), ("Eb", eb), ("Epart", ep)]:
            if not math.isfinite(val):
                return False, f"Non-finite {name} at cycle {cyc}"
        if etot <= 0.0 or eb <= 0.0 or ep <= 0.0:
            return False, f"Non-positive energy at cycle {cyc}: Etot={etot}, Eb={eb}, Epart={ep}"
        drift = abs(etot - e0) / e0
        if drift > MAX_ETOT_DRIFT:
            return False, (
                f"Total energy drift {drift * 100:.2f}% at cycle {cyc} "
                f"(limit {MAX_ETOT_DRIFT * 100:.0f}%)"
            )

    stdout_path = os.path.join(RUN_DIR, "fleks_stdout.log")
    if not os.path.isfile(stdout_path):
        return False, f"Missing stdout log at {stdout_path}"
    with open(stdout_path, "r", encoding="utf-8", errors="replace") as f:
        stdout_text = f.read()
    if WARNING_NEEDLE not in stdout_text:
        return False, f"Expected warning not found in stdout: {WARNING_NEEDLE}"

    return True, "Log and SetTargetPPC warning verified"


def validate_plot(test_name=None):
    """Verify upstream L0 collar PPC, boundary density, and L0/L1 smoothness."""
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "z=0*.out")))
    if not out_files:
        return False, f"No z=0*.out plot files found in {plots_dir}"

    frame = _load_frame(out_files[-1])
    if frame is None:
        return False, f"Failed to parse plot frame {out_files[-1]}"

    required = ("X", "Y", "RHOS0", "RHOS1", "QC", "PPCS0", "PPCS1")
    missing = [k for k in required if k not in frame]
    if missing:
        return False, f"Missing required plot columns: {', '.join(missing)}"

    x = frame["X"]
    y = frame["Y"]
    ne = frame["RHOS0"] / ME
    ni = frame["RHOS1"] / MI
    ppc0 = frame["PPCS0"]
    ppc1 = frame["PPCS1"]

    strip = np.abs(y) < 0.75
    coarse_mask = strip & (x > 1.0)
    fine_mask = strip & (np.abs(x) < 1.0)
    if not np.any(coarse_mask) or not np.any(fine_mask):
        return False, "Missing expected coarse (x > 1.0) or fine (|x| < 1.0) points"

    collar = strip & (x >= 1.0) & (x <= 1.5)
    if not np.any(collar):
        return False, "No points in upstream L0 collar [1.0, 1.5]"

    mean_ppc0 = float(np.mean(ppc0[collar]))
    mean_ppc1 = float(np.mean(ppc1[collar]))
    logger.debug(
        "  Upstream collar mean PPC: ppcS0=%.1f, ppcS1=%.1f", mean_ppc0, mean_ppc1
    )
    if mean_ppc0 < 45.0 or mean_ppc1 < 45.0:
        return False, (
            f"Upstream collar pre-split PPC collapsed: "
            f"mean(ppcS0)={mean_ppc0:.1f}, mean(ppcS1)={mean_ppc1:.1f} (min 45)"
        )

    _, ne_c, ni_c = _per_x_means(x, ne, ni, collar)
    dn_c = float(np.max(np.abs(ni_c - ne_c)))
    if (
        np.min(ne_c) < 0.75
        or np.max(ne_c) > 1.25
        or np.min(ni_c) < 0.75
        or np.max(ni_c) > 1.25
        or dn_c > 0.20
    ):
        return False, (
            f"Upstream collar density/charge artifact: "
            f"ne=[{np.min(ne_c):.3f}, {np.max(ne_c):.3f}], "
            f"ni=[{np.min(ni_c):.3f}, {np.max(ni_c):.3f}], max|ni-ne|={dn_c:.3f}"
        )

    trans = (np.abs(y) < 0.95) & (x >= 0.5) & (x <= 1.25)
    if not np.any(trans):
        return False, "No points in L0/L1 transition [0.5, 1.25]"

    _, ne_t, ni_t = _per_x_means(x, ne, ni, trans)
    dn_t = float(np.max(np.abs(ni_t - ne_t)))
    logger.debug(
        "  L0/L1 transition [0.5, 1.25] x-means: ne=[%.3f, %.3f], ni=[%.3f, %.3f], max|ni-ne|=%.3f",
        np.min(ne_t), np.max(ne_t), np.min(ni_t), np.max(ni_t), dn_t,
    )
    if (
        np.min(ne_t) < 0.80
        or np.max(ne_t) > 1.20
        or np.min(ni_t) < 0.80
        or np.max(ni_t) > 1.20
        or dn_t > 0.25
    ):
        return False, (
            f"L0/L1 transition density/charge artifact: "
            f"ne=[{np.min(ne_t):.3f}, {np.max(ne_t):.3f}], "
            f"ni=[{np.min(ni_t):.3f}, {np.max(ni_t):.3f}], max|ni-ne|={dn_t:.3f}"
        )

    return True, (
        f"Verified: collar PPC=({mean_ppc0:.1f}, {mean_ppc1:.1f}), "
        f"max|ni-ne| collar={dn_c:.3f} transition={dn_t:.3f}"
    )
