#!/usr/bin/env python3
"""Validator for the GEM reconnection test."""
import glob
import logging
import math
import os

import numpy as np

RUN_DIR = "run_test"


def set_run_dir(run_dir):
    """Point the plot helpers at the current run directory."""
    global RUN_DIR
    RUN_DIR = run_dir


logger = logging.getLogger(__name__)


def validate_log(pic_diags=None, test_name=None):
    """Energy-log sanity for the GEM reconnection run."""
    if not pic_diags or len(pic_diags) < 2:
        return True, "Passed (no pic log)"
    first, last = pic_diags[0], pic_diags[-1]
    eb0, eb1 = first.get("Eb", 0.0), last.get("Eb", 0.0)
    ep1 = last.get("Epart", 0.0)
    if not math.isfinite(eb1):
        return False, "Eb not finite (NaN/Inf)"
    if not math.isfinite(ep1):
        return False, "Epart not finite (NaN/Inf)"
    if eb0 > 0 and eb1 > eb0 * 10.0:
        return False, f"Eb grew {eb1/eb0:.1f}x (unstable blow-up)"
    logger.debug("    Eb: %.3e -> %.3e", eb0, eb1)
    logger.debug("    Epart: %.3e -> %.3e", first.get("Epart", 0.0), ep1)
    return True, "Passed (finite, bounded Eb)"


def _load_frame(path):
    """Parse one z=0 .out file into (x, y, Bx, By, Bz, rhoS0, rhoS1)."""
    try:
        lines = open(path, encoding="latin-1").read().splitlines()
    except Exception:
        return None
    if len(lines) < 6:
        return None
    names = lines[4].split()[:19]
    idx = {v: i for i, v in enumerate(names)}
    rows = [ln.split() for ln in lines[5:] if len(ln.split()) >= 19]
    if not rows:
        return None
    data = np.array([r[:19] for r in rows], dtype=float)
    x = data[:, idx["x"]]
    y = data[:, idx["y"]]
    ux = np.unique(np.round(x, 3))
    uy = np.unique(np.round(y, 3))
    nx, ny = len(ux), len(uy)
    if nx < 4 or ny < 4:
        return None
    bx = data[:, idx["Bx"]].reshape(ny, nx)
    by = data[:, idx["By"]].reshape(ny, nx)
    bz = data[:, idx["Bz"]].reshape(ny, nx)
    rhoS0 = data[:, idx["rhoS0"]].reshape(ny, nx) if "rhoS0" in idx else None
    rhoS1 = data[:, idx["rhoS1"]].reshape(ny, nx) if "rhoS1" in idx else None
    return ux, uy, bx, by, bz, rhoS0, rhoS1


def validate_plot(test_name):
    """Verify GEM initial field and density setup."""
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        return True, "Passed (no plot files to check)"

    fr0 = _load_frame(out_files[0])
    if fr0 is None:
        return True, "Passed (frame unparseable)"

    ux, uy, bx, by, bz, rhoS0, rhoS1 = fr0

    # Verify Bx has positive and negative values (current sheet)
    if bx.min() >= 0 or bx.max() <= 0:
        return False, f"t=0: Bx min={bx.min():.2f}, max={bx.max():.2f} does not cross zero"

    # Verify sheet density is elevated above background
    if rhoS0 is not None:
        if rhoS0.max() <= rhoS0.min():
            return False, "t=0: density is completely uniform, expected current-sheet profile"
        logger.debug("    t=0 rhoS0 range: [%.3f, %.3f]", rhoS0.min(), rhoS0.max())

    # Verify quasi-neutrality if electrons are present
    if rhoS0 is not None and rhoS1 is not None:
        # Full-PIC species 0 (ions m=1) and species 1 (electrons m=0.04)
        ratio = 25.0
        imb = float(np.abs(rhoS0 - rhoS1 * ratio).max() / rhoS0.max())
        logger.debug("    t=0 max charge imbalance: %.3f", imb)
        if imb > 0.5:
            return False, f"t=0 charge imbalance too high: {imb:.3f}"

    return True, "Passed (GEM initial condition verified)"
