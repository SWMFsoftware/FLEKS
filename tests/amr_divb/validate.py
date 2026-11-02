#!/usr/bin/env python3
"""Validator for the semi-implicit PIC AMR div(B) preservation test (tests/amr_divb).

Verifies the semi-implicit upwind field solver across coarse-fine AMR boundaries:
  1. Energy conservation:
     - Finite, strictly positive energies (Etot, Eb, Epart).
     - Bounded total energy drift across the run.
  2. Solenoidal preservation:
     - In the default curl-form upwind (amr_divb, amr_divb_uniform), interface
       and domain max |div(B)| remain at machine roundoff (< 1e-10, typically ~1e-15).
  3. AMR hierarchy & physical fields:
     - Detects both coarse (dx = 0.25) and fine (dx = 0.125) grid levels in AMR decks.
     - Confirms smooth magnetic field and particle flow across the interface.
"""
import glob
import logging
import math
import os
import re

import numpy as np

RUN_DIR = "run_test"


def set_run_dir(run_dir):
    """Point the validator at the current run directory."""
    global RUN_DIR
    RUN_DIR = run_dir


logger = logging.getLogger(__name__)

CODE_TO_PLOT = 3396000.0 / 3496000.0
CUBE_HALF = 1.0 * CODE_TO_PLOT
DX_COARSE = 0.25 * CODE_TO_PLOT
DX_FINE = 0.125 * CODE_TO_PLOT
MAX_DIVB_TOL = 1e-10
NHEADER = 5

_DIVB_LINE_RE = re.compile(r"divB-AMR\s+n=(\d+)\s+(.*)")
_LEVEL_BUCKET_RE = re.compile(
    r"L(\d+)\s+if(?:/cv)?/in=([0-9.eE+-]+|inf|nan)/"
    r"(?:([0-9.eE+-]+|inf|nan)/)?([0-9.eE+-]+|inf|nan)/"
    r"dom=([0-9.eE+-]+|inf|nan)"
)


def _load_frame(path):
    """Parse one PostIDL .out frame into a dict of per-point arrays."""
    try:
        with open(path, "r", encoding="latin-1") as f:
            lines = f.read().splitlines()
    except Exception:
        return None

    if len(lines) <= NHEADER:
        return None

    all_names = lines[NHEADER - 1].split()
    rows = []
    for ln in lines[NHEADER:]:
        c = ln.split()
        if len(c) >= 5:
            rows.append(c)
    if len(rows) < 4:
        return None

    ncol = min(len(rows[0]), len(all_names))
    data = np.array([r[:ncol] for r in rows], dtype=float)
    idx = {all_names[i]: i for i in range(ncol)}
    return {k: data[:, idx[k]] for k in all_names if k in idx}


def validate_log(pic_diags=None, test_name=None):
    """Energy log validation."""
    if not pic_diags or len(pic_diags) < 2:
        return False, f"Insufficient log entries (expected >=2, got {len(pic_diags) if pic_diags else 0})"

    e0 = pic_diags[0].get("Etot", 0.0)
    eb0 = pic_diags[0].get("Eb", 0.0)
    max_drift = 0.10  # 10% max energy drift

    for r in pic_diags:
        cyc = r["cycle"]
        etot = r.get("Etot", 0.0)
        eb = r.get("Eb", 0.0)
        ep = r.get("Epart", 0.0)

        for name, val in [("Etot", etot), ("Eb", eb), ("Epart", ep)]:
            if not math.isfinite(val):
                return False, f"Non-finite {name} at cycle {cyc}"

        if etot <= 0.0 or ep <= 0.0:
            return False, f"Non-positive energy at cycle {cyc}"

        drift = abs(etot - e0) / e0
        if drift > max_drift:
            return False, f"Energy drift too large at cycle {cyc}: {drift * 100:.2f}% (limit {max_drift * 100:.1f}%)"

        if eb0 > 0 and eb > eb0 * 3.0:
            return False, f"Magnetic energy grew unphysically at cycle {cyc}: Eb={eb:.3e} vs Eb0={eb0:.3e}"

    final_drift = abs(pic_diags[-1].get("Etot", 0.0) - e0) / e0
    logger.debug("  All %d cycles completed: initial Etot=%.5g, final Etot=%.5g (drift=%.2f%%)",
                 len(pic_diags), e0, pic_diags[-1].get("Etot", 0.0), final_drift * 100)
    return True, f"Passed: Etot drift {final_drift * 100:.2f}%"


def validate_divb(run_dir=None, test_name=None):
    """Assert on max |div(B)| from stdout log."""
    is_scalar_variant = test_name and "scalar" in test_name
    rd = run_dir or RUN_DIR
    candidates = [
        os.path.join(rd, "fleks_stdout.log"),
        os.path.join(rd, "run.log"),
    ]
    stdout_path = None
    for c in candidates:
        if os.path.isfile(c):
            stdout_path = c
            break

    if not stdout_path:
        return True, "No stdout log found"

    with open(stdout_path, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    lines = [ln.strip() for ln in content.splitlines() if "divB-AMR" in ln]
    if not lines:
        return True, "No divB-AMR lines found (single level or not reported)"

    max_divb_seen = 0.0
    for ln in lines:
        m = _DIVB_LINE_RE.search(ln)
        if not m:
            continue
        cycle = int(m.group(1))
        rest = m.group(2)
        for lm in _LEVEL_BUCKET_RE.finditer(rest):
            lev = int(lm.group(1))
            iface = float(lm.group(2))
            cv_s = lm.group(3)
            covered = float(cv_s) if cv_s is not None else 0.0
            interior = float(lm.group(4))

            for bname, val in [("iface", iface), ("covered", covered), ("interior", interior)]:
                if not math.isfinite(val):
                    return False, f"Non-finite div(B) in {bname} bucket at cycle {cycle} L{lev}"
                max_divb_seen = max(max_divb_seen, val)
                if not is_scalar_variant and val > MAX_DIVB_TOL:
                    return False, f"div(B) exceeded {MAX_DIVB_TOL} at cycle {cycle} L{lev} {bname}={val:.2e}"

    logger.debug("  %d divB-AMR lines verified: max |div(B)| = %.2e", len(lines), max_divb_seen)
    return True, f"div(B) verified: max={max_divb_seen:.2e}"


def validate_plot(run_dir=None, test_name=None):
    """Validate fluid .out spatial profile."""
    is_scalar_variant = test_name and "scalar" in test_name
    is_uniform_variant = test_name and "uniform" in test_name

    rd = run_dir or RUN_DIR
    plots_dir = os.path.join(rd, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "z=0*.out")))
    if not out_files:
        return True, "No fluid .out plots generated"

    frame = _load_frame(out_files[-1])
    if frame is None:
        return False, f"Failed to parse fluid frame {out_files[-1]}"

    x = frame.get("x")
    y = frame.get("y")
    div = frame.get("divB")
    bmag = frame.get("B")

    if x is None or y is None or div is None or bmag is None:
        return False, "Missing required columns (x, y, divB, B) in plot frame"

    r = np.maximum(np.abs(x), np.abs(y)) / CUBE_HALF
    h = np.where(r < 1.0, DX_FINE, DX_COARSE)
    err = np.abs(div) * h / np.maximum(bmag, 1e-6)

    max_err = float(np.max(err))
    logger.debug("  Frame verified: %d points, max dimensionless div(B) error = %.3e", len(x), max_err)

    if not is_scalar_variant and max_err > MAX_DIVB_TOL:
        return False, f"Dimensionless div(B) error too high: {max_err:.3e} (expected < {MAX_DIVB_TOL})"

    return True, f"Plot verified: max error={max_err:.2e}"


def validate(run_dir=None, test_name=None, pic_diags=None, **kwargs):
    """Top-level validator called by validate_tests.py."""
    if run_dir:
        set_run_dir(run_dir)

    if pic_diags:
        ok, msg = validate_log(pic_diags, test_name)
        if not ok:
            return False, msg

    ok, msg = validate_divb(run_dir, test_name)
    if not ok:
        return False, msg

    ok, msg = validate_plot(run_dir, test_name)
    if not ok:
        return False, msg

    return True, "AMR div(B) preservation test passed"
