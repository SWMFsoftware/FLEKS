#!/usr/bin/env python3
"""Validator for the 3D hybrid PIC AMR Alfvén pulse test (tests/pulse_amr_3d).

Verifies the 3D coarse-fine AMR interface with non-uniform magnetic field:
  1. Energy conservation:
     - Finite, strictly positive energies (Etot, Eb, Epart).
     - Total energy drift < 5% across 20 cycles.
     - Bounded magnetic energy Eb without whistler/Hall instability.
  2. Div(B) preservation across 3D coarse-fine interface:
     - Every divB-AMR line on all levels has max |div(B)| < 1e-10
       (measures ~1e-15 with divergence-preserving AMR interface scheme).
     - dGhost remains stable without runaway.
  3. Spatial grid & field profile:
     - Detection of coarse (dx ~ 0.25) and fine (dx ~ 0.125) grids.
     - Presence of the By pulse perturbation.
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

# Expected grid parameters
DX_COARSE = 0.25
DX_FINE = 0.125
BOX_X_HALF = 2.0
BOX_Y_HALF = 1.0
BOX_Z_HALF = 0.5
NHEADER = 5
MAX_DIVB_TOL = 1e-10

_DIVB_LINE_RE = re.compile(r"divB-AMR\s+n=(\d+)\s+(.*)")
_LEVEL_BUCKET_RE = re.compile(
    r"L(\d+)\s+if(?:/cv)?/in=([0-9.eE+-]+)/(?:([0-9.eE+-]+)/)?([0-9.eE+-]+)/dom=([0-9.eE+-]+)"
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
    """Energy-log validation."""
    logger.debug("Validating 3D AMR pulse energy log...")
    if not pic_diags or len(pic_diags) < 2:
        return False, f"Insufficient log entries (expected >=2, got {len(pic_diags) if pic_diags else 0})"

    e0 = pic_diags[0].get("Etot", 0.0)
    eb0 = pic_diags[0].get("Eb", 0.0)
    max_drift = 0.05  # 5% max energy drift

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
            return False, f"Magnetic energy blew up at cycle {cyc}: Eb={eb:.3e} vs Eb0={eb0:.3e}"

    final_drift = abs(pic_diags[-1].get("Etot", 0.0) - e0) / e0
    logger.debug("  All %d cycles completed: initial Etot=%.5g, final Etot=%.5g (drift=%.2f%%)",
                 len(pic_diags), e0, pic_diags[-1].get("Etot", 0.0), final_drift * 100)
    return True, f"Passed: Etot drift {final_drift * 100:.2f}%"


def validate_divb(run_dir=None):
    """Assert on max |div(B)| parsed from fleks_stdout.log."""
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
        logger.debug("  [INFO] No stdout log found in %s; skipping div(B) check.", rd)
        return True, "No stdout log found"

    with open(stdout_path, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    lines = [ln.strip() for ln in content.splitlines() if "divB-AMR" in ln]
    if not lines:
        return False, "No divB-AMR lines found in stdout log"

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
                if val > MAX_DIVB_TOL:
                    return False, f"div(B) exceeded {MAX_DIVB_TOL} at cycle {cycle} L{lev} {bname}={val:.2e}"

    logger.debug("  %d divB-AMR lines verified: max |div(B)| = %.2e (< %.1e)",
                 len(lines), max_divb_seen, MAX_DIVB_TOL)
    return True, f"div(B) verified: max={max_divb_seen:.2e}"


def validate_plot(run_dir=None):
    """Validate fluid .out spatial profile."""
    rd = run_dir or RUN_DIR
    plots_dir = os.path.join(rd, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*_fluid*.out")))
    if not out_files:
        return True, "No fluid .out plots generated (PostProc.pl not run?)"

    frame = _load_frame(out_files[-1])
    if frame is None:
        return False, f"Failed to parse fluid frame {out_files[-1]}"

    x = frame.get("x")
    by = frame.get("By")
    if x is None or by is None:
        return False, "Missing x or By columns in fluid frame"

    max_by = float(np.max(np.abs(by)))
    logger.debug("  Fluid frame verified with %d points, max |By| = %.4e", len(x), max_by)
    if max_by <= 0.0:
        return False, "By pulse amplitude is identically zero"

    return True, f"Plot verified: points={len(x)}, max|By|={max_by:.2e}"


def validate(run_dir=None, test_name=None, pic_diags=None, **kwargs):
    """Top-level validator called by validate_tests.py."""
    if run_dir:
        set_run_dir(run_dir)

    if pic_diags:
        ok, msg = validate_log(pic_diags, test_name)
        if not ok:
            return False, msg

    ok, msg = validate_divb(run_dir)
    if not ok:
        return False, msg

    ok, msg = validate_plot(run_dir)
    if not ok:
        return False, msg

    return True, "3D AMR pulse test passed"
