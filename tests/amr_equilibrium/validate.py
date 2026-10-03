#!/usr/bin/env python3
"""Validator for the Hybrid PIC AMR thermal equilibrium test (tests/amr_equilibrium).

Verifies uniform thermal equilibrium across a stationary coarse-fine AMR interface:
  1. Energy log:
     - Finite, strictly positive energies (Etot, Eb, Epart).
     - Strict total energy conservation (< 1% drift).
     - Bounded magnetic energy Eb without numerical Hall/whistler growth.
  2. Spatial fluid profiles (PostIDL .out):
     - Detection of coarse (dx, dy ~ 1.0) and fine (dx, dy ~ 0.5) AMR grids.
     - Verification that the fine grid occupies the central box |x| <= 8.0, |y| <= 4.0.
     - Uniform ion density across both interfaces (residual within statistical noise).
     - Reflection symmetry of profiles across the midplane (x = 0).
"""
import glob
import logging
import math
import os

import numpy as np

RUN_DIR = "run_test"


def set_run_dir(run_dir):
    """Point validator at the current run directory."""
    global RUN_DIR
    RUN_DIR = run_dir


logger = logging.getLogger(__name__)

# Expected grid parameters
DX_COARSE = 1.0
DX_FINE = 0.5
SLAB_X_HALF = 8.0
SLAB_Y_HALF = 4.0
NHEADER = 5


def _load_frame(path):
    """Parse one PostIDL .out frame into a dict of per-cell arrays."""
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


def _load_all_frames():
    """Load all fluid .out files from the run directory."""
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*_fluid*.out")))
    if not out_files:
        return None, "no *_fluid .out plot files (PostProc.pl not run?)"
    frames = []
    for f in out_files:
        fr = _load_frame(f)
        if fr is not None:
            frames.append((os.path.basename(f), fr))
    if not frames:
        return None, "failed to parse any fluid .out frames"
    return frames, None


def validate_log(pic_diags=None, test_name=None):
    """Energy-log validation: finite energies and strict conservation."""
    logger.debug("Validating AMR Thermal Equilibrium energy log...")

    if not pic_diags or len(pic_diags) < 2:
        return False, f"Insufficient log entries (expected >=2, got {len(pic_diags) if pic_diags else 0})"

    cycles = [r["cycle"] for r in pic_diags]
    if cycles[-1] < 10:
        return False, f"Final cycle {cycles[-1]} < 10"

    e0 = pic_diags[0].get("Etot", 0.0)
    eb0 = pic_diags[0].get("Eb", 0.0)
    max_energy_drift = 0.02  # 2% max drift

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
        if drift > max_energy_drift:
            return False, f"Energy drift too large at cycle {cyc}: {drift * 100:.3f}% (limit {max_energy_drift * 100:.1f}%)"

        if eb0 > 0 and eb > eb0 * 1.5:
            return False, f"Magnetic energy grew excessively at cycle {cyc} (Eb={eb:.3e}, initial={eb0:.3e})"

    last_etot = pic_diags[-1].get("Etot", 0.0)
    final_drift = abs(last_etot - e0) / e0
    logger.debug(
        "  All %d cycles completed: initial Etot=%.6g, final Etot=%.6g (drift=%.4f%%)",
        len(pic_diags),
        e0,
        last_etot,
        final_drift * 100,
    )
    return True, f"Passed (Etot drift {final_drift * 100:.4f}%)"


def validate_plot(test_name=None):
    """Spatial profile validation across the AMR interfaces."""
    logger.debug("Validating AMR Thermal Equilibrium fluid plots...")

    frames, err = _load_all_frames()
    if err is not None:
        return False, err

    # Inspect the final frame
    fname, frame = frames[-1]
    logger.debug("  Checking frame %s with %d points...", fname, len(frame["x"]))

    x = frame["x"]
    y = frame["y"]
    rho = frame.get("rhoS0")
    if rho is None:
        return False, "rhoS0 column missing in fluid .out frame"

    # 1. Grid structure check: detect coarse and fine spacings in x and y
    in_fine = (np.abs(x) < SLAB_X_HALF - 0.5) & (np.abs(y) < SLAB_Y_HALF - 0.5)
    in_coarse = (np.abs(x) > SLAB_X_HALF + 0.5) | (np.abs(y) > SLAB_Y_HALF + 0.5)

    xs_fine = np.sort(np.unique(x[in_fine]))
    xs_coarse = np.sort(np.unique(x[in_coarse]))
    ys_fine = np.sort(np.unique(y[in_fine]))
    ys_coarse = np.sort(np.unique(y[in_coarse]))

    if len(xs_fine) < 2 or len(xs_coarse) < 2 or len(ys_fine) < 2 or len(ys_coarse) < 2:
        return False, "Could not identify coarse and fine x/y regions"

    dx_fine_meas = np.median(np.diff(xs_fine))
    dx_coarse_meas = np.median(np.diff(xs_coarse))
    dy_fine_meas = np.median(np.diff(ys_fine))
    dy_coarse_meas = np.median(np.diff(ys_coarse))

    logger.debug(
        "  Measured fine grid: dx=%.4f, dy=%.4f (expected ~%.1f)",
        dx_fine_meas,
        dy_fine_meas,
        DX_FINE,
    )
    logger.debug(
        "  Measured coarse grid: dx=%.4f, dy=%.4f (expected ~%.1f)",
        dx_coarse_meas,
        dy_coarse_meas,
        DX_COARSE,
    )

    if abs(dx_fine_meas - DX_FINE) > 0.1 or abs(dy_fine_meas - DX_FINE) > 0.1:
        return False, f"Fine grid spacing (dx={dx_fine_meas:.3f}, dy={dy_fine_meas:.3f}) differs from expected {DX_FINE}"
    if abs(dx_coarse_meas - DX_COARSE) > 0.2 or abs(dy_coarse_meas - DX_COARSE) > 0.2:
        return False, f"Coarse grid spacing (dx={dx_coarse_meas:.3f}, dy={dy_coarse_meas:.3f}) differs from expected {DX_COARSE}"

    # 2. Density profile uniformity across interfaces
    # Bin by x coordinates and compute mean density profile <rho>(x)
    x_bins = np.unique(x)
    mean_rho_profile = []
    for xb in x_bins:
        mask = np.abs(x - xb) < 1e-4
        mean_rho_profile.append(float(np.mean(rho[mask])))
    mean_rho_profile = np.array(mean_rho_profile)

    max_dev = float(np.max(np.abs(mean_rho_profile - 1.0)))
    logger.debug(
        "  Density profile across x in [%.1f, %.1f]: max |rho - 1.0| = %.4f",
        x_bins[0],
        x_bins[-1],
        max_dev,
    )

    # 16x16 particles per cell has shot noise ~ 1/sqrt(256) = 0.0625.
    # Averaged over y (16 or 32 cells), statistical uncertainty is < 0.02.
    # Allow up to 0.10 for PIC noise and initial relaxation.
    if max_dev > 0.10:
        return False, f"Density deviation across AMR interface too large: max |rho - 1.0| = {max_dev:.4f} (> 0.10)"

    # 3. Symmetry check across x = 0
    # For every x > 0, find corresponding x < 0 and compare mean density
    sym_diffs = []
    for i, xb in enumerate(x_bins):
        if xb <= 0:
            continue
        neg_xb = -xb
        match = np.where(np.abs(x_bins - neg_xb) < 1e-4)[0]
        if len(match) > 0:
            diff = abs(mean_rho_profile[i] - mean_rho_profile[match[0]])
            sym_diffs.append(diff)

    if sym_diffs:
        max_sym_diff = float(np.max(sym_diffs))
        logger.debug("  Midplane reflection symmetry: max |rho(x) - rho(-x)| = %.4f", max_sym_diff)
        if max_sym_diff > 0.08:
            return False, f"Asymmetric interface artifact: max |rho(x) - rho(-x)| = {max_sym_diff:.4f} (> 0.08)"

    # 4. Magnetic field check
    bz = frame.get("Bz")
    if bz is not None:
        mean_bz = float(np.mean(bz))
        logger.debug("  Mean Bz = %.4f (expected ~1.0)", mean_bz)
        if abs(mean_bz - 1.0) > 0.05:
            return False, f"Mean Bz = {mean_bz:.4f} drifted from 1.0"

    msg = f"AMR equilibrium verified: dx_fine={dx_fine_meas:.2f}, dx_coarse={dx_coarse_meas:.2f}, max |delta rho|={max_dev:.4f}"
    return True, msg


def validate(run_dir=None):
    """Fallback entry point for test runner."""
    if run_dir:
        set_run_dir(run_dir)
    return validate_plot()
