#!/usr/bin/env python3
"""Validator for tests/amr_dive_inflow (AMR + div(E) + inflow + co-moving B0)."""
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
COMOVING_MAX_ETOT_DRIFT = 0.05
COMOVING_MAX_EE_RATIO = 10.0
WARNING_PREFIX = "WARNING [Pic::SetTargetPPC]"
WARNING_NEEDLE = (
    f"{WARNING_PREFIX}: L1 refinement is only 2 L0 cell(s) "
    "from the inflow face (x+)"
)

logger = logging.getLogger(__name__)


def set_run_dir(run_dir):
    """Point the validator at the current run directory."""
    global RUN_DIR
    RUN_DIR = run_dir


def _is_comoving_b0(test_name):
    return bool(test_name and "comoving_b0" in test_name)


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

    is_comoving_b0 = _is_comoving_b0(test_name)
    drift_limit = COMOVING_MAX_ETOT_DRIFT if is_comoving_b0 else MAX_ETOT_DRIFT

    e0 = pic_diags[0].get("Etot", 0.0)
    ee0 = pic_diags[0].get("Ee", 0.0)
    if not math.isfinite(e0) or e0 <= 0.0:
        return False, f"Invalid initial Etot={e0}"

    max_ee = 0.0
    for r in pic_diags:
        cyc = r["cycle"]
        etot = r.get("Etot", 0.0)
        ee = r.get("Ee", 0.0)
        eb = r.get("Eb", 0.0)
        ep = r.get("Epart", 0.0)
        for name, val in [("Etot", etot), ("Ee", ee), ("Eb", eb), ("Epart", ep)]:
            if not math.isfinite(val):
                return False, f"Non-finite {name} at cycle {cyc}"
        if etot <= 0.0 or eb <= 0.0 or ep <= 0.0 or ee < 0.0:
            return False, (
                f"Invalid energy at cycle {cyc}: "
                f"Etot={etot}, Ee={ee}, Eb={eb}, Epart={ep}"
            )
        max_ee = max(max_ee, ee)
        drift = abs(etot - e0) / e0
        if drift > drift_limit:
            return False, (
                f"Total energy drift {drift * 100:.2f}% at cycle {cyc} "
                f"(limit {drift_limit * 100:.0f}%)"
            )

    if is_comoving_b0:
        ee_limit = COMOVING_MAX_EE_RATIO * max(ee0, 1e-6)
        if max_ee > ee_limit:
            return False, (
                f"Electric energy grew to {max_ee:.3e} (limit {ee_limit:.3e})"
            )

    stdout_path = os.path.join(RUN_DIR, "fleks_stdout.log")
    if not os.path.isfile(stdout_path):
        return False, f"Missing stdout log at {stdout_path}"
    with open(stdout_path, "r", encoding="utf-8", errors="replace") as f:
        stdout_text = f.read()

    if is_comoving_b0:
        if WARNING_PREFIX in stdout_text:
            return False, f"Unexpected clearance warning in stdout: {WARNING_PREFIX}"
        return True, "Co-moving B0 log verified (bounded Etot/Ee, no clearance warning)"

    if WARNING_NEEDLE not in stdout_text:
        return False, f"Expected warning not found in stdout: {WARNING_NEEDLE}"

    return True, "Log and SetTargetPPC warning verified"


def _validate_plot_comoving_b0(frame):
    """Validate non-periodic coarse-fine div(E), #BODY, and co-moving B0 fields."""
    required = (
        "X",
        "Y",
        "RHOS0",
        "RHOS1",
        "EX",
        "EY",
        "EZ",
        "B0Z",
        "QC",
        "BODY",
        "PPCS0",
        "PPCS1",
    )
    missing = [k for k in required if k not in frame]
    if missing:
        return False, f"Missing required plot columns: {', '.join(missing)}"

    x = frame["X"]
    y = frame["Y"]
    r = np.hypot(x, y)
    ne = frame["RHOS0"] / ME
    ni = frame["RHOS1"] / MI
    qc = frame["QC"]
    body = frame["BODY"]
    emag = np.sqrt(frame["EX"] ** 2 + frame["EY"] ** 2 + frame["EZ"] ** 2)

    inside = body > 0.5
    if not np.any(inside):
        return False, "No points inside #BODY mask (BODY > 0.5)"
    max_ne_in = float(np.max(np.abs(ne[inside])))
    max_ni_in = float(np.max(np.abs(ni[inside])))
    max_e_in = float(np.max(emag[inside]))
    if max_ne_in > 1e-10 or max_ni_in > 1e-10 or max_e_in > 1e-8:
        return False, (
            f"Non-zero fields/plasma inside #BODY: "
            f"max|ne|={max_ne_in:.2e}, max|ni|={max_ni_in:.2e}, max|E|={max_e_in:.2e}"
        )

    ext_inner = (r > 0.65) & (np.abs(x) < 1.85) & (np.abs(y) < 1.85)
    if not np.any(ext_inner):
        return False, "No exterior interior points (r > 0.65, |x|,|y| < 1.85)"
    max_qc_all = float(np.max(np.abs(qc)))
    max_qc_ext = float(np.max(np.abs(qc[ext_inner])))
    if max_qc_ext >= 8.0e4 or max_qc_all >= 3.0e6:
        return False, (
            f"Excessive net charge |QC|: exterior={max_qc_ext:.2e} (limit 8.0e4), "
            f"overall={max_qc_all:.2e} (limit 3.0e6)"
        )

    strips = [
        ("side", (np.abs(y) >= 0.80) & (np.abs(y) <= 1.25) & (np.abs(x) < 0.70), 1.6, 2.4),
        ("upstream", (x >= 0.80) & (x <= 1.35) & (np.abs(y) < 0.70), 1.6, 2.4),
        ("tail", (x >= -1.25) & (x <= -0.80) & (np.abs(y) < 0.60), 1.4, 2.2),
    ]
    strip_msgs = []
    for label, mask, n_lo, n_hi in strips:
        if not np.any(mask):
            return False, f"No points in {label} coarse-fine transition strip"
        mean_ne = float(np.mean(ne[mask]))
        mean_ni = float(np.mean(ni[mask]))
        mean_dn = float(np.mean(np.abs(ni[mask] - ne[mask])))
        max_qc_strip = float(np.max(np.abs(qc[mask])))
        logger.debug(
            "  %s strip: mean_ne=%.3f, mean_ni=%.3f, mean|ni-ne|=%.3f, max|QC|=%.1f",
            label,
            mean_ne,
            mean_ni,
            mean_dn,
            max_qc_strip,
        )
        if not (n_lo <= mean_ne <= n_hi and n_lo <= mean_ni <= n_hi):
            return False, (
                f"{label} strip density out of range [{n_lo}, {n_hi}]: "
                f"mean_ne={mean_ne:.3f}, mean_ni={mean_ni:.3f}"
            )
        if mean_dn >= 0.45 or max_qc_strip >= 2.0e4:
            return False, (
                f"{label} coarse-fine charge artifact: "
                f"mean|ni-ne|={mean_dn:.3f}, max|QC|={max_qc_strip:.1f}"
            )
        strip_msgs.append(f"{label} dn={mean_dn:.2f}")

    dayside = (x >= 0.60) & (x <= 1.05) & (np.abs(y) <= 0.25)
    if not np.any(dayside):
        return False, "No points in dayside subsolar band"
    mean_b0z = float(np.mean(frame["B0Z"][dayside]))
    mean_ey = float(np.mean(frame["EY"][dayside]))
    if mean_b0z >= -3.0 or mean_ey <= 100.0:
        return False, (
            f"Co-moving B0 motional Ey mismatch on dayside: "
            f"mean(B0Z)={mean_b0z:.2f} nT (< -3.0), mean(EY)={mean_ey:.1f} uV/m (> 100)"
        )

    up_ez = (x >= 1.15) & (x <= 1.75) & (np.abs(y) <= 0.50)
    if not np.any(up_ez):
        return False, "No points in upstream convection band"
    mean_ez = float(np.mean(frame["EZ"][up_ez]))
    if not (350.0 <= mean_ez <= 550.0):
        return False, (
            f"Upstream convection EZ={mean_ez:.1f} uV/m outside [350, 550]"
        )

    return True, (
        f"Verified: ext max|QC|={max_qc_ext:.1e}, {', '.join(strip_msgs)}, "
        f"dayside EY={mean_ey:.1f} uV/m, upstream EZ={mean_ez:.1f} uV/m"
    )


def validate_plot(test_name=None):
    """Verify upstream L0 collar PPC, boundary density, and L0/L1 smoothness."""
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "z=0*.out")))
    if not out_files:
        return False, f"No z=0*.out plot files found in {plots_dir}"

    frame = _load_frame(out_files[-1])
    if frame is None:
        return False, f"Failed to parse plot frame {out_files[-1]}"

    if _is_comoving_b0(test_name):
        return _validate_plot_comoving_b0(frame)

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

    max_qc_fine = float(np.max(np.abs(frame["QC"][fine_mask])))
    logger.debug("  Interior L1 strip max|QC| = %.1f", max_qc_fine)
    if max_qc_fine >= 2.0e4:
        return False, f"Interior L1 max|QC|={max_qc_fine:.1f} exceeds 2.0e4"

    collar = strip & (x > 1.0) & (x <= 1.5)
    if not np.any(collar):
        return False, "No points in upstream L0 collar (1.0, 1.5]"

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
        f"max|ni-ne| collar={dn_c:.3f} transition={dn_t:.3f}, fine max|QC|={max_qc_fine:.1f}"
    )
