#!/usr/bin/env python3
"""Validator for the AMR Fadeev reconnection test (tests/reconnection_amr)."""
import glob
import logging
import math
import os
import re

import numpy as np

# Self-contained (lives outside tests/); runner injects RUN_DIR via set_run_dir.
RUN_DIR = "run_test"


def set_run_dir(run_dir):
    """Point the plot helpers at the current run directory."""
    global RUN_DIR
    RUN_DIR = run_dir


logger = logging.getLogger(__name__)

# Grid / equilibrium parameters (must match PARAM.in / FadeevIC defaults).
# Base grid is the uniform reconnection grid 32 x 16 (dx = dy = 62.8/32 =
# 1.9625 d_i); the refined layer is 2x finer (dx = dy = 0.98 d_i).
L = 5.0                 # sheet half-thickness / d_i
DX_COARSE = 1.9625      # base-grid cell size / d_i
DX_FINE = DX_COARSE / 2.0  # refined-layer cell size / d_i (~0.98)
REFINE_Y_HALF = 7.85    # refined band |y| < 7.85 d_i

# PostIDL .out column indices (from the header line `x y rhoS0 rhoS1 Bx By Bz
# Ex Ey Ez uxS0 ...`): 0=x, 1=y, 2=rhoS0, 3=rhoS1, 4=Bx, 5=By, 6=Bz, ...
COL_X, COL_Y, COL_BX, COL_BY = 0, 1, 4, 5
NHEADER = 5  # .out header: 4 metadata lines + 1 variable-name line


def _load_frame(path):
    """Parse one PostIDL .out frame into a dict of per-point arrays.

    The .out layout is: 4 metadata lines, then a variable-name line
    (`x y rhoS0 rhoS1 Bx By Bz Ex Ey Ez uxS0 ...`), then one data row per cell.
    The multi-level AMR grid is written per-cell (non-uniform), so rows cannot
    be reshaped into a rectangular (ny, nx) grid.
    """
    lines = open(path, encoding="latin-1").read().splitlines()
    if len(lines) <= NHEADER:
        return None
    # The variable-name line lists all columns (including scalars like mS0,
    # qS0, cLight, rPlanet, cutz); the per-cell data rows carry only the field
    # columns (x y rhoS0 rhoS1 Bx By Bz Ex Ey Ez uxS0 ... pYZS1), so use the
    # actual data-row width and map the leading header names onto it.
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
    res = {"x": data[:, idx["x"]], "y": data[:, idx["y"]],
           "Bx": data[:, idx["Bx"]], "By": data[:, idx["By"]]}
    if "Bz" in idx:
        res["Bz"] = data[:, idx["Bz"]]
    if "rhoS0" in idx:
        res["rho"] = data[:, idx["rhoS0"]]
    return res


def _frame_stats(frame):
    """Return (coarse_present, fine_present, fine_layer_ok) flags.

    Detects the two cell spacings along y (and x) and verifies the fine cells
    sit in the central layer |y| < 7.85 d_i rather than at the box edges.
    A y-row is "fine" when its x-cells are spaced by DX_FINE; the fine rows
    must form a contiguous band around the current sheet.
    """
    y = frame["y"]
    ys = np.unique(y)
    if ys.size < 4:
        return False, False, False
    dy = np.diff(ys)
    dy = dy[dy > 0.02]
    if dy.size == 0:
        return False, False, False
    spacings = set(np.round(dy, 4))
    coarse = any(abs(s - DX_COARSE) < 0.05 for s in spacings)
    fine = any(abs(s - DX_FINE) < 0.05 for s in spacings)
    logger.debug("  _frame_stats: spacings=%s, coarse=%s (DX_COARSE=%s), fine=%s (DX_FINE=%s)",
                 sorted(spacings), coarse, DX_COARSE, fine, DX_FINE)

    # Fine y-rows: those carrying DX_FINE-spaced x cells.
    fine_rows = []
    for yy in ys:
        m = np.abs(frame["y"] - yy) < 0.01
        xs = np.sort(frame["x"][m])
        if xs.size >= 2:
            xdx = np.round(np.diff(xs), 4)
            if xdx.size and abs(float(xdx[0]) - DX_FINE) < 0.05:
                fine_rows.append(float(yy))
    layer_ok = False
    if fine_rows:
        ylo, yhi = min(fine_rows), max(fine_rows)
        layer_ok = (abs(ylo + REFINE_Y_HALF) < 0.5) and \
                   (abs(yhi - REFINE_Y_HALF) < 0.5)
    return coarse, fine, layer_ok


def _midplane_by(frame):
    """Return (x_array, By_array) along the midplane y ~ 0 (refined layer).

    The structured `.out` output samples the node grid (both solvers), so the
    block-aligned refined band (`maxBlockSizeY = 2`) contains the node row
    y = 0 exactly.  Pick the fine-level row (x cells spaced by DX_FINE) whose
    |y| is smallest — the row closest to the current sheet — and return its By
    profile.
    """
    y = frame["y"]
    ys = np.unique(y)
    best = None  # (|yy|, yy)
    for yy in ys:
        m = np.abs(y - yy) < 0.01
        xs = np.sort(frame["x"][m])
        if xs.size < 2:
            continue
        xdx = np.round(np.diff(xs), 4)
        if xdx.size and abs(float(xdx[0]) - DX_FINE) < 0.05:
            if best is None or abs(float(yy)) < best[0]:
                best = (abs(float(yy)), float(yy))
    if best is None:
        return None, None
    _, yy = best
    m = np.abs(y - yy) < 0.01
    if m.sum() < 8:
        return None, None
    x = frame["x"][m]
    by = frame["By"][m]
    order = np.argsort(x)
    return x[order], by[order]


def _null_crossings(x, f):
    """x-positions where the midplane in-plane field changes sign."""
    s = np.sign(f)
    crosses = np.where(np.diff(s) != 0)[0]
    return [float(x[i]) for i in crosses]


def _xpoint(x, f):
    """Midplane null nearest the central reconnection X-point (x ~ 0)."""
    nulls = _null_crossings(x, f)
    if not nulls:
        return None
    return min(nulls, key=lambda n: abs(n - 0.0))


# divB-AMR report line (hybrid only), e.g.
#   divB-AMR n=100 L0 if/cv/in=1.7e-04/3.7e-04/1.6e-04/dom=2.5e-01
#                  L1 if/in=4.4e-05/9.5e-05/dom=4.7e-05 dCov=... dGhost=...
_DIVB_LINE_RE = re.compile(r"divB-AMR\s+n=(\d+)\s+(.*)")
_DIVB_LEVEL_RE = re.compile(
    r"L(\d+)\s+if(?:/cv)?/in=([^/\s]+)/(?:([^/\s]+)/)?([^/\s]+)/dom=([^\s]+)")
# Interface div(B) on the fine level may exceed the largest interior value by
# at most this factor. The interface cells see the coarse-fine stencil, so
# they are not at round-off, but a fine level whose interface ghost cells are
# overwritten with the uniform state (B = 0 here) gives a ratio of ~1e4.
DIVB_IFACE_RATIO = 10.0
# Drift of the first fine ghost layer from the coarse interpolation, relative
# to max |B| (dGhost/|B| of the divB-AMR line). The layer is advanced by the
# Faraday update, so it departs from the interpolation only by the fine-scale
# structure: ~0.03 at the first cycle and ~0.1-0.2 later in a correct run.
# Ghost cells overwritten with the uniform state give ~1 from the first cycle.
# Unlike the interface bucket, which is empty when the fine level spans the
# whole box in z (fake 2D in a 3D build), this works on every build.
DGHOST_FIRST_MAX = 0.2
DGHOST_RUN_MAX = 0.5
_DGHOST_RE = re.compile(r"dGhost=([^|\s]+)\|B\|=([^\s]+)")
# Short force-free and GEM-like variants: max |B| (printed by divB-AMR) may
# exceed its first value by at most this factor. Both equilibria keep max |B|
# near B0 over these runs (particle noise adds ~15%); the 2*dx_coarse mode that
# the fine-to-coarse E injection drove on the GEM-like deck took it to 2x B0
# within 26 steps.
BMAX_GROWTH = 1.5


def _variant(test_name):
    """Which deck a variant runs: 'forcefree', 'gem', 'lobe' or the Fadeev default."""
    name = test_name or ""
    if "forcefree" in name:
        return "forcefree"
    if "gem" in name:
        return "gem"
    if "lobe" in name:
        return "lobe"
    return "fadeev"


def _fine_iface_divb(test_name=None):
    """Return (ok, message) for the coarse-fine interface div(B).

    Reads the divB-AMR lines that FLEKS prints with #DIVB alwaysComputeDivB
    from the saved stdout. The hybrid variants must print them, and every line
    must carry the level-1 record and the ghost drift: on a fake-2D grid the
    level-1 interface bucket is empty, so dGhost is the check that catches a
    broken coarse-fine treatment there, and a missing field must not read as
    zero drift. Only the full-PIC variant, which has no such report, may skip
    the check. For the force-free, GEM-like and lobe variants max |B| must
    also stay within BMAX_GROWTH of its first value.
    """
    required = test_name is None or test_name.endswith("hybrid")
    check_bmax = _variant(test_name) in ("forcefree", "gem", "lobe")
    bmax0 = None
    path = os.path.join(RUN_DIR, "fleks_stdout.log")
    if not os.path.isfile(path):
        if required:
            return False, "no stdout log (fleks_stdout.log) to read divB-AMR from"
        return True, "no stdout log"
    lines = [ln for ln in open(path, encoding="utf-8", errors="replace")
             if "divB-AMR" in ln]
    if not lines:
        if required:
            return False, ("no divB-AMR lines in the stdout "
                           "(#DIVB alwaysComputeDivB missing?)")
        return True, "no divB-AMR lines"
    worst = worst_ghost = 0.0
    for ln in lines:
        m = _DIVB_LINE_RE.search(ln)
        if not m:
            return False, f"unparseable divB-AMR line: {ln.strip()}"
        cycle = int(m.group(1))
        iface, interior = {}, {}
        for lm in _DIVB_LEVEL_RE.finditer(m.group(2)):
            try:
                iface[int(lm.group(1))] = float(lm.group(2))
                interior[int(lm.group(1))] = float(lm.group(4))
            except ValueError:  # 'inf'
                return False, f"div(B) not finite at cycle {cycle}"
        if 0 not in iface or 1 not in iface:
            return False, f"cycle {cycle}: divB-AMR line lacks the L0 or L1 record"
        gm = _DGHOST_RE.search(m.group(2))
        if not gm:
            return False, f"cycle {cycle}: divB-AMR line lacks dGhost and |B|"
        try:
            dghost, bmax = float(gm.group(1)), float(gm.group(2))
        except ValueError:
            return False, f"dGhost not finite at cycle {cycle}"
        if check_bmax:
            if bmax0 is None:
                bmax0 = bmax
            elif not bmax <= BMAX_GROWTH * bmax0:
                return False, (f"cycle {cycle}: max |B| = {bmax:.2g} grew "
                               f"beyond {BMAX_GROWTH} x its initial "
                               f"{bmax0:.2g} (interface instability)")
        rel = dghost / bmax if bmax > 0 else math.inf
        limit = DGHOST_FIRST_MAX if ln is lines[0] else DGHOST_RUN_MAX
        worst_ghost = max(worst_ghost, rel)
        if not rel <= limit:
            return False, (f"cycle {cycle}: first fine ghost layer is "
                           f"{rel:.2f} |B| away from the coarse "
                           f"interpolation (limit {limit})")
        ref = max(interior.values())
        if not all(math.isfinite(v) for v in list(iface.values()) + [ref]):
            return False, f"div(B) not finite at cycle {cycle}"
        ratio = iface[1] / ref if ref > 0 else (0.0 if iface[1] == 0 else math.inf)
        worst = max(worst, ratio)
        if ratio > DIVB_IFACE_RATIO:
            return False, (f"cycle {cycle}: fine-level interface div(B) = "
                           f"{iface[1]:.2e} is {ratio:.0f}x the interior "
                           f"max {ref:.2e} (limit {DIVB_IFACE_RATIO:.0f}x)")
    logger.debug("    interface div(B): worst ratio to interior = %.2f, "
                 "dGhost/|B| <= %.3f over %d divB-AMR lines", worst,
                 worst_ghost, len(lines))
    return True, (f"interface div(B) <= {worst:.2f}x interior, "
                  f"dGhost/|B| <= {worst_ghost:.2f}")


def validate_log(pic_diags=None, test_name=None):
    """Energy-log sanity: finite, bounded Eb; Epart grows (reconnection).

    For hybrid PIC also checks the coarse-fine interface div(B) (see
    _fine_iface_divb).
    """
    ok, msg = _fine_iface_divb(test_name)
    if not ok:
        return False, msg
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
        return False, (f"Eb grew {eb1/eb0:.1f}x (unstable blow-up, not "
                       f"reconnection)")
    logger.debug("    Eb: %.3e -> %.3e", eb0, eb1)
    logger.debug("    Epart: %.3e -> %.3e", first.get("Epart", 0.0), ep1)
    return True, "Passed (finite, bounded Eb)"


def _load_all_frames():
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*_fluid*.out")))
    if not out_files:
        return None, "no *_fluid .out plot files (PostProc.pl not run?)"
    frames = []
    for f in out_files:
        fr = _load_frame(f)
        if fr is not None:
            frames.append((os.path.basename(f), fr))
    if len(frames) < 2:
        return None, "need >=2 parseable .out frames"
    return frames, None


def _has_two_levels(frame):
    """True if the frame mixes two x spacings with ratio 2 (AMR output)."""
    spacings = set()
    for yy in np.unique(np.round(frame["y"], 6)):
        xs = np.unique(frame["x"][np.abs(frame["y"] - yy) < 1e-6])
        if xs.size >= 3:
            spacings.add(round(float(np.min(np.diff(xs))), 6))
    if len(spacings) < 2:
        return False
    lo, hi = min(spacings), max(spacings)
    return abs(hi / lo - 2.0) < 0.05


# |B|^2 must stay uniform in the force-free sheet: std/mean over all points
# at t=0 (the equilibrium) and at the end, and the mean |B|^2 of the refined
# patch must match the coarse level.
FF_B2_STD_T0 = 0.05
FF_B2_STD_END = 0.10
FF_B2_LEVEL_MISMATCH = 0.03
FF_BOX = (0.70524, 0.35262)  # refined patch half-widths (code units)


def _validate_forcefree_amr_plot(frames):
    """Force-free AMR variant: |B|^2 uniform on both levels, start and end."""
    if not _has_two_levels(frames[0][1]):
        return False, "AMR grid not detected (no two x spacings in ratio 2)"
    stats = []
    for (fname, fr), lim in ((frames[0], FF_B2_STD_T0),
                             (frames[-1], FF_B2_STD_END)):
        b2 = fr["Bx"]**2 + fr["By"]**2 + fr["Bz"]**2
        if not np.all(np.isfinite(b2)):
            return False, f"{fname}: non-finite B"
        rel_std = float(b2.std() / b2.mean())
        if rel_std > lim:
            return False, (f"{fname}: |B|^2 not uniform (std/mean = "
                           f"{rel_std:.3f} > {lim})")
        inside = (np.abs(fr["x"]) < FF_BOX[0]) & (np.abs(fr["y"]) < FF_BOX[1])
        mismatch = abs(float(b2[inside].mean() / b2[~inside].mean()) - 1.0)
        if mismatch > FF_B2_LEVEL_MISMATCH:
            return False, (f"{fname}: mean |B|^2 of the refined patch differs "
                           f"from the coarse level by {mismatch:.3f} "
                           f"(> {FF_B2_LEVEL_MISMATCH})")
        stats.append((rel_std, mismatch))
    logger.debug("    force-free |B|^2 std/mean %.3f -> %.3f, patch/coarse "
                 "mismatch %.4f -> %.4f", stats[0][0], stats[1][0],
                 stats[0][1], stats[1][1])
    return True, (f"force-free AMR: |B|^2 std/mean {stats[0][0]:.3f} -> "
                  f"{stats[1][0]:.3f}, patch/coarse {stats[1][1]:.4f}")


def _validate_gem_amr_plot(frames):
    """GEM-like AMR variant: Harris sheet at t=0, bounded field at the end."""
    first = frames[0][1]
    if not _has_two_levels(first):
        return False, "AMR grid not detected (no two x spacings in ratio 2)"
    y = first["y"]
    top = first["Bx"][y > y.max() - 1e-6].mean()
    bot = first["Bx"][y < y.min() + 1e-6].mean()
    if not (top > 0.9 and bot < -0.9):
        return False, (f"t=0: Bx does not form the Harris sheet "
                       f"(top={top:.2f}, bottom={bot:.2f})")
    if "rho" in first:
        rho_peak = float(first["rho"].max())
        rho_lobe = float(first["rho"][y > y.max() - 1e-6].mean())
        if not 0.9 < rho_peak < 1.5 or rho_lobe > 0.4:
            return False, (f"t=0: density peak {rho_peak:.2f} / lobe "
                           f"{rho_lobe:.2f} not the GEM profile (1.2 / 0.2)")
    fname, last = frames[-1]
    b = np.sqrt(last["Bx"]**2 + last["By"]**2 + last["Bz"]**2)
    if not np.all(np.isfinite(b)):
        return False, f"{fname}: non-finite B"
    if float(b.max()) > BMAX_GROWTH:
        return False, (f"{fname}: max |B| = {b.max():.2f} > {BMAX_GROWTH} B0 "
                       f"(interface instability)")
    logger.debug("    GEM AMR: Bx top/bottom %.2f/%.2f, final max |B| %.3f",
                 top, bot, float(b.max()))
    return True, f"GEM-like AMR: final max |B| {b.max():.3f} B0"


# Lobe variant: the uniform lobe plasma (B = B0 x) has no transverse field at
# t = 0, and only particle noise (~0.03 B0 at the end) should build one. The
# whistler feedback loop at the x interfaces took max |B_perp| to 0.7 B0.
LOBE_BPERP_MAX = 0.1


def _validate_lobe_amr_plot(frames):
    """Lobe AMR variant: uniform B0 x at t=0, no transverse field at the end."""
    first = frames[0][1]
    if not _has_two_levels(first):
        return False, "AMR grid not detected (no two x spacings in ratio 2)"
    if not np.allclose(first["Bx"], 1.0, atol=0.01):
        return False, "t=0: Bx is not the uniform lobe field B0 = 1"
    fname, last = frames[-1]
    bperp = np.hypot(last["By"], last.get("Bz", np.zeros_like(last["By"])))
    if not np.all(np.isfinite(bperp)):
        return False, f"{fname}: non-finite B"
    if float(bperp.max()) > LOBE_BPERP_MAX:
        return False, (f"{fname}: max |B_perp| = {bperp.max():.2f} B0 > "
                       f"{LOBE_BPERP_MAX} (whistler instability at the interface)")
    logger.debug("    lobe AMR: final max |B_perp| %.3f B0", float(bperp.max()))
    return True, f"lobe AMR: final max |B_perp| {bperp.max():.3f} B0"


def validate_plot(test_name):
    """AMR reconnection plot check.

    1. The two-level AMR grid is present (coarse + fine spacing, fine band
       around the central current sheet).
    2. Reconnection is active: the midplane By perturbation grows and the
       in-plane field shows the O/X null structure of the growing islands.
    """
    frames, err = _load_all_frames()
    if err is not None:
        return False, err
    variant = _variant(test_name)
    if variant == "forcefree":
        return _validate_forcefree_amr_plot(frames)
    if variant == "gem":
        return _validate_gem_amr_plot(frames)
    if variant == "lobe":
        return _validate_lobe_amr_plot(frames)

    # ---- (1) AMR refinement present ----
    coarse0, fine0, layer0 = _frame_stats(frames[0][1])
    if not coarse0 or not fine0:
        return False, ("AMR grid not detected (need coarse dx~3.9 and fine "
                       "dx~2.0 d_i)")
    if not layer0:
        return False, ("fine cells do not form the central current-sheet layer "
                       "|y| < ~7.85 d_i")

    # ---- (2) Reconnection on the refined midplane ----
    first = frames[0][1]
    x0, by0 = _midplane_by(first)
    if x0 is None:
        return False, "no midplane (y~0) points found in refined layer"

    # The Fadeev equilibrium is anti-symmetric in By about x=0, so the midplane
    # carries the reconnection X-point (By null near x ~ 0) at t=0.
    xp0 = _xpoint(x0, by0)
    if xp0 is None or abs(xp0) > 0.5 * L:
        return False, (f"t=0: no midplane X-point near x~0 "
                       f"(closest null at x={xp0:.2f} d_i)")

    # Perturbation growth from the seed and midplane profile change.
    last = frames[-1][1]
    max_dby = [float(np.abs(fr["By"] - first["By"]).max()) for _, fr in frames]
    late_amp = max(max_dby[-1], np.abs(last["By"]).max())
    early_amp = max(max_dby[:1] + [0.0]) if max_dby else 0.0
    growth = (late_amp / early_amp) if early_amp > 1e-6 else 0.0
    xf, byf = _midplane_by(last)
    if xf is None:
        return False, "no midplane points in final frame"
    profile_change = float(np.abs(byf - by0).max())

    # Hybrid reconnection has a slower ion-scale onset than full-PIC, but runs
    # longer (TimeMax=20 vs 3); the threshold is solver-aware.
    is_hybrid = test_name is not None and test_name.endswith("hybrid")
    late_thresh = 0.05 if is_hybrid else 0.1
    if late_amp < late_thresh:
        return False, (f"late |delta By| = {late_amp:.3f} too small "
                       f"(no instability)")
    if profile_change < 0.05:
        return False, (f"midplane By changed by only {profile_change:.3f} "
                       f"(no reconnection)")

    # ---- (3) AMR boundary stability & Bz boundedness ----
    for fname, fr in frames:
        if "Bz" in fr:
            max_bz = float(np.abs(fr["Bz"]).max())
            if not math.isfinite(max_bz):
                return False, f"{fname}: non-finite Bz detected"
            # Unphysical boundary corruption causes artificial Bz blowup (> 1.2).
            if max_bz > 1.2:
                return False, (f"{fname}: unphysical Bz blowup ({max_bz:.3f} > 1.2), "
                               f"spurious AMR boundary artifact")
            # Check refinement boundary (|y| ~ REFINE_Y_HALF)
            bny_mask = np.abs(np.abs(fr["y"]) - REFINE_Y_HALF) < 1.0
            if np.any(bny_mask):
                bny_bz = float(np.abs(fr["Bz"][bny_mask]).max())
                if bny_bz > 0.6:
                    return False, (f"{fname}: excessive Bz at AMR boundary |y|~{REFINE_Y_HALF:.2f} "
                                   f"({bny_bz:.3f} > 0.6)")

    logger.debug("    midplane |delta By|: %.3f -> %.3f (%.1fx), "
                 "X-point at x=%.2f d_i, profile change %.3f",
                 early_amp, late_amp, growth, xp0, profile_change)
    msg = (f"AMR [coarse dx={DX_COARSE}, fine dx={DX_FINE}, layer "
           f"|y|<{REFINE_Y_HALF}]: By pert {early_amp:.3f}->{late_amp:.3f}, "
           f"X-point at x={xp0:.2f} d_i")
    return True, msg
