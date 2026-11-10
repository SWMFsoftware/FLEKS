#!/usr/bin/env python3
"""Validator for the planar solar wind discontinuity BC tests (tests/bc_discontinuity/).

Variants:
  - PARAM.in    -> base_name "bc_discontinuity"    (2D tilted TD across x-lo and y-lo)
  - PARAM.in.rd -> base_name "bc_discontinuity_rd" (1D/2D RD with circular B_t rotation and Walen jump)
"""
import glob
import logging
import math
import os
import subprocess
import tempfile

from tests._shared import run_dir as _run_dir

logger = logging.getLogger(__name__)

set_run_dir = _run_dir.set_run_dir

MP_SI = 1.67262192e-27  # proton mass [kg]


def validate_log(pic_diags=None, test_name=None):
    """Verify finite, positive, bounded energies throughout the run."""
    if not pic_diags or len(pic_diags) < 2:
        return True, "Passed (no pic log)"

    first, last = pic_diags[0], pic_diags[-1]
    for key in ("Etot", "Ee", "Eb", "Epart"):
        v0, v1 = first.get(key, 0.0), last.get(key, 0.0)
        if not (math.isfinite(v0) and math.isfinite(v1)):
            return False, f"{key} is non-finite (NaN/Inf)"

    et0, et1 = first.get("Etot", 0.0), last.get("Etot", 0.0)
    if et0 <= 0.0 or et1 <= 0.0:
        return False, f"Non-positive total energy: {et0} -> {et1}"
    if et1 / et0 < 0.2 or et1 / et0 > 5.0:
        return False, f"Etot ratio {et1 / et0:.3f} outside [0.2, 5.0]"

    return True, "Passed (energies finite and bounded)"


def _load_out(out_file):
    """Load an ASCII .out file and return ({VAR: col_idx}, float_rows)."""
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


def _mean(vals):
    return sum(vals) / len(vals) if vals else 0.0


def _check_negative_startup_abort(param_path, bad_old, bad_new, expected_substr):
    """Verify FLEKS.exe aborts at startup when a jump condition is violated."""
    fleks_exe = os.path.abspath(os.path.join(_run_dir.RUN_DIR, "FLEKS.exe"))
    if not os.path.exists(fleks_exe):
        fleks_exe = os.path.abspath("bin/FLEKS.exe")
    if not os.path.exists(fleks_exe):
        return True, "FLEKS.exe not found for negative check"

    with open(param_path, "r", encoding="latin-1") as f:
        text = f.read()
    if bad_old not in text:
        return False, f"Substring '{bad_old}' not found in {param_path}"
    bad_text = text.replace(bad_old, bad_new, 1)

    with tempfile.TemporaryDirectory() as tmpdir:
        with open(os.path.join(tmpdir, "PARAM.in"), "w", encoding="latin-1") as f:
            f.write(bad_text)
        proc = subprocess.run(
            [fleks_exe], cwd=tmpdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=30
        )
        out = (proc.stdout or b"").decode("utf-8", errors="replace") + \
              (proc.stderr or b"").decode("utf-8", errors="replace")
        if proc.returncode == 0:
            return False, f"Expected FLEKS.exe to abort on invalid jump ({bad_new}), but it exited 0"
        if expected_substr.lower() not in out.lower():
            return False, (f"Expected abort message containing '{expected_substr}', "
                           f"got: {out[-300:]}")
    return True, "Startup validation correctly rejected invalid jump condition"


def _validate_td(out_files):
    vidx0, rows0 = _load_out(out_files[0])
    vidxl, rowsl = _load_out(out_files[-1])
    if not vidx0 or not rows0 or not vidxl or not rowsl:
        return False, "Failed to parse initial or final .out file"

    xi = vidxl.get("X")
    yi = vidxl.get("Y")
    bzi = vidxl.get("BZ")
    rhoi = vidxl.get("RHOS0")
    if None in (xi, yi, bzi, rhoi):
        return False, f"Missing X/Y/BZ/RHOS0 in .out columns: {list(vidxl.keys())}"

    # 1. At t = 0, the discontinuity is padded outside the domain, so the
    # entire domain interior must be in State 1 (Bz ~= +5 nT, n ~= 5 /cc).
    bz0_mean = _mean([r[vidx0["BZ"]] for r in rows0])
    if abs(bz0_mean - 5.0e-9) > 0.5e-9:
        return False, f"Initial frame Bz mean {bz0_mean:.3e} T != +5.0e-9 T (IC polluted by BC)"

    # 2. At t = 6 s, the 45-degree tilted TD (normal = (1,1,0)/sqrt(2)) has
    # entered across both x-lo and y-lo.
    # Upstream corner (x + y <= 2.5e5 m): State 2 (Bz ~= -3 nT, n ~= 10 /cc)
    # Downstream corner (x + y >= 12.5e5 m): State 1 (Bz ~= +5 nT, n ~= 5 /cc)
    up_rows = [r for r in rowsl if (r[xi] + r[yi]) <= 2.5e5]
    dn_rows = [r for r in rowsl if (r[xi] + r[yi]) >= 12.5e5]
    if not up_rows or not dn_rows:
        return False, "Empty upstream or downstream sample region in 2D TD plot"

    bz_up = _mean([r[bzi] for r in up_rows])
    bz_dn = _mean([r[bzi] for r in dn_rows])
    n_up = _mean([r[rhoi] / (1.0e6 * MP_SI) for r in up_rows])
    n_dn = _mean([r[rhoi] / (1.0e6 * MP_SI) for r in dn_rows])

    logger.debug("  [TD] upstream: Bz = %.3e T, n = %.2f /cc; downstream: Bz = %.3e T, n = %.2f /cc",
                 bz_up, n_up, bz_dn, n_dn)

    if bz_up > -2.0e-9:
        return False, f"Post-discontinuity Bz = {bz_up:.3e} T (expected ~ -3.0e-9 T)"
    if bz_dn < 4.0e-9:
        return False, f"Pre-discontinuity Bz = {bz_dn:.3e} T (expected ~ +5.0e-9 T)"
    if n_up < 7.5:
        return False, f"Post-discontinuity density n = {n_up:.2f} /cc (expected ~ 10.0 /cc)"
    if abs(n_dn - 5.0) > 1.5:
        return False, f"Pre-discontinuity density n = {n_dn:.2f} /cc (expected ~ 5.0 /cc)"

    # 3. Check 45-degree symmetry between the x-lo inflow wall and y-lo inflow wall
    x_wall_rows = [r for r in rowsl if r[xi] <= 1.5e5 and 2.5e5 <= r[yi] <= 5.5e5]
    y_wall_rows = [r for r in rowsl if r[yi] <= 1.5e5 and 2.5e5 <= r[xi] <= 5.5e5]
    if x_wall_rows and y_wall_rows:
        bz_xwall = _mean([r[bzi] for r in x_wall_rows])
        bz_ywall = _mean([r[bzi] for r in y_wall_rows])
        if abs(bz_xwall - bz_ywall) > 1.0e-9:
            return False, (f"Tilted TD x-lo vs y-lo asymmetry too large: "
                           f"Bz(x-wall)={bz_xwall:.3e} vs Bz(y-wall)={bz_ywall:.3e}")

    # 4. Negative startup check: unbalanced TD pressure must abort with hint
    ok, msg = _check_negative_startup_abort(
        os.path.join("tests", "bc_discontinuity", "PARAM.in"),
        "96110.182               T [K]",
        "20000.0                 T [K]",
        "total pressure is not balanced"
    )
    if not ok:
        return False, msg

    return True, (f"Passed (2D tilted TD: Bz {bz_dn*1e9:+.2f} -> {bz_up*1e9:+.2f} nT, "
                  f"n {n_dn:.1f} -> {n_up:.1f} /cc, startup check ok)")


def _validate_rd(out_files):
    vidx0, rows0 = _load_out(out_files[0])
    vidxl, rowsl = _load_out(out_files[-1])
    if not vidx0 or not rows0 or not vidxl or not rowsl:
        return False, "Failed to parse initial or final .out file"

    xi = vidxl.get("X")
    bxi = vidxl.get("BX")
    byi = vidxl.get("BY")
    bzi = vidxl.get("BZ")
    uyi = vidxl.get("UYS0")
    if None in (xi, bxi, byi, bzi, uyi):
        return False, f"Missing X/BX/BY/BZ/UYS0 in .out columns: {list(vidxl.keys())}"

    # 1. Initial frame is pure State 1 (By ~= +4 nT)
    by0_mean = _mean([r[vidx0["BY"]] for r in rows0])
    if abs(by0_mean - 4.0e-9) > 0.5e-9:
        return False, f"Initial frame By mean {by0_mean:.3e} T != +4.0e-9 T"

    # 2. Final frame: upstream (x <= 1.5e5 m) is State 2 (By ~= -4 nT, uy < 0),
    # downstream (x >= 8.5e5 m) is State 1 (By ~= +4 nT, uy > 0).
    up_rows = [r for r in rowsl if r[xi] <= 1.5e5]
    dn_rows = [r for r in rowsl if r[xi] >= 8.5e5]
    if not up_rows or not dn_rows:
        return False, "Empty upstream or downstream sample region in RD plot"

    by_up = _mean([r[byi] for r in up_rows])
    by_dn = _mean([r[byi] for r in dn_rows])
    uy_up = _mean([r[uyi] for r in up_rows])
    uy_dn = _mean([r[uyi] for r in dn_rows])

    if by_up > -2.5e-9:
        return False, f"RD State 2 By = {by_up:.3e} T (expected ~ -4.0e-9 T)"
    if by_dn < 2.5e-9:
        return False, f"RD State 1 By = {by_dn:.3e} T (expected ~ +4.0e-9 T)"
    if uy_up > -15.0e3:
        return False, f"RD State 2 uy = {uy_up:.3e} m/s (expected ~ -3.90e4 m/s)"
    if uy_dn < 15.0e3:
        return False, f"RD State 1 uy = {uy_dn:.3e} m/s (expected ~ +3.90e4 m/s)"

    # 3. Circular rotation inside the RD layer: where |By| < 2 nT, Bz must
    # rotate out-of-plane so |B| stays near 5 nT instead of collapsing to 3 nT.
    mid_rows = [r for r in rowsl if abs(r[byi]) < 2.0e-9]
    if not mid_rows:
        return False, "No cells found inside RD transition layer (|By| < 2 nT)"
    max_bz_mid = max(abs(r[bzi]) for r in mid_rows)
    mean_bmag = _mean([math.sqrt(r[bxi]**2 + r[byi]**2 + r[bzi]**2) for r in rowsl])
    if max_bz_mid < 2.0e-9:
        return False, (f"RD transition layer did not rotate B out-of-plane: "
                       f"max |Bz| = {max_bz_mid:.3e} T (expected ~ 4.0e-9 T)")
    if abs(mean_bmag - 5.0e-9) > 0.6e-9:
        return False, f"RD mean |B| = {mean_bmag:.3e} T drifted from 5.0e-9 T"

    # 4. Negative startup check: violating Walen relation must abort
    ok, msg = _check_negative_startup_abort(
        os.path.join("tests", "bc_discontinuity", "PARAM.in.rd"),
        "-39.01855               uy [km/s]",
        "39.01855                uy [km/s]",
        "Walen relation"
    )
    if not ok:
        return False, msg

    return True, (f"Passed (RD: By {by_dn*1e9:+.2f} -> {by_up*1e9:+.2f} nT, "
                  f"layer max |Bz| = {max_bz_mid*1e9:.2f} nT, mean |B| = {mean_bmag*1e9:.2f} nT)")


def validate_plot(test_name=None):
    plots_dir = os.path.join(_run_dir.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if len(out_files) < 2:
        return False, f"Expected at least 2 .out files in {plots_dir}, found {len(out_files)}"

    if test_name and "rd" in test_name.lower():
        return _validate_rd(out_files)
    return _validate_td(out_files)
