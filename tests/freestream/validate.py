#!/usr/bin/env python3
"""Strict validator for the free-stream test.

The free-stream is an exact uniform steady state:
- Fast bulk flow ux across a 3D oblique magnetic field (Bx, By, Bz)
- Convective electric field E = -u x B
- Plasma beta ~ 1 (T = 314000 K) with thermal pressure p = nkT balancing magnetic pressure
- Full PIC includes both kinetic ions (Species 0) and electrons (Species 1)
- Hybrid PIC includes kinetic ions (Species 0) and fluid electrons

The test validates that the uniform equilibrium remains steady:
- Energy-log checks (validate_log): total particle kinetic energy Epart,
  per-species energies Epart0/Epart1, magnetic energy Eb, and electric energy Ee
  are conserved.
- Plot checks (validate_plot):
  * Bulk velocities (uxS0, and uxS1 if present) conserved; transverse velocities negligible.
  * Particle temperatures/pressures (pS0, and pS1 if present) preserved without numerical heating/cooling.
  * Normal magnetic field Bx strictly uniform and preserved.
  * Mean magnetic and electric fields preserved.
  * For Full PIC: comoving-upwind scheme ensures spatial spread of B and E stays small (<10%).
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

# Strict free-stream tolerances.
EPART_TOL = 0.02         # particle kinetic energy conserved to +/-2%
EB_TOL = 0.05            # magnetic energy conserved to +/-5%
EE_TOL = 0.05            # electric field energy conserved to +/-5%
PRESS_TOL = 0.05         # mean particle pressure/temperature conserved to +/-5%
VEL_TOL = 0.05           # bulk velocity preserved to +/-5%
BX_SPREAD_TOL = 0.05     # normal field Bx spatial spread <= 5% of mean
FULL_PIC_SPREAD_TOL = 0.10 # Full PIC transverse field spread <= 10% of magnitude
TRANSVERSE_VEL_MAX = 0.05  # transverse spurious velocity <= 5% of bulk flow


def set_run_dir(run_dir):
    """Mirror the runner's RUN_DIR into the shared hybrid helper."""
    import tests._shared.hybrid as _hyb
    _hyb.set_run_dir(run_dir)


def validate_log(pic_diags=None, test_name=None):
    """Strict energy-log checks: Epart, Eb, Ee, and species energies conserved."""
    if not pic_diags or len(pic_diags) < 2:
        return True, "Passed (no pic log)"

    first, last = pic_diags[0], pic_diags[-1]
    passed = True
    reasons = []

    # Total particle kinetic energy conservation.
    ep0, ep1 = first.get("Epart", 0.0), last.get("Epart", 0.0)
    if not (math.isfinite(ep0) and math.isfinite(ep1)):
        passed = False
        reasons.append("Epart not finite (NaN/Inf)")
    elif ep0 > 0:
        ratio = ep1 / ep0
        logger.debug("    Epart: %s -> %s (ratio %.5f)",
                     f"{ep0:.6e}", f"{ep1:.6e}", ratio)
        if abs(ratio - 1.0) > EPART_TOL:
            passed = False
            reasons.append(
                f"Epart ratio {ratio:.4f} not within [{1-EPART_TOL:.3f}, "
                f"{1+EPART_TOL:.3f}] (particle kinetic energy not conserved)")

    # Per-species kinetic energy conservation (e.g. Epart0 for ions, Epart1 for electrons).
    for s_idx in (0, 1):
        s_key = f"Epart{s_idx}"
        if s_key in first and s_key in last:
            es0, es1 = first[s_key], last[s_key]
            if es0 > 0 and math.isfinite(es1):
                s_ratio = es1 / es0
                logger.debug("    %s: %s -> %s (ratio %.5f)",
                             s_key, f"{es0:.6e}", f"{es1:.6e}", s_ratio)
                if abs(s_ratio - 1.0) > EPART_TOL:
                    passed = False
                    reasons.append(
                        f"{s_key} ratio {s_ratio:.4f} not within [{1-EPART_TOL:.3f}, "
                        f"{1+EPART_TOL:.3f}] (species {s_idx} energy drift)")

    # Magnetic energy conservation (uniform state => Eb constant).
    eb0, eb1 = first.get("Eb", 0.0), last.get("Eb", 0.0)
    if not math.isfinite(eb1):
        passed = False
        reasons.append("Eb not finite (NaN/Inf)")
    elif eb0 > 0:
        ratio = eb1 / eb0
        logger.debug("    Eb: %s -> %s (ratio %.5f)",
                     f"{eb0:.6e}", f"{eb1:.6e}", ratio)
        if abs(ratio - 1.0) > EB_TOL:
            passed = False
            reasons.append(
                f"Eb ratio {ratio:.4f} not within [{1-EB_TOL:.3f}, "
                f"{1+EB_TOL:.3f}] (magnetic field energy drift)")

    # Electric field energy conservation (convective field E = -u x B => Ee constant).
    ee0, ee1 = first.get("Ee", 0.0), last.get("Ee", 0.0)
    if not (math.isfinite(ee0) and math.isfinite(ee1)):
        passed = False
        reasons.append("Ee not finite (NaN/Inf)")
    elif ee0 > 0:
        ratio = ee1 / ee0
        logger.debug("    Ee: %s -> %s (ratio %.5f)",
                     f"{ee0:.6e}", f"{ee1:.6e}", ratio)
        if abs(ratio - 1.0) > EE_TOL:
            passed = False
            reasons.append(
                f"Ee ratio {ratio:.4f} not within [{1-EE_TOL:.3f}, "
                f"{1+EE_TOL:.3f}] (electric field energy drift)")

    if passed:
        return True, "Passed (strict: Epart, Eb, Ee conserved)"
    return False, "; ".join(reasons)


def _fs_load_out(out_file):
    """Load a .out frame: return ({VAR: col_idx}, float rows)."""
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


def _fs_col(vidx, rows, name):
    """Return the column array for *name* (or None if absent)."""
    i = vidx.get(name)
    if i is None or not rows or i >= len(rows[0]):
        return None
    return [r[i] for r in rows]


def validate_plot(test_name):
    """Strict plot checks: bulk velocities, temperatures, B and E fields."""
    import tests._shared.hybrid as _hyb
    plots_dir = os.path.join(_hyb.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [FS] No .out files (PostProc.pl not run?) -- skipping.")
        return True, "No .out files (skipped)"

    vidx0, rows0 = _fs_load_out(out_files[0])
    vidxl, rowsl = _fs_load_out(out_files[-1])
    if vidx0 is None or vidxl is None or not rows0 or not rowsl:
        return True, "Could not parse .out frames (skipped)"

    passed = True
    reasons = []
    is_full_pic = "FULL PIC" in (test_name or "").upper()

    # 1. Bulk velocity ux conservation for all present species (UXS0, UXS1).
    for s_name in ("UXS0", "UXS1"):
        ux0 = _fs_col(vidx0, rows0, s_name)
        uxl = _fs_col(vidxl, rowsl, s_name)
        if ux0 and uxl:
            mean0, meanl = sum(ux0) / len(ux0), sum(uxl) / len(uxl)
            logger.debug("    [FS] <%s>: %s -> %s", s_name, f"{mean0:.5f}", f"{meanl:.5f}")
            if abs(mean0) > 1e-12 and abs(meanl / mean0 - 1.0) > VEL_TOL:
                passed = False
                reasons.append(
                    f"bulk velocity <{s_name}> {mean0:.4f} -> {meanl:.4f} "
                    f"(>{VEL_TOL*100:.0f}% drift)")

    # 2. Transverse velocities remain negligible noise (< 5% of bulk flow).
    ux_ref = _fs_col(vidx0, rows0, "UXS0")
    u_scale = abs(sum(ux_ref) / len(ux_ref)) if ux_ref else 1.0
    for u_trans in ("UYS0", "UZS0", "UYS1", "UZS1"):
        ut = _fs_col(vidxl, rowsl, u_trans)
        if ut:
            mean_ut = sum(ut) / len(ut)
            if abs(mean_ut) > TRANSVERSE_VEL_MAX * u_scale:
                passed = False
                reasons.append(
                    f"spurious transverse velocity <{u_trans}> = {mean_ut:.4f} "
                    f"(>{TRANSVERSE_VEL_MAX*100:.0f}% of bulk flow)")

    # 3. Particle temperature / pressure conservation (pS0, pS1).
    # Verifies the mean pressure/temperature is preserved without numerical heating/cooling.
    for p_name in ("PS0", "PS1"):
        p0 = _fs_col(vidx0, rows0, p_name)
        pl = _fs_col(vidxl, rowsl, p_name)
        if p0 and pl:
            mean0, meanl = sum(p0) / len(p0), sum(pl) / len(pl)
            spread0 = max(p0) - min(p0)
            spreadl = max(pl) - min(pl)
            logger.debug("    [FS] <%s>: %s -> %s (spread %.4e -> %.4e)",
                         p_name, f"{mean0:.5e}", f"{meanl:.5e}", spread0, spreadl)
            if abs(mean0) > 1e-12 and abs(meanl / mean0 - 1.0) > PRESS_TOL:
                passed = False
                reasons.append(
                    f"particle pressure <{p_name}> {mean0:.4e} -> {meanl:.4e} "
                    f"(>{PRESS_TOL*100:.0f}% change: numerical heating/cooling)")
            # Spread check: verify shot noise doesn't blow up relative to initial sampling.
            max_allowed_spread = 2.0 * max(spread0, 0.05 * abs(mean0))
            if spreadl > max_allowed_spread:
                passed = False
                reasons.append(
                    f"particle pressure <{p_name}> spread grew excessively "
                    f"({spread0:.4e} -> {spreadl:.4e} > {max_allowed_spread:.4e})")

    # 4. Magnetic field components (BX, BY, BZ).
    bx0 = _fs_col(vidx0, rows0, "BX")
    b_scale = abs(sum(bx0) / len(bx0)) if bx0 else 1.0
    for b_comp in ("BX", "BY", "BZ"):
        b0 = _fs_col(vidx0, rows0, b_comp)
        bl = _fs_col(vidxl, rowsl, b_comp)
        if b0 and bl:
            mb0, mbl = sum(b0) / len(b0), sum(bl) / len(bl)
            spreadl = max(bl) - min(bl)
            logger.debug("    [FS] <%s>: %s -> %s (spread %.4f)",
                         b_comp, f"{mb0:.5f}", f"{mbl:.5f}", spreadl)
            if abs(mb0) > 0.05 * b_scale:
                if abs(mbl / mb0 - 1.0) > EB_TOL:
                    passed = False
                    reasons.append(f"mean magnetic field {b_comp} changed {mb0:.4f} -> {mbl:.4f}")
                if b_comp == "BX" and spreadl > BX_SPREAD_TOL * abs(mbl):
                    passed = False
                    reasons.append(
                        f"normal field {b_comp} not uniform (spread {spreadl:.4f} > "
                        f"{BX_SPREAD_TOL*100:.0f}% of mean)")
                elif is_full_pic and spreadl > FULL_PIC_SPREAD_TOL * abs(mbl):
                    passed = False
                    reasons.append(
                        f"{b_comp} not uniform (spread {spreadl:.4f} > "
                        f"{FULL_PIC_SPREAD_TOL*100:.0f}% of mean)")
            else:
                if abs(mbl) > 0.05 * b_scale:
                    passed = False
                    reasons.append(f"spurious magnetic field component <{b_comp}> = {mbl:.4f}")

    # 5. Electric field components (EX, EY, EZ).
    ez0 = _fs_col(vidx0, rows0, "EZ")
    e_scale = abs(sum(ez0) / len(ez0)) if ez0 else 1.0
    for e_comp in ("EX", "EY", "EZ"):
        e0 = _fs_col(vidx0, rows0, e_comp)
        el = _fs_col(vidxl, rowsl, e_comp)
        if e0 and el:
            me0, mel = sum(e0) / len(e0), sum(el) / len(el)
            spreadl = max(el) - min(el)
            logger.debug("    [FS] <%s>: %s -> %s (spread %.4f)",
                         e_comp, f"{me0:.5f}", f"{mel:.5f}", spreadl)
            if abs(me0) > 0.05 * e_scale:
                if abs(mel / me0 - 1.0) > EE_TOL:
                    passed = False
                    reasons.append(f"mean electric field {e_comp} changed {me0:.4f} -> {mel:.4f}")
                if is_full_pic and spreadl > FULL_PIC_SPREAD_TOL * abs(mel):
                    passed = False
                    reasons.append(
                        f"{e_comp} not uniform (spread {spreadl:.4f} > "
                        f"{FULL_PIC_SPREAD_TOL*100:.0f}% of mean)")
            else:
                if abs(mel) > 0.05 * e_scale:
                    passed = False
                    reasons.append(f"spurious electric field component <{e_comp}> = {mel:.4f}")

    if passed:
        return True, "Passed (strict: state stays uniform, bulk flow, pressure, B & E preserved)"
    return False, "; ".join(reasons)
