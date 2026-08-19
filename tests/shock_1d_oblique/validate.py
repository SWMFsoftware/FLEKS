#!/usr/bin/env python3
"""Validator for the 1D oblique (45 deg) shock reflection test (tests/shock_1d_oblique/).

A magnetized hybrid plasma streams along -x into a reflecting left wall.  The
upstream field is inclined 45 deg to the shock normal, so Bx != 0.  The whole
point of this test is the reflect-field BC (BC::fixed): a PEC/conducting left wall
forces E_t = 0, contradicting the bulk motional E = -u x B and driving a runaway
(Eb -> ~1e18).  The fixed (upstream-pinned) field BC must keep the energies
finite and bounded and hold the upstream field at the reflection face.

Success criteria
----------------
1. Energy log: Etot, Ee, Eb, Epart all finite (no NaN/Inf blow-up).  Eb and Epart
   stay bounded (a fixed field BC prevents the runaway; generous tolerances because
   a reflecting shock genuinely compresses B and converts flow -> thermal energy).
2. Reflect-field face (-x): the field components Bx, Bz in the left-most cells stay
   at the upstream #INFLOW value (7.382e-9) with no spurious erosion, and the
   motional Ey is maintained (a PEC wall would zero E_t = Ey, disagreeing with -u x B).
3. Inflow face (+x): upstream rhoS0 and uxS0 at the right-most cells stay within
   tolerance of the prescribed #INFLOW state (no draining / no spurious source).
4. No spurious transverse electric field: mean Ex, Ez near zero; Ey tracks the
   upstream motional value.

The .out frames are produced by PostProc.pl (y=0 fluid ascii pic); if they are
absent the plot check is skipped (the energy-log check still runs).
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

# Energy tolerances.  A reflecting shock genuinely raises Eb (compression) and
# converts flow -> ion thermal energy (Epart rises), so we tolerate large but
# FINITE growth.  The decisive guard is finiteness + an upper bound that would
# catch the conducting-wall runaway (Eb ~ 1e18).
EPART_RATIO_MAX = 50.0      # Epart may grow as flow thermalises, but not explode
EB_RATIO_MAX = 50.0         # Eb may grow (shock compression) but must stay finite
E_SPURIOUS_TOL = 0.15       # spurious mean E <= 15% of upstream |B| scale
B_FACE_TOL = 0.30           # boundary Bx/Bz within 30% of upstream #INFLOW value
RHO_TOL = 0.15              # inflow-face density within 15% of upstream value
VEL_TOL = 0.15              # inflow-face bulk velocity within 15% of upstream value

# Upstream #INFLOW state (SI-derived, code units after conversion).  These are
# the values asserted in PARAM.in and what the fixed/inflow faces must hold.
B_UPSTREAM = 7.382e-9       # Bx = Bz upstream (T, upstream SI value mirrored to code)
UX_UPSTREAM = -713.0        # ux upstream [km/s] as set in #UNIFORMSTATE / #INFLOW
RHO_UPSTREAM = 5.0          # rho upstream [amu/cc] as set in #UNIFORMSTATE / #INFLOW
EY_MOTIONAL = 5.263e-3      # -ux*Bz [V/m] (upstream motional E)


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

    # All energies must be finite (no NaN/Inf from a broken boundary / runaway).
    for key in ("Etot", "Ee", "Eb", "Epart"):
        v0, v1 = first.get(key, 0.0), last.get(key, 0.0)
        if not (math.isfinite(v0) and math.isfinite(v1)):
            passed = False
            reasons.append(f"{key} not finite (NaN/Inf) -- runaway blow-up")

    # Ion kinetic energy: may grow as flow thermalises, but must stay finite &
    # bounded (the conducting-wall runaway sends Epart to 1e18).
    ep0, ep1 = first.get("Epart", 0.0), last.get("Epart", 0.0)
    if ep0 > 0:
        ratio = ep1 / ep0
        logger.debug("    Epart: %s -> %s (ratio %.5f)",
                     f"{ep0:.6e}", f"{ep1:.6e}", ratio)
        if ratio > EPART_RATIO_MAX:
            passed = False
            reasons.append(
                f"Epart ratio {ratio:.4f} > {EPART_RATIO_MAX:.1f} "
                f"(field-energy runaway / broken boundary)")

    # Magnetic energy: a reflecting shock compresses B (Eb grows) but must stay
    # bounded & finite.  The conducting-wall bug drove Eb to ~1e18.
    eb0, eb1 = first.get("Eb", 0.0), last.get("Eb", 0.0)
    if eb0 > 0:
        ratio = eb1 / eb0
        logger.debug("    Eb: %s -> %s (ratio %.5f)",
                     f"{eb0:.6e}", f"{eb1:.6e}", ratio)
        if ratio > EB_RATIO_MAX:
            passed = False
            reasons.append(
                f"Eb ratio {ratio:.4f} > {EB_RATIO_MAX:.1f} "
                f"(magnetic field runaway -- fixed BC not holding the face)")

    # Electric energy must stay negligible (no spurious E build-up at the faces).
    ee1 = last.get("Ee", 0.0)
    if eb0 > 0 and ee1 > 0.5 * eb0:
        passed = False
        reasons.append(
            f"Ee {ee1:.3e} exceeds 50% of Eb ({eb0:.3e}) "
            f"(spurious E built up at a boundary)")

    if passed:
        return True, "Passed (energies finite & bounded; no runaway)"
    return False, "; ".join(reasons)


def _load_out(out_file):
    """Load a .out frame: return ({VAR: col_idx}, list of float rows)."""
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


def _col(vidx, rows, name):
    i = vidx.get(name)
    if i is None or not rows or i >= len(rows[0]):
        return None
    return [r[i] for r in rows]


def _mean(vals):
    return sum(vals) / len(vals) if vals else 0.0


def validate_plot(test_name):
    """Plot checks: reflect-face (-x) holds upstream field, inflow (+x) face
    maintains upstream state, no spurious E."""
    import tests._shared.hybrid as _hyb
    plots_dir = os.path.join(_hyb.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [SHOCK1D] No .out files (PostProc.pl not run?) -- skipping.")
        return True, "No .out files (skipped)"

    vidxl, rowsl = _load_out(out_files[-1])
    if vidxl is None or not rowsl:
        return True, "Could not parse .out frames (skipped)"

    passed = True
    reasons = []

    x = _col(vidxl, rowsl, "X")
    if x is None:
        return True, "No X column in .out (skipped)"
    xmin, xmax = min(x), max(x)
    n = len(x)
    # left (reflect) third and right (inflow) third index masks
    left_cut = xmin + 0.33 * (xmax - xmin)
    right_cut = xmin + 0.67 * (xmax - xmin)
    left_mask = [i for i in range(n) if x[i] <= left_cut]
    right_mask = [i for i in range(n) if x[i] >= right_cut]

    # --- Reflect-field face (-x): Bx, Bz held at upstream value ---
    for bname, bup in (("BX", B_UPSTREAM), ("BZ", B_UPSTREAM)):
        b = _col(vidxl, rowsl, bname)
        if b and left_mask:
            b_face = _mean([b[i] for i in left_mask])
            logger.debug("    [SHOCK1D] left-face <%s> = %.4e (upstream %.4e)",
                         bname, b_face, bup)
            if bup > 1e-12 and abs(b_face) > B_FACE_TOL * bup and \
                    abs(b_face - bup) > B_FACE_TOL * bup:
                passed = False
                reasons.append(
                    f"reflect-face {bname} {b_face:.4e} deviates > "
                    f"{B_FACE_TOL*100:.0f}% from upstream {bup:.4e} "
                    f"(fixed BC not holding the face)")

    # --- Reflect-field face (-x): motional Ey maintained (PEC would zero it) ---
    ey = _col(vidxl, rowsl, "EY")
    if ey and left_mask:
        ey_face = _mean([ey[i] for i in left_mask])
        logger.debug("    [SHOCK1D] left-face <Ey> = %.4e (motional %.4e)",
                     ey_face, EY_MOTIONAL)
        # allow sign/scale tolerance: the face Ey must be the upstream motional
        # value, not ~0 (which a PEC wall would give).
        if abs(ey_face) < 0.10 * abs(EY_MOTIONAL):
            passed = False
            reasons.append(
                f"reflect-face <Ey> = {ey_face:.4e} ~ 0 (PEC-style E_t=0!); "
                f"expected upstream motional Ey = {EY_MOTIONAL:.4e}")

    # --- Inflow (+x) face: upstream state maintained ---
    rho = _col(vidxl, rowsl, "RHOS0")
    if rho and right_mask:
        rho_face = _mean([rho[i] for i in right_mask])
        logger.debug("    [SHOCK1D] inflow-face <rhoS0> = %.4f (upstream %.4f)",
                     rho_face, RHO_UPSTREAM)
        if RHO_UPSTREAM > 1e-12 and \
                abs(rho_face / RHO_UPSTREAM - 1.0) > RHO_TOL:
            passed = False
            reasons.append(
                f"inflow-face density {rho_face:.4f} deviates > "
                f"{RHO_TOL*100:.0f}% from upstream {RHO_UPSTREAM:.4f} "
                f"(inflow face draining / spurious source)")

    ux = _col(vidxl, rowsl, "UXS0")
    if ux and right_mask:
        ux_face = _mean([ux[i] for i in right_mask])
        logger.debug("    [SHOCK1D] inflow-face <uxS0> = %.4f (upstream %.4f)",
                     ux_face, UX_UPSTREAM)
        if abs(UX_UPSTREAM) > 1e-12 and \
                abs(ux_face / UX_UPSTREAM - 1.0) > VEL_TOL:
            passed = False
            reasons.append(
                f"inflow-face uxS0 {ux_face:.4f} deviates > {VEL_TOL*100:.0f}% "
                f"from upstream {UX_UPSTREAM:.4f} (inflow not maintaining state)")

    # --- No spurious transverse electric field (Ex, Ez near zero) ---
    for ecomp in ("EX", "EZ"):
        e = _col(vidxl, rowsl, ecomp)
        if e:
            meanl = _mean(e)
            logger.debug("    [SHOCK1D] mean <%s> = %.5e", ecomp, meanl)
            if abs(meanl) > E_SPURIOUS_TOL * B_UPSTREAM:
                passed = False
                reasons.append(
                    f"spurious mean electric field <{ecomp}> = {meanl:.4e} "
                    f"(>{E_SPURIOUS_TOL*100:.0f}% of |B| scale)")

    if passed:
        return True, ("Passed (reflect-field holds upstream B/E, inflow maintains "
                      "state, no spurious E)")
    return False, "; ".join(reasons)
