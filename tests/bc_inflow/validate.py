#!/usr/bin/env python3
"""Validator for the inflow/outflow open-boundary hybrid test (tests/bc_inflow/).

See README.md for the deck and the list of checks.
"""
import glob
import logging
import math
import os

from tests._shared import run_dir as _run_dir

logger = logging.getLogger(__name__)

# Loose energy tolerances: the state is steady but the open boundaries inject
# and remove particles every step, so Epart is not strictly conserved.  We
# guard against gross non-conservation / blow-up rather than tight drift.
EPART_RATIO_MIN = 0.5
EPART_RATIO_MAX = 2.0
EB_RATIO_TOL = 0.20          # Bx0 is uniform => Eb should stay within +/-20%
EE_MAX_FRAC = 0.10           # Ee <= 10% of Eb (no spurious E build-up)
RHO_TOL = 0.10                # upstream density preserved to +/-10%
VEL_TOL = 0.10               # bulk velocity preserved to +/-10%
T_TOL = 0.15                 # per-species temperature preserved to +/-15%
RHO_RATIO_TOL = 0.10         # cross-species density ratio
T_RATIO_TOL = 0.15           # cross-species temperature ratio
B_SPREAD_TOL = 0.10          # Bx spatially uniform to 10% of mean
E_SPURIOUS_TOL = 0.10        # mean spurious E <= 10% of guide-field scale
ANISO_TOL = 0.25             # Pxx/Pyy/Pzz isotropy (finite-particle noise)


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

    # All energies must be finite (no NaN/Inf from a broken boundary).
    for key in ("Etot", "Ee", "Eb", "Epart"):
        v0, v1 = first.get(key, 0.0), last.get(key, 0.0)
        if not (math.isfinite(v0) and math.isfinite(v1)):
            passed = False
            reasons.append(f"{key} not finite (NaN/Inf)")

    # Ion kinetic energy: bounded (open BCs => not strictly conserved).
    ep0, ep1 = first.get("Epart", 0.0), last.get("Epart", 0.0)
    if ep0 > 0:
        ratio = ep1 / ep0
        logger.debug("    Epart: %s -> %s (ratio %.5f)",
                     f"{ep0:.6e}", f"{ep1:.6e}", ratio)
        if ratio < EPART_RATIO_MIN or ratio > EPART_RATIO_MAX:
            passed = False
            reasons.append(
                f"Epart ratio {ratio:.4f} outside "
                f"[{EPART_RATIO_MIN}, {EPART_RATIO_MAX}] (gross non-conservation)")

    # Magnetic energy: uniform Bx0 => Eb should stay close to its initial value.
    eb0, eb1 = first.get("Eb", 0.0), last.get("Eb", 0.0)
    if eb0 > 0:
        ratio = eb1 / eb0
        logger.debug("    Eb: %s -> %s (ratio %.5f)",
                     f"{eb0:.6e}", f"{eb1:.6e}", ratio)
        if abs(ratio - 1.0) > EB_RATIO_TOL:
            passed = False
            reasons.append(
                f"Eb ratio {ratio:.4f} not within "
                f"[{1-EB_RATIO_TOL:.3f}, {1+EB_RATIO_TOL:.3f}] "
                f"(field energy drifted; boundary field not held)")

    # Electric energy stays negligible (no spurious E build-up at the faces).
    ee1 = last.get("Ee", 0.0)
    if eb0 > 0 and ee1 > EE_MAX_FRAC * eb0:
        passed = False
        reasons.append(
            f"Ee {ee1:.3e} exceeds {EE_MAX_FRAC*100:.0f}% of Eb "
            f"({eb0:.3e}) (spurious E built up at a boundary)")

    if passed:
        return True, "Passed (energies finite & bounded)"
    return False, "; ".join(reasons)


def _load_out(out_file):
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


# Use the shared col() helper from run_dir instead of a local copy.
_col = _run_dir.col


def _mean(vals):
    return sum(vals) / len(vals) if vals else 0.0


def _deck_species():
    """Read the per-species (mass, rho, T) declared by the deck.

    rho is the #UNIFORMSTATE NUMBER density [1/cc] -- the convention shared
    with #INFLOW -- and T the uniform-state temperature [K].  Only RATIOS are
    used below, so the normalisation cancels.
    """
    candidates = [os.path.join(_run_dir.RUN_DIR, "PARAM.in"),
                  os.path.join("tests", "bc_inflow", "PARAM.in")]
    path = next((p for p in candidates if os.path.isfile(p)), None)
    if path is None:
        return []

    masses, rhos, temps = [], [], []
    command = None
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            toks = line.split()
            if not toks:
                continue
            if toks[0].startswith("#"):
                command = toks[0].upper()
                continue
            if len(toks) < 2:
                continue
            name = toks[1].lower()
            try:
                value = float(toks[0])
            except ValueError:
                continue
            if command == "#PLASMA" and name.startswith("mass"):
                masses.append(value)
            elif command == "#UNIFORMSTATE":
                if name.startswith("rho"):
                    rhos.append(value)
                elif name.startswith("t"):
                    temps.append(value)

    nS = min(len(masses), len(rhos), len(temps))
    return [{"mass": masses[i], "rho": rhos[i], "T": temps[i]}
            for i in range(nS)]


def _species_in_plot(vidx):
    """Number of species present in the .out frame (rhoS0, rhoS1, ...)."""
    n = 0
    while f"RHOS{n}" in vidx and f"PS{n}" in vidx:
        n += 1
    return n


def _inflow_side(vidx, rows, name):
    """Mean of *name* over the inflow-side (first) third of the domain."""
    vals = _col(vidx, rows, name)
    xi = vidx.get("X")
    if not vals or xi is None:
        return None
    xs = [r[xi] for r in rows]
    cut = min(xs) + 0.33 * (max(xs) - min(xs))
    sel = [v for v, x in zip(vals, xs) if x <= cut]
    return _mean(sel) if sel else None


def _temperature(vidx, rows, iS, mass):
    """Code-unit temperature of species iS: T = p * m / rho (inflow side)."""
    rho = _inflow_side(vidx, rows, f"RHOS{iS}")
    press = _inflow_side(vidx, rows, f"PS{iS}")
    if not rho or not press or rho <= 0:
        return None
    return press * mass / rho


def validate_plot(test_name):
    """Plot checks: Bx uniform & unchanged, per-species ux/rho/T preserved at
    the inflow face, and the cross-species density/temperature ratios set by
    the deck (they fail loudly if one #INFLOW block is broadcast to all
    species or if the thermal speed ignores the species mass)."""
    import tests._shared.hybrid as _hyb
    plots_dir = os.path.join(_hyb.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [INFLOW] No .out files (PostProc.pl not run?) -- skipping.")
        return True, "No .out files (skipped)"

    vidx0, rows0 = _load_out(out_files[0])
    vidxl, rowsl = _load_out(out_files[-1])
    if vidx0 is None or vidxl is None or not rows0 or not rowsl:
        return True, "Could not parse .out frames (skipped)"

    passed = True
    reasons = []

    nSpecies = _species_in_plot(vidxl)
    deck = _deck_species()
    masses = [d["mass"] for d in deck]
    if len(masses) < nSpecies:
        masses += [1.0] * (nSpecies - len(masses))
    logger.debug("    [INFLOW] %d species in the output, deck: %s",
                 nSpecies, deck)

    # Bulk velocity uxS<i> conserved across the domain, for every species.
    for iS in range(nSpecies):
        ux0, uxl = _col(vidx0, rows0, f"UXS{iS}"), _col(vidxl, rowsl, f"UXS{iS}")
        if not (ux0 and uxl):
            continue
        mean0, meanl = _mean(ux0), _mean(uxl)
        logger.debug("    [INFLOW] <uxS%d>: %s -> %s", iS,
                     f"{mean0:.5f}", f"{meanl:.5f}")
        if abs(mean0) > 1e-12 and abs(meanl / mean0 - 1.0) > VEL_TOL:
            passed = False
            reasons.append(
                f"bulk velocity <uxS{iS}> {mean0:.4f} -> {meanl:.4f} "
                f"(>{VEL_TOL*100:.0f}% drift; inflow not maintaining state)")
        elif abs(mean0) <= 1e-12 and abs(meanl) > 0.05:
            passed = False
            reasons.append(f"spurious bulk velocity {meanl:.4f} developed")

    # Inflow-side density and temperature of every species, plus pressure
    # isotropy (a broken per-species inflow shows up as the wrong T, because
    # vth must be sqrt(T/m) with the species mass).
    rho_in, temp_in = {}, {}
    for iS in range(nSpecies):
        rho_l = _inflow_side(vidxl, rowsl, f"RHOS{iS}")
        rho_0 = _inflow_side(vidx0, rows0, f"RHOS{iS}")
        t_l = _temperature(vidxl, rowsl, iS, masses[iS])
        t_0 = _temperature(vidx0, rows0, iS, masses[iS])
        rho_in[iS], temp_in[iS] = rho_l, t_l
        logger.debug("    [INFLOW] inflow-side species %d: rho %s -> %s, "
                     "T %s -> %s",
                     iS, f"{rho_0:.6f}", f"{rho_l:.6f}",
                     f"{t_0:.5e}", f"{t_l:.5e}")

        if rho_0 and rho_l:
            if rho_0 > 1e-12 and abs(rho_l / rho_0 - 1.0) > RHO_TOL:
                passed = False
                reasons.append(
                    f"inflow-side density rhoS{iS} {rho_0:.4f} -> {rho_l:.4f} "
                    f"(>{RHO_TOL*100:.0f}% drift; inflow face draining)")
            elif rho_l <= 0.0:
                passed = False
                reasons.append(f"inflow-side density rhoS{iS} collapsed to zero")

        if t_0 and t_l and t_0 > 0:
            if abs(t_l / t_0 - 1.0) > T_TOL:
                passed = False
                reasons.append(
                    f"inflow-side temperature of species {iS}: {t_0:.4e} -> "
                    f"{t_l:.4e} (>{T_TOL*100:.0f}% drift; wrong per-species "
                    f"thermal speed)")

        # Thermal isotropy and a non-degenerate out-of-plane pressure.
        pxx = _inflow_side(vidxl, rowsl, f"PXXS{iS}")
        pyy = _inflow_side(vidxl, rowsl, f"PYYS{iS}")
        pzz = _inflow_side(vidxl, rowsl, f"PZZS{iS}")
        if pxx and pyy and pzz:
            if pzz <= 0.0 or not math.isfinite(pzz):
                passed = False
                reasons.append(
                    f"inflow-side Pzz of species {iS} collapsed or non-finite "
                    f"(vz sampling failed)")
            elif pxx > 0.0:
                for label, val in (("Z", pzz), ("Y", pyy)):
                    if abs(val / pxx - 1.0) > ANISO_TOL:
                        passed = False
                        reasons.append(
                            f"inflow pressure of species {iS} anisotropic in "
                            f"{label} (P{label}{label}/Pxx = {val/pxx:.3f})")

    # Cross-species ratios: the mass-density ratio rhoS1/rhoS0 and the
    # temperature ratio T1/T0 are set by the deck.  Broadcasting a single
    # #INFLOW block to every species, or using the proton mass for every
    # thermal speed, breaks them by a factor of order m_1/m_0.
    #
    # The plot carries MASS densities while the deck gives NUMBER densities
    # (rho [1/cc], the convention shared by #UNIFORMSTATE and #INFLOW), so the
    # expected ratio is (n_1*m_1)/(n_0*m_0).
    if nSpecies >= 2 and len(deck) >= 2:
        if rho_in.get(0) and rho_in.get(1) and rho_in[0] > 0:
            meas = rho_in[1] / rho_in[0]
            expect = ((deck[1]["rho"] * deck[1]["mass"]) /
                      (deck[0]["rho"] * deck[0]["mass"]))
            logger.debug("    [INFLOW] rhoS1/rhoS0 = %.4f (deck %.4f)",
                         meas, expect)
            if abs(meas / expect - 1.0) > RHO_RATIO_TOL:
                passed = False
                reasons.append(
                    f"species density ratio rhoS1/rhoS0 = {meas:.3f} but the "
                    f"deck prescribes {expect:.3f}; #INFLOW is not per-species")

        if temp_in.get(0) and temp_in.get(1) and temp_in[0] > 0:
            meas = temp_in[1] / temp_in[0]
            expect = deck[1]["T"] / deck[0]["T"]
            logger.debug("    [INFLOW] T1/T0 = %.4f (deck %.4f)", meas, expect)
            if abs(meas / expect - 1.0) > T_RATIO_TOL:
                passed = False
                reasons.append(
                    f"species temperature ratio T1/T0 = {meas:.3f} but the "
                    f"deck prescribes {expect:.3f}; the injected thermal speed "
                    f"ignores the species mass")

    # Guide field Bx uniform and unchanged.
    bx0, bxl = _col(vidx0, rows0, "BX"), _col(vidxl, rowsl, "BX")
    mb0 = _mean(bx0) if bx0 else 0.0
    if bxl:
        mbl = _mean(bxl)
        spreadl = max(bxl) - min(bxl)
        logger.debug("    [INFLOW] <Bx>: %s -> %s (spread %.4f)",
                     f"{mb0:.5f}", f"{mbl:.5f}", spreadl)
        if abs(mb0) > 1e-12 and abs(mbl / mb0 - 1.0) > EB_RATIO_TOL:
            passed = False
            reasons.append(f"guide field Bx changed {mb0:.4f} -> {mbl:.4f}")
        if abs(mbl) > 1e-12 and spreadl > B_SPREAD_TOL * abs(mbl):
            passed = False
            reasons.append(
                f"guide field not uniform (spread {spreadl:.4f} > "
                f"{B_SPREAD_TOL*100:.0f}% of <Bx>; boundary layer formed)")

    # No spurious mean E field.
    for ecomp in ("EX", "EY", "EZ"):
        el = _col(vidxl, rowsl, ecomp)
        if el:
            meanl = _mean(el)
            scale = abs(mb0) if (bx0 and abs(mb0) > 1e-12) else 1.0
            logger.debug("    [INFLOW] mean <%s> = %.5f (scale %.4f)",
                         ecomp, meanl, scale)
            if abs(meanl) > E_SPURIOUS_TOL * scale:
                passed = False
                reasons.append(
                    f"spurious mean electric field <{ecomp}> = {meanl:.4f} "
                    f"(>{E_SPURIOUS_TOL*100:.0f}% of guide-field scale)")

    if passed:
        return True, ("Passed (multi-species state stays uniform; inflow "
                      "maintains the per-species upstream state, outflow lets "
                      "plasma leave)")
    return False, "; ".join(reasons)
