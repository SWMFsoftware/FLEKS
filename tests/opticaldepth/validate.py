#!/usr/bin/env python3
"""Validator for the optical-depth attenuation test (tests/opticaldepth).

The deck produces exospheric ions from two neutral components (H and O) whose
EUV attenuation is

    S_i(r) = n_i(r) * nu0_i * A(r)

with the BATSRUS optical-depth model

    tau(r) = sum_c n_c(r) sigma_c H_c(r),   mu = max(x.solarDir/r, cosSzaFloor)
    A(r)   = exp(-tau/mu)                  (plain), or
    tau    = n_c* (r) sigma_c* H_c*        (Chapman: only one component)
    A(r)   = exp(-tau chap(r/H_c*, cosSZA))

This module rebuilds that model from PARAM.in and compares it with the ion
density produced by the run.  The source contribution is measured as the
difference between the final and the initial frame, which cancels the uniform
background plasma and needs no knowledge of the plot density unit.

Checks performed (all on the final frame):
  1. the O+ source profile along the subsolar line follows exp(-tau/mu),
  2. the H+ source profile does too, even though H has sigma = 0: this proves
     that tau sums over *all* neutral components and not just the parent one,
  3. the nightside production is suppressed with respect to the dayside,
  4. the frame contains only finite values.
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"

# Particle-tracking tolerance passed to validate_test_particles() by the runner.
# This test does not enable #PARTICLETRACKER.
PARTICLE_TOL = None

# Radial samples along the subsolar line, in units of rPlanet.
SAMPLE_RADII = (1.2, 1.5, 1.8, 2.2, 2.6)
REFERENCE_RADIUS = SAMPLE_RADII[0]

# |y| accepted for a row to count as "on the subsolar line", in units of rPlanet.
Y_TOL = 0.15

# Relative agreement required between the measured and modelled profile.
PROFILE_TOL = 0.25


def set_run_dir(run_dir):
    """Point the validator at the current run directory."""
    global RUN_DIR
    RUN_DIR = run_dir


# ---------------------------------------------------------------------------
# PARAM.in parsing
# ---------------------------------------------------------------------------
def _num(token):
    """Parse a float, returning None for free-text lines in the deck."""
    try:
        return float(token)
    except ValueError:
        return None


def _parse_deck():
    """Return the model parameters read from the deck used by the run."""
    param_path = os.path.join(RUN_DIR, "PARAM.in")
    deck = {
        "rPlanet": 3.39e6,
        "lNormSI": 1000.0,
        "xMin": -1.017e7,
        "xMax": 1.017e7,
        "typeProfile": "Exponential",
        "n0": [],
        "H0": [],
        "nu0": [],
        "crossSection": [],
        "scaleHeight": [],
        "solarDir": [1.0, 0.0, 0.0],
        "chapmanFunction": False,
        "chapmanComponent": -1,
        "minProduction": None,
        "tauFloor": None,
        "cosSzaFloor": None,
        "chapmanTauMax": None,
    }
    section = None
    exo_counter = 0
    n_exo = 0
    opt_index = 0
    norm_counter = 0
    geom_counter = 0
    with open(param_path, "r", encoding="latin-1") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                section = s.split()[0].upper()
                if section == "#EXOSPHERE":
                    exo_counter = -1
                elif section == "#OPTICALDEPTH":
                    opt_index = 0
                elif section == "#NORMALIZATION":
                    norm_counter = 0
                elif section == "#GEOMETRY":
                    geom_counter = 0
                continue
            parts = s.split()
            val = _num(parts[0])
            if section == "#EXOSPHERE" and exo_counter == -1:
                deck["typeProfile"] = parts[0]
                exo_counter += 1
                continue
            if section == "#OPTICALDEPTH" and val is None:
                # The T/F flag carries a non-numeric value.
                if "chapmanFunction" in s:
                    deck["chapmanFunction"] = parts[0].upper() == "T"
                    opt_index += 1
                continue
            if val is None:
                continue
            if section == "#BODYSIZE":
                deck["rPlanet"] = val
            elif section == "#NORMALIZATION":
                # First entry is lNormSI [m], the length of one code unit.
                if norm_counter == 0:
                    deck["lNormSI"] = val
                norm_counter += 1
            elif section == "#GEOMETRY":
                # First two entries are xMin and xMax [m].
                if geom_counter == 0:
                    deck["xMin"] = val
                elif geom_counter == 1:
                    deck["xMax"] = val
                geom_counter += 1
            elif section == "#EXOSPHERE":
                if exo_counter == 0:
                    n_exo = int(val)
                elif exo_counter > 0:
                    which = (exo_counter - 1) % 4
                    if which == 0:
                        deck["n0"].append(val)
                    elif which == 1:
                        deck["H0"].append(val)
                exo_counter += 1
            elif section == "#PHOTOIONIZATION":
                deck["nu0"].append(val)
            elif section == "#OPTICALDEPTH":
                if opt_index < n_exo:
                    deck["crossSection"].append(val)
                elif opt_index < 2 * n_exo:
                    deck["scaleHeight"].append(val)
                elif opt_index < 2 * n_exo + 3:
                    deck["solarDir"][opt_index - 2 * n_exo] = val
                elif "chapmanFunction" in s:
                    deck["chapmanFunction"] = parts[0].upper() == "T"
                elif "minProduction" in s:
                    deck["minProduction"] = val
                elif "tauFloor" in s:
                    deck["tauFloor"] = val
                elif "cosSzaFloor" in s:
                    deck["cosSzaFloor"] = val
                elif "chapmanTauMax" in s:
                    deck["chapmanTauMax"] = val
                elif "chapmanComponent" in s:
                    deck["chapmanComponent"] = int(val)
                opt_index += 1

    # BATSRUS defaults for anything the deck left out.
    if deck["minProduction"] is None:
        deck["minProduction"] = 1.0e-6 if deck["chapmanFunction"] else 1.0e-5
    if deck["tauFloor"] is None:
        deck["tauFloor"] = 6.0e-3
    if deck["cosSzaFloor"] is None:
        deck["cosSzaFloor"] = 5.0e-4
    if deck["chapmanTauMax"] is None:
        deck["chapmanTauMax"] = 13.8
    if deck["chapmanComponent"] < 0 and deck["chapmanFunction"]:
        # Component with the largest vertical optical depth at the reference
        # radius, mirroring the automatic selection in ExoSource.
        best, best_tau = 0, -1.0
        for iC in range(len(deck["crossSection"])):
            tau = (_neutral(deck, iC, deck["rPlanet"]) *
                   deck["crossSection"][iC] * _scale_height(deck, iC,
                                                            deck["rPlanet"]))
            if tau > best_tau:
                best, best_tau = iC, tau
        deck["chapmanComponent"] = best
    return deck


def _scale_height(deck, iC, r):
    """Local scale height [m] of component iC at radius r."""
    if 0 <= iC < len(deck["scaleHeight"]) and deck["scaleHeight"][iC] > 0.0:
        return deck["scaleHeight"][iC]
    h0 = deck["H0"][iC]
    if h0 <= 0.0:
        return 0.0
    if deck["typeProfile"] == "Exponential":
        return h0
    if deck["typeProfile"] == "Chamberlain":
        return r * r / h0
    if deck["typeProfile"] == "Power-Law":
        k0 = deck.get("k0", [0.0] * len(deck["H0"]))[iC]
        return r / k0 if k0 > 0.0 else 0.0
    return 0.0


def _neutral(deck, iC, r):
    """Neutral number density [m^-3] of component iC at radius r."""
    if r < deck["rPlanet"]:
        return 0.0
    h0 = deck["H0"][iC]
    if h0 <= 0.0:
        return 0.0
    if deck["typeProfile"] == "Exponential":
        return deck["n0"][iC] * math.exp(-(r - deck["rPlanet"]) / h0)
    if deck["typeProfile"] == "Chamberlain":
        return deck["n0"][iC] * math.exp(-h0 * (1.0 / deck["rPlanet"] - 1.0 / r))
    if deck["typeProfile"] == "Power-Law":
        k0 = deck.get("k0", [0.0] * len(deck["H0"]))[iC]
        return deck["n0"][iC] * (deck["rPlanet"] / r) ** k0
    return 0.0


def _chapman(Xp, cosSZA):
    """Chapman function (Smith & Smith 1972), as in BATSRUS ModUserMars."""
    y = math.sqrt(0.5 * Xp) * abs(cosSZA)
    if cosSZA > 0.0:
        if y < 8.0:
            return (math.sqrt(0.5 * math.pi * Xp) *
                    (1.0606963 + 0.5564383 * y) /
                    (1.0619896 + 1.7245609 * y + y * y))
        if y < 100.0:
            return math.sqrt(0.5 * math.pi * Xp) * 0.56498823 / (0.6651874 + y)
        return 0.0
    sin_sza = math.sqrt(max(1.0 - cosSZA * cosSZA, 0.0))
    if y < 8.0:
        return (math.sqrt(2.0 * math.pi * Xp) *
                (math.sqrt(sin_sza) * math.exp(Xp * (1.0 - sin_sza)) -
                 0.5 * (1.0606963 + 0.5564383 * y) /
                 (1.0619896 + 1.7245609 * y + y * y)))
    if y < 100.0:
        return (math.sqrt(2.0 * math.pi * Xp) *
                (math.sqrt(sin_sza) *
                 math.exp(min(100.0, Xp * (1.0 - sin_sza))) -
                 0.5 * 0.56498823 / (0.6651874 + y)))
    return 0.0


def _tau_vertical(deck, r, only_comp):
    total = 0.0
    for iC in range(len(deck["crossSection"])):
        if only_comp >= 0 and iC != only_comp:
            continue
        sigma = deck["crossSection"][iC]
        h = _scale_height(deck, iC, r)
        if sigma <= 0.0 or h <= 0.0:
            continue
        total += _neutral(deck, iC, r) * sigma * h
    return total


def _attenuation(deck, x, y, r):
    """Model attenuation A(r) at a point (x, y) [m]."""
    r_safe = max(r, 1.0e-3)
    proj = (x * deck["solarDir"][0] + y * deck["solarDir"][1] +
            deck["solarDir"][2] * 0.0)
    cos_sza = proj / r_safe
    mu = max(cos_sza, deck["cosSzaFloor"])
    if deck["chapmanFunction"]:
        iC = deck["chapmanComponent"]
        tau_v = _tau_vertical(deck, r, iC)
        if tau_v > deck["chapmanTauMax"]:
            return deck["minProduction"] * max(mu, deck["minProduction"])
        h = _scale_height(deck, iC, r)
        if h <= 0.0:
            return deck["minProduction"]
        chap = _chapman(r / h, cos_sza)
        if chap <= 0.0:
            # Negative fit value deep on the nightside: optically thick.
            return deck["minProduction"] * max(mu, deck["minProduction"])
        return max(math.exp(-tau_v * chap), deck["minProduction"])
    tau = max(_tau_vertical(deck, r, -1), deck["tauFloor"]) / mu
    if tau >= 1.0e30 or proj <= 0.0:
        return deck["minProduction"]
    return max(math.exp(-tau), deck["minProduction"])


def _production(deck, iC, x, y, r):
    """Model source strength n_i(r) * nu0_i * A(r) at (x, y) [m], r = |x, y|."""
    if iC >= len(deck["nu0"]):
        return 0.0
    return _neutral(deck, iC, r) * deck["nu0"][iC] * _attenuation(
        deck, x, y, r)


# ---------------------------------------------------------------------------
# Plot reading
# ---------------------------------------------------------------------------
def _read_out(out_file):
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


def _plot_scale(rows, vidx, deck):
    """Return the length of one plot unit in metres.

    The plot coordinates span exactly the #GEOMETRY box, so the box length
    from the deck and the coordinate span in the data give the unit directly.
    No assumption about the output normalisation is needed.
    """
    xs = [r[0] for r in rows]
    span_plot = max(xs) - min(xs)
    span_si = deck["xMax"] - deck["xMin"]
    if span_plot > 0.0 and span_si > 0.0:
        return span_si / span_plot
    return deck["lNormSI"]


def _sample(rows, vidx, deck, target_rp, nightside=False):
    """Pick the row closest to (x = +/-target_rp, y = 0) in rPlanet units."""
    scale = _plot_scale(rows, vidx, deck)
    rp = deck["rPlanet"]
    best = None
    for row in rows:
        x_si = row[0] * scale
        y_si = row[1] * scale
        if (x_si > 0.0) == nightside or x_si == 0.0:
            continue
        if abs(y_si) / rp > Y_TOL:
            continue
        d = abs(abs(x_si) / rp - target_rp)
        if best is None or d < best[0]:
            best = (d, row, x_si, y_si)
    if best is None:
        return None
    _, row, x_si, y_si = best
    return row, x_si, y_si, math.hypot(x_si, y_si)


# ---------------------------------------------------------------------------
# Validators
# ---------------------------------------------------------------------------
def validate_log(pic_diags=None, test_name=None):
    """Require a completed run with a usable PIC energy log."""
    if not pic_diags or len(pic_diags) < 2:
        return False, "No PIC energy log produced"
    last = pic_diags[-1]
    for key in ("Epart", "Eb", "Ee", "Etot"):
        v = last.get(key)
        if v is not None and not math.isfinite(v):
            return False, f"Non-finite {key} in the PIC log: {v}"
    logger.debug("    PIC log frames: %d, final Etot=%s",
                 len(pic_diags), last.get("Etot"))
    return True, "Passed"


def validate_plot(test_name=None):
    """Compare the produced ion profile with the analytic optical-depth model."""
    deck = _parse_deck()
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if len(out_files) < 2:
        return False, f"Expected at least 2 plot frames, found {len(out_files)}"

    vidx0, rows0 = _read_out(out_files[0])
    vidx1, rows1 = _read_out(out_files[-1])
    if not vidx0 or not vidx1:
        return False, "Could not parse the plot frames"

    # 1. the frame must contain only finite values
    for row in rows1:
        for v in row:
            if not math.isfinite(v):
                return False, "Non-finite value in the final plot frame"

    model = "Chapman" if deck["chapmanFunction"] else "exp(-tau/mu)"
    logger.debug("    [OD] model=%s rPlanet=%.3e type=%s tau(rP)=%.4g "
                 "plot unit=%.1f m",
                 model, deck["rPlanet"], deck["typeProfile"],
                 _tau_vertical(deck, deck["rPlanet"], -1),
                 _plot_scale(rows1, vidx1, deck))

    # ExoSource maps exosphere component iC to plasma species iC + 1, so
    # component 1 (O) feeds rhoS2 and component 0 (H) feeds rhoS1.
    # 2. the radial source profile follows the attenuation model
    for iC, sp in ((1, "O+"), (0, "H+")):
        key = f"RHOS{iC + 1}"
        s0 = _sample(rows0, vidx0, deck, REFERENCE_RADIUS)
        s1 = _sample(rows1, vidx1, deck, REFERENCE_RADIUS)
        if s0 is None or s1 is None:
            return False, f"No sample row at r = {REFERENCE_RADIUS} rPlanet"
        row0, row1 = s0[0], s1[0]
        ref_measured = row1[vidx1[key]] - row0[vidx0[key]]
        if ref_measured <= 0.0:
            return False, (f"{sp}: no source signal at r = "
                           f"{REFERENCE_RADIUS} rPlanet "
                           f"(delta = {ref_measured:.3e})")
        ref_model = _production(deck, iC, s1[1], s1[2], s1[3])

        for r_units in SAMPLE_RADII:
            a = _sample(rows0, vidx0, deck, r_units)
            b = _sample(rows1, vidx1, deck, r_units)
            if a is None or b is None:
                continue
            measured = b[0][vidx1[key]] - a[0][vidx0[key]]
            if ref_measured <= 0.0:
                continue
            m_ratio = measured / ref_measured
            p_ratio = _production(deck, iC, b[1], b[2], b[3]) / ref_model
            logger.debug(
                "    [OD] %s r=%.2f rP: measured=%.4e model=%.4e "
                "ratio meas/model=%.3f",
                sp, r_units, measured, p_ratio,
                m_ratio / p_ratio if p_ratio else float("nan"))

            if p_ratio <= 0.0:
                continue
            err = abs(m_ratio / p_ratio - 1.0)
            if err > PROFILE_TOL:
                return False, (
                    f"{sp} radial profile at r = {r_units} rPlanet does not "
                    f"follow the {model} model: measured/model "
                    f"= {m_ratio / p_ratio:.3f} (tolerance "
                    f"{PROFILE_TOL:.0%})")

    # 3. nightside suppression
    day = _sample(rows1, vidx1, deck, 1.5)
    night = _sample(rows1, vidx1, deck, 1.5, nightside=True)
    if day is not None and night is not None:
        key = "RHOS2"
        day_v = day[0][vidx1[key]]
        night_v = night[0][vidx1[key]]
        logger.debug("    [OD] nightside rhoS2=%.4e vs dayside=%.4e",
                     night_v, day_v)
        if day_v > 0.0 and night_v > 0.1 * day_v:
            return False, (f"nightside production not suppressed: "
                           f"rhoS2(-1.5 rP) = {night_v:.3e} vs "
                           f"rhoS2(+1.5 rP) = {day_v:.3e}")

    return True, "Passed"
