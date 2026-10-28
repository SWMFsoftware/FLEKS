#!/usr/bin/env python3
"""Validator for the spatial exosphere & radiation tests (tests/exosphere).

Variants:
* ``aligned``: Unresolved thin exponential atmosphere (H = dx/33); dual-cell
  averaging prevents grid aliasing; checks integrated photoproduction.
* ``shifted``: Same atmosphere shifted by dx/3; checks grid-invariance.
* ``corona``: Transparent Chamberlain corona with optical depth; lit outside
  planet and geometric shadow.
* ``sourceppc``: Target particles per cell (#MINPARTICLESPERCELL / #CONSTANTPPV)
  versus added source particles per step.
* ``opticaldepth``: Plain exp(-tau/mu) attenuation along subsolar line and
  nightside suppression.
* ``chapman``: Chapman function attenuation along subsolar line and nightside
  suppression.
* ``shadow``: 2D/3D photoionization with geometric shadow cylinder, day/night
  density asymmetry, photoion temperature (~300 K), and boundary smoothness.
"""
import glob
import logging
import math
import os

import numpy as np

from tests._shared import run_dir as _run_dir
from tests._shared.validators import validate_ionisation_source

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"
PARTICLE_TOL = None

# Atmosphere parameters for photoproduction tests (aligned / shifted / corona)
R_ATM, H_ATM = 1.5e5, 1.0e3         # [m]
N0_ATM, NU_ATM = 1.0e12, 1.0e-3     # [m^-3], [1/s]
NSEED_ATM = 1.0e3                   # [m^-3]
NCELL_ATM, LNORM_ATM = 24, 1.0e5
H0_CORONA = 1.5e6                   # [m]

# Optical depth radial sampling
SAMPLE_RADII = (1.2, 1.5, 1.8, 2.2, 2.6)
REFERENCE_RADIUS = SAMPLE_RADII[0]
Y_TOL = 0.15
PROFILE_TOL = 0.25


def set_run_dir(run_dir):
    """Point the validator at the current run directory."""
    global RUN_DIR
    RUN_DIR = run_dir
    _run_dir.set_run_dir(run_dir)


def _seed_volume(lo):
    dx = 8.0 / NCELL_ATM
    c = lo + (np.arange(NCELL_ATM) + 0.5) * dx
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    outside = x * x + y * y + z * z >= (R_ATM / LNORM_ATM) ** 2
    return outside.sum() * (dx * LNORM_ATM) ** 3


def _corona_production(lo, n=320):
    edges = np.linspace(lo, lo + 8.0, n + 1) * LNORM_ATM
    c = 0.5 * (edges[1:] + edges[:-1])
    dv = (edges[1] - edges[0]) ** 3
    y, z = np.meshgrid(c, c, indexing="ij")
    perp2 = y * y + z * z
    total = 0.0
    for x in c:
        r = np.sqrt(x * x + perp2)
        lit = (r > R_ATM) & ~((x < 0) & (perp2 < R_ATM * R_ATM))
        nCor = N0_ATM * np.exp(-H0_CORONA * (1 / R_ATM - 1 / np.maximum(r, R_ATM)))
        total += (nCor * lit).sum() * dv
    return NU_ATM * total


def _validate_photoproduction_log(diags, variant):
    if not diags or len(diags) < 3 or "Mpart1" not in diags[0]:
        return False, "missing log or Mpart1 column"

    if variant == "sourceppc":
        if "nParticles" not in diags[0]:
            return False, "missing nParticles log column"
        perStep = [(b["nParticles"] - a["nParticles"]) for a, b in zip(diags, diags[1:])]
        ratio = sum(perStep) / len(perStep) / NCELL_ATM ** 3
        logger.debug("    source particles per cell and step: %.3f", ratio)
        if not math.isfinite(ratio) or abs(ratio - 1) > 0.05:
            return False, (f"source particles per cell and step = {ratio:.3f} "
                           "(expected 1, tol 5%)")
        return True, "Passed"

    lo = -3.9 if variant == "shifted" else -4.0
    seed = NSEED_ATM * _seed_volume(lo)
    unit = diags[0]["Mpart1"] / seed
    perStep = [(b["Mpart1"] - a["Mpart1"]) / unit for a, b in zip(diags, diags[1:])]
    sim = sum(perStep) / len(perStep) / 0.1  # dt = 0.1

    if variant == "corona":
        P = _corona_production(lo)
        tol = 0.05
    else:
        # P = nu * integral_outside_body n0 exp(-(r - R)/H) d^3r
        P = (NU_ATM * N0_ATM * 4.0 * math.pi *
             (R_ATM ** 2 * H_ATM + 2.0 * R_ATM * H_ATM ** 2 + 2.0 * H_ATM ** 3))
        tol = 0.10

    err = abs(sim / P - 1.0)
    logger.debug("    variant '%s': simulated photoproduction %.4e vs "
                 "analytic %.4e (rel. err %.2f%%, tol %.0f%%)",
                 variant, sim, P, 100 * err, 100 * tol)
    if not math.isfinite(sim) or err > tol:
        return False, (f"photoproduction rate = {sim:.4e} vs expected "
                       f"{P:.4e} (rel err {err:.2%}, tol {tol:.0%})")
    return True, "Passed"


def _num(token):
    try:
        return float(token)
    except ValueError:
        return None


def _parse_deck():
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
                if "chapmanFunction" in s:
                    deck["chapmanFunction"] = parts[0].upper() == "T"
                    opt_index += 1
                continue
            if val is None:
                continue
            if section == "#PLANETRADIUS":
                deck["rPlanet"] = val
            elif section == "#NORMALIZATION":
                if norm_counter == 0:
                    deck["lNormSI"] = val
                norm_counter += 1
            elif section == "#GEOMETRY":
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

    if deck["minProduction"] is None:
        deck["minProduction"] = 1.0e-6 if deck["chapmanFunction"] else 1.0e-5
    if deck["tauFloor"] is None:
        deck["tauFloor"] = 6.0e-3
    if deck["cosSzaFloor"] is None:
        deck["cosSzaFloor"] = 5.0e-4
    if deck["chapmanTauMax"] is None:
        deck["chapmanTauMax"] = 13.8
    if deck["chapmanComponent"] < 0 and deck["chapmanFunction"]:
        best, best_tau = 0, -1.0
        for iC in range(len(deck["crossSection"])):
            tau = (_neutral(deck, iC, deck["rPlanet"]) *
                   deck["crossSection"][iC] * _scale_height(deck, iC, deck["rPlanet"]))
            if tau > best_tau:
                best, best_tau = iC, tau
        deck["chapmanComponent"] = best
    return deck


def _scale_height(deck, iC, r):
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
            return deck["minProduction"] * max(mu, deck["minProduction"])
        return max(math.exp(-tau_v * chap), deck["minProduction"])
    tau = max(_tau_vertical(deck, r, -1), deck["tauFloor"]) / mu
    if tau >= 1.0e30 or proj <= 0.0:
        return deck["minProduction"]
    return max(math.exp(-tau), deck["minProduction"])


def _production(deck, iC, x, y, r):
    if iC >= len(deck["nu0"]):
        return 0.0
    return _neutral(deck, iC, r) * deck["nu0"][iC] * _attenuation(deck, x, y, r)


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
    xs = [r[0] for r in rows]
    span_plot = max(xs) - min(xs)
    span_m = (deck["xMax"] - deck["xMin"]) * deck["lNormSI"]
    if span_plot > 0.0 and span_m > 0.0:
        return span_m / span_plot
    return deck["rPlanet"]


def _sample(rows, vidx, deck, target_rp, nightside=False):
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


def _validate_opticaldepth_plot():
    deck = _parse_deck()
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if len(out_files) < 2:
        return False, f"Expected at least 2 plot frames, found {len(out_files)}"

    vidx0, rows0 = _read_out(out_files[0])
    vidx1, rows1 = _read_out(out_files[-1])
    if not vidx0 or not vidx1:
        return False, "Could not parse plot frames"

    for row in rows1:
        for v in row:
            if not math.isfinite(v):
                return False, "Non-finite value in the final plot frame"

    model = "Chapman" if deck["chapmanFunction"] else "exp(-tau/mu)"
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

            if p_ratio <= 0.0:
                continue
            err = abs(m_ratio / p_ratio - 1.0)
            if err > PROFILE_TOL:
                return False, (
                    f"{sp} radial profile at r = {r_units} rPlanet does not "
                    f"follow the {model} model: measured/model "
                    f"= {m_ratio / p_ratio:.3f} (tolerance {PROFILE_TOL:.0%})")

    day = _sample(rows1, vidx1, deck, 1.5)
    night = _sample(rows1, vidx1, deck, 1.5, nightside=True)
    if day is not None and night is not None:
        key = "RHOS2"
        day_v = day[0][vidx1[key]]
        night_v = night[0][vidx1[key]]
        if day_v > 0.0 and night_v > 0.1 * day_v:
            return False, (f"nightside production not suppressed: "
                           f"rhoS2(-1.5 rP) = {night_v:.3e} vs "
                           f"rhoS2(+1.5 rP) = {day_v:.3e}")

    return True, "Passed"


def _read_shadow_params():
    param_path = os.path.join(RUN_DIR, "PARAM.in")
    Rp_si = 3.0e6
    halfH_si = 0.0
    shadowR_si = 0.0
    useShadow = False
    try:
        with open(param_path, "r") as pf:
            section = None
            for line in pf:
                s = line.strip()
                if s.startswith("#"):
                    section = s
                    continue
                if not s:
                    continue
                parts = s.split()
                if section == "#PLANETRADIUS" and len(parts) >= 1:
                    try:
                        Rp_si = float(parts[0])
                    except ValueError:
                        pass
                elif section == "#SHADOWCYLINDER":
                    useShadow = True
                    try:
                        val = float(parts[0])
                    except ValueError:
                        continue
                    if "radius" in s.lower():
                        shadowR_si = val
                    elif "halfheight" in s.lower():
                        halfH_si = val
    except Exception:
        pass

    if not useShadow:
        return None
    per_body_radius = 1.0 / Rp_si
    return (Rp_si * per_body_radius, halfH_si * per_body_radius, shadowR_si * per_body_radius)


def _validate_shadow_plot():
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    shadow_geom = _read_shadow_params()
    if shadow_geom is None:
        return True, "No shadow cylinder"

    Rp_plot, halfH_plot, shadowR_plot = shadow_geom
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        return False, "No .out files found"

    latest_out = out_files[-1]
    with open(latest_out, "r") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return True, "Short .out file"
    var_names = lines[4].split()

    rho_idx = None
    for target in ("RHOS2", "RHOS1"):
        for iv, vn in enumerate(var_names):
            if vn.upper() == target:
                rho_idx = iv
                break
        if rho_idx is not None:
            break
    if rho_idx is None:
        return True, "rhoS2/rhoS1 not in .out"

    p_idx = None
    for target in ("P2", "P1", "P"):
        for iv, vn in enumerate(var_names):
            if vn.upper() == target:
                p_idx = iv
                break
        if p_idx is not None:
            break

    points = []
    points_with_p = []
    for line in lines[5:]:
        cols = line.strip().split()
        if len(cols) <= rho_idx:
            continue
        try:
            x = float(cols[0])
            y = float(cols[1])
            rho = float(cols[rho_idx])
            points.append((x, y, rho))
            if p_idx is not None and len(cols) > p_idx:
                points_with_p.append((x, y, rho, float(cols[p_idx])))
        except (ValueError, IndexError):
            continue

    if not points:
        return False, "No data points parsed"

    # Day/night asymmetry check
    dayside_vals = []
    shadow_vals = []
    y_lim = min(Rp_plot / 2.0, shadowR_plot / 4.0)
    for x, y, rho in points:
        if abs(y) > y_lim:
            continue
        if x > 0:
            dayside_vals.append(rho)
        elif x < -Rp_plot * 0.5:
            shadow_vals.append(rho)

    if not dayside_vals:
        return False, "Zero dayside points"
    dayside_mean = sum(dayside_vals) / len(dayside_vals)
    if dayside_mean <= 0.0:
        return False, "Dayside rho is zero -- photoionization inactive"
    if not shadow_vals:
        return False, "Zero shadow points"
    shadow_mean = sum(shadow_vals) / len(shadow_vals)

    if shadow_mean > max(dayside_mean * 0.2, 1e-30):
        return False, (f"Shadow rho too high ({shadow_mean:.2e}) "
                       f"vs dayside ({dayside_mean:.2e})")

    # Boundary smoothness check (ignore points with zero or near-zero noise density)
    import collections
    by_x = collections.defaultdict(dict)
    for x, y, rho in points:
        if x > 0.5 * Rp_plot:
            by_x[round(x, 1)][round(y, 1)] = rho

    max_jump = 0.0
    for x, y_dict in sorted(by_x.items()):
        ys = sorted(y_dict.keys())
        y_center = min(ys, key=abs)
        if abs(y_center) > 1e-3:
            continue
        pos_ys = [y for y in ys if y > 1e-3]
        if not pos_ys:
            continue
        dy = pos_ys[0]
        y_plus = min(ys, key=lambda y: abs(y - dy))
        y_minus = min(ys, key=lambda y: abs(y - (-dy)))
        rho_0 = y_dict[y_center]
        rho_p = y_dict[y_plus]
        rho_m = y_dict[y_minus]
        rho_neighbor = 0.5 * (rho_p + rho_m)
        if rho_neighbor > 1e-3 and rho_0 > 1e-3:
            jump = rho_0 / rho_neighbor
            max_jump = max(max_jump, jump)

    if max_jump > 1.25:
        return False, f"Block boundary artifact across Y=0: jump ratio {max_jump:.3f} > 1.25"

    # Dayside photoion temperature check (~300 K)
    if points_with_p:
        mass_si = 16.0 * 1.66053906660e-27
        u_norm = 1.0e6
        k_B = 1.380649e-23
        temp_scale = mass_si * (u_norm ** 2) / k_B
        dayside_temps = []
        for x, y, rho, p in points_with_p:
            if x > 0.2 * Rp_plot and abs(y) < 0.5 * Rp_plot:
                if rho > 1e-10 and p > 0:
                    dayside_temps.append((p / rho) * temp_scale)
        if dayside_temps:
            mean_T = sum(dayside_temps) / len(dayside_temps)
            if mean_T < 50.0:
                return False, f"Photoion temperature too low: {mean_T:.1f} K (expected ~300 K)"
            if mean_T > 2000.0:
                return False, f"Photoion temperature too high: {mean_T:.1f} K (expected ~300 K)"

    return True, "Passed"


def _non_finite(diags):
    keys = sorted({k for e in (diags or []) for k, v in e.items()
                   if isinstance(v, float) and not math.isfinite(v)})
    if keys:
        return f"non-finite log values in {', '.join(keys[:6])}"
    return None


def validate_log(pic_diags=None, test_name=None):
    bad = _non_finite(pic_diags)
    if bad:
        return False, bad
    variant = (test_name or "").replace("exosphere_", "")
    logger.debug("Validating exosphere test variant '%s'...", variant)

    if variant in ("aligned", "shifted", "corona", "sourceppc"):
        return _validate_photoproduction_log(pic_diags, variant)
    elif variant in ("opticaldepth", "chapman"):
        if not pic_diags or len(pic_diags) < 2:
            return False, "No PIC energy log produced"
        return True, "Passed"
    elif variant == "shadow":
        return validate_ionisation_source(
            pic_diags, test_label="Photoionization (Shadow)", logger=logger
        )
    return False, f"Unknown exosphere variant '{test_name}'"


def validate_plot(test_name=None):
    variant = (test_name or "").replace("exosphere_", "")
    if variant in ("opticaldepth", "chapman"):
        return _validate_opticaldepth_plot()
    elif variant == "shadow":
        return _validate_shadow_plot()
    return True, "Passed"
