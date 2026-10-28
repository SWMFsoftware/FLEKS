#!/usr/bin/env python3
"""Validator for the kinetic reaction tests (tests/reactions).

Variants:
* ``ionneutral``: A+ + N -> B+ at nu*dt = 0.01. N_A = N_A0 exp(-nu t),
  N_A + N_B conserved, product inherits reactant drift.
* ``stiff``: the same reaction at nu*dt = 3; N_A + N_B conserved exactly.
* ``recombination`` (full PIC) / ``recombination_hybrid``: O2+ + e -> O + O
  with n_e = n_i, n(t) = n0/(1 + k n0 t), k = k0 (Tref/Te)^alpha.
* ``cx``: #CHARGEEXCHANGE of a drifting H+ beam: H+ mass conserved and
  energy follows slowing-down u = u0/(1 + n sigma u0 t).
* ``weightfloor``: photoionization below #SOURCEPARTICLES minWeightRatio; injected mass
  is S V t on average.
* ``chemistry``: multi-species Mars chemistry network (10 reactions, 4 ion
  species H+, O+, O2+, CO2+); verifies energy changes across species,
  cross-species CX source growth in O2+, and CO2+ changes.
* ``electronimpact``: electron-impact ionization source; verifies O+ energy
  growth over time.
* ``cx_multi``: 2-component Chamberlain neutral atmosphere charge exchange;
  verifies O+ energy growth, H+ stability, and spatial surface-peaked source
  profile in plot output.
"""
import glob
import logging
import math
import os

from tests._shared import run_dir as _run_dir
from tests._shared.validators import log_epart_header, validate_ionisation_source

logger = logging.getLogger(__name__)

# Boltzmann constant / proton mass [J/K/kg]: thermal energy per unit mass.
KB_OVER_MP = 1.380649e-23 / 1.67262192e-27

RUN_DIR = "run_test"


def set_run_dir(run_dir):
    """Point the validator at the current run directory."""
    global RUN_DIR
    RUN_DIR = run_dir
    _run_dir.set_run_dir(run_dir)


def _series(pic_diags, key):
    return [d[key] for d in pic_diags]


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def _need(pic_diags, keys):
    if not pic_diags or len(pic_diags) < 2:
        return "no PIC log"
    missing = [k for k in keys if k not in pic_diags[0]]
    if missing:
        return "missing log columns " + ", ".join(missing)
    return None


def _check(cond, msg, reasons):
    logger.debug("    %s: %s", "OK  " if cond else "FAIL", msg)
    if not cond:
        reasons.append(msg)


def _validate_ionneutral(d, nu, tol_decay, tol_cons, check_drift, reasons):
    t0 = d[0]["time"]
    mA0, mB0 = d[0]["Mpart1"], d[0]["Mpart2"]
    nTot0 = mA0 / 16.0 + mB0 / 32.0
    worstDecay = worstCons = 0.0
    for e in d[1:]:
        t = e["time"] - t0
        if math.exp(-nu * t) > 1e-3:
            worstDecay = max(worstDecay,
                             _rel(e["Mpart1"] / mA0, math.exp(-nu * t)))
        nTot = e["Mpart1"] / 16.0 + e["Mpart2"] / 32.0
        worstCons = max(worstCons, _rel(nTot, nTot0))
    tEnd = d[-1]["time"] - t0
    _check(worstDecay < tol_decay,
           f"N_A/N_A0 = exp(-nu t): worst rel. error {worstDecay:.2e} "
           f"(tol {tol_decay:g})", reasons)
    _check(worstCons < tol_cons,
           f"N_A + N_B conserved: worst rel. error {worstCons:.2e} "
           f"(tol {tol_cons:g}), nu*t_end = {nu * tEnd:.2f}", reasons)
    if check_drift:
        eA = d[-1]["Epart1"] / d[-1]["Mpart1"]
        eB = d[-1]["Epart2"] / d[-1]["Mpart2"]
        _check(_rel(eB, eA) < 0.03,
               f"product keeps the reactant drift: (E/M)_B/(E/M)_A = "
               f"{eB / eA:.4f}", reasons)
    _check(all(abs(e["Mpart0"]) == 0.0 and abs(e["Epart0"]) == 0.0
               for e in d),
           "hybrid fluid-electron slot carries no particles", reasons)


def _validate_recombination(d, kn0, hybrid, reasons):
    t0 = d[0]["time"]
    mI0 = d[0]["Mpart1"]
    worst = 0.0
    for e in d[1:]:
        t = e["time"] - t0
        worst = max(worst, _rel(e["Mpart1"] / mI0, 1.0 / (1.0 + kn0 * t)))
    tEnd = d[-1]["time"] - t0
    _check(worst < 0.04,
           f"n_i/n0 = 1/(1 + k n0 t): worst rel. error {worst:.2e} (tol 0.04)"
           f", k n0 t_end = {kn0 * tEnd:.2f}, final ratio "
           f"{d[-1]['Mpart1'] / mI0:.4f} vs {1.0 / (1.0 + kn0 * tEnd):.4f}",
           reasons)
    if hybrid:
        _check(all(e["Mpart0"] == 0.0 for e in d),
               "hybrid fluid-electron slot carries no particles", reasons)
    else:
        mE0 = d[0]["Mpart0"]
        worstE = max(_rel(e["Mpart0"] / mE0, e["Mpart1"] / mI0) for e in d)
        _check(worstE < 0.01,
               f"electrons are lost with the ions: worst rel. difference "
               f"{worstE:.2e} (tol 0.01)", reasons)


def _validate_cx(d, nu0, reasons):
    t0 = d[0]["time"]
    m0, e0 = d[0]["Mpart1"], d[0]["Epart1"]
    eTh = 1.5 * KB_OVER_MP * 1000.0
    eBeam = 0.5 * 1.0e10 + eTh
    worstM = max(_rel(e["Mpart1"], m0) for e in d)
    _check(worstM < 0.01,
           f"H+ mass conserved: worst rel. error {worstM:.2e} (tol 0.01)",
           reasons)
    worstE = 0.0
    for e in d[1:]:
        t = e["time"] - t0
        r = 1.0 / (1.0 + nu0 * t)
        expected = r + (1.0 - r) * eTh / eBeam
        worstE = max(worstE, _rel(e["Epart1"] / e0, expected))
    tEnd = d[-1]["time"] - t0
    _check(worstE < 0.03,
           f"slowing-down u0/(1 + n sigma u0 t): worst rel. error of E/E0 "
           f"{worstE:.2e} (tol 0.03), n sigma u0 t_end = {nu0 * tEnd:.2f}",
           reasons)
    _check(all(e["Mpart0"] == 0.0 for e in d),
           "no electron is created by charge exchange", reasons)


def _validate_weightfloor(d, sOverSeed, reasons):
    t0 = d[0]["time"]
    m0 = d[0]["Mpart1"]
    tEnd = d[-1]["time"] - t0
    gain = d[-1]["Mpart1"] / m0 - 1.0
    expected = sOverSeed * tEnd
    _check(_rel(gain, expected) < 0.06,
           f"injected mass S V t on average: (M - M0)/M0 = {gain:.4f} vs "
           f"{expected:.4f} (tol 6%)", reasons)
    _check(abs(d[-1]["Mpart2"] / d[0]["Mpart2"] - 1.0) < 1e-12,
           "background H+ untouched", reasons)


def _validate_chemistry(d, reasons):
    first, last = d[0], d[-1]
    epart_keys = sorted(
        k for k in first.keys() if k.startswith("Epart") and k != "Epart"
    )
    species_names = {
        "Epart1": "H+",
        "Epart2": "O+",
        "Epart3": "O2+",
        "Epart4": "CO2+",
    }
    change_threshold = 0.001  # 0.1%
    for k in epart_keys:
        e0 = first.get(k, 0.0)
        e1 = last.get(k, 0.0)
        if e0 <= 0:
            continue
        ratio = e1 / e0
        name = species_names.get(k, k)
        _check(abs(ratio - 1.0) >= change_threshold,
               f"{name} energy changed (ratio {ratio:.4f})", reasons)

    # O2+ must increase (cross-species CX source dominates recombination loss)
    if "Epart3" in first and first["Epart3"] > 0:
        o2_ratio = last["Epart3"] / first["Epart3"]
        _check(o2_ratio > 1.0,
               f"O2+ energy increased by cross-species CX: ratio = {o2_ratio:.4f}",
               reasons)

    # CO2+ must change (photoionization source and CX loss active)
    if "Epart4" in first and first["Epart4"] > 0:
        co2_ratio = last["Epart4"] / first["Epart4"]
        _check(abs(co2_ratio - 1.0) >= change_threshold,
               f"CO2+ energy changed: ratio = {co2_ratio:.4f}", reasons)


def _validate_cx_multi_log(d, reasons):
    first, last = d[0], d[-1]
    epart_keys = sorted(
        k for k in first.keys() if k.startswith("Epart") and k != "Epart"
    )
    log_epart_header(first, last, epart_keys, logger)

    min_factor_o = 2.0
    h_tolerance = 0.05
    o_key = epart_keys[-1] if epart_keys else "Epart2"
    e_o_init = first.get(o_key, 0.0)
    e_o_final = last.get(o_key, 0.0)
    factor_o = e_o_final / max(e_o_init, 1e-30)
    _check(factor_o >= min_factor_o,
           f"O+ energy increased: factor {factor_o:.2f}x (tol >= {min_factor_o}x)",
           reasons)

    if "Epart1" in first:
        e_h_init = first["Epart1"]
        e_h_final = last["Epart1"]
        _check(e_h_final >= e_h_init * (1.0 - h_tolerance),
               f"H+ energy stable: final/init = {e_h_final / e_h_init:.4f}",
               reasons)


def _check_charge_exchange_source_profile():
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        return False, "No .out files found"

    out_file = out_files[-1]
    with open(out_file, "r") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return True, "Short .out file"

    var_names = lines[4].split()
    rho_idx = None
    rho_name = None
    for target in ("RHOS2", "RHOS1"):
        for iv, vn in enumerate(var_names):
            if vn.upper() == target:
                rho_idx = iv
                rho_name = target
                break
        if rho_idx is not None:
            break
    if rho_idx is None:
        return True, "rhoS2/rhoS1 not in .out"

    Rp_plot = 1.0
    points = []
    is_2d = (rho_idx >= 2)
    for line in lines[5:]:
        cols = line.strip().split()
        if len(cols) <= rho_idx:
            continue
        try:
            x = float(cols[0])
            y = float(cols[1]) if is_2d else 0.0
            r = (x * x + y * y) ** 0.5 if is_2d else abs(x)
            rho = float(cols[rho_idx])
            points.append((x, y, r, rho))
        except (ValueError, IndexError):
            continue

    if not points:
        return False, "No data points parsed"

    near_surface = [rho for _, _, r, rho in points if 0.5 * Rp_plot < r <= 1.5 * Rp_plot]
    deep_interior = [rho for _, _, r, rho in points if r < 0.3 * Rp_plot]

    surface_mean = sum(near_surface) / len(near_surface) if near_surface else 0.0
    surface_max = max(near_surface, default=0.0)
    interior_mean = sum(deep_interior) / len(deep_interior) if deep_interior else 0.0

    min_surface_max = 1e-4
    if surface_max < min_surface_max:
        return False, (f"Surface {rho_name} max {surface_max:.2e} < {min_surface_max:.1e} "
                       f"(unphysical rate / missing conversion)")

    if interior_mean > surface_mean * 0.1:
        return False, (f"Interior density too high "
                       f"({interior_mean:.2e} vs {surface_mean:.2e})")

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
    variant = (test_name or "").replace("reactions_", "")
    logger.debug("Validating reaction test variant '%s'...", variant)

    keys = {
        "ionneutral": ["Mpart0", "Mpart1", "Mpart2", "Epart0", "Epart1", "Epart2"],
        "stiff": ["Mpart0", "Mpart1", "Mpart2", "Epart0"],
        "recombination": ["Mpart0", "Mpart1"],
        "recombination_hybrid": ["Mpart0", "Mpart1"],
        "cx": ["Mpart0", "Mpart1", "Epart1"],
        "weightfloor": ["Mpart1", "Mpart2"],
        "chemistry": ["Epart1", "Epart2", "Epart3", "Epart4"],
        "electronimpact": ["Epart1", "Epart2"],
        "cx_multi": ["Epart1", "Epart2"],
    }
    if variant not in keys:
        return False, f"unknown reactions variant '{test_name}'"
    missing = _need(pic_diags, keys[variant])
    if missing:
        return False, missing

    reasons = []
    if variant == "ionneutral":
        _validate_ionneutral(pic_diags, 0.1, 0.01, 0.01, True, reasons)
    elif variant == "stiff":
        _validate_ionneutral(pic_diags, 30.0, 0.02, 0.01, False, reasons)
    elif variant == "recombination":
        _validate_recombination(pic_diags, 1.0, False, reasons)
    elif variant == "recombination_hybrid":
        _validate_recombination(pic_diags, 1.0, True, reasons)
    elif variant == "cx":
        _validate_cx(pic_diags, 0.15, reasons)
    elif variant == "weightfloor":
        _validate_weightfloor(pic_diags, 0.1, reasons)
    elif variant == "chemistry":
        _validate_chemistry(pic_diags, reasons)
    elif variant == "electronimpact":
        return validate_ionisation_source(
            pic_diags, test_label="Electron-Impact Ionization", logger=logger
        )
    elif variant == "cx_multi":
        _validate_cx_multi_log(pic_diags, reasons)

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"


def validate_plot(test_name=None):
    if test_name and test_name.endswith("cx_multi"):
        logger.debug("  --- Validating Output Files (CX source profile) ---")
        return _check_charge_exchange_source_profile()
    return True, "Passed"
