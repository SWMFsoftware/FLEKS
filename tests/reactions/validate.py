#!/usr/bin/env python3
"""Validator for the reaction tests (tests/reactions).

Every variant is a uniform, periodic 1D plasma in which one reaction runs
against its analytic solution. The checks use the per-species particle mass
(``MpartN``) and energy (``EpartN``) columns of the PIC log, so they measure
what the particles really carry:

* ``ionneutral``: A+ + N -> B+ at nu*dt = 0.01. N_A = N_A0 exp(-nu t),
  N_A + N_B conserved, the product inherits the reactant drift, and the hybrid
  fluid-electron slot (species 0) has no particles.
* ``stiff``: the same reaction at nu*dt = 3; N_A + N_B is conserved exactly.
* ``recombination`` (full PIC) / ``recombination_hybrid``: O2+ + e -> O + O
  with n_e = n_i, n(t) = n0/(1 + k n0 t), k = k0 (Tref/Te)^alpha.
* ``cx``: #CHARGEEXCHANGE of a drifting H+ beam: the H+ mass is constant and
  the energy follows the slowing-down u = u0/(1 + n sigma u0 t).
* ``weightfloor``: photoionization below the #SOURCEWEIGHTFLOOR floor; the
  injected mass is S V t on average.

All of them read the rates from the plasma carried by the particles, so they
fail if the source terms use the frozen #UNIFORMSTATE instead.
"""
import logging
import math

logger = logging.getLogger(__name__)

# Boltzmann constant / proton mass [J/K/kg]: thermal energy per unit mass.
KB_OVER_MP = 1.380649e-23 / 1.67262192e-27


def _series(pic_diags, key):
    return [d[key] for d in pic_diags]


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def _need(pic_diags, keys):
    if not pic_diags or len(pic_diags) < 3:
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
        # Below ~1e-6 of a reference particle the decayed weights play Russian
        # roulette (apply_loss), so the decay is only compared well above it.
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
        # Energy per mass = u^2/2 + 3/2 kT/m: the products must carry the
        # 100 km/s drift of the reactant (thermal part ~1e-4 of the total).
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
        # The electrons recombine with the ions one for one.
        mE0 = d[0]["Mpart0"]
        worstE = max(_rel(e["Mpart0"] / mE0, e["Mpart1"] / mI0) for e in d)
        _check(worstE < 0.01,
               f"electrons are lost with the ions: worst rel. difference "
               f"{worstE:.2e} (tol 0.01)", reasons)


def _validate_cx(d, nu0, reasons):
    t0 = d[0]["time"]
    m0, e0 = d[0]["Mpart1"], d[0]["Epart1"]
    # Energy per mass of the beam and of the new ions (T0 = 1000 K, m = 1).
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


def validate_log(pic_diags=None, test_name=None):
    variant = (test_name or "").replace("reactions_", "")
    logger.debug("Validating reaction test variant '%s'...", variant)

    keys = {
        "ionneutral": ["Mpart0", "Mpart1", "Mpart2", "Epart0", "Epart1",
                       "Epart2"],
        "stiff": ["Mpart0", "Mpart1", "Mpart2", "Epart0"],
        "recombination": ["Mpart0", "Mpart1"],
        "recombination_hybrid": ["Mpart0", "Mpart1"],
        "cx": ["Mpart0", "Mpart1", "Epart1"],
        "weightfloor": ["Mpart1", "Mpart2"],
    }
    if variant not in keys:
        return False, f"unknown reactions variant '{test_name}'"
    missing = _need(pic_diags, keys[variant])
    if missing:
        return False, missing

    reasons = []
    if variant == "ionneutral":
        # nu = k n_N = 1e-16 m^3/s * 1e15 m^-3
        _validate_ionneutral(pic_diags, 0.1, 0.01, 0.01, True, reasons)
    elif variant == "stiff":
        # nu = 3e-14 m^3/s * 1e15 m^-3, nu*dt = 3
        _validate_ionneutral(pic_diags, 30.0, 0.02, 0.01, False, reasons)
    elif variant == "recombination":
        # k = 0.4*(2500/1e4)^0.5 = 0.2 cm^3/s, n0 = 5 /cc
        _validate_recombination(pic_diags, 1.0, False, reasons)
    elif variant == "recombination_hybrid":
        _validate_recombination(pic_diags, 1.0, True, reasons)
    elif variant == "cx":
        # n sigma u0 = 1e15 m^-3 * 1.5e-21 m^2 * 1e5 m/s
        _validate_cx(pic_diags, 0.15, reasons)
    elif variant == "weightfloor":
        # S/n_seed = 1e3 m^-3 s^-1 / 1e4 m^-3
        _validate_weightfloor(pic_diags, 0.1, reasons)

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"
