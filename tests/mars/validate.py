#!/usr/bin/env python3
"""Validator for the coarse 3D Mars smoke test (tests/mars).

The run is far from a steady state after 20 s, so the checks are sanity
bounds that the setup bugs fixed on the mars-standalone-fixes branch violated:

* no heavy ion enters through the upstream face (#INFLOW parsing): the O2+
  mass, which no reaction produces from the trace seeds, stays put in step 1;
* the photochemistry runs on the plasma carried by the particles: O2+, the
  product of CO2+ + O and O+ + CO2, outnumbers CO2+ by the end;
* the plasma stays quasi-neutral (|sum q n| < 0.2% of the solar-wind H+;
  8e-4 with the outflowthermal faces, 4e-3 when the electrons drain through
  plain outflow faces);
* the field energies stay bounded.
"""
import glob
import logging
import math
import os

import numpy as np

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"


def set_run_dir(run_dir):
    global RUN_DIR
    RUN_DIR = run_dir

MASS = [0.0625, 1.0, 32.0, 16.0, 44.0]


def _non_finite(diags):
    """Reason string if a log value is NaN/inf (comparisons with NaN are false,
    so a tolerance check alone would pass a broken run), else None."""
    keys = sorted({k for e in (diags or []) for k, v in e.items()
                   if isinstance(v, float) and not math.isfinite(v)})
    if keys:
        return f"non-finite log values in {', '.join(keys[:6])}"
    return None


def validate_log(pic_diags=None, test_name=None):
    bad = _non_finite(pic_diags)
    if bad:
        return False, bad
    if (test_name or "").endswith("emission"):
        return _validate_emission_log(pic_diags)
    d = pic_diags
    if not d or len(d) < 10:
        return False, "no PIC log"
    if "Mpart4" not in d[0]:
        return False, "missing MpartN log columns"

    reasons = []

    def check(cond, msg):
        logger.debug("    %s: %s", "OK  " if cond else "FAIL", msg)
        if not cond:
            reasons.append(msg)

    n = lambda e, i: e[f"Mpart{i}"] / MASS[i]
    first, step1, last = d[0], d[1], d[-1]

    dO2 = step1["Mpart2"] - first["Mpart2"]
    check(abs(dO2) < 0.05 * first["Mpart2"] + 1.0,
          f"no O2+ injected at the inflow face in step 1 (dM = {dO2:.3g})")

    check(n(last, 2) > n(last, 4),
          f"O2+ outnumbers CO2+ at t = {last['time']:.0f} s "
          f"({n(last, 2):.4g} vs {n(last, 4):.4g})")

    worstQ = 0.0
    for e in d:
        q = -n(e, 0) + sum(n(e, i) for i in range(1, 5))
        worstQ = max(worstQ, abs(q) / n(e, 1))
    check(worstQ < 2e-3, f"quasi-neutral: max |net charge|/N(H+) = {worstQ:.2e}")

    for k in ("Eb", "Ee"):
        r = last[k] / first[k]
        check(0.5 < r < 1.6, f"{k} bounded: {k}/{k}0 = {r:.3f}")

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"


# ---------------------------------------------------------------------------
# fullpic_emission: the final boundary setup (run_mars/final_coarse), 40 s
# ---------------------------------------------------------------------------
FACES = ["xm", "xp", "ym", "yp", "zm", "zp"]
RM = 3396e3          # Mars radius [m]
MP = 1.67262192e-27  # [kg]
NSW = 3.0e6          # solar-wind H+ [m^-3]


def _validate_emission_log(d):
    """Log checks of the emission-shell / floating-body / open-face setup.

    * quasi-neutral: max |net charge| < 0.5% of N(H+) (measured 0.16%);
    * floating body: the re-emitted charge equals the net absorbed charge
      (reemit returns every step's electron excess);
    * H+ budget: dM - (face in - face out - body) is the net H+ source, the
      corona (+1.2 code mass/s), constant; checked in [0.8, 1.6];
    * the planetary ions accumulate: emitted O2+ minus what fell back;
    * the field and particle energies stay bounded (Ee peaks near 18 s while
      the shock forms).
    """
    if not d or len(d) < 100 or "mIn_xm_S1" not in d[0]:
        return False, "missing log or face budget columns"
    reasons = []

    def check(cond, msg):
        logger.debug("    %s: %s", "OK  " if cond else "FAIL", msg)
        if not cond:
            reasons.append(msg)

    n = lambda e, i: e[f"Mpart{i}"] / MASS[i]
    first, last = d[0], d[-1]

    worstQ = max(abs(sum(n(e, i) for i in range(1, 5)) - n(e, 0)) / n(e, 1)
                 for e in d)
    check(worstQ < 5e-3,
          f"quasi-neutral: max |net charge|/N(H+) = {worstQ:.2e} (< 5e-3)")

    netAbs = -last["qBodyAbsorb"]
    check(abs(last["qBodyReemit"] - netAbs) < 1e-3 * abs(netAbs),
          f"re-emitted = net absorbed charge ({last['qBodyReemit']:.4g} vs "
          f"{netAbs:.4g})")

    # H+ budget over the second half (after the start-up absorption).
    a = next(e for e in d if e["time"] >= last["time"] / 2)
    dt = last["time"] - a["time"]
    faces = sum((last[f"mIn_{f}_S1"] - a[f"mIn_{f}_S1"]) -
                (last[f"mOut_{f}_S1"] - a[f"mOut_{f}_S1"]) for f in FACES)
    body = last["qBodyS1"] - a["qBodyS1"]
    src = (last["Mpart1"] - a["Mpart1"] - faces + body) / dt
    check(0.8 < src < 1.6,
          f"H+ net source (corona) = {src:.3f} code mass/s (0.8 - 1.6)")

    o2 = n(last, 2) - n(first, 2)
    check(100 < o2 < 360,
          f"O2+ in the plasma = {o2:.1f} (emitted 360, less what fell back)")

    rEb = last["Eb"] / first["Eb"]
    check(1.5 < rEb < 6, f"Eb/Eb0 = {rEb:.2f} (pileup: 1.5 - 6)")
    rEe = max(e["Ee"] for e in d) / first["Ee"]
    check(rEe < 50, f"max Ee/Ee0 = {rEe:.1f} (< 50)")
    rEp = last["Epart"] / first["Epart"]
    check(0.95 < rEp < 1.05, f"Epart/Epart0 = {rEp:.3f} (0.95 - 1.05)")

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"


def _read_out(path):
    with open(path, "r", encoding="latin-1") as f:
        lines = f.readlines()
    names = lines[4].split()
    data = np.array([[float(v) for v in l.split()] for l in lines[5:]
                     if l.split()])
    return names, data


def validate_plot(test_name=None):
    """Bow shock and pileup at the end of fullpic_emission (z = 0 cut, SI).

    * bow shock: going sunward along the subsolar line, the last point where
      n(H+) > 1.5 n_sw lies between 1.3 and 2.1 R_M (1.75 measured at 40 s,
      steady to 200 s);
    * pileup: max |B| on the dayside within 2 R_M exceeds 15 nT (30 - 40 nT;
      the absorbing body alone gives ~5 nT).
    """
    if not (test_name or "").endswith("emission"):
        return True, "no plot check"
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots",
                                          "z=0*.out")))
    if not plots:
        return False, "no z=0 plot"
    names, d = _read_out(plots[-1])
    x, y = d[:, 0] / RM, d[:, 1] / RM
    nH = d[:, names.index("rhoS1")] / MP
    B = np.sqrt(sum(d[:, names.index(c)]**2 for c in ("Bx", "By", "Bz")))
    reasons = []

    line = (np.abs(y) < 0.2) & (x > 0)
    xs, ns = x[line], nH[line]
    o = np.argsort(-xs)
    shock = next((xx for xx, nn in zip(xs[o], ns[o]) if nn > 1.5 * NSW),
                 np.nan)
    logger.debug("    bow shock at %.2f R_M", shock)
    if not 1.3 < shock < 2.1:
        reasons.append(f"bow shock at {shock:.2f} R_M (1.3 - 2.1)")

    day = (x > 0) & (np.hypot(x, y) < 2.0)
    bMax = B[day].max() * 1e9
    logger.debug("    dayside max |B| = %.1f nT", bMax)
    if not bMax >= 15:
        reasons.append(f"dayside max |B| = {bMax:.1f} nT (> 15)")

    if reasons:
        return False, "; ".join(reasons)
    return True, f"Passed (shock {shock:.2f} R_M, pileup {bMax:.0f} nT)"
