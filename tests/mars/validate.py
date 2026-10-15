#!/usr/bin/env python3
"""Validator for the coarse 3D Mars smoke test (tests/mars).

The run is far from a steady state after 20 s, so the checks are sanity
bounds that the setup bugs fixed on the mars-standalone-fixes branch violated:

* no heavy ion enters through the upstream face (#INFLOW parsing): the O2+
  mass, which no reaction produces from the trace seeds, stays put in step 1;
* the photochemistry runs on the plasma carried by the particles: O2+, the
  product of CO2+ + O and O+ + CO2, outnumbers CO2+ by the end;
* the plasma stays quasi-neutral (|sum q n| < 1% of the solar-wind H+);
* the field energies stay bounded.
"""
import logging

logger = logging.getLogger(__name__)

MASS = [0.0625, 1.0, 32.0, 16.0, 44.0]


def validate_log(pic_diags=None, test_name=None):
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
    check(worstQ < 0.01, f"quasi-neutral: max |net charge|/N(H+) = {worstQ:.2e}")

    for k in ("Eb", "Ee"):
        r = last[k] / first[k]
        check(0.5 < r < 1.6, f"{k} bounded: {k}/{k}0 = {r:.3f}")

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"
