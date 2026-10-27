#!/usr/bin/env python3
"""Validator for the open particle boundary benchmark (tests/openboundary).

Evaluates plasma preservation across open boundary faces:
1. Mass of each species (Mpart_s(t)/Mpart_s(0)) within 3% tolerance.
2. At rest: Temperature (Epart_s/Mpart_s) relative to initial value within 3%.
3. Face mass budget closure (dM - (inflow - outflow) < 1e-9).
"""
import logging
import math

logger = logging.getLogger(__name__)


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
    d = pic_diags
    if not d or len(d) < 50:
        return False, "missing or incomplete log"

    reasons = []

    # 1. Check species mass preservation within 3%
    for iS in (0, 1):
        key = f"Mpart{iS}" if f"Mpart{iS}" in d[0] else f"Epart{iS}"
        if key in d[0] and d[0][key] > 0:
            ratios = [e[key] / d[0][key] for e in d]
            worst = max(ratios, key=lambda r: abs(r - 1))
            logger.debug("    species %d (%s): ratio in [%.4f, %.4f], end %.4f",
                         iS, key, min(ratios), max(ratios), ratios[-1])
            if abs(worst - 1) > 0.03:
                reasons.append(f"species {iS} ({key}): ratio = {worst:.4f} (tol 3%)")

    # 2. At rest, check temperature (energy per mass) within 3%
    isDrift = "drift" in (test_name or "")
    if not isDrift and "Mpart0" in d[0] and "Epart0" in d[0]:
        for iS in (0, 1):
            e0 = d[0][f"Epart{iS}"] / d[0][f"Mpart{iS}"]
            ratios = [e[f"Epart{iS}"] / e[f"Mpart{iS}"] / e0 for e in d]
            worst = max(ratios, key=lambda r: abs(r - 1))
            logger.debug("    species %d: T/T0 in [%.4f, %.4f], end %.4f",
                         iS, min(ratios), max(ratios), ratios[-1])
            if abs(worst - 1) > 0.03:
                reasons.append(f"species {iS}: T/T0 = {worst:.4f} (tol 3%)")

    # 3. Check face budget closure if face columns are logged
    faces = ["xm", "xp", "ym", "yp", "zm", "zp"]
    if "mIn_ym_S0" in d[0]:
        for iS in (0, 1):
            net = sum(d[-1][f"mIn_{f}_S{iS}"] - d[-1][f"mOut_{f}_S{iS}"]
                      for f in faces if f"mIn_{f}_S{iS}" in d[-1])
            dM = d[-1][f"Mpart{iS}"] - d[0][f"Mpart{iS}"]
            inflow = sum(d[-1][f"mIn_{f}_S{iS}"] for f in faces if f"mIn_{f}_S{iS}" in d[-1])
            if inflow > 0:
                err = abs(dM - net) / inflow
                logger.debug("    species %d: budget residual %.2e", iS, err)
                if err > 1e-9:
                    reasons.append(f"species {iS}: dM - (in - out) = {err:.2e} "
                                   f"of injected mass (tol 1e-9)")
        if (test_name or "").endswith("rest"):
            xFlux = sum(d[-1].get(f"m{io}_{f}_S{iS}", 0.0) for io in ("In", "Out")
                        for f in ("xm", "xp") for iS in (0, 1))
            if xFlux != 0:
                reasons.append(f"periodic x faces carry a flux ({xFlux:.3g})")

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"
