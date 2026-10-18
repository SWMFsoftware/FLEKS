#!/usr/bin/env python3
"""Validator for the floating body (tests/floatingbody).

An absorbing sphere (R = 1.5) sits in a closed box (reflecting walls) of
electron-proton plasma at rest, T = 1e5 K (kT/e = 8.62 V), mi/me = 16, full
PIC with the fields off, so the particles fly ballistically. The electrons
reach the sphere sigma_e/sigma_i = 4 times faster than the ions. In the
closed box the net charge of the plasma is exactly the opposite of the net
charge absorbed by the body (and re-emitted).

* absorb (reference): the body takes 3x more electron than ion charge and the
  plasma charges up (> +5% of N(H+) in 5 s, measured 10%); the bookkeeping
  holds (1%).
* sheath (#FLOATINGBODY sheath, 8 s): after the potential has relaxed (5 s)
  the net absorbed current vanishes (< 10% of the ion current) and the
  plasma charge stays flat (< 0.5% of N(H+)); the potential is of order
  kT/e ln(sigma_e/sigma_i) = 1.4 kT/e (accepted 0.7 - 1.6 kT/e: the closed box
  depletes and the electrons near the sphere are not exactly Maxwellian).
* reemit (#FLOATINGBODY reemit): the re-emitted charge equals the absorbed
  electron minus ion charge (1e-6) and the plasma stays neutral (1e-3).
"""
import logging
import math

logger = logging.getLogger(__name__)

ME = 0.0625
KTE = 8.617  # [V]


def _net(e):
    """Net charge of the plasma (ions - electrons), code units."""
    return e["Mpart1"] - e["Mpart0"] / ME


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
    d = pic_diags
    if not d or len(d) < 40 or "qBodyS0" not in d[0]:
        return False, "missing log or qBodyS columns"
    variant = (test_name or "").rsplit("_", 1)[-1]
    reasons = []

    def check(cond, msg):
        logger.debug("    %s: %s", "OK  " if cond else "FAIL", msg)
        if not cond:
            reasons.append(msg)

    end = d[-1]
    nH = end["Mpart1"]
    absorbed = -(end["qBodyS0"] + end["qBodyS1"])  # e - i, > 0 if e excess

    if variant == "absorb":
        check(_net(end) / nH > 0.05,
              f"plasma net charge/N(H+) = {_net(end) / nH:+.3f} (> 0.05)")
        check(abs(_net(end) - absorbed) < 0.01 * abs(absorbed),
              "plasma net charge = absorbed (e - i) charge "
              f"({_net(end):.4f} vs {absorbed:.4f})")
    elif variant == "sheath":
        mid = next(e for e in d if e["time"] >= 5.0)
        dE = -(end["qBodyS0"] - mid["qBodyS0"])
        dI = end["qBodyS1"] - mid["qBodyS1"]
        check(abs(dE - dI) < 0.1 * dI,
              f"net absorbed current after 5 s = {(dE - dI) / dI:+.3f} of the "
              "ion current (tol 10%)")
        drift = (_net(end) - _net(mid)) / nH
        check(abs(drift) < 0.005,
              f"plasma charge change after 5 s = {drift:+.4f} of N(H+) "
              "(tol 0.5%)")
        phi = [e["phiBody"] / KTE for e in d if e["time"] >= 5.0]
        check(0.7 < min(phi) and max(phi) < 1.6,
              f"phi in [{min(phi):.2f}, {max(phi):.2f}] kT/e (0.7 - 1.6)")
    elif variant == "reemit":
        check(abs(end["qBodyReemit"] - absorbed) < 1e-6 * abs(absorbed),
              f"re-emitted = absorbed (e - i) charge ({end['qBodyReemit']:.6f}"
              f" vs {absorbed:.6f})")
        worst = max(abs(_net(e)) / e["Mpart1"] for e in d)
        check(worst < 1e-3, f"max |plasma net charge|/N(H+) = {worst:.2e}")
    else:
        return False, f"unknown variant {variant}"

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"
