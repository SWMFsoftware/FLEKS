#!/usr/bin/env python3
"""Validator for the explicit-PIC thermal-plasma test (tests/thermal_plasma).

A uniform, unmagnetized, thermal electron-ion plasma has nothing to drive it,
so the explicit solver must conserve the total (field + particle) energy. Two
variants, discovered from this directory:
  - PARAM.in.leapfrog -> leapfrog field integrator: |Etot/Etot0 - 1| < 1%
  - PARAM.in.leapfrog_mag -> the same in a uniform magnetic field at 45 degrees
                         (electron w_ce = 0.3 w_pe; B0 is ~16% of Etot): < 1%
  - PARAM.in.pc       -> predictor-corrector integrator (slightly dissipative):
                         |Etot/Etot0 - 1| < 2%
The bound applies at every logged step. Before the fix of FLEKS issue 170 the
leapfrog variant gained a factor of ~560 in 500 steps (the explicit Maxwell
solver used the mass flux instead of the charge current, centered at t^n).
"""
import logging
import math

from tests._shared import run_dir as _run_dir

logger = logging.getLogger(__name__)

set_run_dir = _run_dir.set_run_dir


def validate_log(pic_diags=None, test_name=None):
    tol = 0.02 if test_name and test_name.endswith("pc") else 0.01
    if not pic_diags or len(pic_diags) < 10:
        n = len(pic_diags) if pic_diags else 0
        return False, f"Expected >= 10 energy log entries, got {n}"

    e0 = pic_diags[0].get("Etot", 0.0)
    if not (math.isfinite(e0) and e0 > 0.0):
        return False, f"Initial Etot = {e0} is not finite and positive"

    worst = 0.0
    for d in pic_diags:
        e = d.get("Etot", float("nan"))
        if not math.isfinite(e):
            return False, f"Non-finite Etot at cycle {d.get('cycle')}"
        worst = max(worst, abs(e / e0 - 1.0))
    last = pic_diags[-1]["Etot"] / e0
    logger.debug("    Etot_final/Etot0 = %.5f, max |Etot/Etot0 - 1| = %.5f",
                 last, worst)
    if worst > tol:
        return False, (f"Energy not conserved: max |Etot/Etot0 - 1| = "
                       f"{worst:.4f} > {tol} (final ratio {last:.4f})")
    return True, f"Passed (max |Etot/Etot0 - 1| = {worst:.4f} < {tol})"


def validate_plot(test_name):
    return True, "No plot check (energy log only)"
