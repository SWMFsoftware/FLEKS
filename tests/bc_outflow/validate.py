#!/usr/bin/env python3
"""Validator for the outflow field-boundary tests (tests/bc_outflow/).

  - PARAM.in.largedt.fields -> base_name "bc_outflow_largedt_fields"
      Vacuum GEM field between outflow y-faces, semi-implicit solver at
      c*dt/dx = 7, theta = 0.5. The zero-gradient outflow faces may raise the
      field energy once (by about 13x here), but it must not keep growing.
"""
import logging
import math

logger = logging.getLogger(__name__)

# Last-step energy over the maximum of the first half of the run. The fixed
# solver gives 0.95; zero ghost nodes in the GMRES operator give ~1e10.
GROWTH_MAX = 1.05
# Overall bound on Etot / Etot_initial (13.3 with the fixed solver).
ETOT_MAX = 100.0


def validate_log(pic_diags=None, test_name=None):
    logger.debug("Validating %s (log)...", test_name)
    if not pic_diags or len(pic_diags) < 10:
        return False, "Energy log missing or too short"

    etot = [d.get("Etot", 0.0) for d in pic_diags]
    if not all(math.isfinite(e) for e in etot):
        return False, "Non-finite field energy (NaN/Inf)"
    if etot[0] <= 0:
        return False, "Initial Etot is zero (field not seeded)"

    half = len(etot) // 2
    early = max(etot[:half])
    growth = etot[-1] / early
    peak = max(etot) / etot[0]
    logger.debug("    Etot_last / max(Etot first half) = %.3e, max Etot/Etot0 = %.3e",
                 growth, peak)

    if growth > GROWTH_MAX or peak > ETOT_MAX:
        return False, (
            f"Field energy keeps growing (last/max(first half) = {growth:.3e}, "
            f"max/initial = {peak:.3e}) -- the GMRES operator does not use the "
            f"zero-gradient ghosts of the outflow faces")
    return True, (f"Passed (energy bounded: max/initial = {peak:.1f}, "
                  f"last/max(first half) = {growth:.2f})")
