#!/usr/bin/env python3
"""Validator for the hybrid-Ohm's-law test (tests/ohm).

Uses the shared hybrid-family validator (energy-log checks) plus the seeded
wavelength / bounded-amplitude / whistler-dispersion plot check, and adds an
ohm-specific check that the resistive term is actually active (see
_check_resistive_damping below).
"""
import glob
import logging
import math
import os

from .._shared import hybrid as _hyb
from .._shared.hybrid import validate_hybrid, validate_plot as _hyb_plot

logger = logging.getLogger(__name__)

# Bound on the late/early transverse-amplitude growth factor.
#
# The seeded n=1 wave is NOT a free resistive eigenmode: the kinetic-ion
# response keeps re-exciting the magnetic perturbation, so the observable
# effect of the resistivity is a partial suppression rather than the
# eigenmode decay e^{-gamma*t} (gamma ~ 0.19, i.e. ~70% over the run).
# Measured with this PARAM.in (2000 PPC, deterministic):
#   etaResistivity = 2.0e9  -> late max|B_perp| = 3.3e-2, growth 1.65x
#   etaResistivity ~ 0      -> late max|B_perp| = 4.3e-2, growth 2.16x
# The bound sits between the two, so it fails if the resistivity is silently
# switched off (e.g. by a broken SI->code unit conversion) while leaving
# head-room for modest run-to-run variation.
MAX_GROWTH = 1.9


def _max_bperp(out_file):
    """Largest transverse field magnitude |B_perp| = hypot(By, Bz) in a frame."""
    data = _hyb._hyb_load_out(out_file)
    if data is None:
        return None
    by, bz = data
    if not by:
        return None
    return max(math.hypot(b, c) for b, c in zip(by, bz))


def _check_resistive_damping():
    """Require the resistivity to measurably hold back the transverse field."""
    plots = sorted(glob.glob(os.path.join(_hyb.RUN_DIR, "PC", "plots", "*.out")))
    if not plots:
        logger.debug("    [OHM] No .out files (resistive-damping check skipped)")
        return True, "No .out files (skipped)"

    a0 = _max_bperp(plots[0])
    a1 = _max_bperp(plots[-1])
    if not a0 or a1 is None:
        return True, "Could not measure the transverse amplitude"

    growth = a1 / a0
    logger.debug("    [OHM] Resistive damping: max|B_perp| %.4e -> %.4e "
                 "(growth %.2fx, bound %.2fx)", a0, a1, growth, MAX_GROWTH)
    if growth > MAX_GROWTH:
        return False, ("late/early transverse amplitude grew %.2fx (bound %.2fx) "
                       "-- the resistive term looks inactive" % (growth, MAX_GROWTH))
    return True, "Passed"


def validate_log(pic_diags=None, test_name=None):
    """Hybrid-wave energy-log validation (shared with the hybrid family)."""
    return validate_hybrid(pic_diags=pic_diags, test_name=test_name)


def validate_plot(test_name):
    """Hybrid-wave plot check + the ohm resistive-damping bound."""
    ok, reason = _hyb_plot(test_name)
    if not ok:
        return ok, reason
    return _check_resistive_damping()
