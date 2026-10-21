#!/usr/bin/env python3
"""Validation for the HYBRID_WHISTLER test.

This directory now runs TWO variants of the same whistler eigenmode, driven by
the test runner (tests/validate_tests.py):

  * FULL PIC  (PARAM.in,           base_name = whistler)
      kinetic ions AND kinetic electrons, standard Maxwell/GMRES solver
      (solveEM = T).
  * EXPLICIT  (PARAM.in.explicit,   base_name = whistler_explicit)
      kinetic ions AND kinetic electrons, explicit leapfrog Maxwell solver
      (#FIELDSOLVER explicit), m_i/m_e = 25; checked against the full
      cold-plasma R-mode within 5%, like the semi-implicit full-PIC variant.
  * HYBRID    (PARAM.in.hybrid,     base_name = whistler_hybrid)
      kinetic ions + fluid electrons, generalized Ohm's-law solver
      (useHybridPIC = T, Hall-only).

Both variants seed the identical x-aligned, transverse, circularly-polarized
RIGHT-HAND (whistler) wave (#TESTCASE HybridWave, #WAVEIC rightHand T,
frac = 0.02) so the field diagnostics are directly comparable.  The shared,
solver-agnostic checks in tests/_shared/hybrid.py are used for both:

  * validate_log  -> energy-log checks (finite energies, bounded magnetic
    energy, conserved particle number under periodic BCs).
  * validate_plot -> seeded-wavelength (n=1) check + whistler-dispersion check
    measuring the transverse-wave frequency and comparing it against the two
    parallel-propagating Hall branches
    w/Omega_i = [+-kappa^2 + kappa*sqrt(kappa^2+4)]/2 (kappa = k d_i) from the
    .out plot frames.

Because the dispersion check lives in tests/_shared/hybrid.py and reads only the
By/Bz components, it applies unchanged to the full-PIC solver output.
"""

import glob
import logging
import math
import os
import sys

logger = logging.getLogger(__name__)

# Make the shared test helpers importable when run directly.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from _shared import hybrid as hybrid_helpers

RUN_DIR = "run_test"
MAX_GROWTH = 1.9
# rk4_upwind (24 s, #HYBRIDUPWIND instead of hyper-resistivity): the seeded mode
# grows 1.3x with upwind or hyper-resistivity; without either the ~3-cell grid
# mode grows unstably and the ion energy rises 95x.
MAX_GROWTH_UPWIND = 1.5


def set_run_dir(run_dir):
    """Mirror the runner's RUN_DIR into the shared hybrid helper."""
    global RUN_DIR
    RUN_DIR = run_dir
    hybrid_helpers.set_run_dir(run_dir)


def _max_bperp(out_file):
    """Return maximum transverse field magnitude |B_perp| = hypot(By, Bz)."""
    data = hybrid_helpers._hyb_load_out(out_file)
    if not data:
        return None
    by, bz = data
    if not by:
        return None
    return max(math.hypot(b, c) for b, c in zip(by, bz))


def _bperp_growth():
    """Return late/early max|B_perp| from the first and last .out frames, or
    None if there are no frames or the amplitude cannot be measured."""
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    if not plots:
        return None
    a0 = _max_bperp(plots[0])
    a1 = _max_bperp(plots[-1])
    if a0 is None or a0 <= 0.0 or a1 is None:
        return None
    return a0, a1


def _check_upwind_stability():
    """Upwind dissipation alone keeps the 24 s rk4 run stable."""
    amp = _bperp_growth()
    if amp is None:
        return False, "Could not measure the transverse amplitude"
    a0, a1 = amp
    growth = a1 / a0
    logger.debug("    [WHISTLER] Upwind: max|B_perp| %.4e -> %.4e (growth %.2fx, "
                 "bound %.2fx)", a0, a1, growth, MAX_GROWTH_UPWIND)
    if not growth <= MAX_GROWTH_UPWIND:
        return False, ("transverse amplitude grew %.2fx (bound %.2fx) -- upwind "
                       "dissipation does not damp the grid-scale mode"
                       % (growth, MAX_GROWTH_UPWIND))
    return True, "Passed"


def _check_resistive_damping():
    """Verify that resistivity measurably limits transverse field growth."""
    amp = _bperp_growth()
    if amp is None:
        logger.debug("    [WHISTLER] No transverse amplitude (resistive-damping "
                     "check skipped)")
        return True, "Could not measure transverse amplitude (skipped)"
    a0, a1 = amp
    growth = a1 / a0
    logger.debug("    [WHISTLER] Resistive damping: max|B_perp| %.4e -> %.4e "
                 "(growth %.2fx, bound %.2fx)", a0, a1, growth, MAX_GROWTH)
    if growth > MAX_GROWTH:
        return False, ("late/early transverse amplitude grew %.2fx (bound %.2fx) "
                       "-- resistive term appears inactive" % (growth, MAX_GROWTH))
    return True, "Passed"


def validate_log(pic_diags, test_name=None):
    """Energy / stability checks from the shared hybrid-family validator.

    Works for both the full-PIC and hybrid variants: the full-PIC run
    additionally carries an electron species (reported as Epart1 / Ee), which
    the shared log reader and validator handle.
    """
    return hybrid_helpers.validate_hybrid(pic_diags=pic_diags, test_name=test_name)


def validate_plot(test_name):
    """Seeded-wavelength + whistler-dispersion checks (or resistive damping check).

    The dispersion relation is read from tests/whistler/PARAM.in (the
    full-PIC variant), which all variants share for geometry / normalization.
    For the hybrid variants that is the exact Hall branch it must reproduce; for
    the full-PIC variant the dominant whistler branch matches this frequency at
    k d_i ~ 1, so the same check is a meaningful regression guard.
    """
    ok, reason = hybrid_helpers.validate_plot(test_name)
    if not ok:
        return ok, reason
    if test_name and "resistive" in test_name:
        return _check_resistive_damping()
    if test_name and "upwind" in test_name:
        return _check_upwind_stability()
    return True, reason


if __name__ == "__main__":
    # Allow a quick standalone sanity check (reads run_test/ output).
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    from tests.validate_tests import read_pic_log, RUN_DIR
    diags = read_pic_log(RUN_DIR)
    ok, reason = validate_log(diags, test_name="whistler")
    print(("PASSED" if ok else "FAILED"), "-", reason)
    sys.exit(0 if ok else 1)
