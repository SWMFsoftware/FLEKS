#!/usr/bin/env python3
"""Shared primitive validator helpers for the FLEKS standalone test suite.

This module collects small, reusable building-blocks that were previously
duplicated verbatim across several per-test ``validate.py`` modules.  Every
helper returns ``(passed: bool, reason: str)`` so callers can forward the
result directly to the runner.

Exported helpers
----------------
log_epart_header(first, last, epart_keys, logger)
    Emit the standard per-species energy-diagnostics DEBUG block.

validate_ionisation_source(pic_diags, *, test_label, source_key, logger)
    Check that the heaviest ion species (source_key) energy grew over time.
    Shared by ``electronimpact`` and ``photoionization``.

validate_frozen_field(pic_diags, *, eb_tol, ee_rel_max, test_label, logger)
    Check Eb conserved + Ee ≈ 0.
    Shared by ``singlecell`` (ee_rel_max=1e-2) and ``zerocurrent``
    (ee_rel_max=1e-6), with different ``eb_tol`` windows.

validate_hybrid_energy_bounded(pic_diags, *, growth_max, test_label, logger)
    Check all log entries finite, initial/final energy positive, and total
    energy growth bounded.  Shared by ``bc_absorb`` and ``bc_reflecting``
    hybrid-field variants.
"""
import logging
import math

_module_logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# log_epart_header
# ---------------------------------------------------------------------------
def log_epart_header(first, last, epart_keys, logger=None):
    """Emit the standard per-species energy-diagnostics DEBUG block.

    Logs one line per species showing initial→final energy, plus the total
    Epart initial and final values.  Call this near the top of a
    ``validate_log`` that works with per-species ``EpartN`` columns.

    Parameters
    ----------
    first, last : dict
        First and last rows of ``pic_diags`` (list of dicts from
        ``read_pic_log``).
    epart_keys : list[str]
        Sorted list of per-species keys, e.g. ``["Epart0", "Epart1", ...]``.
    logger : logging.Logger, optional
        Logger to use; defaults to this module's logger.
    """
    lg = logger or _module_logger
    lg.debug("  --- Energy Diagnostics (from log_pic log) ---")
    for k in epart_keys:
        lg.debug("    %s: %s -> %s",
                 k, f"{first.get(k, 0):.6e}", f"{last.get(k, 0):.6e}")
    lg.debug("    Initial total Epart: %s", f"{first.get('Epart', 0):.6e}")
    lg.debug("    Final total Epart:   %s", f"{last.get('Epart', 0):.6e}")


# ---------------------------------------------------------------------------
# validate_ionisation_source
# ---------------------------------------------------------------------------
def validate_ionisation_source(pic_diags, *, test_label="Ionization Source",
                                source_key=None, logger=None):
    """Check that the heaviest ion species energy increased over time.

    Used by ``electronimpact`` and ``photoionization``: both seed a heavy
    ion (O+) population that should gain energy from the ionization source
    term.  The last ``EpartN`` key in the log is assumed to be the heaviest
    species unless ``source_key`` is given explicitly.

    Parameters
    ----------
    pic_diags : list[dict]
        Output of ``read_pic_log``.
    test_label : str
        Short name used in log messages (e.g. ``"Electron-Impact Ionization"``).
    source_key : str or None
        Explicit key such as ``"Epart2"``; if None the last ``EpartN`` key
        found in the first row is used.
    logger : logging.Logger, optional

    Returns
    -------
    (bool, str)
        ``(True, "Passed")`` or ``(False, reason)``.
    """
    lg = logger or _module_logger
    lg.debug("Validating %s Test...", test_label)

    if not pic_diags or len(pic_diags) < 2:
        lg.debug("  [INFO] No PIC energy log found; skipping energy checks.")
        return True, "Passed (no pic log)"

    first = pic_diags[0]
    last = pic_diags[-1]

    # Determine the source (heaviest ion) species key from EpartN columns.
    epart_keys = sorted(
        k for k in first.keys() if k.startswith("Epart") and k != "Epart"
    )
    if not epart_keys:
        lg.debug("  [INFO] No per-species energy columns; skipping.")
        return True, "Passed (no Epart columns)"

    if source_key is None:
        source_key = epart_keys[-1]   # last EpartN = heaviest ion
    source_idx = source_key.replace("Epart", "")

    log_epart_header(first, last, epart_keys, lg)

    e_src_initial = first.get(source_key, 0.0)
    e_src_final = last.get(source_key, 0.0)
    lg.debug("    Initial %s (species %s, O+): %s",
             source_key, source_idx, f"{e_src_initial:.6e}")
    lg.debug("    Final   %s (species %s, O+): %s",
             source_key, source_idx, f"{e_src_final:.6e}")
    lg.debug("    Growth factor: %.3fx",
             e_src_final / max(e_src_initial, 1e-30))

    if e_src_final <= e_src_initial:
        lg.debug("    FAIL: %s energy did not increase.", source_key)
        lg.debug("    Ionization source may not be working correctly.")
        return False, (
            f"{source_key} energy did not increase "
            f"(initial={e_src_initial:.2e}, final={e_src_final:.2e})"
        )
    lg.debug("    SUCCESS: %s energy increased (ionization source active).",
             source_key)
    return True, "Passed"


# ---------------------------------------------------------------------------
# validate_frozen_field
# ---------------------------------------------------------------------------
def validate_frozen_field(pic_diags, *, eb_tol=(0.999, 1.001),
                           ee_rel_max=1e-6,
                           test_label="Frozen Field", logger=None):
    """Check that the magnetic field is truly frozen (Eb conserved, Ee ≈ 0).

    Used by ``singlecell`` and ``zerocurrent``.  The two tests share the same
    structure but differ in tolerance:

    * singlecell  : ``eb_tol=(0.9999, 1.0001)``, ``ee_rel_max=1e-2``
    * zerocurrent : ``eb_tol=(0.999,  1.001 )``, ``ee_rel_max=1e-6``

    Parameters
    ----------
    pic_diags : list[dict]
        Output of ``read_pic_log``.
    eb_tol : (float, float)
        Acceptable range ``(lo, hi)`` for the Eb_final / Eb_initial ratio.
    ee_rel_max : float
        Maximum acceptable Ee / Eb_initial ratio (must be ≈ 0).
    test_label : str
        Short name used in log messages.
    logger : logging.Logger, optional

    Returns
    -------
    (bool, str)
    """
    lg = logger or _module_logger
    lg.debug("=== Validating %s Test ===", test_label)

    if not pic_diags:
        return False, "No diagnostics found"

    first = pic_diags[0]
    last = pic_diags[-1]
    passed = True
    reasons = []

    eb0 = first.get("Eb", 0.0)
    eb1 = last.get("Eb", 0.0)
    lg.debug("    Eb (magnetic): %s -> %s", f"{eb0:.6e}", f"{eb1:.6e}")

    if not math.isfinite(eb1):
        passed = False
        reasons.append("Eb not finite (NaN/Inf)")
    if eb0 > 0:
        ratio = eb1 / eb0
        lg.debug("    Eb ratio: %.6f", ratio)
        lo, hi = eb_tol
        if ratio < lo or ratio > hi:
            passed = False
            reasons.append(
                f"Eb ratio {ratio:.6f} not in [{lo}, {hi}] "
                f"(field should be frozen)")

    eemax = max((d.get("Ee", 0.0) for d in pic_diags), default=0.0)
    lg.debug("    Ee (electric, max): %s", f"{eemax:.6e}")
    if not math.isfinite(eemax):
        passed = False
        reasons.append("Ee not finite (NaN/Inf)")
    if eb0 > 0 and eemax > ee_rel_max * eb0:
        passed = False
        reasons.append(
            f"Ee {eemax:.3e} not ~0 vs Eb {eb0:.3e} "
            f"(field is evolving / propagating, not frozen)")

    if passed:
        lg.debug("%s Test: PASSED", test_label)
        return True, "Passed"
    return False, "; ".join(reasons)


# ---------------------------------------------------------------------------
# validate_hybrid_energy_bounded
# ---------------------------------------------------------------------------
def validate_hybrid_energy_bounded(pic_diags, *, growth_max=1e3,
                                    test_label="Hybrid Energy", logger=None):
    """Check that the hybrid-PIC run stays finite, positive, and bounded.

    Used by the hybrid-field variants of ``bc_absorb`` and ``bc_reflecting``:
    the field wall is only required to be stable (no blow-up), not to conserve
    energy precisely.

    Parameters
    ----------
    pic_diags : list[dict]
        Output of ``read_pic_log``.
    growth_max : float
        Maximum allowed Etot_final / Etot_initial ratio.
    test_label : str
        Short name for log messages (e.g. ``"Absorbing Field Wall"``).
    logger : logging.Logger, optional

    Returns
    -------
    (bool, str)
    """
    lg = logger or _module_logger

    e0 = pic_diags[0].get("Etot", 0.0)
    e1 = pic_diags[-1].get("Etot", 0.0)
    finite = all(
        math.isfinite(d.get("Etot", 0.0)) and
        math.isfinite(d.get("Epart", 0.0)) for d in pic_diags
    )
    lg.debug("    Etot: %.4e -> %.4e (%.1f growth)", e0, e1,
             1.0 if e0 == 0 else e1 / e0)

    if not finite:
        return False, "Non-finite energy (NaN/Inf) in energy log"
    if e0 <= 0 or e1 <= 0:
        return False, "Non-positive total energy (plasma not initialised / drained)"
    if e0 > 0 and e1 > growth_max * e0:
        return False, (
            f"Etot grew from {e0:.3e} to {e1:.3e} "
            f"(>{growth_max:.0f}x) -- {test_label} instability"
        )
    return True, f"Passed (finite, bounded energy => {test_label} is stable)"
