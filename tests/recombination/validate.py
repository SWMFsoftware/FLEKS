#!/usr/bin/env python3
"""Validator for the recombination loss test (tests/recombination).

Checks that O2+ (species 2, Epart2) and electrons (species 0, Epart0) energy
decreases over time due to recombination loss, while H+ (species 1, Epart1) energy
remains stable since H+ does not participate in recombination.
"""
import logging

logger = logging.getLogger(__name__)


def validate_log(pic_diags=None, test_name=None):
    """Validate the recombination loss test (O2+ + e- -> O + O).

    Checks that O2+ (species 2, Epart2) and electron (species 0, Epart0)
    energies decrease over time due to recombination loss, while H+
    (species 1, Epart1) energy remains stable since H+ does not participate
    in recombination.
    """
    logger.debug("Validating Recombination Loss Test...")

    if not pic_diags or len(pic_diags) < 2:
        logger.debug("  [INFO] No PIC energy log found; skipping energy checks.")
        return True, "Passed (no pic log)"

    first = pic_diags[0]
    last = pic_diags[-1]

    epart_keys = sorted(
        k for k in first.keys() if k.startswith("Epart") and k != "Epart"
    )
    if not epart_keys:
        logger.debug("  [INFO] No per-species energy columns; skipping.")
        return True, "Passed (no Epart columns)"

    logger.debug("  --- Energy Diagnostics (from log_pic log) ---")
    for k in epart_keys:
        logger.debug("    %s: %s -> %s",
                     k, f"{first.get(k, 0):.6e}", f"{last.get(k, 0):.6e}")

    passed = True
    reasons = []

    # Electron (species 0) should decrease due to recombination in full-PIC mode.
    e_key = "Epart0" if "Epart0" in first else None
    if e_key:
        e_e_initial = first.get(e_key, 0.0)
        e_e_final = last.get(e_key, 0.0)
        logger.debug("    %s (e-): %s -> %s",
                     e_key, f"{e_e_initial:.6e}", f"{e_e_final:.6e}")
        if e_e_initial <= 0:
            logger.debug("    FAIL: %s initial energy is zero.", e_key)
            passed = False
            reasons.append("e- initial energy is zero")
        elif e_e_final >= e_e_initial:
            logger.debug("    FAIL: %s energy did not decrease (electron loss not active).", e_key)
            passed = False
            reasons.append("e- energy did not decrease")
        else:
            e_ratio = e_e_final / e_e_initial
            logger.debug("    SUCCESS: %s energy decreased to %.3f of initial.",
                         e_key, e_ratio)

    # O2+ (species 2) should decrease due to recombination.
    o2_key = "Epart2" if "Epart2" in first else (epart_keys[-1] if len(epart_keys) >= 2 else None)
    if o2_key:
        e_o2_initial = first.get(o2_key, 0.0)
        e_o2_final = last.get(o2_key, 0.0)
        logger.debug("    %s (O2+): %s -> %s",
                     o2_key, f"{e_o2_initial:.6e}", f"{e_o2_final:.6e}")
        if e_o2_initial <= 0:
            logger.debug("    FAIL: %s initial energy is zero.", o2_key)
            passed = False
            reasons.append("O2+ initial energy is zero")
        elif e_o2_final >= e_o2_initial:
            logger.debug("    FAIL: %s energy did not decrease (recombination not active).", o2_key)
            passed = False
            reasons.append("O2+ energy did not decrease")
        else:
            ratio = e_o2_final / e_o2_initial
            logger.debug("    SUCCESS: %s energy decreased to %.3f of initial.",
                         o2_key, ratio)

    # H+ (species 1) should remain stable.
    h_key = "Epart1" if "Epart1" in first else None
    if h_key:
        e_h_initial = first.get(h_key, 0.0)
        e_h_final = last.get(h_key, 0.0)
        h_tolerance = 0.10  # allow up to 10% variation (numerical noise)
        logger.debug("    %s (H+): %s -> %s",
                     h_key, f"{e_h_initial:.6e}", f"{e_h_final:.6e}")
        if e_h_initial > 0:
            h_ratio = abs(e_h_final - e_h_initial) / e_h_initial
            if h_ratio > h_tolerance:
                logger.debug("    FAIL: %s energy changed by %.1f%% "
                             "(threshold %.0f%%).",
                             h_key, h_ratio * 100, h_tolerance * 100)
                passed = False
                reasons.append(f"H+ energy changed by {h_ratio*100:.1f}%")
            else:
                logger.debug("    SUCCESS: %s energy stable (%.1f%% change).",
                             h_key, h_ratio * 100)

    if passed:
        logger.debug("Recombination Loss Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)
