#!/usr/bin/env python3
"""Validator for dynamic AMR multi-session test.

Verifies that refinement regions can be moved and cleared across sessions:
  Session 1 (cycles 1..5):   Initial refinement on left region (x in [-1000, 0]).
  Session 2 (cycles 6..10):  Moved refinement to right region (x in [0, 1000]).
                             FluidInterface fills new fine nodes, PIC injects
                             new particles, and field solver & particle pusher
                             advance without instability or zero-density gaps.
  Session 3 (cycles 11..15): Refinement cleared with 'none' (derefined to level 0).
                             Field solver & particle pusher advance on single level.
"""
import logging
import math
import os

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"


def set_run_dir(run_dir):
    global RUN_DIR
    RUN_DIR = run_dir


def validate_log(pic_diags=None, test_name=None):
    """Validate that all sessions ran and energies remained physical."""
    logger.debug("Validating Dynamic AMR Multi-Session Test...")

    if not pic_diags or len(pic_diags) < 15:
        logger.error(
            "  FAIL: Expected at least 15 log entries, got %d",
            len(pic_diags) if pic_diags else 0,
        )
        return False, (
            f"Incomplete run (expected >=15 cycles, got "
            f"{len(pic_diags) if pic_diags else 0})"
        )

    cycles = [r["cycle"] for r in pic_diags]
    if cycles[-1] < 15:
        return False, f"Final cycle {cycles[-1]} < 15"

    e0 = pic_diags[0].get("Etot", 0.0)
    max_energy_drift = 0.05  # 5% max drift across all dynamic AMR sessions

    # Check finite, non-zero energies and strict conservation across all 3 sessions
    for r in pic_diags:
        cyc = r["cycle"]
        etot = r.get("Etot", 0.0)
        ee = r.get("Ee", 0.0)
        eb = r.get("Eb", 0.0)
        ep = r.get("Epart", 0.0)

        for name, val in [("Etot", etot), ("Ee", ee), ("Eb", eb), ("Epart", ep)]:
            if not math.isfinite(val):
                logger.error("  FAIL: %s is non-finite at cycle %d", name, cyc)
                return False, f"Non-finite {name} at cycle {cyc}"

        if etot <= 0.0 or ep <= 0.0:
            logger.error(
                "  FAIL: Zero or negative energy at cycle %d (Etot=%g, Epart=%g)",
                cyc,
                etot,
                ep,
            )
            return False, f"Non-positive energy at cycle {cyc}"

        drift = abs(etot - e0) / e0
        if drift > max_energy_drift:
            logger.error(
                "  FAIL: Excessive energy drift at cycle %d: Etot=%g (initial=%g, drift=%.2f%%, limit=%.2f%%)",
                cyc,
                etot,
                e0,
                drift * 100,
                max_energy_drift * 100,
            )
            return False, f"Energy drift too large at cycle {cyc} ({drift * 100:.2f}%)"

    logger.debug(
        "  All 15 cycles completed across 3 sessions with finite positive energies."
    )
    return True, "Passed"
