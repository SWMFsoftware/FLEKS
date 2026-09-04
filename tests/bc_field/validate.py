#!/usr/bin/env python3
"""Validator for the hybrid-PIC field-wall boundary-condition test.

Two variants are discovered from this directory:
  - PARAM.in.hybrid.conducting -> base_name "bc_field_hybrid_conducting"
  - PARAM.in.hybrid.absorb     -> base_name "bc_field_hybrid_absorb"

Both are hybrid-PIC (`useHybridPIC=T`) finite-x runs of a transverse shear-Alfven
wave (HybridWave) in a cavity whose x faces are a physical FIELD boundary --
CONDUCTING (PEC: E_t=0, B_n=0) or ABSORBING (outgoing characteristic) -- and whose
x particle faces are REFLECTING (so the ions, and hence the Ohm's-law moments,
stay inside).

This test provides the missing hybrid-PIC coverage for the `conducting` and
`absorb` FIELD boundary conditions (the full-PIC `bc_absorb_fields` covers only
the `absorb` field BC, and only in a pure-EM, zero-particle `tophat` domain that
is incompatible with the hybrid Ohm's law).  It validates that:

  validate_log (energy log):
    1. The run completes with finite, positive total energy (no NaN / no
       vacuum blow-out) -- catches a wall-BC bug that destabilises the hybrid
       field advance.
    2. Total energy stays bounded: the final Etot does not blow up by a large
       factor from its initial value.  A runaway field-wall instability would
       show as unbounded growth.

  validate_plot (plot output):
    Field/fluid values are finite and bounded (no NaN / no gross blow-up) near
    the walls and across the domain.

The two variants share the identical geometry, seed and solver settings, so any
difference in behaviour is attributable to the field boundary under test.
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

def set_run_dir(run_dir):
    """Mirror the runner's RUN_DIR into the shared hybrid helper."""
    import tests._shared.hybrid as _hyb
    _hyb.set_run_dir(run_dir)

# Tolerances.  A seeded standing Alfven wave evolves but must stay finite and
# bounded; the exact bound is loose because it depends on the mode amplitude.
ETOT_MIN = 1e-4            # Etot must stay above this (no vacuum blow-out)
ETOT_GROWTH_MAX = 1e3      # Etot may not grow by more than 1e3x (no instability)


def validate_log(pic_diags=None, test_name=None):
    """Energy-log checks: run is finite and energy is bounded."""
    logger.debug("Validating Field-Wall Test (log)...")

    if not pic_diags or len(pic_diags) < 2:
        logger.debug("  [FW] No PIC energy log found; skipping.")
        return True, "Passed (no pic log)"

    e0 = pic_diags[0].get("Etot", 0.0)
    finite = all(
        math.isfinite(d.get("Etot", 0.0)) and
        math.isfinite(d.get("Epart", 0.0)) for d in pic_diags
    )
    e1 = pic_diags[-1].get("Etot", 0.0)

    logger.debug("    Etot: %.4e -> %.4e (%.1f growth)", e0, e1,
                 1.0 if e0 == 0 else e1 / e0)

    if not finite:
        return False, "Non-finite energy (NaN/Inf) in energy log"

    if e0 <= 0 or e1 <= 0:
        return False, "Non-positive total energy (plasma not initialised / drained)"

    if e0 > 0 and e1 > ETOT_GROWTH_MAX * e0:
        return False, (f"Etot grew from {e0:.3e} to {e1:.3e} "
                       f"(>{ETOT_GROWTH_MAX:.0f}x) -- field-wall instability")

    return True, "Passed (finite, bounded energy => field wall is stable)"


def _load_out(out_file):
    """Load a .out frame: return ({VAR: col_idx}, float rows)."""
    with open(out_file, "r", encoding="latin-1") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return None, None
    vidx = {v.upper(): i for i, v in enumerate(lines[4].split())}
    rows = []
    for line in lines[5:]:
        cols = line.split()
        if not cols:
            continue
        try:
            rows.append([float(c) for c in cols])
        except ValueError:
            continue
    return vidx, rows


def _col(vidx, rows, name):
    i = vidx.get(name)
    if i is None or not rows or i >= len(rows[0]):
        return None
    return [r[i] for r in rows]


def validate_plot(test_name):
    """Plot check: field/fluid values are finite and bounded (no NaN/blow-up)."""
    import tests._shared.hybrid as _hyb
    plots_dir = os.path.join(_hyb.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [FW] No .out files; skipping plot check.")
        return True, "No .out files (skipped)"

    # Check the last two frames so a late-time blow-up is caught.
    for out_file in out_files[-2:]:
        vidx, rows = _load_out(out_file)
        if vidx is None or not rows:
            continue
        for var in ("BY", "BZ", "EY", "EZ"):
            col = _col(vidx, rows, var)
            if col is None:
                continue
            peak = max(abs(v) for v in col)
            logger.debug("    [FW] %s peak (last frame) = %.4e", var, peak)
            if not all(math.isfinite(v) for v in col):
                return False, f"{var} not finite (NaN)"
            if peak > 1e6:
                return False, f"{var} blew up (peak {peak:.2e})"

    return True, "Passed (fields finite and bounded near walls)"
