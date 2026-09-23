#!/usr/bin/env python3
"""Validator for the reflecting / conducting wall boundary tests (tests/bc_reflecting/).

Four variants are discovered from this directory:
  - PARAM.in.fields          -> base_name "bc_reflecting_fields"
      Full-PIC conducting field wall (tophat EM pulse bounces off PEC walls,
      conserving total electromagnetic energy).
  - PARAM.in.hybrid.fields   -> base_name "bc_reflecting_hybrid_fields"
      Hybrid-PIC conducting field wall (shear-Alfvén wave in PEC cavity;
      finite and stable field advance).
  - PARAM.in.particles       -> base_name "bc_reflecting_particles"
      Full-PIC specular particle reflection (pure particle push, solveEM=F).
  - PARAM.in.hybrid.particles-> base_name "bc_reflecting_hybrid_particles"
      Hybrid-PIC specular particle reflection (hybrid mover, useHybridPIC=T).

This validator branches on test_name to verify the appropriate physics.
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

import tests._shared.hybrid as _hyb

RUN_DIR = "run_test"


def set_run_dir(run_dir):
    """Mirror the runner's RUN_DIR into the shared hybrid helper."""
    global RUN_DIR
    RUN_DIR = run_dir
    _hyb.set_run_dir(run_dir)


# Tolerances
EPART_MIN_FRAC = 0.5       # particles: Epart may not fall below 50% of initial
RHO_MIN = 1e-3             # particles: deposited rhoS0 must stay > this fraction of peak
ETOT_CONDUCTING_MIN = 0.8  # full-PIC fields: Etot_final / Etot_initial >= 0.8 (conserved)
ETOT_CONDUCTING_MAX = 1.25 # full-PIC fields: Etot_final / Etot_initial <= 1.25 (conserved)
ETOT_HYBRID_GROWTH = 1e3   # hybrid fields: Etot may not grow by more than 1e3x (stable)


def validate_log(pic_diags=None, test_name=None):
    """Energy-log checks, branching on the variant."""
    logger.debug("Validating %s (log)...", test_name)

    if not pic_diags or len(pic_diags) < 2:
        logger.debug("  No PIC energy log found; skipping.")
        return True, "Passed (no pic log)"

    first, last = pic_diags[0], pic_diags[-1]
    e0 = first.get("Etot", 0.0)
    e1 = last.get("Etot", 0.0)

    if test_name == "bc_reflecting_fields":
        return _validate_log_fields(e0, e1, first, last)
    if test_name == "bc_reflecting_hybrid_fields":
        return _validate_log_hybrid_fields(pic_diags)

    # Particle reflection (bc_reflecting_particles, bc_reflecting_hybrid_particles,
    # or legacy unadorned names)
    return _validate_log_particles(e0, e1, first, last)


def _validate_log_particles(e0, e1, first, last):
    """Particles: specular reflection conserves kinetic energy (Epart does not drain)."""
    ep0 = first.get("Epart", 0.0)
    ep1 = last.get("Epart", 0.0)
    logger.debug("    Etot: %.4e -> %.4e", e0, e1)
    logger.debug("    Epart: %.4e -> %.4e", ep0, ep1)

    if not (math.isfinite(e0) and math.isfinite(e1) and
            math.isfinite(ep0) and math.isfinite(ep1)):
        return False, "Non-finite energy (NaN/Inf)"

    if e0 <= 0:
        return False, "Initial Etot is zero (plasma not initialised)"

    if ep0 > 0 and ep1 < EPART_MIN_FRAC * ep0:
        return False, (
            f"Epart fell from {ep0:.3e} to {ep1:.3e} "
            f"(< {EPART_MIN_FRAC*100:.0f}% of initial) -- "
            f"particles appear to be lost rather than reflected"
        )

    return True, "Passed (particle energy retained => reflection active)"


def _validate_log_fields(e0, e1, first, last):
    """Full-PIC conducting field wall: EM energy is conserved as pulses reflect."""
    e0_f = first.get("Ee", 0.0) + first.get("Eb", 0.0)
    e1_f = last.get("Ee", 0.0) + last.get("Eb", 0.0)
    logger.debug("    Ee+Eb: %.4e -> %.4e", e0_f, e1_f)

    for v in (e0, e1, e0_f, e1_f):
        if not math.isfinite(v):
            return False, "Non-finite EM energy (NaN/Inf)"
    if e0 <= 0:
        return False, "Initial Etot is zero (pulse not seeded)"

    denom = e0 if e0 > 0 else 1.0
    ratio = e1 / denom
    logger.debug("    Etot_final/Etot_initial = %.4f", ratio)

    if ratio < ETOT_CONDUCTING_MIN:
        return False, (
            f"EM energy lost (final/initial = {ratio:.3f} < {ETOT_CONDUCTING_MIN}) "
            f"-- pulses appear to be absorbed or leaking instead of reflecting"
        )
    if ratio > ETOT_CONDUCTING_MAX:
        return False, (
            f"EM energy grew excessively (final/initial = {ratio:.3f} > {ETOT_CONDUCTING_MAX}) "
            f"-- conducting wall instability"
        )

    return True, "Passed (EM energy conserved across PEC reflections)"


def _validate_log_hybrid_fields(pic_diags):
    """Hybrid-PIC conducting wall: advance stays finite, positive, and bounded."""
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
    if e0 > 0 and e1 > ETOT_HYBRID_GROWTH * e0:
        return False, (
            f"Etot grew from {e0:.3e} to {e1:.3e} "
            f"(>{ETOT_HYBRID_GROWTH:.0f}x) -- field-wall instability"
        )

    return True, "Passed (finite, bounded energy => PEC field wall is stable)"


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
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


def _find_boundary_rows(vidx, rows, n_cells=1):
    """Return (lo_rows, hi_rows) -- the rows nearest the low and high x boundaries.

    Uses the 'X' column when present to find the actual min/max-x rows, so the
    check is correct regardless of whether AMReX writes rows in sorted order.
    Falls back to rows[0:n_cells] / rows[-n_cells:] when no X column exists.

    Parameters
    ----------
    n_cells : int
        How many boundary-adjacent rows to include on each side.
    """
    xi = vidx.get("X")
    if xi is None or not rows:
        # No X column -- fall back to positional indexing.
        return rows[:n_cells], rows[-n_cells:]

    xs = [r[xi] for r in rows]
    xmin, xmax = min(xs), max(xs)
    dx_est = (xmax - xmin) / max(1, len(set(xs)) - 1) if len(set(xs)) > 1 else 1.0
    half = max(0.5, n_cells) * dx_est

    lo_rows = [r for r in rows if r[xi] <= xmin + half]
    hi_rows = [r for r in rows if r[xi] >= xmax - half]
    if not lo_rows:
        lo_rows = rows[:n_cells]
    if not hi_rows:
        hi_rows = rows[-n_cells:]
    return lo_rows, hi_rows


def validate_plot(test_name):
    """Plot checks branching on variant."""
    logger.debug("Validating %s (plot)...", test_name)
    run_dir = getattr(_hyb, "RUN_DIR", RUN_DIR)
    plots_dir = os.path.join(run_dir, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [RW] No .out files; skipping plot check.")
        return True, "No .out files (skipped)"

    if test_name == "bc_reflecting_fields":
        return _check_fields_plot(out_files)
    if test_name == "bc_reflecting_hybrid_fields":
        return _check_hybrid_fields_plot(out_files)

    return _check_particles_plot(out_files[-1])


def _check_fields_plot(out_files):
    """Full-PIC conducting fields: verify PEC reflection via finiteness,
    boundedness, and interior pulse detection.

    Note: the .out file writes volume-sampled values (not wall-face values), so
    we cannot directly check Et=0 / Bn=0 from plot output -- the PEC BC
    constrains the *face* value, not the adjacent output point.  Instead we rely
    on:
      * All field components finite (no NaN/Inf).
      * No blow-up (|field| < 1e6).
      * The pulse is still present in the interior (not absorbed by the wall).
    Energy conservation across reflections is verified in validate_log
    (_validate_log_fields), which robustly catches the anti-symmetric-reflection
    bug (absorbing wall would drain EM energy below the ETOT_CONDUCTING_MIN
    threshold).
    """
    if isinstance(out_files, str):
        out_files = [out_files]

    has_interior_pulse = False
    for out_file in out_files:
        vidx, rows = _load_out(out_file)
        if vidx is None or not rows:
            continue

        for var in ("EY", "BZ", "EX", "EZ", "BX", "BY"):
            col = _col(vidx, rows, var)
            if col is None:
                continue
            if not all(math.isfinite(v) for v in col):
                return False, f"{var} not finite (NaN)"
            peak = max(abs(v) for v in col)
            if peak > 1e6:
                return False, f"{var} blew up (peak {peak:.2e})"

        col_ey = _col(vidx, rows, "EY")
        if col_ey and max(abs(v) for v in col_ey) > 0.05:
            has_interior_pulse = True

    if not has_interior_pulse:
        return False, "Pulse not detected in interior (EM fields stayed near 0)"

    return True, "Passed (fields finite, bounded, and interior pulse persists => PEC reflection active)"


# Soft tolerance for near-wall |Et| in the hybrid pulse test.  The cell-
# centred hybrid advance applies the PEC BC on the wall face; the nearest
# plotted cell (at the cell centre, ~dx/2 from the face) can carry a small
# residual Ey/Ez from the Ohm's-law evaluation before the BC overwrites it.
# NOTE: we only check tangential E (Ey, Ez) here, NOT normal B (Bx), because
# the test carries a background guide field Bx ≈ 1.0 throughout the domain.
_HYBRID_WALL_ET_TOL = 0.5

# Minimum fraction of the initial pulse peak |By| that must still be present
# at the late frame.  A fully absorbing wall drains all By energy out of the
# domain; a properly reflecting wall conserves it.  The threshold 0.3 is loose
# enough to tolerate phase-mixing and numerical dispersion over many bounces
# while still failing if the wall is absorbing most of the energy.
_HYBRID_PULSE_AMPLITUDE_RATIO_MIN = 0.3


def _check_hybrid_fields_plot(out_files):
    """Hybrid-PIC Gaussian-pulse conducting-wall test.

    Checks:
      * All field components finite (no NaN/Inf) at every frame.
      * No blow-up (|field| < 1e6) at every frame.
      * Guide field Bx stays in (0.5, 2.0) in interior rows.
      * Near-wall tangential E (Ey, Ez) does not spike at the x walls.
      * Pulse-energy conservation: max|By| at the last frame is at least
        _HYBRID_PULSE_AMPLITUDE_RATIO_MIN × max|By| at the first frame.
        A properly reflecting wall conserves the Alfven-pulse amplitude;
        an absorbing wall drains it below the threshold.

    NOTE on ghost-cell sensitivity: the wall-node Et=0 BC is enforced by
    section 1 of apply_conducting_wall (correct in both old and new code),
    so the 1D single-process pulse test validates section 1.  Section 2
    (ghost-cell anti-symmetry for MPI halos / AMR boundaries) requires a
    multi-level or multi-MPI test to distinguish.
    """
    if isinstance(out_files, str):
        out_files = [out_files]

    # Sort frames by filename (they embed the step number).
    sorted_files = sorted(out_files)
    if not sorted_files:
        return True, "No .out files to check (skipped)"

    early_peak_by = None
    late_peak_by = None

    for frame_idx, out_file in enumerate(sorted_files):
        vidx, rows = _load_out(out_file)
        if vidx is None or not rows:
            continue

        # Global finiteness / blow-up check.
        for var in ("BX", "BY", "BZ", "EX", "EY", "EZ"):
            col = _col(vidx, rows, var)
            if col is None:
                continue
            if not all(math.isfinite(v) for v in col):
                return False, f"[frame {frame_idx}] {var} not finite (NaN)"
            peak = max(abs(v) for v in col)
            if peak > 1e6:
                return False, f"[frame {frame_idx}] {var} blew up (peak {peak:.2e})"

        # Near-wall tangential-E check (Ey, Ez only; Bx has guide-field offset).
        lo_rows, hi_rows = _find_boundary_rows(vidx, rows, n_cells=1)
        for side, brows in (("low", lo_rows), ("high", hi_rows)):
            for var in ("EY", "EZ"):
                idx = vidx.get(var)
                if idx is None:
                    continue
                worst = max(abs(r[idx]) for r in brows)
                if worst > _HYBRID_WALL_ET_TOL:
                    return False, (
                        f"[frame {frame_idx}] Near-wall {var} spike at {side}-x "
                        f"wall: {worst:.3e} > {_HYBRID_WALL_ET_TOL:.2g}"
                    )

        # Guide field Bx: interior rows only (exclude potential ghost effects).
        boundary_set = set(id(r) for r in lo_rows + hi_rows)
        interior_rows = [r for r in rows if id(r) not in boundary_set] or rows
        col_bx_int = _col(vidx, interior_rows, "BX") if interior_rows else None
        if col_bx_int:
            min_bx = min(col_bx_int)
            max_bx = max(col_bx_int)
            if min_bx < 0.5 or max_bx > 2.0:
                return False, (
                    f"[frame {frame_idx}] Guide field Bx abnormal in interior: "
                    f"min={min_bx:.3f}, max={max_bx:.3f} (expected ~1.0)"
                )

        # Track peak |By| for pulse-amplitude conservation check.
        col_by = _col(vidx, rows, "BY")
        if col_by:
            peak_by = max(abs(v) for v in col_by)
            if early_peak_by is None:
                early_peak_by = peak_by
            late_peak_by = peak_by

    # Pulse-amplitude conservation: the reflecting wall must preserve the By peak.
    if early_peak_by is not None and late_peak_by is not None:
        if early_peak_by < 1e-12:
            return False, "Initial By pulse is essentially zero (not seeded?)"
        ratio = late_peak_by / early_peak_by
        logger.debug("    [HYB] pulse amplitude ratio late/early |By| = %.3f", ratio)
        if ratio < _HYBRID_PULSE_AMPLITUDE_RATIO_MIN:
            return False, (
                f"Alfven pulse was absorbed (not reflected) by the conducting "
                f"wall: |By| ratio = {ratio:.3f} < {_HYBRID_PULSE_AMPLITUDE_RATIO_MIN}"
            )

    return True, (
        f"Passed (hybrid pulse reflected: amplitude ratio "
        f"{late_peak_by / early_peak_by:.2f}, "
        f"fields finite, near-wall |Et| OK, guide field stable)"
    )



def _check_particles_plot(out_file):
    """Particles: verify density remains confined and positive at wall cells."""
    vidx, rows = _load_out(out_file)
    if vidx is None or not rows:
        return True, "Could not parse .out (skipped)"

    rho = _col(vidx, rows, "RHOS0")
    if rho is None:
        return True, "No rhoS0 column (skipped)"

    peak = max(abs(v) for v in rho)
    logger.debug("    [RW] rhoS0 peak = %.4e", peak)

    # No-solver mode does not deposit output moments (zero array): skip.
    if peak <= 0:
        return True, "rhoS0 not deposited in no-solver mode (skipped)"

    finite = all(math.isfinite(v) for v in rho)
    if not finite:
        return False, "rhoS0 not finite (NaN)"

    # Ion density should remain confined (not drained to vacuum)
    min_frac = min(abs(v) for v in rho) / peak
    if min_frac < RHO_MIN:
        logger.debug("    [RW] min/peak rhoS0 = %.4e (below %.4e)",
                     min_frac, RHO_MIN)
        return False, "rhoS0 drained toward vacuum (particle loss at walls)"

    # Explicit boundary cell verification
    xi = vidx.get("X")
    ux = _col(vidx, rows, "UXS0")
    if xi is not None and ux is not None:
        xs = [r[xi] for r in rows]
        xmin, xmax = min(xs), max(xs)
        dx_est = (xmax - xmin) / max(1, len(set(xs)) - 1)
        lo_cells = [r for r in rows if r[xi] <= xmin + 1.5 * dx_est]
        hi_cells = [r for r in rows if r[xi] >= xmax - 1.5 * dx_est]

        rho_lo = [r[vidx["RHOS0"]] for r in lo_cells]
        rho_hi = [r[vidx["RHOS0"]] for r in hi_cells]
        ux_lo = [r[vidx["UXS0"]] for r in lo_cells]
        ux_hi = [r[vidx["UXS0"]] for r in hi_cells]

        for r_val in rho_lo + rho_hi:
            if not math.isfinite(r_val) or r_val < 0.0:
                return False, "Non-finite or negative density at boundary cells"

        for u_val in ux_lo + ux_hi:
            if not math.isfinite(u_val):
                return False, "Non-finite normal velocity at boundary cells"

    return True, "Passed (particles confined, boundary cells verified, rhoS0 positive and finite)"
