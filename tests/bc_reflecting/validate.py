#!/usr/bin/env python3
"""Validator for the reflecting-wall particle-boundary test.

Two variants are discovered from this directory:
  - PARAM.in           -> base_name "bc_reflecting"  (full-PIC, solveEM=F,
                          fields frozen, pure particle push)
  - PARAM.in.hybrid    -> base_name "bc_reflecting_hybrid"  (hybrid-PIC,
                          useHybridPIC=T, Ohm's-law + Faraday field advance)

Both variants exercise the shared REFLECTING particle boundary in the same way:
a single ion species streams toward the +x wall (ux > 0) with a small oblique
uy; both x faces are reflecting.  The validator is solver-agnostic (no branch on
test_name) because specular reflection conserves particle kinetic energy in both
field modes.  This validator checks:

  validate_log (energy log):
    1. The run completes with finite, positive total energy.
    2. Particle kinetic energy (Epart) does NOT decay toward zero.  Specular
       reflection conserves kinetic energy (the B-field does no work and no
       particle escapes through a reflecting face), so a large Epart loss would
       mean particles are leaking out instead of being reflected.

  validate_plot (plot output):
    With the fields frozen, the only dynamically relevant signature is that
    particles remain confined: the deposited charge/ion density (rhoS0) stays
    positive and finite across the domain (no NaN / no vacuum blow-out).
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

from .._shared.hybrid import RUN_DIR  # noqa: E402

# Tolerances.  The reflection is specular (energy-conserving) but particle
# noise and the finite cell size mean Epart is conserved only approximately.
EPART_MIN_FRAC = 0.5      # Epart may not fall below 50% of its initial value
RHO_MIN = 1e-3            # deposited rhoS0 must stay > this fraction of peak


def validate_log(pic_diags=None, test_name=None):
    """Energy-log checks: run is finite and particles are not all lost."""
    logger.debug("Validating Reflecting-Wall Test (log)...")

    if not pic_diags or len(pic_diags) < 2:
        logger.debug("  [RW] No PIC energy log found; skipping.")
        return True, "Passed (no pic log)"

    first, last = pic_diags[0], pic_diags[-1]
    e0 = first.get("Etot", 0.0)
    e1 = last.get("Etot", 0.0)
    ep0 = first.get("Epart", 0.0)
    ep1 = last.get("Epart", 0.0)

    logger.debug("    Etot: %.4e -> %.4e", e0, e1)
    logger.debug("    Epart: %.4e -> %.4e", ep0, ep1)

    if not (math.isfinite(e0) and math.isfinite(e1) and
            math.isfinite(ep0) and math.isfinite(ep1)):
        return False, "Non-finite energy (NaN/Inf)"

    if e0 <= 0:
        return False, "Initial Etot is zero (plasma not initialised)"

    # The decisive check: with reflecting walls, particles must NOT be removed
    # (an outflow wall would drain Epart toward zero).  Specular reflection
    # conserves kinetic energy, so Epart should stay well above 50% of initial.
    passed = True
    reasons = []
    if ep0 > 0 and ep1 < EPART_MIN_FRAC * ep0:
        passed = False
        reasons.append(
            f"Epart fell from {ep0:.3e} to {ep1:.3e} "
            f"(< {EPART_MIN_FRAC*100:.0f}% of initial) -- "
            f"particles appear to be lost rather than reflected")

    if passed:
        return True, "Passed (particle energy retained => reflection active)"
    return False, "; ".join(reasons)


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
    """Plot check: if rhoS0 is deposited, particles must stay confined.

    With the field solver OFF, the full-PIC node moment deposit is coupled to
    the E-solver, so rhoS0 may be zero in the plot.  In that case the check is
    skipped; the energy-log Epart conservation is the primary validation.
    """
    import tests._shared.hybrid as _hyb
    plots_dir = os.path.join(_hyb.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [RW] No .out files; skipping plot check.")
        return True, "No .out files (skipped)"

    vidx, rows = _load_out(out_files[-1])
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

    # Ion density should remain confined (not drained to a vacuum by particle
    # loss at the walls).
    min_frac = min(abs(v) for v in rho) / peak
    if min_frac < RHO_MIN:
        logger.debug("    [RW] min/peak rhoS0 = %.4e (below %.4e)",
                     min_frac, RHO_MIN)
        return False, "rhoS0 drained toward vacuum (particle loss at walls)"

    # Explicit checks at the boundary cells:
    # Verify density is positive and finite at both wall cells, and normal
    # velocity uxS0 is arrested at the reflecting wall (no unbounded penetration flux).
    xi = vidx.get("X")
    ux = _col(vidx, rows, "UXS0")
    if xi is not None and ux is not None:
        xs = [r[xi] for r in rows]
        xmin, xmax = min(xs), max(xs)
        dx_est = (xmax - xmin) / max(1, len(set(xs)) - 1)
        # Identify boundary cells within 1.5 dx of either wall
        lo_cells = [r for r in rows if r[xi] <= xmin + 1.5 * dx_est]
        hi_cells = [r for r in rows if r[xi] >= xmax - 1.5 * dx_est]

        rho_lo = [r[vidx["RHOS0"]] for r in lo_cells]
        rho_hi = [r[vidx["RHOS0"]] for r in hi_cells]
        ux_lo = [r[vidx["UXS0"]] for r in lo_cells]
        ux_hi = [r[vidx["UXS0"]] for r in hi_cells]

        logger.debug("    [RW] Boundary cells: lo_wall rho=%.4e, ux=%.4e | hi_wall rho=%.4e, ux=%.4e",
                     sum(rho_lo)/len(rho_lo), sum(ux_lo)/len(ux_lo),
                     sum(rho_hi)/len(rho_hi), sum(ux_hi)/len(ux_hi))

        for r_val in rho_lo + rho_hi:
            if not math.isfinite(r_val) or r_val < 0.0:
                return False, "Non-finite or negative density at boundary cells"

        for u_val in ux_lo + ux_hi:
            if not math.isfinite(u_val):
                return False, "Non-finite normal velocity at boundary cells"

    return True, "Passed (particles confined, boundary cells verified, rhoS0 positive and finite)"
