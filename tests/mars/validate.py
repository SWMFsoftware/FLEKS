#!/usr/bin/env python3
"""Validator for the 3D Mars test (tests/mars): the full boundary setup
(run_mars/final_coarse) over 10 s, see _validate_emission_log.
"""
import logging
import math

logger = logging.getLogger(__name__)

MASS = [0.0625, 1.0, 32.0, 16.0, 44.0]


def _non_finite(diags):
    """Reason string if a log value is NaN/inf (comparisons with NaN are false,
    so a tolerance check alone would pass a broken run), else None."""
    keys = sorted({k for e in (diags or []) for k, v in e.items()
                   if isinstance(v, float) and not math.isfinite(v)})
    if keys:
        return f"non-finite log values in {', '.join(keys[:6])}"
    return None


def validate_log(pic_diags=None, test_name=None):
    bad = _non_finite(pic_diags)
    if bad:
        return False, bad
    return _validate_emission_log(pic_diags)


# ---------------------------------------------------------------------------
# fullpic_emission: the final boundary setup (run_mars/final_coarse), 10 s
# ---------------------------------------------------------------------------
FACES = ["xm", "xp", "ym", "yp", "zm", "zp"]
O2_RATE = 9.0  # emitted O2+, 2e25 /s, in units of the O2+ particle count


def _validate_emission_log(d):
    """Log checks of the emission-shell / floating-body / open-face setup.

    * quasi-neutral: max |net charge| < 0.5% of N(H+) (measured 0.13%);
    * floating body: the re-emitted charge equals the sum of the per-step
      electron excesses of the absorbed charge (steps in which the ions
      dominate re-emit nothing), to round-off;
    * H+ budget: dM - (face in - face out - body) is the net H+ source, the
      corona (+1.2 code mass/s), constant; checked in [0.8, 1.6];
    * the planetary ions accumulate: emitted O2+ minus what fell back;
    * the field and particle energies stay bounded.

    The bow shock (1.75 R_M) and the pileup (30 - 40 nT) form only after
    30 - 40 s; the 200 s run_mars/final_coarse checks them, not this test.
    """
    if not d or len(d) < 100 or "mIn_xm_S1" not in d[0]:
        return False, "missing log or face budget columns"
    reasons = []

    def check(cond, msg):
        logger.debug("    %s: %s", "OK  " if cond else "FAIL", msg)
        if not cond:
            reasons.append(msg)

    n = lambda e, i: e[f"Mpart{i}"] / MASS[i]
    first, last = d[0], d[-1]

    worstQ = max(abs(sum(n(e, i) for i in range(1, 5)) - n(e, 0)) / n(e, 1)
                 for e in d)
    check(worstQ < 5e-3,
          f"quasi-neutral: max |net charge|/N(H+) = {worstQ:.2e} (< 5e-3)")

    excess = sum(max(0.0, -sum(b[f"qBodyS{i}"] - a[f"qBodyS{i}"]
                               for i in range(5)))
                 for a, b in zip(d, d[1:]))
    check(abs(last["qBodyReemit"] - excess) < 1e-6 * excess,
          f"re-emitted = absorbed electron excess ({last['qBodyReemit']:.6g}"
          f" vs {excess:.6g})")

    # H+ budget over the second half (after the start-up absorption).
    a = next(e for e in d if e["time"] >= last["time"] / 2)
    dt = last["time"] - a["time"]
    faces = sum((last[f"mIn_{f}_S1"] - a[f"mIn_{f}_S1"]) -
                (last[f"mOut_{f}_S1"] - a[f"mOut_{f}_S1"]) for f in FACES)
    body = last["qBodyS1"] - a["qBodyS1"]
    src = (last["Mpart1"] - a["Mpart1"] - faces + body) / dt
    check(0.8 < src < 1.6,
          f"H+ net source (corona) = {src:.3f} code mass/s (0.8 - 1.6)")

    o2 = n(last, 2) - n(first, 2)
    emitted = O2_RATE * (last["time"] - first["time"])
    check(0.5 * emitted < o2 < 1.02 * emitted,
          f"O2+ in the plasma = {o2:.1f} (emitted {emitted:.0f}, less what "
          "fell back)")

    rEb = last["Eb"] / first["Eb"]
    check(0.8 < rEb < 6, f"Eb/Eb0 = {rEb:.2f} (0.8 - 6)")
    rEe = max(e["Ee"] for e in d) / first["Ee"]
    check(rEe < 20, f"max Ee/Ee0 = {rEe:.1f} (< 20)")
    rEp = last["Epart"] / first["Epart"]
    check(0.95 < rEp < 1.05, f"Epart/Epart0 = {rEp:.3f} (0.95 - 1.05)")

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"
