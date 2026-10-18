#!/usr/bin/env python3
"""Validator for the open particle boundary (tests/openboundary).

A uniform plasma of H+ and O+ (T = 1e5 K, thermal speeds 28.7 and 7.2 km/s)
fills a box whose six faces are outflowthermal, without fields. Each face
re-injects the inward half-flux n vtherm g(vd) of the local Maxwellian of each
species, with the normal drift limited to its outward part (the re-entered
particles bias the measured drift inward, which would feed back), so each
species keeps its mass: at rest (rest), and drifting along +x at 2 sigma(H+)
= 8 sigma(O+) (drift), entering through an inflow face (#INFLOW) at x- and
leaving through the others; the re-entry at x+ is reduced by the outward
drift, and the flanks see a flow parallel to them. (The inflow face injects
the convective flux n u only, so a subsonic drift would lose the thermal
back-streaming through it.) With plain outflow H+ would lose about half of
its mass in the 10 s of the run.

Checks over the whole run, both species: Mpart_s(t)/Mpart_s(0) within 3%,
and (rest only) the temperature, Epart_s/Mpart_s relative to its initial
value, within 3%: the re-entered particles must carry the energy flux 2 kT
per particle (a flux-weighted speed sampler that over-sampled its tail heated
H+ by 23%; a clamp of the noisy single-step drift cooled it by 7%). In the
drift variant Epart includes the bulk flow, and the inflow face, which
injects |u + sigma g| at the convective rate, lowers it by 17% with plain
outflow at x+ as well, so the energy is not checked there.

Budget: without body or sources the mass change of each species equals the
mass injected minus removed at the faces (mIn_/mOut_ log columns), 1e-9.
drift_pc runs drift with the predictor-corrector step (rk4-pc), whose trial
push is discarded: its fluxes must not stay in the tallies. periodic is rest
with a periodic x: nothing may be booked to the x faces.
"""
import logging
import math

logger = logging.getLogger(__name__)


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
    d = pic_diags
    if not d or len(d) < 50 or "Mpart1" not in d[0]:
        return False, "missing log or Mpart columns"
    reasons = []
    for iS in (0, 1):
        ratios = [e[f"Mpart{iS}"] / d[0][f"Mpart{iS}"] for e in d]
        worst = max(ratios, key=lambda r: abs(r - 1))
        logger.debug("    species %d: mass/initial in [%.4f, %.4f], end %.4f",
                     iS, min(ratios), max(ratios), ratios[-1])
        if abs(worst - 1) > 0.03:
            reasons.append(f"species {iS}: mass/initial = {worst:.4f} "
                           f"(tol 3%)")
    isDrift = "drift" in (test_name or "")
    for iS in ((0, 1) if not isDrift else ()):
        e0 = d[0][f"Epart{iS}"] / d[0][f"Mpart{iS}"]
        ratios = [e[f"Epart{iS}"] / e[f"Mpart{iS}"] / e0 for e in d]
        worst = max(ratios, key=lambda r: abs(r - 1))
        logger.debug("    species %d: T/T0 in [%.4f, %.4f], end %.4f",
                     iS, min(ratios), max(ratios), ratios[-1])
        if abs(worst - 1) > 0.03:
            reasons.append(f"species {iS}: T/T0 = {worst:.4f} (tol 3%)")
    faces = ["xm", "xp", "ym", "yp", "zm", "zp"]
    if "mIn_xm_S0" in d[0]:
        for iS in (0, 1):
            net = sum(d[-1][f"mIn_{f}_S{iS}"] - d[-1][f"mOut_{f}_S{iS}"]
                      for f in faces)
            dM = d[-1][f"Mpart{iS}"] - d[0][f"Mpart{iS}"]
            inflow = sum(d[-1][f"mIn_{f}_S{iS}"] for f in faces)
            err = abs(dM - net) / inflow
            logger.debug("    species %d: budget residual %.2e", iS, err)
            if err > 1e-9:
                reasons.append(f"species {iS}: dM - (in - out) = {err:.2e} "
                               f"of the injected mass (tol 1e-9)")
        if (test_name or "").endswith("periodic"):
            # A particle that wraps through x is no loss; one that leaves a y
            # or z face across a box edge must be booked there, not at x.
            xFlux = sum(d[-1][f"m{io}_{f}_S{iS}"] for io in ("In", "Out")
                        for f in ("xm", "xp") for iS in (0, 1))
            if xFlux != 0:
                reasons.append(f"periodic x faces carry a flux ({xFlux:.3g})")
    else:
        reasons.append("missing mIn_/mOut_ face budget columns")
    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"
