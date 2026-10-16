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
"""
import logging

logger = logging.getLogger(__name__)


def validate_log(pic_diags=None, test_name=None):
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
    for iS in ((0, 1) if not (test_name or "").endswith("drift") else ()):
        e0 = d[0][f"Epart{iS}"] / d[0][f"Mpart{iS}"]
        ratios = [e[f"Epart{iS}"] / e[f"Mpart{iS}"] / e0 for e in d]
        worst = max(ratios, key=lambda r: abs(r - 1))
        logger.debug("    species %d: T/T0 in [%.4f, %.4f], end %.4f",
                     iS, min(ratios), max(ratios), ratios[-1])
        if abs(worst - 1) > 0.03:
            reasons.append(f"species {iS}: T/T0 = {worst:.4f} (tol 3%)")
    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"
