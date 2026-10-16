#!/usr/bin/env python3
"""Validator for the ionospheric emission shell (tests/ionosphereemission).

H+ (1e16 /s) and O+ (2e15 /s) are emitted at r = 2 above an absorbing body
(R = 1.5) into vacuum, with E = B = 0. In code units one proton has the mass
m1 = e^2 / (1e7 lNormSI m_p) (FluidInterface: Mnorm = 1e3 mNormSI [g]), so

* the injected mass grows at (1e16 * 1 + 2e15 * 16) m1 per second (1e-6);
* the emission is outward, so the body absorbs nothing (< 1e-3);
* within 2 s no ion reaches the walls: the mass of each species in the box
  grows at its own rate times its mass (1e-3; the vacuum seed is negligible).
"""
import logging

logger = logging.getLogger(__name__)

E_SI = 1.602176634e-19
MP_SI = 1.67262192369e-27
L_NORM = 1.0e5
M1 = E_SI**2 / (1.0e7 * L_NORM * MP_SI)
RATES = [1.0e16, 2.0e15]
MASSES = [1.0, 16.0]


def validate_log(pic_diags=None, test_name=None):
    d = pic_diags
    if not d or len(d) < 10 or "mBodyInject" not in d[0]:
        return False, "missing log or mBodyInject column"

    reasons = []

    def check(cond, msg):
        logger.debug("    %s: %s", "OK  " if cond else "FAIL", msg)
        if not cond:
            reasons.append(msg)

    a, b = d[1], d[-1]
    dt = b["time"] - a["time"]
    expected = sum(r * m for r, m in zip(RATES, MASSES)) * M1
    rate = (b["mBodyInject"] - a["mBodyInject"]) / dt
    check(abs(rate / expected - 1) < 1e-6,
          f"injected mass rate = {rate / expected:.8f} x expected (tol 1e-6)")

    absorbed = b["mBodyAbsorb"] / max(b["mBodyInject"], 1e-300)
    check(absorbed < 1e-3,
          f"absorbed / injected mass = {absorbed:.2e} (outward emission, "
          f"tol 1e-3)")

    for i, (r, m) in enumerate(zip(RATES, MASSES)):
        grow = (b[f"Mpart{i}"] - a[f"Mpart{i}"]) / dt / (r * m * M1)
        check(abs(grow - 1) < 1e-3,
              f"species {i}: box mass grows at {grow:.6f} x its emission rate "
              f"(tol 1e-3)")

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"
