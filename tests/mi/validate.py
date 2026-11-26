#!/usr/bin/env python3
"""Validator for the 3D Hybrid-PIC MI coupling test (tests/mi)."""
import logging
import math
import os
import re

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"


def set_run_dir(run_dir):
    global RUN_DIR
    RUN_DIR = run_dir


def _non_finite(diags):
    keys = sorted({k for e in (diags or []) for k, v in e.items()
                   if isinstance(v, float) and not math.isfinite(v)})
    if keys:
        return f"non-finite log values in {', '.join(keys[:6])}"
    return None


def validate_log(pic_diags=None, test_name=None):
    bad = _non_finite(pic_diags)
    if bad:
        return False, bad
    if not pic_diags or len(pic_diags) < 2:
        return False, "missing or insufficient PIC log rows"
    return True, f"Passed ({len(pic_diags)} log steps)"


def validate_plot(test_name=None):
    diag_path = os.path.join(RUN_DIR, "PC", "plots", "ionosphere_n0000020.dat")
    if not os.path.isfile(diag_path):
        return False, f"missing MI diagnostic file {diag_path}"

    hdr = {}
    rows = []
    with open(diag_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                for m in re.finditer(r"(\w+)\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", s):
                    hdr[m.group(1)] = float(m.group(2))
                continue
            vals = [float(tok) for tok in s.split()]
            if len(vals) < 6 or any(not math.isfinite(v) for v in vals):
                return False, "non-finite or malformed row in ionosphere diagnostic output"
            rows.append(vals)

    n_theta = int(hdr.get("nTheta", -1))
    n_phi = int(hdr.get("nPhi", -1))
    if n_theta <= 0 or n_phi <= 0 or len(rows) != n_theta * n_phi:
        return False, f"unexpected grid count {len(rows)} vs {n_theta}x{n_phi}"

    cpcp_n = hdr.get("cpcpN_kV", float("nan"))
    cpcp_s = hdr.get("cpcpS_kV", float("nan"))
    if not (math.isfinite(cpcp_n) and math.isfinite(cpcp_s)):
        return False, "non-finite CPCP header values"

    reasons = []
    if not (0.01 < cpcp_n < 1e4) or not (0.01 < cpcp_s < 1e4):
        reasons.append(f"CPCP out of physical range: N={cpcp_n:.4f} kV, S={cpcp_s:.4f} kV")
    sym = abs(cpcp_n - cpcp_s) / max(cpcp_n, cpcp_s, 1e-30)
    if sym >= 0.35:
        reasons.append(f"north/south CPCP asymmetry {sym:.3f} >= 0.35 "
                       f"(N={cpcp_n:.4f}, S={cpcp_s:.4f})")

    max_jr = max(abs(r[3]) for r in rows)
    max_sp = max(r[4] for r in rows)
    if max_jr <= 1e-4:
        reasons.append(f"max(|Jr|)={max_jr:.3e} uA/m2 too small")
    if max_sp <= 3.01:
        reasons.append(f"max(SigmaP)={max_sp:.3f} S shows no auroral/EUV enhancement above 3.0 S")

    if reasons:
        return False, "; ".join(reasons)

    return True, (f"Passed (CPCP_N={cpcp_n:.3f} kV, CPCP_S={cpcp_s:.3f} kV, "
                  f"max|Jr|={max_jr:.3e} uA/m2, maxSigmaP={max_sp:.2f} S)")
