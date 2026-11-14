#!/usr/bin/env python3
"""Validation for the solver_coupling test suite.

Covers five coupled multi-solver variants:
  * solver_coupling                  : Implicit (x < 0) + RK4 Hybrid w/ Pe (x >= 0)
  * solver_coupling_explicit_hybrid  : Explicit PC (x < 0) + RK4 Hybrid (x >= 0)
  * solver_coupling_implicit_explicit: Implicit (x < 0) + Explicit Leapfrog (x >= 0)
  * solver_coupling_blended_amr      : 2D AMR straddling seam + nBlendCell = 2
  * solver_coupling_three_way        : Implicit + Explicit Leapfrog + RK4 Hybrid w/ Pe
"""

import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"

# Solver codes matching PicSolverMode in include/Grid.h
SOLVER_IMPLICIT = 0
SOLVER_EXPLICIT = 1
SOLVER_HYBRID = 2


def set_run_dir(run_dir):
    global RUN_DIR
    RUN_DIR = run_dir


def _subdomain_specs(test_name):
    """Return (label, interior_x_bounds, full_x_bounds, solver_code) tuples."""
    name = (test_name or "").lower()
    if "three_way" in name:
        return [
            ("left", (-9.1, -3.7), (-math.inf, -3.2), SOLVER_IMPLICIT),
            ("middle", (-2.7, 2.7), (-3.2, 3.2), SOLVER_EXPLICIT),
            ("right", (3.7, 9.1), (3.2, math.inf), SOLVER_HYBRID),
        ]
    if "explicit_hybrid" in name:
        left_s, right_s = SOLVER_EXPLICIT, SOLVER_HYBRID
    elif "implicit_explicit" in name:
        left_s, right_s = SOLVER_IMPLICIT, SOLVER_EXPLICIT
    else:
        left_s, right_s = SOLVER_IMPLICIT, SOLVER_HYBRID
    return [
        ("left", (-5.9, -0.5), (-math.inf, 0.0), left_s),
        ("right", (0.5, 5.9), (0.0, math.inf), right_s),
    ]


def validate_log(pic_diags, test_name=None):
    """Verify finite energies and conservative bounds across the coupled seam."""
    if not pic_diags or len(pic_diags) < 5:
        return False, "pic log is empty or has fewer than 5 records"

    for idx, row in enumerate(pic_diags):
        for key in ("Etot", "Ee", "Eb", "Epart"):
            val = row.get(key, math.nan)
            if not math.isfinite(val) or val < 0.0:
                return False, f"non-finite or negative {key}={val} at log entry {idx}"

    etot0 = pic_diags[0]["Etot"]
    etot1 = pic_diags[-1]["Etot"]
    if etot0 <= 0.0:
        return False, f"initial Etot must be positive, got {etot0}"
    drift = abs(etot1 - etot0) / etot0
    if drift >= 0.15:
        return False, f"total energy relative drift {drift:.4f} exceeds 0.15"

    eb0 = pic_diags[0]["Eb"]
    eb1 = pic_diags[-1]["Eb"]
    if eb0 <= 0.0:
        return False, f"initial Eb must be positive, got {eb0}"
    eb_ratio = eb1 / eb0
    if not (0.80 <= eb_ratio <= 1.20):
        return False, f"magnetic energy ratio Eb[-1]/Eb[0]={eb_ratio:.4f} outside [0.80, 1.20]"

    logger.debug(
        "  [%s] log ok: rows=%d dEtot/Etot=%.4e Eb_ratio=%.4f",
        test_name or "solver_coupling",
        len(pic_diags),
        drift,
        eb_ratio,
    )
    return True, "Passed"


def _read_out(path):
    """Parse an SWMF ASCII .out file into (col_map, rows)."""
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    if len(lines) < 6:
        raise ValueError(f"Malformed .out file (too few lines): {path}")

    # Header line 2 stores: iStep time nDim nParam nPlotVar (negative nDim for AMR).
    # Header line 3 stores the grid point/cell counts (first token is nPts when nDim < 0).
    # Header line 5 lists (|nDim| + nPlotVar) column names followed by nParam
    # scalar parameter names, while each data row contains (|nDim| + nPlotVar) values.
    hdr = lines[1].split()
    if len(hdr) < 5:
        raise ValueError(f"Malformed .out header line 2 in {path}: {lines[1]!r}")
    raw_ndim = int(hdr[2])
    n_dim = abs(raw_ndim)
    n_plotvar = int(hdr[4])
    n_cols = n_dim + n_plotvar
    grid_dims = [int(t) for t in lines[2].split()]
    if len(grid_dims) < (1 if raw_ndim < 0 else n_dim):
        raise ValueError(f"Malformed .out grid line 3 in {path}: {lines[2]!r}")
    n_pts = grid_dims[0] if raw_ndim < 0 else math.prod(grid_dims[:n_dim])

    all_names = lines[4].split()
    if n_cols <= 0 or len(all_names) < n_cols:
        raise ValueError(
            f"Malformed .out column header in {path}: need {n_cols}, got {len(all_names)}"
        )
    cols = all_names[:n_cols]
    col_map = {c.lower(): idx for idx, c in enumerate(cols)}

    tokens = []
    for line in lines[5:]:
        stripped = line.strip()
        if not stripped:
            continue
        tokens.extend(stripped.split())

    expected_tokens = n_pts * n_cols
    if len(tokens) < expected_tokens:
        raise ValueError(
            f"Truncated .out data in {path}: need {expected_tokens}, got {len(tokens)}"
        )
    tokens = tokens[:expected_tokens]

    rows = []
    for r in range(n_pts):
        base = r * n_cols
        rows.append([float(tokens[base + c]) for c in range(n_cols)])
    return col_map, rows


def _group_1d_x_profile(col_map, rows):
    """Average fields onto unique sorted x positions for 1D spectral checks."""
    ix = col_map["x"]
    iby = col_map["by"]
    ibz = col_map["bz"]
    buckets = {}
    for row in rows:
        xk = round(row[ix], 8)
        entry = buckets.get(xk)
        if entry is None:
            buckets[xk] = [row[ix], row[iby], row[ibz], 1]
        else:
            entry[1] += row[iby]
            entry[2] += row[ibz]
            entry[3] += 1
    profile = []
    for xk in sorted(buckets):
        x, by_sum, bz_sum, cnt = buckets[xk]
        profile.append((x, by_sum / cnt, bz_sum / cnt))
    return profile


def _mode1_power_and_coeff(profile):
    """Compute (mode1_fraction, c_plus, c_minus) from a uniform 1D (x, By, Bz) profile."""
    n = len(profile)
    if n < 8:
        raise ValueError(f"Too few x points ({n}) for mode-1 Fourier check")
    xs = [p[0] for p in profile]
    dx = (xs[-1] - xs[0]) / (n - 1)
    lx = dx * n
    k1 = 2.0 * math.pi / lx

    mean_by = sum(p[1] for p in profile) / n
    mean_bz = sum(p[2] for p in profile) / n
    by_ac = [p[1] - mean_by for p in profile]
    bz_ac = [p[2] - mean_bz for p in profile]

    total_ac_power = 0.0
    for by, bz in zip(by_ac, bz_ac):
        total_ac_power += by * by + bz * bz
    if total_ac_power <= 0.0:
        raise ValueError("Zero transverse AC magnetic energy in final frame")

    # Compute real-space Fourier power at harmonic m=1 summed over By and Bz.
    # For real signals f_j, Parseval with unitary norm gives total_ac_power = sum |F_k|^2 / n.
    mode1_power = 0.0
    for sig in (by_ac, bz_ac):
        re = sum(val * math.cos(k1 * x) for x, val in zip(xs, sig))
        im = sum(-val * math.sin(k1 * x) for x, val in zip(xs, sig))
        # Negative frequency -k1 carries equal energy for real signal, so multiply by 2/n.
        mode1_power += 2.0 * (re * re + im * im) / n

    mode1_frac = min(1.0, mode1_power / total_ac_power)

    # Complex transverse wave field W(x) = By(x) + i Bz(x) projected onto exp(+/- i k1 x)
    c_plus = 0j
    c_minus = 0j
    for x, by, bz in zip(xs, by_ac, bz_ac):
        w = complex(by, bz)
        phase = k1 * x
        ep = complex(math.cos(phase), math.sin(phase))
        em = complex(math.cos(phase), -math.sin(phase))
        c_plus += w * ep
        c_minus += w * em
    return mode1_frac, c_plus / n, c_minus / n


def _unwrap_phase(angles):
    """Unwrap a sequence of angles (radians) in chronological order."""
    if not angles:
        return []
    out = [angles[0]]
    for a in angles[1:]:
        diff = a - out[-1]
        while diff > math.pi:
            diff -= 2.0 * math.pi
        while diff < -math.pi:
            diff += 2.0 * math.pi
        out.append(out[-1] + diff)
    return out


def validate_plot(test_name):
    """Check solver partition, divB/divBc, moment continuity, Pe, and wave propagation."""
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    if len(plots) < 2:
        return False, f"expected at least 2 .out frames in {RUN_DIR}/PC/plots, found {len(plots)}"

    specs = _subdomain_specs(test_name)
    name_lower = (test_name or "").lower()
    is_blended_amr = "blended_amr" in name_lower
    is_three_way = "three_way" in name_lower

    max_divb = 0.0
    max_divbc = 0.0
    frames = []
    for path in plots:
        col_map, rows = _read_out(path)
        for req in ("x", "by", "bz", "divb", "divbc", "solver", "rhos0", "rhos1", "ps0", "ps1"):
            if req not in col_map:
                return False, f"missing required plot column '{req}' in {path}"

        idivb = col_map["divb"]
        idivbc = col_map["divbc"]
        for row in rows:
            vb = abs(row[idivb])
            vbc = abs(row[idivbc])
            if not math.isfinite(vb):
                return False, f"non-finite divB in {path}"
            if not math.isfinite(vbc):
                return False, f"non-finite divBc in {path}"
            if vb > max_divb:
                max_divb = vb
            if vbc > max_divbc:
                max_divbc = vbc
        frames.append((path, col_map, rows))

    if max_divb >= 1e-9:
        return False, f"max |divB|={max_divb:.3e} >= 1e-9"
    if not is_blended_amr and max_divbc >= 1e-9:
        return False, f"max |divBc|={max_divbc:.3e} >= 1e-9"

    # Check solver partition, species densities, pressures, and Pe in first and last frames.
    # Exclude periodic wrap boundaries and immediate seam cells via interior_x_bounds.
    for tag, (path, col_map, rows) in (("first", frames[0]), ("last", frames[-1])):
        ix = col_map["x"]
        isol = col_map["solver"]
        ir0, ir1 = col_map["rhos0"], col_map["rhos1"]
        ip0, ip1 = col_map["ps0"], col_map["ps1"]
        ipe = col_map.get("pe")

        pe_means = []
        for label, (ilo, ihi), _, exp_code in specs:
            sub_rows = [r for r in rows if ilo < r[ix] < ihi]
            if not sub_rows:
                return False, f"missing {label} subdomain points in {tag} frame ({path})"
            for r in sub_rows:
                code = int(round(r[isol]))
                if code != exp_code:
                    return False, (
                        f"{tag} frame {path}: expected {label} solver={exp_code} "
                        f"at x={r[ix]:.3f}, got {r[isol]}"
                    )
                for vname, idx in (("rhoS0", ir0), ("rhoS1", ir1), ("pS0", ip0), ("pS1", ip1)):
                    if not math.isfinite(r[idx]) or r[idx] <= 0.0:
                        return False, (
                            f"{tag} frame {path}: non-positive {vname}={r[idx]} at x={r[ix]:.3f}"
                        )
            if ipe is not None:
                pe_means.append(sum(r[ipe] for r in sub_rows) / len(sub_rows))

        if ipe is not None:
            for r in rows:
                if not math.isfinite(r[ipe]) or r[ipe] <= 0.0:
                    return False, f"{tag} frame {path}: non-positive Pe={r[ipe]} at x={r[ix]:.3f}"
            pe_max = max(pe_means)
            rel_pe = (pe_max - min(pe_means)) / pe_max
            if rel_pe > 0.30:
                return False, (
                    f"{tag} frame {path}: subdomain mean Pe spread {rel_pe:.3f} > 0.30 "
                    f"(means={[round(p, 4) for p in pe_means]})"
                )

    # Transverse wave balance and Fourier coherence across the seam(s).
    _, col_last, rows_last = frames[-1]
    ix, iby, ibz = col_last["x"], col_last["by"], col_last["bz"]
    rms_vals = []
    for label, _, (flo, fhi), _ in specs:
        bperp = [r[iby] ** 2 + r[ibz] ** 2 for r in rows_last if flo <= r[ix] < fhi]
        if not bperp:
            return False, f"no {label} subdomain points in final frame"
        rms = math.sqrt(sum(bperp) / len(bperp))
        if rms <= 0.0:
            return False, f"zero {label}-subdomain rms |B_perp| in final frame"
        rms_vals.append(rms)

    if is_blended_amr:
        rms_ratio = rms_vals[0] / rms_vals[-1]
        if not (0.50 <= rms_ratio <= 2.00):
            return False, (
                f"blended_amr rms_Bperp_left/right={rms_ratio:.3f} outside [0.50, 2.00]"
            )
        logger.debug(
            "  [%s] plot ok: max|divB|=%.2e max|divBc|=%.2e rmsL/R=%.3f",
            test_name,
            max_divb,
            max_divbc,
            rms_ratio,
        )
        return True, "Passed"

    if is_three_way:
        rms_spread = max(rms_vals) / min(rms_vals)
        if rms_spread > 1.65:
            return False, f"final three_way rms_Bperp spread {rms_spread:.3f} > 1.65"
        rms_diag = f"rms_spread={rms_spread:.3f}"
    else:
        rms_ratio = rms_vals[0] / rms_vals[-1]
        if not (0.65 <= rms_ratio <= 1.55):
            return False, f"final rms_Bperp_left/right={rms_ratio:.3f} outside [0.65, 1.55]"
        rms_diag = f"rmsL/R={rms_ratio:.3f}"

    prof_last = _group_1d_x_profile(col_last, rows_last)
    mode1_frac, _, _ = _mode1_power_and_coeff(prof_last)
    if mode1_frac < 0.75:
        return False, f"final mode-1 power fraction={mode1_frac:.3f} < 0.75"

    # Measure unwrapped phase advance on the dominant seeded circular polarization.
    _, col_first, rows_first = frames[0]
    prof_first = _group_1d_x_profile(col_first, rows_first)
    _, cp0, cm0 = _mode1_power_and_coeff(prof_first)
    use_plus = abs(cp0) >= abs(cm0)

    angles = []
    for _, col_m, rows_f in frames:
        prof = _group_1d_x_profile(col_m, rows_f)
        _, cp, cm = _mode1_power_and_coeff(prof)
        c = cp if use_plus else cm
        angles.append(math.atan2(c.imag, c.real))
    unwrapped = _unwrap_phase(angles)
    dphase = unwrapped[-1] - unwrapped[0]
    if dphase <= 0.1:
        return False, f"insufficient right-hand wave phase advance dphase={dphase:.4f} rad"

    logger.debug(
        "  [%s] plot ok: max|divB|=%.2e max|divBc|=%.2e %s mode1=%.3f dphase=%.3f rad",
        test_name or "solver_coupling",
        max_divb,
        max_divbc,
        rms_diag,
        mode1_frac,
        dphase,
    )
    return True, "Passed"
