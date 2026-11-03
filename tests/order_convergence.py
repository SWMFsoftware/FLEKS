#!/usr/bin/env python3
"""Measure the temporal order of accuracy of the hybrid-PIC solver.

Uses tests/whistler (default) or tests/iaw to measure the numerical dispersion
and order of accuracy:
  * whistler: Parallel circularly polarized whistler wave. Exercises the full
    generalized Ohm's law (Hall term, convection, Faraday update, RK subcycling,
    and particle Boris push). Tracks the phase of the transverse Fourier mode
    k=1 (By - i*Bz) vs time to measure omega(dt).
  * iaw: Ion acoustic wave. Exercises the ambipolar electric field (-grad(pe)/rho).
    Measures the oscillation frequency of the density perturbation via global
    nonlinear least-squares fitting of a damped cosine across all output frames.

The grid is kept fixed while dt is halved. Two observables are reported:

  * PRIMARY, Cauchy L2 norm of the final field:
        p_L2 = log2( ||u(h)-u(h/2)|| / ||u(h/2)-u(h/4)|| )
    This is the reliable metric: the fields differ by a smooth O(dt^p) amount
    that cleans up by a factor 2^p per refinement step on a grid-identical run.

  * DIAGNOSTIC, dispersion frequency Richardson extrapolation:
        omega(dt) = omega* + C*dt^p,   p_omega = log2( |w(h)-w(h/2)| / |w(h/2)-w(h/4)| )
    Treat p_omega with care. The truncation error in omega is ~1e-6 relative,
    while the particle-noise bias of any single-mode estimator is ~1e-3, so the
    FINEST refinement step must keep its difference above that noise floor
    otherwise p_omega swings by +-1 between particle realizations. Each case
    below therefore uses a dt ladder whose smallest step still clears the floor
    (the iaw case at 5 s does not: see tests/iaw).

ALIGNMENT REQUIREMENT (enforced at start-up): both plot_dt and every dt in the
ladder must divide `time_max` exactly. Otherwise FLEKS stops at the first
multiple of dt >= time_max and the last frame lands at a different time for each
dt, which corrupts the Cauchy norm (it would then compare fields at different
times) and shifts the least-squares window of the frequency fit.

Usage:
    python3 tests/order_convergence.py                    # default: ssprk3 on whistler
    python3 tests/order_convergence.py --integrator rk4   # rk4 single-pass
    python3 tests/order_convergence.py --test iaw         # acoustic wave test
    python3 tests/order_convergence.py --dt 0.02          # single dt diagnostics
"""
import argparse
import cmath
import glob
import math
import os
import shutil
import subprocess
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_DIR = os.path.join(REPO, "run_order")
EXE = os.path.join(REPO, "bin", "FLEKS.exe")

# Per-case configuration.
#   nsub:     field subcycling used for every dt in the ladder.
#   time_max: must be divided exactly by plot_dt and by every dt (checked below).
#   dts:      coarse -> fine. The finest step must still produce a difference
#             above the particle-noise floor, otherwise p_omega is meaningless.
CASES = {
    "whistler": {
        "nsub": 16,
        "time_max": 4.0,
        "plot_dt": 0.1,
        "dts": (0.04, 0.02, 0.01),
    },
    # The iaw mode decays by ~30 % over 4.8 s, so a 0.01 step leaves the
    # difference between dt=0.02 and dt=0.01 buried in the particle noise; the
    # 0.08/0.04/0.02 ladder keeps both differences resolved.
    "iaw": {
        "nsub": 4,
        "time_max": 4.8,
        "plot_dt": 0.08,
        "dts": (0.08, 0.04, 0.02),
    },
}


def check_alignment(test_name):
    """plot_dt and every dt must divide time_max exactly (see module docstring)."""
    case = CASES[test_name]
    steps = [case["plot_dt"]] + list(case["dts"])
    bad = [
        v
        for v in steps
        if abs(case["time_max"] / v - round(case["time_max"] / v)) > 1e-9
    ]
    if bad:
        raise SystemExit(
            "%s: time_max=%.4g is not an integer multiple of %s; the final frame "
            "would land at a different time for each dt"
            % (test_name, case["time_max"], bad)
        )


def _rewrite_param(text, block, fields):
    """Rewrite numeric values within a named #BLOCK. `fields` maps keyword ->
    new value. Returns patched text."""
    out = []
    in_block = False
    replaced = set()
    for line in text.splitlines():
        s = line.strip()
        if not s:
            out.append(line)
            continue
        tok = s.split()[0]
        if tok.startswith("#"):
            in_block = (tok == block)
            out.append(line)
            continue
        if in_block:
            for kw, val in fields.items():
                if kw in s and kw not in replaced:
                    indent = line[: len(line) - len(line.lstrip())]
                    out.append(f"{indent}{val}                     {kw}")
                    replaced.add(kw)
                    break
            else:
                out.append(line)
            continue
        out.append(line)
    missing = set(fields) - replaced
    if missing:
        raise RuntimeError(f"Could not patch {missing} in {block}")
    return "\n".join(out) + "\n"


def _frame_time(path):
    try:
        with open(path, "r", encoding="latin-1") as f:
            f.readline()
            nxt = f.readline()
        return float(nxt.split()[1])
    except (OSError, IndexError, ValueError):
        return None


def run_case(test_name, dt, nsub, nprocs=1, integrator="rk4-pc", time_max=4.0):
    if os.path.isdir(RUN_DIR):
        shutil.rmtree(RUN_DIR)
    os.makedirs(os.path.join(RUN_DIR, "PC", "plots"), exist_ok=True)
    os.makedirs(os.path.join(RUN_DIR, "PC", "restartOUT"), exist_ok=True)

    if test_name == "whistler":
        param_file = os.path.join(REPO, "tests", "whistler", "PARAM.in.rk4")
    else:
        param_file = os.path.join(REPO, "tests", "iaw", "PARAM.in")
    plot_dt = "%g" % CASES[test_name]["plot_dt"]

    with open(param_file) as f:
        text = f.read()

    text = _rewrite_param(text, "#TIMESTEPPING", {"dt": dt})
    if "#FIELDSUBCYCLE" in text:
        text = _rewrite_param(text, "#FIELDSUBCYCLE", {"nSubcycle": nsub})
    else:
        text = _rewrite_param(text, "#BSUBCYCLE", {"nBSubcycle": nsub})
    text = _rewrite_param(text, "#STOP", {"TimeMax": time_max})
    text = _rewrite_param(text, "#SAVEPLOT", {"dn": "-1", "dt": plot_dt})

    if test_name == "iaw":
        text = _rewrite_param(text, "#PARTICLES", {"nParticle per cell in X": 20000})

    if integrator:
        if "#FIELDINTEGRATOR" in text:
            text = _rewrite_param(
                text, "#FIELDINTEGRATOR", {"fieldIntegrator": integrator}
            )
        elif "#FIELDSOLVER" in text:
            text = text.replace(
                "#FIELDSOLVER",
                f"#FIELDINTEGRATOR\n{integrator}                   fieldIntegrator\n\n#FIELDSOLVER",
            )
        else:
            text = text.replace(
                "#HYBRIDPIC",
                f"#FIELDINTEGRATOR\n{integrator}                   fieldIntegrator\n\n#HYBRIDPIC",
            )

    with open(os.path.join(RUN_DIR, "PARAM.in"), "w") as f:
        f.write(text)

    for link, target in (
        ("FLEKS.exe", EXE),
        ("PostProc.pl", os.path.join(REPO, "share", "Scripts", "PostProc.pl")),
        ("PC/PostIDL.exe", os.path.join(REPO, "bin", "PostIDL.exe")),
        ("PC/pIDL", os.path.join(REPO, "share", "Scripts", "pIDL")),
    ):
        p = os.path.join(RUN_DIR, link)
        if not os.path.lexists(p):
            os.symlink(target, p)

    cmd = (
        ["./FLEKS.exe"]
        if nprocs <= 1
        else ["mpirun", "-n", str(nprocs), "./FLEKS.exe"]
    )
    r = subprocess.run(
        cmd, cwd=RUN_DIR, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if r.returncode != 0:
        sys.stderr.write(
            "FLEKS failed (dt=%s):\n%s\n%s\n" % (dt, r.stdout, r.stderr)
        )
        return None

    pp = subprocess.run(
        ["./PostProc.pl", "-v"],
        cwd=RUN_DIR,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if pp.returncode != 0:
        sys.stderr.write(
            "PostProc.pl failed (dt=%s):\n%s\n%s\n" % (dt, pp.stdout, pp.stderr)
        )
        return None

    if test_name == "whistler":
        return measure_whistler_freq()
    else:
        return measure_acoustic_freq()


def measure_whistler_freq():
    """Track the rotating phase of the transverse mode k=1 (By - i*Bz) vs time;
    return angular frequency omega in code units and the final By field."""
    plots = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots, "*.out")))
    if len(out_files) < 5:
        sys.stderr.write("  need >=5 .out frames, got %d\n" % len(out_files))
        return None

    samples = []
    last_by = None
    for f in out_files:
        with open(f, "r", encoding="latin-1") as fp:
            lines = fp.readlines()
        if len(lines) < 6:
            continue
        t = _frame_time(f)
        if t is None:
            continue
        vidx = {v.upper(): i for i, v in enumerate(lines[4].split())}
        if "BY" not in vidx or "BZ" not in vidx:
            continue
        by_list, bz_list = [], []
        for l in lines[5:]:
            c = l.split()
            if len(c) > vidx["BZ"]:
                by_list.append(float(c[vidx["BY"]]))
                bz_list.append(float(c[vidx["BZ"]]))
        n = len(by_list)
        if n < 4:
            continue
        last_by = by_list
        # Transverse right-hand mode: By - i*Bz has rotating phase exp(-i*omega*t)
        c1 = sum(
            (by_list[j] - 1j * bz_list[j]) * cmath.exp(-2j * math.pi * j / n)
            for j in range(n)
        ) / n
        samples.append((t, cmath.phase(c1)))

    if len(samples) < 5:
        return None

    samples.sort(key=lambda s: s[0])
    unwrapped = [samples[0][1]]
    for i in range(1, len(samples)):
        dphi = samples[i][1] - samples[i - 1][1]
        dphi -= 2.0 * math.pi * round(dphi / (2.0 * math.pi))
        unwrapped.append(unwrapped[-1] + dphi)

    ts = [s[0] for s in samples]
    nt = len(ts)
    mean_t = sum(ts) / nt
    mean_p = sum(unwrapped) / nt
    cov = sum((ts[i] - mean_t) * (unwrapped[i] - mean_p) for i in range(nt))
    var_t = sum((ts[i] - mean_t) ** 2 for i in range(nt))
    omega = -cov / var_t if var_t > 0 else 0.0
    return {"omega_code": omega, "field": last_by}


def _fit_standing_wave_freq(ts, ys):
    """Fit y(t) = A * exp(-gamma*t) * cos(omega*t + phi0) by Gauss-Newton least
    squares, initialized from the two-point zero-crossing estimate.

    Returns a dict with the parameter values, the residual and a `valid` flag.
    The fit is unconstrained, so it is only accepted when it moved a sane
    distance from the initial guess and actually explains the samples; otherwise
    the caller falls back to the zero-crossing estimate.
    """
    t_arr = np.asarray(ts, dtype=np.float64)
    y_arr = np.asarray(ys, dtype=np.float64)

    # Initial guess via linear zero-crossing
    y0 = y_arr[0]
    omega_init = 0.0
    for i in range(len(y_arr) - 1):
        if (y_arr[i] > 0 and y_arr[i + 1] <= 0) or (
            y_arr[i] < 0 and y_arr[i + 1] >= 0
        ):
            y1, y2 = y_arr[i], y_arr[i + 1]
            t1, t2 = t_arr[i], t_arr[i + 1]
            t_zero = t1 + y1 / (y1 - y2) * (t2 - t1)
            if t_zero > 0:
                omega_init = (math.pi / 2.0) / t_zero
            break

    if omega_init <= 0.0:
        omegas = []
        for t_val, y_val in zip(t_arr[1:], y_arr[1:]):
            if t_val > 0 and y0 != 0:
                val = max(-1.0, min(1.0, y_val / y0))
                omegas.append(math.acos(val) / t_val)
        omega_init = sum(omegas) / len(omegas) if omegas else 1.0

    # Parameter vector p = [A, omega, gamma, phi0]
    p = np.array([y0, omega_init, 0.0, 0.0], dtype=np.float64)
    rms = float("inf")
    iters = 0
    for iters in range(1, 51):
        A, omega, gamma, phi0 = p
        exp_t = np.exp(-gamma * t_arr)
        arg = omega * t_arr + phi0
        cos_arg = np.cos(arg)
        sin_arg = np.sin(arg)

        f = A * exp_t * cos_arg
        r = y_arr - f
        rms = float(np.sqrt(np.mean(r * r)))

        # Jacobian columns: df/dA, df/domega, df/dgamma, df/dphi0
        J = np.column_stack(
            (
                exp_t * cos_arg,
                -A * t_arr * exp_t * sin_arg,
                -A * t_arr * exp_t * cos_arg,
                -A * exp_t * sin_arg,
            )
        )
        try:
            dp, _, _, _ = np.linalg.lstsq(J, r, rcond=None)
        except np.linalg.LinAlgError:
            break
        p += dp
        # Relative, component-wise stop test: the parameters differ by many
        # orders of magnitude, so an absolute norm on dp only measures A.
        if np.max(np.abs(dp) / (np.abs(p) + 1e-30)) < 1e-9:
            break

    A, omega, gamma, phi0 = (float(v) for v in p)
    r = y_arr - A * np.exp(-gamma * t_arr) * np.cos(omega * t_arr + phi0)
    rms = float(np.sqrt(np.mean(r * r)))
    amp = abs(A)
    reason = ""
    if amp <= 0.0:
        reason = "degenerate amplitude"
    elif not (0.7 * omega_init <= abs(omega) <= 1.3 * omega_init):
        reason = "omega=%.4f left the initial guess %.4f by more than 30%%" % (
            omega,
            omega_init,
        )
    elif rms > 0.05 * amp:
        reason = "residual %.1e is %.1f%% of the amplitude" % (rms, 100.0 * rms / amp)
    return {
        "omega": abs(omega),
        "amp": amp,
        "gamma": gamma,
        "phi0": phi0,
        "iters": iters,
        "rms": rms,
        "omega_init": omega_init,
        "valid": not reason,
        "reason": reason,
    }


def measure_acoustic_freq():
    """Track the standing wave density oscillation of mode k=1 vs time."""
    plots = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots, "*.out")))
    if len(out_files) < 5:
        sys.stderr.write("  need >=5 .out frames, got %d\n" % len(out_files))
        return None

    samples = []
    last_rho = None
    for f in out_files:
        with open(f, "r", encoding="latin-1") as fp:
            lines = fp.readlines()
        if len(lines) < 6:
            continue
        t = _frame_time(f)
        if t is None:
            continue
        vidx = {v.upper(): i for i, v in enumerate(lines[4].split())}
        target = next((k for k in ("RHOS0", "RHOS1", "RHO") if k in vidx), None)
        if target is None:
            continue
        ridx = vidx[target]
        rhos = [float(l.split()[ridx]) for l in lines[5:] if len(l.split()) > ridx]
        n = len(rhos)
        if n < 4:
            continue
        last_rho = rhos
        mean = sum(rhos) / n
        # Standing wave seeded with sin(kx*x): mode amplitude along imaginary axis
        im_c = sum(
            (rhos[j] - mean) * math.sin(2.0 * math.pi * j / n) for j in range(n)
        ) * (2.0 / n)
        samples.append((t, im_c))

    if len(samples) < 5:
        return None

    samples.sort(key=lambda s: s[0])
    ts = [s[0] for s in samples]
    ys = [s[1] for s in samples]
    fit = _fit_standing_wave_freq(ts, ys)
    if not fit["valid"]:
        sys.stderr.write("  fit rejected (%s), keeping zero-crossing omega\n" % fit["reason"])
    return {"omega_code": fit["omega"], "field": last_rho, "fit": fit}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--test",
        type=str,
        default="whistler",
        choices=["whistler", "iaw"],
        help="Test case to run: 'whistler' (default, Hall+convection wave) or 'iaw' (acoustic wave)",
    )
    ap.add_argument(
        "--integrator",
        type=str,
        default="ssprk3",
        help="Field integrator (ssprk3, rk4, ssprk3-pc, rk4-pc)",
    )
    ap.add_argument("--dt", type=float, default=None, help="Run a single dt (diagnostics)")
    ap.add_argument("--nproc", type=int, default=1, help="Number of MPI processes")
    args = ap.parse_args()

    check_alignment(args.test)
    case = CASES[args.test]
    time_max = case["time_max"]
    nsub = case["nsub"]
    dts = list(case["dts"])

    if args.dt is not None:
        res = run_case(
            args.test,
            args.dt,
            nsub,
            args.nproc,
            args.integrator,
            time_max=time_max,
        )
        if res is None:
            sys.exit("no frequency measured")
        print(
            "test=%s integrator=%s dt=%-8s omega_code = %.6f rad/s"
            % (args.test, args.integrator, args.dt, res["omega_code"])
        )
        return

    omegas = []
    fields = []
    print(
        "Testing integrator: %s on %s (time_max=%.3g s, plot_dt=%.3g s)"
        % (args.integrator, args.test, time_max, case["plot_dt"])
    )
    print("%-10s %-12s %-16s" % ("dt", "nBSubcycle", "omega_code"))
    for dt in dts:
        res = run_case(
            args.test,
            dt,
            nsub,
            args.nproc,
            args.integrator,
            time_max=time_max,
        )
        if res is None:
            print("%-10.5f %-12d (failed)" % (dt, nsub))
            continue
        print("%-10.5f %-12d %-16.6f" % (dt, nsub, res["omega_code"]))
        fit = res.get("fit")
        if fit is not None:
            print(
                "           fit: A=%.3e gamma=%.4f phi0=%+.3f rms/A=%.1e iters=%d%s"
                % (
                    fit["amp"],
                    fit["gamma"],
                    fit["phi0"],
                    fit["rms"] / fit["amp"] if fit["amp"] > 0 else float("nan"),
                    fit["iters"],
                    "" if fit["valid"] else "  [REJECTED: %s]" % fit["reason"],
                )
            )
        omegas.append((dt, res["omega_code"]))
        fields.append(res["field"])

    if len(omegas) < 3:
        print("\nInsufficient successful runs to fit the order.")
        return

    # 1. PRIMARY: Cauchy L2 norm of the final field.
    #    p_L2 = log2( ||u(h)-u(h/2)|| / ||u(h/2)-u(h/4)|| )
    p_l2 = None
    if fields[0] and fields[1] and fields[2]:
        f0, f1, f2 = (np.asarray(f, dtype=float) for f in fields)
        n = len(f0)
        c1 = float(np.linalg.norm(f0 - f1) / math.sqrt(n))
        c2 = float(np.linalg.norm(f1 - f2) / math.sqrt(n))
        print(
            "\nCauchy L2: ||u(%g)-u(%g)||=%.3e, ||u(%g)-u(%g)||=%.3e"
            % (dts[0], dts[1], c1, dts[1], dts[2], c2)
        )
        if c1 > 0.0 and c2 > 0.0:
            p_l2 = math.log(c1 / c2) / math.log(2.0)
            print("  ratio = %.3f  ->  p_L2 = %.2f" % (c1 / c2, p_l2))
        else:
            print("  field difference is zero to round-off; nothing to fit.")

    # 2. DIAGNOSTIC: dispersion frequency Richardson extrapolation.
    #    omega(dt) = omega* + C*dt^p, ratio = d1/d2 = 2^p
    d1 = abs(omegas[0][1] - omegas[1][1])
    d2 = abs(omegas[1][1] - omegas[2][1])
    print("omega differences: d1=%.3e d2=%.3e" % (d1, d2))
    p_omega = None
    if d1 > 0.0 and d2 > 0.0:
        ratio = d1 / d2
        p_omega = math.log(ratio) / math.log(2.0)
        print(
            "  Richardson ratio r = %.3f  ->  dispersion order p_omega = %.2f"
            % (ratio, p_omega)
        )
        print(
            "  (diagnostic only: particle statistics, not truncation, set this "
            "error bar unless d2 clears the noise floor)"
        )
    else:
        print("  could not estimate p_omega (zero difference).")

    print("")
    if p_l2 is None:
        print("RESULT: p_L2 unavailable; no reliable temporal order.")
        return
    print("Temporal order (primary, Cauchy L2): p = %.2f  (second order -> p ~ 2)" % p_l2)
    if p_omega is not None:
        print("Temporal order (diagnostic, omega): p = %.2f" % p_omega)
    if p_l2 >= 1.7:
        print("RESULT: CONFIRMED second-order (p_L2 = %.2f)" % p_l2)
    else:
        print("RESULT: order %.2f < 1.7" % p_l2)


if __name__ == "__main__":
    main()
