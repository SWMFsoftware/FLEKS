#!/usr/bin/env python3
"""Measure the temporal order of accuracy of the hybrid-PIC solver.

Uses tests/whistler (default) or tests/iaw to measure the numerical dispersion
and order of accuracy:
  * whistler: Parallel circularly polarized whistler wave. Exercises the full
    generalized Ohm's law (Hall term, convection, Faraday update, RK subcycling,
    and particle Boris push). Tracks the phase of the transverse Fourier mode
    k=1 (By - i*Bz) vs time to measure omega(dt).
  * iaw: Ion acoustic wave. Exercises the ambipolar electric field (-grad(pe)/rho).
    Measures the oscillation frequency of the density perturbation.

The grid is kept fixed while dt is halved. The measured frequency omega(dt) is
evaluated via Richardson extrapolation:
    omega(dt) = omega* + C * dt^p
    p = log( |w(h) - w(h/2)| / |w(h/2) - w(h/4)| ) / log(2)
Effective temporal order p is reported (second order -> p ~ 2).

Usage:
    python3 tests/order_convergence.py                    # default: rk4-pc on whistler
    python3 tests/order_convergence.py --integrator rk4   # legacy rk4
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

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_DIR = os.path.join(REPO, "run_order")
EXE = os.path.join(REPO, "bin", "FLEKS.exe")


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
        plot_dt = "0.1"
    else:
        param_file = os.path.join(REPO, "tests", "iaw", "PARAM.in")
        plot_dt = "0.2"

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
    # Estimate frequency from zero crossing of cos(omega*t): omega*t0 = pi/2
    y0 = samples[0][1]
    t_zero = None
    for i in range(len(samples) - 1):
        if (samples[i][1] > 0 and samples[i + 1][1] <= 0) or (
            samples[i][1] < 0 and samples[i + 1][1] >= 0
        ):
            y1, y2 = samples[i][1], samples[i + 1][1]
            t1, t2 = samples[i][0], samples[i + 1][0]
            t_zero = t1 + y1 / (y1 - y2) * (t2 - t1)
            break

    if t_zero and t_zero > 0:
        omega = (math.pi / 2.0) / t_zero
    else:
        # Fallback to linear arccos approximation
        omegas = []
        for t, y in samples[1:]:
            val = max(-1.0, min(1.0, y / y0))
            if t > 0:
                omegas.append(math.acos(val) / t)
        omega = sum(omegas) / len(omegas) if omegas else 0.0

    return {"omega_code": omega, "field": last_rho}


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
        default="rk4-pc",
        help="Field integrator (rk4-pc, rk4, ssprk3-pc, ssprk3)",
    )
    ap.add_argument("--dt", type=float, default=None, help="Run a single dt (diagnostics)")
    ap.add_argument("--nproc", type=int, default=1, help="Number of MPI processes")
    args = ap.parse_args()

    time_max = 4.0 if args.test == "whistler" else 5.0
    nsub = 16 if args.test == "whistler" else 4

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

    # Refinement sequence for testing temporal order
    dts = [0.04, 0.02, 0.01] if args.test == "whistler" else [0.04, 0.02, 0.01]
    omegas = []
    fields = []
    print("Testing integrator: %s on %s" % (args.integrator, args.test))
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
        omegas.append((dt, res["omega_code"]))
        fields.append(res["field"])

    if len(omegas) >= 3:
        # Richardson extrapolation:
        # omega(dt) = omega* + C * dt^p
        # ratio = (omega(h) - omega(h/2)) / (omega(h/2) - omega(h/4)) = 2^p
        d1 = abs(omegas[0][1] - omegas[1][1])
        d2 = abs(omegas[1][1] - omegas[2][1])
        print("\nomega differences: ['%.3e', '%.3e']" % (d1, d2))

        if d2 > 0:
            ratio = d1 / d2
            p = math.log(ratio) / math.log(2.0)
            print("  Richardson ratio r = %.3f  ->  effective order p = %.2f" % (ratio, p))
            print("\nMean effective order p = %.2f  (second order -> p ~ 2)" % p)
            if p >= 1.7:
                print("RESULT: CONFIRMED second-order (p = %.2f)" % p)
            else:
                print("RESULT: order %.2f < 1.7" % p)
        else:
            print("\nCould not estimate the order (differences too small).")

        # Cauchy L2 error check if fields are available
        if len(fields) >= 3 and fields[0] and fields[1] and fields[2]:
            import numpy as np

            f0, f1, f2 = (
                np.array(fields[0]),
                np.array(fields[1]),
                np.array(fields[2]),
            )
            c1 = np.linalg.norm(f0 - f1) / np.sqrt(len(f0))
            c2 = np.linalg.norm(f1 - f2) / np.sqrt(len(f0))
            c_p = math.log(c1 / c2) / math.log(2.0) if c2 > 0 else 0.0
            print(
                "Cauchy L2 field differences: ||u(0.04)-u(0.02)||=%.3e, ||u(0.02)-u(0.01)||=%.3e -> p=%.2f"
                % (c1, c2, c_p)
            )
    else:
        print("\nInsufficient successful runs to fit the order.")


if __name__ == "__main__":
    main()
