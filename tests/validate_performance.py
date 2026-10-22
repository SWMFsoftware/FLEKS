#!/usr/bin/env python3
"""Performance regression runner for standalone FLEKS.

Runs the full-PIC beam, hybrid-PIC whistler, particle tracker and 2D
reconnection benchmarks on 1 and 2 MPI processes (3 runs each), parses AMReX
TinyProfiler output to extract particle-mover and field-solver timings, and
reports particle-step rates (μs/part-step) and the 2-process speedup.

Two modes:

* With --base-exe (the CI): the candidate and a reference executable (the
  merge base) run on the same machine, alternating run by run, and a metric
  fails only if the candidate is more than --tolerance slower than the
  reference. Absolute timings are not comparable across CI runners, whose
  hardware differs by up to a factor of two, so this is the regression check.
  Each executable runs the benchmark decks of its own commit (--base-decks).
* Without it (a local run): the candidate is compared with the absolute
  targets in BASELINES, which describe the slowest GitHub runner class.
"""
import argparse
import os
import shutil
import subprocess
import re
import sys
import time

MOVERS = [
    "Pts::charged_particle_mover",
    "TestParticles::move_charged_particles",
]

FULLPIC_SOLVERS = [
    "Pic::update_E_impl",
    "Pic::E_iterate",
    "Pic::E_matvec",
]

HYBRID_SOLVERS = [
    "Pic::update_B_hybrid",
    "Pic::assemble_ohm_E",
]

PT_SOLVERS = []

_REQUIRED_MOVER = {
    "fullpic": "Pts::charged_particle_mover",
    "hybrid":  "Pts::charged_particle_mover",
    "pt":      "TestParticles::move_charged_particles",
    "reconnection2d": "Pts::charged_particle_mover",
}

# Absolute targets (μs/part-step) and 2-core speedup floor, one set per solver,
# used only without a reference executable.
BASELINES = {
    "fullpic": {
        "total_pps": 6.2,    # total wall-clock
        "mover_pps": 0.105,  # isolated particle mover
        "solver_pps": 4.7,   # isolated implicit field solver
        "speedup": 1.5,      # 2-core scaling floor
    },
    "hybrid": {
        "total_pps": 0.31,   # total wall-clock
        "mover_pps": 0.052,  # isolated particle mover
        "solver_pps": 0.105, # explicit hybrid field advance is cheaper
        "speedup": 1.45,     # 2-core scaling floor
    },
    "pt": {
        "total_pps": 3.2,    # total wall-clock rate
        "mover_pps": 0.14,   # isolated test particle mover rate
        "solver_pps": 0.001, # EM solver disabled
        "speedup": 1.9,      # 2-core mover speedup floor
    },
    "reconnection2d": {
        "total_pps": 0.97,   # total wall-clock rate
        "mover_pps": 0.056,  # isolated particle mover rate
        "solver_pps": 0.02,  # isolated implicit field solver rate
        "speedup": 1.85,     # 2-core scaling floor
    },
}

def safe_symlink(src, dst):
    if os.path.lexists(dst):
        if os.path.islink(dst) or os.path.isfile(dst):
            os.remove(dst)
        elif os.path.isdir(dst):
            shutil.rmtree(dst)
    os.symlink(src, dst)

def prepare_run_dir(run_dir, exe):
    os.makedirs(run_dir, exist_ok=True)
    safe_symlink(os.path.abspath(exe), os.path.join(run_dir, "FLEKS.exe"))

    # Component plots and restart directories
    pc_dir = os.path.join(run_dir, "PC")
    os.makedirs(pc_dir, exist_ok=True)
    os.makedirs(os.path.join(pc_dir, "plots"), exist_ok=True)
    os.makedirs(os.path.join(pc_dir, "restartOUT"), exist_ok=True)

def _param_float(block, idx):
    """Read the idx-th numeric token of a PARAM block (list of stripped lines)."""
    try:
        return float(block[idx].split()[0])
    except (IndexError, ValueError):
        return None


def _read_param_block(param_file, command):
    """Return the (stripped) lines belonging to a #COMMAND block of a PARAM file."""
    block = []
    capture = False
    with open(param_file) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                capture = (s.split()[0].upper() == command.upper())
                continue
            if capture:
                block.append(s)
    return block


def count_particles_from_param(param_file):
    """Total seeded macroparticles = prod(nCell_d * nPartPerCell_d) * nSpecies
    from the #NCELL, #PARTICLES (or #TPPARTICLES), and #PLASMA blocks.  FLEKS seeds
    nPartPerCell particles per cell for each species, so this is exact for benchmarks.
    Returns None if the blocks are missing/malformed."""
    plasma = _read_param_block(param_file, "#PLASMA")
    ns = int(_param_float(plasma, 0)) if plasma else 1
    ncell = _read_param_block(param_file, "#NCELL")
    ppc = _read_param_block(param_file, "#PARTICLES")
    if not ppc:
        ppc = _read_param_block(param_file, "#TPPARTICLES")
    if len(ncell) < 3 or len(ppc) < 3:
        return None
    total = ns
    for d in range(3):
        nc = _param_float(ncell, d)
        np = _param_float(ppc, d)
        if nc is None or np is None:
            return None
        total *= int(nc) * int(np)
    return total


def cycles_from_param(param_file):
    """Number of timesteps for a fixed-dt run, from #TIMESTEPPING dt and
    #STOP (MaxIter if > 0, else round(TimeMax / dt)).  Returns None if the
    values cannot be determined."""
    dt_block = _read_param_block(param_file, "#TIMESTEPPING")
    stop_block = _read_param_block(param_file, "#STOP")
    use_fixed = dt_block[0].split()[0] if len(dt_block) >= 1 else None
    dt = _param_float(dt_block, 1) if len(dt_block) >= 2 else None
    max_iter = _param_float(stop_block, 0) if len(stop_block) >= 1 else None
    time_max = _param_float(stop_block, 1) if len(stop_block) >= 2 else None
    if max_iter is not None and int(max_iter) > 0:
        return int(max_iter)
    if use_fixed and use_fixed.upper() == "T" and dt and dt > 0 and time_max:
        return int(round(time_max / dt))
    return None

def parse_tiny_profiler(stdout):
    """Parse AMReX TinyProfiler total time, exclusive and inclusive timing
    tables from FLEKS stdout into nested dicts."""
    profiler_data = {}

    total_time_match = re.search(r"TinyProfiler total time across processes\s+\[min\.\.\.avg\.\.\.max\]:\s+([\d.]+)", stdout)
    if total_time_match:
        profiler_data["total_time"] = float(total_time_match.group(1))
    else:
        total_time_match2 = re.search(r"TinyProfiler total time across processes.*?:\s+([\d.]+)", stdout)
        if total_time_match2:
            profiler_data["total_time"] = float(total_time_match2.group(1))

    excl_section_match = re.search(
        r"Name\s+NCalls\s+Excl\.\s+Min\s+Excl\.\s+Avg\s+Excl\.\s+Max\s+Max\s+%\n-+([\s\S]+?)-+",
        stdout
    )
    if excl_section_match:
        excl_text = excl_section_match.group(1)
        excl_entries = {}
        for line in excl_text.splitlines():
            line = line.strip()
            if not line:
                continue
            match = re.match(r"^(.+?)\s+(\d+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)%", line)
            if match:
                name = match.group(1).strip()
                excl_entries[name] = {
                    "calls": int(match.group(2)),
                    "min": float(match.group(3)),
                    "avg": float(match.group(4)),
                    "max": float(match.group(5)),
                    "percent": float(match.group(6))
                }
        profiler_data["exclusive"] = excl_entries

    incl_section_match = re.search(
        r"Name\s+NCalls\s+Incl\.\s+Min\s+Incl\.\s+Avg\s+Incl\.\s+Max\s+Max\s+%\n-+([\s\S]+?)-+",
        stdout
    )
    if incl_section_match:
        incl_text = incl_section_match.group(1)
        incl_entries = {}
        for line in incl_text.splitlines():
            line = line.strip()
            if not line:
                continue
            match = re.match(r"^(.+?)\s+(\d+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)%", line)
            if match:
                name = match.group(1).strip()
                incl_entries[name] = {
                    "calls": int(match.group(2)),
                    "min": float(match.group(3)),
                    "avg": float(match.group(4)),
                    "max": float(match.group(5)),
                    "percent": float(match.group(6))
                }
        profiler_data["inclusive"] = incl_entries

    return profiler_data


SUITES = ("fullpic", "hybrid", "pt", "reconnection2d")

LABELS = {
    "fullpic": "FULL PIC",
    "hybrid": "HYBRID PIC",
    "pt": "PARTICLE TRACKER",
    "reconnection2d": "2D RECONNECTION (FULL PIC)",
}

DECKS = {
    "fullpic": "PARAM.in",
    "hybrid": "PARAM.in.hybrid",
    "pt": "PARAM.in.pt",
    "reconnection2d": "PARAM.in.reconnection2d",
}

# (key, report label, unit) of the compared metrics; a rate passes when it is
# low, the speedup when it is high.
METRICS = (
    ("total", "Total Wall-Clock Rate", "μs/pt"),
    ("mover", "Particle Mover Rate", "μs/pt"),
    ("solver", "Field Solver Rate", "μs/pt"),
    ("speedup", "Parallel Speedup (2 Cores)", "x"),
)

# Below this a rate is a disabled component (the tracker has no field solver),
# not something to compare.
RATE_FLOOR = 1e-6


def run_once(exe, n_proc, run_dir, param_file, solver_kind):
    """Run one benchmark and return its wall-clock, mover and solver times
    (seconds). Raises RuntimeError if FLEKS fails."""
    if solver_kind == "hybrid":
        solver_names = HYBRID_SOLVERS
    elif solver_kind == "pt":
        solver_names = PT_SOLVERS
    else:
        solver_names = FULLPIC_SOLVERS

    prepare_run_dir(run_dir, exe)
    shutil.copy(param_file, os.path.join(run_dir, "PARAM.in"))

    start_time = time.perf_counter()
    result = subprocess.run(["mpirun", "-n", str(n_proc), "./FLEKS.exe"],
                            cwd=run_dir, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    duration = time.perf_counter() - start_time
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"{exe} failed with exit code {result.returncode}")

    prof = parse_tiny_profiler(result.stdout)

    # Mover time (exclusive); warn if the required mover region is absent.
    required_mover = _REQUIRED_MOVER[solver_kind]
    mover_time = 0.0
    for name in MOVERS:
        if name in prof.get("exclusive", {}):
            mover_time = prof["exclusive"][name]["avg"]
            if name == required_mover:
                break
    if required_mover not in prof.get("exclusive", {}):
        print(f"  ⚠ WARNING: Required mover region '{required_mover}' NOT "
              f"FOUND in the TinyProfiler exclusive table of {exe}.")

    # Solver time (inclusive); warn if no solver region was found.
    solver_time = 0.0
    for name in solver_names:
        if name in prof.get("inclusive", {}):
            solver_time = prof["inclusive"][name]["avg"]
            break
    else:
        if solver_names:
            print(f"  ⚠ WARNING: No solver region in {solver_names} found in "
                  f"the TinyProfiler inclusive table of {exe}.")

    return {"duration": duration, "mover_time": mover_time,
            "solver_time": solver_time}


def median(values):
    s = sorted(values)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def summarize(runs):
    durations = [r["duration"] for r in runs]
    return {
        "runs": runs,
        "duration_median": median(durations),
        "duration_min": min(durations),
        "mover_median": median([r["mover_time"] for r in runs]),
        "solver_median": median([r["solver_time"] for r in runs]),
    }


def evaluate_suite(solver_kind, targets, count):
    """Benchmark one suite for every target {name: (exe, param_file)} on 1 and
    2 processes, alternating the targets run by run so drift of the machine
    affects them alike. Returns {name: metrics} or None on a failed run."""
    label = LABELS[solver_kind]
    steps = {}
    for name, (exe, param_file) in targets.items():
        particles = count_particles_from_param(param_file)
        cycles = cycles_from_param(param_file)
        if particles is None or cycles is None:
            print(f"Error: could not derive particle/cycle counts from "
                  f"{param_file}.")
            return None
        steps[name] = (particles, cycles)

    runs = {(name, n_proc): [] for name in targets for n_proc in (1, 2)}
    for n_proc in (1, 2):
        print(f"Running {label} on {n_proc} MPI process(es), {count} runs "
              f"per executable ({', '.join(targets)})...")
        for i in range(count):
            for name, (exe, param_file) in targets.items():
                print(f"  Run {i + 1}/{count} [{name}]...")
                try:
                    runs[(name, n_proc)].append(run_once(
                        exe, n_proc, f"run_perf_{name}_{n_proc}", param_file,
                        solver_kind))
                except RuntimeError as exc:
                    print(f"FAIL: {label}: {exc}")
                    return None

    metrics = {}
    for name in targets:
        particles, cycles = steps[name]
        total_steps = particles * cycles
        serial = summarize(runs[(name, 1)])
        parallel = summarize(runs[(name, 2)])
        # For test particles the mover is the computation; the total
        # wall-clock is dominated by redistributing the many attributes and
        # by trajectory output, so the speedup is that of the mover.
        if solver_kind == "pt" and parallel["mover_median"] > 0:
            speedup = serial["mover_median"] / parallel["mover_median"]
        else:
            speedup = serial["duration_median"] / parallel["duration_median"]
        metrics[name] = {
            "particles": particles,
            "cycles": cycles,
            "total_steps": total_steps,
            "serial": serial,
            "parallel": parallel,
            "total": serial["duration_median"] * 1e6 / total_steps,
            "mover": serial["mover_median"] * 1e6 / total_steps,
            "solver": serial["solver_median"] * 1e6 / total_steps,
            "speedup": speedup,
        }
    return metrics


def judge_absolute(solver_kind, m):
    """Pass flags against the absolute targets (no reference executable)."""
    base = BASELINES[solver_kind]
    return {
        "total": m["total"] <= base["total_pps"],
        "mover": m["mover"] <= base["mover_pps"],
        "solver": m["solver"] <= base["solver_pps"],
        "speedup": m["speedup"] >= base["speedup"],
    }


def judge_relative(cand, ref, tolerance):
    """Pass flags and candidate/reference ratios. A rate fails if it is more
    than `tolerance` above the reference; the speedup fails if it is more than
    `tolerance` below it. A rate the reference does not have (a disabled
    component) is not compared."""
    passed, ratio = {}, {}
    for key, _, _ in METRICS:
        if key != "speedup" and ref[key] < RATE_FLOOR:
            passed[key], ratio[key] = True, None
            continue
        ratio[key] = cand[key] / ref[key]
        passed[key] = (ratio[key] >= 1.0 - tolerance if key == "speedup"
                       else ratio[key] <= 1.0 + tolerance)
    return passed, ratio


def _status_markdown(flag):
    return "🟢 **PASSED**" if flag else "🔴 **FAILED**"


def _status_console(flag):
    return "PASSED" if flag else "FAILED"


def _speedup_label(solver_kind):
    return ("Parallel Speedup (Mover, 2 Cores)" if solver_kind == "pt"
            else "Parallel Speedup (2 Cores)")


def _fmt(key, value):
    return f"{value:.2f}x" if key == "speedup" else f"{value:.3f} μs/pt"


def _runs_rows(m, tag):
    rows = []
    for title, s in ((f"1 MPI Process (Serial){tag}", m["serial"]),
                     (f"2 MPI Processes (Parallel){tag}", m["parallel"])):
        cells = " | ".join(f"{r['duration']:.3f}s" for r in s["runs"])
        rows.append(f"| {title} | {cells} | {s['duration_median']:.3f}s | "
                    f"{s['duration_min']:.3f}s |")
    return rows


def write_report(results, relative, tolerance, path="performance_summary.md"):
    with open(path, "w") as f:
        f.write("### ⚡ Standalone FLEKS Performance Report\n\n")
        if relative:
            f.write(f"The PR build and the merge-base build ran on the same "
                    f"runner, alternating run by run (median of the runs per "
                    f"executable). A metric fails if the PR is more than "
                    f"{tolerance:.0%} slower than the base.\n\n")
        else:
            f.write("No reference build: compared with the absolute targets "
                    "(median of the runs).\n\n")
        for solver_kind in SUITES:
            if solver_kind not in results:
                continue
            res = results[solver_kind]
            cand = res["metrics"]["pr"]
            f.write(f"\n#### {LABELS[solver_kind]}\n\n")
            if res["relative"]:
                ref = res["metrics"]["base"]
                f.write("| Performance Metric | PR (Median) | Base (Median) | "
                        "PR / Base | Status |\n")
                f.write("| :--- | :--- | :--- | :--- | :--- |\n")
                for key, name, _ in METRICS:
                    if key == "speedup":
                        name = _speedup_label(solver_kind)
                    r = res["ratio"][key]
                    f.write(f"| {name} | {_fmt(key, cand[key])} | "
                            f"{_fmt(key, ref[key])} | "
                            f"{'n/a' if r is None else f'{r:.3f}'} | "
                            f"{_status_markdown(res['passed'][key])} |\n")
            else:
                base = BASELINES[solver_kind]
                targets = {"total": base["total_pps"], "mover": base["mover_pps"],
                           "solver": base["solver_pps"], "speedup": base["speedup"]}
                if relative:
                    f.write("*No base deck for this benchmark: compared with "
                            "the absolute targets.*\n\n")
                f.write("| Performance Metric | Measured (Median) | Target | "
                        "Status |\n")
                f.write("| :--- | :--- | :--- | :--- |\n")
                for key, name, _ in METRICS:
                    if key == "speedup":
                        name = _speedup_label(solver_kind)
                        target = f">= {targets[key]:.2f}x"
                    else:
                        target = f"<= {targets[key]:.3f} μs/pt"
                    f.write(f"| {name} | {_fmt(key, cand[key])} | {target} | "
                            f"{_status_markdown(res['passed'][key])} |\n")
            part_label = "Test particles" if solver_kind == "pt" else "Macroparticles"
            f.write(f"*{part_label}: {cand['particles']}, cycles: "
                    f"{cand['cycles']} (total steps: {cand['total_steps']}).*\n")
            f.write("\n**Detailed Runs (Wall-Clock Runtime)**\n")
            n_runs = len(cand["serial"]["runs"])
            f.write("| Process Count | " +
                    " | ".join(f"Run {k + 1}" for k in range(n_runs)) +
                    " | Median | Minimum |\n")
            f.write("| :--- | " + " | ".join(":---" for _ in range(n_runs)) +
                    " | :--- | :--- |\n")
            if res["relative"]:
                rows = (_runs_rows(cand, ", PR") +
                        _runs_rows(res["metrics"]["base"], ", base"))
            else:
                rows = _runs_rows(cand, "")
            f.write("\n".join(rows) + "\n")
        f.write("\n*Note: Benchmark ran on virtualized runner.*\n")


def print_profile(solver_kind, res, tolerance):
    cand = res["metrics"]["pr"]
    print("\n" + "=" * 85)
    print(" " * 26 + f"PERFORMANCE PROFILE ({LABELS[solver_kind]})")
    print("=" * 85)
    if res["relative"]:
        ref = res["metrics"]["base"]
        print(f" {'Metric / Component':<30} | {'PR':<14} | {'Base':<14} | "
              f"{'PR/Base':<7} | Status")
        print("-" * 85)
        for key, name, _ in METRICS:
            if key == "speedup":
                name = _speedup_label(solver_kind)
            r = res["ratio"][key]
            print(f" {name:<30} | {_fmt(key, cand[key]):<14} | "
                  f"{_fmt(key, ref[key]):<14} | "
                  f"{'n/a' if r is None else f'{r:.3f}':<7} | "
                  f"{_status_console(res['passed'][key])}")
        print(f" (fails beyond {tolerance:.0%} slower than the base)")
    else:
        base = BASELINES[solver_kind]
        targets = {"total": base["total_pps"], "mover": base["mover_pps"],
                   "solver": base["solver_pps"], "speedup": base["speedup"]}
        print(f" {'Metric / Component':<30} | {'Measured':<16} | "
              f"{'Target':<18} | Status")
        print("-" * 85)
        for key, name, _ in METRICS:
            if key == "speedup":
                name = _speedup_label(solver_kind)
                target = f">= {targets[key]:.2f}x"
            else:
                target = f"<= {targets[key]:.3f} μs/pt"
            print(f" {name:<30} | {_fmt(key, cand[key]):<16} | {target:<18} | "
                  f"{_status_console(res['passed'][key])}")
    print("=" * 85)


def main():
    parser = argparse.ArgumentParser(
        description="Performance regression runner for standalone FLEKS.")
    parser.add_argument("--suite", choices=list(SUITES) + ["all"],
                        default="all",
                        help="Benchmark suite to run (default: all)")
    parser.add_argument("--count", type=int, default=3,
                        help="Runs per process count and executable "
                             "(default: 3)")
    parser.add_argument("--pt", action="store_true",
                        help="Run the particle tracker benchmark suite")
    parser.add_argument("--exe", default="../bin/FLEKS.exe",
                        help="Candidate executable (default: ../bin/FLEKS.exe)")
    parser.add_argument("--base-exe",
                        help="Reference executable (e.g. the merge base). "
                             "Enables the relative comparison.")
    parser.add_argument("--base-decks",
                        help="Directory with the reference commit's "
                             "benchmark decks (default: the candidate's "
                             "performance/ decks)")
    parser.add_argument("--tolerance", type=float, default=0.15,
                        help="Allowed slowdown relative to the reference "
                             "(default: 0.15)")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    exe = os.path.abspath(args.exe)
    base_exe = os.path.abspath(args.base_exe) if args.base_exe else None
    base_decks = os.path.abspath(args.base_decks) if args.base_decks else None
    os.chdir(script_dir)

    for path in [exe] + ([base_exe] if base_exe else []):
        if not os.path.exists(path):
            print(f"Error: executable {path} not found. Please compile it "
                  f"first.")
            sys.exit(1)

    print("=" * 80)
    print("                 FLEKS ROBUST PERFORMANCE REGRESSION TEST")
    print("=" * 80)
    if base_exe:
        print(f"Candidate: {exe}\nReference: {base_exe}\n"
              f"Tolerance: {args.tolerance:.0%}")

    if args.pt or args.suite == "pt":
        suites = ["pt"]
    elif args.suite == "all":
        suites = list(SUITES)
    else:
        suites = [args.suite]

    results = {}
    all_passed = True
    for solver_kind in suites:
        deck = os.path.join("performance", DECKS[solver_kind])
        if not os.path.exists(deck):
            print(f"Error: benchmark deck {deck} not found.")
            sys.exit(1)
        targets = {"pr": (exe, deck)}
        if base_exe:
            base_deck = os.path.join(base_decks or "performance",
                                     DECKS[solver_kind])
            if os.path.exists(base_deck):
                # Reference first, so neither side always runs on a warm cache.
                targets = {"base": (base_exe, base_deck), "pr": (exe, deck)}
            else:
                print(f"  No reference deck {base_deck}: {solver_kind} is "
                      f"compared with the absolute targets.")

        print("\n" + "#" * 85)
        print(f"#  SOLVER: {LABELS[solver_kind]}")
        print("#" * 85)
        metrics = evaluate_suite(solver_kind, targets, args.count)
        if metrics is None:
            all_passed = False
            continue
        res = {"metrics": metrics, "relative": "base" in metrics}
        if res["relative"]:
            res["passed"], res["ratio"] = judge_relative(
                metrics["pr"], metrics["base"], args.tolerance)
        else:
            res["passed"] = judge_absolute(solver_kind, metrics["pr"])
        results[solver_kind] = res
        print_profile(solver_kind, res, args.tolerance)
        all_passed = all_passed and all(res["passed"].values())

    try:
        write_report(results, base_exe is not None, args.tolerance)
    except Exception as e:
        print(f"Warning: Could not write performance_summary.md: {e}")

    for d in os.listdir("."):
        if d.startswith("run_perf_") and os.path.isdir(d):
            shutil.rmtree(d)

    if all_passed:
        msg = "ALL PERFORMANCE TESTS PASSED SUCCESSFULLY!"
        print(f"\033[92;1m\n{msg}\033[0m\n" if sys.stdout.isatty()
              else f"\n{msg}\n")
        sys.exit(0)
    msg = "PERFORMANCE REGRESSION DETECTED."
    print(f"\033[91;1m\n{msg}\033[0m\n" if sys.stdout.isatty() else f"\n{msg}\n")
    sys.exit(1)


if __name__ == "__main__":
    main()
