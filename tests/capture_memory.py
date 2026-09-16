#!/usr/bin/env python3
"""Memory regression, step 1 of 2: capture a memory profile.

To check whether your changes regress memory, run this and then
``tests/compare_memory.py`` (step 2 of 2), which does the comparison::

    python3 tests/capture_memory.py --out mine.json    # where you are now
    python3 tests/capture_memory.py --out base.json --ref master   # reference
    python3 tests/compare_memory.py base.json mine.json

In CI ``.github/workflows/memory_test.yml`` runs both for you on the same
runner and posts the result on the PR.

Every standalone FLEKS run already ends with the AMReX TinyProfiler report and
the FLEKS load-balance report; this script just runs a curated selection of the
standalone tests with the profiler switched fully on and writes one JSON file
holding the per-region allocation counts / peak bytes plus the RSS.  The
selection is deliberately small: it has to run twice (reference + candidate)
inside one CI job.

Other options::

    --verify    run twice and report how reproducible memory is here, which is
                what --rss-tol in step 2 should be tuned from
    --list      show the captured selection
    --test X    restrict to one entry of the selection
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from tests import profiler  # noqa: E402
from tests import validate_tests  # noqa: E402

# The profiling selection: (test directory, PARAM.in suffix or None, nprocs).
#
# Chosen to cover the dominant cost centres of both field solvers while
# staying short: full-PIC implicit E solve + particle mover (beam), hybrid
# Ohm assembly + Faraday advance (performance/PARAM.in.hybrid), 2D moment
# deposition and current calculation (reconnection), boundary injection
# (shock), and one 2-rank entry so MPI-related allocations are covered too.
PROFILE_TESTS = [
    ("beam", None, 1),
    ("performance", "hybrid", 1),
    ("reconnection", None, 1),
    ("shock", None, 1),
    ("beam", None, 2),
]

DEFAULT_RUN_DIR = "run_test_prof"

# Force the profiler to report every region instead of folding the small ones
# into "Other" (AMReX default print_threshold is 1 %), and to write the
# report to a file so it does not have to be scraped out of the physics log.
PROFILER_ARGS = ["tiny_profiler.output_file=prof.txt",
                 "tiny_profiler.print_threshold=0"]

# Which of the periodic RSS reports to keep. A run emits one table per report
# step (63 of them for the beam deck); storing the whole series per test would
# bloat the profile document for no benefit, and the two ends are what carry
# the signal: the first shows the settled baseline, the last the end state.
LOAD_BALANCE_KEEP = ("first", "last")


def git(*args):
    """Return git output from the repository root, or '' on failure."""
    try:
        out = subprocess.run(("git",) + args, cwd=REPO_ROOT,
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        return out.stdout.decode("utf-8", "replace").strip()
    except OSError:
        return ""


def param_path(test, variant):
    name = "PARAM.in" + ("." + variant if variant else "")
    return os.path.join(REPO_ROOT, "tests", test, name)


def test_key(test, variant, nprocs):
    return f"{test}{'.' + variant if variant else ''}.n{nprocs}"


def prepare_run_dir(run_dir, exe):
    """Minimal run-directory setup: no PostIDL / PostProc needed for profiling."""
    os.makedirs(run_dir, exist_ok=True)
    link = os.path.join(run_dir, "FLEKS.exe")
    if os.path.lexists(link):
        os.remove(link)
    os.symlink(os.path.abspath(exe), link)


def clean_output(run_dir):
    """Drop plot/restart output between runs, keeping the directory skeleton."""
    for sub in (os.path.join(run_dir, "PC", "plots"),
                os.path.join(run_dir, "PC", "restartOUT")):
        if not os.path.isdir(sub):
            continue
        for entry in os.listdir(sub):
            path = os.path.join(sub, entry)
            try:
                shutil.rmtree(path) if os.path.isdir(path) and not \
                    os.path.islink(path) else os.remove(path)
            except OSError:
                pass


def run_one(test, variant, nprocs, run_dir, exe, keep_prof=False,
            keep_log=False):
    """Run one entry and return its profile record."""
    prepare_run_dir(run_dir, exe)
    clean_output(run_dir)

    source = param_path(test, variant)
    with open(source) as handle:
        param_text = handle.read()
    with open(os.path.join(run_dir, "PARAM.in"), "w") as handle:
        handle.write(param_text)

    if nprocs > 1:
        cmd = ["mpirun", "-n", str(nprocs), "./FLEKS.exe"] + PROFILER_ARGS
    else:
        cmd = ["./FLEKS.exe"] + PROFILER_ARGS

    env = dict(os.environ, OMP_NUM_THREADS="1")
    start = time.monotonic()
    result = subprocess.run(cmd, cwd=run_dir, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    wall_s = time.monotonic() - start
    log = (result.stdout or b"").decode("utf-8", "replace")
    result.stdout = log

    prof_file = os.path.join(run_dir, "prof.txt")
    record = {
        "test": test,
        "variant": variant,
        "nprocs": nprocs,
        "param": os.path.relpath(source, REPO_ROOT),
        "wall_s": round(wall_s, 4),
        "exit_code": result.returncode,
        "profile": {},
        "load_balance": {},
    }

    if keep_log:
        with open(os.path.join(run_dir, "run.log"), "w") as handle:
            handle.write(log)

    if result.returncode != 0:
        record["error"] = log[-4000:]
        return record

    if os.path.isfile(prof_file):
        with open(prof_file) as handle:
            record["profile"] = profiler.parse_tinyprofiler(handle.read())
        if not keep_prof:
            os.remove(prof_file)
    else:
        record["error"] = "no TinyProfiler output; is AMREX_TINY_PROFILING on?"

    # RSS. Pic::report_load_balance() prints one table per report step; keeping
    # every one would bloat the document for no benefit, so only the first and
    # the last are retained (see LOAD_BALANCE_KEEP).
    tables = profiler.parse_load_balance(log)
    if tables:
        record["load_balance"] = {"n_tables": len(tables)}
        for name in LOAD_BALANCE_KEEP:
            index = 0 if name == "first" else len(tables) - 1
            record["load_balance"][name] = tables[index]

    return record


def capture(selection, run_dir, exe, verbose=False, keep_prof=False,
            keep_log=False):
    """Run every entry in *selection* and return the full profile document."""
    if not os.path.isfile(exe):
        print(f"error: executable not found: {exe}\n"
              f"  build it first, or pass --exe with the path to one.",
              file=sys.stderr)
        sys.exit(2)

    os.makedirs(run_dir, exist_ok=True)
    tests = {}
    for test, variant, nprocs in selection:
        ok, reason, _skip = validate_tests.preflight_check(test)
        if not ok:
            print(f"  SKIP {test_key(test, variant, nprocs)}: {reason}")
            continue
        if not os.path.isfile(param_path(test, variant)):
            print(f"  SKIP {test_key(test, variant, nprocs)}: "
                  f"missing {param_path(test, variant)}")
            continue

        key = test_key(test, variant, nprocs)
        print(f"  RUN  {key} ...", flush=True)
        record = run_one(test, variant, nprocs, run_dir, exe,
                         keep_prof=keep_prof, keep_log=keep_log)
        tests[key] = record
        if "error" in record:
            print(f"       FAILED (exit {record['exit_code']})")
        elif verbose:
            print(profiler.format_summary(record["profile"], top=5))
            rss = record.get("load_balance", {}).get("last", {}).get("Memory(MB)")
            if rss:
                print(f"       RSS last report (MB): min {rss['min']:.1f} "
                      f"avg {rss['avg']:.1f} max {rss['max']:.1f}")

    return {
        "meta": {
            "ref": git("rev-parse", "--abbrev-ref", "HEAD"),
            "commit": git("rev-parse", "HEAD"),
            "commit_subject": git("log", "-1", "--format=%s"),
            "dirty": bool(git("status", "--porcelain")),
            "amrex_dim": validate_tests.configured_amrex_dim(),
            "nlevmax": validate_tests.configured_nlevmax(),
            "exo_source": validate_tests.configured_user_source_is_exo(),
            "omp_num_threads": os.environ.get("OMP_NUM_THREADS", "1"),
            "date": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "profiler_args": PROFILER_ARGS,
            "exe": os.path.relpath(exe, REPO_ROOT),
        },
        "tests": tests,
    }


# ---------------------------------------------------------------------------
# Determinism check
# ---------------------------------------------------------------------------
def _memory_index(profile):
    """Flat {arena.region: {nalloc, maxmem_max, curmem_max}} for comparisons."""
    out = {}
    for arena, regions in (profile.get("memory") or {}).items():
        for region, values in regions.items():
            out[f"{arena}::{region}"] = {
                k: values.get(k) for k in ("nalloc", "maxmem_max", "curmem_max")
            }
    return out


def _rss_index(record):
    """Flat {rss.<first|last>.<min|avg|max>: MB} for the determinism check."""
    out = {}
    for name, table in (record.get("load_balance") or {}).items():
        if name == "n_tables" or not isinstance(table, dict):
            continue
        rss = table.get("Memory(MB)")
        if rss:
            out[f"rss.{name}"] = (rss["min"], rss["avg"], rss["max"])
    return out


def verify(selection, run_dir, exe, repeats=2):
    """Run the selection twice and report non-deterministic memory values.

    Allocation counts and peak bytes are exact integers for a fixed problem,
    rank count and RNG seed.  If they are not reproducible the memory checks
    in compare_memory.py can only be used as warnings, not as a gate.
    """
    runs = []
    for i in range(repeats):
        print(f"run {i + 1}/{repeats}")
        runs.append(capture(selection, run_dir, exe))

    baseline = runs[0]
    problems, rss_spread = [], []
    for key, record in baseline["tests"].items():
        if "error" in record:
            problems.append(f"{key}: run 1 failed")
            continue
        reference = _memory_index(record["profile"])
        rss_reference = _rss_index(record)
        for i, other in enumerate(runs[1:], start=2):
            candidate = other.get("tests", {}).get(key, {})
            if "error" in candidate:
                problems.append(f"{key}: run {i} failed")
                continue
            current = _memory_index(candidate["profile"])
            for name, values in reference.items():
                got = current.get(name)
                if got != values:
                    problems.append(
                        f"{key}: {name} differs between run 1 and run {i}: "
                        f"{values} vs {got}")
            for name, values in rss_reference.items():
                got = _rss_index(candidate).get(name)
                if got is None:
                    continue
                for column, base, cand in zip(("min", "avg", "max"), values, got):
                    if base != cand:
                        rss_spread.append((key, f"{name}.{column}", base, cand,
                                           abs(cand - base),
                                           abs(cand - base) / base if base else 0))

    if problems:
        print("\nARENA MEMORY IS NOT DETERMINISTIC:")
        for problem in problems:
            print(f"  {problem}")
        print("\nThe allocation checks in compare_memory.py must be treated "
              "as warnings on this setup.")
        return 1

    print("\nArena memory (nalloc / maxmem_max / curmem_max) is bit-identical "
          f"across {repeats} runs -- it can be gated strictly.")

    if not rss_spread:
        print("RSS is bit-identical too.")
        return 0

    worst_abs = max(r[4] for r in rss_spread)
    worst_rel = max(r[5] for r in rss_spread)
    worst = max(rss_spread, key=lambda r: r[4])
    print(f"RSS is not bit-identical: {len(rss_spread)} of "
          f"{len(rss_spread) + 0} compared values moved between runs.")
    print(f"  largest absolute change: {worst_abs:.1f} MB "
          f"({worst[0]} {worst[1]}: {worst[2]:.1f} -> {worst[3]:.1f})")
    print(f"  largest relative change: {100 * worst_rel:.1f} %")
    print(f"  -> gate with a tolerance well above that, e.g. "
          f"--rss-tol {max(1.0, round(3 * worst_abs, 1)):g}")
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Memory regression step 1 of 2: capture a memory profile.")
    parser.add_argument("--out", help="write the profile document as JSON")
    parser.add_argument("--ref", help="label recorded in meta['ref']")
    parser.add_argument("--run-dir", default=DEFAULT_RUN_DIR,
                        help=f"run directory (default: {DEFAULT_RUN_DIR})")
    parser.add_argument("--exe", default=os.path.join("bin", "FLEKS.exe"),
                        help="FLEKS executable to measure (default bin/FLEKS.exe); "
                             "use this to compare against a binary built from a "
                             "different commit without rebuilding it")
    parser.add_argument("--test", action="append",
                        help="restrict to test[.variant][:nprocs]; repeatable")
    parser.add_argument("--list", action="store_true", help="show the selection")
    parser.add_argument("--verify", action="store_true",
                        help="run twice and check memory determinism")
    parser.add_argument("--repeats", type=int, default=2, help="--verify runs")
    parser.add_argument("--keep-prof", action="store_true",
                        help="keep the raw prof.txt of the last test")
    parser.add_argument("--keep-log", action="store_true",
                        help="keep the raw run.log of the last test")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    os.chdir(REPO_ROOT)

    if args.list:
        for test, variant, nprocs in PROFILE_TESTS:
            print(f"  {test_key(test, variant, nprocs):<28} "
                  f"{os.path.relpath(param_path(test, variant), REPO_ROOT)}")
        return 0

    selection = PROFILE_TESTS
    if args.test:
        wanted = set()
        for item in args.test:
            name, _, nproc = item.partition(":")
            wanted.add((name, int(nproc) if nproc else None))
        selection = [e for e in PROFILE_TESTS
                     if (e[0] if not e[1] else f"{e[0]}.{e[1]}", None) in wanted
                     or any(w[0] in (e[0], f"{e[0]}.{e[1]}") and
                            (w[1] is None or w[1] == e[2]) for w in wanted)]
        if not selection:
            print(f"no selection matches {args.test}; see --list")
            return 1

    exe = args.exe if os.path.isabs(args.exe) else os.path.join(REPO_ROOT, args.exe)

    if args.verify:
        return verify(selection, args.run_dir, exe, repeats=args.repeats)

    print("Capturing profiles:")
    document = capture(selection, args.run_dir, exe, verbose=args.verbose,
                       keep_prof=args.keep_prof, keep_log=args.keep_log)
    if args.ref:
        document["meta"]["ref"] = args.ref

    if args.out:
        with open(args.out, "w") as handle:
            json.dump(document, handle, indent=2, sort_keys=True)
        print(f"wrote {args.out}")
    else:
        print(json.dumps(document, indent=2, sort_keys=True))

    failed = [k for k, v in document["tests"].items() if "error" in v]
    if failed:
        print(f"\n{len(failed)} test(s) produced no profile: {', '.join(failed)}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
