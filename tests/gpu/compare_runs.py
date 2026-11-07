#!/usr/bin/env python3
"""Run FLEKS test decks and compare the outputs of two runs.

The GPU port keeps one code path for the CPU and the GPU. This script checks
the two reproducibility promises of the port:

  * a CPU build gives bitwise-identical results before and after a change
    that is not meant to change them (compare two CPU runs), and
  * a GPU build gives bitwise-identical results run to run (run the GPU
    executable twice and compare), and agrees with the CPU build to round-off
    (compare a GPU run with a CPU run using --rtol).

Usage:
  compare_runs.py run --exe bin/FLEKS.exe --out /tmp/runA [--np 1]
                      [--fleks-args "amrex.the_arena_is_managed=1"]
                      beam/PARAM.in.hybrid whistler/PARAM.in ...
  compare_runs.py compare /tmp/runA /tmp/runB [--rtol 0]

A deck is given relative to the tests directory. 'run' keeps every output
(validate_tests.py deletes them) in OUT/<deck>/. 'compare' compares all the
files in the PC directories of the two runs: identical files pass; text files
that differ are compared number by number with the relative tolerance rtol
(0 = bitwise); binary files must be identical unless rtol > 0, in which case
they are reported as differing.
"""

import argparse
import filecmp
import math
import os
import re
import shutil
import subprocess
import sys
import time

TESTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(TESTS)
NUMBER = re.compile(
    r"[-+]?(?:(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?"
    r"|(?<![A-Za-z_])(?:nan|inf)(?![A-Za-z_]))",
    re.IGNORECASE,
)


def deck_name(deck):
    return deck.replace("/", "__")


def prepare_run_dir(deck, out):
    """Create OUT/<deck>/ with the deck as PARAM.in and the data files next
    to it, and return its path."""
    src = os.path.join(TESTS, deck)
    run_dir = os.path.join(out, deck_name(deck))
    if os.path.isdir(run_dir):
        shutil.rmtree(run_dir)
    os.makedirs(os.path.join(run_dir, "PC", "plots"))
    os.makedirs(os.path.join(run_dir, "PC", "restartOUT"))
    shutil.copy(src, os.path.join(run_dir, "PARAM.in"))
    # Data files next to the deck (e.g. field coefficient tables).
    deck_dir = os.path.dirname(src)
    for f in os.listdir(deck_dir):
        path = os.path.join(deck_dir, f)
        if (os.path.isfile(path) and not f.startswith("PARAM.in")
                and not f.endswith((".py", ".md"))):
            shutil.copy(path, run_dir)
    return run_dir


def run(args):
    exe = os.path.abspath(args.exe)
    os.makedirs(args.out, exist_ok=True)
    status = 0
    for deck in args.decks:
        run_dir = prepare_run_dir(deck, args.out)
        cmd = [exe] + args.fleks_args.split()
        if args.np > 1:
            cmd = ["mpirun", "-n", str(args.np)] + cmd
        t0 = time.time()
        with open(os.path.join(run_dir, "stdout.log"), "w") as log:
            ret = subprocess.call(cmd, cwd=run_dir, stdout=log,
                                  stderr=subprocess.STDOUT)
        print(f"{deck:45s} exit={ret} {time.time() - t0:7.1f}s")
        status |= ret != 0
    return status


def numbers(path):
    with open(path, errors="replace") as f:
        return [float(x.replace("d", "e").replace("D", "e"))
                for x in NUMBER.findall(f.read())]


def is_text(path):
    with open(path, "rb") as f:
        return b"\0" not in f.read(4096)


def compare_files(a, b, rtol):
    """Return (ok, message) for two files."""
    if filecmp.cmp(a, b, shallow=False):
        return True, ""
    if not (is_text(a) and is_text(b)):
        return False, "binary files differ"
    na, nb = numbers(a), numbers(b)
    if len(na) != len(nb):
        return False, f"{len(na)} vs {len(nb)} numbers"
    worst = 0.0
    for x, y in zip(na, nb):
        if not (math.isfinite(x) and math.isfinite(y)):
            if not ((math.isnan(x) and math.isnan(y)) or x == y):
                return False, f"non-finite mismatch ({x} vs {y})"
            continue
        scale = max(abs(x), abs(y))
        if scale > 0:
            worst = max(worst, abs(x - y) / scale)
    if worst <= rtol and rtol > 0:
        return True, f"max rel diff {worst:.2e}"
    return False, f"max rel diff {worst:.2e}"


def list_output_files(pc_dir):
    out = set()
    for root, _, files in os.walk(pc_dir):
        for f in files:
            if f.endswith((".lock",)) or "timing" in f.lower():
                continue
            out.add(os.path.relpath(os.path.join(root, f), pc_dir))
    return out


def compare(args):
    status = 0
    decks_a = {n for n in os.listdir(args.a)
               if os.path.isdir(os.path.join(args.a, n, "PC"))}
    decks_b = {n for n in os.listdir(args.b)
               if os.path.isdir(os.path.join(args.b, n, "PC"))}
    for name in sorted(decks_a | decks_b):
        if name not in decks_a:
            print(f"{name}: missing in {args.a}")
            status = 1
            continue
        if name not in decks_b:
            print(f"{name}: missing in {args.b}")
            status = 1
            continue
        dir_a = os.path.join(args.a, name, "PC")
        dir_b = os.path.join(args.b, name, "PC")
        files_a = list_output_files(dir_a)
        files_b = list_output_files(dir_b)
        all_files = sorted(files_a | files_b)
        bad = []
        for rel in all_files:
            if rel not in files_a:
                bad.append(f"{rel}: missing in {args.a}")
                continue
            if rel not in files_b:
                bad.append(f"{rel}: missing in {args.b}")
                continue
            fa = os.path.join(dir_a, rel)
            fb = os.path.join(dir_b, rel)
            ok, msg = compare_files(fa, fb, args.rtol)
            if not ok:
                bad.append(f"{rel}: {msg}")
        n = len(all_files)
        verdict = "IDENTICAL" if not bad else "DIFFERENT"
        print(f"{name:45s} {verdict} ({n} files)")
        for line in bad[:args.max_report]:
            print("    " + line)
        status |= bool(bad)
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("--exe", default=os.path.join(REPO, "bin", "FLEKS.exe"))
    p.add_argument("--out", required=True)
    p.add_argument("--np", type=int, default=1)
    p.add_argument("--fleks-args", default="")
    p.add_argument("decks", nargs="+")
    p.set_defaults(func=run)
    p = sub.add_parser("compare")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--rtol", type=float, default=0.0)
    p.add_argument("--max-report", type=int, default=5)
    p.set_defaults(func=compare)
    args = parser.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
