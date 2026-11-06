#!/usr/bin/env python3
"""Profile FLEKS runs with Nsight Systems and count launches and syncs.

Small FLEKS decks (and FLX-coupled runs with 8^3 boxes) are bound by kernel
launches and stream synchronizations, not by kernel throughput. This script
runs a deck under `nsys`, then attributes every CUDA runtime call to the
innermost FLEKS timer region (timing_func / BL_PROFILE, which the AMReX
TinyProfiler also emits as NVTX ranges), so a change can be judged by how
many launches and syncs it removes from each region.

Usage:
  profile_run.py run --exe FLEKS.exe --out /tmp/profA [--fleks-args ...]
                     beam/PARAM.in.implicit whistler/PARAM.in ...
  profile_run.py cmd --dir RUN_DIR -- mpirun -n 1 ./FLX.exe   (any command,
                     e.g. an FLX coupled test in its prepared run directory)
  profile_run.py show /tmp/profA/<deck> [--top 15]
  profile_run.py analyze /tmp/profA/<deck>   (redo the summary from prof.sqlite)
  profile_run.py compare /tmp/profA /tmp/profB [--top 15]

'run' prepares OUT/<deck>/ as compare_runs.py does. Each profiled directory
gets prof.nsys-rep, prof.sqlite, stdout.log and profile.json (the summary that
'show' and 'compare' read). 'compare' takes two profiled directories, or two
OUT directories, in which case it compares the decks they share.

On a shared workstation, nsys runs inside a memory-capped systemd scope
(--mem, default 10G; --no-scope to disable).

Counts per step are taken inside the --step-region NVTX range (default
Domain::update, one per PIC step), so initialization is left out.
"""

import argparse
import json
import os
import re
import sqlite3
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_runs import REPO, prepare_run_dir  # noqa: E402

PUSH_POP = 59  # NVTX_EVENTS.eventType of a push/pop range
CATEGORIES = ("launch", "sync", "memcpy", "memset", "alloc")


def category(api):
    """The category of a CUDA runtime/driver API name, or None."""
    name = re.sub(r"_v\d+$", "", api)
    if "LaunchKernel" in name or "LaunchCooperativeKernel" in name:
        return "launch"
    if name.endswith("Synchronize"):
        return "sync"
    if "Memcpy" in name:
        return "memcpy"
    if "Memset" in name:
        return "memset"
    if re.search(r"Malloc|Free|HostAlloc", name):
        return "alloc"
    return None


# ---------------------------------------------------------------- profiling
def nsys_cmd(args, cmd, report):
    prof = ["nsys", "profile", "-t", "cuda,nvtx", "-s", "none",
            "--cpuctxsw=none", "-f", "true", "-o", report,
            "--export=sqlite"] + cmd
    if not args.no_scope:
        prof = ["systemd-run", "--user", "--scope", "-q", "-p",
                f"MemoryMax={args.mem}", "-p", "MemorySwapMax=0"] + prof
    return prof


def profile_dir(args, run_dir, cmd):
    report = os.path.join(run_dir, "prof")
    t0 = time.time()
    with open(os.path.join(run_dir, "stdout.log"), "w") as log:
        ret = subprocess.call(nsys_cmd(args, cmd, report), cwd=run_dir,
                              stdout=log, stderr=subprocess.STDOUT)
    wall = time.time() - t0
    if ret != 0 or not os.path.exists(report + ".sqlite"):
        print(f"{run_dir}: exit={ret}, no profile")
        return 1
    summary = analyze(report + ".sqlite", os.path.join(run_dir, "stdout.log"),
                      args.step_region)
    summary["wall_nsys"] = wall
    summary["cmd"] = cmd
    with open(os.path.join(run_dir, "profile.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print_summary(os.path.basename(run_dir.rstrip("/")), summary, args.top)
    return 0


def run(args):
    exe = os.path.abspath(args.exe)
    os.makedirs(args.out, exist_ok=True)
    status = 0
    for deck in args.decks:
        run_dir = prepare_run_dir(deck, args.out)
        cmd = [exe] + args.fleks_args.split()
        if args.np > 1:
            cmd = ["mpirun", "-n", str(args.np)] + cmd
        status |= profile_dir(args, run_dir, cmd)
    return status


def run_cmd(args):
    cmd = args.command[1:] if args.command[:1] == ["--"] else args.command
    return profile_dir(args, os.path.abspath(args.dir), cmd)


# ---------------------------------------------------------------- analysis
def kernel_label(short, demangled):
    """A readable kernel name: AMReX kernels (launch_global) are named after
    the functions that enclose their lambdas, e.g. sum_moments/deposit_to_nodes.
    """
    if "launch_global" not in short or not demangled:
        return short
    names = []
    pos = demangled.find(")::[lambda")
    while pos >= 0:
        # Back to the matching '(' of the argument list, then the name.
        depth, i = 0, pos
        while i >= 0:
            depth += {")": 1, "(": -1}.get(demangled[i], 0)
            if depth == 0:
                break
            i -= 1
        j, angle = i - 1, 0
        while j >= 0:
            c = demangled[j]
            if c == ">":
                angle += 1
            elif c == "<":
                if angle == 0:
                    break
                angle -= 1
            elif angle == 0 and not (c.isalnum() or c in "_:"):
                break
            j -= 1
        name = re.sub(r"<[^<>]*>", "", demangled[j + 1:i])
        while "<" in name:
            name = re.sub(r"<[^<>]*>", "", name)
        name = name.split("::")[-1]
        if name and not demangled[j + 1:i].startswith("amrex::") and \
                name not in names:
            names.append(name)
        pos = demangled.find(")::[lambda", pos + 1)
    return "/".join(names[:3]) or short


def table_exists(db, name):
    return db.execute("SELECT 1 FROM sqlite_master WHERE name=?",
                      (name,)).fetchone() is not None


def tiny_profiler(stdout):
    """Exclusive times of the first TinyProfiler table: name -> [n, sec]."""
    out = {}
    try:
        lines = open(stdout, errors="replace").read().splitlines()
    except OSError:
        return out
    started = False
    for line in lines:
        if line.startswith("Name") and "Excl." in line:
            started = True
            continue
        if not started or line.startswith("---"):
            continue
        parts = line.split()
        if len(parts) < 6:
            if out:
                break
            continue
        try:
            out[" ".join(parts[:-5])] = [int(parts[-5]), float(parts[-3])]
        except ValueError:
            break
    return out


def analyze(sqlite_path, stdout, step_region):
    db = sqlite3.connect(sqlite_path)
    strings = dict(db.execute("SELECT id, value FROM StringIds"))

    # NVTX push/pop ranges per thread, sorted for the nesting sweep.
    ranges = {}
    if table_exists(db, "NVTX_EVENTS"):
        for start, end, text, text_id, tid in db.execute(
                "SELECT start, end, text, textId, globalTid FROM NVTX_EVENTS "
                "WHERE eventType=? AND end IS NOT NULL", (PUSH_POP,)):
            name = text if text is not None else strings.get(text_id, "?")
            ranges.setdefault(tid, []).append((start, -end, name))
    for r in ranges.values():
        r.sort()

    calls = {}
    if table_exists(db, "CUPTI_ACTIVITY_KIND_RUNTIME"):
        for start, end, name_id, tid, corr in db.execute(
                "SELECT start, end, nameId, globalTid, correlationId FROM "
                "CUPTI_ACTIVITY_KIND_RUNTIME"):
            cat = category(strings.get(name_id, ""))
            if cat:
                calls.setdefault(tid, []).append((start, end, cat, corr))

    zero = lambda: {c: 0 for c in CATEGORIES} | {  # noqa: E731
        "sync_s": 0.0, "gpu_s": 0.0}
    total, in_step = zero(), zero()
    regions = {}
    region_of_launch = {}  # correlationId -> region
    steps = 0
    for tid, evs in calls.items():
        evs.sort()
        rs = ranges.get(tid, [])
        stack = []  # open ranges: (end, name)
        ir = 0
        for start, end, cat, corr in evs:
            while ir < len(rs) and rs[ir][0] <= start:
                while stack and stack[-1][0] < rs[ir][0]:
                    stack.pop()
                stack.append((-rs[ir][1], rs[ir][2]))
                ir += 1
            while stack and stack[-1][0] < start:
                stack.pop()
            region = stack[-1][1] if stack else "(none)"
            if cat == "launch":
                region_of_launch[corr] = region
            dt = (end - start) * 1e-9 if cat == "sync" else 0.0
            for d, inc in ((total, True),
                           (in_step, any(n == step_region for _, n in stack)),
                           (regions.setdefault(region, zero()), True)):
                if inc:
                    d[cat] += 1
                    d["sync_s"] += dt
    for rs in ranges.values():
        steps += sum(1 for r in rs if r[2] == step_region)

    kernels = {}
    if table_exists(db, "CUPTI_ACTIVITY_KIND_KERNEL"):
        labels = {}
        for short, dem, corr, t in db.execute(
                "SELECT shortName, demangledName, correlationId, end - start "
                "FROM CUPTI_ACTIVITY_KIND_KERNEL"):
            key = (short, dem)
            if key not in labels:
                labels[key] = kernel_label(strings.get(short, "?"),
                                           strings.get(dem, ""))
            k = kernels.setdefault(labels[key], [0, 0.0])
            k[0] += 1
            k[1] += t * 1e-9
            region = region_of_launch.get(corr)
            if region is not None:
                regions[region]["gpu_s"] += t * 1e-9
    gpu_s = sum(t for _, t in kernels.values())

    tiny = tiny_profiler(stdout)
    m = re.search(r"TinyProfiler total time across processes \[min\.\.\.avg"
                  r"\.\.\.max\]: \S+ \.\.\. (\S+)",
                  open(stdout, errors="replace").read())
    return {"tiny_total": float(m.group(1)) if m else None,
            "steps": steps, "gpu_kernel_s": gpu_s, "total": total,
            "in_step": in_step, "regions": regions, "kernels": kernels,
            "tiny": tiny}


# ---------------------------------------------------------------- reports
def per_step(s, key):
    return s["in_step"][key] / s["steps"] if s["steps"] else float("nan")


def print_summary(name, s, top):
    print(f"== {name}")
    print(f"   TinyProfiler total {s['tiny_total']} s, GPU kernel time "
          f"{s['gpu_kernel_s']:.3f} s, {s['steps']} steps")
    t = s["total"]
    print("   total:    " + ", ".join(f"{c} {t[c]}" for c in CATEGORIES) +
          f", sync wait {t['sync_s']:.3f} s")
    if s["steps"]:
        print("   per step: " + ", ".join(
            f"{c} {per_step(s, c):.0f}" for c in CATEGORIES))
    if top <= 0:
        return
    for key, label in (("sync", "syncs"), ("launch", "launches"),
                       ("gpu_s", "GPU time")):
        print(f"   regions by {label} (launch/sync/memcpy/alloc, sync wait,"
              " GPU time):")
        rows = sorted(s["regions"].items(), key=lambda kv: -kv[1].get(key, 0))
        for region, r in rows[:top]:
            if r.get(key, 0) == 0:
                break
            print(f"     {region:45s} {r['launch']:8d} {r['sync']:7d} "
                  f"{r['memcpy']:7d} {r['alloc']:6d} {r['sync_s']:8.3f} s "
                  f"{r.get('gpu_s', 0):8.3f} s")
    print("   kernels by time (count, s):")
    for k, (n, t) in sorted(s["kernels"].items(), key=lambda kv: -kv[1][1])[
            :top]:
        print(f"     {k[:70]:70s} {n:8d} {t:8.3f}")


def load(d):
    with open(os.path.join(d, "profile.json")) as f:
        return json.load(f)


def compare_one(name, a, b, top):
    print(f"== {name}   (A -> B)")
    print(f"   TinyProfiler total {a['tiny_total']} -> {b['tiny_total']} s, "
          f"GPU kernel {a['gpu_kernel_s']:.3f} -> {b['gpu_kernel_s']:.3f} s")
    for c in CATEGORIES:
        print(f"   {c:7s} total {a['total'][c]:9d} -> {b['total'][c]:9d}"
              f"   per step {per_step(a, c):9.1f} -> {per_step(b, c):9.1f}")
    print(f"   sync wait {a['total']['sync_s']:.3f} -> "
          f"{b['total']['sync_s']:.3f} s")
    if top <= 0:
        return
    names = set(a["regions"]) | set(b["regions"])
    z = {c: 0 for c in CATEGORIES} | {"sync_s": 0.0}
    rows = []
    for n in names:
        ra, rb = a["regions"].get(n, z), b["regions"].get(n, z)
        rows.append((abs(ra["sync"] - rb["sync"]) +
                     abs(ra["launch"] - rb["launch"]) / 4, n, ra, rb))
    rows.sort(reverse=True)
    print("   largest region changes (launch A->B, sync A->B, "
          "excl. time A->B):")
    for _, n, ra, rb in rows[:top]:
        ta = a["tiny"].get(n, [0, 0.0])[1]
        tb = b["tiny"].get(n, [0, 0.0])[1]
        print(f"     {n:42s} {ra['launch']:7d} -> {rb['launch']:7d}  "
              f"{ra['sync']:6d} -> {rb['sync']:6d}  {ta:7.3f} -> {tb:7.3f}")


def compare(args):
    if os.path.exists(os.path.join(args.a, "profile.json")):
        compare_one(os.path.basename(args.a.rstrip("/")), load(args.a),
                    load(args.b), args.top)
        return 0
    for name in sorted(os.listdir(args.a)):
        da, db = os.path.join(args.a, name), os.path.join(args.b, name)
        if (os.path.exists(os.path.join(da, "profile.json")) and
                os.path.exists(os.path.join(db, "profile.json"))):
            compare_one(name, load(da), load(db), args.top)
    return 0


def reanalyze(args):
    for d in args.dirs:
        summary = analyze(os.path.join(d, "prof.sqlite"),
                          os.path.join(d, "stdout.log"), args.step_region)
        old = load(d) if os.path.exists(os.path.join(d, "profile.json")) \
            else {}
        summary["wall_nsys"] = old.get("wall_nsys")
        summary["cmd"] = old.get("cmd")
        with open(os.path.join(d, "profile.json"), "w") as f:
            json.dump(summary, f, indent=1)
        print_summary(os.path.basename(d.rstrip("/")), summary, args.top)
    return 0


def show(args):
    for d in args.dirs:
        print_summary(os.path.basename(d.rstrip("/")), load(d), args.top)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)

    def profiling_options(p):
        p.add_argument("--mem", default="10G")
        p.add_argument("--no-scope", action="store_true")
        p.add_argument("--step-region", default="Domain::update")
        p.add_argument("--top", type=int, default=12)

    p = sub.add_parser("run")
    p.add_argument("--exe", default=os.path.join(REPO, "bin", "FLEKS.exe"))
    p.add_argument("--out", required=True)
    p.add_argument("--np", type=int, default=1)
    p.add_argument("--fleks-args", default="")
    p.add_argument("decks", nargs="+")
    profiling_options(p)
    p.set_defaults(func=run)
    p = sub.add_parser("cmd")
    p.add_argument("--dir", required=True)
    p.add_argument("command", nargs=argparse.REMAINDER)
    profiling_options(p)
    p.set_defaults(func=run_cmd)
    p = sub.add_parser("analyze")
    p.add_argument("dirs", nargs="+")
    p.add_argument("--step-region", default="Domain::update")
    p.add_argument("--top", type=int, default=12)
    p.set_defaults(func=reanalyze)
    p = sub.add_parser("show")
    p.add_argument("dirs", nargs="+")
    p.add_argument("--top", type=int, default=12)
    p.set_defaults(func=show)
    p = sub.add_parser("compare")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--top", type=int, default=12)
    p.set_defaults(func=compare)
    args = parser.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
