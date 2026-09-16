#!/usr/bin/env python3
"""Diff two TinyProfiler captures (master vs PR) and flag regressions.

Reads two documents produced by ``tests/profile_tests.py`` and compares them
region by region.  Three families are reported:

* **arena allocation counts and peak bytes** -- exact integers for a fixed
  problem, rank count and RNG seed (``profile_tests.py --verify`` checks this),
  so they are gated strictly;
* **RSS** from the FLEKS load-balance report -- covers the whole process
  including ``std::vector`` / ``operator new`` traffic the arena profiler
  cannot see, but is not bit-identical between runs, so it is gated with a
  tolerance band (``--rss-tol`` / ``--rss-abs-tol``);
* **wall-clock timings** -- not reproducible on shared machines, so reported
  as warnings unless ``--gate-timing`` is passed.

Usage::

    python3 tests/compare_profiles.py master.json pr.json --out diff.md
    python3 tests/compare_profiles.py master.json pr.json --warn-only
"""
import argparse
import json
import sys

FAIL = "FAIL"
WARN = "WARN"
INFO = "INFO"

# Metrics compared per region.
MEMORY_COUNT = "nalloc"
MEMORY_PEAK = "maxmem_max"
MEMORY_LEAK = "curmem_max"
TIMING_METRICS = ("excl_max", "incl_max")


def _fmt(value, unit=""):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:,.4g}{unit}"
    return f"{value:,}{unit}"


def _fmt_bytes(value):
    if value is None:
        return "-"
    scaled = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(scaled) < 1024 or unit == "TiB":
            return f"{scaled:,.4g} {unit}"
        scaled /= 1024.0
    return f"{scaled:,.4g} TiB"


class Findings:
    """Ordered collection of comparison results."""

    def __init__(self):
        self.items = []

    def add(self, severity, test, region, metric, base, cand, note="",
            arena="", unit=""):
        delta = None
        if base is not None and cand is not None:
            delta = cand - base
        rel = None
        if base not in (None, 0) and delta is not None:
            rel = delta / abs(base)
        self.items.append({
            "severity": severity, "test": test, "arena": arena,
            "region": region, "metric": metric, "base": base, "cand": cand,
            "delta": delta, "delta_pct": None if rel is None else 100.0 * rel,
            "note": note, "unit": unit,
        })

    def by_severity(self, severity):
        return [i for i in self.items if i["severity"] == severity]

    def failures(self):
        return self.by_severity(FAIL)


def _ncalls(profile, region):
    return (profile.get("timing") or {}).get(region, {}).get("ncalls")


def _compare_load_balance(base, cand, key, cfg, findings):
    """Compare the RSS reports recorded in the test records.

    RSS is gated with a tolerance band rather than for equality: unlike the
    arena counters it is not bit-identical between runs (page granularity and
    allocator behaviour), by up to ~0.5 MB on the decks used here.  Measure the
    spread with ``profile_tests.py --verify`` and tune --rss-tol /
    --rss-abs-tol from it.  The allowed change is
    ``max(abs_tol, rel_tol * baseline)``, so the absolute floor keeps the
    tolerance from collapsing to nothing on small baselines.
    """
    bl_b = base.get("load_balance") or {}
    bl_c = cand.get("load_balance") or {}
    if not bl_b and not bl_c:
        return

    if bl_b.get("n_tables") != bl_c.get("n_tables"):
        findings.add(INFO, key, "-", "reports", bl_b.get("n_tables"),
                     bl_c.get("n_tables"), "number of RSS reports changed",
                     arena="RSS")

    for name in ("first", "last"):
        tb, tc = bl_b.get(name), bl_c.get(name)
        if not isinstance(tb, dict) or not isinstance(tc, dict):
            continue

        rss_b, rss_c = tb.get("Memory(MB)") or {}, tc.get("Memory(MB)") or {}
        for column in ("min", "avg", "max"):
            vb, vc = rss_b.get(column), rss_c.get(column)
            if vb is None or vc is None:
                continue
            allowed = max(cfg.rss_abs_tol, cfg.rss_tol * abs(vb))
            # Only the per-rank maximum is gated; min/avg are context.
            severity = FAIL if column == "max" else INFO
            if vc > vb + allowed:
                findings.add(severity, key, f"Memory(MB) [{name}]", f"rss_{column}",
                             vb, vc, "RSS grew", arena="RSS", unit="MB")
            elif vc < vb - allowed:
                findings.add(INFO, key, f"Memory(MB) [{name}]", f"rss_{column}",
                             vb, vc, "RSS shrank", arena="RSS", unit="MB")

        # Problem size, as context: a grid/particle change explains an RSS move
        # and means the two runs are not measuring the same thing.
        for label in ("Cells # of all levs", "Parts # of all levs"):
            vb = (tb.get(label) or {}).get("max")
            vc = (tc.get(label) or {}).get("max")
            if vb is not None and vc is not None and vb != vc:
                findings.add(INFO, key, label, "problem size", vb, vc,
                             "problem size changed", arena="RSS")


def _compare_timing(base, cand, key, cfg, findings):
    b = base.get("timing") or {}
    c = cand.get("timing") or {}
    for region in sorted(set(b) | set(c)):
        rb, rc = b.get(region), c.get(region)
        if rb is None:
            findings.add(INFO, key, region, "timing", None, None,
                         "region added", unit="s")
            continue
        if rc is None:
            findings.add(INFO, key, region, "timing", None, None,
                         "region removed", unit="s")
            continue

        nb, nc = rb.get("ncalls"), rc.get("ncalls")
        if nb != nc:
            findings.add(INFO, key, region, "ncalls", nb, nc,
                         "call count changed; timings not comparable")

        for metric in TIMING_METRICS:
            vb, vc = rb.get(metric), rc.get(metric)
            if vb is None or vc is None:
                continue
            # Ignore regions below the noise floor; on a shared runner a
            # sub-millisecond region can move by 50 % without meaning anything.
            if max(abs(vb), abs(vc)) < cfg.min_time_s:
                continue
            if vb <= 0:
                continue
            rel = (vc - vb) / vb
            if rel > cfg.time_tol:
                severity = FAIL if cfg.gate_timing else WARN
                findings.add(severity, key, region, metric, vb, vc,
                             f"+{100 * rel:.1f}% slower", unit="s")
            elif rel < -cfg.time_tol:
                findings.add(INFO, key, region, metric, vb, vc,
                             f"{100 * rel:.1f}% faster", unit="s")


def _compare_memory(base, cand, key, cfg, findings):
    bm = base.get("memory") or {}
    cm = cand.get("memory") or {}
    for arena in sorted(set(bm) | set(cm)):
        rb_a, rc_a = bm.get(arena, {}), cm.get(arena, {})
        for region in sorted(set(rb_a) | set(rc_a)):
            rb, rc = rb_a.get(region), rc_a.get(region)
            if rb is None:
                findings.add(INFO, key, region, MEMORY_PEAK, None, None,
                             "region added", arena=arena)
                continue
            if rc is None:
                findings.add(INFO, key, region, MEMORY_PEAK, None, None,
                             "region removed", arena=arena)
                continue

            nb, nc = _ncalls(base, region), _ncalls(cand, region)
            if nb is not None and nc is not None and nb > 0 and nc > 0 and nb != nc:
                # More calls legitimately means more allocations; normalise so
                # that "the refactor allocated something extra per call" is what
                # gets flagged.
                metric, vb, vc = "nalloc/call", rb.get("nalloc", 0) / nb, \
                    rc.get("nalloc", 0) / nc
                unit = ""
            else:
                metric, vb, vc = MEMORY_COUNT, rb.get("nalloc"), rc.get("nalloc")
                unit = ""

            if vb is not None and vc is not None:
                allowed = vb * (1.0 + cfg.alloc_tol)
                if vc > allowed:
                    findings.add(FAIL, key, region, metric, vb, vc,
                                 "more allocations per call", arena=arena,
                                 unit=unit)
                elif vc < vb * (1.0 - cfg.alloc_tol):
                    findings.add(INFO, key, region, metric, vb, vc,
                                 "fewer allocations per call", arena=arena,
                                 unit=unit)

            vb, vc = rb.get(MEMORY_PEAK), rc.get(MEMORY_PEAK)
            if vb is not None and vc is not None and vc > vb * (1.0 + cfg.mem_tol):
                findings.add(FAIL, key, region, MEMORY_PEAK, vb, vc,
                             "peak bytes grew", arena=arena)

            vb, vc = rb.get(MEMORY_LEAK), rc.get(MEMORY_LEAK)
            if vb is not None and vc is not None and vc > max(vb, 0):
                findings.add(FAIL, key, region, MEMORY_LEAK, vb, vc,
                             "memory not freed at finalize", arena=arena)


def compare(baseline, candidate, cfg):
    """Compare two profile documents and return a Findings collection."""
    findings = Findings()
    bmeta, cmeta = baseline.get("meta", {}), candidate.get("meta", {})
    btests, ctests = baseline.get("tests", {}), candidate.get("tests", {})

    if bmeta.get("amrex_dim") != cmeta.get("amrex_dim"):
        findings.add(FAIL, "-", "-", "amrex_dim", bmeta.get("amrex_dim"),
                     cmeta.get("amrex_dim"), "profiles are not comparable")

    for key in sorted(set(btests) | set(ctests)):
        rb, rc = btests.get(key), ctests.get(key)
        if rb is None:
            findings.add(INFO, key, "-", "test", None, None, "test added")
            continue
        if rc is None:
            findings.add(WARN, key, "-", "test", None, None, "test not run")
            continue
        if "error" in rb or "error" in rc:
            findings.add(WARN, key, "-", "test", None, None,
                         "one side produced no profile")
            continue
        if rb.get("nprocs") != rc.get("nprocs"):
            findings.add(WARN, key, "-", "nprocs", rb.get("nprocs"),
                         rc.get("nprocs"), "rank counts differ")
            continue

        _compare_timing(rb["profile"], rc["profile"], key, cfg, findings)
        _compare_memory(rb["profile"], rc["profile"], key, cfg, findings)
        _compare_load_balance(rb, rc, key, cfg, findings)

    return findings


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def _table(rows, show_arena=False):
    header = ["| Test", "Region", "Metric", "master", "PR", "Δ", "Note |"]
    if show_arena:
        header.insert(1, "Arena")
    lines = [" | ".join(header), "| :--- | :--- | :--- | ---: | ---: | ---: | :--- |"
             if not show_arena else
             "| :--- | :--- | :--- | :--- | ---: | ---: | ---: | :--- |"]
    for row in rows:
        cells = [f"| `{row['test']}`"]
        if show_arena:
            cells.append(f"`{row['arena']}`" if row["arena"] else "")
        cells.append(f"`{row['region']}`")
        cells.append(f"`{row['metric']}`")
        if row["metric"] in (MEMORY_PEAK, MEMORY_LEAK):
            cells.append(_fmt_bytes(row["base"]))
            cells.append(_fmt_bytes(row["cand"]))
            cells.append(_fmt_bytes(row["delta"]) if row["delta"] else "-")
        elif row.get("unit") == "MB":
            cells.append(_fmt(row["base"], " MB"))
            cells.append(_fmt(row["cand"], " MB"))
            cells.append(f"{row['delta']:+,.1f} MB" if row["delta"] else "-")
        else:
            unit = row.get("unit", "")
            cells.append(_fmt(row["base"], unit))
            cells.append(_fmt(row["cand"], unit))
            cells.append(_fmt(row["delta"], unit) if row["delta"] is not None else "-")
        cells.append(f"{row['note']} |")
        lines.append(" | ".join(cells))
    return "\n".join(lines)


def format_markdown(baseline, candidate, findings, title):
    bmeta, cmeta = baseline.get("meta", {}), candidate.get("meta", {})
    out = [f"### {title}", ""]
    out.append(f"* master: `{bmeta.get('commit', '?')[:12]}` "
               f"({bmeta.get('ref', '?')}) — {bmeta.get('commit_subject', '')}")
    out.append(f"* PR: `{cmeta.get('commit', '?')[:12]}` "
               f"({cmeta.get('ref', '?')}) — {cmeta.get('commit_subject', '')}")
    out.append("")

    failures = findings.failures()
    warnings = findings.by_severity(WARN)
    infos = findings.by_severity(INFO)

    if not failures and not warnings:
        out.append("✅ **No profiler regressions.**")
    else:
        out.append(f"**{len(failures)} regression(s), {len(warnings)} warning(s)** "
                   f"across {len(set(i['test'] for i in findings.items))} tests.")

    if failures:
        out += ["", "#### 🔴 Regressions (allocations, peak bytes, RSS)",
                "", _table(failures, show_arena=True)]
    if warnings:
        out += ["", "#### 🟡 Warnings (timing, not gated)", "",
                _table(warnings)]
    improvements = [i for i in infos if "faster" in i["note"]
                    or "fewer" in i["note"]]
    if improvements:
        out += ["", "#### 🟢 Improvements", "",
                _table(sorted(improvements, key=lambda r: r["delta"] or 0)[:10])]
    structural = [i for i in infos if i["note"] in
                  ("region added", "region removed", "test added",
                   "test not run", "call count changed; timings not comparable",
                   "number of RSS reports changed", "problem size changed")]
    if structural:
        out += ["", "#### ⚪ Added / removed / not comparable", "",
                _table(structural[:40], show_arena=True)]
    out.append("")
    return "\n".join(out)


def main():
    parser = argparse.ArgumentParser(
        description="Diff two TinyProfiler captures and flag regressions.")
    parser.add_argument("baseline", help="reference profile JSON (master)")
    parser.add_argument("candidate", help="candidate profile JSON (PR)")
    parser.add_argument("--out", help="write the markdown report here")
    parser.add_argument("--json", dest="json_out", help="write findings as JSON")
    parser.add_argument("--mem-tol", type=float, default=0.02,
                        help="allowed relative growth of peak bytes (default 0.02)")
    parser.add_argument("--rss-tol", type=float, default=0.02,
                        help="allowed relative growth of RSS (default 0.02)")
    parser.add_argument("--rss-abs-tol", type=float, default=2.0,
                        help="allowed absolute RSS growth in MB; overrides "
                             "--rss-tol for small baselines (default 2.0)")
    parser.add_argument("--alloc-tol", type=float, default=0.0,
                        help="allowed relative growth of allocations per call "
                             "(default 0.0 = any increase fails)")
    parser.add_argument("--time-tol", type=float, default=0.25,
                        help="relative timing change to report (default 0.25)")
    parser.add_argument("--min-time-s", type=float, default=0.01,
                        help="ignore timing regions smaller than this (default 0.01)")
    parser.add_argument("--gate-timing", action="store_true",
                        help="make timing regressions fail instead of warn")
    parser.add_argument("--warn-only", action="store_true",
                        help="never exit non-zero")
    parser.add_argument("--title", default="🔬 Profiler regression report")
    args = parser.parse_args()

    with open(args.baseline) as handle:
        baseline = json.load(handle)
    with open(args.candidate) as handle:
        candidate = json.load(handle)

    findings = compare(baseline, candidate, args)
    report = format_markdown(baseline, candidate, findings, args.title)

    if args.out:
        with open(args.out, "w") as handle:
            handle.write(report)
        print(f"wrote {args.out}")
    else:
        print(report)

    if args.json_out:
        with open(args.json_out, "w") as handle:
            json.dump({"baseline": baseline.get("meta", {}),
                       "candidate": candidate.get("meta", {}),
                       "findings": findings.items}, handle, indent=2)
        print(f"wrote {args.json_out}")

    if args.warn_only:
        return 0
    return 1 if findings.failures() else 0


if __name__ == "__main__":
    sys.exit(main())
