#!/usr/bin/env python3
"""Diff the memory captured by two ``tests/profile_tests.py`` runs.

Two families are compared, both from the same profile document:

* **Arena allocations**, per ``BL_PROFILE`` region: the allocation count, the
  peak bytes held and anything left allocated at finalize.  These are exact
  integers for a fixed problem, rank count and RNG seed, so any increase is a
  real regression (``profile_tests.py --verify`` checks that they are
  bit-identical across runs).
* **RSS** from the FLEKS load-balance report: covers the whole process,
  including the ``std::vector`` / ``operator new`` traffic the arena tables
  cannot attribute, but it moves by a few tenths of a MB between runs of the
  same binary, so it gets a tolerance (``--rss-tol``).

Timing is deliberately **not** compared here.  It is captured in the same
document (``tests/profiler.py`` can show it) but two runs of identical code
already move individual regions by tens of percent, so it is not a useful
gate; ``tests/validate_performance.py`` is the tool for tracking speed.

Usage::

    python3 tests/compare_profiles.py master.json pr.json
    python3 tests/compare_profiles.py master.json pr.json --out diff.md
    python3 tests/compare_profiles.py master.json pr.json --rss-tol 5
"""
import argparse
import json
import sys

FAIL = "FAIL"
INFO = "INFO"

REPORT_TITLE = "🔬 Profiler memory report"


def _fmt_bytes(value):
    if value is None:
        return "-"
    scaled = abs(float(value))
    sign = "-" if value < 0 else ""
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if scaled < 1024 or unit == "TiB":
            return f"{sign}{scaled:,.4g} {unit}"
        scaled /= 1024.0
    return f"{sign}{scaled:,.4g} TiB"


def _cell(row, value):
    """Render one numeric cell according to the metric's kind."""
    if value is None:
        return "-"
    if row["kind"] == "bytes":
        return _fmt_bytes(value)
    if row["kind"] == "mb":
        return f"{value:,.1f} MB"
    if row["kind"] == "ratio":
        return f"{value:,.4g}"
    if row["kind"] == "number":
        return f"{value:,.0f}"
    return f"{value:,}"


class Findings:
    """Ordered collection of comparison results."""

    def __init__(self):
        self.items = []

    def add(self, severity, test, region, metric, base, cand, note="",
            arena="", kind="count"):
        delta = None
        if base is not None and cand is not None:
            delta = cand - base
        self.items.append({
            "severity": severity, "test": test, "arena": arena,
            "region": region, "metric": metric, "base": base, "cand": cand,
            "delta": delta, "note": note, "kind": kind,
        })

    def by_severity(self, severity):
        return [i for i in self.items if i["severity"] == severity]

    def failures(self):
        return self.by_severity(FAIL)


def _ncalls(document, region):
    return (document.get("timing") or {}).get(region, {}).get("ncalls")


def _compare_arena(base, cand, key, findings):
    """Compare the per-region arena tables.  Strict: the values are exact."""
    bm = base.get("memory") or {}
    cm = cand.get("memory") or {}
    for arena in sorted(set(bm) | set(cm)):
        rb_a, rc_a = bm.get(arena, {}), cm.get(arena, {})
        for region in sorted(set(rb_a) | set(rc_a)):
            rb, rc = rb_a.get(region), rc_a.get(region)
            if rb is None:
                findings.add(INFO, key, region, "nalloc", None, None,
                             "region added", arena=arena)
                continue
            if rc is None:
                findings.add(INFO, key, region, "nalloc", None, None,
                             "region removed", arena=arena)
                continue

            nb, nc = _ncalls(base, region), _ncalls(cand, region)
            vb, vc = rb.get("nalloc"), rc.get("nalloc")
            if vb is None or vc is None:
                continue
            if nb and nc and nb != nc:
                # More calls legitimately means more allocations; normalise so
                # that "the refactor allocated something extra per call" is what
                # gets flagged.
                findings.add(INFO, key, region, "ncalls", nb, nc,
                             "call count changed; allocations per call",
                             arena=arena, kind="number")
                vb, vc = vb / nb, vc / nc
                metric = "nalloc/call"
                kind = "ratio"
            else:
                metric = "nalloc"
                kind = "count"

            if vc > vb:
                findings.add(FAIL, key, region, metric, vb, vc,
                             "more allocations", arena=arena, kind=kind)
            elif vc < vb:
                findings.add(INFO, key, region, metric, vb, vc,
                             "fewer allocations", arena=arena, kind=kind)

            vb, vc = rb.get("maxmem_max"), rc.get("maxmem_max")
            if vb is not None and vc is not None:
                if vc > vb:
                    findings.add(FAIL, key, region, "maxmem_max", vb, vc,
                                 "peak bytes grew", arena=arena, kind="bytes")
                elif vc < vb:
                    findings.add(INFO, key, region, "maxmem_max", vb, vc,
                                 "peak bytes shrank", arena=arena, kind="bytes")

            # A missing column means "nothing was left allocated", so default
            # to 0: otherwise a leak newly introduced by the PR would be
            # skipped precisely because the baseline has no column for it.
            vb, vc = rb.get("curmem_max", 0), rc.get("curmem_max", 0)
            if vc > vb:
                findings.add(FAIL, key, region, "curmem_max", vb, vc,
                             "memory not freed at finalize", arena=arena,
                             kind="bytes")


def _compare_rss(base, cand, key, cfg, findings):
    """Compare the RSS reports recorded in the test records.

    RSS is gated with an absolute tolerance in MB (``--rss-tol``) rather than an
    equality check: it is not bit-identical between runs (page granularity and
    allocator behaviour).  Measure the spread with ``profile_tests.py --verify``
    and keep the tolerance comfortably above it.  Only the per-rank maximum is
    gated; min/avg are reported as context.
    """
    bl_b = base.get("load_balance") or {}
    bl_c = cand.get("load_balance") or {}
    if not bl_b and not bl_c:
        return

    if bl_b.get("n_tables") != bl_c.get("n_tables"):
        findings.add(INFO, key, "reports", "n_tables", bl_b.get("n_tables"),
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
            severity = FAIL if column == "max" else INFO
            if vc > vb + cfg.rss_tol:
                findings.add(severity, key, f"Memory(MB) [{name}]",
                             f"rss_{column}", vb, vc, "RSS grew", arena="RSS",
                             kind="mb")
            elif vc < vb - cfg.rss_tol:
                findings.add(INFO, key, f"Memory(MB) [{name}]", f"rss_{column}",
                             vb, vc, "RSS shrank", arena="RSS", kind="mb")

        # Problem size, as context: a grid/particle change explains an RSS move
        # and means the two runs are not measuring the same thing.
        for label in ("Cells # of all levs", "Parts # of all levs"):
            vb = (tb.get(label) or {}).get("max")
            vc = (tc.get(label) or {}).get("max")
            if vb is not None and vc is not None and vb != vc:
                findings.add(INFO, key, label, "problem size", vb, vc,
                             "problem size changed", arena="RSS", kind="number")


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
            findings.add(INFO, key, "-", "test", None, None, "test not run")
            continue
        if "error" in rb or "error" in rc:
            findings.add(INFO, key, "-", "test", None, None,
                         "one side produced no profile")
            continue
        if rb.get("nprocs") != rc.get("nprocs"):
            findings.add(INFO, key, "-", "nprocs", rb.get("nprocs"),
                         rc.get("nprocs"), "rank counts differ")
            continue

        _compare_arena(rb["profile"], rc["profile"], key, findings)
        _compare_rss(rb, rc, key, cfg, findings)

    return findings


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def _table(rows, show_arena=False):
    header = ["| Test", "Region", "Metric", "master", "PR", "Δ", "Note |"]
    if show_arena:
        header.insert(1, "Arena")
    separator = ("| :--- | :--- | :--- | :--- | ---: | ---: | ---: | :--- |"
                 if show_arena
                 else "| :--- | :--- | :--- | ---: | ---: | ---: | :--- |")
    lines = [" | ".join(header), separator]

    for row in rows:
        cells = [f"| `{row['test']}`"]
        if show_arena:
            cells.append(f"`{row['arena']}`" if row["arena"] else "")
        cells.append(f"`{row['region']}`")
        cells.append(f"`{row['metric']}`")

        if row["base"] is None and row["cand"] is None:
            cells += ["-", "-", "-"]
        else:
            delta = row["delta"]
            if delta is not None and row["kind"] in ("mb", "bytes", "ratio"):
                delta = round(delta, 3)
            cells.append(_cell(row, row["base"]))
            cells.append(_cell(row, row["cand"]))
            cells.append("-" if delta is None else _cell(row, delta))
        cells.append(f"{row['note']} |")
        lines.append(" | ".join(cells))

    return "\n".join(lines)


def format_markdown(baseline, candidate, findings):
    bmeta, cmeta = baseline.get("meta", {}), candidate.get("meta", {})
    out = [f"### {REPORT_TITLE}", ""]
    out.append(f"* master: `{str(bmeta.get('commit', '?'))[:12]}` "
               f"({bmeta.get('ref', '?')}) — {bmeta.get('commit_subject', '')}")
    out.append(f"* PR: `{str(cmeta.get('commit', '?'))[:12]}` "
               f"({cmeta.get('ref', '?')}) — {cmeta.get('commit_subject', '')}")
    out.append("")

    failures = findings.failures()
    others = findings.by_severity(INFO)

    if not failures:
        out.append("✅ **No memory regressions.**")
    else:
        out.append(f"**{len(failures)} memory regression(s)** "
                   f"across {len(set(i['test'] for i in failures))} test(s).")

    if failures:
        out += ["", "#### 🔴 Regressions (allocations, peak bytes, RSS)", "",
                _table(failures, show_arena=True)]
    if others:
        out += ["", "#### ⚪ Other changes", "",
                _table(sorted(others, key=lambda r: (r["test"], r["region"]))[:40],
                       show_arena=True)]
    out.append("")
    return "\n".join(out)


def main():
    parser = argparse.ArgumentParser(
        description="Diff the memory captured by two profile_tests.py runs.")
    parser.add_argument("baseline", help="reference profile JSON (master)")
    parser.add_argument("candidate", help="candidate profile JSON (PR)")
    parser.add_argument("--out", help="write the markdown report to this file")
    parser.add_argument("--rss-tol", type=float, default=2.0,
                        help="allowed RSS growth in MB (default 2.0); re-tune "
                             "with 'profile_tests.py --verify'")
    args = parser.parse_args()

    with open(args.baseline) as handle:
        baseline = json.load(handle)
    with open(args.candidate) as handle:
        candidate = json.load(handle)

    findings = compare(baseline, candidate, args)
    report = format_markdown(baseline, candidate, findings)

    if args.out:
        with open(args.out, "w") as handle:
            handle.write(report)
        print(f"wrote {args.out}")
    else:
        print(report)

    return 1 if findings.failures() else 0


if __name__ == "__main__":
    sys.exit(main())
