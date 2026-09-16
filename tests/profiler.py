#!/usr/bin/env python3
"""Parsers for the FLEKS standalone run report.

Two independent reports are parsed:

* the **AMReX TinyProfiler** report (``parse_tinyprofiler``), printed at
  ``amrex::Finalize()``: per-region timings, plus allocation count and peak
  bytes per region for each profiled arena;
* the **FLEKS load-balance report** (``parse_load_balance``), printed by
  ``Pic::report_load_balance()``: blocks / cells / particles per level and the
  resident set size (RSS) in MB as a min/avg/max across ranks.

Tables are parsed by header rather than by column position, and values are
consumed right-to-left: the column set varies with the run (single-rank runs
print one value per metric, ``Nfree``/``CurrentMem`` only appear when memory is
still allocated at finalize) and memory values are two tokens (``4755 KiB``),
so region names must be allowed to contain spaces.

Usage::

    python3 tests/profiler.py prof.txt              # timing + arena memory
    python3 tests/profiler.py prof.txt -o prof.json # write JSON
    python3 tests/profiler.py run.log --load-balance --step 10   # RSS series
    python3 tests/profiler.py --self-test           # check the parsers
"""
import argparse
import json
import os
import re
import sys

# Byte units used by amrex::TinyProfiler's mem_to_string().  Note that it
# divides by 1024 repeatedly, so "8192 KiB" really is 8 MiB of memory.
_BYTE_UNITS = {"B": 1, "KiB": 1024, "MiB": 1024**2, "GiB": 1024**3, "TiB": 1024**4}

_HLINE_RE = re.compile(r"^-{10,}\s*$")
_TOTAL_RE = re.compile(
    r"TinyProfiler total time across processes.*?:\s*"
    r"(?P<min>[\d.eE+-]+)\s*\.\.\.\s*"
    r"(?P<avg>[\d.eE+-]+)\s*\.\.\.\s*"
    r"(?P<max>[\d.eE+-]+)")
_USAGE_RE = re.compile(r"^(?P<name>.+?)\s+Usage:\s*$")

# Column labels emitted by AMReX, longest first so that "AvgMem min" is
# preferred over "AvgMem" when the header is split.
_LABELS = [
    "CurrentMem max",
    "AvgMem min", "AvgMem avg", "AvgMem max",
    "MaxMem min", "MaxMem avg", "MaxMem max",
    "Excl. Min", "Excl. Avg", "Excl. Max",
    "Incl. Min", "Incl. Avg", "Incl. Max",
    "CurrentMem", "AvgMem", "MaxMem", "NCalls", "Max %", "Nfree", "Nalloc",
    "Name",
]

# Columns holding a "<number> <unit>" memory value rather than a plain number.
_MEMORY_COLUMNS = {c for c in _LABELS if "Mem" in c}

# Header label -> canonical key.  With a single MPI rank AMReX prints one value
# per metric instead of min/avg/max; that value is the per-rank maximum, so it
# is mapped onto the "_max" key to keep 1-rank and N-rank profiles comparable.
_KEY = {
    "NCalls": "ncalls",
    "Excl. Min": "excl_min", "Excl. Avg": "excl_avg", "Excl. Max": "excl_max",
    "Incl. Min": "incl_min", "Incl. Avg": "incl_avg", "Incl. Max": "incl_max",
    "Max %": "max_pct",
    "Nalloc": "nalloc", "Nfree": "nfree",
    "AvgMem": "avgmem_max",
    "AvgMem min": "avgmem_min", "AvgMem avg": "avgmem_avg",
    "AvgMem max": "avgmem_max",
    "MaxMem": "maxmem_max",
    "MaxMem min": "maxmem_min", "MaxMem avg": "maxmem_avg",
    "MaxMem max": "maxmem_max",
    "CurrentMem": "curmem_max", "CurrentMem max": "curmem_max",
}


def _split_header(line):
    """Split a TinyProfiler header line into its (possibly multi-word) labels."""
    columns, i, n = [], 0, len(line)
    while i < n:
        if line[i].isspace():
            i += 1
            continue
        for label in _LABELS:
            if line.startswith(label, i):
                columns.append(label)
                i += len(label)
                break
        else:
            match = re.match(r"\S+", line[i:])
            columns.append(match.group(0))
            i += match.end()
    return columns


def _to_number(token):
    return float(token.rstrip("%"))


def _to_bytes(number, unit):
    return int(float(number) * _BYTE_UNITS[unit])


def _parse_row(line, columns):
    """Return (region name, {canonical key: value}) for one data row."""
    tokens = line.split()
    values, i = {}, len(tokens) - 1
    for column in reversed(columns[1:]):
        if column in _MEMORY_COLUMNS:
            if i < 1:
                return None
            values[_KEY[column]] = _to_bytes(tokens[i - 1], tokens[i])
            i -= 2
        else:
            if i < 0:
                return None
            values[_KEY[column]] = _to_number(tokens[i])
            i -= 1
    return " ".join(tokens[:i + 1]), values


def _iter_tables(lines):
    """Yield (title, columns, data rows) for each hline-delimited table.

    AMReX prints ``hline / header / hline / rows... / hline``.  *title* is the
    nearest preceding "<name> Usage:" line, or "" for the timing tables.
    """
    i, n = 0, len(lines)
    while i < n:
        if not _HLINE_RE.match(lines[i]):
            i += 1
            continue
        header_idx = i + 1
        if header_idx + 1 >= n or not _HLINE_RE.match(lines[header_idx + 1]):
            i += 1
            continue
        columns = _split_header(lines[header_idx])
        if columns[0] != "Name":
            i += 1
            continue

        rows = []
        j = header_idx + 2
        while j < n and not _HLINE_RE.match(lines[j]):
            rows.append(lines[j])
            j += 1

        title = ""
        for back in range(i - 1, max(-1, i - 4), -1):
            match = _USAGE_RE.match(lines[back])
            if match:
                title = match.group("name")
                break
            if lines[back].strip():
                break

        yield title, columns, rows
        i = j + 1 if j < n else n


def parse_tinyprofiler(text):
    """Parse TinyProfiler output into a comparable dictionary.

    Returns::

        {"total_time_s": {"min", "avg", "max"},
         "timing": {region: {"ncalls", "excl_*", "incl_*", "max_pct"}},
         "memory": {arena: {region: {"nalloc", "nfree", "*mem_*"}}}}

    ``timing`` and ``memory`` are empty (not an error) when profiling was
    disabled, e.g. with ``tiny_profiler.enabled=0``.
    """
    lines = text.splitlines()
    result = {"total_time_s": {}, "timing": {}, "memory": {}}

    match = _TOTAL_RE.search(text)
    if match:
        result["total_time_s"] = {k: float(v) for k, v in match.groupdict().items()}

    for title, columns, rows in _iter_tables(lines):
        if "Nalloc" in columns:
            arena = result["memory"].setdefault(title or "Unknown", {})
            target = arena
        else:
            target = result["timing"]

        for row in rows:
            if not row.strip():
                continue
            parsed = _parse_row(row, columns)
            if parsed is None:
                continue
            name, values = parsed
            target.setdefault(name, {}).update(values)

    return result


def parse_tinyprofiler_file(path):
    with open(path, "r") as handle:
        return parse_tinyprofiler(handle.read())


# ---------------------------------------------------------------------------
# FLEKS load-balance report (RSS)
# ---------------------------------------------------------------------------
_LOAD_BALANCE_START = re.compile(r"^=+\s*Load balance report\s*=+\s*$")
_RULE_RE = re.compile(r"^[-=]+\s*$")


def parse_load_balance(text):
    """Parse every FLEKS load-balance table in a run log.

    ``Pic::report_load_balance()`` prints one table per report step::

        ===============================Load balance report=============================
        |     Value          |      Min      |     Avg      |      Max     |where(max)|
        |Cells  # of all levs|          64.0 |         64.0 |         64.0 |         0|
        |Memory(MB)          |          54.4 |         54.4 |         54.4 |         0|
        ===============================================================================

    Returns a list of ``{label: {"min", "avg", "max", "where"}}`` in order of
    appearance; ``Memory(MB)`` is the resident set size of each rank in MB and
    the three columns are min/avg/max across MPI ranks.

    Note that the RSS series must **not** be obtained by enabling ``#MEMORY``:
    that command frees arena memory, flushes the tile cache, calls
    ``ShrinkToFit()`` on every particle container and ``malloc_trim(0)`` every
    ``dnMemory`` cycles (``Pic::free_memory()``), which perturbs both the RSS
    trajectory and the allocation counts we are measuring.  The report is
    printed anyway whenever ``doReport`` is set, so no parameter is needed.
    """
    tables, current = [], None
    for line in text.splitlines():
        if _LOAD_BALANCE_START.match(line):
            current = {}
            continue
        if current is None:
            continue
        if _RULE_RE.match(line) and "=" in line:
            # The closing "===" line terminates the table; the "-" rules
            # between sections are just visual separators.
            tables.append(current)
            current = None
            continue
        if _RULE_RE.match(line) or "|" not in line:
            continue

        fields = [f.strip() for f in line.strip().strip("|").split("|")]
        if len(fields) < 5 or fields[0] == "Value":
            continue
        label = re.sub(r"\s+", " ", fields[0])
        try:
            current[label] = {"min": float(fields[1]), "avg": float(fields[2]),
                              "max": float(fields[3]), "where": int(fields[4])}
        except ValueError:
            continue

    return tables


def parse_load_balance_file(path):
    with open(path, "r") as handle:
        return parse_load_balance(handle.read())


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def _fmt_bytes(nbytes):
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if nbytes < 1024 or unit == "TiB":
            return f"{nbytes:.0f} {unit}" if unit == "B" else f"{nbytes:.1f} {unit}"
        nbytes /= 1024.0
    return f"{nbytes:.1f} TiB"


def format_summary(profile, top=None, memory_only=False):
    """Render a profile as a plain-text report."""
    out = []
    total = profile.get("total_time_s") or {}
    if total:
        out.append("Total wall time [min...avg...max]: "
                   f"{total.get('min', 0):.4f} ... {total.get('avg', 0):.4f} ... "
                   f"{total.get('max', 0):.4f} s")

    if not memory_only:
        timing = profile.get("timing") or {}
        out.append("")
        out.append(f"Timing ({len(timing)} regions), top by exclusive max:")
        out.append(f"  {'region':<45} {'ncalls':>9} {'excl_max':>10} {'incl_max':>10}")
        rows = sorted(timing.items(), key=lambda kv: kv[1].get("excl_max", 0),
                      reverse=True)
        for name, values in rows[:top] if top else rows:
            out.append(f"  {name:<45} {int(values.get('ncalls', 0)):>9} "
                       f"{values.get('excl_max', 0):>10.5f} "
                       f"{values.get('incl_max', 0):>10.5f}")

    for arena, regions in (profile.get("memory") or {}).items():
        out.append("")
        out.append(f"{arena} ({len(regions)} regions), top by peak bytes:")
        out.append(f"  {'region':<45} {'nalloc':>9} {'nfree':>9} {'maxmem_max':>12}")
        rows = sorted(regions.items(), key=lambda kv: kv[1].get("maxmem_max", 0),
                      reverse=True)
        for name, values in rows[:top] if top else rows:
            # nfree is only reported when memory is still allocated at
            # finalize; do not render a missing value as a real zero.
            nfree = values.get("nfree")
            out.append(f"  {name:<45} {int(values.get('nalloc', 0)):>9} "
                       f"{'-' if nfree is None else int(nfree):>9} "
                       f"{_fmt_bytes(values.get('maxmem_max', 0)):>12}")

    return "\n".join(out)


def format_load_balance(tables, step=1):
    """Render a list of load-balance tables as a plain-text report."""
    if not tables:
        return "no load-balance report found"
    out = [f"{len(tables)} load-balance report(s)",
           f"  {'report':>6} {'RSS min':>9} {'RSS avg':>9} {'RSS max':>9} "
           f"{'cells':>10} {'parts':>12}"]
    for i, table in enumerate(tables):
        if i % step and i != len(tables) - 1:
            continue
        rss = table.get("Memory(MB)", {})
        cells = table.get("Cells # of all levs", {})
        parts = table.get("Parts # of all levs", {})
        out.append(f"  {i:>6} {rss.get('min', 0):>9.1f} {rss.get('avg', 0):>9.1f} "
                   f"{rss.get('max', 0):>9.1f} {cells.get('max', 0):>10.1f} "
                   f"{parts.get('max', 0):>12.1f}")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Self-test
#
# A miniature run report covering the report variants that matter: the two
# timing tables, a memory table with unit suffixes, a single-rank load-balance
# table, a multi-rank one where min/avg/max differ, and one with the optional
# Nfree/CurrentMem columns.
# ---------------------------------------------------------------------------
_SELF_TEST_REPORT = """\
TinyProfiler total time across processes [min...avg...max]: 0.5 ... 0.5 ... 0.5

-------------------------------------------------------------------------------------------
Name                                        NCalls  Excl. Min  Excl. Avg  Excl. Max   Max %
-------------------------------------------------------------------------------------------
Pts::calc_mass_matrix                          126    0.07161    0.07161    0.07161  18.43%
Pts::charged_particle_mover                    126    0.01967    0.01967    0.01967   5.06%
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
Name                                        NCalls  Incl. Min  Incl. Avg  Incl. Max   Max %
-------------------------------------------------------------------------------------------
Pts::calc_mass_matrix                          126     0.0717     0.0717     0.0717  18.43%
-------------------------------------------------------------------------------------------

Cpu Memory Usage:
---------------------------------------------------
Name                     Nalloc    AvgMem    MaxMem
---------------------------------------------------
Grid::regrid_base            50  4755 KiB  4859 KiB
Pic::calculate_phi           63     0   B   512   B
---------------------------------------------------

Async Memory Usage:
-------------------------------------------------------------------------
Name      Nalloc  Nfree  AvgMem min  AvgMem avg  AvgMem max  CurrentMem max
-------------------------------------------------------------------------
Pts::merge    4      3       1 MiB       2 MiB       3 MiB          64 KiB
-------------------------------------------------------------------------

==== FLEKS1:  Cycle 1 ====

===============================Load balance report=============================
|     Value          |      Min      |     Avg      |      Max     |where(max)|
|Cells  # of all levs|          64.0 |         64.0 |         64.0 |         0|
|Memory(MB)          |          54.4 |         54.4 |         54.4 |         0|
===============================================================================

===============================Load balance report=============================
|     Value          |      Min      |     Avg      |      Max     |where(max)|
|Cells  # of all levs|           0.0 |         32.0 |         64.0 |         0|
|Memory(MB)          |          42.2 |         48.5 |         54.7 |         0|
===============================================================================
"""


def self_test():
    """Parse _SELF_TEST_REPORT and check the values it encodes."""
    checks, failures = [], []

    def check(label, ok, detail=""):
        checks.append((label, ok, detail))
        if not ok:
            failures.append(f"{label}: {detail}")

    profile = parse_tinyprofiler(_SELF_TEST_REPORT)
    total = profile.get("total_time_s", {})
    check("total time", total.get("max") == 0.5, f"got {total}")

    timing = profile.get("timing", {})
    mover = timing.get("Pts::charged_particle_mover", {})
    check("timing ncalls", mover.get("ncalls") == 126, f"got {mover.get('ncalls')}")
    check("timing excl and incl merged",
          mover.get("excl_max") == 0.01967 and mover.get("incl_max") is None,
          f"got {mover}")
    check("max_pct without the % sign",
          timing.get("Pts::calc_mass_matrix", {}).get("max_pct") == 18.43,
          f"got {timing.get('Pts::calc_mass_matrix')}")
    check("inclusive table merged into the same region",
          timing.get("Pts::calc_mass_matrix", {}).get("incl_max") == 0.0717,
          f"got {timing.get('Pts::calc_mass_matrix')}")

    cpu = profile.get("memory", {}).get("Cpu Memory", {})
    regrid = cpu.get("Grid::regrid_base", {})
    # 4859 KiB, as AMReX prints it.  Name and unit are two tokens.
    check("KiB converted to bytes",
          regrid.get("maxmem_max") == 4859 * 1024, f"got {regrid}")
    phi = cpu.get("Pic::calculate_phi", {})
    # Columns are assigned right-to-left, so the trailing "512   B" must land
    # on MaxMem and the "0   B" before it on AvgMem.
    check("trailing columns assigned in order",
          (phi.get("nalloc"), phi.get("avgmem_max"), phi.get("maxmem_max"))
          == (63, 0, 512), f"got {phi}")
    check("Nfree absent unless reported", "nfree" not in phi, f"got {phi}")

    async_ = profile.get("memory", {}).get("Async Memory", {})
    merge = async_.get("Pts::merge", {})
    check("multi-column memory table",
          (merge.get("nalloc"), merge.get("nfree"), merge.get("curmem_max"))
          == (4, 3, 64 * 1024), f"got {merge}")
    check("AvgMem max column", merge.get("avgmem_max") == 3 * 1024**2, f"got {merge}")

    tables = parse_load_balance(_SELF_TEST_REPORT)
    check("load-balance table count", len(tables) == 2, f"got {len(tables)}")
    if len(tables) == 2:
        check("1-rank RSS", tables[0].get("Memory(MB)", {}).get("max") == 54.4,
              f"got {tables[0].get('Memory(MB)')}")
        check("label normalised",
              tables[0].get("Cells # of all levs", {}).get("max") == 64.0,
              f"got {sorted(tables[0])}")
        rss = tables[1].get("Memory(MB)", {})
        check("2-rank RSS min/avg/max distinct",
              (rss.get("min"), rss.get("avg"), rss.get("max")) == (42.2, 48.5, 54.7),
              f"got {rss}")

    for label, ok, detail in checks:
        print(f"  [{'ok' if ok else 'FAIL'}] {label}" + (f" ({detail})" if not ok else ""))
    print(f"\n{len(checks) - len(failures)}/{len(checks)} checks passed")
    return 1 if failures else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", nargs="?", help="run report / prof.txt to parse")
    parser.add_argument("-o", "--output", help="write the parsed profile as JSON")
    parser.add_argument("--top", type=int, help="show only the N largest entries")
    parser.add_argument("--memory-only", action="store_true",
                        help="report only the memory tables")
    parser.add_argument("--load-balance", action="store_true",
                        help="parse the FLEKS load-balance (RSS) reports instead")
    parser.add_argument("--step", type=int, default=1,
                        help="with --load-balance, print every Nth report")
    parser.add_argument("--self-test", action="store_true",
                        help="parse the bundled sample and verify known values")
    args = parser.parse_args(argv)

    if args.self_test:
        return self_test()

    if not args.input:
        parser.error("an input file is required (or use --self-test)")

    if args.load_balance:
        print(format_load_balance(parse_load_balance_file(args.input),
                                  step=max(1, args.step)))
        return 0

    profile = parse_tinyprofiler_file(args.input)
    if args.output:
        with open(args.output, "w") as handle:
            json.dump(profile, handle, indent=2, sort_keys=True)
        print(f"wrote {args.output}")
    else:
        print(format_summary(profile, top=args.top, memory_only=args.memory_only))
    return 0


if __name__ == "__main__":
    sys.exit(main())
