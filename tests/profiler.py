#!/usr/bin/env python3
"""Parser for AMReX TinyProfiler output (timing + per-region memory).

FLEKS prints the AMReX TinyProfiler report at the end of every standalone run
(``amrex::Finalize()`` -> ``BL_TINY_PROFILE_FINALIZE()`` and
``BL_TINY_PROFILE_MEMORYFINALIZE()``).  It contains

* a total wall-clock line,
* two timing tables (exclusive, then inclusive) with columns
  ``Name NCalls Excl. Min Excl. Avg Excl. Max Max %``,
* one ``<Arena> Usage:`` table per profiled AMReX arena -- ``Cpu Memory`` on a
  CPU build -- with columns
  ``Name Nalloc [Nfree] AvgMem[...] MaxMem[...] [CurrentMem[...]]``
  where ``Name`` is the *enclosing* ``BL_PROFILE`` / ``timing_func`` region,
  i.e. the FLEKS function that made the allocation.

Two properties make positional parsing unsafe:

* the column set depends on the run -- single-rank runs print one value per
  metric instead of ``min``/``avg``/``max``, and ``Nfree`` / ``CurrentMem``
  only appear when memory is still allocated at finalize;
* memory values carry a unit suffix (``42   B``, ``4755 KiB``, ``19 MiB``)
  and are therefore two whitespace-separated tokens.

So every table is parsed **by header**, and values are consumed
right-to-left so that region names may contain spaces.

Usage::

    python3 tests/profiler.py prof.txt              # human-readable summary
    python3 tests/profiler.py prof.txt -o prof.json # write JSON
    python3 tests/profiler.py --self-test           # parse the bundled sample
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


# ---------------------------------------------------------------------------
# Self-test against the bundled sample
# ---------------------------------------------------------------------------
_SAMPLE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "profiler_samples", "tinyprofiler_beam.txt")


def self_test():
    """Parse the bundled beam sample and check a few known values."""
    if not os.path.isfile(_SAMPLE):
        print(f"sample not found: {_SAMPLE}")
        return 1

    profile = parse_tinyprofiler_file(_SAMPLE)
    checks, failures = [], []

    def check(label, ok, detail=""):
        checks.append((label, ok, detail))
        if not ok:
            failures.append(f"{label}: {detail}")

    total = profile.get("total_time_s", {})
    check("total time parsed", total.get("max", 0) > 0, f"got {total}")

    timing = profile.get("timing", {})
    check("timing regions parsed", len(timing) > 50, f"got {len(timing)}")
    mover = timing.get("Pts::charged_particle_mover")
    check("mover region present", mover is not None, "missing")
    if mover:
        check("mover ncalls == 126", int(mover.get("ncalls", -1)) == 126,
              f"got {mover.get('ncalls')}")
        check("mover has excl and incl",
              "excl_max" in mover and "incl_max" in mover, f"got {mover}")

    memory = profile.get("memory", {})
    check("Cpu Memory arena present", "Cpu Memory" in memory,
          f"got {sorted(memory)}")
    cpu = memory.get("Cpu Memory", {})
    regrid = cpu.get("Grid::regrid_base")
    check("Grid::regrid_base present", regrid is not None, "missing")
    if regrid:
        check("regrid nalloc == 50", int(regrid.get("nalloc", -1)) == 50,
              f"got {regrid.get('nalloc')}")
        # 4859 KiB, as printed by AMReX.
        check("regrid maxmem in bytes",
              regrid.get("maxmem_max", 0) == 4859 * 1024,
              f"got {regrid.get('maxmem_max')}")

    for label, ok, detail in checks:
        print(f"  [{'ok' if ok else 'FAIL'}] {label}" + (f" ({detail})" if not ok else ""))
    print(f"\n{len(checks) - len(failures)}/{len(checks)} checks passed")
    return 1 if failures else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", nargs="?", default=_SAMPLE,
                        help="TinyProfiler output file (default: bundled sample)")
    parser.add_argument("-o", "--output", help="write the parsed profile as JSON")
    parser.add_argument("--top", type=int, help="show only the N largest entries")
    parser.add_argument("--memory-only", action="store_true",
                        help="report only the memory tables")
    parser.add_argument("--self-test", action="store_true",
                        help="parse the bundled sample and verify known values")
    args = parser.parse_args(argv)

    if args.self_test:
        return self_test()

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
