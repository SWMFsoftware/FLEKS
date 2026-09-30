#!/usr/bin/env python3
"""Require the two coupled test sessions to rebuild different fine grids."""

import argparse
import pathlib
import re
import sys


SESSION = re.compile(r"Starting Session\s+(\d+)")
FINE_GRID = re.compile(
    r"FLEKS0 fi:\s+iLev = 1\s+# of boxes =\s+(\d+)"
    r"\s+# of cells =\s+(\d+)"
)


def check_mesh_change(log):
    """Return both fine-grid sizes or raise for a missing mesh transition."""
    session = None
    grids = {}
    regridded = set()
    for line in log.splitlines():
        match = SESSION.search(line)
        if match:
            session = int(match.group(1))
        if session is None:
            continue
        if "FLEKS0: Domain::regrid is called" in line:
            regridded.add(session)
        match = FINE_GRID.search(line)
        if match:
            grids[session] = tuple(map(int, match.groups()))

    if 1 not in regridded or 2 not in regridded:
        raise ValueError("both sessions must call Domain::regrid")
    if 1 not in grids or 2 not in grids:
        raise ValueError("both sessions must report a level-1 fluid grid")
    if grids[1] == grids[2]:
        raise ValueError("level-1 fine-grid size did not change")
    return (*grids[1], *grids[2])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runlog", type=pathlib.Path)
    args = parser.parse_args()
    try:
        old_boxes, old_cells, new_boxes, new_cells = check_mesh_change(
            args.runlog.read_text()
        )
    except (OSError, ValueError) as error:
        print(f"Dynamic AMR mesh check failed: {error}", file=sys.stderr)
        return 1
    print(
        f"Dynamic AMR level 1: {old_boxes} boxes/{old_cells} cells -> "
        f"{new_boxes} boxes/{new_cells} cells"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
