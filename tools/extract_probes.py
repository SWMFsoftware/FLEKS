#!/usr/bin/env python3
"""Generic virtual-probe post-processor for FLEKS/SWMF .out and plotfile outputs.

Supports ASCII, binary ``real4``, and binary ``real8`` ``.out`` frames (plus
optional AMReX plotfiles via ``flekspy`` when available) and samples physical
variables onto:

* explicit point coordinates (``--point``, ``--points-file``),
* 1D line rakes (``--line x0,y0,z0 x1,y1,z1 N``),
* 2D planar arrays (``--plane ox,oy,oz ux,uy,uz vx,vy,vz Nu Nv``).

Results can be exported to tidy CSV (``.csv``) or compressed NumPy (``.npz``).
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import struct
import sys
from typing import Iterable, Sequence

import numpy as np


def _build_rectilinear_stencil(
    axes: Sequence[np.ndarray],
    query: np.ndarray,
) -> tuple[np.ndarray, list[tuple[tuple[np.ndarray, ...], np.ndarray]]]:
    """Precompute N-D multilinear stencil indices and corner weights once."""
    n_pts = query.shape[0]
    ndim = len(axes)
    if n_pts == 0:
        return np.zeros(0, dtype=bool), []

    lo_indices: list[np.ndarray] = []
    hi_indices: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    out_of_bounds = np.zeros(n_pts, dtype=bool)

    for dim, ax in enumerate(axes):
        q = query[:, dim]
        invalid = ~np.isfinite(q) | (q < ax[0]) | (q > ax[-1])
        out_of_bounds |= invalid

        if len(ax) == 1:
            zeros = np.zeros(n_pts, dtype=np.intp)
            lo_indices.append(zeros)
            hi_indices.append(zeros)
            weights.append(np.zeros(n_pts, dtype=np.float64))
            continue

        q_safe = np.where(invalid, ax[0], q)
        idx = np.clip(np.searchsorted(ax, q_safe, side="right") - 1, 0, len(ax) - 2)
        idx_hi = idx + 1
        span = ax[idx_hi] - ax[idx]
        w = np.where(span > 0.0, (q_safe - ax[idx]) / span, 0.0)
        lo_indices.append(idx)
        hi_indices.append(idx_hi)
        weights.append(w)

    corners: list[tuple[tuple[np.ndarray, ...], np.ndarray]] = []
    for mask in range(1 << ndim):
        corner_weight = np.ones(n_pts, dtype=np.float64)
        corner_idx: list[np.ndarray] = []
        for dim in range(ndim):
            if (mask >> dim) & 1:
                corner_idx.append(hi_indices[dim])
                corner_weight *= weights[dim]
            else:
                corner_idx.append(lo_indices[dim])
                corner_weight *= 1.0 - weights[dim]
        corners.append((tuple(corner_idx), corner_weight))

    return out_of_bounds, corners


def _apply_rectilinear_stencil(
    stencil: tuple[np.ndarray, list[tuple[tuple[np.ndarray, ...], np.ndarray]]],
    grid_vals: np.ndarray,
) -> np.ndarray:
    """Evaluate a precomputed rectilinear stencil on one variable array."""
    out_of_bounds, corners = stencil
    n_pts = out_of_bounds.shape[0]
    if n_pts == 0:
        return np.empty(0, dtype=np.float64)

    out = np.zeros(n_pts, dtype=np.float64)
    for corner_idx, corner_weight in corners:
        out += corner_weight * grid_vals[corner_idx]
    out[out_of_bounds] = np.nan
    return out


def _interp_rectilinear_nd(
    axes: Sequence[np.ndarray],
    grid_vals: np.ndarray,
    query: np.ndarray,
) -> np.ndarray:
    """Multilinear interpolation on an N-D rectilinear grid (NaN outside bounds)."""
    return _apply_rectilinear_stencil(
        _build_rectilinear_stencil(axes, query), grid_vals
    )


def _interp_idw_nd(
    coords: np.ndarray,
    values: np.ndarray,
    query: np.ndarray,
    k: int = 8,
) -> np.ndarray:
    """Pure-NumPy inverse-distance-weighted fallback when SciPy is unavailable."""
    vals = np.asarray(values, dtype=np.float64)
    squeeze_out = vals.ndim == 1
    vals_2d = vals[:, np.newaxis] if squeeze_out else vals
    n_query = query.shape[0]
    n_comp = vals_2d.shape[1]

    if n_query == 0:
        empty = np.empty((0, n_comp), dtype=np.float64)
        return empty[:, 0] if squeeze_out else empty
    if coords.shape[0] == 0:
        full = np.full((n_query, n_comp), np.nan, dtype=np.float64)
        return full[:, 0] if squeeze_out else full

    mins = np.min(coords, axis=0)
    maxs = np.max(coords, axis=0)
    tol = 1.0e-12 * np.maximum(maxs - mins, 1.0)
    out = np.full((n_query, n_comp), np.nan, dtype=np.float64)
    k_eff = min(k, coords.shape[0])

    for i in range(n_query):
        q = query[i]
        if not np.all(np.isfinite(q)):
            continue
        if np.any(q < mins - tol) or np.any(q > maxs + tol):
            continue
        d2 = np.sum((coords - q) ** 2, axis=1)
        if k_eff < coords.shape[0]:
            nn = np.argpartition(d2, k_eff - 1)[:k_eff]
        else:
            nn = np.arange(coords.shape[0])
        d2_nn = d2[nn]
        zero_mask = d2_nn == 0.0
        if np.any(zero_mask):
            out[i] = np.mean(vals_2d[nn[zero_mask]], axis=0)
        else:
            w = 1.0 / d2_nn
            out[i] = (w @ vals_2d[nn]) / np.sum(w)

    return out[:, 0] if squeeze_out else out


# ---------------------------------------------------------------------------
# Frame readers (.out ASCII / binary real4 / binary real8 / AMReX plotfile)
# ---------------------------------------------------------------------------
def _is_binary_out(path: str) -> bool:
    """Return True if *path* begins with a valid Fortran sequential record."""
    try:
        file_size = os.path.getsize(path)
    except OSError:
        return False
    if file_size < 12:
        return False
    with open(path, "rb") as f:
        hdr = f.read(4)
        if len(hdr) < 4:
            return False
        for endian in ("<", ">"):
            (nbytes,) = struct.unpack(endian + "i", hdr)
            if 1 <= nbytes <= min(10000, file_size - 8):
                f.seek(4 + nbytes)
                ftr = f.read(4)
                if len(ftr) == 4 and struct.unpack(endian + "i", ftr)[0] == nbytes:
                    return True
    return False


def _read_fortran_record(f, endian: str) -> bytes:
    """Read one sequential unformatted Fortran record."""
    hdr = f.read(4)
    if len(hdr) < 4:
        raise ValueError("Unexpected EOF reading Fortran record header")
    (nbytes,) = struct.unpack(endian + "i", hdr)
    payload = f.read(nbytes)
    if len(payload) < nbytes:
        raise ValueError(
            f"Truncated Fortran record payload: expected {nbytes}, got {len(payload)}"
        )
    ftr = f.read(4)
    if len(ftr) < 4:
        raise ValueError("Unexpected EOF reading Fortran record trailer")
    (nbytes2,) = struct.unpack(endian + "i", ftr)
    if nbytes2 != nbytes:
        raise ValueError(f"Fortran record marker mismatch: {nbytes} != {nbytes2}")
    return payload


def _skip_fortran_record(f, endian: str) -> None:
    """Skip over one sequential unformatted Fortran record without loading payload."""
    hdr = f.read(4)
    if len(hdr) < 4:
        raise ValueError("Unexpected EOF skipping Fortran record header")
    (nbytes,) = struct.unpack(endian + "i", hdr)
    f.seek(nbytes + 4, os.SEEK_CUR)


def _read_binary_out(
    path: str, var_names: Sequence[str] | None = None
) -> dict:
    """Parse a binary (real4 or real8) IDL/SWMF .out frame."""
    wanted_upper = (
        {v.strip().upper() for v in var_names} if var_names is not None else None
    )
    with open(path, "rb") as f:
        hdr = f.read(4)
        endian = "<"
        (nbytes_le,) = struct.unpack("<i", hdr)
        f.seek(0)
        if not (1 <= nbytes_le <= 10000):
            endian = ">"

        header_str = _read_fortran_record(f, endian).decode(
            "ascii", errors="replace"
        ).strip()
        meta_rec = _read_fortran_record(f, endian)
        if len(meta_rec) == 20:
            n_step, sim_time, n_dim_out, n_param, n_var = struct.unpack(
                endian + "ifiii", meta_rec
            )
            fdtype = endian + "f4"
        elif len(meta_rec) == 24:
            n_step, sim_time, n_dim_out, n_param, n_var = struct.unpack(
                endian + "idiii", meta_rec
            )
            fdtype = endian + "f8"
        else:
            raise ValueError(
                f"Unsupported metadata record size ({len(meta_rec)} bytes) in {path}"
            )

        n_dim = abs(n_dim_out)
        n_sizes = struct.unpack(
            endian + f"{n_dim}i", _read_fortran_record(f, endian)
        )
        n_point = n_sizes[0] if n_dim_out < 0 else math.prod(n_sizes)

        params = (
            np.frombuffer(_read_fortran_record(f, endian), dtype=fdtype)
            .astype(np.float64)
            .copy()
            if n_param > 0
            else np.array([], dtype=np.float64)
        )
        all_names = (
            _read_fortran_record(f, endian).decode("ascii", errors="replace").split()
        )
        coord_names = all_names[:n_dim]
        file_var_names = all_names[n_dim : n_dim + n_var]

        raw_coords = np.frombuffer(
            _read_fortran_record(f, endian), dtype=fdtype
        ).astype(np.float64)
        coords = raw_coords.reshape(n_dim, n_point).T.copy()

        fields = {}
        for vname in file_var_names:
            if wanted_upper is not None and vname.upper() not in wanted_upper:
                _skip_fortran_record(f, endian)
                continue
            arr = np.frombuffer(
                _read_fortran_record(f, endian), dtype=fdtype
            ).astype(np.float64)
            fields[vname] = arr.copy()

    return {
        "path": path,
        "format": "real4" if "4" in fdtype else "real8",
        "header": header_str,
        "nStep": int(n_step),
        "time": float(sim_time),
        "nDim": n_dim,
        "structured": n_dim_out > 0,
        "dims": tuple(int(s) for s in n_sizes),
        "params": params,
        "coord_names": coord_names,
        "var_names": file_var_names,
        "coords": coords,
        "fields": fields,
    }


def _read_ascii_out(path: str) -> dict:
    """Parse an ASCII SWMF/FLEKS .out frame."""
    with open(path, "r", encoding="latin-1") as f:
        lines = f.readlines()
    if len(lines) < 5:
        raise ValueError(f"ASCII .out file {path} has fewer than 5 lines")

    header_str = lines[0].strip()
    meta_toks = lines[1].split()
    n_step = int(meta_toks[0])
    sim_time = float(meta_toks[1])
    n_dim_out = int(meta_toks[2])
    n_param = int(meta_toks[3])
    n_var = int(meta_toks[4])
    n_dim = abs(n_dim_out)

    dim_toks = lines[2].split()
    n_sizes = tuple(int(tok) for tok in dim_toks[:n_dim])

    idx = 3
    params = np.array([], dtype=np.float64)
    if n_param > 0 and idx < len(lines):
        params = np.array(
            [float(tok) for tok in lines[idx].split()[:n_param]], dtype=np.float64
        )
        idx += 1

    # Skip blank or purely numeric parameter continuation lines until we reach
    # the alphabetic column-header line.
    while idx < len(lines):
        toks = lines[idx].split()
        if not toks:
            idx += 1
            continue
        is_numeric = True
        for tok in toks:
            try:
                float(tok)
            except ValueError:
                is_numeric = False
                break
        if not is_numeric:
            break
        idx += 1

    if idx >= len(lines):
        raise ValueError(f"Could not locate column-name header in {path}")

    all_names = lines[idx].split()
    coord_names = all_names[:n_dim]
    var_names = all_names[n_dim : n_dim + n_var]
    data_start = idx + 1

    data = np.loadtxt(lines[data_start:], dtype=np.float64)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    if data.size == 0:
        raise ValueError(f"ASCII .out file {path} contains no data rows")

    coords = data[:, :n_dim].copy()
    fields = {
        vname: data[:, n_dim + k].copy() for k, vname in enumerate(var_names)
    }

    return {
        "path": path,
        "format": "ascii",
        "header": header_str,
        "nStep": n_step,
        "time": sim_time,
        "nDim": n_dim,
        "structured": n_dim_out > 0,
        "dims": n_sizes,
        "params": params,
        "coord_names": coord_names,
        "var_names": var_names,
        "coords": coords,
        "fields": fields,
    }


def _read_amrex_plotfile(
    path: str, var_names: Sequence[str] | None = None
) -> dict:
    """Parse an AMReX plotfile directory via flekspy when available."""
    try:
        import flekspy
    except ImportError as exc:
        raise ImportError(
            f"Reading AMReX plotfile '{path}' requires flekspy, which is not "
            "installed in this Python environment."
        ) from exc

    ds = flekspy.load(path)
    sim_time = float(getattr(ds, "time", 0.0))
    n_step = int(getattr(ds, "nStep", 0))
    # Extract cell centers and field arrays from finest covering or raw grids.
    ad = ds.all_data() if hasattr(ds, "all_data") else ds
    coord_names = []
    coord_cols = []
    for axis in ("x", "y", "z"):
        if hasattr(ad, axis) or (isinstance(ad, dict) and axis in ad):
            coord_names.append(axis.upper())
            coord_cols.append(np.asarray(ad[axis], dtype=np.float64).ravel())
    if not coord_cols:
        raise ValueError(f"Could not extract coordinates from AMReX plotfile {path}")

    coords = np.column_stack(coord_cols)
    file_var_names = [
        v for v in getattr(ds, "field_list", []) if v.upper() not in ("X", "Y", "Z")
    ]
    wanted_upper = (
        {v.strip().upper() for v in var_names} if var_names is not None else None
    )
    fields = {
        v: np.asarray(ad[v], dtype=np.float64).ravel()
        for v in file_var_names
        if wanted_upper is None or v.upper() in wanted_upper
    }
    return {
        "path": path,
        "format": "amrex",
        "header": os.path.basename(path),
        "nStep": n_step,
        "time": sim_time,
        "nDim": len(coord_names),
        "structured": False,
        "dims": (coords.shape[0],),
        "params": np.array([], dtype=np.float64),
        "coord_names": coord_names,
        "var_names": file_var_names,
        "coords": coords,
        "fields": fields,
    }


def read_frame(path: str, var_names: Sequence[str] | None = None) -> dict:
    """Load one FLEKS output frame (.out ASCII/real4/real8 or AMReX plotfile)."""
    if os.path.isdir(path) and os.path.isfile(os.path.join(path, "Header")):
        return _read_amrex_plotfile(path, var_names=var_names)
    if _is_binary_out(path):
        return _read_binary_out(path, var_names=var_names)
    return _read_ascii_out(path)


# ---------------------------------------------------------------------------
# Probe geometry builders
# ---------------------------------------------------------------------------
def _parse_vec3(spec: str) -> np.ndarray:
    """Parse 'x,y' or 'x,y,z' into a float64 3-vector (z defaults to 0)."""
    parts = [float(p.strip()) for p in spec.split(",") if p.strip()]
    if len(parts) == 2:
        return np.array([parts[0], parts[1], 0.0], dtype=np.float64)
    if len(parts) == 3:
        return np.array(parts, dtype=np.float64)
    raise ValueError(f"Expected 2 or 3 comma-separated floats, got '{spec}'")


def make_line_probes(
    start: Sequence[float], end: Sequence[float], count: int
) -> np.ndarray:
    """Return *count* linearly spaced 3D probe coordinates along [start, end]."""
    if count < 1:
        raise ValueError(f"Line probe count must be >= 1, got {count}")
    p0 = np.asarray(start, dtype=np.float64).reshape(3)
    p1 = np.asarray(end, dtype=np.float64).reshape(3)
    t = np.linspace(0.0, 1.0, count)[:, np.newaxis]
    return p0[np.newaxis, :] + t * (p1 - p0)[np.newaxis, :]


def make_plane_probes(
    origin: Sequence[float],
    vec_u: Sequence[float],
    vec_v: Sequence[float],
    nu: int,
    nv: int,
) -> np.ndarray:
    """Return (nu*nv, 3) probe coordinates over origin + s*u + t*v."""
    if nu < 1 or nv < 1:
        raise ValueError(f"Plane grid counts must be >= 1, got ({nu}, {nv})")
    o = np.asarray(origin, dtype=np.float64).reshape(3)
    u = np.asarray(vec_u, dtype=np.float64).reshape(3)
    v = np.asarray(vec_v, dtype=np.float64).reshape(3)
    su = np.linspace(0.0, 1.0, nu)
    tv = np.linspace(0.0, 1.0, nv)
    uu, vv = np.meshgrid(su, tv, indexing="ij")
    pts = (
        o[np.newaxis, np.newaxis, :]
        + uu[:, :, np.newaxis] * u[np.newaxis, np.newaxis, :]
        + vv[:, :, np.newaxis] * v[np.newaxis, np.newaxis, :]
    )
    return pts.reshape(-1, 3)


def load_points_file(path: str) -> np.ndarray:
    """Load probe coordinates from an ASCII/CSV points file."""
    pts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            tokens = line.replace(",", " ").split()
            vals = [float(t) for t in tokens]
            if len(vals) == 2:
                vals.append(0.0)
            elif len(vals) < 2 or len(vals) > 3:
                raise ValueError(
                    f"Expected 2 or 3 coordinates per line in {path}, got '{line}'"
                )
            pts.append(vals)
    if not pts:
        raise ValueError(f"No probe coordinates found in {path}")
    return np.asarray(pts, dtype=np.float64)


# ---------------------------------------------------------------------------
# Field interpolation onto arbitrary probe coordinates
# ---------------------------------------------------------------------------
def _resolve_axis_indices(coord_names: Sequence[str]) -> list[int]:
    """Map frame coordinate column names to Cartesian indices (0=X, 1=Y, 2=Z)."""
    axis_map = {"X": 0, "Y": 1, "Z": 2}
    indices = []
    for cname in coord_names:
        key = cname.strip().upper()[:1]
        if key not in axis_map:
            raise ValueError(f"Unsupported coordinate axis name '{cname}'")
        indices.append(axis_map[key])
    return indices


def _try_structured_grid(
    coords: np.ndarray, tol: float = 1e-10
) -> tuple[list[np.ndarray], np.ndarray, list[int]] | None:
    """Detect whether *coords* forms a complete Cartesian tensor-product grid.

    Returns ``(active_axes, sort_order, active_dims)`` when *coords* lies on a
    rectilinear product grid, or ``None`` otherwise. Axes of length 1 are
    squeezed into *active_axes* so fake-2D 3D slices interpolate cleanly in
    plane without out-of-bounds failures along the invariant coordinate.
    """
    n_point, n_dim = coords.shape
    raw_axes = []
    # Deduplicate each coordinate column within relative spacing tolerance.
    for d in range(n_dim):
        vals = np.unique(coords[:, d])
        # Merge nearly identical floating-point coordinates (e.g. real4 noise).
        if len(vals) > 1:
            span = max(abs(vals[-1] - vals[0]), 1.0)
            merged = [vals[0]]
            for v in vals[1:]:
                if abs(v - merged[-1]) > tol * span:
                    merged.append(v)
            vals = np.asarray(merged, dtype=np.float64)
        raw_axes.append(vals)

    grid_sizes = [len(a) for a in raw_axes]
    if math.prod(grid_sizes) != n_point:
        return None

    # Sort points in lexicographic C-order (axis 0 slowest, axis n_dim-1 fastest)
    # after 1D searchsorted nearest-tick snapping (avoids N_pts x N_axis memory).
    snapped = []
    for d in range(n_dim):
        ax = raw_axes[d]
        col = coords[:, d]
        if len(ax) == 1:
            snapped.append(np.full(n_point, ax[0], dtype=np.float64))
        else:
            pos = np.searchsorted(ax, col, side="left")
            right = np.minimum(pos, len(ax) - 1)
            left = np.maximum(pos - 1, 0)
            choose_left = (col - ax[left]) <= (ax[right] - col)
            snapped.append(np.where(choose_left, ax[left], ax[right]))
    order = np.lexsort([snapped[d] for d in range(n_dim - 1, -1, -1)])
    for d in range(n_dim):
        expected = np.repeat(
            np.tile(raw_axes[d], math.prod(grid_sizes[:d])),
            math.prod(grid_sizes[d + 1 :]),
        )
        span = max(abs(raw_axes[d][-1] - raw_axes[d][0]), 1.0)
        if not np.allclose(coords[order, d], expected, atol=tol * span, rtol=0.0):
            return None

    active_dims = [d for d in range(n_dim) if len(raw_axes[d]) > 1]
    active_axes = [raw_axes[d] for d in active_dims]
    return active_axes, order, active_dims


def _prepare_sampling_plan(
    coords: np.ndarray,
    query: np.ndarray,
) -> dict:
    """Detect structured layout and precompute interpolation stencil for a mesh."""
    structured = _try_structured_grid(coords)
    if structured is None:
        return {"structured": False}

    active_axes, sort_order, active_dims = structured
    if not active_dims:
        return {
            "structured": True,
            "active_dims": active_dims,
        }

    shape = tuple(len(ax) for ax in active_axes)
    stencil = _build_rectilinear_stencil(active_axes, query[:, active_dims])
    return {
        "structured": True,
        "active_dims": active_dims,
        "sort_order": sort_order,
        "shape": shape,
        "stencil": stencil,
    }


def sample_frame(
    frame: dict,
    probe_coords: np.ndarray,
    var_names: Sequence[str] | None = None,
    _cache: dict | None = None,
) -> tuple[list[str], dict[str, np.ndarray]]:
    """Interpolate *var_names* from *frame* onto 3D *probe_coords*."""
    probe_coords = np.asarray(probe_coords, dtype=np.float64).reshape(-1, 3)

    # Resolve requested variable names case-insensitively while preserving order.
    file_vars_upper = {v.upper(): v for v in frame["var_names"]}
    if var_names is None:
        selected = list(frame["var_names"])
        out_names = list(frame["var_names"])
    else:
        selected = []
        out_names = []
        for req in var_names:
            key = req.strip().upper()
            if key not in file_vars_upper:
                raise KeyError(
                    f"Variable '{req}' not found in {frame['path']}. "
                    f"Available: {', '.join(frame['var_names'])}"
                )
            selected.append(file_vars_upper[key])
            out_names.append(req.strip())

    axis_indices = _resolve_axis_indices(frame["coord_names"])
    coords = frame["coords"]
    query = probe_coords[:, axis_indices]

    plan = None
    if _cache is not None:
        prev_axes = _cache.get("axis_indices")
        prev_coords = _cache.get("coords")
        if (
            prev_axes == axis_indices
            and prev_coords is not None
            and coords.shape == prev_coords.shape
            and np.array_equal(coords, prev_coords)
        ):
            plan = _cache["plan"]

    if plan is None:
        plan = _prepare_sampling_plan(coords, query)
        if _cache is not None:
            _cache["axis_indices"] = axis_indices
            _cache["coords"] = coords
            _cache["plan"] = plan

    results: dict[str, np.ndarray] = {}

    if plan["structured"]:
        if not plan["active_dims"]:
            for out_name, f_name in zip(out_names, selected):
                results[out_name] = np.full(
                    probe_coords.shape[0],
                    float(frame["fields"][f_name][0]),
                    dtype=np.float64,
                )
            return out_names, results

        sort_order = plan["sort_order"]
        shape = plan["shape"]
        stencil = plan["stencil"]
        for out_name, f_name in zip(out_names, selected):
            grid_vals = frame["fields"][f_name][sort_order].reshape(shape)
            results[out_name] = _apply_rectilinear_stencil(stencil, grid_vals)
        return out_names, results

    # Fallback for unstructured / AMR point clouds: build triangulation / k-NN
    # once across all selected variables.
    if not selected:
        return out_names, results

    val_matrix = np.column_stack([frame["fields"][f_name] for f_name in selected])
    try:
        from scipy.interpolate import LinearNDInterpolator
    except ImportError:
        LinearNDInterpolator = None

    if LinearNDInterpolator is not None:
        interp = LinearNDInterpolator(coords, val_matrix, fill_value=np.nan)
        sampled_mat = interp(query).astype(np.float64)
        if sampled_mat.ndim == 1:
            sampled_mat = sampled_mat[:, np.newaxis]
    else:
        sampled_mat = _interp_idw_nd(coords, val_matrix, query)

    for idx, out_name in enumerate(out_names):
        results[out_name] = sampled_mat[:, idx]
    return out_names, results


def sample_frames(
    paths: Sequence[str],
    probe_coords: np.ndarray,
    var_names: Sequence[str] | None = None,
) -> dict:
    """Sample a sorted sequence of frames at fixed 3D *probe_coords*."""
    if not paths:
        raise ValueError("No input frame paths provided")

    probe_coords = np.asarray(probe_coords, dtype=np.float64).reshape(-1, 3)
    frames_data = [read_frame(p, var_names=var_names) for p in paths]
    frames_data.sort(key=lambda d: (d["time"], d["nStep"], d["path"]))

    times = np.array([d["time"] for d in frames_data], dtype=np.float64)
    steps = np.array([d["nStep"] for d in frames_data], dtype=np.int64)

    resolved_names = None
    series: dict[str, list[np.ndarray]] = {}
    grid_cache: dict = {}
    for frame in frames_data:
        names, sampled = sample_frame(
            frame, probe_coords, var_names=var_names, _cache=grid_cache
        )
        if resolved_names is None:
            resolved_names = names
            series = {n: [] for n in resolved_names}
        elif names != resolved_names:
            raise ValueError(
                f"Inconsistent variable list across frames (saw {names}, "
                f"expected {resolved_names})"
            )
        for n in resolved_names:
            series[n].append(sampled[n])

    arrays = {n: np.vstack(vals) for n, vals in series.items()}
    return {
        "paths": [d["path"] for d in frames_data],
        "times": times,
        "steps": steps,
        "coords": probe_coords,
        "var_names": resolved_names or [],
        "data": arrays,
    }


# ---------------------------------------------------------------------------
# Writers (.csv and .npz)
# ---------------------------------------------------------------------------
def write_csv(path: str, result: dict) -> None:
    """Write sampled probe time series to a tidy CSV file."""
    var_names = result["var_names"]
    coords = result["coords"]
    n_frames = len(result["times"])
    n_probes = coords.shape[0]

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["frame", "step", "time", "probe_id", "x", "y", "z"] + var_names)
        for i_frame in range(n_frames):
            t_val = result["times"][i_frame]
            step_val = int(result["steps"][i_frame])
            for i_probe in range(n_probes):
                row = [
                    i_frame,
                    step_val,
                    f"{t_val:.10e}",
                    i_probe,
                    f"{coords[i_probe, 0]:.10e}",
                    f"{coords[i_probe, 1]:.10e}",
                    f"{coords[i_probe, 2]:.10e}",
                ]
                for vname in var_names:
                    row.append(f"{result['data'][vname][i_frame, i_probe]:.10e}")
                writer.writerow(row)


def write_npz(
    path: str, result: dict, plane_shape: tuple[int, int] | None = None
) -> None:
    """Write sampled probe time series to a compressed NumPy archive."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    payload = {
        "times": result["times"],
        "steps": result["steps"],
        "coords": result["coords"],
        "var_names": np.asarray(result["var_names"]),
    }
    if plane_shape is not None:
        payload["plane_shape"] = np.asarray(plane_shape, dtype=np.int64)
    for vname, arr in result["data"].items():
        payload[vname] = arr
    np.savez_compressed(path, **payload)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def _expand_inputs(inputs: Iterable[str]) -> list[str]:
    paths: list[str] = []
    for item in inputs:
        matched = sorted(glob.glob(item))
        if matched:
            paths.extend(matched)
        elif os.path.exists(item):
            paths.append(item)
        else:
            raise FileNotFoundError(f"Input file or pattern not found: {item}")
    return paths


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Extract virtual probe time series or spatial rakes from FLEKS outputs."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input .out file(s), glob pattern(s), or AMReX plotfile director(y/ies).",
    )
    parser.add_argument(
        "--point",
        action="append",
        default=[],
        metavar="X,Y[,Z]",
        help="Single probe coordinate (repeatable).",
    )
    parser.add_argument(
        "--points-file",
        help="ASCII/CSV file containing probe coordinates (2 or 3 values per row).",
    )
    parser.add_argument(
        "--line",
        nargs=3,
        metavar=("START", "END", "N"),
        help="1D linear probe rake: x0,y0,z0 x1,y1,z1 N.",
    )
    parser.add_argument(
        "--plane",
        nargs=5,
        metavar=("ORIGIN", "U", "V", "NU", "NV"),
        help="2D planar probe array: ox,oy,oz ux,uy,uz vx,vy,vz Nu Nv.",
    )
    parser.add_argument(
        "--vars",
        help="Comma-separated list of variables to extract (default: all).",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output path ending in .csv or .npz.",
    )

    args = parser.parse_args(argv)
    probe_blocks = []
    plane_shape = None

    for p_str in args.point:
        probe_blocks.append(_parse_vec3(p_str)[np.newaxis, :])
    if args.points_file:
        probe_blocks.append(load_points_file(args.points_file))
    if args.line:
        p0 = _parse_vec3(args.line[0])
        p1 = _parse_vec3(args.line[1])
        n_pts = int(args.line[2])
        probe_blocks.append(make_line_probes(p0, p1, n_pts))
    if args.plane:
        o = _parse_vec3(args.plane[0])
        u = _parse_vec3(args.plane[1])
        v = _parse_vec3(args.plane[2])
        nu = int(args.plane[3])
        nv = int(args.plane[4])
        probe_blocks.append(make_plane_probes(o, u, v, nu, nv))
        plane_shape = (nu, nv)

    if not probe_blocks:
        parser.error("Specify at least one of --point, --points-file, --line, or --plane.")

    probe_coords = np.vstack(probe_blocks)
    var_list = [v.strip() for v in args.vars.split(",") if v.strip()] if args.vars else None
    input_paths = _expand_inputs(args.inputs)
    result = sample_frames(input_paths, probe_coords, var_names=var_list)

    out_lower = args.output.lower()
    if out_lower.endswith(".csv"):
        write_csv(args.output, result)
    elif out_lower.endswith(".npz"):
        write_npz(args.output, result, plane_shape=plane_shape)
    else:
        parser.error("Output filename must end in .csv or .npz")

    return 0


if __name__ == "__main__":
    sys.exit(main())
