#!/usr/bin/env python3
"""Validator for the inner-body tests (tests/body/).

Four variants are discovered from this directory:

  - PARAM.in              -> "body"             (default: absorb + linetied)
  - PARAM.in.conducting   -> "body_conducting"  (absorb + conducting: PEC,
                             tangential E = 0 and radial B = 0)
  - PARAM.in.insulating   -> "body_insulating"  (absorb + insulating: no field
                             constraint, the magnetic field passes through)
  - PARAM.in.reflect      -> "body_reflect"     (reflect + linetied: specular
                             reflection on the sphere, nothing is absorbed)

All of them use the same setup: a uniform plasma streams in +x through a
periodic 2D box past an absorbing sphere of radius R_BODY at the origin.
"""
import logging
import math

from tests._shared import run_dir as _run_dir

logger = logging.getLogger(__name__)

set_run_dir = _run_dir.set_run_dir

R_BODY = 1.2          # #BODY radius in code units (see PARAM.in)
WAKE_MAX_FRAC = 0.5   # wake density must stay below 50% of the upstream value
ZERO_TOL = 1e-12      # "exactly zero" threshold inside the body
ETOT_GROWTH_MAX = 10.0  # the body must not inject energy
CONSTRAINT_TOL = 1e-6   # relative tolerance of the conducting constraint
Z_HALF = 0.05           # half of the z extent of the fake-2D decks (one cell)
PASS_THROUGH_MIN = 0.5  # insulating: |B| and |E| inside vs outside
# Density next to the staircase surface relative to the far field. The initial
# frame is uniform; the allowed band covers the particle noise (4 ppc).
SURFACE_RHO_MIN = 0.7
SURFACE_RHO_MAX = 1.35
EPART_KEEP_MIN = 0.8    # reflect: elastic reflection keeps the particle energy


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _columns(vidx, rows):
    """Return the columns needed by the field checks, or None if missing."""
    out = {}
    for name in ("BODY", "X", "Y", "Z", "EX", "EY", "EZ", "BX", "BY", "BZ",
                 "RHOS0"):
        out[name] = _run_dir.col(vidx, rows, name)
    return out


def _inside_indices(cols):
    body = cols["BODY"]
    return [i for i in range(len(body)) if body[i] > 0.5]


def _radial(cols, i):
    """Outward radial unit vector (and radius) of point i about the origin."""
    x, y, z = cols["X"][i], cols["Y"][i], cols["Z"][i]
    r = math.sqrt(x * x + y * y + z * z)
    if r <= 0.0:
        return None, 0.0
    return (x / r, y / r, z / r), r


def _radial_in_plane(cols, i):
    """In-plane (x, y) radial unit vector of point i about the origin."""
    x, y = cols["X"][i], cols["Y"][i]
    r = math.sqrt(x * x + y * y)
    if r <= 0.0:
        return None, 0.0
    return (x / r, y / r), r


def _vec(cols, prefix, i):
    return (cols[prefix + "X"][i], cols[prefix + "Y"][i], cols[prefix + "Z"][i])


def _grid_spacing(cols):
    """Grid spacing of the output frame (the .out files are node centred)."""
    xs = sorted(set(cols["X"]))
    if len(xs) < 2:
        return 0.0
    return min(b - a for a, b in zip(xs, xs[1:]))


def _median(values):
    values = sorted(values)
    n = len(values)
    if n == 0:
        return None
    mid = n // 2
    return values[mid] if n % 2 else 0.5 * (values[mid - 1] + values[mid])


def _norm(vec):
    return math.sqrt(sum(v * v for v in vec))


def _magnitude_inside_outside(cols, inside, prefix):
    """Max |field| inside the body and in an annulus just outside it."""
    inside_max = 0.0
    for i in inside:
        inside_max = max(inside_max, _norm(_vec(cols, prefix, i)))

    outside_max = 0.0
    for i in range(len(cols["BODY"])):
        n, r = _radial(cols, i)
        if n is None or cols["BODY"][i] > 0.5:
            continue
        if r < R_BODY + 0.1 or r > R_BODY + 0.6:
            continue
        outside_max = max(outside_max, _norm(_vec(cols, prefix, i)))

    return inside_max, outside_max


# ---------------------------------------------------------------------------
# Log checks
# ---------------------------------------------------------------------------
def validate_log(pic_diags=None, test_name=None):
    """Check the absorption tallies and that the run stays finite."""
    logger.debug("Validating %s (log)...", test_name)

    if not pic_diags or len(pic_diags) < 2:
        return True, "Passed (no pic log)"

    for entry in pic_diags:
        for key in ("Etot", "Ee", "Eb", "Epart"):
            if not math.isfinite(entry.get(key, 0.0)):
                return False, f"Non-finite {key} (NaN/Inf)"

    first, last = pic_diags[0], pic_diags[-1]
    e0 = first.get("Etot", 0.0)
    e1 = last.get("Etot", 0.0)
    logger.debug("    Etot: %.4e -> %.4e", e0, e1)
    if e0 > 0 and e1 > ETOT_GROWTH_MAX * e0:
        return False, (f"Energy blew up: Etot {e0:.3e} -> {e1:.3e} "
                       f"(>{ETOT_GROWTH_MAX:g}x)")

    absorbed = [entry.get("nBodyAbsorb") for entry in pic_diags]
    if any(a is None for a in absorbed):
        return False, ("nBodyAbsorb column missing from log_pic -- the body "
                       "tallies are not written")

    logger.debug("    nBodyAbsorb: %.4e -> %.4e", absorbed[0], absorbed[-1])
    for prev, cur in zip(absorbed[:-1], absorbed[1:]):
        if cur < prev - 1e-9:
            return False, "nBodyAbsorb decreased (tallies are cumulative)"

    if test_name == "body_reflect":
        # Reflection keeps every particle: nothing may be absorbed.
        if absorbed[-1] > 0:
            return False, (f"{absorbed[-1]:.0f} particles were absorbed by a "
                           "reflecting body")
        ep0 = first.get("Epart", 0.0)
        ep1 = last.get("Epart", 0.0)
        if ep0 > 0 and ep1 < EPART_KEEP_MIN * ep0:
            return False, (f"Reflecting body lost particle energy: Epart "
                           f"{ep0:.3e} -> {ep1:.3e}")
        return True, "Passed (nothing absorbed, particle energy kept)"

    if absorbed[-1] <= 0:
        return False, ("No particle was absorbed by the body (nBodyAbsorb = 0) "
                       "-- the inner boundary is inactive")

    return True, (f"Passed (body absorbed {absorbed[-1]:.0f} particles, "
                  "tallies cumulative)")


# ---------------------------------------------------------------------------
# Plot checks
# ---------------------------------------------------------------------------
def validate_plot(test_name):
    """Check the body interior and, for the default variant, the wake."""
    logger.debug("Validating %s (plot)...", test_name)

    vidx, rows = _run_dir.load_last_out()
    if vidx is None or not rows:
        return True, "No .out frames (skipped)"

    cols = _columns(vidx, rows)
    if cols["BODY"] is None:
        return False, "BODY column missing from the output"

    # A collapsed 2D frame may not carry every component (e.g. Z): treat the
    # missing ones as zero instead of failing on the column lookup.
    for name, values in cols.items():
        if name == "BODY":
            continue
        if values is None:
            cols[name] = [0.0] * len(rows)
            continue
        if not all(math.isfinite(v) for v in values):
            return False, f"Non-finite {name} in the final plot frame"

    inside = _inside_indices(cols)
    if not inside:
        return False, "No point is marked as inside the body (mask is empty)"

    # The initial frame is uniform, so it shows whether the moments next to the
    # staircase surface are diluted by the empty body cells.
    first_vidx, first_rows = _run_dir.load_first_out()
    first_cols = None
    if first_vidx is not None and first_rows:
        first_cols = _columns(first_vidx, first_rows)
        ok, msg = _check_surface_density(first_cols)
        if not ok:
            return False, msg
        logger.debug("    initial frame: %s", msg)

    if test_name == "body_conducting":
        return _check_conducting(cols, inside, first_cols)
    if test_name == "body_insulating":
        return _check_insulating(cols, inside)
    if test_name == "body_reflect":
        return _check_reflect(cols, inside)

    return _check_linetied(cols, inside, rows)


def _check_surface_density(cols):
    """The density next to the body must be the density of the plasma there.

    A node-centred CIC moment averages over the cells around the node and the
    cells inside the body are empty, so without a correction the nodes on the
    staircase surface report only a half or three quarters of the plasma
    density. The moments are rescaled by the fraction of the surrounding cells
    that are outside the body, so the first frame, which is uniform, has to
    come out uniform right up to the surface.
    """
    rho = cols.get("RHOS0")
    if rho is None:
        return True, "No rhoS0 column (surface density not checked)"

    dx = _grid_spacing(cols)
    if dx <= 0:
        return True, "Cannot determine the grid spacing (not checked)"

    far = _median([rho[i] for i in range(len(rho))
                   if _radial_in_plane(cols, i)[1] > R_BODY + 0.6])
    if not far:
        return True, "No far-field sample (surface density not checked)"

    ratios = [rho[i] / far for i in range(len(rho))
              if cols["BODY"][i] < 0.5 and
              R_BODY - 1.5 * dx <= _radial_in_plane(cols, i)[1] <= R_BODY + 1.5 * dx]
    if not ratios:
        return True, "No node next to the body (surface density not checked)"

    logger.debug("    rho/rho_far next to the body: %.3f .. %.3f (n = %d)",
                 min(ratios), max(ratios), len(ratios))
    if min(ratios) < SURFACE_RHO_MIN or max(ratios) > SURFACE_RHO_MAX:
        return False, (f"The density next to the body is not the plasma density "
                       f"(rho/rho_far in [{min(ratios):.3f}, {max(ratios):.3f}], "
                       f"expected [{SURFACE_RHO_MIN}, {SURFACE_RHO_MAX}])")

    return True, (f"surface density uniform (rho/rho_far in "
                  f"[{min(ratios):.3f}, {max(ratios):.3f}])")


def _check_linetied(cols, inside, rows):
    """Default: the plasma moments and the electric field vanish inside."""
    for name in ("RHOS0", "RHOS1", "EX", "EY", "EZ"):
        values = cols.get(name)
        if values is None:
            continue
        peak = max(abs(values[i]) for i in inside)
        logger.debug("    max |%s| inside the body = %.3e", name, peak)
        if peak > ZERO_TOL:
            return False, (f"{name} is not zero inside the body "
                           f"(max |{name}| = {peak:.3e})")

    rho = cols["RHOS0"]
    x, y = cols["X"], cols["Y"]

    def mean_rho(x_lo, x_hi):
        sel = [rho[i] for i in range(len(rows))
               if x_lo <= x[i] <= x_hi and abs(y[i]) <= 0.5 * R_BODY]
        return (sum(sel) / len(sel)) if sel else None

    upstream = mean_rho(-2.4, -1.3)
    wake = mean_rho(R_BODY + 0.1, 2.4)

    if upstream is None or wake is None:
        return True, "Passed (wake sample regions empty)"
    if upstream <= 0:
        return False, "Upstream density is zero (plasma not initialized)"

    logger.debug("    rho: upstream %.4e, wake %.4e", upstream, wake)
    if wake >= WAKE_MAX_FRAC * upstream:
        return False, (f"No wake behind the body (wake/upstream = "
                       f"{wake / upstream:.3f} >= {WAKE_MAX_FRAC}) -- particles "
                       "are not absorbed")

    return True, (f"Passed (body interior empty, wake/upstream = "
                  f"{wake / upstream:.3f})")


def _check_conducting(cols, inside, first_cols=None):
    """conducting: E is purely radial (E_t = 0) and B is purely tangential.

    The conditions are surface conditions: they hold on the one-cell-thick
    surface layer of the body, while the interior is a cavity with E = 0 and B
    frozen at its initial value. The two layers are told apart by the radius,
    with a band of ambiguous nodes in between that is not checked (the
    staircase surface sits at R_BODY - 1.5 dx to R_BODY - 2.5 dx).

    The check uses the in-plane (x, y) components: the run is fake 2D (one
    cell in z), the plot output carries no z coordinate, and the radial
    direction of a body node has a z component of order dz/R that cannot be
    reconstructed from the output.
    """
    n_pt = len(cols["BODY"])
    scale_e = max(math.hypot(cols["EX"][i], cols["EY"][i]) for i in range(n_pt))
    scale_b = max(math.hypot(cols["BX"][i], cols["BY"][i]) for i in range(n_pt))
    if scale_e <= 0 or scale_b <= 0:
        return False, "The ambient E or B field is zero (test is vacuous)"

    dx = _grid_spacing(cols)
    r_surf = R_BODY - 1.5 * dx   # solidly in the surface layer
    r_int = R_BODY - 2.5 * dx    # solidly in the frozen interior

    max_et = 0.0
    max_br = 0.0
    max_e_int = 0.0
    max_db_int = 0.0
    n_surf = n_int = 0
    for i in inside:
        n2, r2 = _radial_in_plane(cols, i)
        if n2 is None:
            continue

        if r2 >= r_surf:
            n_surf += 1

            # E x n = 0 gives Ex*y - Ey*x = 0, which does not involve z, so
            # the in-plane tangential field must vanish exactly.
            ex, ey = cols["EX"][i], cols["EY"][i]
            er = ex * n2[0] + ey * n2[1]
            max_et = max(max_et, math.hypot(ex - er * n2[0], ey - er * n2[1]))

            # B . n = 0 reads Bx*x + By*y + Bz*z = 0. The output plane carries
            # no z and a body node sits at z = +-dz/2, so the in-plane part may
            # keep a residual of |Bz| * dz/2 / r.
            bx, by = cols["BX"][i], cols["BY"][i]
            br = abs(bx * n2[0] + by * n2[1])
            allowed = abs(cols["BZ"][i]) * Z_HALF / r2
            max_br = max(max_br, max(0.0, br - allowed))
        elif r2 <= r_int:
            n_int += 1

            # The interior is shielded: no electric field and the magnetic
            # field keeps its initial value.
            max_e_int = max(max_e_int, math.hypot(cols["EX"][i], cols["EY"][i],
                                                  cols["EZ"][i]))
            if first_cols is not None:
                db = math.sqrt(sum((cols["B" + d][i] - first_cols["B" + d][i]) ** 2
                                   for d in ("X", "Y", "Z")))
                max_db_int = max(max_db_int, db)

    logger.debug("    surface: %d nodes, max |E_t| = %.3e (|E| = %.3e), "
                 "max |B_r| beyond the fake-2D residual = %.3e (|B| = %.3e)",
                 n_surf, max_et, scale_e, max_br, scale_b)
    logger.debug("    interior: %d nodes, max |E| = %.3e, max |B - B(t=0)| = %.3e",
                 n_int, max_e_int, max_db_int)

    if n_surf == 0 or n_int == 0:
        return False, ("The conducting body is too small to separate the "
                       "surface layer from the interior")

    if max_et > CONSTRAINT_TOL * scale_e:
        return False, (f"Tangential E does not vanish on the conducting surface "
                       f"(max |E_t| = {max_et:.3e} > "
                       f"{CONSTRAINT_TOL} * |E| = {CONSTRAINT_TOL * scale_e:.3e})")
    if max_br > CONSTRAINT_TOL * scale_b:
        return False, (f"Radial B does not vanish on the conducting surface "
                       f"(max |B_r| = {max_br:.3e} > "
                       f"{CONSTRAINT_TOL} * |B| = {CONSTRAINT_TOL * scale_b:.3e})")

    if max_e_int > CONSTRAINT_TOL * scale_e:
        return False, (f"The electric field is not zero inside the conducting "
                       f"body (max |E| = {max_e_int:.3e})")
    if first_cols is not None and max_db_int > CONSTRAINT_TOL * scale_b:
        return False, (f"The magnetic field is not frozen inside the conducting "
                       f"body (max |B - B(t=0)| = {max_db_int:.3e})")

    return True, (f"Passed (E radial and B tangential on the {n_surf} surface "
                  f"nodes; E = 0 and B frozen on the {n_int} interior nodes)")


def _check_insulating(cols, inside):
    """insulating: no field constraint, so |B| and |E| pass through the body."""
    e_in, e_out = _magnitude_inside_outside(cols, inside, "E")
    b_in, b_out = _magnitude_inside_outside(cols, inside, "B")
    logger.debug("    |E|: inside %.3e, outside %.3e; |B|: inside %.3e, "
                 "outside %.3e", e_in, e_out, b_in, b_out)

    if b_out <= 0:
        return False, "No magnetic field outside the body (test is vacuous)"

    ratio = b_in / b_out
    if not (PASS_THROUGH_MIN <= ratio <= 1.0 / PASS_THROUGH_MIN):
        return False, (f"The magnetic field does not pass through the "
                       f"insulating body (|B| inside/outside = {ratio:.3f})")

    if e_in <= 0:
        return False, ("The electric field inside the insulating body is zero "
                       "-- it should not be constrained")

    return True, (f"Passed (fields pass through the insulating body, "
                  f"|B| inside/outside = {ratio:.3f})")


def _check_reflect(cols, inside):
    """reflect: nothing is inside the body, but the fields are untouched."""
    rho = cols["RHOS0"]
    if rho is None:
        return True, "No rhoS0 column (skipped)"
    peak = max(abs(rho[i]) for i in inside)
    logger.debug("    max |rhoS0| inside the body = %.3e", peak)
    if peak > ZERO_TOL:
        return False, (f"Particles inside a reflecting body "
                       f"(max |rhoS0| = {peak:.3e})")
    return True, "Passed (no particle inside the reflecting body)"
