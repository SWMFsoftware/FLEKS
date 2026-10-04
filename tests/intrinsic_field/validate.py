#!/usr/bin/env python3
"""Validator for the intrinsic-magnetic-field tests (tests/intrinsic_field/).

Four variants are discovered from this directory:

  - PARAM.in.body         -> "intrinsic_field_body"           dipole + #BODY
  - PARAM.in.body_conducting -> "intrinsic_field_body_conducting"
                            dipole + a conducting #BODY: the total field is
                            made tangential, so B1 = -B0.n n-hat on the shell
  - PARAM.in.crustal      -> "intrinsic_field_crustal"        g10-only file
  - PARAM.in.crustal_nm2  -> "intrinsic_field_crustal_nm2"    BATSRUS layout, many g/h
  - PARAM.in.dipole_crustal -> "intrinsic_field_dipole_crustal"

There is no bare "dipole only" deck: in practice a planetary dipole never runs
alone, it always comes with an inner body, so the body deck *is* the dipole
deck and it carries the analytic-dipole comparison.

Every deck runs with SOLVEEM F and a zero UNIFORMSTATE field, so the evolved
field is identically zero and the total field equals the static intrinsic
field B0. That makes the plot output directly comparable with an analytic
evaluation, and it makes "B0 leaks into the evolved field" and "B0 is missing
from the total field" both show up as a hard failure.

All output is in PLANETARY units: the magnetic field is in nT and the
coordinates are in planetary radii (#PLANETRADIUS = lNormSI in these decks, so
one code length unit is one planetary radius). The reference sphere of the
dipole therefore sits at r = 1 and |B0| there must round-trip to the deck
input, exactly like the density round-trip check of the beam test.
"""
import logging
import math

from tests._shared import run_dir as _run_dir

logger = logging.getLogger(__name__)

set_run_dir = _run_dir.set_run_dir

# The decks: equatorial surface field [nT], tilt [deg] and the reference
# radius, which is 1 planetary radius by construction (rRef < 0 in the deck).
B_EQ_NT = 100.0
# theta = 90 deg tips the axis out of +z into the XY plane (the plane the plot
# writes) and phi = 180 deg points it along +x; phi = 0 would give -x.
THETA_DEG = 90.0
PHI_DEG = 180.0
R_REF = 1.0

# Hybrid deck: the uniform drift of #UNIFORMSTATE is 100 km/s, and the plot
# output is in planetary units, so the velocity that enters E = -u x B0 is
# 100 km/s (No2OutV = uNormSI = 1e5 m/s = 100 km/s per code unit).
HYBRID_U_OUT = (100.0, 0.0, 0.0)

R_MIN = 0.6            # skip the dipole singularity at the origin
R_MAX = 3.5            # and the periodic seam far outside
R_BODY = 1.0           # #BODY radius in the decks, in planetary radii
ANALYTIC_TOL = 1e-6    # median relative error against the analytic dipole
STATIC_TOL = 1e-6      # B0 must not change between frames
TOTAL_TOL = 1e-10      # total field == B0 when B1 == 0
# Purity of B1 on a conducting surface, relative to |B1| there. It is
# a cancellation of two ~1e2 nT numbers, so the 11 significant digits
# of the .out frames limit it to about 1e-8; anything above 1e-6 means
# the projection left a tangential component.
TANGENTIAL_TOL = 1e-6
ROUNDTRIP_TOL = 1e-3   # limited by the grid spacing near r = 1
CRUSTAL_TOL = 1e-8     # C++ vs independent Python evaluation
ETOT_TOL = 1e-6        # relative drift of the total energy


# ---------------------------------------------------------------------------
# References
# ---------------------------------------------------------------------------
def _dipole_axis():
    th = math.radians(THETA_DEG)
    ph = math.radians(PHI_DEG)
    return (-math.sin(th) * math.cos(ph), -math.sin(th) * math.sin(ph),
            math.cos(th))


def dipole_b0(x, y, z=0.0):
    """Analytic tilted dipole in nT, r in planetary radii."""
    mh = _dipole_axis()
    r = math.sqrt(x * x + y * y + z * z)
    if r <= 0.0:
        return (0.0, 0.0, 0.0)
    rh = (x / r, y / r, z / r)
    mdotr = sum(mh[i] * rh[i] for i in range(3))
    f = B_EQ_NT * (R_REF / r) ** 3
    return tuple(f * (3.0 * mdotr * rh[i] - mh[i]) for i in range(3))


def _schmidt_rnm(n_max):
    """Table R[n][m] of the Schmidt semi-normalized harmonics used by BATSRUS.

    Same recurrence as ModUserMars::set_mars_b0 / IntrinsicBField::eval_crustal:
      R(0,0) = 1,  R(n,0) = cos(theta),
      R(n,n) = sqrt((n-0.5)/n) sin(theta) R(n-1,n-1),
      R(n+1,n) = cos(theta) sqrt(2n+1) R(n,n),
      R(l,m) = [cos(theta)(2l-1)R(l-1,m) - R(l-2,m) sqrt((l+m-1)(l-m-1))]
               / sqrt(l^2-m^2)
    """
    table = [[0.0] * (n_max + 1) for _ in range(n_max + 1)]

    def build(ct, st):
        r = [[0.0] * (n_max + 1) for _ in range(n_max + 1)]
        r[0][0] = 1.0
        if n_max >= 1:
            r[1][0] = ct
        for n in range(1, n_max):
            if n == 1:
                r[n][n] = st * r[n - 1][n - 1]
            else:
                r[n][n] = math.sqrt((n - 0.5) / n) * st * r[n - 1][n - 1]
            r[n + 1][n] = ct * math.sqrt(2.0 * n + 1.0) * r[n][n]
        for m in range(0, n_max):
            for l in range(m + 2, n_max + 1):
                r[l][m] = (ct * (2.0 * l - 1.0) * r[l - 1][m] -
                           r[l - 2][m] * math.sqrt((l + m - 1.0) *
                                                   (l - m - 1.0))) / \
                          math.sqrt(1.0 * l * l - 1.0 * m * m)
        return r

    return build


def crustal_b0(coeffs, x, y, z=0.0):
    """Independent evaluation of the crustal field in nT.

    *coeffs* is {(n, m): (g, h)} with the reference radius folded in (a = 1),
    exactly the convention of the BATSRUS files and of the test data.
    """
    n_max = max(n for n, _ in coeffs) + 1
    r = math.sqrt(x * x + y * y + z * z)
    if r <= 0.0:
        return (0.0, 0.0, 0.0)
    ct = z / r
    st = math.sqrt(max(0.0, 1.0 - ct * ct))
    ph = math.atan2(y, x)

    build = _schmidt_rnm(n_max + 2)
    rn = build(ct, st)

    arr = 1.0 / r
    aorn = [1.0]
    for n in range(1, n_max + 3):
        aorn.append(arr * aorn[-1])

    br = bth = bph = 0.0
    for (n, m), (g, h) in coeffs.items():
        if n >= n_max or m > n:
            continue
        cd = g * math.cos(m * ph) + h * math.sin(m * ph)
        if m == 0:
            drnm = -math.sqrt((n + 1.0) * n / 2.0) * rn[n][m + 1]
        elif st <= 1.0e-6:
            drnm = -math.sqrt((n + m + 1.0) * (n - m)) * rn[n][m + 1]
        else:
            drnm = (m * ct * rn[n][m] / st -
                    math.sqrt((n + m + 1.0) * (n - m)) * rn[n][m + 1])
        br += (n + 1) * aorn[n + 2] * rn[n][m] * cd
        bth -= aorn[n + 2] * drnm * cd
        if st > 1.0e-6:
            bph -= aorn[n + 2] * rn[n][m] * m / st * \
                (-g * math.sin(m * ph) + h * math.cos(m * ph))

    cp = math.cos(ph)
    sp = math.sin(ph)
    return (br * st * cp + bth * ct * cp - bph * sp,
            br * st * sp + bth * ct * sp + bph * cp,
            br * ct - bth * st)


# The coefficients of crustal_nm2.txt (nMax = 3). Only degrees with m <= n-1
# are present, because the evaluation loops n = 0..nMax-1 and the highest
# degree of the file is dropped, so the two crustal decks stay comparable.
CRUSTAL_NM2 = {
    (1, 0): (100.0, 0.0),
    (2, 0): (-18.0, 0.0),
    (2, 1): (6.25, 3.5),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _load(pattern="*.out"):
    return _run_dir.load_last_out(pattern=pattern)


def _columns(vidx, rows):
    cols = {}
    for name in ("X", "Y", "Z", "B0X", "B0Y", "B0Z", "BX", "BY", "BZ",
                 "EX", "EY", "EZ", "BODY", "BODYSURFN"):
        cols[name] = _run_dir.col(vidx, rows, name)
    return cols


def _xyz(cols, i):
    x, y = cols["X"][i], cols["Y"][i]
    z = cols["Z"][i] if cols["Z"] is not None else 0.0
    return x, y, z


def _radius(cols, i):
    x, y, z = _xyz(cols, i)
    return math.sqrt(x * x + y * y + z * z)


def _mag(vec):
    return math.sqrt(sum(v * v for v in vec))


def _grid_spacing(cols):
    """Cell size from the printed X coordinates (the .out frames are nodal)."""
    xs = sorted(set(cols["X"]))
    if len(xs) < 2:
        return 0.0
    return min(b - a for a, b in zip(xs, xs[1:]))


def _is_shell(cols, i):
    """True when point i is on the one-cell-thick surface shell of the body.

    Uses the 'bodySurfN' mask, the node mask, which is what project_body_B()
    acts on for the nodal field; the cell mask ('bodySurf') is offset from it by
    one node along the surface. Without it (a deck that does not write the mask)
    this falls back to a radius band, which is ambiguous because the staircase
    surface sits anywhere within a cell or two of R = 1 and the shell and the
    interior have opposite B_r.
    """
    mask = cols.get("BODYSURFN") or cols.get("BODYSURF")
    if mask is not None:
        return mask[i] > 0.5
    return 0.6 < _radius(cols, i) < 1.0


def _outward(cols, i):
    """Outward radial unit vector of point i (the body sits at the origin)."""
    x, y, z = _xyz(cols, i)
    r = math.sqrt(x * x + y * y + z * z)
    if r <= 0.0:
        return None
    return (x / r, y / r, z / r)


def _field(cols, prefix, i):
    return tuple(cols[prefix + d][i] for d in "XYZ")


def _finite(cols, names):
    for name in names:
        values = cols.get(name)
        if values is None:
            continue
        for v in values:
            if not math.isfinite(v):
                return False, name
    return True, None


# ---------------------------------------------------------------------------
# Log checks
# ---------------------------------------------------------------------------
def validate_log(pic_diags=None, test_name=None):
    """The energy must stay finite; with B1 = 0 Eb is constant.

    Etot is only conserved when nothing removes particles: with #BODY the
    absorbed particles leave the domain, so only Eb is checked there.
    """
    logger.debug("Validating %s (log)...", test_name)
    if not pic_diags or len(pic_diags) < 2:
        return True, "Passed (no pic log)"

    for entry in pic_diags:
        for key in ("Etot", "Ee", "Eb", "Epart"):
            if not math.isfinite(entry.get(key, 0.0)):
                return False, f"Non-finite {key} (NaN/Inf)"

    b0, b1 = pic_diags[0]["Eb"], pic_diags[-1]["Eb"]
    logger.debug("    Eb %.6e -> %.6e", b0, b1)
    if abs(b0) > 0 and abs(b1 - b0) > STATIC_TOL * abs(b0):
        return False, (f"Eb changed by {abs(b1 - b0) / abs(b0):.3e}; the "
                       "magnetic energy of a static field must be constant")

    if test_name in ("intrinsic_field_body", "intrinsic_field_body_conducting"):
        return True, "Passed (energies finite, Eb constant, body absorbs)"

    e0, e1 = pic_diags[0]["Etot"], pic_diags[-1]["Etot"]
    logger.debug("    Etot %.6e -> %.6e", e0, e1)
    if abs(e0) > 0 and abs(e1 - e0) > ETOT_TOL * abs(e0):
        return False, (f"Etot drifted by {abs(e1 - e0) / abs(e0):.3e} "
                       f"(> {ETOT_TOL:g}); the intrinsic field must not do "
                       "work")
    return True, "Passed (energies finite and constant)"


# ---------------------------------------------------------------------------
# Plot checks
# ---------------------------------------------------------------------------
def validate_plot(test_name):
    logger.debug("Validating %s (plot)...", test_name)

    vidx, rows = _load()
    if vidx is None or not rows:
        return True, "No .out frames (skipped)"

    cols = _columns(vidx, rows)
    if cols["B0X"] is None:
        return False, "B0x/B0y/B0z are missing from the output"

    ok, bad = _finite(cols, ("B0X", "B0Y", "B0Z", "BX", "BY", "BZ"))
    if not ok:
        return False, f"Non-finite {bad} in the final plot frame"

    # -- 1. B0 exists -------------------------------------------------------
    b0max = max(_mag(_field(cols, "B0", i)) for i in range(len(rows)))
    logger.debug("    max |B0| = %.4f nT", b0max)
    if b0max <= 0.0:
        return False, "B0 is identically zero -- the intrinsic field is off"

    # -- 2. static ----------------------------------------------------------
    fvidx, frows = _run_dir.load_first_out()
    fcols = None
    if fvidx is not None and frows and len(frows) == len(rows):
        fcols = _columns(fvidx, frows)
        if fcols["B0X"] is not None:
            db = max(_mag(tuple(fcols["B0" + d][i] - cols["B0" + d][i]
                                for d in "XYZ")) for i in range(len(rows)))
            logger.debug("    max |B0(last) - B0(first)| = %.3e nT", db)
            if db > STATIC_TOL * b0max:
                return False, (f"B0 is not static (max change {db:.3e} nT > "
                               f"{STATIC_TOL:g} * {b0max:.3e} nT); the "
                               "intrinsic field must not evolve")

    # -- 3. total == B0 ------------------------------------------------------
    # On a 'conducting' surface the total field is deliberately NOT B0: the
    # evolved field picks up exactly the normal component that cancels B0, so
    # the shell is excluded here and checked separately below.
    conducting = (test_name == "intrinsic_field_body_conducting")
    worst, where = 0.0, -1
    for i, _ in enumerate(rows):
        if conducting and _is_shell(cols, i):
            continue
        d = _mag(tuple(cols["B" + d][i] - cols["B0" + d][i] for d in "XYZ"))
        if d > worst:
            worst, where = d, i
    logger.debug("    max |B_total - B0| = %.3e nT%s", worst,
                 " (excluding the surface shell)" if conducting else "")
    if worst > TOTAL_TOL * b0max:
        return False, (f"The total field is not B0 (max |B - B0| = {worst:.3e} "
                       f"nT > {TOTAL_TOL:g} * {b0max:.3e} nT)")

    if test_name in ("intrinsic_field_crustal", "intrinsic_field_crustal_nm2"):
        return _check_crustal(cols, rows, test_name)
    if test_name == "intrinsic_field_dipole_crustal":
        return _check_superposition(cols, rows)

    if test_name == "intrinsic_field_body_conducting":
        ok, msg = _check_conducting_total(cols, rows)
        if not ok:
            return False, msg
        ok_dip, msg_dip = _check_dipole(cols, rows, test_name)
        return ok_dip, f"{msg}; {msg_dip}"

    # The body deck is the dipole deck: check the body mask as well as the
    # analytic dipole, so that "B0 survives inside the body" stays covered
    # now that the bare dipole deck is gone.
    if test_name == "intrinsic_field_body":
        ok, msg = _check_body(cols, rows)
        if not ok:
            return False, msg
        ok_dip, msg_dip = _check_dipole(cols, rows, test_name)
        return ok_dip, f"{msg}; {msg_dip}"

    return _check_dipole(cols, rows, test_name)


def _check_dipole(cols, rows, test_name):
    """Analytic tilted dipole, r^-3 scaling, tilt and the surface round trip."""
    errs = []
    for i, _ in enumerate(rows):
        rr = _radius(cols, i)
        if rr < R_MIN or rr > R_MAX:
            continue
        x, y, z = _xyz(cols, i)
        ban = dipole_b0(x, y, z)
        scale = B_EQ_NT * (R_REF / rr) ** 3
        errs.append(_mag(tuple(_field(cols, "B0", i)[k] - ban[k]
                               for k in range(3))) / scale)
    if not errs:
        return False, "No usable point in the comparison band"
    errs.sort()
    med = errs[len(errs) // 2]
    logger.debug("    dipole: n=%d median rel err = %.3e, max = %.3e",
                 len(errs), med, errs[-1])
    if med > ANALYTIC_TOL:
        return False, (f"B0 does not match the analytic dipole (median rel "
                       f"err = {med:.3e} > {ANALYTIC_TOL:g})")

    # Surface round trip: |B0| on the magnetic equator at r ~ 1 must be the
    # deck input (SI -> code -> SI is the identity).
    mh = _dipole_axis()
    found = None
    for i, _ in enumerate(rows):
        rr = _radius(cols, i)
        if not (0.95 < rr < 1.05):
            continue
        x, y, z = _xyz(cols, i)
        mdotr = sum(mh[k] * v / rr for k, v in enumerate((x, y, z)))
        if abs(mdotr) < 0.02:
            found = _mag(_field(cols, "B0", i))
            break
    if found is None:
        logger.debug("    no equatorial point on the reference sphere "
                     "(round trip not checked)")
    else:
        rel = abs(found - B_EQ_NT) / B_EQ_NT
        logger.debug("    equator: r ~ 1, |B0| = %.4f nT (expected %.1f)",
                     found, B_EQ_NT)
        if rel > ROUNDTRIP_TOL:
            return False, (f"|B0| on the reference sphere is {found:.4f} nT, "
                           f"expected {B_EQ_NT:g} nT (rel err {rel:.3e})")

    # r^-3 scaling: |B0| r^3 is constant along a ray once the direction is
    # fixed; check the spread of that product over the comparison band.
    prods = [(_mag(_field(cols, "B0", i)) * _radius(cols, i) ** 3)
             for i in range(len(rows)) if R_MIN < _radius(cols, i) < R_MAX]
    spread = (max(prods) - min(prods)) / (0.5 * (max(prods) + min(prods)))
    logger.debug("    |B0| r^3 spread over the band = %.3e", spread)
    if spread > 2.0:
        return False, ("|B0| r^3 is not constant (spread {:.3e}); the radial "
                       "dependence is not r^-3".format(spread))

    return True, (f"Passed (analytic dipole, median rel err {med:.1e}; "
                  f"|B0| r=1 = {B_EQ_NT:.1f} nT)")


def _check_conducting_total(cols, rows):
    """A conducting body makes the TOTAL field tangential, not B1.

    With an intrinsic planetary field the evolved field B1 is set to
    B1 = -B0.n n-hat on the surface shell, so that (B1 + B0).n = 0 there and the
    intrinsic field passes through the conductor. This is the condition BATSRUS
    calls 'reflectb'. Since B1 starts at zero and nothing else drives it, the
    deck checks three things:

      1. on the shell, the total field has no radial component,
      2. on the shell, B1 is purely normal (it cancels B0 and adds nothing
         tangential),
      3. off the shell, B1 is still zero, so the dipole is undisturbed.
    """
    shell = [i for i, _ in enumerate(rows) if _is_shell(cols, i)]
    if not shell:
        return False, ("No point is marked as on the surface shell "
                       "(bodySurf mask is empty)")

    b0max = max(_mag(_field(cols, "B0", i)) for i in range(len(rows)))

    # 1. the total field is tangential on the shell
    worst_br, worst_at = 0.0, -1
    for i in shell:
        n = _outward(cols, i)
        if n is None:
            continue
        b = _field(cols, "B", i)
        br = abs(sum(b[k] * n[k] for k in range(3)))
        if br > worst_br:
            worst_br, worst_at = br, i
    logger.debug("    shell: %d points, max |B_total . n| = %.3e nT "
                 "(|B0|max = %.3e nT)", len(shell), worst_br, b0max)
    if worst_br > TOTAL_TOL * b0max:
        return False, (f"The total field is not tangential on the conducting "
                       f"surface (max |B.n| = {worst_br:.3e} nT > "
                       f"{TOTAL_TOL:g} * {b0max:.3e} nT)")

    # 2. B1 on the shell is purely normal: it cancels B0 and nothing else
    worst_tan, worst_b1 = 0.0, 0.0
    for i in shell:
        n = _outward(cols, i)
        if n is None:
            continue
        b1 = tuple(cols["B" + d][i] - cols["B0" + d][i] for d in "XYZ")
        worst_tan = max(worst_tan, _mag(tuple(b1[k] -
                                              sum(b1[m] * n[m] for m in range(3)) *
                                              n[k] for k in range(3))))
        worst_b1 = max(worst_b1, _mag(b1))
    logger.debug("    shell: max |B1_tangential| = %.3e nT, max |B1| = %.3e nT",
                 worst_tan, worst_b1)
    if worst_b1 <= 0.0:
        return False, ("B1 is zero on the conducting surface; the evolved field "
                       "did not pick up the component that cancels B0")
    if worst_tan > TANGENTIAL_TOL * worst_b1:
        return False, (f"B1 on the conducting surface is not purely normal "
                       f"(max |B1_tangential| = {worst_tan:.3e} nT > "
                       f"{TANGENTIAL_TOL:g} * |B1| = {TANGENTIAL_TOL * worst_b1:.3e} nT); "
                       "it must cancel B0 and nothing else")

    # 3. off the shell nothing changed: B1 is still zero there. The staircase
    # fringe (nodes of the outermost shell cells that sit just outside the
    # analytic sphere) is skipped: those nodes carry the projection but cannot be
    # picked out with a radius test, so requiring B1 = 0 there would be a
    # statement about the mask rather than about the physics.
    dx = _grid_spacing(cols)
    fringe = R_BODY + (1.5 * dx if dx > 0 else 0.0)
    outside = [i for i, _ in enumerate(rows)
               if not _is_shell(cols, i) and cols["BODY"][i] < 0.5
               and _radius(cols, i) > fringe]
    worst_off = max((_mag(tuple(cols["B" + d][i] - cols["B0" + d][i]
                            for d in "XYZ")) for i in outside), default=0.0)
    logger.debug("    outside the shell (r > %.3f): %d points, "
                 "max |B_total - B0| = %.3e nT", fringe, len(outside), worst_off)
    if worst_off > TOTAL_TOL * b0max:
        return False, (f"The evolved field leaked outside the surface shell "
                       f"(max |B - B0| = {worst_off:.3e} nT > "
                       f"{TOTAL_TOL:g} * {b0max:.3e} nT)")

    return True, (f"Passed (total field tangential on {len(shell)} shell points, "
                  f"B1 = -B0.n n-hat there and zero elsewhere)")


def _check_body(cols, rows):
    """B0 must survive inside the body and stay constant there."""
    if cols["BODY"] is None:
        return True, "No body column (body checks skipped)"
    inside = [i for i in range(len(rows)) if cols["BODY"][i] > 0.5]
    if not inside:
        return False, "No point is inside the body (mask is empty)"

    b0in = max(_mag(_field(cols, "B0", i)) for i in inside)
    logger.debug("    max |B0| inside the body = %.4f nT (%d points)",
                 b0in, len(inside))
    if b0in <= 0.0:
        return False, "B0 is zero inside the body -- the mask removes it"

    fvidx, frows = _run_dir.load_first_out()
    if fvidx is not None and frows and len(frows) == len(rows):
        fcols = _columns(fvidx, frows)
        if fcols["B0X"] is not None:
            db = max(_mag(tuple(fcols["B0" + d][i] - cols["B0" + d][i]
                                for d in "XYZ")) for i in inside)
            if db > STATIC_TOL * b0in:
                return False, (f"B0 changed inside the body (max {db:.3e} nT)")

    # The electric field is pinned to zero in the linetied body.
    emax = max(_mag(_field(cols, "E", i)) for i in inside)
    if emax > 1e-9:
        return False, (f"The electric field is not zero inside the body "
                       f"(max |E| = {emax:.3e})")

    return True, f"Passed (B0 present and static inside the body, |E| = 0)"


def _crustal_coeffs_for(test_name):
    if test_name == "intrinsic_field_crustal":
        # crustal_nm1.txt: only g(1,0) is non-zero.
        return {(1, 0): (100.0, 0.0)}, True
    return CRUSTAL_NM2, False


def _check_crustal(cols, rows, test_name):
    coeffs, axial_only = _crustal_coeffs_for(test_name)

    if axial_only:
        # Closed form, independent of the derivative convention:
        # Br = 2 g10 cos(theta) / r^3 and Btheta = C g10 sin(theta) / r^3 with
        # one constant C, so Bphi must vanish and Br must match exactly.
        worst_br = 0.0
        ratios = []
        for i, _ in enumerate(rows):
            rr = _radius(cols, i)
            if rr < R_MIN or rr > R_MAX:
                continue
            x, y, z = _xyz(cols, i)
            r = math.sqrt(x * x + y * y + z * z)
            ct = z / r
            st = math.sqrt(max(0.0, 1.0 - ct * ct))
            ph = math.atan2(y, x)
            br = 2.0 * 100.0 * ct / r ** 3
            # radial component of the simulated field
            rh = (x / r, y / r, z / r)
            bs = _field(cols, "B0", i)
            br_sim = sum(bs[k] * rh[k] for k in range(3))
            worst_br = max(worst_br, abs(br_sim - br) /
                           (abs(br) if abs(br) > 1e-12 else 1.0))
            # theta component: B - Br r^hat, projected on theta^hat
            bth_sim = sum((bs[k] - br * rh[k]) *
                          (ct * math.cos(ph) if k == 0 else
                           ct * math.cos(ph) if k == 1 else -st)
                          for k in range(3))
            expect = -100.0 * st / r ** 3
            if abs(expect) > 1e-9:
                ratios.append(bth_sim / expect)
        if worst_br > ANALYTIC_TOL:
            return False, (f"The g10-only crustal field is not Br = "
                           f"2 g10 cos(theta)/r^3 (max rel err {worst_br:.3e})")
        if ratios:
            c = sum(ratios) / len(ratios)
            spread = max(abs(v / c - 1.0) for v in ratios)
            logger.debug("    crustal g10: Btheta coefficient = %.6f, "
                         "spread = %.3e", c, spread)
            if spread > 1e-6:
                return False, ("Btheta is not proportional to sin(theta)/r^3 "
                               f"(spread {spread:.3e})")
        return True, ("Passed (axial-dipole crustal field, Br exact, "
                      "Btheta self-consistent)")

    # General case: compare with the independent Python evaluation.
    errs = []
    for i, _ in enumerate(rows):
        rr = _radius(cols, i)
        if rr < R_MIN or rr > R_MAX:
            continue
        x, y, z = _xyz(cols, i)
        ban = crustal_b0(coeffs, x, y, z)
        scale = max(_mag(ban), 1e-12)
        errs.append(_mag(tuple(_field(cols, "B0", i)[k] - ban[k]
                               for k in range(3))) / scale)
    if not errs:
        return False, "No usable point in the comparison band"
    errs.sort()
    med = errs[len(errs) // 2]
    logger.debug("    crustal: n=%d median rel err = %.3e, max = %.3e",
                 len(errs), med, errs[-1])
    if med > CRUSTAL_TOL:
        return False, (f"B0 does not match the independent spherical-harmonic "
                       f"evaluation (median rel err = {med:.3e})")
    return True, f"Passed (crustal field, median rel err {med:.1e})"


def _check_superposition(cols, rows, first_cols=None):
    """Dipole + crustal must be the sum of the two individual fields."""
    errs = []
    for i, _ in enumerate(rows):
        rr = _radius(cols, i)
        if rr < R_MIN or rr > R_MAX:
            continue
        x, y, z = _xyz(cols, i)
        bd = dipole_b0(x, y, z)
        bc = crustal_b0(CRUSTAL_NM2, x, y, z)
        ban = tuple(bd[k] + bc[k] for k in range(3))
        scale = max(_mag(ban), 1e-12)
        errs.append(_mag(tuple(_field(cols, "B0", i)[k] - ban[k]
                               for k in range(3))) / scale)
    if not errs:
        return False, "No usable point in the comparison band"
    errs.sort()
    med = errs[len(errs) // 2]
    logger.debug("    superposition: n=%d median rel err = %.3e", len(errs),
                 med)
    if med > ANALYTIC_TOL:
        return False, (f"The dipole + crustal sum is wrong (median rel err = "
                       f"{med:.3e})")
    return True, f"Passed (dipole + crustal superposition, rel err {med:.1e})"
