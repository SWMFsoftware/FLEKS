# Intrinsic Magnetic Field

A static planetary magnetic field declared with `#DIPOLE` (analytic dipole)
and/or `#CRUSTALFIELD` (spherical harmonic expansion read from a file, in the
two BATSRUS layouts). The field is frozen: it never enters Faraday's law and it
is not written to the restart files. It is added to the evolved field B1
wherever a *total* magnetic field is needed, B = B1 + B0: the particle push,
the mass matrix, the generalized Ohm's law (convective and Hall terms), the
magnetic energy and the plot output.

The two are stored apart so that the div(B) cleaning and the upwind correction,
which act on `nodeB`/`centerB`, never touch B0.

## Variants

| File | Checks |
|------|--------|
| `PARAM.in` | Tilted dipole, field solver off, `B1 = 0` |
| `PARAM.in.body` | Dipole + absorbing `#BODY`: B0 must survive inside the body |
| `PARAM.in.crustal` | Crustal field, new BATSRUS layout |
| `PARAM.in.crustal_nm2` | Crustal field with several degrees and both g and h |
| `PARAM.in.crustal_old` | Same coefficients, legacy `marsmgsp` layout |
| `PARAM.in.dipole_crustal` | Dipole + crustal must be their superposition |

All variants switch the field solver off (`#SOLVEEM F`) and set the
`#UNIFORMSTATE` field to zero, so `B1 = 0` and the plotted field *is* B0. That
makes the output comparable point by point with an independent evaluation of
the model, and any leakage of B0 into the evolved field (or a failure to add it
to the total) shows up immediately.

The tests are cheap: the whole suite runs in about a minute.

## Checks

1. **Static**: B0 in the first and the last frame agree to `1e-6` relative.
2. **Total = B0**: because `B1 = 0`, `B` and `B0` must agree to `1e-10`.
3. **Surface round trip**: at `r = 1` (planetary radii) on the magnetic
   equator, `|B0|` equals the `#DIPOLE` strength in nT. The plot unit is
   `planet`, so B comes out in nT and this is an identity check of the
   SI -> code -> SI conversion.
4. **Dipole shape**: `r^-3` scaling, `Br/Btheta = 2 cot(theta)`, and the
   recovered dipole axis matches the requested tilt.
5. **Crustal**: compared against an independent spherical harmonic evaluation
   in `validate.py`; `crustal_nm2` and `crustal_old` must be identical.
6. **Superposition**: dipole + crustal is the sum of the two.
7. **Body**: B0 is non-zero and unchanged inside the body mask.
8. **Energies**: `Etot` and `Eb` stay constant (the test particles gyrate, and
   a magnetic field does no work). With `#BODY` the absorbed particles leave
   the domain, so only `Eb` is checked.

## Units and reference radius

The two length conventions of the standalone decks apply unchanged:

- `#BODY` radius/center and the `#GEOMETRY` box are in **code units**
  (multiples of `#NORMALIZATION lNormSI`);
- `#PLANETRADIUS` and the `#DIPOLE` strength/`rRef` are **SI**.

The dipole reference radius resolves as: an explicit `rRef > 0` in `#DIPOLE`,
else the `#BODY` radius, else `#PLANETRADIUS`. Every deck here sets
`#PLANETRADIUS` to the reference radius, so one code unit of length is one
planetary radius and the reference sphere sits at `r = 1` in the plot output.

## Known limitations

- `#BODYBOUNDARY fieldBoundary = conducting` cannot be combined with an
  intrinsic field: it pins the radial field of the *evolved* part to zero on
  the body surface, which contradicts the radial field of a planet. FLEKS
  aborts on the combination; use `linetied` or `insulating`.
- The hybrid solver uses the total field only for the convective and Hall
  terms of the Ohm's law. The current is computed from B1 alone, because the
  intrinsic field is current-free and the discrete curl of `B1 + B0` would
  inject its truncation error into J. A dedicated hybrid test is left for
  future work.
