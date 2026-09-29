# Intrinsic Magnetic Field (`#DIPOLE`, `#CRUSTALFIELD`)

## Description

This test verifies the static intrinsic magnetic field of a planet, declared
with `#DIPOLE` (analytic dipole) and/or `#CRUSTALFIELD` (spherical harmonic
expansion read from a file, matching the BATSRUS `ModUserMars` layout). The
field is frozen: it never enters Faraday's law and it is not written to the
restart files. It is added to the evolved field B1 wherever a *total* magnetic
field is needed, B = B1 + B0: the particle push, the mass matrix, the
generalized Ohm's law (convective and Hall terms), the magnetic energy and the
plot output.

The two are stored apart so that the div(B) cleaning and the upwind correction,
which act on `nodeB`/`centerB`, never touch B0.

Every variant switches the field solver off (`#SOLVEEM F`) and sets the
`#UNIFORMSTATE` field to zero, so `B1 = 0` and the plotted field *is* B0. That
makes the output comparable point by point with an independent evaluation of
the model, and any leakage of B0 into the evolved field (or a failure to add it
to the total) shows up immediately.

## Variants

| Variant | File | Physics | Validation |
|---------|------|---------|------------|
| `intrinsic_field` | `PARAM.in` | Tilted dipole: 100 nT on the reference sphere `r = 1`, axis 30° from `+z` towards `-x`; `B1 = 0` | Median relative error against the analytic dipole `< 1e-6`; `\|B0\| r^3` constant; `\|B0\| = 100 nT` at the magnetic equator on `r = 1` |
| `intrinsic_field_body` | `PARAM.in.body` | Same dipole plus an absorbing `#BODY` (radius 1, `absorb` + `linetied`) | `B0` non-zero and unchanged inside the body mask, `\|E\| = 0` there; only `Eb` checked in the log (absorbed particles leave the domain) |
| `intrinsic_field_crustal` | `PARAM.in.crustal` | Crustal field `crustal_nm1.txt`, `nMax = 3`, only `g(1,0) = 100 nT` non-zero | Closed-form axial dipole: `Br = 2 g10 cos(theta)/r^3` exact, `Btheta ∝ sin(theta)/r^3`, no `Bphi` |
| `intrinsic_field_crustal_nm2` | `PARAM.in.crustal_nm2` | Crustal field `crustal_nm2.txt`, `nMax = 3`, several degrees with both `g` and `h` | Matches an independent Schmidt semi-normalized evaluation in `validate.py` (median relative error `< 1e-8`) |
| `intrinsic_field_dipole_crustal` | `PARAM.in.dipole_crustal` | Dipole + crustal field together | Linear superposition: `B0` equals the sum of the dipole-only and crustal-only fields (median relative error `< 1e-6`) |

All variants run in about a minute in total.

## Physics & Solver Setup

- **Geometry & Boundaries**: `32 × 32 × 4` cells over `[-4, 4] × [-4, 4] ×
  [-2, 2]` code units (`dx = 0.25`), periodic in every direction, blocks of
  `16 × 16 × 4`. Four cells in z put `z = 0` exactly on a node plane, which is
  the plane the plot writes; under a 2D AMReX build (`-amrex2d`) the grid
  collapses to a single cell in z and the same checks apply.

- **Reference radius**: `#PLANETRADIUS` is `1.0e5 m` and `lNormSI` is the same
  value, so one code unit of length is one planetary radius and the reference
  sphere sits at `r = 1` in the plot output. The dipole reference radius
  resolves as: an explicit `rRef > 0` in `#DIPOLE`, else the `#BODY` radius,
  else `#PLANETRADIUS`; every deck here leaves `rRef < 0`.

- **Plasma Species**:

  | Species | Mass [amu] | Charge [e] | n [amu/cc] | u [km/s] | T [K] | Role |
  |---------|------------|------------|------------|----------|-------|------|
  | 0 | 1.0 | +1 | 5.0 | 0 (20 in y for `PARAM.in`) | 1.0e4 | Ions |
  | 1 | 0.04 | −1 | 0.2 | 0 (20 in y for `PARAM.in`) | 1.0e4 | Electrons (n_e = n_i) |

- **Electromagnetic Fields**: `#SOLVEEM F`, so the evolved field stays at its
  initial value (zero) and `E = 0`: the only field acting on the particles is
  the intrinsic one. The dipole deck gives the plasma a 20 km/s drift in y so
  that the test particles gyrate.

- **Output**: `planet` units, so the magnetic field comes out in nT and the
  coordinates in planetary radii. Two frames are written (`t = 0` and
  `t = 10 s`) carrying `B0x B0y B0z Bx By Bz Ex Ey Ez`, plus `body` in the body
  deck.

- **Time Stepping**: adaptive (`useFixedDt = F`, `cfl = 0.2`), `TimeMax = 10 s`,
  `2 × 2 × 1` particles per cell.

## Validation

From the pic-log history:

1. all energies stay finite, and `Eb` is constant — a static field carries a
   constant magnetic energy. With `#BODY` the absorbed particles leave the
   domain, so only `Eb` is checked there; elsewhere `Etot` must also stay
   constant to `1e-6` (a magnetic field does no work).

From the plot frames:

2. **Static**: `B0` in the first and the last frame agree to `1e-6` relative.
3. **Total = B0**: because `B1 = 0`, `B` and `B0` must agree to `1e-10`.
4. **Surface round trip**: at `r ≈ 1` on the magnetic equator, `|B0|` equals
   the `#DIPOLE` strength in nT — an identity check of the
   SI → code → SI conversion.
5. **Dipole shape**: `r^-3` scaling, `Br/Btheta = 2 cot(theta)`, and the
   recovered dipole axis matches the requested tilt.
6. **Crustal**: compared against an independent spherical harmonic evaluation
   in `validate.py`.
7. **Superposition**: dipole + crustal is the sum of the two.
8. **Body**: `B0` is non-zero and unchanged inside the body mask, and `E = 0`
   there.

## Known Limitations

- `#BODYBOUNDARY fieldBoundary = conducting` cannot be combined with an
  intrinsic field: it pins the radial field of the *evolved* part to zero on
  the body surface, which contradicts the radial field of a planet. FLEKS
  aborts on the combination; use `linetied` or `insulating`.
- The hybrid solver uses the total field (B1 + B0) for the convective and Hall
  terms of the Ohm's law, while the current J is computed from B1 alone because
  the intrinsic field is current-free. Both the hybrid PIC solver and full PIC
  solver support intrinsic magnetic fields and inner body boundaries.

## Running

From the FLEKS root directory (requires compiled `bin/FLEKS.exe`):

```bash
python3 tests/validate_tests.py --test=intrinsic_field                  # all variants
python3 tests/validate_tests.py --test=intrinsic_field.crustal          # single variant
python3 tests/validate_tests.py --test=intrinsic_field.dipole_crustal
```
