# AMR Plasma Thermal Equilibrium Test — True 3D

Same uniform-plasma thermal equilibrium as `tests/amr_equilibrium`, but with a
**real z extent** (`nCellZ = 3`), so the 3D branch of the hybrid coarse-fine
treatment actually runs.

## Why this test exists

`relax_covered_B_to_fine` used to start with

```cpp
if (nDim != 2) { MultiFab::Copy(B, target, ...); return; }   // plain average_down
```

and `nDim` is the compile-time `amrex::SpaceDim`. In the default 3D build the
whole divergence-free relaxation was therefore dead code, and the covered
coarse cells were overwritten with a plain `average_down`, which reintroduced
an O(1e-4) div(B) error at the coarse-fine interface.

Nothing in the suite caught this, because every other AMR hybrid deck is
**fake-2D** (`nCellZ == 1`) — `isFake2D` routes those into the same scalar-
potential branch a true-2D build would use, so the 3D path was never executed.

With `nCellZ = 3` this deck gets `isFake2D == false` and `nDim == 3`, hence
`is_bz_div_free() == false` and the full nodal **vector-potential** relaxation
runs.

## Physical setup

Identical to `tests/amr_equilibrium/PARAM.in.hybrid`: kinetic ions
($q = m = 1$, $n_0 = 1\ \mathrm{amu/cc}$, $T_i = 10\ \mathrm{eV}$), massless
isothermal electron fluid ($T_e = 10\ \mathrm{eV}$, $\gamma_e = 1$), uniform
guide field $B_{z0} = 1\ \mathrm{nT}$, and a uniform drift $u_x = 100\ \mathrm{km/s}$.

## Mesh

- Domain $x \in [-16, 16]$, $y \in [-8, 8]$, $z \in [-1.5, 1.5]$, periodic in
  all three directions.
- Level 0: $32 \times 16 \times 3$ cells, $\Delta x = \Delta y = \Delta z = 1$.
- Level 1: central slab $|x| \le 8$, $|y| \le 4$ (full z) refined by 2, so
  $\Delta_1 = 0.5$.
- `nCellZ` is odd so that $z = 0$ is a cell centre and the `z=-0.5` output plane
  is non-empty; `dz = 1` keeps the mesh isotropic.
- Particle count is reduced to $8 \times 8 \times 2$ per cell to keep the 3D
  run cheap (the z plane doubles the cell count versus the fake-2D deck).

## Validation

`validate.py` reuses `tests/amr_equilibrium/validate.py` for the equilibrium
checks (energy conservation, coarse/fine grid detection, flat density across
the interfaces, reflection symmetry, B and u drift) and adds the gate that was
missing:

- **`#DIVB`** is on with `alwaysComputeDivB = T`, so every cycle prints a
  `divB-AMR` line with the max |div(B)| per level, split into interface /
  covered / interior / domain-edge buckets.
- The validator parses those lines and requires the interface, covered and
  interior buckets on every level to stay below `1e-10`. With the
  coarse-fine treatment working they sit at ~1.6e-16; if the relaxation ever
  becomes unreachable again, div(B) jumps to ~1e-4 and CI goes red.

The domain is fully periodic, so the `dom` (physical-boundary) bucket is zero
and cannot mask the reading.

Two checks are calibrated on the measured 3D behaviour rather than inherited
verbatim from the 2D sibling:

- **Energy** is measured from cycle 1. Step 1 seeds the electron pressure and
  the ambipolar E, which moves `Etot` once by ~0.09% and then leaves it flat to
  ~0.007%. The 2D sibling's 0.1% budget is almost entirely consumed by that
  one-off transient; skipping it allows a 10x tighter bound.
- **Density** bins within `INTERFACE_MARGIN = 2.0` of the slab edge are
  skipped, because a node on the coarse-fine interface averages in covered
  coarse cells that hold no particles and reports ~0.32 instead of ~1. Away
  from the interface the profile is flat to ~0.06, well inside the 0.10
  tolerance (that floor is a uniform pattern, not shot noise — quadrupling the
  particle count does not reduce it).

## Running

```bash
python3 tests/validate_tests.py --test=amr_equilibrium_3d -v
```

Requires a true-3D AMReX library; the test is in `AMREX3D_TESTS`, so it is
SKIPPED (not failed) under a `-amrex2d` build.
