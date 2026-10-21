# Magnetic Reconnection Standalone Tests

This directory contains standalone magnetic reconnection test suites in FLEKS across multiple physical configurations and field solvers:

1. **`PARAM.in.fadeev_pic`** — **Fadeev full-PIC**: Force-free Fadeev current-sheet island equilibrium with kinetic ions and kinetic electrons ($m_i/m_e = 25$, Maxwell/GMRES solver).
2. **`PARAM.in.fadeev_hybrid`** — **Fadeev hybrid-PIC**: Kinetic ions + massless fluid electrons with generalized Ohm's law.
3. **`PARAM.in.gem_pic`** — **Classic GEM Challenge full-PIC**: Standard GEM reconnection benchmark (Birn et al. 2001) with a Harris current sheet, conducting walls in $y$, and a central magnetic perturbation.
4. **`PARAM.in.gem_hybrid`** — **Classic GEM Challenge hybrid-PIC**: Kinetic ions + isothermal fluid electrons with generalized Ohm's law, stationary ion background (`useUniformIonPressure = T`), and electron fluid carrying the diamagnetic current.
5. **`PARAM.in.asym_pic`** — **Asymmetric full-PIC**: Double current sheet reconnection with asymmetric magnetic fields ($B_1 = 1.0, B_2 = 2.0$) and temperatures ($T_1 = 1.33, T_2 = 3.33$) in a periodic domain.
6. **`PARAM.in.forcefree_hybrid`** — **Force-Free Sheet hybrid-PIC**: Force-free current sheet reconnection (Le et al. 2016, WarpX benchmark) with uniform plasma density, uniform total magnetic pressure, and fixed 40 subcycles.
7. **`PARAM.in.forcefree_hybrid_adaptive`** — **Force-Free Sheet hybrid-PIC (Adaptive)**: Same force-free setup with CFL-driven adaptive macro time-stepping (`#TIMESTEPPING F 0.25`) and automatic magnetic subcycling (`#BSUBCYCLE 1 T 40`). Fast CI benchmark.
8. **`PARAM.in.forcefree_pic`** — **Force-Free Sheet full-PIC**: Force-free current sheet reconnection with kinetic ions and electrons ($m_i/m_e = 25$, Maxwell/GMRES solver).
9. **`PARAM.in.fadeev_hybrid_upwind_long`** — **Fadeev hybrid-PIC through coalescence**: `fadeev_hybrid` run to $t = 120$ with the upwind dissipation `#HYBRIDUPWIND` in place of hyper-resistivity (see [Fadeev](#1-fadeev-equilibrium-paraminfadeev_pic-paraminfadeev_hybrid)).

## Coordinate Mapping

FLEKS uses a fake-2D convention with 1 cell along $z$:
- $x$: Reconnection outflow / periodic drive direction
- $y$: Current-sheet normal direction
- $z$: Out-of-plane / current / guide-field direction

## Physical Configurations

Physical parameters, grid dimensions, boundary conditions, and timestepping details are documented in comments inside each `PARAM.in.*` file.

### 1. Fadeev Equilibrium (`PARAM.in.fadeev_pic`, `PARAM.in.fadeev_hybrid`)
Uses `#TESTCASE fadeev` (`FadeevIC`): two-island equilibrium with magnetic field nulls and density concentrated in the current sheet. The sheet is bounded by outflow field boundaries in $y$.

Reconnection is measured by the island flux $\psi = |A_z(X) - A_z(O)|$, with
$A_z = -\int B_y\,dx$ on the midplane and the X- and O-points at the midplane
$B_y$ nulls nearest $x = 0$ and $x = \pm\pi L$. $\psi \approx 3$ at $t = 0$;
it must change by more than 0.02 (0.04–0.06 in the hybrid runs by $t = 20$,
0.16 in the full-PIC run by $t = 3$). By $t = 20$ the islands have hardly
started to coalesce, so `fadeev_hybrid` cannot tell dissipation schemes apart.

`fadeev_hybrid_upwind_long` runs through the coalescence to $t = 120$ (10–15
minutes on 2 ranks) and must lose more than 0.2 of $\psi$ (0.49 measured).
Without dissipation the run fails at $t = 73$, and with the `fadeev_hybrid`
hyper-resistivity (`etaHyperCh = 0.001`) at $t = 116$: plasma drains through
the open $y$ boundaries (the ion energy falls by 95 %), density holes form
(6×6 ions per cell) and $|B|/\rho$ runs away. Grid-mode hyper-resistivity
with `etaHyperCh = 0.0035`, which damps resolved modes about as much as the
upwind term at this resolution, also survives, with the island flux within
7 % of the upwind run; the deck's 0.001 is simply too weak for this length.

### 2. Classic GEM Challenge (`PARAM.in.gem_pic`, `PARAM.in.gem_hybrid`)
Uses `#TESTCASE gem` with `useStandardGem = T`: Harris current sheet with periodic boundaries in $x$ and conducting/reflecting walls in $y$, following the Birn et al. 2001 benchmark specifications.

### 3. Asymmetric Double Current Sheet (`PARAM.in.asym_pic`)
Uses `#TESTCASE gem` with `isAsymmetryReconnection = T`: two current sheets at $y = \pm 7\,d_i$ with asymmetric lobe fields and temperatures, fully periodic in $x$ and $y$.

### 4. Force-Free Current Sheet (`PARAM.in.forcefree_hybrid`, `PARAM.in.forcefree_pic`)
Uses `#TESTCASE forcefree` (`ForceFreeIC`), ported from the WarpX benchmark (Le et al. 2016). Uniform total magnetic pressure ($B_x^2 + B_z^2 = \text{const}$) so no plasma pressure gradient is needed for mechanical equilibrium ($\mathbf{J} \times \mathbf{B} = 0$).

## Running

Run all reconnection test variants together:
```bash
# Serial (default)
python3 tests/validate_tests.py --test=reconnection

# MPI (e.g. 2 ranks)
python3 tests/validate_tests.py --test=reconnection -n 2
```

Run a single variant:
```bash
python3 tests/validate_tests.py --test=reconnection.fadeev_pic        # Fadeev full-PIC
python3 tests/validate_tests.py --test=reconnection.fadeev_hybrid     # Fadeev hybrid-PIC
python3 tests/validate_tests.py --test=reconnection.gem_pic           # Classic GEM challenge full-PIC
python3 tests/validate_tests.py --test=reconnection.gem_hybrid        # Classic GEM challenge hybrid-PIC
python3 tests/validate_tests.py --test=reconnection.asym_pic          # Asymmetric reconnection full-PIC
python3 tests/validate_tests.py --test=reconnection.forcefree_hybrid  # Force-free sheet hybrid-PIC (fixed 40 subcycles)
python3 tests/validate_tests.py --test=reconnection.forcefree_hybrid_adaptive # Force-free sheet hybrid-PIC (adaptive)
python3 tests/validate_tests.py --test=reconnection.forcefree_pic     # Force-free sheet full-PIC
python3 tests/validate_tests.py --test=reconnection.fadeev_hybrid_upwind_long -n 2  # Fadeev hybrid-PIC to t = 120
```
*(Note: only the two force-free variants (`forcefree_hybrid`, `forcefree_pic`) run in the
default full suite. The Fadeev, GEM and asymmetric variants are marked expensive and skipped;
run them by explicitly specifying e.g. `--test=reconnection.gem_pic` or by adding
`--include-expensive` / `--all`.)*

## Validation

Validation logic and per-check documentation live in [`validate.py`](validate.py). Each `_validate_*` function docstring describes the exact checks and thresholds applied to that configuration.
