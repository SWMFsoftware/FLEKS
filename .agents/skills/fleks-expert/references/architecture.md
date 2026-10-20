# FLEKS Architecture Reference

> Canonical: `doc/Algorithm.tex` (mathematics); the class hierarchy is mirrored
> in `references/file-layout.md`. This page adds the solver comparison that
> decides which path a run takes.

## Class hierarchy

```
Domain                        — Top-level simulation manager
 ├── DomainGrid               — Grid information container
 ├── Pic : Grid : AmrCore     — PIC solver (fields + particle push on AMR grid)
 ├── ParticleTracker : Grid   — Test particle tracker
 ├── FluidInterface : Grid    — MHD/fluid state on grid (coupling data)
 ├── SourceInterface          — Source terms for coupling
 │    └── UserSource          — User-defined source (selected from userfiles/)
 ├── OHInterface              — Outer-heliosphere coupling data
 └── TimeCtr                  — Time stepping, event control, plot scheduling
      ├── EventCtr            — Periodic event trigger (dn or dt based)
      └── PlotCtr             — Plot scheduling (combines EventCtr + PlotWriter)
```

Design rule: ionization parameters live in **SourceInterface** (physics in
`UserSource`), never in `FluidInterface`, so the MHD coupling layer stays
uncluttered.

## Core classes

| Class | Header | Purpose |
|---|---|---|
| `Domain` | `Domain.h` | Top-level orchestrator: owns `Pic`, `FluidInterface`, `TimeCtr` |
| `Pic` | `Pic.h` | PIC solver: field solve, particle push, moments, I/O |
| `Grid` | `Grid.h` | AMR grid management (inherits `AmrCore`) |
| `Particles` | `Particles.h` | Templated particle container (`PicParticles`, `PTParticles`) |
| `TestParticles` | `TestParticles.h` | Test particle class (trajectory recording, I/O) |
| `FluidInterface` | `FluidInterface.h` | Fluid/MHD state variables on the PIC grid |
| `LinearSolver` | `LinearSolver.h` | GMRES Krylov solver for the implicit E field |
| `TimeCtr` | `TimeCtr.h` | Time step management, CFL, event scheduling |
| `PlotWriter` | `PlotWriter.h` | Output formatting (IDL, AMReX, HDF5, VTK, Tecplot) |
| `DataContainer` | `Converter/DataContainer.h` | Data reading (IDL, AMReX formats) — converter only, not used by the solver |
| `GridUtility` | `GridUtility.h` | Discrete operators (curl, div, grad, averaging) |
| `BC` | `BC.h` | Boundary condition types |

## Field solvers

FLEKS provides three simulation field solver modes selected via `#FIELDSOLVER` (plus static field mode `none`). By default, `#FIELDSOLVER implicit` is active:

### Full (implicit) PIC (default: `#FIELDSOLVER implicit`)

- Kinetic ions *and* electrons, semi-implicit θ-scheme (default θ = 0.51) with
  a Boris mover.
- The implicit electric field is solved every step by the GMRES
  `LinearSolver`.
- Most robust for strongly kinetic phenomena: reconnection, shock surfing,
  light waves with dt larger than explicit CFL.
- Key files: `src/PicFieldSolver.cpp`, `src/PicDivE.cpp`.
- Configured by `#DISCRETIZATION` and `#EFIELDSOLVER`.

### Hybrid PIC (`#FIELDSOLVER hybrid`)

- **Kinetic ions + massless fluid electrons**: electrons are a
  charge-neutralizing fluid governed by a generalized Ohm's law (Hall term,
  electron pressure, resistivity, hyper-resistivity), removing the electron
  time-scale.
- Fields advance explicitly with an **RK4/SSPRK3 Faraday advance** — no GMRES
  solve.
- Enabled by `#FIELDSOLVER hybrid` (legacy `#HYBRIDPIC` deprecated); tuned with
  `#HALLTERM`, `#RESISTIVITY`, `#HYPERRESISTIVITY`, `#ELECTRON*`,
  `#FIELDINTEGRATOR` (`rk4`, `ssprk3`, `rk4-pc`, `ssprk3-pc`),
  `#FIELDSUBCYCLE` (legacy `#BSUBCYCLE` deprecated), `#MINIMUMDENSITY`,
  `#SMOOTHMOMENTS`, `#DIVE`/`#DIVB`.
- Key file: `src/PicHybrid.cpp`. Standalone tests ship a `PARAM.in.hybrid`
  variant; the runner reports all solver variants when present.

### Explicit Full PIC (`#FIELDSOLVER explicit`)

- **Kinetic ions and kinetic electrons**: resolves full kinetic physics including
  electron scales without requiring matrix inversions or Krylov solvers.
- **Explicit Maxwell evolution**: discrete dual staggering where `nodeE` (vertices)
  and `centerB` (cell centers) are evolved explicitly via discrete curl operators
  `curl_node_to_center` (Faraday) and `curl_center_to_node` (Ampere).
- **Exact discrete divergence preservation**: $\nabla \cdot \mathbf{B} = 0$ is
  preserved to machine precision by construction across all AMR levels.
- **Integrator options**: `#FIELDINTEGRATOR`
  * `leapfrog` (default): Subcycled, time-staggered Boris-Maxwell leapfrog.
  * `pc` / `rk2`: 2nd-order predictor-corrector with intermediate particle push.
- **Subcycling**: `#FIELDSUBCYCLE` sets the number of explicit Maxwell subcycles
  per particle time step dt.
- **Current**: Ampere's law takes the charge current sum_s q_s n_s u_s
  (`nodeJexpl`, `deposit_explicit_current`), not the mass-flux moments in
  `nodePlasma`. Leapfrog: kick v^{n+1/2} with (E^n, B^n), deposit at x^n,
  drift, deposit at x^{n+1}; J^{n+1/2} is the average. Predictor-corrector:
  deposits at t^n and t^{n+1}, interpolated in time across subcycles.
  `tests/thermal_plasma` checks energy conservation.
- **Zero regression / minimal footprint**: bypasses allocation of GMRES scratch
  MultiFabs and mass matrices when active, leaving semi-implicit and hybrid
  allocations and performance untouched.
- Enabled by `#FIELDSOLVER explicit`; key file: `src/PicExplicit.cpp`. Standalone
  tests ship a `PARAM.in.explicit` variant.

## Algorithm summary

- **Time stepping:** CFL-based or fixed Δt (`#TIMESTEPPING`, typical CFL 0.1–0.4).
- **Divergence cleaning:** accurate div(E) correction via particle position
  adjustment; optional hyperbolic div(B) cleaning.
- **Particle management:** splitting, merging, fast merge with Lagrange
  multipliers (`#RESAMPLING`, `#FASTMERGE`).
- **AMR:** block-structured AMR via AMReX, configurable levels and regions
  (`./Config.pl -lev=N`, `#REGION`, `#AMR`); dynamic load balancing through
  `FleksDistributionMap`.
- **Coupling:** bi-directional GM↔PC exchange (MHD-AEPIC); see
  `references/coupling.md`.
- Mathematical foundations for both solvers: `doc/Algorithm.tex`.
