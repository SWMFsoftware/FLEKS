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

### Full (implicit) PIC

- Kinetic ions *and* electrons, semi-implicit θ-scheme (default θ = 0.51) with
  a Boris mover.
- The implicit electric field is solved every step by the GMRES
  `LinearSolver`.
- Most robust for strongly kinetic phenomena: reconnection, shock surfing,
  light waves.
- Key files: `src/PicFieldSolver.cpp`, `src/PicDivE.cpp`.

### Hybrid PIC

- **Kinetic ions + massless fluid electrons**: electrons are a
  charge-neutralizing fluid governed by a generalized Ohm's law (Hall term,
  electron pressure, resistivity, hyper-resistivity), removing the electron
  time-scale.
- Fields advance explicitly with an **RK4/SSPRK3 Faraday advance** — no GMRES
  solve.
- Enabled by `#HYBRIDPIC`; tuned with `#HALLTERM`, `#RESISTIVITY`,
  `#HYPERRESISTIVITY`, `#ELECTRON*`, `#FIELDINTEGRATOR`, `#BSUBCYCLE`,
  `#MINIMUMDENSITY`, `#AVGFIELDB`, `#SMOOTHMOMENTS`, `#DIVE`/`#DIVB`.
- Key file: `src/PicHybrid.cpp`. Standalone tests ship a `PARAM.in.hybrid`
  variant; the runner reports both solver variants when one exists.

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
