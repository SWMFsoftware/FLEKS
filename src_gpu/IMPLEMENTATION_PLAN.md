# FLEKS GPU PIC demo implementation plan

Status: implemented and verified on local branch `gpu`. See `PROGRESS.md` for
completed stages and validation evidence, and `README.md` for runnable commands.

## Purpose and scope

Build a standalone AMReX PIC demo that exposes the field, particle, communication,
and solver work required to port FLEKS to GPUs. Keep the implementation small
enough to understand and profile while preserving the important FLEKS data
layouts and algorithm boundaries.

In the planned layout, all demo C++ sources and headers belong in `src_gpu/`.
The executable will remain
`bin/main_gpu.exe`, built with `make GPUTEST` from the repository root or
`make -C src_gpu GPUTEST`. The demo links to AMReX and the configured C++/CUDA
libraries independently of `libFLEKS.a` and the SWMF coupling libraries.

The initial physical configuration is a uniform, periodic 3D Cartesian mesh,
one mesh level, configurable box size and particles per cell, and normalized
units. Prescribed-field tests use one species. The coupled PIC milestone uses
kinetic electrons and ions with a neutral initial state. AMR, physical walls,
SWMF coupling, reactions, resampling, and load balancing are later extensions.

## Stage 0: planned build separation

1. Commit the uncommitted changes in `src/main_gpu.cpp`, then move it with
   `git mv` to `src_gpu/main_gpu.cpp` so its history follows the file.
   Preserve the existing test behavior.
2. Add `src_gpu/Makefile` and route root `GPUTEST` to it. It includes
   `../Makefile.def` and `../Makefile.conf` like `src/Makefile` so that the
   nvcc compiler and flags, `AMREXLIB`, and `CPPLIB` match. Compile demo
   objects in `src_gpu/`; future demo `.cpp` and `.h` files stay there too.
   Generate header dependencies (`-MMD -MP` or `depend.pl`): the current
   `GPUTEST` rule tracks only `main_gpu.cpp`, which stops being sufficient once
   the demo has its own headers.
3. Remove the `GPUTEST` rule from `src/Makefile`, drop `main_gpu.cpp` from its
   `filter-out` list and `main_gpu.o` from its `clean` rule, and delete the
   stale `src/main_gpu.o`. Keep the production source list restricted to the
   production source directories.
4. Include `src_gpu/` in root clean operations, the compilation database
   (directory loop in `tools/generate_compile_commands.py`; `-I../include`
   remains correct one level below the root), and the default C++ formatter
   scope (`repo_dirs` in `tools/format_all.py`).
5. Verify build routing, incremental rebuilds after touching a demo header,
   and executable behavior where the available compiler, MPI runtime, and GPU
   permit.

Stage 0 is complete. The original `FabArrayGpuTest` checks are preserved in
`DemoChecks`; the stages below describe the implementation and acceptance.

## Planned source layout

Create files as their stages are implemented; do not add empty placeholder
classes. Device helpers and templates must be visible from their callers.

| File | Responsibility |
| --- | --- |
| `main_gpu.cpp` | AMReX initialization/finalization, mode dispatch, exit status |
| `DemoConfig.h`, `DemoConfig.cpp` | Read and validate runtime parameters |
| `DemoGrid.h`, `DemoGrid.cpp` | Own geometry, cell/nodal box arrays, distribution map |
| `DemoFields.h`, `DemoFields.cpp` | Own fields, moments, matrix storage, scratch arrays |
| `DemoFieldKernels.h` | Device-callable interpolation and discrete operators |
| `DemoParticles.h`, `DemoParticles.cpp` | One AMReX particle container per species |
| `DemoParticleKernels.h` | CIC weights, gather, Boris push, atomic deposit helpers |
| `DemoPic.h`, `DemoPic.cpp` | Species ownership and timestep orchestration |
| `DemoLinearSolver.h`, `DemoLinearSolver.cpp` | GPU-resident Krylov vector operations and solve |
| `DemoChecks.h`, `DemoChecks.cpp` | Field/particle reference checks and diagnostics |

Extract the current field-test class into `DemoChecks` during stage 1. Keep
`main_gpu.cpp` as the driver rather than growing it into the complete solver.

## Architecture and ownership

`DemoPic` owns a single `DemoGrid`, `DemoFields`, and a collection of
`DemoParticles` species. Consumers hold references to the grid; they do not
duplicate or mutate mesh metadata. Construct the grid before its consumers and
destroy it after them. Destroy every AMReX allocation before finalization.

Start with a single-level AMReX particle container. Keep mesh queries behind
`DemoGrid` so a later AMR version can adopt the production hierarchy interface
without changing particle physics kernels.

Keep the host as the control plane: parameters, allocation, iteration, MPI
orchestration, and reporting. Fields, particle data, mass matrices, scratch
vectors, and repeated numerical operations stay on the device. Copy only small
diagnostic scalars and explicit reference-check snapshots to the host.

Follow FLEKS naming and header conventions: PascalCase classes/files,
snake_case methods, camelCase variables, include guards, and no namespace-wide
using directives in headers. Match the existing nvcc restrictions: kernel
launching members are public where needed, and device lambdas capture views and
scalars by value rather than capturing `this`.

## Stage 1: preserve the baseline and introduce classes

Implement configuration, grid ownership, the field/check harness, and mode
dispatch. Preserve the existing field-only checks as the default until the
particle mode is ready.

Planned modes are `fields`, `particles`, `deposit`, `implicit`, and `benchmark`.
Existing arguments `n_cell`, `max_grid_size`, and `n_iter` remain usable.
Introduce `n_step`, `dt`, `particles_per_cell`, `seed`, and diagnostic intervals
only when their stages need them. Validate positive dimensions and iteration
counts and reject box-size choices that cannot support a requested check.

Acceptance:

- All current fill, reduction, ghost, arithmetic, stencil, copy, and host
  reference checks retain their PASS/FAIL behavior.
- The executable returns nonzero when a check fails on any rank.
- Both root and direct-directory builds use the new source directory.
- Production FLEKS targets do not compile or link demo objects.

## Stage 2: implement FLEKS-style field layouts and kernels

Allocate three-component nodal `nodeE` and `nodeB`, and cell-centered
`centerB`. Add nodal deposited quantities and persistent work arrays as needed.
Choose ghost widths from each operator's actual stencil radius, including
composed operations; document those requirements at the method boundary.
Production FLEKS allocates one uniform ghost width (`Grid::get_n_ghost()`).
Narrower per-operator widths reduce communication volume, so stage 6 reports
ghost-exchange cost for both the demo widths and the production width.

Implement manufactured initialization, cell-to-node B interpolation,
node-to-cell curl, cell-to-node curl, and field vector arithmetic. Translate
FLEKS's discrete operators rather than substituting a Yee-grid discretization.
Reference `include/GridUtility.h` and `src/PicFieldSolver.cpp`.

Make three distinct communication operations explicit: filling ghost values,
adding deposited contributions, and reconciling duplicated valid nodal values.
An additive deposit must be summed once before its result is used or copied to
ghost cells. A field computed independently on nodes must not be additively
summed merely because nodes overlap.

Acceptance:

- Manufactured discrete results agree with an independent host reference.
- Checks cover multiple boxes, periodic seams, and MPI boundaries.
- Nodal reductions use ownership masks, including periodic duplicates.
- Results remain consistent when only box size or rank count changes.
- Communication smoke test for the stage 5 mass-matrix type: allocate a small
  `NodeMMFab` (`FabArray<BaseFab<RealMM>>`, `include/Array1D.h`), add values
  to individual `Real` components with device atomics, call `SumBoundary` and
  `FillBoundary`, and compare with a host reference. This confirms early that
  AMReX device communication compiles and works for the custom value type; if
  it does not, decide the storage layout before stage 5.

## Stage 3: implement particles, gather, push, and migration

Retain an AoS baseline with the production particle type: FLEKS uses
`AmrParticleContainer<nPicPartReal, nPicPartInt>` = `<4, 1>`
(`include/Particles.h`, `include/Constants.h`), with velocity at
`iup_ = 0..2`, macro-particle weight/charge at `iqp_ = 3`, and one integer
attribute. Use the same template parameters and attribute indices, and use
`AmrParticleContainer` with a single level, so that gather, push, and deposit
kernels transfer to production unchanged. Keep species charge and mass
separate from weight and document the conversion to deposited charge and mass.
Use AMReX particle storage and migration rather than a custom device particle
vector. Before starting, confirm that the installed AMReX version still
supports AoS particle containers without deprecation or performance penalties;
the stage 6 SoA comparison may require `ParticleContainerPureSoA`, which
changes kernel access patterns.

Initialize particles from global cell indices, species, and within-cell sample
indices. Initialization must describe the same physical sample after changing
box decomposition. Device-safe random generation, if introduced, must retain
that property; otherwise begin with deterministic quiet-start samples.

Implement CIC gathering from the eight surrounding nodes, the Boris velocity
update, position advance, periodic wrapping, and `Redistribute()`. Initially
perform gather and push separately for diagnosis, storing the gathered E and B
in six runtime particle attributes or a per-tile scratch array; the production
mover fuses them, so stage 6 compares against a fused kernel. Before gathering, guarantee
that particles are assigned to the appropriate boxes and field ghosts are valid.
Redistribute after drift before any new-position operation uses box-local data.

Acceptance:

- Free streaming has the expected position and unchanged velocity.
- Uniform E gives the expected acceleration for positive and negative charge.
- Uniform B preserves speed to roundoff; orbit error converges with timestep.
- Particles deliberately cross box, rank, and periodic boundaries with no
  change in total count, weight, or identifying attributes.
- Tiny-case particle results agree with a serial reference.

## Stage 4: implement particle-to-grid moments

Zero deposition storage, accumulate CIC density/current contributions with
device atomics, perform additive boundary synchronization, normalize by cell
volume, and populate the required ghosts. Begin with density and three current
components; add the six independent second moments after this path is correct.
Document whether each component is number, charge, mass, or velocity weighted.

Reuse one definition of CIC weights for gather and deposition. The serial
reference should independently calculate expected values for a few particles
at cell centers, nodes, arbitrary positions, and periodic boundaries.

Ordinary `q*v` current moments test particle scatter and support diagnostics.
They do not establish discrete charge continuity for an explicit Maxwell
integrator. An explicit branch would require a trajectory-based
charge-conserving current deposit; port the existing FLEKS Esirkepov deposit
(declared in `include/Particles.h`, used by `src/PicExplicit.cpp`) and
validate it separately rather than writing a new one.

Acceptance:

- Integrated charge and current agree with direct particle sums using unique
  nodal ownership and the chosen normalization.
- Multiple particles depositing to the same nodes produce no lost updates.
- Changing decomposition changes results only within accumulation tolerances.
- CPU/GPU differences use precision- and problem-size-aware tolerances rather
  than requiring bitwise equality from atomic additions.

Milestone A is complete after stages 1-4: fields, particles, gather/push/deposit,
and migration work on the GPU with independent correctness checks. This remains
a prescribed-field particle demo until stage 5 closes the physical coupling.

## Stage 5: implement the semi-implicit FLEKS timestep

Audit the exact position/velocity time levels and theta-stage field conventions
in `src/Pic.cpp`, `src/ParticlesMover.cpp`, `src/ParticlesMassMatrix.cpp`, and
`src/PicFieldSolver.cpp`, and `src/PicDivE.cpp` before coding the coupled
loop. Record the resulting stage table next to the driver. A Boris test by
itself does not validate the time centering of the coupled FLEKS algorithm.

The default single-level implicit step in `Pic::update()` (`src/Pic.cpp`) is:
`calc_mass_matrix()` → `update_E()` → `advance_particles_and_sources(0)` →
`update_B()` → `divE_correction()` (on by default, `doCorrectDivE = true`) →
`sum_moments()`. Co-moving frame, upwind, hybrid, and AMR branches are out of
scope.

The mass matrix dominates memory and scatter cost. `RealMM` holds
`nMMComponents = 243` reals per node in 3D (27 neighbor nodes × 3 × 3), about
1.9 KB per node in double precision, or roughly 500 MB for 64^3 nodes before
ghosts. The full CIC scatter performs 8 destination nodes × 8 source nodes ×
9 tensor entries = 576 atomic adds per particle; production's symmetry
optimization reduces this work further. Size benchmark problems
against device memory, and time step 1 on its own as soon as it is correct
rather than waiting for stage 6.

Implement this stage in separately reviewable parts:

1. Accumulate `jHat` and the nodal particle-response mass matrix using the
   single-level FLEKS formulas. Validate contributions from a single particle
   and from particles sharing nodes against a serial reference.
2. Apply the mass matrix to a nodal field, including inter-node stencil
   couplings and ghosts. Start with the production custom FabArray/RealMM
   representation verified in stage 2; measure a MultiFab component
   representation later.
3. Assemble the electric-field RHS and complete matrix-vector operator with
   the actual FLEKS coefficients and discrete operators. Verify each term
   separately and use a manufactured solution for the assembled system.
4. Implement restarted GMRES with vectors and vector operations on the device.
   Host-side convergence decisions and a small Hessenberg system are allowed.
   Represent each Krylov vector as a three-component nodal MultiFab, and use
   an owner mask for unique-node dot products with MPI global reductions. Do
   not copy full Krylov vectors through the host each iteration. This is a
   deliberate departure from production: `LinearSolver`
   (`include/LinearSolver.h`) operates on flattened host `double*` vectors
   through a `MATVEC` callback. Record what the device version needs from the
   interface; changing it is the main solver finding for the production port.
5. Close the timestep: accumulate response, solve the theta-stage E, push and
   redistribute particles, update center B and node B. Follow the audited
   FLEKS sequencing and stage conventions.
6. Port the divE correction: for `nDivECorrection` iterations, sum charge to
   cell centers, solve for the correction potential with the same device
   GMRES (`divESolver`, loose tolerance and small iteration limit), and
   correct particle positions (`divE_correct_particle_position()`), then
   redistribute. Deposit diagnostic moments after the correction. Until this
   part lands, run with the correction disabled and state that in every
   reported result.

Allocate matrix and Krylov scratch storage once. Begin without preconditioning
on small cases, expose restart length, iteration limit, and residual tolerance,
and report nonconvergence as a failure. Verify the true residual independently
of the solver's recurrence before interpreting the physical result.

Acceptance:

- Mass-matrix contributions and operator outputs match a serial reference.
- Manufactured systems converge to the requested relative/absolute residual.
- A small neutral plasma-oscillation problem with omega_p dt well below 1 has
  the expected frequency and converges under timestep refinement. FLEKS uses
  `theta = 0.51` (`FieldSolver` in `include/Pic.h`), which damps and shifts
  the oscillation as omega_p dt grows; state theta and omega_p dt for every
  run, and add a `theta = 0.5` run as an energy-conservation check.
- Record kinetic/field energy, total charge, Gauss residual, and divergence of
  B; define acceptance bounds from the discretization, solver tolerances,
  theta dissipation, and whether the divE correction is enabled.
- A short CPU/GPU run agrees within documented floating-point tolerances.

Milestone B is a self-consistent single-level semi-implicit PIC demo.

## Stage 6: profile and compare porting choices

Profile initialization, field kernels, gather, push, deposition, matrix
assembly/application, the full solve, ghost exchange, and redistribution.
Exclude startup and reference checks from steady-state measurements. Warm up
kernels, synchronize at timing boundaries, and report the maximum elapsed time
over ranks. Avoid per-kernel synchronization in the production timestep.

Report total timestep time alongside particles/s, cells/s, solver iterations,
memory use, and an explicitly stated definition of any effective bandwidth.
Benchmark uniform and clustered particle distributions and several box sizes.
Begin with host-staged MPI as in the current demo; test GPU-aware MPI only on a
runtime that supports device buffers.

After recording the AoS/atomic baseline, change one choice at a time: AoS vs
SoA, particle sorting frequency, atomic vs buffered deposit, mass-matrix storage,
and separate vs fused kernels. Mass-matrix accumulation is the expected
hotspot, so evaluate atomic vs buffered deposit there first. Retain a change only when checks pass and the
measurement identifies its effect on complete timestep cost.

## Verification and build workflow

Planned build and smoke commands after stage 0 and CUDA AMReX configuration:

```bash
make GPUTEST -j8
make -C src_gpu GPUTEST
bin/main_gpu.exe n_cell=16 max_grid_size=8 n_iter=2
mpirun -np 2 bin/main_gpu.exe n_cell=16 max_grid_size=8 n_iter=2
python3 tools/format_all.py src_gpu
```

Keep the executable invocation above compatible while moving the sources. Add
new mode examples when those modes actually exist. Run CPU builds with a
separately configured CPU AMReX installation and clean demo objects when
switching compilers/backends. CUDA compilation alone does not verify GPU
execution. If no GPU is visible, report build verification and the missing
runtime checks separately.

Use tiny deterministic reference cases for correctness and larger cases only
for profiling. Independent serial references can repeat a coefficient or sign
error from the same reading of FLEKS, so also compare stage 2 operators and
stage 5 single-step fields against production FLEKS on CPU for one matching
tiny configuration. Exchange data through output files so the demo build stays
independent of `libFLEKS.a`. GPU memory/race checking is appropriate after adding particle
scatter and matrix assembly. Future AMR work starts only after the single-level
multi-box and multi-rank checks pass.
