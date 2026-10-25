# FLEKS GPU PIC demo

A standalone, single-level 3D AMReX PIC application for evaluating the FLEKS GPU
port. It includes nodal E/interpolated B, cell-centered evolved B, production
AoS particle attributes, CIC gather/scatter, Boris motion, periodic MPI
migration, the 243-component particle response, device GMRES, and particle
position correction for div(E). It links independently of `libFLEKS.a`.

See [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) for scope and acceptance,
[PROGRESS.md](PROGRESS.md) for the implementation record,
[FIELD_NOTES.md](FIELD_NOTES.md) for equations and solver findings, and
[PARTICLE_NOTES.md](PARTICLE_NOTES.md) for particle conventions and benchmarks.
Measured initial costs and their limits are recorded in
[BENCHMARKS.md](BENCHMARKS.md).
The current performance review and prioritized next steps are in
[PERFORMANCE_PLAN.md](PERFORMANCE_PLAN.md).

## Build and run

Use the repository's configured compiler and AMReX installation:

```bash
make GPUTEST -j4
mpirun -np 1 bin/main_gpu.exe mode=fields n_cell=8 max_grid_size=4 n_iter=2
mpirun -np 2 bin/main_gpu.exe mode=particles n_cell=8 max_grid_size=4
mpirun -np 1 bin/main_gpu.exe mode=implicit n_cell=8 n_step=100
```

For a new CUDA configuration, the repository recipe is
`./Config.pl -install -compiler=nvfortran,nvcc`. Do not reconfigure an existing
working build unless switching toolchains intentionally. For compilation on a
host where `nvidia-smi` is unavailable, supply the intended architecture, for
example `make GPUTEST CUDA_ARCH=89` for the RTX 4090 used in validation.

`make -C src_gpu` is equivalent to building `GPUTEST` directly. Sources and
objects stay in this directory; compiler-generated `.d` files track headers.
Root clean targets include this directory. Clean demo objects when switching
AMReX installations, compilers, or CPU/GPU backends.

This host's relocated HPC-X installation needs the following environment for
MPI launch (the paths are machine-specific):

```bash
export OPAL_PREFIX=/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/comm_libs/12.9/hpcx/hpcx-2.22.1/ompi
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader
```

The demo defaults to host-staged MPI buffers. Enable
`amrex.use_gpu_aware_mpi=1` only on a transport supporting device pointers.
Its AMReX device arena starts at 64 MiB and grows as needed; override
`amrex.the_arena_init_size` to study allocation behavior.

## Modes and parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `mode` | `fields` | `fields`, `particles`, `deposit`, `implicit`, `benchmark`, or `all` |
| `n_cell` | 16 | Cells per dimension in a periodic unit cube |
| `max_grid_size` | 8 | Largest box edge |
| `n_iter` | 10 | Repetitions for microbenchmarks |
| `particles_per_cell` | 8 | Per species; cube counts use a quiet subcell lattice |
| `seed` | 1 | Deterministic within-cell sample permutation |
| `n_step` | 20 | Coupled timesteps |
| `dt` | 0.01 | Timestep in normalized units, c=1 |
| `theta` | 0.51 | Electric-field time centering, 0.5 through 1 |
| `mass_ratio` | 25 | Ion mass; electron mass is 1 |
| `initial_condition` | `uniform` | Uniform electron drift or `langmuir` sinusoidal velocity |
| `perturbation` | 0.001 | Initial electron velocity amplitude |
| `b0` | 0 | Uniform initial Bz |
| `coef_diff` | 0.1 | Production electric operator's divergence mixture |
| `gmres_restart` | 20 | Restart length for persistent device Krylov vectors |
| `gmres_max_iter` | 200 | Electric solve iteration limit |
| `solver_tolerance` | 1e-9 | Electric relative RHS-norm tolerance |
| `correct_div_e` | 1 | Enable particle-position div(E) correction |
| `n_div_e_correction` | 3 | Correction passes per timestep |
| `div_e_tolerance` | 0.01 | Scalar correction relative tolerance |
| `div_e_max_iter` | 20 | Scalar solve iteration limit |
| `diagnostic_interval` | 1 | DIAG output interval; always include initial/final states |
| `fused_push` | 0 | Fuse field gather with the particle push |
| `response` | atomic | Mass-matrix assembly: `atomic` (per particle), `cell` (GPU block reduction per cell), or `buffered` |
| `buffered_response` | 0 | Alias for `response=buffered` |
| `deposit` | atomic | div(E) charge/matrix and moment deposits: `atomic` (per particle) or `cell` (GPU block reduction per cell); the charge-only deposit always stays per particle |
| `sort_interval` | 0 | Additional periodic particle sorting; zero disables |
| `benchmark_timing` | stages | `stages` adds synchronized phase timers; `throughput` times one complete measured window |
| `warmup_steps` | 1 | Coupled steps excluded before benchmark measurements |
| `microbenchmarks` | 1 | Run layout/deposit/ghost microbenchmarks after benchmark evolution |
| `pair_symmetry` | 0 | Experimental cell-kernel reuse of unordered node-pair arithmetic |
| `gmres_batched` | 1 | Twice classical Gram–Schmidt with one host wait per Arnoldi step; 0 selects the reference twice modified Gram–Schmidt |
| `solver_sections` | 0 | In `benchmark_timing=stages`, time solver sections (operator, orthogonalization, updates, residual); adds synchronization |
| `div_e_method` | gmres | div(E) correction solve: `gmres` (production), `jacobi_gmres`, or `pcg`; the last two change the applied correction |
| `cell_precision` | double | `single` stages and sums the response matrix per cell in float (and the div(E) matrix with `center_kernel=pairs`; GPU cell kernels only); charge, jHat, and moments stay double. Faster on GPUs with slow FP64, but theta=0.5 energy drift rises from roundoff to about 1e-8 |
| `center_kernel` | moments | div(E) matrix deposit with `deposit=cell`: `moments` (separable per-cell sums, double precision) or `pairs` (node-pair kernel) |
| `fold_charge` | 0 | Experimental cell correction deposit staging of charge-weighted gradients |
| `cell_chunk_size` | 0 | Experimental maximum particles per cell block; zero is unlimited, otherwise at least 64 |
| `clustered` | 0 | Compress quiet-start positions to test scatter contention |
| `production_ghost` | 3 | Comparison width in the ghost-exchange benchmark |

`fields` preserves the original array checks and adds dual-grid operators,
custom-value communication, unique-node reductions, and manufactured GMRES
systems. `particles` and `deposit` run the same comprehensive positive/negative
species reference checks, including all moments and response matrices. `all`
also advances the coupled plasma. A failed check or solve gives a nonzero exit.

The default uniform electron drift excites a spatially uniform electric plasma
oscillation with omega_p=sqrt(4*pi*(1+1/mass_ratio)). This isolates time
centering from interpolation dispersion and initially satisfies Gauss's law.
`initial_condition=langmuir` exercises nonuniform gather, response, and div(E)
correction; its finite-grid frequency differs from that uniform limit.

Particles start with x at the half stage and v at integer time. Each step
assembles response, solves the theta-stage E, kicks v and drifts with v_new,
advances B, corrects particle positions, and deposits diagnostics. The scalar
charge time-centering parameter is 0.51 independently of field theta, matching
the current production source. Theta=0.51 permits numerical dissipation; use
theta=0.5 for the energy-conservation comparison.

`DIAG` lines report kinetic/field energy, charge, particle count, Gauss/divB
norms, absolute true L2 solver residuals and iterations, and the correction/theta
settings. `electric_x` is the global minimum of Ex: it is the signed uniform
field for the default test, not a mean for nonuniform runs. Gauss is evaluated
against production's time-centered old/new cell charge before saving new charge.

## Validation

```bash
make GPUTEST GPUTEST_REFERENCE -j4
python3 src_gpu/validate_demo.py --mpi --physics \
  --reference bin/gpu_reference.exe --output-dir /tmp/fleks-gpu-validation
```

The standard-library-only runner keeps logs, diagnostics, paired reference CSVs,
and `results.json` in the requested output directory. It checks single/two-rank
runs, decomposition invariance, plasma frequency, timestep refinement, and
energy behavior. `gpu_reference.exe` executes current production operator
expressions on CPU pinned memory and compares them with the demo GPU kernels.
Its generated header is stamped with the production source hash. The reference
fixture described in [VALIDATION.md](VALIDATION.md) also compares a matched
one-step plasma result with the existing production CPU FLEKS executable.

## Profiling and porting decisions

```bash
mpirun -np 1 bin/main_gpu.exe mode=benchmark n_cell=16 max_grid_size=8 \
  n_step=10 n_iter=20 correct_div_e=0
mpirun -np 1 bin/main_gpu.exe mode=benchmark n_cell=16 max_grid_size=8 \
  n_step=10 n_iter=20 correct_div_e=0 fused_push=1 buffered_response=1 clustered=1
```

Benchmark mode excludes `warmup_steps` coupled steps (default one). Stage timings use GPU
synchronization at measured boundaries and the maximum rank elapsed time. They
include response assembly/exchange, electric solve, gather/push/migration,
magnetic advance, charge/moment deposition, and correction when enabled. Normal
implicit mode does not insert those per-stage timing synchronizations.

For complete evolution throughput without phase-timer barriers, use
`benchmark_timing=throughput microbenchmarks=0`. This synchronizes the GPU
only at the endpoints of the measured window and reduces elapsed time across
ranks once at the end. Algorithm-required synchronization remains in the
solvers and diagnostics. Window time includes diagnostics, reported separately
as `TIMING stage=diagnostics_included`; this measures wall time in diagnostic
calls, including waits for pending evolution work. It is a subset of timestep time and
must not be added to it. Stage-mode timestep time excludes diagnostics. Use
matching `diagnostic_interval` in comparisons. Both modes divide by completed
measured steps; aborted runs are marked and do not report throughput. The
initial diagnostic and warmup are outside both measurement windows.

For the measured 64^3, 100-particle-per-cell RTX 4090 workload, the best
validated settings are `response=cell deposit=cell pair_symmetry=1 fused_push=1`.
Five alternating trials reduced median timestep time by 14.7% and reduced
peak device-arena usage by about 2,406 AMReX-reported MB. See
[BENCHMARKS.md](BENCHMARKS.md) for spread, solver work, and workload limits.

The optional GPU cell options can be evaluated independently:
`pair_symmetry=1` computes 36 unordered spatial pairs and writes both
destinations, `fold_charge=1` stages charge-weighted correction gradients,
and `cell_chunk_size=256` splits populated bins into bounded particle chunks
for response, correction-matrix, and moment deposits. Chunk scheduling adds
a prefix scan and a host synchronization per tile; zero preserves one block
per cell. Folded gradients use an additional 12 KiB of shared memory per block.
These options preserve atomic output accumulation and CPU fallback paths.
`gmres_batched=1` (the default) uses twice classical Gram–Schmidt. On
CUDA/HIP both projection passes, their vector updates, and the norm run on
the device without host waits; the host waits once per Arnoldi step for the
Hessenberg column (once per pass with several MPI ranks). The basis is stored
once, in a packed array that the per-vector MultiFabs alias. CPU builds use
local dot products. `gmres_batched=0` keeps twice modified Gram–Schmidt as the
reference. Benchmark mode prints one `SOLVER` line per solver with
iterations, host waits per iteration, and, with `solver_sections=1`, section
times. `DIVE_SOLVER_OPTIMIZATION_PLAN.md` records the measurements.

For repeated comparisons, with logs, diagnostics, source/build provenance,
CPU affinity/load, arena peaks, and median/spread in `results.json`, run:

```bash
python3 src_gpu/measure_performance.py --output-dir /tmp/fleks-perf-trials \
  --cases baseline pairs folded chunks batched fused combined --repeats 5
```

The runner alternates case order, disables microbenchmarks, checks physical
diagnostics against the baseline, and uses the documented 64^3 reference
problem by default. Use an empty output directory. Set `CUDA_VISIBLE_DEVICES`
and CPU affinity for an idle device/core before launch; when pinning the
parent, `--launcher 'mpirun --bind-to none'` preserves inherited affinity.
Use `--set n_cell=32 max_grid_size=16 particles_per_cell=27` for common
problem overrides, `--ranks 2` for MPI, or `--timing stages` for attribution.
AMReX profiling ranges distinguish bin construction, migration, response
communication operations, and phase submissions. A GPU timeline is needed
to distinguish asynchronous kernel execution from host submission time.

Microbenchmarks compare AoS with explicit component-SoA scratch, atomic with
cell-buffered response assembly, custom with component matrix storage, and
ghost widths. The SoA test freezes particle positions and excludes conversion
from kernel time while reporting setup separately; it does not measure a
PureSoA container or migration. The buffered method includes sorting/range
construction and temporary allocations. Fused push and sorting choices must
also be compared on complete timestep timings using matching initial states.

The custom mass matrix uses 243 reals per node plus ghosts. Startup prints a
conservative per-rank allocation estimate and rejects configurations above
70% of device capacity. These small cases demonstrate correctness and expose
costs; they do not establish production-scale speedups. AMR, wall boundaries,
coupling, and trajectory-based explicit current deposition remain outside this
single-level implicit demo.
