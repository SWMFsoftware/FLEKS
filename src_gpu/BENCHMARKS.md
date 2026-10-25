# GPU measurements — 2026-10-24 to 2026-10-25

Sections are chronological; later sections supersede earlier numbers. The
final state is summarized in `GPU_DEVELOPMENT_SUMMARY.md`.

These are short exploratory measurements on NVIDIA GeForce RTX 4090 GPUs
(24 GiB each), AMReX `45a1fc5e1e5a`, nvcc sm89, double precision, and host-staged
MPI. They are correctness-preserving comparisons, not production scaling
results or statistically established speedups. No other demo or sanitizer
process was running during this series.

All eight cases passed their conservation and storage/layout comparisons.
Common arguments were:

```text
mode=benchmark n_cell=16 max_grid_size=8 particles_per_cell=8
n_step=8 n_iter=10 correct_div_e=0 diagnostic_interval=8
```

Each case has 65,536 particles, two species, and a uniform cold electron-drift
initial condition, except that the clustered cases compress positions into
the central region. Coupled timings exclude the first step and diagnostics,
include charge-history bookkeeping, and synchronize at measured phase
boundaries. Two-rank times use the maximum rank elapsed time.

| Case | Changed arguments | Mean timestep (ms) |
| --- | --- | ---: |
| Uniform baseline | none | 6.072 |
| Uniform fused gather/push | `fused_push=1` | 5.706 |
| Uniform sorting each step | `sort_interval=1` | 6.520 |
| Uniform buffered response | `buffered_response=1` | 14.242 |
| Clustered atomic | `clustered=1` | 7.653 |
| Clustered buffered | `clustered=1 buffered_response=1` | 388.876 |
| Smaller boxes | `max_grid_size=4` | 98.612 |
| Two GPUs | `mpirun -np 2` | 68.085 |

Baseline phase means were 0.768 ms for response assembly, 3.520 ms for matrix
exchange, 0.958 ms for the electric solve, and 0.364 ms for separate particle
gather plus push/migration. The tiny problem is dominated by communication and
launch/control overhead, so the two-GPU case provides a communication test,
not evidence of useful strong scaling.

For the baseline state's microbenchmarks:

| Operation | Mean time (ms) |
| --- | ---: |
| Atomic response assembly including clearing | 0.764 |
| Cell-buffered response including sorting/setup/clearing | 10.017 |
| Custom 243-real matrix application | 0.0570 |
| Component MultiFab matrix application | 0.0434 |
| Ghost exchange, width 1 | 0.00749 |
| Ghost exchange, width 3 | 0.00770 |
| AoS frozen-position gather/push, electron species | 0.0257 |
| Component-SoA scratch gather/push, electron species | 0.0279 |

The matrix layouts produced matching results; the particle layout test's
velocity difference was zero. SoA setup/transpose was measured separately
(0.0736 ms for that species) and excluded from kernel time. This is a scratch
layout test, not an AMReX PureSoA container/migration benchmark.

The cell-buffered implementation serially accumulates all particles in a cell
before flushing 64 tensor blocks. It reduces global atomics but introduces
large thread-local storage and poor parallelism in densely occupied cells.
The clustered result is a concrete reason to keep the atomic path as default;
this buffering strategy is an experimental comparison, not an optimization to
transfer directly into FLEKS. Each thread holds 600 doubles (4.8 KB) of
private buffers, which spill from registers to local memory, and each call
allocates range arrays and synchronizes the stream once per box. The result
therefore measures this implementation, not buffering in general; it does not
rule out block- or warp-level reduction in shared memory. Similarly, the small fused-kernel difference
needs longer repeated trials before attributing a speedup to fusion.

Logs and structured results from this run are in
`/tmp/fleks-gpu-benchmarks/`; rerun the documented commands to reproduce them.
The most useful next profiling target is custom matrix communication at larger
per-box workloads, followed by cooperative buffering with bounded work per
thread rather than a serial per-cell accumulation.

## Cell-reduced mass-matrix assembly — 2026-10-24

`response=cell` (`DemoParticles::accumulate_response_cell`) bins particles by
cell with an AMReX `DenseBins` counting sort (an index permutation; particles
are not moved) and launches one 64-thread block per cell. Thread `t` owns the
node pair `(t/8, t%8)` and accumulates its 3x3 tensor in registers over every
particle in the cell; particle weights, `alpha` and `jHat` terms are staged in
shared memory in chunks of 64. Each cell then issues 576 atomics, with
distinct addresses inside a warp, instead of 576 atomics per particle with
all lanes of a warp hitting the same addresses. The CPU build falls back to
the per-particle path.

One RTX 4090, 64^3 cells, 100 particles per cell per species (52.4M
particles), 8 boxes of 32^3, Langmuir wave, `b0=0.5`, `dt=0.1`, div(E)
correction on, `solver_tolerance=1e-6`; steps 2-6 timed:

| Assembly | Seconds per assembly | Timestep (s) | Particles/s |
| --- | ---: | ---: | ---: |
| `atomic` | 0.536 | 1.064 | 49.3M |
| `cell` | 0.073 | 0.603 | 87.0M |
| `buffered` | 0.672 | (not run) | |
| 32 CPU cores, `atomic` | 0.443 | 2.163 | 24.2M |

All three GPU assemblies agree entry by entry within 7e-17, and the coupled
diagnostics after six steps match the atomic run to 2e-15. Clustered particles
(`clustered=1`, 10 particles per cell, about 1,250 per occupied cell) remain
favorable: 0.050 s for `cell` versus 0.112 s for `atomic`.

Nsight Systems attributes about 0.5 ms per box and species to the counting
sort and about 5.9 ms to the block kernel. Hardware counters are restricted to
administrators on this host (`RmProfilingAdminOnly=1`), so Nsight Compute
could not run. A simple throughput model (FP64 on this consumer GPU is 1/64 of
FP32) estimates about 5 ms: roughly 3.6 ms for the ten FP64 instructions per
particle and node pair, plus 1.5 ms of per-particle staging. The kernel is
therefore likely FP64-bound. Exploiting the symmetry of the pair weights
(36 unique pairs instead of 64) could remove about 40% of the pair work; at
this point, however, the divE charge/matrix deposit (0.22 s) and the moment
deposit (0.12 s) cost more per step and would benefit from the same per-cell
reduction.

## Cell-reduced particle deposits — 2026-10-25

`deposit=cell` applies the per-cell block reduction to the div(E)
charge/correction-matrix deposit (`deposit_center_charge_cell`, bins over the
half-cell-shifted base cells including one low ghost layer) and to the
ten-moment deposit (`deposit_cell`). Same reference case as above, with
`response=cell`, the process pinned to one idle core:

| Stage | `deposit=atomic` (s) | `deposit=cell` (s) | Speedup |
| --- | ---: | ---: | ---: |
| div(E) charge + matrix, 3 passes | 0.219 | 0.111 | 2.0x |
| Moment deposit | 0.122 | 0.016 | 7.5x |
| Final charge deposit (per particle in both) | 0.013 | 0.013 | 1.0x |
| **Timestep** | **0.599** | **0.387** | **1.55x** |

Throughput is 135.5M particles per second, 5.6x the 32-core CPU run. The
div(E) kernel is FP64-bound after the atomic reduction, and clustered
particles do not benefit; `DEPOSIT_OPTIMIZATION_PLAN.md` ("Outcome") gives
the profile and the remaining options.

A competing 16-rank CPU job started during the first measurements and
doubled the timestep (1.28 s, with the electric solve at 0.195 s instead of
0.047 s): the host drives every Krylov iteration, so GPU timings on this
shared machine are sensitive to CPU load. The numbers above are from a quiet
period and were reproduced twice.

## Spatial pair reuse and throughput timing — 2026-10-25

Implemented the next performance plan with opt-in pair reuse, charge folding,
bounded cell chunks, and batched CGS2. The timing harness now supports a full
throughput window with configurable warmup and no per-phase timer barriers.
`measure_performance.py` alternates trial order and records logs, diagnostics,
configuration, source/executable hashes, hardware, affinity, CPU load, and
median/spread in JSON. It rejects failed solves, diagnostic disagreement,
and additional correction-limit warnings.

On idle GPU 1 (RTX 4090), CPU affinity 8 with MPI binding disabled, the same
64^3/100-ppc reference configuration was measured for six steps, excluding
one warmup step. `diagnostic_interval=6`, `microbenchmarks=0`, and
`benchmark_timing=throughput`; diagnostics are included in this window.
Each case ran five times in alternating order. All cases used
`response=cell deposit=cell`. GPU 0 had another user's workload, and host
load was recorded; these are shared-host results, not isolated-machine results.

| Additional flags | Median timestep (s) | Min–max (s) | Speedup |
| --- | ---: | ---: | ---: |
| None | 0.41661 | 0.41645–0.41801 | 1.000x |
| `pair_symmetry=1` | 0.36040 | 0.36000–0.36111 | 1.156x |
| `pair_symmetry=1 fold_charge=1` | 0.37508 | 0.37455–0.37531 | 1.111x |
| `fused_push=1` | 0.41196 | 0.41191–0.41312 | 1.011x |

Pair reuse reduces timestep time by **13.5%**, raising throughput from about
125.8M to 145.5M particles/s. Physical diagnostics agree with the baseline.
Every trial had the same measured solver work: 70 electric iterations,
287 correction iterations, and six best-effort correction-limit warnings.
The warning count is a property of this reference run and did not increase
with the candidates. AMReX reported peak device-arena usage of about 8,551 MB
for the baseline; this is arena usage, not total process/device memory.

The first pair implementation assigned 36 active lanes out of 64 and showed
no improvement (0.42257 versus 0.42064 s in a separate ten-step series).
The final implementation keeps the second logical warp useful: its 32 lanes
share particle accumulation for the last four unordered pairs, followed by
shuffle reductions. It writes both ordered destinations without transposing
the magnetized response tensor. Charge folding adds 12 KiB of shared memory
and staging work; it is slower than pair reuse alone on this GPU.

The ten-step exploratory series also found that chunk size 256 on the uniform
reference added overhead, and an MPI-only CGS2 prototype was slower. The
CGS2 implementation was subsequently revised to batch device-local projections
as well. The discarded lane mapping and MPI-only prototype are not retained
as additional runtime options. Do not attribute their measurements to the
final implementations.

Raw logs and structured results for the retained pair mapping are in
`/tmp/fleks-performance-measure-revised/`; the exploratory series is in
`/tmp/fleks-performance-measure-v2/`. The earlier 0.387 s stage-timed result
used a different measurement series and excludes diagnostics; compare
candidates within the alternating series above rather than treating those
two baselines as equivalent.

### Final combined settings and clustered cells

After batching GPU-local CGS2 projections as well, a fresh alternating
five-trial series with the same six-step protocol gave:

| Additional flags | Median timestep (s) | Min–max (s) | Speedup |
| --- | ---: | ---: | ---: |
| None | 0.41599 | 0.41491–0.41631 | 1.000x |
| `pair_symmetry=1 fused_push=1` | **0.35493** | 0.35391–0.35507 | **1.172x** |
| `gmres_batched=1` | 0.42389 | 0.42290–0.42993 | 0.981x |
| Pair reuse + fused push + batched GMRES | 0.36301 | 0.36155–0.37039 | 1.146x |

Recommended flags for this workload are
`response=cell deposit=cell pair_symmetry=1 fused_push=1`, with the other
experimental options disabled. Timestep time fell **14.7%**, and throughput
rose from about 126.0M to **147.7M particles/s**. Solver iteration and warning
counts remained identical to the preceding table. Peak device-arena usage
fell from 8,551 to 6,145 AMReX-reported MB, consistent with eliminating the
six-double gathered-field scratch per particle. CGS2 adds about 180 reported
MB for the packed basis. Batching removes synchronization but its extra
packing, updates, and reductions still cost more here; retain it as an
optional comparison for other workloads/transports.

For `clustered=1 particles_per_cell=10` (5.24M particles), otherwise the same
configuration and protocol, another five alternating trials measured:

| Additional flags | Median timestep (s) | Min–max (s) |
| --- | ---: | ---: |
| None | 0.24425 | 0.24404–0.24429 |
| `pair_symmetry=1` | 0.23792 | 0.23787–0.23840 |
| `pair_symmetry=1 cell_chunk_size=256` | 0.23616 | 0.23614–0.23636 |

Chunking added a small gain for clustered cells, unlike the uniform case.
All clustered trials had 227 electric iterations, 300 correction iterations,
and 15 correction-limit warnings; diagnostics matched. These are best-effort
correction runs, not evidence that the scalar solve converged. Keep chunking
explicit rather than selecting it automatically from this limited sample.

Logs and JSON: `/tmp/fleks-performance-measure-final/` and
`/tmp/fleks-performance-clustered/`. Single- and two-rank GPU checks use the
one selected device; they validate MPI correctness, not two-GPU scaling.

To repeat the final comparison after selecting an idle GPU and configuring
the MPI environment from README.md:

```bash
python3 src_gpu/measure_performance.py --output-dir /tmp/fleks-perf-final \
  --cases baseline pairs-fused batched combined --repeats 5 \
  --steps 6 --warmup-steps 1
```

Add `--set clustered=1 particles_per_cell=10 --cases baseline pairs chunks`
with a new output directory for the clustered comparison. Preserve the same
GPU assignment, CPU affinity, and MPI binding policy across candidates.

### Profiling and validation

An Nsight Systems CUDA timeline for the recommended settings is saved as
`/tmp/fleks-performance-profile.nsys-rep`, with kernel/API summaries in
`/tmp/fleks-performance-profile-stats.csv`. The profiled, stage-synchronized
run took 0.392 s/step; instrumentation makes it unsuitable as a replacement
for the unprofiled throughput result. Its largest stages were center deposit
(106 ms), electric solve (64 ms), correction plus migration (57 ms), response
assembly (54 ms), and scalar solves (39 ms). These guide the next review.

Across the entire trace, including initialization and the first step, the
correction-deposit kernel used 22.4% of GPU kernel time, electric response
matvec 15.0%, response assembly 12.3%, position correction 8.1%, and particle
destination partitioning 7.0%. CUDA stream waits were 85.1% of CUDA API wall
time; waits include GPU execution and must not be added to kernel time.
Similarly, inclusive host bin-building scopes can absorb preceding GPU work;
use the CUDA timeline rather than interpreting their wall time as sort cost.
Neither this trace nor the throughput tests establishes an FP64 hardware
counter diagnosis or useful two-GPU scaling. GPU 0 remained occupied by
another user's workload, so cross-device scaling measurements were deferred.

The final CUDA build passed all 45 cases in `validate_demo.py --mpi --physics
--reference ...`, including magnetized pair reuse, folded/split deposits,
fused push, difficult batched solves, timing-mode equivalence, and production
operator/physics checks. Results: `/tmp/fleks-performance-validation-final/`.
The CPU build passed the 45-case suite and 12 timing edge cases; refreshed
final-solver checks also passed `mode=all`, a two-rank magnetized large-dt
comparison, and the measurement-runner smoke test. CPU artifacts are under
`/tmp/fleks-performance-cpu/`.

Compute Sanitizer reported zero memory errors and zero race hazards for
clustered cooperative kernels and again after the final device-batched
projection change. Logs are `/tmp/fleks-performance-{memcheck,racecheck}.log`
and `/tmp/fleks-performance-final-{memcheck,racecheck}.log`.

## Synchronization-free CGS2 solver — 2026-10-25

`gmres_batched=1` was rewritten to remove host waits and is now the default
(`DIVE_SOLVER_OPTIMIZATION_PLAN.md`, "Outcome"). With the recommended flags
on the reference case, three alternating trials per configuration:

| Layout | Scalar div(E) solves, MGS2 -> CGS2 (ms/step) | Throughput step, MGS2 -> CGS2 (s) |
| --- | ---: | ---: |
| 8 boxes of 32^3 | 33.80 -> 13.79 | 0.3449 -> 0.3203 |
| 64 boxes of 16^3 | 35.58 -> 15.30 | 0.3999 -> 0.3728 |

Host waits per scalar iteration fell from 43.8 to 1.10 with identical
iteration counts and warnings. The scalar solve is now 2.0x faster than 64
CPU cores. The earlier `gmres_batched=1` measurements above (0.981x) refer
to the previous implementation.

## Single-precision cell matrix sums — 2026-10-25

`cell_precision=single` stages the pair inputs of the two matrix kernels in
float and sums each cell's response tensor and correction-matrix entries in
float; each per-cell sum is still added into the double-precision global
matrix. Charge density, jHat, and moments cancel between species and stay in
double precision. The RTX 4090 runs FP64 at 1/64 of its FP32 rate, and the
pair loops were FP64-bound.

Reference case, recommended flags, 8 boxes of 32^3, one MPI rank pinned to an
idle core, three alternating trials per mode:

| Stage (ms/step, stage mode, GPU 0) | Double | Single | Speedup |
| --- | ---: | ---: | ---: |
| Response matrix assembly | 50.56 | 25.24 | 2.00x |
| div(E) charge + matrix, 3 passes | 95.80 | 70.53 | 1.36x |
| All other stages | unchanged | unchanged | 1.00x |
| **Stage-timed step** | **318.5** | **267.9** | **1.19x** |

Throughput mode on GPU 1 (another user's job occupied GPU 0 during the
first throughput series, which was discarded): 0.3321 s -> 0.2794 s per step
(all trials within 0.0004 s), 157.9M -> 187.7M particles per second. GPU 1
is about 3% slower than GPU 0 in absolute terms.

The speedup is smaller than the up-to-5x estimated for these kernels before
measuring: per-particle staging (CIC weights, the nodal-B gather, alpha,
gradients), the double-precision charge and jHat sums, counting sorts, and
the boundary exchanges still run in double precision or are not
arithmetic-bound. Moving the staging arithmetic to float is the next
candidate.

Accuracy, double versus single:

| Check | Double | Single |
| --- | ---: | ---: |
| Reference case: solver iterations and warnings | 70 / 287 / 6 | 70 / 287 / 6 |
| Reference case: kinetic energy after 6 steps, relative difference | | 2.5e-9 |
| Matrix entries versus double, relative to the largest | | about 1.3e-7 |
| theta=0.5 energy drift, 600 steps | 2.2e-15 | 8.6e-9 |
| Plasma-frequency error | 1.09e-4 | 1.09e-4 |
| Nonuniform div(E) correction, final Gauss residual | 9.4e-15 | 7.9e-15 |
| Large-dt maximum corrected Gauss residual | 7.5e-9 | 1.5e-8 |
| Large-dt electric / last-pass div(E) iterations, 20 steps | 462 / 134 | 467 / 140 |

Energy conservation at theta=0.5 depends on the mass matrix being exactly
consistent with the particle push; float per-cell sums break that at the
1e-8 level. The option therefore stays opt-in. On data-center GPUs with
half-rate FP64 the arithmetic advantage largely disappears.

Tensor cores (assessment, not measured): the RTX 4090 has 512
fourth-generation Tensor Cores for FP16, BF16, TF32, FP8, INT8 and INT4, but
none for FP64. Each cell's matrix sum is a small matrix product (about
72 x 100 times 100 x 8 per species), so it could map to tensor-core
multiply-accumulate, but TF32 and lower precisions carry about three decimal
digits or fewer, far below the accuracy that already costs energy
conservation in single precision; TF32 tensor throughput on this card equals
plain FP32; and the per-cell products are too small to amortize data layout
changes. The solver work is bandwidth-bound. FP64 tensor cores on A100/H100
could run the same per-cell sums in full double precision and are the
relevant target for data-center GPUs.

## Separable-moment div(E) deposit — 2026-10-25

`center_kernel=moments` (default) rebuilds the div(E) correction matrix from
27 per-cell sums of non-negative terms (CIC weights factor by dimension),
plus 8 charge sums, instead of 36 node-pair gradient dot products per
particle. See `DIVE_DEPOSIT_OPTIMIZATION_PLAN.md`.

Reference case, stage mode, GPU 0, two undisturbed trials each:

| | Pairs | Moments | Gain |
| --- | ---: | ---: | ---: |
| div(E) deposit, double | 96.1 ms | 38.9 ms | 2.47x |
| Whole step, double | 319.2 ms | 262.0 ms | 1.22x |
| Whole step, `cell_precision=single` | 268.4 ms | 236.6 ms | 1.13x |

Throughput mode: 0.3233 -> 0.2660 s. 64 boxes of 16^3: 371.5 -> 309.7 ms per
step. Solver work identical. The kernel is double precision only, so
`cell_precision=single` now speeds up only the response kernel.
