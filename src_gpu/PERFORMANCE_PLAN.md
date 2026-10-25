# GPU demo performance review and next-step plan

> Historical. This review and its implementation (commit `422512c`) were
> followed by a rewritten synchronization-free GMRES (now the default
> `gmres_batched=1`, which supersedes the CGS2 description below), optional
> single-precision cell sums, and the separable-moment div(E) deposit.
> Development concluded on 2026-10-25; see `GPU_DEVELOPMENT_SUMMARY.md` for
> the final state.

Reviewed 2026-10-25. Scope: the current `src_gpu` coupled PIC timestep,
including the already implemented `response=cell` and `deposit=cell` paths.
The original review synthesized recorded measurements. Implementation and
new validation followed on the same date after host execution exposed the
GPU devices hidden by the sandbox. The sections below retain the review's
baseline, rationale, and acceptance criteria.

## Implementation status

- Phase 0: throughput-window timing, warmup control, optional microbenchmarks,
  completed-step/abort accounting, iteration and correction-warning totals,
  profiling ranges, and an alternating repeated-trial runner are implemented.
- Phases 1–2: optional `pair_symmetry=1` reduces both cell matrix kernels to
  36 pairs. The final four pairs share particle work across eight lane groups
  so the second warp remains useful. Reverse response entries retain the
  original nonsymmetric tensor. Optional `fold_charge=1` is measured separately.
- Phase 3: operator helpers have explicit halo contracts; duplicate fills and
  unused output ghosts were removed. Existing fused push can be tested by the
  runner. Custom response exchange has separate profiling ranges for each
  accumulation, ownership synchronization, and fill operation.
- Phase 4: basis normalization is fused. Optional `gmres_batched=1` uses CGS2
  projection batches and fused basis/solution updates; true residual and
  ownership-mask checks remain mandatory.
- Phase 5: optional `cell_chunk_size` bounds particles per block for response,
  correction matrix, and moments. Prefix-scan scheduling overhead is explicit.
  Migration is profiled separately. Deferred migration, preconditioning,
  automatic selection, and matrix-storage changes remain conditional on
  profiling and correctness evidence, as specified below.

See [BENCHMARKS.md](BENCHMARKS.md) for measured outcomes and recommendations,
and [README.md](README.md) for the flags and reproducible measurement command.

## Baseline and interpretation

[BENCHMARKS.md](BENCHMARKS.md) and the outcome section of
[DEPOSIT_OPTIMIZATION_PLAN.md](DEPOSIT_OPTIMIZATION_PLAN.md) report the
following reference: one RTX 4090, double precision, 64^3 cells, two species
with 100 particles per cell each (52.4 million particles), eight 32^3 boxes,
Langmuir initial condition, `b0=0.5`, `dt=0.1`, three div(E) corrections,
and electric solver tolerance 1e-6.

With `response=cell deposit=cell`, the recorded timestep is **0.387 s**
(135.5 million particles/s), versus 0.599 s with atomic deposits and cell
response. Those deposit optimizations are complete; do not repeat the old
plan as new work.

| Cost | Recorded time | Interpretation |
| --- | ---: | --- |
| Three div(E) charge/matrix deposits | 111 ms | Largest quantified remaining deposit cost; about 29% of 387 ms |
| Response assembly | 73 ms in the preceding reference series | Next arithmetic target; remeasure with the current configuration |
| Ten-moment deposit | 16 ms | Already reduced from 122 ms |
| Final charge deposit | 13 ms | Atomic path beat the tested cell path (16.4 ms) |
| Binning across the step | About 18 ms | Included in deposit/assembly costs, not an additional stage |

The 73 ms response timing and the latest 387 ms total are from different
measurement series. Do not present their difference as a measured breakdown
of all remaining stages. Obtain fresh solve, exchange, push, and migration
times before ranking those stages against each other.

The tiny 16^3 benchmark tells a different story: matrix exchange took
3.52 ms out of 6.07 ms, smaller boxes raised the timestep to 98.6 ms, and
two GPUs took 68.1 ms. These establish overhead sensitivity, not useful GPU
scaling. The recorded FP64 bottleneck diagnosis for cell kernels is based
on a throughput model and Nsight Systems timing; restricted hardware
counters prevented direct confirmation.

## Findings at the original review

1. **Pair arithmetic is duplicated in the two largest cell kernels.**
   `DemoParticles::deposit_center_charge_cell` assigns all 64 ordered node
   pairs to 64 threads. Every thread loops over every particle in the cell,
   recomputing a symmetric gradient dot product for the reverse pair.
   `accumulate_response_cell` similarly computes 64 weight products and
   nine tensor accumulations per pair. Only 36 unordered pairs are unique.
   Reducing arithmetic is promising after the existing atomic reduction.

2. **Dense cells serialize work inside too few blocks.** Both kernels launch
   one block per bin and traverse its particles in chunks. The recorded
   clustered correction deposit is slightly slower than atomic deposition
   (33 versus 30 ms). Empty-bin launches and long occupied-bin loops can
   coexist. The same strategy need not win at low particles per cell either.

3. **GMRES has a growing synchronization cost.** In
   `DemoLinearSolver::solve`, twice modified Gram–Schmidt performs
   `2*(j+1)` host-returning dot products and separate SAXPY operations at
   Arnoldi column `j`, plus a norm. A full restart of length `m` therefore
   has `m*(m+1)` projection reductions (420 at `m=20`), excluding norms and
   true-residual checks. Device vectors do not make the iteration
   device-controlled. Existing measurements also show strong sensitivity
   to competing CPU work.

4. **Field operator helpers repeat halo work.** In
   `DemoFields::apply_e_operator`, `curl_node_to_center` fills `curlCenter`,
   then `curl_center_to_node` fills it again. The input is synchronized and
   filled by both the first curl helper and `divergence` without an
   intervening input update. Some intermediate output ghost fills also
   deserve a consumer-by-consumer audit. `finish_response` separately sums,
   synchronizes, and fills a 243-real custom matrix. Preserve the distinct
   accumulation and ownership semantics when optimizing these operations.

5. **Migration and gather traffic are material candidates, not yet ranked.**
   `push` and every `correct_positions` call unconditionally call
   `Redistribute`: eight calls per reference timestep across two species.
   Separate gather writes six doubles per particle that push reads back;
   at 52.4M particles that is about 5.0 GB of logical scratch traffic per
   step, excluding allocation and other particle traffic. `fused_push=1`
   already removes it, but its large-case benefit has not been established.

6. **The timing harness changes execution.** `DemoPic::measure` adds device
   synchronization and an MPI maximum reduction at every stage. The outer
   timestep includes this instrumentation overhead while excluding
   diagnostics. It is useful for attribution but does not measure the
   normal execution path. Summed per-stage rank maxima need not equal the
   maximum rank's total step. Benchmark mode also runs expensive layout
   microbenchmarks after evolution, so process wall time is not timestep time.

## Prioritized implementation plan

### 0. Establish a repeatable current baseline

Before changing kernels, add a coupled-throughput timing option with no
per-stage barriers or reductions: synchronize at the endpoints of a measured
step window and reduce the elapsed time across ranks once. Keep the existing
stage mode for attribution. Account for diagnostics explicitly, provide a
way to skip post-run microbenchmarks, and divide by completed measured steps
if a solve fails. Record all solver iterations and correction warnings.

Start by reproducing the documented stage baseline:

```bash
mpirun -np 1 bin/main_gpu.exe mode=benchmark n_cell=64 max_grid_size=32 \
  particles_per_cell=100 n_step=6 n_iter=3 initial_condition=langmuir \
  perturbation=0.01 b0=0.5 dt=0.1 correct_div_e=1 \
  solver_tolerance=1e-6 response=cell deposit=cell diagnostic_interval=6
```

Use the MPI environment in [README.md](README.md). Reserve an idle GPU and
CPU core, record affinity, CPU load, clocks, compiler/build flags, AMReX and
source revisions, and peak device memory. Collect at least five independent
baseline/candidate runs in alternating order; report median and spread.
Extend warmup and measurement windows until results stabilize. Keep seed,
physics, tolerances, diagnostic cadence, and decomposition identical.

Capture an Nsight Systems timeline with ranges for bins, cell kernels,
halo operations, orthogonalization, and migration. Use hardware counters
only if available; otherwise label occupancy, bandwidth, and FP64 diagnoses
as hypotheses. First gate: reproduce the baseline or explain the difference
before assigning a speedup to a code change.

### 1. Reduce div(E) matrix arithmetic

Implement an optional 36-pair kernel in `deposit_center_charge_cell`.
Compute each unordered pair once and scatter off-diagonal results to both
ordered destinations; write diagonal pairs once. Retain charge accumulation
exactly once per node. Test particle-split thread groups and block sizes so
that reducing the pair count does not simply leave lanes idle.

Separately test factoring `abs(q)/V` into one staged gradient operand.
Account for added staging arithmetic and shared-memory footprint; do not
assume removing an inner-loop multiply is a net gain. Preserve signed
charge for rho and absolute charge for the matrix.

There are 43.75% fewer pair evaluations, but staging, binning, and output
atomics remain. Treat **20–30% lower total center-deposit time** as an
experimental target, not a prediction. On the recorded baseline that would
save about 22–33 ms per timestep. Accept only with a repeatable total-step
gain and unchanged validation results.

### 2. Apply pair reuse to response assembly

In `accumulate_response_cell`, reuse the nine accumulated tensor components
for each reverse spatial pair. The weights satisfy `w_a*w_c = w_c*w_a`;
the magnetized alpha tensor is generally nonsymmetric. **Copy the same
tensor to the reverse spatial destination; do not transpose alpha or assume
the full matrix is symmetric.** Keep jHat accumulation complete and unique.

Benchmark register pressure, staging, and total assembly time. Retain the
full 243-component output representation initially, isolating arithmetic
changes from communication/layout changes. A 20–30% assembly improvement
would save roughly 15–22 ms if the 73 ms baseline reproduces. Combined with
phase 1, this suggests an exploratory 0.33–0.35 s timestep target; it is not
a measured result or a guaranteed outcome.

### 3. Remove avoidable communication and test existing fusion

Give operator helpers explicit input-halo and output-halo contracts. Remove
the duplicate `curlCenter` fill and repeated input preparation within one
operator application first. Fill an intermediate's ghosts only when its next
consumer reads them. Validate box interfaces, periodic ownership, and
production operator equivalence after each change.

Split response exchange timings into SumBoundary, OverrideSync, and
FillBoundary. Compare custom and component matrix storage end to end,
including deposit writes, packing if any, memory footprint, exchange, and
matvec. The old component-layout microbenchmark alone does not justify a
storage rewrite. Sweep box sizes 16 and 32 at fixed 64^3 domain; try 64 only
if the allocation estimate fits. Compare one and two GPUs using meaningful
per-rank work, reporting strong and weak scaling separately. Enable GPU-aware
MPI only with a verified supporting transport.

Run `fused_push=0/1` on the large case now, as a low-effort experiment. Retain
fusion as a recommended option only if timestep and memory results support it.

### 4. Reduce Krylov overhead if the new profile warrants it

First fuse vector copy/scale and multi-vector updates where possible without
changing the algorithm. Then prototype batched projections using a numerically
appropriate two-pass classical Gram–Schmidt variant and fused basis updates.
The current modified Gram–Schmidt projections depend on preceding updates;
they cannot simply be issued concurrently with identical semantics.

Retain owner-mask weighting, true-residual checks, and the current solver as
a comparison path. Measure reduction counts, iterations, time per iteration,
and complete solve time for both electric and div(E) solves. Reject apparent
speedups caused by looser convergence or extra correction-limit failures.
Consider preconditioning only after measuring iteration counts; a nodal
block-diagonal electric preconditioner and a scalar correction preconditioner
need separate mathematical and convergence validation.

### 5. Handle clustered occupancy and migration

Use bin occupancy statistics to prototype bounded particle chunks per cell,
with multiple blocks for oversized bins. Reduce each chunk locally and
atomically accumulate its partial outputs. Tune the threshold against extra
atomics and scheduling overhead, including low-occupancy uniform cases.
Keep explicit atomic/cell selection for reproducibility before considering
an adaptive policy. Response and correction deposits can require different
thresholds.

Profile `Redistribute` separately from push/correction kernels, including
cross-box particle counts. Investigate a proven local migration path or a
device-computed no-crossing fast path only if migration is substantial.
Do not defer redistribution across correction passes without establishing
that every particle remains associated with valid field and deposit storage.
Include the cost of any crossing detection and collective decision.

## Work to defer

- Bin reuse is currently below 5% of timestep time. Revisit only after larger
  costs fall; cached permutations must be invalidated on position changes,
  sorting, migration, initialization, and particle replacement. Center-shape
  and nodal bins are different and cannot be interchanged.
- Keep the atomic final-charge path and ten-moment deposit every step.
  Removing moments would evade production's timestep-update workload.
- Do not promote the serial buffered response path: its large private arrays
  and dense-cell serialization already performed poorly.
- Do not infer a PureSoA advantage from the frozen-position scratch test or
  sort every step without amortizing sorting and migration costs.
- Mixed precision, relaxed tolerances, fewer correction passes, and disabled
  physics require separate accuracy studies and are outside this plan.

## Correctness and acceptance gates

For each kernel or communication change, run per-species comparisons against
the atomic paths, including all matrix entries, charge, and ten moments.
Use the existing checks' tolerances; neutral-species cancellation makes a
relative error on summed charge a poor sole metric. Include magnetized and
unmagnetized cases, uniform/clustered distributions, low/high occupancy,
empty cells, and particles at periodic corners and box faces.

Run the existing suite after each accepted implementation phase:

```bash
python3 src_gpu/validate_demo.py --mpi --physics \
  --reference bin/gpu_reference.exe --output-dir /tmp/fleks-gpu-perf-validation
```

Also check the CPU fallback and, for new cooperative kernels, memory/race
sanitizers where available. Preserve conservation, field/operator agreement,
plasma evolution, true residuals, and decomposition invariance. Solver changes
need difficult nonuniform cases, not only uniform oscillations.

Publish complete timestep results alongside phase timings, iteration counts,
peak memory, and uncertainty in `BENCHMARKS.md`. A candidate must exceed run
variability and improve its target workload without unexplained regressions;
retain an opt-in path when a tradeoff is workload-dependent. Re-rank the
remaining phases after each accepted improvement.
