# Proposal: batched Gram-Schmidt (CGS2) in the GMRES Arnoldi step

Status: not scheduled. It changes CPU outputs, so it is kept out of the GPU
port (`feat/gpu-port`) and needs its own branch and new reference outputs.

## Problem

`linear_solver_gmres` (`src/LinearSolver.cpp`) orthogonalizes each new Krylov
vector with modified Gram-Schmidt (MGS). In Arnoldi step `i`, MGS computes
`i + 1` dot products one after the other, because each dot uses the vector
already updated by the previous projection:

```
for j = 0..i:
    h[j]  = <v_j, w>          # one reduction
    w    -= h[j] * v_j        # one kernel
h[i+1] = ||w||; w /= h[i+1]   # one reduction, two kernels
```

Each step therefore costs about `2(i + 1) + 3` kernel launches and `i + 2`
reductions. With `nKrylov` vectors per restart, that is O(nKrylov²) launches
and reductions per restart.

- **GPU, one rank:** the coefficients stay on the device, so there is no host
  round trip per dot. But every reduction is still a separate launch, and on
  the small grids of coupled runs (8³ boxes, few cells per rank) the solve is
  bound by launch latency, not bandwidth.
- **Several ranks (CPU or GPU):** each dot is a separate `MPI_Allreduce`, so
  a step pays `i + 2` network latencies.

## Proposal

Use classical Gram-Schmidt with one reorthogonalization pass (CGS2):

```
h   = V_i^T w           # all i+1 dots in one batched reduction
w  -= V_i h             # one kernel
h2  = V_i^T w           # second pass, one batched reduction
w  -= V_i h2; h += h2   # one kernel
h[i+1] = ||w||; w /= h[i+1]
```

This uses 2 batched reductions plus one norm per step, independent of `i`. With
MPI, the per-rank partial sums of a batch go in one `MPI_Allreduce` of `i + 1`
values. CGS2 is as stable as MGS in practice (orthogonality to working
precision), which is why it is the usual choice in GPU and distributed Krylov
solvers.

The batched dot product must stay deterministic: a fixed reduction tree per
coefficient, as in `deterministic_sum_to` (`include/GpuBackend.h`), so GPU runs
stay bitwise repeatable.

## Consequences

- **CPU outputs change.** The coefficients are computed in a different order,
  so the round-off differs. GMRES may also take a different number of
  iterations, so solutions can differ by up to the solver tolerance, not just
  round-off. All CPU reference outputs (standalone tests, SWMF tests, FLX
  coupled tests) need regenerating.
- **GPU vs CPU agreement** stays at round-off only if both use CGS2.
- **Cost:** CGS2 does about twice the flops of MGS (two passes). That is
  negligible next to the matvec, but it is not free on large, bandwidth-bound
  grids. Measure both small (coupled 8³) and large standalone decks.

## Options

1. Use CGS2 on all platforms and regenerate the references (one code path, as
   the GPU port's design rule asks).
2. Add a parameter (e.g. in `#EFIELDSOLVER`) to choose MGS or CGS2, with MGS as
   the default until the references are updated. This keeps the old outputs
   available, but leaves two methods to maintain.

## Acceptance

- Standalone and coupled tests pass with the new references; the field-solver
  convergence (iterations per solve, final residual) is no worse than MGS on
  the test decks.
- GPU runs are bitwise repeatable; GPU and CPU agree at round-off.
- Timing: launches and reductions per GMRES iteration no longer grow with the
  Krylov index; report `E_iterate` / `E_matvec` timings on test17 and AEPIC
  (coupled) and one large standalone deck, on the RTX 5070 and the A100.
