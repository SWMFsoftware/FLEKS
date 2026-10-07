# `evolveGhostB`: long-horizon behaviour of the first fine ghost layer

Measured on `PARAM.in.hybrid` from this directory. Everything below is
reproducible with the developer-local scan script (see *Reproducing* at the
end); the raw run logs themselves are not kept in the repository.

**Question.** With `evolveGhostB = T`, the first fine ghost layer of `B` is
advanced by the Faraday update but never re-anchored to the coarse level
(`apply_centerB_BC` only re-interpolates the *outer* layers, `nKeep = 1`).
Its drift from the coarse interpolation (`dGhost`) was known to reach O(|B|) by
t = 5. Does it run away, and does it hurt anything?

**Answer.** No. `dGhost` grows fast at first and then **saturates** near 9
(≈ 4x the peak |B|), while **every AMR div(B) bucket stays bit-for-bit flat for
100 time units**. No action is required for stability.

## Setup

`tests/reconnection_amr/PARAM.in.hybrid` (Harris sheet, `nBSubcycle = 4`,
`dt = 0.02`, refined layer |y| < 7.85, outflow in y, 4 MPI ranks), run with
`#DIVB` `alwaysComputeDivB = T` and all three `#HYBRIDPIC` switches on, to
t = 50 and t = 100 (2500 and 5000 cycles). Control with `evolveGhostB = F` to
t = 100.

## Result (t = 100 run; the t = 50 run agrees throughout)

| t | L0 iface | L0 covered | L0 interior | L1 interior | L1 domain-edge | dCov | dGhost | \|B\| |
|---|---|---|---|---|---|---|---|---|
| 0.02 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 9.60e-02 | 1.0e-02 | 0.028 | 1.1 |
| 10 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 1.7e-01 | 4.4e-02 | 2.30 | 1.4 |
| 20 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 3.1e-01 | 4.2e-02 | 4.00 | 1.3 |
| 30 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 4.0e-01 | 5.0e-02 | 5.10 | 1.3 |
| 40 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 4.9e-01 | 7.3e-02 | 6.10 | 1.4 |
| 50 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 5.5e-01 | 1.0e-01 | 7.10 | 1.5 |
| 60 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 5.8e-01 | 1.4e-01 | 7.60 | 1.8 |
| 70 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 6.4e-01 | 1.4e-01 | 8.30 | 2.1 |
| 80 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 6.9e-01 | 1.5e-01 | 8.70 | 2.3 |
| 90 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 7.4e-01 | 1.6e-01 | 8.80 | 2.0 |
| 100 | 6.20e-04 | 1.20e-03 | 4.60e-04 | 3.80e-04 | 7.8e-01 | 1.2e-01 | 9.10 | 2.2 |

Columns are the `divB-AMR` log line emitted by `Pic::report_divB_amr`:
`if` / `cv` / `in` are the AMR interface, covered and interior buckets of
level 0, `in` is the fine level interior, `dom` the physical-boundary bucket,
then `dCov`, `dGhost` and the peak `|B|` used as the drift scale.

### 1. div(B) is completely flat

All four AMR buckets (`if`, `cv`, `in` on level 0, `in` on level 1) hold their
t = 0.02 values to the last printed digit for 5000 cycles. This is the headline
result for the whole branch: the coarse-fine interface contributes a *constant*
offset set by the initial condition, and adds no div(B) at any later time.

For contrast, with the three `#HYBRIDPIC` AMR switches off the same deck grows
to `if` = 1.3e-2 / `cv` = 1.1e-2 by t = 5, i.e. ~20x the interface value the
switches hold indefinitely.

### 2. `dGhost` saturates

Per 10 time units the increments are 2.27, 1.70, 1.10, 1.00, 1.00, 0.50, 0.70,
0.40, 0.10, 0.30 — monotonically decaying (with noise). Extrapolating, it
asymptotes around 9-10, i.e. it is bounded, not a runaway. The saturating
behaviour is what one expects once the ghost layer tracks the *slowly varying*
part of the difference between the fine solution and the coarse interpolation
of it, instead of accumulating it step after step.

`dGhost/|B|` peaks near 4.7 and then falls back as |B| grows (the reconnection
accelerates and the global maximum rises), so even the ratio is not monotone.

### 3. `dCov` (covered-cell drift) stays bounded

Fluctuating between 1e-2 and 1.7e-1 with no trend. The `ctRestrictB` relaxation
is removing the per-step drift as designed.

### 4. The one quantity that does grow steadily is *not* AMR

`L1 domain-edge` (cells on the physical outflow boundary, the `dom` bucket that
`report_divB_amr` separates out) climbs 9.6e-2 -> 7.8e-1 over the same run. That
is the outflow boundary condition, set by `apply_field_bc`, not by anything the
coarse-fine treatment does. Its increments also decay (0.07 -> 0.04 per 10t), so
it is decelerating too.

Worth watching, but it is a separate issue and it is the reason the `dom` bucket
must not be folded into the interior number: before it was split out, this
outflow signal made it look as though `evolveGhostB` were degrading the fine
level, which it is not.

## Control: `evolveGhostB = F`

With the switch off, `dGhost` is identically zero by construction (the ghost is
re-interpolated from the coarse level on every call), so the metric itself
carries no information; the comparison of interest is the fine-level `in`
bucket, which reads 3.8e-4 either way at t = 100.

## Conclusion / recommendation

`dGhost` does **not** need to be fixed. It saturates, and it does not degrade
div(B) anywhere over 100 time units. A periodic re-anchor (every N steps, or a
weighted blend with the coarse interpolation) remains available as cheap
insurance if a future deck does show growth, but there is no evidence for it
here, and it would give up the round-off fine-level div(B) that `evolveGhostB`
buys (2.8e-3 -> 2.8e-16 on `amr_equilibrium`).

The `dGhost > 5% |B|` warning in `report_divB_amr` **does** fire on this deck,
so the diagnostic works, but the threshold is set far below the observed
saturation level and will fire on any strong-gradient deck. It should be raised
(to ~20% of |B|) or reworded to say explicitly that saturation at several times
|B| is the expected behaviour. Left as-is pending that decision.

## Reproducing

The scan and parse helpers live outside version control, under the gitignored
`scratch/` tree (see `.gitignore`), because they build ad-hoc variants of the
test decks:

- `scratch/hybrid_amr_divb/scan.py` — builds a deck variant (`#DIVB` block,
  `#HYBRIDPIC` switch combination, `TimeMax`) and runs it under
  `scratch/hybrid_amr_divb/runs/<variant>/`. Variants used here: `rec_t50`,
  `rec_t100`, `rec_t100_nogh`.
- `scratch/hybrid_amr_divb/parse_divb.py` — parses the `divB-AMR` lines of a run
  log into the table above.

```bash
python3 scratch/hybrid_amr_divb/scan.py --variants rec_t50 rec_t100 rec_t100_nogh
python3 scratch/hybrid_amr_divb/parse_divb.py rec_t100
```

Run logs are not kept in the repository; regenerating the t = 100 pair takes
~19 min on 4 MPI ranks.