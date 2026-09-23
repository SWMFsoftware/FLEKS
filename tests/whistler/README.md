# HYBRID_WHISTLER

Parallel whistler–Alfvén wave test in FLEKS.  This directory runs **two
variants of the same Helicon-type transverse wave**:

| Variant | Config file | Field solver | Species | Expected dispersion |
| --- | --- | --- | --- | --- |
| **Full PIC** | `PARAM.in` | Maxwell / θ-scheme (`solveEM=T`) | kinetic ions **and** kinetic electrons | full cold plasma (electron inertia + displacement current) |
| **Hybrid** | `PARAM.in.hybrid` | generalized Ohm's law (`useHybridPIC=T`) | kinetic ions + massless fluid electrons | Hall / massless-electron limit |

Both variants are seeded with the same x-aligned, transverse, circularly
polarized **right-hand (whistler)** wave (`#TESTCASE HybridWave`,
`#WAVEIC rightHand T`, `frac = 0.02`).

> **The two variants are *not* directly comparable at this normalization.**
> The decks set `cLight = uNormSI`, and with `d_i ~ 1` in code units this gives
> `c / v_A ≈ 0.98` — i.e. the displacement current is *not* negligible, so the
> full-PIC run measures the full electromagnetic R-mode, not the Hall
> (massless-electron, no displacement current) whistler.  See
> [Dispersion relations](#dispersion-relations) below.  The hybrid variant is
> the one that isolates the Hall term; the full-PIC variant verifies the
> Maxwell/θ-scheme against the full cold-plasma dispersion.

## Physics / normalization

All quantities are in normalized CGS/code units chosen so ion scales are O(1):

* Box length `Lx = 6.4` code units (`nCellX = 64`, `maxBlockSizeX = 64`).
* `lNormSI = uNormSI = 1.0e5`, so `tNorm = 1 s` and `Ω_i = 1` in code units.
* `ρ = 5 amu/cc` → `ρ_code = 0.07673`, `v_A = 1.0183`, and
  `d_i = 1/ω_pi = 1.0183` code units (FLEKS takes `c = 1` in code units, so
  `d_i = 1/√(4πρ)`).
* Guide field `B = 1.0439e-8 T` → `B_code ≈ 1 ≈ Ω_i` (with `q/m = c = 1`).
* Box mode `m` has `κ = k d_i = 2π m d_i / Lx = 1.000 m`, so mode 1 sits
  almost exactly at `k d_i = 1`.
* Seeded amplitude `frac = 0.02`, `walenFactor` chosen for the exact eigenmode
  (see below) — the seed is a **single** circular branch, so the mode amplitude
  is constant in time instead of beating.

### Dispersion relations

For a parallel-propagating circularly polarized wave with `κ = k d_i` the
Hall (massless-electron, no displacement current) branches — what the hybrid
solver must reproduce — are

```
omega/Omega_i = [+kappa^2 + kappa*sqrt(kappa^2+4)] / 2     right-hand / whistler
omega/Omega_i = [-kappa^2 + kappa*sqrt(kappa^2+4)] / 2     left-hand / ion-cyclotron
```

with limits `ω = k v_A` (κ→0) and `ω → Ω_i κ²` (whistler, i.e. `ω ∝ k²`) or
`ω → Ω_i` (ion-cyclotron resonance).  Keeping the displacement current and the
electron inertia, and using FLEKS code units (`c = Ω_i = d_i = 1`,
`μ = m_i/m_e`, `n_e/n_i`):

```
kappa^2/w^2 = 1 + (n_e/n_i)/(w (1 - w/mu)) - 1/(w (w+1))     right-hand
kappa^2/w^2 = 1 - (n_e/n_i)/(w (1 + w/mu)) - 1/(w (w-1))     left-hand
```

`tests/_shared/hybrid.py` picks the appropriate pair from the deck
(`#SOLVEEM T` → full cold plasma, otherwise Hall) and bisects for the root.

## Seeding

`#WAVEIC` (see `PARAM.XML`) builds

```
dB_y(x,0) = B1 cos(kx),   dB_z(x,0) = hand * B1 sin(kx),   B1 = frac * Bx0
u_perp    = -walenFactor * B_perp / B0        (code units, B0 ~ 1)
```

* `rightHand = T` sets `hand = -1`, i.e. the transverse field rotates
  `y -> z` about `+B0` (the electron-gyration sense) for a wave propagating
  along `+B0`.  This is the whistler helicity.  The default (`F`) is the
  left-hand/ion-cyclotron seed.
* `walenFactor = 1` is the incompressible (Alfvén) relation, which is the
  exact eigenmode only for `κ → 0`; at `κ ~ 1` it excites **both** branches
  (≈84 %/16 % amplitude for the whistler).  The exact single-branch seed needs

  ```
  walenFactor = kappa * (v_A/B0) / (omega/Omega_i)
  ```

  `dispersion.py` computes this automatically (`--seed eigenmode`, default).
  With it, `|C(t)|` is constant to ~1 % and the frequency fit residual drops
  from ~11 % to ~1 %.

## Running

```
python3 tests/validate_tests.py --test=whistler          # both variants
python3 tests/validate_tests.py --test=whistler.full     # full PIC only (PARAM.in)
python3 tests/validate_tests.py --test=whistler.hybrid   # hybrid only (PARAM.in.hybrid)
```

### Dispersion-relation scan

`dispersion.py` is the physics-level check: it runs a scan over wave modes
`m = 1 … 6` (i.e. `κ = 1 … 6`), tracks `δB_y(x0, t)` at a probe point, measures
`ω(k)`, and compares against the analytic whistler:

```bash
python3 tests/whistler/dispersion.py                          # hybrid, m = 1..4
python3 tests/whistler/dispersion.py --modes 1 2 3 4 5 6      # wider k range
python3 tests/whistler/dispersion.py --variant full           # full-PIC deck
python3 tests/whistler/dispersion.py --no-run                 # re-analyse runs
```

It writes `whistler_dispersion.png` (ω–k with both analytic branches, `dω/dk`,
`v_p`) and `whistler_waveform.png` (probe time series, `(δB_y, δB_z)` hodogram,
seeded-mode spectrum, measured-vs-analytic) and prints a PASS/FAIL table for:

1. **Right-hand polarization** — the transverse field at the probes rotates
   `y → z` about `+B0` (hodogram circulation positive), selecting the whistler
   branch.  Checked both on the seeded spatial harmonic and on the raw probe.
2. **`ω ∝ k²`** — the measured local `d ln ω / d ln k` follows the analytic
   value, which tends to 2 (`Ω_i κ²` asymptote) for large `κ`.
3. **Phase velocity** — `v_p = ω/k` agrees with the analytic whistler within
   the discrete truncation error of the compact curl/gradient stencils.
4. **Probe vs mode frequency** — the raw probe estimate agrees with the
   seeded-harmonic estimate (reported as noise-limited when the probe is
   buried in grid noise).

Measured (hybrid deck, `--modes 1 2 3 4 5 6 --periods 6`, `dn` per mode; all
four criteria PASS):

| m | κ | measured | analytic | error | local `d ln ω / d ln k` (meas / anal.) | hand |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.9998 | 1.6159 | 1.6174 | −0.09 % | – | R |
| 2 | 1.9995 | 4.8055 | 4.8262 | −0.43 % | 1.572 / 1.577 | R |
| 3 | 2.9993 | 9.8276 | 9.9034 | −0.77 % | 1.764 / 1.773 | R |
| 4 | 3.9991 | 16.7110 | 16.9357 | −1.33 % | 1.845 / 1.865 | R |
| 5 | 4.9988 | 25.4178 | 25.9495 | −2.05 % | 1.879 / 1.912 | R |
| 6 | 5.9986 | 35.8891 | 36.9544 | −2.88 % | 1.892 / 1.939 | R |

The residual error grows with `k dx` as expected for the compact discrete
operators (`k dx = 0.098 m`), and the measured curve sits consistently *below*
the continuum one, i.e. the discrete stencil under-estimates the whistler
frequency.  The local log-log slope climbs towards the `ω ∝ k²` asymptote
(2) from below, tracking the analytic slope to ≲ 2.5 %.  Phase velocity agrees
to ≤ 1.74 % and the raw probe-point frequency to ≤ 0.11 %, and all 6 seeded
periods are usable as one clean fit window (before the fix in gotcha 4 the
window had to be trimmed to 2.4–3.2 periods at `m = 1…3`).  The full-PIC
variant measures `ω/Ω_i = 0.790` against the full cold-plasma root `0.786`
(0.5 %).

## Validation

`validate.py` runs the shared checks from `tests/_shared/hybrid.py` for **both**
variants: seeded-wavelength (`n = 1`, ≥ 50 % of the non-DC power) and
bounded-amplitude checks, plus the whistler-dispersion check above.  For the
hybrid variant it additionally runs the shared hybrid energy-log checks (see
[hybrid README](../_shared/hybrid.md)).

## Gotchas (all verified on this test)

1. **Do not reintroduce any EMA time-averaging of B into the Ohm's law.**  The
   (now removed) `#AVGFIELDB` filter lagged the convection term `-u_i × B` by
   `τ ~ (nAvgFieldB-1)/2 · dt` and thus *anti-damped* the wave at `γ ~ ω²τ/2`.
   Measured on this deck with `nAvgFieldB = 20`: `γ = 0.20 /Ω_i` versus
   `0.014 /Ω_i` without it — the wave grew 13× in 13 s instead of holding its
   seed amplitude.
2. **The hybrid seed must match the branch.**  `wf_neg` (wrong sign of
   `walenFactor`) excites the ion-cyclotron branch and the probe amplitude
   beats to zero; `rightHand`/`walenFactor` are the two knobs.
3. **Output cadence.**  The whistler period at `κ = 1` is `T = 3.89 s`; with
   `dn = 100` (the old deck) the phase step per frame was `ω dt = 3.3 rad > π`
   and the phase fit aliased.  The decks now use `dn = 10`
   (`dt_frame = 0.2 s`, phase step 0.32 rad).
4. **The hybrid Hall advance needs grid-scale hyper-resistivity.**  The explicit
   Hall advance of this cold, low-β plasma is unstable to a non-propagating,
   symmetric `±n` mode at the few-cell scale.  Measured on this deck (`κ ≈ 1`,
   24 s, 64 cells): the strongest growth is at `λ ≈ 3` cells (`n ≈ ±19`) with
   `γ ≈ 0.65 /Ω_i`, while the physical range `k d_i ≲ 6` is stable.  It is pinned
   to the Hall-term discretization by two checks — it is independent of the
   particle count (20 000 ppc: unchanged) and it disappears with the Hall term
   off (`useHallTerm = F`: `γ = −0.007`) — while more B sub-cycling does not
   help (`nBSubcycle = 32`: `γ = +0.50`).  The deck therefore enables grid-mode
   hyper-resistivity (`#HYPERRESISTIVITY`, `etaHyperMode = grid`,
   `etaHyperCh = 1e-2`), which damps at `γ_h = etaHyperCh (k dx)^4 / dt` per
   step: ≈ `6.0 /Ω_i` at the unstable band (versus its `0.65`), but only
   `4.6e-5` at `k d_i = 1` and `0.06` at `k d_i = 6` (`< 1 %` amplitude loss per
   period).  With it, `max|B_perp|` grows by only 1.08–1.42× over the whole
   6-period run (the seed's own beat) instead of 6.9–24×, and the seed mode is
   untouched (`γ = 0.000` versus `−0.002`).  Anything `>= 1e-2` works equally
   well; `1e-3` is too weak (`γ` only falls to `0.42`).  Before this was
   understood the hybrid deck stopped at `TimeMax = 10 s` (2.6 periods) and both
   the validator and `dispersion.py` trimmed the fit window once `max|B_perp|`
   exceeded ~1.5× its seeded value; neither is necessary now.  (This option only
   exists in the hybrid solver — the full-PIC deck has no Ohm's law to damp.)
5. **Full-PIC electron density** must be quasi-neutral (`ρ_e = 5/1836
   amu/cc`).  The previous value (`2.5552e-5`) violated `n_e = n_i` by 106×,
   which moved the measured frequency from 0.79 to 1.24 and would have made the
   dispersion check meaningless.
6. **Never write a command keyword inside a *comment* line of a PARAM.in file.**
   The reader scans for the command keyword *anywhere* in the text, including
   mid-line, and then takes the parameter value from the first token of the
   **next** line.  A comment that referred to the Hall-term block by name
   (`#HALLTERM`, followed on the next line by the word "not" …) silently built a
   real `#HALLTERM` command with `useHallTerm = "not"` → false, which disabled
   the Hall term; the seeded mode then appeared to grow 21× in 24 s instead of
   holding its amplitude, and the k-scan returned `ω ≈ 0.38` for `m = 1`
   instead of 1.62.  Nothing warns you.  The quickest checks:
   * `grep -n '#' PARAM.in | grep -v '^[0-9]*:#'` — mid-line hashes;
   * the startup echo in the run log (`PC: #COMMAND` followed by each
     `<value> <name>`) lists every command the reader actually built.
   Spell out command names in prose ("the HALLTERM block", "nBSubcycle"),
   never as hash-prefixed tokens.  All other decks under `tests/` are clean.
