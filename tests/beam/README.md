# Ion Beam-Beam Instability Test

## Description
This is a standalone test designed to verify the simulation of the ion beam-beam instability in FLEKS. The test initializes a uniform background plasma moving in the $-x$ direction and a beam of ions of a smaller density moving in the $+x$ direction, solving for electromagnetic field dynamics. The test directory holds these parameter files:

- **`PARAM.in.implicit`** — semi-implicit full PIC: kinetic electrons and ions (`#FIELDSOLVER implicit`).
- **`PARAM.in.explicit`** — explicit full PIC: kinetic electrons and ions with explicit Maxwell advance (`#FIELDSOLVER explicit`).
- **`PARAM.in.hybrid`** — hybrid PIC: kinetic ions / massless fluid electrons (`#FIELDSOLVER hybrid`).
- **`PARAM.in.implicit.30s.waves`** — the same instability in a long thin box (2000 x 2 cells) with a warm plasma and the semi-implicit solver, run long enough for the transverse wave to grow from noise and saturate (see [Long-thin-box waves variant](#long-thin-box-waves-variant)).

The `implicit`, `explicit` and `hybrid` variants share the identical grid, guide field, beam ratio/speeds, density, normalization and output, so the same FFT-based validation applies.

## Physics & Solver Setup
- **Geometry & Boundaries**: 1D periodic grid (64 cells in $X$, 1 in $Y$ and $Z$); the waves variant uses a 2000 x 2 grid instead, see [below](#long-thin-box-waves-variant).
- **Physical Model**: A single charged species ($m = 1.0$, $q = 1.0$) with uniform background density $5.0\,\text{amu/cc}$ and background velocity $-400\,\text{km/s}$ ($u_x = -0.4$ in normalized units).
- **Beam setup**: $1\%$ of the particles have a beam velocity of $+400\,\text{km/s}$ ($u_x = +0.4$ in normalized units).
- **Electromagnetic Fields**:
  - Semi-implicit full PIC: `#FIELDSOLVER implicit` with a guiding magnetic field $B_x = 5\,\text{nT}$ ($5.0 \times 10^{-9}$ T).
  - Explicit full PIC: `#FIELDSOLVER explicit`, Boris-Maxwell leapfrog advance.
  - Hybrid PIC: `#FIELDSOLVER hybrid`, Hall term ON, resistivity and electron-pressure-gradient OFF so the instability grows undamped; the guide field is dynamically weak ($v_A \ll u_{\mathrm{beam}}$) so the explicit hybrid advance is not stiff at this resolution.
- **Purpose**: Verifies that the particle count and mean velocity are conserved, and tracks the growth of cyclotron waves in the transverse magnetic fields $B_y$ and $B_z$ from seed noise.

## Long-thin-box waves variant

`PARAM.in.implicit.30s.waves` keeps the beam setup but changes the configuration
so the transverse wave becomes a resolved, box-filling mode instead of a
box-limited one (see `_check_beam_transverse_wave`):

- 2000 x 2 cells over a 64 km box ($\Delta x = 32$ m), 100 particles/cell,
  periodic in all directions, plasma at $T = 72\,500$ K with a $0.48\,u_{\text{norm}}$
  ($960$ km/s) beam and the semi-implicit solver with `#UPWINDE`, `#UPWINDB`,
  `#SMOOTHJ`, `#DIVE` and `#COMOVING` enabled.
- The initial transverse field is exactly zero, so the wave is driven purely by
  the instability: `max|B_perp|` rises to about $5\%$ of the $5$ nT guide field
  by $t \approx 2$ (code time) and then stays saturated, with its power in
  modes $n = 130$-$440$ ($\ge 80\%$ below $n = N/4$).
- `TimeMax = 4` covers the growth phase plus two saturated units. The deck used
  to run to $t = 121$ (about 22 min on 2 MPI ranks); raising `TimeMax` restores
  that late-time study, which is a research run, not a CI test.

## Validation

The validator checks that the total particle count is conserved, the mean
velocity `MeanVx` is conserved within 1%, and that transverse `By`/`Bz`
waves grow from noise (power concentrated in low-order box-compatible modes). The hybrid variant additionally runs the shared hybrid energy checks.

The waves variant saturates (so the generic monotonic-growth check does not
apply) and is validated by `_validate_waves_plot` instead: the first frame must
be exactly field free, the wave must reach at least $1\%$ of the guide field
with at least $40\%$ of its power below $n = N/4$ (a coherent mode, not
grid-scale noise), the amplitude must stay below $20\%$ of the guide field and
must not collapse in the last frame, $B_x$ must stay uniform to $10^{-4}$
relative, and the ion density must round-trip to the `#UNIFORMSTATE` value.
The energy log is also guarded for finiteness and a steady `Eb`.

## Running

All variants run together under the test name:

```bash
python3 tests/validate_tests.py --test=beam          # all non-expensive variants
```

To run a **single variant**, append its token to the test name:

```bash
python3 tests/validate_tests.py --test=beam.implicit  # semi-implicit full PIC
python3 tests/validate_tests.py --test=beam.explicit  # explicit full PIC
python3 tests/validate_tests.py --test=beam.hybrid    # hybrid PIC
python3 tests/validate_tests.py --test=beam.implicit.30s.waves  # long thin box
```

The output table lists `BEAM (IMPLICIT)`, `BEAM (EXPLICIT)`,
`BEAM (HYBRID)` and `BEAM (IMPLICIT.30S.WAVES)` results. The long-thin-box
variant is the slowest of the four (about 55 s on 2 MPI ranks; the other three
take 1-2 s).
