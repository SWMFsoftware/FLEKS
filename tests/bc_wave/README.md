# Wave-Injection Boundary-Condition Tests

This directory groups several wave-injection tests (`#WAVEBC`).  Each variant
is a `PARAM.in.<suffix>` file; the runner discovers all of them and runs each
once.  A single `validate.py` branches on the variant name.

## Variants

| Variant | File | Physics | Validation |
|---------|------|---------|------------|
| `bc_inject_mono` | `PARAM.in.mono` | Monochromatic Bz wave injected at the `x-lo` face into a pure-EM domain (`#SOLVEEM T`, `lightwave`) | Bz perturbation reaches the interior |
| `bc_inject_alfven` | `PARAM.in.alfven` | Shear Alfvén wave (By + matching velocity kick) injected at the `-x` face into an initially uniform magnetized hybrid plasma, propagating along +x | By reaches the interior; By/Bx < 10% (linear) |
| `bc_inject_window` | `PARAM.in.window` | One out-of-plane shear Alfvén wave (Bz, consistent E, Walén kick; lambda = 16 d_i) injected at the `+x` face of a flowing solar wind (full PIC, 2D), with `#WAVEBCWINDOW` limiting the source to `0 < y < 10` d_i | Bz rms outside the window band < 0.2 of inside (outside is particle noise) |
| `bc_inject_multiface` | `PARAM.in.multiface` | Two `wave` faces (`y:hi` listed before `x:hi` in `#WAVEBC`) with opposite-sign `Bz` and `uz` mode tables on a 2x2-block non-periodic 2D grid | Non-wave block bounds stay zero; `y:hi` and `x:hi` get their distinct kicks; `(x:hi, y:hi)` corner follows `#WAVEBC` face order (`y:hi`) for both `Bz` and `uzS0` |
| `bc_inject_turbulent` | `PARAM.in.turbulent` | Broadband out-of-plane shear-Alfvén spectrum (`#TURBULENTINFLOW`, 8 modes, `lambda = [8, 16]` $d_i$, `dBrms/B0 = 0.25`) injected at the `+x` face of a flowing solar wind with periodic `y` | `wave_modes.json` written with periodic-snapped $k_y$; interior $B_z$ RMS matches spectrum; Walén correlation $r \ge 0.85$ and slope $\approx 1$ |

## Validation

Both variants exercise the **wave-injection boundary condition** (`BC::wave` +
`#WAVEBC`): a wave seeded at a boundary face must propagate into the interior.
`validate.py` checks (branching on the variant name):
1. **Energy log** (both): `Etot` finite, positive, approximately conserved.
2. **`bc_inject_mono` plot**: a monochromatic transverse `Bz` wave injected at a
   boundary face propagates into a pure-EM domain; `Bz` non-zero in the interior.
3. **`bc_inject_alfven` plot**: a self-consistent shear Alfvén wave (`δBy` +
   `δv_y = -δB_y/sqrt(μ0 ρ)`, about -1.02 km/s; `#WAVEBC` amplitudes are SI)
   injected at the `-x` face
   into an initially uniform magnetized hybrid plasma propagates along `+x`;
   `By` non-zero in the interior and `By/Bx < 10%` (linear wave).

## Running

All variants run together under the test name:

```bash
python3 tests/validate_tests.py --test=bc_wave          # all variants
```

To run a **single variant**, append its token to the test name:

```bash
python3 tests/validate_tests.py --test=bc_wave.mono     # monochromatic EM wave only
python3 tests/validate_tests.py --test=bc_wave.alfven   # shear Alfvén wave only
```
