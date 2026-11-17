# External Coil, Toroidal Field, and Inductive Waveform Tests

## What Is Tested

The `coil_field` standalone suite validates the generic external magnetic field
primitives (`#COIL` and `#TOROIDALFIELD`), time-dependent current waveforms
(`#COILWAVEFORM`), inductive electric fields `E0 = -dA0/dt`, and the virtual
probe post-processor (`tools/extract_probes.py`) without using planetary
defaults (`#PLANETRADIUS` omitted):

1. **`PARAM.in.helmholtz_toroidal`** — Coaxial multi-turn `#COIL` solenoids in
   Helmholtz spacing combined with a regularized `1/R` `#TOROIDALFIELD` under
   3D `#FIELDSOLVER none`. Verifies analytic elliptic-integral and toroidal
   superposition in `real8` SI `.out` files.
2. **`PARAM.in.segment_loop`** — Closed square current loop built from four
   regularized `#COIL` `segment` conductors in ASCII `.out` format.
3. **`PARAM.in.filament2d`** — Fake-2D (`nCellZ = 1`) opposing `filament2d`
   pair driven by matching `linear_ramp` and `table` (`waveform_table.txt`)
   waveforms. Verifies analytic `(Bx, By)`, out-of-plane inductive `Ez`, and
   discrete 2D Faraday consistency `curl(E0) = -dB0/dt`.
4. **`PARAM.in.waveform_static`** — 3D tilted and axial `loop` coils driven by
   `cosine_ramp` and `sinusoid` waveforms in `#FIELDSOLVER none`. Verifies 3D
   vector potential inductive electric fields, discrete 3D Faraday consistency,
   and subprocess CLI extraction (`--line` CSV and `--plane` NPZ) using
   `tools/extract_probes.py`.
5. **`PARAM.in.waveform_hybrid`** — Time-dependent `filament2d` coils driven by a
   `crowbar_rlc` waveform inside a resistive `#FIELDSOLVER hybrid` plasma.
   Verifies hybrid `E0` / `Eind` diagnostics, crowbar rise-to-decay waveform
   evaluation, and inductive eddy response (`B1` opposing rising `B0`).

## Running Locally

```bash
python3 tests/validate_tests.py --test=coil_field -n 2
```
