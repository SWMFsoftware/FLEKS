# 2D Photoionization Shadow Cylinder Test (Mars-like Planet)

## Description

This test verifies the `#PHOTOIONIZATION` + `#SHADOWCYLINDER` commands on a 2D
grid with a Mars-like planet at the domain center.  The planet casts a shadow
cylinder on the **nightside (−X)**, creating a clear day/night asymmetry in the
photoionization source.

## Geometry & Grid

```
              Y ↑
                |
       +9e6 m ─ ┼ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─
                |       (dayside +X)
                |          ☀  Sun
                |    ·  ·  ·  ·  ·  ·  ·
                |    ·  ●═══════════●  ·
                |    ·  ║  Planet  ║  ·  ·
                |    ·  ║ R=3e6 m ║  ·  ·
       ─────────┼──────╬───●─────╬────────── X
                |    ·  ║  Shadow  ║  ·
     -9e6 m ─   |    ·  ║ (ν=0)   ║  ·
                |    ·  ╚══════════╝  ·
                |       (nightside −X)
       -9e6 m ─ ┼ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─
                |
```

| Parameter | Value | Note |
|-----------|-------|------|
| X range | [−9×10⁶, 9×10⁶] m | ±3 Rₚ |
| Y range | [−9×10⁶, 9×10⁶] m | ±3 Rₚ |
| Z range | [−0.5, 0.5] m | quasi-2D (1 cell) |
| Grid | 32×32×1 | dx ≈ 562 km |
| Planet radius Rₚ | 3.0×10⁶ m | Mars-like |
| Sun direction | +X | dayside = +X, nightside = −X |
| Shadow radius | 3.0×10⁶ m | = Rₚ |

## Exosphere (Chamberlain Model)

Two neutral components with Mars-like scale heights:

| Component | n₀ [m⁻³] | H₀ [m] | T₀ [K] | Description |
|-----------|----------|--------|--------|-------------|
| H | 1.0×10¹⁰ | 2.0×10⁶ | 300 | Light, large scale height |
| O | 5.0×10¹⁰ | 3.0×10⁵ | 300 | Heavy, small scale height |

Density profile: \( n(r) = n_0 \exp\left[-H_0\left(\frac{1}{R_p} - \frac{1}{r}\right)\right] \) for \( r \geq R_p \).

## Ionization

| Parameter | Value |
|-----------|-------|
| ν_photo₀ (H, O) | 1.0×10⁻⁶ s⁻¹ |
| Shadow cylinder | enabled, R = 3.0×10⁶ m |

Photoionization rate at distance r (outside shadow):

\[ \nu(r) = \nu_0 \left(\frac{R_p}{r}\right)^2 \]

Inside the shadow cylinder (x < 0 and proj·(−solarDir) > 0 and perp < R_shadow):
ν = 0.

## Plasma Species

| Species | Mass [amu] | Charge [e] | Role |
|---------|------------|------------|------|
| 0 (H⁺) | 1.0 | +1 | Background ions |
| 1 (O⁺) | 16.0 | +1 | Heavy ions (receives source) |

## Expected Results

1. **Energy growth**: Species 1 (O⁺) energy increases over time → source active.

2. **Day/night asymmetry** (verified with `--thorough`):
   - Output format: **IDL ASCII** (`z=0` slice, `.out` file) — human-readable,
     easy to load for inspection.
   - **Dayside (+X)**: `ppcS1` > 0 — photoionization produces O⁺ particles.
   - **Nightside shadow (−X, |y| < Rₚ)**: `ppcS1` ≈ 0 — zero photoionization
     inside the shadow cylinder.

3. **No crashes**: Simulation completes all 10 iterations.

## Running

```bash
# Basic (energy growth check only)
python3 tests/validate_tests.py --test=photoionization

# Thorough (includes day/night asymmetry check via flekspy)
python3 tests/validate_tests.py --thorough
```
