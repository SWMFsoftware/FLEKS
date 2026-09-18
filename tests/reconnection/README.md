# Magnetic Reconnection Standalone Tests

This directory contains standalone magnetic reconnection test suites in FLEKS across multiple physical configurations and field solvers:

1. **`PARAM.in`** — **Fadeev full-PIC**: Force-free Fadeev current-sheet island equilibrium with kinetic ions and kinetic electrons ($m_i/m_e = 25$, Maxwell/GMRES solver).
2. **`PARAM.in.hybrid`** — **Fadeev hybrid-PIC**: Kinetic ions + massless fluid electrons with generalized Ohm's law.
3. **`PARAM.in.gem`** — **Classic GEM Challenge full-PIC**: The standard GEM reconnection benchmark (Birn et al. 2001) with a Harris current sheet, conducting walls in $y$, and a central magnetic perturbation.
4. **`PARAM.in.asym`** — **Asymmetric full-PIC**: Double current sheet reconnection with asymmetric magnetic fields ($B_1 = 1.0, B_2 = 2.0$) and temperatures ($T_1 = 1.33, T_2 = 3.33$) in a periodic domain.

## Coordinate Mapping

FLEKS uses a fake-2D convention with 1 cell along $z$:
- $x$: Reconnection outflow / periodic drive direction
- $y$: Current-sheet normal direction
- $z$: Out-of-plane / current / guide-field direction

## Physical Configurations

### 1. Fadeev Equilibrium (`PARAM.in`, `PARAM.in.hybrid`)
Uses `#TESTCASE fadeev` (`FadeevIC`):
```
Bx = b0 * sinh(y/L) / (cosh(y/L) + eps * cos(x/L)) + perturbation
By = b0 * eps * sin(x/L) / (cosh(y/L) + eps * cos(x/L)) + perturbation
Bz = b0 * bg (= 0)
n(x,y) = nb + (1 - nb) * profile(x,y) / profile_max
```
With $L = 5\,d_i$, $\epsilon = 0.4$, 2 islands, $n_b = 0.2$. The sheet is bounded by outflow field boundaries in $y$.

### 2. Classic GEM Challenge (`PARAM.in.gem`)
Uses `#TESTCASE gem` with `useStandardGem = T`:
```
Bx0 = B0 * tanh(y / lambda0)
n(y) = nb + (pB - Bx0^2) / (2 * Tp)
```
- Domain: $[-12.8, 12.8] \times [-6.4, 6.4]\,d_i$ ($64 \times 32 \times 1$ cells)
- Harris sheet thickness $\lambda_0 = 0.5\,d_i$, $B_0 = 1.0$, $T_{tot} = 1.0$ ($T_i/T_e = 5$), $n_b = 0.2$
- Background density $n_b = 0.2$, peak sheet density $n_0 + n_b = 0.5 + 0.2 = 0.7$
- Conducting field walls and reflecting particle walls at $y = \pm 6.4\,d_i$, periodic in $x$
- Perturbation: $\mathbf{B}_1 = \hat{\mathbf{z}} \times \nabla \psi$ with $\psi(x,y) = \psi_0 \cos(2\pi x/L_x) \cos(\pi y/L_y)$

### 3. Asymmetric Double Current Sheet (`PARAM.in.asym`)
Uses `#TESTCASE gem` with `isAsymmetryReconnection = T`:
- Domain: $[-32, 32] \times [-14, 14]\,d_i$ ($64 \times 28 \times 1$ cells), fully periodic in $x$ and $y$
- Two current sheets at $y = \pm 7\,d_i$ ($y = \pm 0.25 W_y$)
- Asymmetric fields: $B_1 = 1.0$ (interior, $y \in [-7, 7]$) and $B_2 = 2.0$ (exterior, $y \notin [-7, 7]$)
- Asymmetric temperatures: $T_1 = 1.33$ and $T_2 = 3.33$
- Localized Gaussian perturbation centered on the sheets

## Running

Run all reconnection test variants together:
```bash
# Serial (default)
python3 tests/validate_tests.py --test=reconnection

# MPI (e.g. 2 ranks)
python3 tests/validate_tests.py --test=reconnection -n 2
```

Run a single variant:
```bash
python3 tests/validate_tests.py --test=reconnection.full     # Fadeev full-PIC
python3 tests/validate_tests.py --test=reconnection.hybrid   # Fadeev hybrid-PIC
python3 tests/validate_tests.py --test=reconnection.gem      # Classic GEM challenge
python3 tests/validate_tests.py --test=reconnection.asym     # Asymmetric reconnection
```

## Validation Checks (`validate.py`)

1. **Energy Log Sanity**:
   - Kinetic ion energy $E_{part}$ and magnetic energy $E_b$ are finite (no NaN/Inf).
   - $E_b$ remains bounded (no numerical blow-up).

2. **Equilibrium Initialization ($t=0$)**:
   - **Fadeev**: In-plane field nulls (O-points) located at $x \approx \pm \pi L \approx \pm 15.7\,d_i$; peak sheet density $\approx 1$, background $\approx 0.2$.
   - **GEM Challenge**: Harris sheet field reversal ($B_x \to \pm 1.0$ at top/bottom boundaries); central X-point null at $x \approx 0$; peak density in $(0.5, 1.0)$, background $< 0.4$.
   - **Asymmetric**: Central region field $B_x \approx +1.0$ ($B_1$), outer boundary field $B_x \approx -2.0$ ($-B_2$); sheet density enhancement in $(0.5, 1.6)$, background $< 0.4$.

3. **Reconnection Dynamics**:
   - Seeded in-plane field perturbation grows nonlinearly ($\delta B_y$ increases).
   - Out-of-plane flux function $A_y$ at the X-point demonstrates active magnetic flux reconnection.

4. **Quasi-Neutrality**:
   - Initial condition verified to satisfy quasi-neutrality $|n_i - n_e| / n_0 < 0.5$ ($\rho_{S0} \approx 25\,\rho_{S1}$), verifying macroparticle charge scaling and electron mass loading.
