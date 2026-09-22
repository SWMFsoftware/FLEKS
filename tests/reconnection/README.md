# Magnetic Reconnection Standalone Tests

This directory contains standalone magnetic reconnection test suites in FLEKS across multiple physical configurations and field solvers:

1. **`PARAM.in.fadeev_pic`** — **Fadeev full-PIC**: Force-free Fadeev current-sheet island equilibrium with kinetic ions and kinetic electrons ($m_i/m_e = 25$, Maxwell/GMRES solver).
2. **`PARAM.in.fadeev_hybrid`** — **Fadeev hybrid-PIC**: Kinetic ions + massless fluid electrons with generalized Ohm's law.
3. **`PARAM.in.gem_pic`** — **Classic GEM Challenge full-PIC**: Standard GEM reconnection benchmark (Birn et al. 2001) with a Harris current sheet, conducting walls in $y$, and a central magnetic perturbation.
4. **`PARAM.in.gem_hybrid`** — **Classic GEM Challenge hybrid-PIC**: Kinetic ions + isothermal fluid electrons with generalized Ohm's law, stationary ion background (`useUniformIonPressure = T`), and electron fluid carrying the diamagnetic current.
5. **`PARAM.in.asym_pic`** — **Asymmetric full-PIC**: Double current sheet reconnection with asymmetric magnetic fields ($B_1 = 1.0, B_2 = 2.0$) and temperatures ($T_1 = 1.33, T_2 = 3.33$) in a periodic domain.
6. **`PARAM.in.forcefree_hybrid`** — **Force-Free Sheet hybrid-PIC**: Force-free current sheet reconnection (Le et al. 2016, WarpX benchmark) with uniform plasma density, uniform total magnetic pressure $B_0^2 + B_g^2$, and $4\pi$-consistent scale normalization.
7. **`PARAM.in.forcefree_pic`** — **Force-Free Sheet full-PIC**: Force-free current sheet reconnection with kinetic ions and electrons ($m_i/m_e = 25$, $10 \times 10$ PPC per species, Maxwell/GMRES solver).

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

### 2. Classic GEM Challenge (`PARAM.in.gem_pic`, `PARAM.in.gem_hybrid`)
Uses `#TESTCASE gem` with `useStandardGem = T`:
```
Bx0 = B0 * tanh(y / lambda0)
n(y) = nb + (pB - Bx0^2) / (2 * Tp)
```
- **Scale Normalization & $4\pi$ Factor**: FLEKS code units are based on Gaussian / CGS units where:
  $$\nabla \times \mathbf{B} = 4\pi \mathbf{J}, \quad \mathbf{J} \times \mathbf{B} = \frac{(\nabla \times \mathbf{B}) \times \mathbf{B}}{4\pi}, \quad P_{\text{mag}} = \frac{B_0^2}{8\pi}$$
  For reference parameters $B_0 = 1.0, \rho_0 = 1.0, q=1, m=1$:
  - $\Omega_{ci} = q B_0 / m = 1.0 \implies \tau_{ci} = 2\pi / \Omega_{ci} \approx 6.2831853$ s
  - $v_A = B_0 / \sqrt{4\pi \rho_0} = 1/\sqrt{4\pi} \approx 0.2820948$ code units ($28.21$ km/s)
  - $d_i = v_A / \Omega_{ci} = 1/\sqrt{4\pi} \approx 0.2820948$ code units ($28.21$ km)
  - $P_{\text{mag}} = \frac{B_0^2}{8\pi} = \frac{1}{8\pi} \approx 0.0397887$ code units
- **GEM Challenge Specifications (Birn et al. 2001)**:
  - **Domain**: $L_x = 25.6\,d_i = 25.6 \times 0.2820948 \approx 7.2216$ code units ($x \in [-3.61081, 3.61081]$), $L_y = 12.8\,d_i = 12.8 \times 0.2820948 \approx 3.6108$ code units ($y \in [-1.80541, 1.80541]$).
  - **Harris sheet half-thickness**: $\lambda_0 = 0.5\,d_i = 0.5 \times 0.2820948 \approx 0.14105$ code units.
  - **Perturbation**: $\psi_0 = 0.1 B_0 d_i = 0.1 \times 1.0 \times 0.282095 \approx 0.02821$ code units (`apert = 0.02821`), yielding dimensionless relative perturbation $\delta B_y / B_0 = \psi_0 \frac{2\pi}{L_x} \approx 0.02454$.
  - **Densities**: Background $n_b = 0.20$ code units ($13.03188$ amu/cc), sheet increment $n_0 = 1.00$ code units (total peak $n(0) = n_b + n_0 = 1.20$), reference $n_{e0} = 1.0$ code unit (`electronDensity0 = 65.15939` amu/cc).
  - **Temperatures**: Magnetic pressure balance requires $n_0 (T_i + T_e) = \frac{B_0^2}{8\pi} = \frac{0.5}{4\pi}$. With $T_i / T_e = 5$:
    - $T_{i,\text{code}} = \frac{5}{6} \times \frac{0.5}{4\pi} = \frac{5}{48\pi} \approx 0.033157$ code units $\implies T_i \approx 40,168$ K.
    - $T_{e,\text{code}} = \frac{1}{6} \times \frac{0.5}{4\pi} = \frac{1}{48\pi} \approx 0.0066315$ code units $\implies T_e \approx 0.6923$ eV (hybrid) or $8,033.6$ K (full-PIC).
  - **Grid**: $64 \times 32 \times 1$ cells ($\Delta x = \Delta y = 0.40\,d_i \approx 0.1128$ code units).
  - **Boundaries**: Periodic in $x$, conducting walls for fields and reflecting walls for particles at $y = \pm 1.80541$ (matching Birn et al.).
  - **Hybrid PIC configuration**: For `PARAM.in.gem_hybrid`, setting `useUniformIonPressure = T` keeps ions stationary ($u_{iz} = 0$), so the fluid electrons carry the diamagnetic current through generalized Ohm's law.

### 3. Asymmetric Double Current Sheet (`PARAM.in.asym`)
Uses `#TESTCASE gem` with `isAsymmetryReconnection = T`:
- Domain: $[-32, 32] \times [-14, 14]\,d_i$ ($64 \times 28 \times 1$ cells), fully periodic in $x$ and $y$
- Two current sheets at $y = \pm 7\,d_i$ ($y = \pm 0.25 W_y$)
- Asymmetric fields: $B_1 = 1.0$ (interior, $y \in [-7, 7]$) and $B_2 = 2.0$ (exterior, $y \notin [-7, 7]$)
- Asymmetric temperatures: $T_1 = 1.33$ and $T_2 = 3.33$
- Localized Gaussian perturbation centered on the sheets

### 4. Force-Free Current Sheet (`PARAM.in.forcefree_hybrid`, `PARAM.in.forcefree_pic`)
Uses `#TESTCASE forcefree` (`ForceFreeIC` / `#FORCEFREEIC`), ported from the WarpX benchmark (Le et al. 2016):
```
Bx(x,y) = b0 * tanh(y/lambda) + deltaBx
By(x,y) = deltaBy
Bz(x,y) = sqrt(bg^2 + b0^2 * sech^2(y/lambda))
n(x,y)  = n0  (uniform)
```
with divergence-free magnetic perturbation:
$$\delta B_x = -\delta B \frac{L_x}{2 L_y} \cos\left(\frac{2\pi x}{L_x}\right) \sin\left(\frac{\pi y}{L_y}\right)$$
$$\delta B_y = \delta B \sin\left(\frac{2\pi x}{L_x}\right) \cos\left(\frac{\pi y}{L_y}\right)$$
- **Scale Normalization & $4\pi$ Factor**: FLEKS code units are based on Gaussian / CGS units where $\nabla \times \mathbf{B} = 4\pi \mathbf{J}$, $\mathbf{J} \times \mathbf{B} = \frac{(\nabla \times \mathbf{B}) \times \mathbf{B}}{4\pi}$, and magnetic pressure $P_{\text{mag}} = \frac{B_0^2}{8\pi}$. For $B_0 = 1.0, \rho_0 = 1.0, q=1, m=1$:
  - $\Omega_{ci} = q B_0 / m = 1.0 \implies \tau_{ci} = 2\pi / \Omega_{ci} \approx 6.2831853$ s
  - $v_A = B_0 / \sqrt{4\pi \rho_0} = 1/\sqrt{4\pi} \approx 0.2820948$ code units ($28.21$ km/s)
  - $d_i = v_A / \Omega_{ci} = 1/\sqrt{4\pi} \approx 0.2820948$ code units ($28.21$ km)
  - $d_i \times v_A = \frac{1}{4\pi} \approx 0.0795775$ code units
  - $P_{\text{mag}} = \frac{B_0^2}{8\pi} = \frac{1}{8\pi} \approx 0.0397887$ code units
- **Force-Free Property**: Because $B_x^2 + B_z^2 = b_0^2 \tanh^2(y/\lambda) + b_g^2 + b_0^2 \operatorname{sech}^2(y/\lambda) \equiv b_0^2 + b_g^2 = \text{const}$, the unperturbed magnetic pressure is uniform everywhere, requiring no plasma pressure gradient for mechanical equilibrium ($\mathbf{J} \times \mathbf{B} = 0$).
- **Domain**:
  - `forcefree_hybrid`: $[-5.6419, 5.6419] \times [-2.82095, 2.82095]$ code units ($40\,d_i \times 20\,d_i$) on a $128 \times 64 \times 1$ grid ($dx = dy = 0.3125\,d_i \approx 0.08815$ code units; production WarpX uses $512 \times 512$).
  - `forcefree_pic`: $[-2.82095, 2.82095] \times [-1.41047, 1.41047]$ code units ($20\,d_i \times 10\,d_i$) on a $64 \times 32 \times 1$ grid ($dx = dy = 0.3125\,d_i$).
- **Physical Parameters**: $\lambda = 1.0\,d_i = 0.2821$ code units, $b_0 = 1.0$, $b_g = 0.3$, $\delta B = 0.01$ (hybrid) or $0.05$ (full-PIC), uniform $n_0 = 1.0$ ($65.15939$ amu/cc), $T_i = 24,191$ K ($\beta_i = 0.50$, $P_{i,\text{code}} = 0.25/(4\pi)$), $T_e = 0.4154$ eV ($\beta_e = 0.10$, isothermal, $P_{e,\text{code}} = 0.05/(4\pi)$), normalized resistivity $\eta = 6 \times 10^{-3}\,d_i v_A$ ($D_m = \frac{6\times 10^{-3}}{4\pi}$, `etaResistivity = 4.77465e6` m$^2$/s), `etaHyperCh = 0.0` (matching WarpX benchmark), `useUniformIonPressure = T` (stationary ion background, matching WarpX).
- **Boundaries**: Periodic in $x$, conducting walls for fields and reflecting walls for particles at $y = \pm y_{\max}$ (matching WarpX Dirichlet and reflecting conditions).
- **Timestepping & Subcycling**: $dt = 10^{-3}\,\tau_{ci} \approx 0.006283$ s, 40 magnetic subcycles (`nBSubcycle = 40`), $T_{\max} = 1.0\,\tau_{ci} \approx 6.283$ s (1,000 steps; full production reconnection runs up to $50\,\tau_{ci}$).

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
python3 tests/validate_tests.py --test=reconnection.fadeev_pic        # Fadeev full-PIC
python3 tests/validate_tests.py --test=reconnection.fadeev_hybrid     # Fadeev hybrid-PIC
python3 tests/validate_tests.py --test=reconnection.gem_pic           # Classic GEM challenge full-PIC
python3 tests/validate_tests.py --test=reconnection.gem_hybrid        # Classic GEM challenge hybrid-PIC
python3 tests/validate_tests.py --test=reconnection.asym_pic          # Asymmetric reconnection full-PIC
python3 tests/validate_tests.py --test=reconnection.forcefree_hybrid  # Force-free sheet hybrid-PIC
python3 tests/validate_tests.py --test=reconnection.forcefree_pic     # Force-free sheet full-PIC
```
*(Note: `reconnection.forcefree_hybrid` is an expensive benchmark and is skipped during default full-suite runs; run it by explicitly specifying `--test=reconnection.forcefree_hybrid` or adding `--include-expensive` / `--all`.)*

## Validation Checks (`validate.py`)

1. **Energy Log Sanity**:
   - Kinetic ion energy $E_{part}$ and magnetic energy $E_b$ are finite (no NaN/Inf).
   - $E_b$ remains bounded (no numerical blow-up).

2. **Equilibrium Initialization ($t=0$)**:
   - **Fadeev**: In-plane field nulls (O-points) located at $x \approx \pm \pi L \approx \pm 15.7\,d_i$; peak sheet density $\approx 1$, background $\approx 0.2$.
   - **GEM Challenge**: Harris sheet field reversal ($B_x \to \pm 1.0$ at top/bottom boundaries); central X-point null at $x \approx 0$; peak density in $(0.5, 1.0)$, background $< 0.4$.
   - **Asymmetric**: Central region field $B_x \approx +1.0$ ($B_1$), outer boundary field $B_x \approx -2.0$ ($-B_2$); sheet density enhancement in $(0.5, 1.6)$, background $< 0.4$.
   - **Force-Free Sheet**: Boundary asymptotic field $B_x \to \pm b_0$; midplane guide field $B_z(0) \approx \sqrt{b_g^2 + b_0^2}$; uniform total magnetic pressure ($|\mathbf{B}|^2 \approx b_0^2 + b_g^2$ with $<5\%$ relative variation); midplane anti-symmetry.

3. **Reconnection Dynamics**:
   - Seeded in-plane field perturbation grows nonlinearly ($\delta B_y$ increases).
   - Out-of-plane flux function $A_y$ at the X-point demonstrates active magnetic flux reconnection.

4. **Quasi-Neutrality**:
   - Initial condition verified to satisfy quasi-neutrality $|n_i - n_e| / n_0 < 0.5$ ($\rho_{S0} \approx 25\,\rho_{S1}$), verifying macroparticle charge scaling and electron mass loading.
