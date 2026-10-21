# Hybrid PIC Upwind (Rusanov) Dissipation Test

Tests the upwind dissipation `#HYBRIDUPWIND` of the hybrid generalized Ohm's
law, which adds a Rusanov (local Lax-Friedrichs) term to the nodal electric
field:

$$
\mathbf{E} \mathrel{+}= \sum_d \frac{a_d}{2}\, \hat{\mathbf{e}}_d \times
\left(\mathbf{B}^R_d - \mathbf{B}^L_d\right),
$$

with $\mathbf{B}^{L,R}_d$ the MUSCL-reconstructed cell-centered field on the
two sides of the node along $d$ and $a_d$ the signal speed.

## Variants

| File | Setup | Check |
|---|---|---|
| `PARAM.in` | Frozen plasma (no particles), 32 cells, circularly polarized mode `waveMode = 8` on a guide field, `upwind1`, fixed speed $a = 0.5$ | Decay rate equals the analytic rate within 0.5 %, single exponential |
| `PARAM.in.mc` | Same with the MC limiter, `waveMode = 4` | Monotone decay, rate between the unlimited second-order and the first-order rate |
| `PARAM.in.divb` | 2D $32\times32$ periodic plasma with 400 ions/cell, Hall term, MC limiter, full whistler speed, no hyper-resistivity, no cleaning | $\max|\nabla\cdot\mathbf{B}(t_{\rm end}) - \nabla\cdot\mathbf{B}(0)| \le 10^{-10}\,\max|\mathbf{B}|/\Delta x$ |
| `PARAM.in.whistler` | The `whistler/PARAM.in.rk4` plasma (2000 ions/cell, Hall term) seeded with the $n = 4$ whistler ($k d_i = 3.9$, 16 cells per wavelength), MC limiter, default signal speed, no hyper-resistivity | Frequency within 4 % of the Hall branch; damping per radian $\gamma/\omega \in [0.0005, 0.005]$ |
| `PARAM.in.whistler_sensor` | Same with `upwindSensorFloor = 0.1` | Damping per radian $\gamma/\omega \in [0.00005, 0.0005]$ (measured 0.00015) |

## Analytic rates

The term adds $(a/2)\,J_{j-1/2}$ to the node between cells $j-1$ and $j$ (for
$E_z$ the jump of $B_y$, for $E_y$ minus the jump of $B_z$), and the
node-to-center curl gives, for each transverse component,

$$
\frac{dB_j}{dt} = \frac{a}{2\Delta x}\left(J_{j+1/2} - J_{j-1/2}\right).
$$

With $B_j \propto e^{ij\theta}$, $\theta = k\Delta x$:

* first order, $J_{j-1/2} = B_j - B_{j-1}$:
  $\gamma = \dfrac{2a}{\Delta x}\sin^2\dfrac{\theta}{2}$;
* unlimited second order (central slopes $s_j = (B_{j+1}-B_{j-1})/2$),
  $J_{j-1/2} = (B_j - s_j/2) - (B_{j-1} + s_{j-1}/2)$, which the MC limiter
  reduces to away from extrema. At the extrema of the mode the limiter falls
  back to first order, so the MC rate lies between the two.

The validator applies the RK4 amplification factor to the complex eigenvalue.
The guide field $B_x$ is normal to the jumps along $x$ and does not decay.

Measured: `upwind1` $\gamma = 2.5000$ (analytic 2.5000, residual
$7\times10^{-13}$); MC $\gamma = 0.187$ (bounds 0.107 and 0.732); div(B)
changes by $5\times10^{-15}\,\max|\mathbf{B}|/\Delta x$.

## Damping of a resolved whistler

`PARAM.in.whistler` guards the cost of the term at resolved scales. The seeded
whistler's complex amplitude $c(t)$ also carries the slow branch of the same
helicity; demodulating $c$ by the whistler frequency and averaging over one
period removes it, and the decay rate is fitted to that envelope. Added
damping per radian, measured at several modes of this deck ($dx = 0.1\,d_i$):

| $k d_i$ (cells per wavelength) | hyper-resistivity, `etaHyperCh` 0.01 | upwind, `upwindCoefWhistler` 0 | 0.25 | 1 |
|---|---|---|---|---|
| 2.0 (32) | 0.0002 | 0.0005 | 0.004 | 0.015 |
| 3.9 (16) | 0.0007 | 0.0022 | 0.020 | 0.070 |
| 7.9 (8) | 0.0026 | 0.0092 | 0.082 | 0.30 |

Both terms damp smooth modes as $k^4$. With the default signal speed (flow
plus fast wave) the upwind term acts like grid-mode hyper-resistivity with
`etaHyperCh` $\approx 0.03$; the whistler term raises that by 10–30 times,
which is why the default is 0. Neither term shifts the frequency.

`upwindSensorFloor` multiplies each component's jump by
$\phi = \min(1, \max(\text{floor}, s))$ with
$s = |\Delta^3 B| / \sum|\Delta B|$ over the same 4 cells. Any linear
correction on this stencil damps as $k^4$, so the sensor is what lets the
damping fall off faster toward resolved scales without a third ghost layer:

| $k d_i$ (cells per wavelength) | sensor floor 0.2 | 0.1 | 0.05 |
|---|---|---|---|
| 2.0 (32) | 0.00009 | 0.00005 | 0.00003 |
| 3.9 (16) | 0.0004 | 0.0002 | 0.00014 |
| 7.9 (8) | 0.002 | 0.002 | 0.002 |

At $k d_i \ge 8$ the sensor opens and the damping rises above
hyper-resistivity's. Grid noise much smaller than a large resolved wave sees
only the floor; a whistler with $B_1 = 0.1 B_0$ and 50 ions per cell stays as
quiet as with the unscaled term.

## div(B)

The term only modifies the nodal $\mathbf{E}$, and
$\mathrm{div}_{\rm node\to center} \circ \mathrm{avg}_{\rm center\to node}
\circ \mathrm{curl}_{\rm node\to center} \equiv 0$, so the update cannot change
the plotted `divB` (the divergence of the node-averaged field). This is unlike
the full-PIC `#UPWINDB`, which diffuses each B component directly and needs
hyperbolic cleaning.

## Running

```bash
python3 tests/validate_tests.py --test=hybrid_upwind
python3 tests/validate_tests.py --test=hybrid_upwind.divb -v
python3 tests/validate_tests.py --test=hybrid_upwind.whistler -v
python3 tests/validate_tests.py --test=hybrid_upwind.whistler_sensor -v
```
