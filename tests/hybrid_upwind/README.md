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
```
