# Ion Acoustic Wave Standalone Test

This test validates the 1D self-consistent low-frequency electrostatic Ion Acoustic Wave (IAW) in a warm plasma with $T_e \gg T_i$.

## Background

In an Ion Acoustic Wave, the massive ions provide the inertia, while the hot, mobile electrons provide the restoring pressure force. The wave speed in the long-wavelength limit is given by:
$$ c_s = \sqrt{\frac{k_B T_e + \gamma_i k_B T_i}{m_i}} $$

In our test setup:
- Electron-to-ion mass ratio: $m_i / m_e = 100$ ($m_e = 0.01, m_i = 1.0$)
- Initial temperatures: $T_e = 1.0$, $T_i = 0.05$ (yielding $T_e / T_i = 20$)
- For $\gamma_i = 3$, this results in an ion acoustic speed of:
  $$ c_s = \sqrt{1.0 + 3 \times 0.05} = \sqrt{1.15} \approx 1.07 $$
- Wave vector: $k_x = 1.0$ (wavelength $\lambda = 2\pi$)
- Wave frequency: $\omega = k_x c_s \approx 1.07$
- Oscillation period: $T = 2\pi / \omega \approx 5.87$

A small initial density perturbation is applied to the electrons (Species 0):
$$ n_e(x, t=0) = 1.0 + 0.01 \cos(x) $$

This density perturbation drives a self-consistent initial electric field:
$$ E_x(x, t=0) = -4\pi \times 0.01 \sin(x) $$

The simulation runs for $T_{\text{max}} = 12.0$, resolving $\approx 2$ full periods of the ion acoustic wave.
