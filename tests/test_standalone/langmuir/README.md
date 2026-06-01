# Langmuir Wave Electrostatic Oscillation Test

This test validates the 1D self-consistent electrostatic Langmuir wave oscillation in a cold, periodic plasma.

## Background

A small initial density perturbation is applied only to electrons (Species 0) using the parameter `waveAmp` (e.g. `0.1`):
$$ n_e(x, t=0) = 1.0 + \text{waveAmp} \cos(x) $$

Under the Gaussian units solved by the PIC engine ($\epsilon_0 = 1 / 4\pi$), this density perturbation drives a self-consistent initial electric field:
$$ E_x(x, t=0) = -4\pi \text{waveAmp} \sin(x) $$

The plasma will oscillate at the electron plasma frequency:
$$ \omega_{pe} = \sqrt{\frac{4\pi n_0 q_e^2}{m_e}} \approx 35.45 $$
for $n_0 = 1.0$, $q_e = -1.0$, $m_e = 0.01$. This yields an oscillation period of:
$$ T = \frac{2\pi}{\omega_{pe}} \approx 0.177 $$

The total simulation time is $T_{\text{max}} = 0.9$, resolving $\approx 5$ full periods of oscillation.
