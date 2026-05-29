# Langmuir Wave Electrostatic Oscillation Test

This test validates the 1D self-consistent electrostatic Langmuir wave oscillation in a cold, periodic plasma.

## Background

A small initial density perturbation is applied only to electrons (Species 0):
$$ n_e(x, t=0) = 1.0 + 0.02 \cos(x) $$

This density perturbation creates a charge imbalance that drives a self-consistent initial electric field:
$$ E_x(x, t=0) = -0.02 \sin(x) $$

The plasma will oscillate at the electron plasma frequency:
$$ \omega_{pe} = \sqrt{\frac{n_0 e^2}{m_e \epsilon_0}} = 10.0 $$
with oscillation period $T = 0.6283$.

The total simulation time is $T_{\text{max}} = 100$, resolving $\approx 159$ periods of oscillation.
