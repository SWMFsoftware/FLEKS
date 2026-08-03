#!/usr/bin/env python3
"""Validator for the hybrid-whistler test (tests/hybrid_whistler).

Uses the shared hybrid-family validator (energy-log checks) plus the seeded
wavelength / bounded-amplitude / whistler-dispersion plot check.
"""
from .._shared.hybrid import validate_hybrid, validate_plot as _hyb_plot


def validate_log(pic_diags=None, test_name=None):
    """Hybrid-wave energy-log validation (shared with the hybrid family)."""
    return validate_hybrid(pic_diags=pic_diags, test_name=test_name)


def validate_plot(test_name):
    """Hybrid-wave plot-output check (seeded wavelength + dispersion)."""
    return _hyb_plot(test_name)
