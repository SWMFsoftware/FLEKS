#!/usr/bin/env python3
"""Validator for the free-stream test (tests/freestream).

The runner exercises the single directory twice (full-PIC and hybrid Hall-OFF
variants), both validated with the shared hybrid-family energy-log checks.  No
dedicated plot-file check is applied.
"""
from .._shared.hybrid import validate_hybrid, validate_plot as _hyb_plot


def validate_log(pic_diags=None, test_name=None):
    """Free-stream energy-log validation (shared with the hybrid family)."""
    return validate_hybrid(pic_diags=pic_diags, test_name=test_name)


def validate_plot(test_name):
    """Free-stream has no dedicated plot check."""
    return _hyb_plot(test_name)
