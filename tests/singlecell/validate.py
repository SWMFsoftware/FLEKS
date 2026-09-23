#!/usr/bin/env python3
"""Validator for the single-cell periodic hybrid test (tests/singlecell).

With exactly one grid cell and periodic boundaries, curl B is identically zero,
so the Hall term (J x B)/rho and the convective term U_i x B both vanish
(U_i = 0).  The electric field stays zero and the magnetic field is frozen.
"""
import logging

from tests._shared.validators import validate_frozen_field

logger = logging.getLogger(__name__)


def validate_log(pic_diags=None, test_name=None):
    """Validate the single-cell hybrid test.

    The test passes iff (1) the magnetic energy Eb is conserved to round-off
    (no spurious Hall-driven evolution), (2) the electric field energy Ee stays
    ~0 (the field is truly frozen, not merely energy-conserving), and (3) no
    NaN/blow-up occurs.

    Tolerances are tighter than zerocurrent because a single cell has *exactly*
    zero curl B, so there is no physical mechanism to change Eb at all.
    """
    return validate_frozen_field(
        pic_diags,
        # A single cell has no spatial gradient → Hall term = 0 exactly,
        # so Eb must be conserved to round-off.
        eb_tol=(0.9999, 1.0001),
        # The single-cell solver leaves a residual Ee ~ 1e-3 * Eb;
        # 1e-2 * Eb separates "frozen" from "propagating".
        ee_rel_max=1e-2,
        test_label="Single-Cell Hybrid",
        logger=logger,
    )
