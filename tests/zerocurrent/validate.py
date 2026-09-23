#!/usr/bin/env python3
"""Validator for the zero-current hybrid wave test (tests/zerocurrent).

No macroparticles are loaded (rho = 0 everywhere), so every cell is left inert:
the generalized Ohm's law is fully short-circuited by the ``if (rho > 0)`` guard
in assemble_ohm_E, so the electric field E = 0.  Faraday's law then gives
dB/dt = -curl E = 0: the seeded sinusoidal B perturbation is FROZEN and does NOT
propagate.
"""
import logging

from tests._shared.validators import validate_frozen_field

logger = logging.getLogger(__name__)


def validate_log(pic_diags=None, test_name=None):
    """Validate the zero-current hybrid wave test.

    The test passes iff (1) the magnetic energy Eb is conserved (the wave
    neither grows, decays, nor travels) AND (2) the electric field energy Ee
    stays ~0 -- the genuine signature that the field is frozen rather than
    merely energy-conserving (a non-dispersive propagating wave would also
    conserve Eb but would have Ee > 0).

    Tolerances are looser than singlecell (which has exactly zero curl B)
    because zerocurrent uses a full spatial grid with a small seeded perturbation
    where minor rounding differences accumulate.
    """
    return validate_frozen_field(
        pic_diags,
        eb_tol=(0.999, 1.001),
        # A propagating wave has Ee > 0 even while conserving Eb; require
        # Ee to be essentially machine-zero relative to Eb.
        ee_rel_max=1e-6,
        test_label="Zero-Current Hybrid",
        logger=logger,
    )
