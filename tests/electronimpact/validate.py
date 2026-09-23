#!/usr/bin/env python3
"""Validator for the electron-impact ionization test (tests/electronimpact).

Checks that the heaviest ion species (O+, which receives the exosphere source)
energy increases over time, confirming the electron-impact ionization source is
active.  Uses the PIC energy log (log_pic_n*.log) as the data source.
"""
import logging

from tests._shared.validators import validate_ionisation_source

logger = logging.getLogger(__name__)


def validate_log(pic_diags=None, test_name=None):
    """Validate that the electron-impact ionization source is active.

    Delegates to the shared ``validate_ionisation_source`` helper which checks
    that the heaviest ion species (last EpartN column, O+) energy increased
    over time.
    """
    return validate_ionisation_source(
        pic_diags,
        test_label="Electron-Impact Ionization",
        logger=logger,
    )
