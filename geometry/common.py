#!/usr/bin/env python3
"""
common.py -- pieces shared by the mesh builders in this directory.

Holds the rig geometry for this case-group. The geometry-independent helpers
(element table, colour modes, sizing) live in cfdtools.mesh and are re-exported
here so the builders keep a single import.
"""

from cfdtools.mesh import (  # noqa: F401  (re-exported)
    COLOR_MODES, ELEM_NAMES,
    cells_for, color_mode, pick_size, resolve_count,
)


# ---------------------------------------------------------------------------
# Rig geometry. Single source of truth for every builder here, so an r-theta
# mesh and an r-z mesh cannot drift into describing different hardware.
# ---------------------------------------------------------------------------
DEF_RI = 0.0055     # inner (rotating bob) radius, m
DEF_RO = 0.00635    # outer (stationary cup) radius, m
