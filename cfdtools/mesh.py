"""
mesh.py -- geometry-independent helpers for the Gmsh mesh builders.

Holds the Gmsh element-type table, the GUI colouring modes, and the sizing
helpers that turn a target cell size into a transfinite node count.

Transfinite meshes ignore Gmsh's own characteristic-length machinery
(Mesh.CharacteristicLengthMax, mesh.setSize, background fields) -- a
transfinite curve takes a node count and nothing else. Size-driven refinement
therefore has to happen before the geometry is built: pick the counts from the
target size, then hand those to Gmsh.
"""

import math


# ---------------------------------------------------------------------------
# Element type codes used by gmsh.model.mesh.getElements
# ---------------------------------------------------------------------------
ELEM_NAMES = {1: "line", 2: "triangle", 3: "quad", 4: "tet", 5: "hex", 6: "prism"}


# ---------------------------------------------------------------------------
# Gmsh GUI mesh-colouring modes (the Mesh.ColorCarousel option).
#
# Display only. The MSH 2 writer emits $PhysicalNames, $Nodes and $Elements and
# nothing else -- there is no colour field in the format -- so whichever mode is
# set, gmshToFoam reads an identical file. Safe to change freely.
# ---------------------------------------------------------------------------
COLOR_MODES = {
    "type": 0,        # by element type: hexes vs prisms
    "entity": 1,      # by elementary entity: gmsh's own default, one per tag
    "physical": 2,    # by physical group: the OpenFOAM patch names
    "partition": 3,   # by mesh partition
}


def color_mode(name):
    """Mesh.ColorCarousel value for a mode name."""
    if name not in COLOR_MODES:
        raise ValueError(
            f"color_by must be one of {sorted(COLOR_MODES)}; got {name!r}"
        )
    return COLOR_MODES[name]


def cells_for(length, max_size, minimum=1):
    """Cells needed to cover `length` with none larger than `max_size`.

    The epsilon keeps a length that divides exactly from tipping into a spare
    cell on floating-point noise: 0.5 mm at 25 um should be 20 cells, not 21.
    """
    if length <= 0.0:
        raise ValueError(f"length must be positive; got {length}")
    if max_size <= 0.0:
        raise ValueError(f"max cell size must be positive; got {max_size}")
    return max(minimum, math.ceil(length / max_size - 1e-9))


def resolve_count(name, explicit, length, max_size, default):
    """Pick the cell count for one edge from whichever control was supplied.

    Precedence is explicit count, then target cell size, then the built-in
    default. Returns (count, source) so the caller can report which control
    actually decided it, which matters when a global --max-size silently
    overrides a default someone was relying on.
    """
    if explicit is not None:
        if explicit < 1:
            raise ValueError(f"{name} must be >= 1; got {explicit}")
        return int(explicit), "explicit"
    if max_size is not None:
        return cells_for(length, max_size), f"<={max_size * 1e6:g}um"
    return int(default), "default"


def pick_size(specific, fallback):
    """Per-direction target size, falling back to the global one."""
    return specific if specific is not None else fallback
