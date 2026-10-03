#!/usr/bin/env python3
"""
wedge.py -- axisymmetric r-z section of a concentric-cylinder (bob-in-cup)
cell for OpenFOAM, via the Gmsh Python API.

Companion to annulus.py. That file slices the cell in r-theta and assumes the
flow is fully developed; this one slices it in r-z, which is what you need as
soon as the bottom gap contributes to the measured torque.

Geometry, with bottom_gap > 0 -- the L-shaped domain:

    r=0      Ri        Ro
     |        |         |
     |        +=========+  z = hb + H     top
     |        |         |
     |        | annulus |     rotorSide <-|-> cupWall
     |        |         |
     +========+=========+  z = hb         rotorBottom
     |                  |
     |    bottom gap    |                 |-> cupWall
     +==================+  z = 0          cupBottom
     ^
    axis

With bottom_gap = 0 the foot disappears and the domain is the plain annular
rectangle Ri <= r <= Ro, 0 <= z <= H.

Bob shape (bottom_gap > 0 only):
  bottom_angle  slope of the bob's bottom face from horizontal, degrees. The
                face is a cone z = hb + r*tan(angle): positive puts the tip
                at the axis pointing down, negative gives a recess, 0 is flat.
                bottom_gap (hb) is the clearance at the axis.
  edge_radius   fillet between the bob's side and bottom face, tangent to
                both. 0 keeps the sharp corner. The fillet surface belongs to
                the rotorBottom patch.
The annulus height H is measured from the sharp-corner height
z_edge = hb + Ri*tan(angle), so changing the fillet does not move the top.

    With a fillet, a fourth block D sits between the upper half of the arc
    and the cup wall. The arc is split at its midpoint M; the lower half is
    part of block A:

            |          |  C      |
            |       T1 +---------+  z = z_t1
            |          )    D    |
            |        M +---------+  z = z_m
            +------T2  |         |
            |   A      |    B    |
            +----------+---------+  z = 0

The cross-section is drawn in the x-z plane (x = r), rotated by -angle/2 about
z, then revolved through +angle. That leaves the two end planes symmetric about
y = 0, which is what OpenFOAM's 'wedge' patch type requires.

Where the domain reaches the axis the wedge collapses to a line, so the cells
in that first column are prisms rather than hexes, and there are no faces on
the axis at all. Both are expected: the 'axis' patch in an OpenFOAM wedge case
is type 'empty' and carries zero faces.

Install:
    conda install -c conda-forge gmsh python-gmsh

Single mesh:
    python wedge.py -o wedge.msh
    python wedge.py -o nogap.msh --bottom-gap 0
    python wedge.py -o wedge.msh --gui       # inspect before writing

Refinement. Counts are explicit by default; a target cell size derives them
instead (rounded up, so the achieved size is at or below what you ask for):
    python wedge.py -o fine.msh  --max-size 50e-6
    python wedge.py -o mixed.msh --max-size-r 250e-6 --max-size-z 100e-6
    python wedge.py -o hand.msh  --max-size-z 100e-6 --nr-gap 40
An explicit --nr-gap / --nz-ann / ... always wins over a size. The core is 20x
the gap, so a single --max-size spends most of its cells under the bob; the
mixed form above keeps the core coarse and resolves the gap by hand. The theta
direction is exempt either way: a wedge is one cell thick by definition.

Sweep (from your own driver script):
    from wedge import build_wedge
    for hb in [0.0, 0.5e-3, 1.0e-3, 2.0e-3]:
        build_wedge(bottom_gap=hb, out_path=f"rz_hb{hb*1e3:.1f}mm.msh")

gmshToFoam writes every patch as type 'patch'. Edit constant/polyMesh/boundary
afterwards and set front and back to 'wedge', or the case is not axisymmetric.
"""

import argparse
import math
import sys

import gmsh

from common import (
    COLOR_MODES, DEF_RI, DEF_RO, ELEM_NAMES,
    color_mode, pick_size, resolve_count,
)


# ---------------------------------------------------------------------------
# Defaults specific to the r-z section. The counts are fallbacks, used when
# neither an explicit count nor a target cell size is given.
# ---------------------------------------------------------------------------
DEF_H = 0.015          # immersed height of the annular section, m
DEF_HB = 0.001         # bottom gap; 0 disables the L and gives a rectangle, m
DEF_ANGLE = 5.0        # wedge opening angle, degrees
DEF_NR_GAP = 20        # radial cells across the annular gap
DEF_NR_CORE = 100      # radial cells from the axis to Ri, bottom gap only
DEF_NZ_ANN = 150       # axial cells up the annulus
DEF_NZ_BOT = 30        # axial cells across the bottom gap
DEF_EDGE_RADIUS = 0.0  # fillet radius at the bob's bottom edge; 0 = sharp, m
DEF_BOTTOM_ANGLE = 0.0 # slope of the bob's bottom face; 0 = flat, degrees


def build_wedge(
    Ri=DEF_RI,
    Ro=DEF_RO,
    H=DEF_H,
    bottom_gap=DEF_HB,
    r_axis=0.0,
    angle=DEF_ANGLE,
    edge_radius=DEF_EDGE_RADIUS,
    bottom_angle=DEF_BOTTOM_ANGLE,
    nr_gap=None,
    nr_core=None,
    nz_ann=None,
    nz_bot=None,
    n_edge=None,
    max_size=None,
    max_size_r=None,
    max_size_z=None,
    bump_r=1.0,
    bump_z=1.0,
    out_path="wedge.msh",
    gui=False,
    color_by="physical",
    verbose=True,
):
    """Generate an axisymmetric wedge of the r-z section, as MSH 2 ASCII.

    Parameters
    ----------
    Ri, Ro : float
        Inner (rotating bob) and outer (stationary cup) radii, metres.
    H : float
        Height of the annular section, from the bob's sharp-corner height
        z_edge = bottom_gap + Ri*tan(bottom_angle) to the top of the fluid,
        metres.
    bottom_gap : float
        Clearance between the bob's bottom face and the cup floor at the axis,
        metres. Zero removes the foot entirely and leaves a plain annular
        rectangle, which reproduces the fully developed case and is the
        baseline to difference the end effect against.
    edge_radius : float
        Fillet radius between the bob's side and bottom face, metres. 0 keeps
        the sharp corner. Requires bottom_gap > 0.
    bottom_angle : float
        Slope of the bob's bottom face from horizontal, degrees. Positive is a
        cone with its tip at the axis pointing down; negative is a recess.
        Requires bottom_gap > 0.
    r_axis : float
        Inner radius of the bottom-gap block. Zero (the default) is the
        physical flat-bottomed bob: the fluid reaches the axis and the wedge
        collapses there. A positive value leaves a solid core on the axis, so
        every cell stays a hex at the cost of a wall that is not real.
    angle : float
        Wedge opening angle in degrees. 5 is the usual OpenFOAM choice; the
        sector is built as +/- angle/2 so the end planes straddle y = 0.
    nr_gap, nr_core : int or None
        Radial cells across Ri..Ro, and across r_axis..Ri. The second is only
        used when bottom_gap > 0. None falls back to max_size, then to the
        DEF_* constant.
    nz_ann, nz_bot : int or None
        Axial cells up the annulus, and across the bottom gap.
    n_edge : int or None
        Cells along the fillet arc (edge_radius > 0 only). None matches the
        radial cell size in the gap.
    max_size : float or None
        Target upper bound on cell size, metres, in every direction at once.
        Counts are rounded up, so the achieved size is at or below this. Note
        the theta direction is exempt: a wedge is one cell thick by definition,
        so its cell size is set by `angle` and Ro alone.
    max_size_r, max_size_z : float or None
        Per-direction overrides of max_size, radial and axial. Radial covers
        both the gap and the core, axial both the annulus and the bottom gap.
        The core is ~20x the gap here, so a single global max_size spends most
        of its cells under the bob; set max_size_r coarse and refine the gap
        with an explicit nr_gap when that matters.
    bump_r : float
        Radial clustering across the annular gap; < 1 clusters toward both
        walls. This is where the wall shear stress is evaluated.
    bump_z : float
        Axial clustering across the bottom gap; < 1 clusters toward the cup
        floor and the bob's bottom face.
    out_path : str
        Output file. Written as MSH 2 ASCII -- gmshToFoam reads nothing else.

    Returns
    -------
    dict
        Mesh statistics.
    """
    if Ro <= Ri:
        raise ValueError(f"outer radius Ro={Ro} must exceed inner radius Ri={Ri}")
    if H <= 0.0:
        raise ValueError(f"annulus height H must be positive; got {H}")
    if bottom_gap < 0.0:
        raise ValueError(f"bottom_gap must be >= 0; got {bottom_gap}")
    if not 0.0 <= r_axis < Ri:
        raise ValueError(f"r_axis must satisfy 0 <= r_axis < Ri={Ri}; got {r_axis}")
    if not 0.0 < angle < 180.0:
        raise ValueError(f"wedge angle must be in (0, 180) degrees; got {angle}")
    if edge_radius < 0.0:
        raise ValueError(f"edge_radius must be >= 0; got {edge_radius}")
    if not -45.0 < bottom_angle < 45.0:
        raise ValueError(f"bottom_angle must be in (-45, 45) degrees; got {bottom_angle}")
    if bottom_gap == 0.0 and (edge_radius > 0.0 or bottom_angle != 0.0):
        raise ValueError("edge_radius and bottom_angle need a bottom gap (bottom_gap > 0)")

    gap = Ro - Ri
    has_bottom = bottom_gap > 0.0

    size_r = pick_size(max_size_r, max_size)
    size_z = pick_size(max_size_z, max_size)

    nr_gap, src_gap = resolve_count("nr_gap", nr_gap, gap, size_r, DEF_NR_GAP)
    nz_ann, src_ann = resolve_count("nz_ann", nz_ann, H, size_z, DEF_NZ_ANN)
    if has_bottom:
        nr_core, src_core = resolve_count(
            "nr_core", nr_core, Ri - r_axis, size_r, DEF_NR_CORE)
        nz_bot, src_bot = resolve_count(
            "nz_bot", nz_bot, bottom_gap, size_z, DEF_NZ_BOT)
    else:
        # Blocks A and B are never built, so these counts go unused. Zero
        # rather than a stale default keeps the cell-count arithmetic honest.
        nr_core = nz_bot = 0
        src_core = src_bot = "unused"

    half = math.radians(angle) / 2.0
    hb = bottom_gap if has_bottom else 0.0

    # Bob bottom: the cone z = hb + r*tan(slope), meeting the side r = Ri at
    # the sharp-corner height z_edge. The fillet (radius rf, centre (rc, zc))
    # replaces that corner. It touches the side at T1 = (Ri, z_t1) and the
    # bottom at T2 = (r_t2, z_t2). With rf = 0 both are the sharp corner.
    slope = math.radians(bottom_angle)
    tan_s, sin_s, cos_s = math.tan(slope), math.sin(slope), math.cos(slope)
    rf = edge_radius
    z_edge = hb + Ri * tan_s
    rc = Ri - rf
    zc = z_edge + rf * (1.0 - sin_s) / cos_s
    r_t2 = rc + rf * sin_s
    z_t2 = zc - rf * cos_s
    z_t1 = zc
    z_top = z_edge + H
    # Midpoint M of the arc: halfway between the directions from the centre
    # to T2 (slope - 90 deg) and to T1 (0 deg).
    theta_m = (slope - math.pi / 2.0) / 2.0
    r_m = rc + rf * math.cos(theta_m)
    z_m = zc + rf * math.sin(theta_m)

    if has_bottom:
        if r_t2 <= r_axis:
            raise ValueError(f"edge_radius={rf} is too large: the fillet reaches r_axis")
        # Lowest point of the bob: the bottom of the fillet circle when it lies
        # on the arc (recess), otherwise T2 (the sharp corner when rf = 0).
        z_low = min(hb + r_axis * tan_s, zc - rf if slope < 0.0 else z_t2)
        if z_low <= 0.0:
            raise ValueError("the bob reaches the cup floor: increase bottom_gap "
                             "or reduce bottom_angle / edge_radius")
        if z_t1 >= z_top:
            raise ValueError(f"edge_radius={rf} is taller than the annulus height H={H}")

    if has_bottom and rf > 0.0:
        arc_len = rf * (math.pi / 2.0 - slope)
        if n_edge is None:
            n_edge = max(2, math.ceil(arc_len / (gap / nr_gap)))
            src_edge = "gap cell size"
        else:
            if n_edge < 2:
                raise ValueError(f"n_edge must be >= 2; got {n_edge}")
            src_edge = "explicit"
        # Lower half of the arc (in block A) and upper half (block D).
        n_lo = n_edge // 2
        n_hi = n_edge - n_lo
    else:
        arc_len = 0.0
        n_edge = n_lo = n_hi = 0
        src_edge = "unused"

    gmsh.initialize()

    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)

        # gmshToFoam reads MSH 2 ASCII only. Gmsh defaults to 4.1.
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

        gmsh.model.add("wedge")
        geo = gmsh.model.geo

        def pt(r, z):
            return geo.addPoint(r, 0.0, z)

        # -------------------------------------------------------------------
        # Cross-section in the x-z plane, x standing in for r.
        #
        # Block C is the annulus and always exists. Blocks A and B are the two
        # halves of the foot and appear only when bottom_gap > 0: A from the
        # axis out to r_t2 under the bob's bottom face, B from r_t2 out to Ro.
        # Block D sits between the fillet arc and the cup wall and appears
        # only when edge_radius > 0. Splitting the domain this way is what
        # lets every block stay a four-sided transfinite patch.
        # -------------------------------------------------------------------
        t1 = pt(Ri, z_t1)       # bob side meets the fillet (the corner if rf = 0)
        q1 = pt(Ro, z_t1)
        p7 = pt(Ri, z_top)
        p8 = pt(Ro, z_top)

        l_int_C = geo.addLine(t1, q1)       # z = z_t1, Ri..Ro
        l_right_C = geo.addLine(q1, p8)     # r = Ro, z_t1..z_top
        l_top_C = geo.addLine(p8, p7)       # z = z_top, Ro..Ri
        l_left_C = geo.addLine(p7, t1)      # r = Ri, z_top..z_t1

        loop_C = geo.addCurveLoop([l_int_C, l_right_C, l_top_C, l_left_C])
        surfaces = [geo.addPlaneSurface([loop_C])]

        radial_gap = [l_int_C, l_top_C]
        axial_ann = [l_right_C, l_left_C]
        radial_core = []
        axial_bot = []
        edge_lo = []
        edge_hi = []
        corners = {}        # surface -> corner points, for blocks with > 4 curves

        if has_bottom:
            p1 = pt(r_axis, 0.0)
            p3 = pt(Ro, 0.0)
            p4 = pt(r_axis, hb + r_axis * tan_s)

            if rf > 0.0:
                # The arc is split at its midpoint M. The lower half continues
                # the bob's bottom face as the top of block A (no kink at T2);
                # the upper half is the left side of block D. The B/D line
                # starts at M, where the arc is at ~45 degrees, so no cell
                # corner is sharper than ~45 degrees.
                t2 = pt(r_t2, z_t2)             # fillet meets the bob's bottom face
                m = pt(r_m, z_m)                # midpoint of the arc
                q2 = pt(Ro, z_m)
                centre = pt(rc, zc)
                f_t2 = pt(r_t2, 0.0)
                f2 = pt(r_m, 0.0)

                l_arc_hi = geo.addCircleArc(t1, centre, m)
                l_arc_lo = geo.addCircleArc(m, centre, t2)
                l_int_BD = geo.addLine(m, q2)       # z = z_m, r_m..Ro
                l_right_D = geo.addLine(q2, q1)     # r = Ro, z_m..z_t1

                loop_D = geo.addCurveLoop([l_int_BD, l_right_D, -l_int_C, l_arc_hi])
                surfaces.append(geo.addPlaneSurface([loop_D]))

                l_bot_A = geo.addLine(p1, f_t2)     # z = 0, r_axis..r_t2
                l_bot_A2 = geo.addLine(f_t2, f2)    # z = 0, r_t2..r_m, under the lower arc
                l_bot_B = geo.addLine(f2, p3)       # z = 0, r_m..Ro
                l_right_B = geo.addLine(p3, q2)     # r = Ro, 0..z_m
                l_int_AB = geo.addLine(f2, m)       # r = r_m, 0..z_m
                l_bob = geo.addLine(t2, p4)         # bob's bottom face, r_t2..r_axis
                l_axis = geo.addLine(p4, p1)        # r = r_axis, down to 0

                loop_A = geo.addCurveLoop(
                    [l_bot_A, l_bot_A2, l_int_AB, l_arc_lo, l_bob, l_axis])
                loop_B = geo.addCurveLoop([l_bot_B, l_right_B, -l_int_BD, -l_int_AB])
                s_A = geo.addPlaneSurface([loop_A])
                surfaces += [s_A, geo.addPlaneSurface([loop_B])]
                # A has six curves; its four corners pair the floor (two
                # curves) with the bottom face plus the lower arc.
                corners[s_A] = [p1, f2, m, p4]

                radial_gap += [l_int_BD, l_bot_B]
                radial_core += [l_bot_A, l_bob]
                axial_bot += [l_right_B, l_axis, l_int_AB]
                edge_hi += [l_arc_hi, l_right_D]
                edge_lo += [l_arc_lo, l_bot_A2]
            else:
                f2 = pt(Ri, 0.0)

                l_bot_A = geo.addLine(p1, f2)       # z = 0, r_axis..Ri
                l_bot_B = geo.addLine(f2, p3)       # z = 0, Ri..Ro
                l_right_B = geo.addLine(p3, q1)     # r = Ro, 0..z_edge
                l_bob = geo.addLine(t1, p4)         # bob's bottom face, Ri..r_axis
                l_axis = geo.addLine(p4, p1)        # r = r_axis, down to 0
                l_int_AB = geo.addLine(f2, t1)      # r = Ri, 0..z_edge

                loop_A = geo.addCurveLoop([l_bot_A, l_int_AB, l_bob, l_axis])
                loop_B = geo.addCurveLoop([l_bot_B, l_right_B, -l_int_C, -l_int_AB])
                surfaces += [geo.addPlaneSurface([loop_A]),
                             geo.addPlaneSurface([loop_B])]

                radial_gap.append(l_bot_B)
                radial_core += [l_bot_A, l_bob]
                axial_bot += [l_right_B, l_axis, l_int_AB]

        # -------------------------------------------------------------------
        # Structured meshing. setTransfiniteCurve takes NODES, hence n + 1.
        # Shared edges get their count from one list only, so opposite sides of
        # every block match by construction:
        #   radial_gap  -> nr_gap    (l_int_C and l_int_BD are shared)
        #   axial_bot   -> nz_bot    (l_int_AB is shared by A and B)
        #   edge_lo     -> n_lo      (lower arc half and the floor under it)
        #   edge_hi     -> n_hi      (upper arc half and the cup wall facing it)
        # -------------------------------------------------------------------
        def transfinite(curves, n, bump):
            for c in curves:
                # Bump is symmetric about the curve midpoint, so a reversed
                # curve in a loop needs no special handling here.
                if abs(bump - 1.0) < 1e-12:
                    geo.mesh.setTransfiniteCurve(c, n + 1)
                else:
                    geo.mesh.setTransfiniteCurve(c, n + 1, "Bump", bump)

        transfinite(radial_gap, nr_gap, bump_r)
        transfinite(axial_ann, nz_ann, 1.0)
        transfinite(radial_core, nr_core, 1.0)
        transfinite(axial_bot, nz_bot, bump_z)
        transfinite(edge_lo, n_lo, 1.0)
        transfinite(edge_hi, n_hi, 1.0)

        for s in surfaces:
            geo.mesh.setTransfiniteSurface(s, cornerTags=corners.get(s, []))
            geo.mesh.setRecombine(2, s)        # triangles -> quads

        # -------------------------------------------------------------------
        # Swing the section to -angle/2, then revolve through +angle. The two
        # end planes land at -angle/2 and +angle/2, straddling y = 0, which is
        # the symmetry OpenFOAM's wedge transform assumes.
        #
        # All blocks go through one revolve call: shared curves are then swept
        # once and the blocks stay conformal across their interfaces. Separate
        # calls would silently produce a split mesh.
        # -------------------------------------------------------------------
        dimtags = [(2, s) for s in surfaces]
        geo.rotate(dimtags, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -half)
        geo.revolve(
            dimtags,
            0.0, 0.0, 0.0,          # a point on the axis
            0.0, 0.0, 1.0,          # the axis direction: z
            2.0 * half,
            numElements=[1],
            recombine=True,
        )

        geo.synchronize()

        # -------------------------------------------------------------------
        # Classify the boundary surfaces geometrically rather than by their
        # position in the revolve return list. With three blocks and a
        # collapsing edge on the axis, index arithmetic is not worth trusting.
        #
        # A surface is on the boundary iff exactly one volume is adjacent to
        # it, which drops the two block interfaces without naming them.
        # -------------------------------------------------------------------
        scale = max(Ro, z_top)
        tol = 1e-6 * scale
        # A wall straddles y = 0; a wedge plane sits wholly on one side of it.
        # The margin has to stay well under the y-halfwidth of the *innermost*
        # wall, which is the core wall when there is one, not Ro.
        r_small = r_axis if r_axis > 0.0 else Ri
        tol_y = 1e-3 * r_small * math.sin(half)

        def rz_points(tag):
            """Corner points of a surface as (r, z, y).

            Not getBoundingBox: on a revolved surface that box is built from a
            tessellation and sits ~1e-7 m inside the true radius, which is far
            too coarse to tell Ri from Ro at a 0.5 mm gap. The corner points
            are exact.
            """
            points = []
            for _, pnt in gmsh.model.getBoundary(
                [(2, tag)], combined=False, oriented=False, recursive=True
            ):
                x, y, z = gmsh.model.getValue(0, pnt, [])
                points.append((math.hypot(x, y), z, y))
            return points

        def on_bob_bottom(r, z):
            """True if (r, z) lies on the bob's bottom face or its fillet."""
            on_face = (r_axis - tol <= r <= r_t2 + tol
                       and abs(z - (hb + r * tan_s)) <= tol)
            on_fillet = (rf > 0.0 and r >= r_t2 - tol and z <= zc + tol
                         and abs(math.hypot(r - rc, z - zc) - rf) <= tol)
            return on_face or on_fillet

        groups = {name: [] for name in (
            "front", "back", "rotorSide", "rotorBottom",
            "cupWall", "cupBottom", "top", "coreWall",
        )}
        collapsed = []
        unclassified = []

        for _, tag in gmsh.model.getEntities(2):
            up, _ = gmsh.model.getAdjacencies(2, tag)
            if len(up) != 1:
                continue                        # interior block interface

            points = rz_points(tag)
            rs = [p[0] for p in points]
            zs = [p[1] for p in points]
            ys = [p[2] for p in points]
            box = (min(rs), max(rs), min(zs), max(zs))

            def every(test):
                return all(test(r, z) for r, z, _ in points)

            # The sliver left behind where the wedge closes on the axis. It has
            # no area, so it must not become a patch: an OpenFOAM wedge simply
            # has no faces there.
            if max(rs) <= tol:
                collapsed.append(tag)
                continue

            if min(ys) >= -tol_y:               # everything at +angle/2
                groups["front"].append(tag)
            elif max(ys) <= tol_y:              # everything at -angle/2
                groups["back"].append(tag)
            elif every(lambda r, z: abs(r - Ro) <= tol):
                groups["cupWall"].append(tag)
            elif every(lambda r, z: abs(r - Ri) <= tol):
                groups["rotorSide"].append(tag)
            elif r_axis > 0.0 and every(lambda r, z: abs(r - r_axis) <= tol):
                groups["coreWall"].append(tag)
            elif every(lambda r, z: abs(z) <= tol):
                groups["cupBottom"].append(tag)
            elif every(lambda r, z: abs(z - z_top) <= tol):
                groups["top"].append(tag)
            elif every(on_bob_bottom):          # bottom face and fillet
                groups["rotorBottom"].append(tag)
            else:
                unclassified.append((tag, *box))

        if unclassified:
            raise RuntimeError(
                "could not classify boundary surfaces "
                f"{[u[0] for u in unclassified]}; (rmin, rmax, zmin, zmax) = "
                f"{[tuple(round(v, 9) for v in u[1:]) for u in unclassified]}"
            )
        for name in ("front", "back", "cupWall", "cupBottom", "top", "rotorSide"):
            if not groups[name]:
                raise RuntimeError(f"no boundary surface was classified as '{name}'")

        # -------------------------------------------------------------------
        # Physical groups. Only entities named here are exported; a boundary
        # surface left out is dropped silently, and the solver then fails at
        # startup about a field boundary with no matching patch.
        #
        # rotorSide and rotorBottom stay separate so their torque can be
        # integrated independently -- which is the whole reason for meshing the
        # r-z section. Merge them in the boundary file if you want one wall.
        # -------------------------------------------------------------------
        def named_group(dim, tags, name):
            tag = gmsh.model.addPhysicalGroup(dim, tags)
            gmsh.model.setPhysicalName(dim, tag, name)
            return tag

        for name, tags in groups.items():
            if tags:
                named_group(2, tags, name)

        vols = [tag for _, tag in gmsh.model.getEntities(3)]
        named_group(3, vols, "internal")

        # -------------------------------------------------------------------
        # Mesh and report.
        # -------------------------------------------------------------------
        gmsh.model.mesh.generate(3)

        etypes, etags, enodes = gmsh.model.mesh.getElements(dim=3)
        counts = {
            ELEM_NAMES.get(etype, str(etype)): len(g)
            for etype, g in zip(etypes, etags)
        }

        # A prism on the axis is fine. A hex carrying a repeated node is a
        # collapsed cell that gmshToFoam cannot make sense of, so count them.
        degenerate = 0
        for etype, nodes in zip(etypes, enodes):
            nn = gmsh.model.mesh.getElementProperties(etype)[3]
            for k in range(len(nodes) // nn):
                chunk = nodes[k * nn:(k + 1) * nn]
                if len(set(chunk)) != nn:
                    degenerate += 1

        cells_ann = nr_gap * nz_ann
        cells_bot = (nr_core + n_lo + nr_gap) * nz_bot if has_bottom else 0
        cells_edge = nr_gap * n_hi
        expected = cells_ann + cells_bot + cells_edge
        # The column against the axis is what the collapse turns into prisms.
        expected_prisms = nz_bot if (has_bottom and r_axis == 0.0) else 0

        stats = {
            "cells": counts,
            "expected_cells": expected,
            "expected_prisms": expected_prisms,
            "degenerate_cells": degenerate,
            "Ri": Ri,
            "Ro": Ro,
            "gap": gap,
            "H": H,
            "bottom_gap": hb,
            "edge_radius": rf,
            "bottom_angle_deg": bottom_angle,
            "z_edge": z_edge,
            "z_top": z_top,
            "wedge_angle_deg": angle,
            "patches": {k: len(v) for k, v in groups.items() if v},
            "collapsed_surfaces": len(collapsed),
        }

        if verbose:
            total = sum(counts.values())
            print()
            print("=" * 62)
            print("MESH  (axisymmetric r-z wedge)")
            print("=" * 62)
            for name, n in sorted(counts.items()):
                print(f"  {name:<10s} {n:8d}")
            print(f"  {'total':<10s} {total:8d}   expected {expected}")
            if total != expected:
                print("  WARNING: cell count does not match the block layout")
            if counts.get("prism", 0) != expected_prisms:
                print(f"  WARNING: {counts.get('prism', 0)} prisms, "
                      f"expected {expected_prisms} on the axis")
            if degenerate:
                print(f"  WARNING: {degenerate} cells carry a repeated node; "
                      "gmshToFoam will not accept these")
            print()
            def line(label, length, n, src):
                print(f"  {label:<10s} {length * 1e3:9.4f} mm   "
                      f"{n:5d} cells  {length / n * 1e6:8.1f} um  [{src}]")

            line("gap", gap, nr_gap, src_gap)
            line("annulus H", H, nz_ann, src_ann)
            if has_bottom:
                line("bottom gap", hb, nz_bot, src_bot)
                line("core", r_t2 - r_axis, nr_core, src_core)
                if rf > 0.0:
                    line("edge arc", arc_len, n_edge, src_edge)
                print(f"  bob bottom {bottom_angle:9.3f} deg slope, "
                      f"edge radius {rf * 1e3:.4f} mm")
            # theta is not a refinement direction: a wedge is one cell thick,
            # so this size follows from the opening angle alone.
            print(f"  {'theta':<10s} {'':12s}     1 cells  "
                  f"{Ro * math.radians(angle) * 1e6:8.1f} um  [angle, at Ro]")
            if has_bottom:
                # Rough end-effect scale, both over a common pi*mu*Omega:
                #   side torque  = 2 pi mu Omega Ri^3 H / gap   (Couette)
                #   base torque  =   pi mu Omega Ri^4 / (2 hb)  (torsional)
                # so the ratio collapses to Ri*gap / (4 hb H). It says how much
                # torque the foot is worth before you have solved anything.
                side = 2.0 * Ri ** 3 * H / gap
                base = Ri ** 4 / (2.0 * hb)
                print(f"  bottom/side torque (rough)  {base / side:6.3f}"
                      f"   [Ri*gap/(4*hb*H)]")
            else:
                print("  bottom gap      none   (plain annular rectangle)")
            print(f"  wedge      {angle:9.3f} deg")
            reaches_axis = has_bottom and r_axis == 0.0
            if reaches_axis:
                print("  axis       wedge collapses; that column is prisms "
                      "and carries no faces")
            if collapsed:
                print(f"             {len(collapsed)} zero-area surface(s) "
                      "excluded from every patch")
            print()
            print("  patches:")
            for name, tags in sorted(groups.items()):
                if tags:
                    print(f"    {name:<12s} {len(tags)} surface(s)")
            print()
            print("  after gmshToFoam, in constant/polyMesh/boundary:")
            print("    front, back -> type wedge;")
            if reaches_axis:
                print("    axis is optional: no faces lie on it. Add")
                print("      axis { type empty; nFaces 0; startFace <nFaces>; }")
                print("      only if a BC file names it.")
            print()

        if gui:
            # Display only -- see COLOR_MODES. Colouring by physical group is
            # what makes the patch assignment visible, which is the reason to
            # open the GUI at all.
            gmsh.option.setNumber("Mesh.ColorCarousel", color_mode(color_by))
            gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
            gmsh.fltk.run()

        gmsh.write(out_path)
        if verbose:
            print(f"wrote {out_path}")

        return stats

    finally:
        gmsh.finalize()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("-o", "--output", default="wedge.msh")
    ap.add_argument("--Ri", type=float, default=DEF_RI, help="bob radius, m")
    ap.add_argument("--Ro", type=float, default=DEF_RO, help="cup radius, m")
    ap.add_argument("-H", "--height", type=float, default=DEF_H,
                    help="annular section height, m")
    ap.add_argument("--bottom-gap", type=float, default=DEF_HB,
                    help="bob-to-floor clearance, m; 0 gives a plain rectangle")
    ap.add_argument("--r-axis", type=float, default=0.0,
                    help="inner radius of the foot, m; 0 reaches the axis")
    ap.add_argument("--angle", type=float, default=DEF_ANGLE,
                    help="wedge opening angle, degrees")
    ap.add_argument("--edge-radius", type=float, default=DEF_EDGE_RADIUS,
                    help="fillet radius at the bob's bottom edge, m; 0 is sharp")
    ap.add_argument("--bottom-angle", type=float, default=DEF_BOTTOM_ANGLE,
                    help="slope of the bob's bottom face, degrees; "
                         "positive = tip down at the axis, 0 = flat")
    ap.add_argument("--nr-gap", type=int, default=None,
                    help=f"radial cells across the gap (default {DEF_NR_GAP})")
    ap.add_argument("--nr-core", type=int, default=None,
                    help=f"radial cells axis to Ri (default {DEF_NR_CORE})")
    ap.add_argument("--nz-ann", type=int, default=None,
                    help=f"axial cells up the annulus (default {DEF_NZ_ANN})")
    ap.add_argument("--nz-bot", type=int, default=None,
                    help=f"axial cells across the foot (default {DEF_NZ_BOT})")
    ap.add_argument("--n-edge", type=int, default=None,
                    help="cells along the fillet arc (default: gap cell size)")
    ap.add_argument("--max-size", type=float, default=None,
                    help="target max cell size in r and z, m")
    ap.add_argument("--max-size-r", type=float, default=None,
                    help="target max radial cell size, m; overrides --max-size")
    ap.add_argument("--max-size-z", type=float, default=None,
                    help="target max axial cell size, m; overrides --max-size")
    ap.add_argument("--bump-r", type=float, default=1.0,
                    help="radial clustering in the gap; <1 clusters at walls")
    ap.add_argument("--bump-z", type=float, default=1.0,
                    help="axial clustering in the bottom gap; <1 at the walls")
    ap.add_argument("--gui", action="store_true", help="open Gmsh before writing")
    ap.add_argument("--color-by", default="physical",
                    choices=sorted(COLOR_MODES),
                    help="GUI mesh colouring; display only, never written")
    args = ap.parse_args()

    build_wedge(
        Ri=args.Ri,
        Ro=args.Ro,
        H=args.height,
        bottom_gap=args.bottom_gap,
        r_axis=args.r_axis,
        angle=args.angle,
        edge_radius=args.edge_radius,
        bottom_angle=args.bottom_angle,
        nr_gap=args.nr_gap,
        nr_core=args.nr_core,
        nz_ann=args.nz_ann,
        nz_bot=args.nz_bot,
        n_edge=args.n_edge,
        max_size=args.max_size,
        max_size_r=args.max_size_r,
        max_size_z=args.max_size_z,
        bump_r=args.bump_r,
        bump_z=args.bump_z,
        out_path=args.output,
        gui=args.gui,
        color_by=args.color_by,
    )


if __name__ == "__main__":
    sys.exit(main())
