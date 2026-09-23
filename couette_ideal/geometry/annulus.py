#!/usr/bin/env python3
"""
annulus.py -- structured hexahedral annulus for OpenFOAM, via the Gmsh
Python API.

Equivalent to annulus.geo, but callable as a function, which is what makes a
parameter sweep possible without editing files. Eccentricity is included as a
parameter because it is a two-line change in this topology: the inner circle
simply gets its own centre.

Topology: four quadrant surfaces between arcs, made transfinite and
recombined to quads, then extruded one layer in z to give hexahedra. The z
direction becomes an 'empty' patch pair in OpenFOAM, making the case 2D.

Install:
    conda install -c conda-forge gmsh python-gmsh

Single mesh:
    python annulus.py -o annulus.msh
    python annulus.py -o ecc.msh --eccentricity 0.0001
    python annulus.py -o annulus.msh --gui      # inspect before writing

Sweep (from your own driver script):
    from annulus import build_annulus
    for e in [0.0, 25e-6, 50e-6, 100e-6]:
        build_annulus(eccentricity=e, out_path=f"mesh_e{e*1e6:.0f}um.msh")
"""

import argparse
import sys

import gmsh


# ---------------------------------------------------------------------------
# Element type codes used by gmsh.model.mesh.getElements
# ---------------------------------------------------------------------------
ELEM_NAMES = {1: "line", 2: "triangle", 3: "quad", 4: "tet", 5: "hex", 6: "prism"}

# ---------------------------------------------------------------------------
# Rig geometry. Single source of truth: both the build_annulus signature and
# the CLI read these, so an import-driven sweep and a command-line run cannot
# silently describe different annuli.
# ---------------------------------------------------------------------------
DEF_RI = 0.010     # inner (rotating) radius, m
DEF_RO = 0.0105    # outer (stationary) radius, m
DEF_NR = 20        # radial cells across the gap
DEF_NA = 20        # azimuthal cells per quadrant
DEF_BUMP = 1.0


def build_annulus(
    Ri=DEF_RI,
    Ro=DEF_RO,
    t=None,
    nr=DEF_NR,
    na=DEF_NA,
    bump=DEF_BUMP,
    eccentricity=0.0,
    out_path="annulus.msh",
    gui=False,
    verbose=True,
):
    """Generate a structured hex annulus and write it as MSH 2 ASCII.

    Parameters
    ----------
    Ri, Ro : float
        Inner (rotating) and outer (stationary) radii, metres.
    t : float or None
        Axial thickness. Arbitrary for a 2D case -- the z direction becomes an
        'empty' patch. None (the default) sets it to the radial cell size,
        (Ro - Ri) / nr, which keeps the aspect ratio near unity so checkMesh
        stays quiet. Pass a number to override.
    nr, na : int
        Radial cells across the gap, and azimuthal cells per quadrant. Total
        cell count is nr * na * 4.
    bump : float
        Radial clustering. 1.0 is uniform; values below 1 cluster nodes toward
        both walls, the analogue of blockMesh's two-sided simpleGrading. Wall
        shear stress is evaluated across the first cell, so this is where
        resolution matters.
    eccentricity : float
        Offset of the inner cylinder axis along x, metres. Positive is +x,
        negative is -x; the mesh is mirror-symmetric in the sign, so only the
        magnitude matters physically. Zero gives the concentric case with a
        closed-form solution to validate against.
    out_path : str
        Output file. Written as MSH 2 ASCII -- gmshToFoam reads nothing else.

    Returns
    -------
    dict
        Mesh statistics.
    """
    # Check the radii first. Every quantity below divides by the gap, so a
    # non-positive gap has to be caught here or it surfaces later as a
    # confusing complaint about some other parameter.
    if Ro <= Ri:
        raise ValueError(
            f"outer radius Ro={Ro} must exceed inner radius Ri={Ri}"
        )

    gap = Ro - Ri

    if nr < 1 or na < 1:
        raise ValueError(f"nr and na must be >= 1; got nr={nr}, na={na}")

    # abs(): a negative eccentricity offsets the inner cylinder along -x, which
    # is just as capable of making the walls intersect as a positive one.
    if abs(eccentricity) >= gap:
        raise ValueError(
            f"|eccentricity| {abs(eccentricity) * 1e6:.1f} um reaches or "
            f"exceeds the gap {gap * 1e6:.1f} um; the cylinders would touch "
            "or intersect"
        )

    if t is None:
        t = gap / nr
    elif t <= 0.0:
        raise ValueError(f"thickness t must be positive; got {t}")

    gmsh.initialize()

    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)

        # gmshToFoam reads MSH 2 ASCII only. Gmsh defaults to 4.1.
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

        gmsh.model.add("annulus")
        geo = gmsh.model.geo

        e = eccentricity
        ae = abs(e)

        # -------------------------------------------------------------------
        # Points. Two centres: the outer cup on the axis, the inner cylinder
        # offset by e. For the concentric case they coincide, which is
        # harmless -- Gmsh tolerates duplicate-location points that are never
        # merged.
        # -------------------------------------------------------------------
        c_out = geo.addPoint(0.0, 0.0, 0.0)
        c_in = geo.addPoint(e, 0.0, 0.0)

        # Inner points, about the offset centre
        p_in = [
            geo.addPoint(e + Ri, 0.0, 0.0),
            geo.addPoint(e, Ri, 0.0),
            geo.addPoint(e - Ri, 0.0, 0.0),
            geo.addPoint(e, -Ri, 0.0),
        ]

        # Outer points, about the axis
        p_out = [
            geo.addPoint(Ro, 0.0, 0.0),
            geo.addPoint(0.0, Ro, 0.0),
            geo.addPoint(-Ro, 0.0, 0.0),
            geo.addPoint(0.0, -Ro, 0.0),
        ]

        # -------------------------------------------------------------------
        # Curves. addCircleArc takes (start, centre, end) and spans < 180 deg,
        # which is why each circle is four quadrant arcs.
        # -------------------------------------------------------------------
        arc_in = [
            geo.addCircleArc(p_in[i], c_in, p_in[(i + 1) % 4]) for i in range(4)
        ]
        arc_out = [
            geo.addCircleArc(p_out[i], c_out, p_out[(i + 1) % 4]) for i in range(4)
        ]
        radial = [geo.addLine(p_in[i], p_out[i]) for i in range(4)]

        # -------------------------------------------------------------------
        # Quadrant surfaces. Each loop runs: radial out, outer arc, radial in
        # (reversed), inner arc (reversed). This ordering fixes the position
        # of each wall in the extrude return list below -- change it and the
        # indices change with it.
        # -------------------------------------------------------------------
        surfaces = []
        for i in range(4):
            j = (i + 1) % 4
            loop = geo.addCurveLoop(
                [radial[i], arc_out[i], -radial[j], -arc_in[i]]
            )
            surfaces.append(geo.addPlaneSurface([loop]))

        # -------------------------------------------------------------------
        # Structured meshing. setTransfiniteCurve takes the number of NODES,
        # hence n + 1 for n cells.
        # -------------------------------------------------------------------
        for c in arc_in + arc_out:
            geo.mesh.setTransfiniteCurve(c, na + 1)

        for c in radial:
            if abs(bump - 1.0) < 1e-12:
                geo.mesh.setTransfiniteCurve(c, nr + 1)
            else:
                geo.mesh.setTransfiniteCurve(c, nr + 1, "Bump", bump)

        for s in surfaces:
            geo.mesh.setTransfiniteSurface(s)
            geo.mesh.setRecombine(2, s)        # triangles -> quads

        # -------------------------------------------------------------------
        # Extrude one layer in z. recombine=True turns the prisms into hexes.
        #
        # The return list has a fixed order per input surface:
        #   [0] top surface   [1] volume   [2:] sides, in curve-loop order
        #
        # With the loop ordering above, within each block of six:
        #   index 3 -> outer wall,  index 5 -> inner wall
        #
        # For geometries where the loop ordering varies, classify the returned
        # surfaces geometrically with gmsh.model.getBoundingBox instead of by
        # index -- more robust, and worth switching to if the topology grows.
        # -------------------------------------------------------------------
        ext = geo.extrude(
            [(2, s) for s in surfaces],
            0.0, 0.0, t,
            numElements=[1],
            recombine=True,
        )

        stride = 6
        if len(ext) != stride * len(surfaces):
            raise RuntimeError(
                f"unexpected extrude return length {len(ext)}; "
                "the curve loops may not all have four curves"
            )

        tops = [ext[stride * i + 0][1] for i in range(4)]
        vols = [ext[stride * i + 1][1] for i in range(4)]
        outers = [ext[stride * i + 3][1] for i in range(4)]
        inners = [ext[stride * i + 5][1] for i in range(4)]

        geo.synchronize()

        # -------------------------------------------------------------------
        # Physical groups. Only entities named here are exported. A boundary
        # surface omitted from a physical group is silently dropped, and the
        # solver then fails at startup complaining about a field boundary with
        # no matching patch -- an error that points nowhere near the cause.
        #
        # The names become the OpenFOAM patch names.
        # -------------------------------------------------------------------
        def named_group(dim, tags, name):
            tag = gmsh.model.addPhysicalGroup(dim, tags)
            gmsh.model.setPhysicalName(dim, tag, name)
            return tag

        named_group(2, inners, "innerWall")
        named_group(2, outers, "outerWall")
        named_group(2, tops + surfaces, "frontAndBack")
        named_group(3, vols, "internal")

        # -------------------------------------------------------------------
        # Mesh and report.
        # -------------------------------------------------------------------
        gmsh.model.mesh.generate(3)

        # etype, not t -- t is the thickness. Comprehension scope keeps them
        # separate in Python 3, but the collision is not worth reading past.
        etypes, etags, _ = gmsh.model.mesh.getElements(dim=3)
        counts = {
            ELEM_NAMES.get(etype, str(etype)): len(g)
            for etype, g in zip(etypes, etags)
        }

        expected = nr * na * 4
        stats = {
            "cells": counts,
            "expected_hexes": expected,
            "Ri": Ri,
            "Ro": Ro,
            "gap": gap,
            "thickness": t,
            "eccentricity": e,
            "eccentricity_ratio": ae / gap,
        }

        if verbose:
            print()
            print("=" * 60)
            print("MESH")
            print("=" * 60)
            for name, n in counts.items():
                print(f"  {name:<10s} {n:8d}")
            print(f"  expected hexes {expected:8d}")
            if counts.get("hex", 0) != expected or "tet" in counts:
                print("  WARNING: not a pure structured hex mesh")
            print(f"  gap        {gap * 1e3:8.3f} mm")
            print(f"  thickness  {t * 1e3:8.4f} mm  "
                  f"(z/radial aspect {t / (gap / nr):.2f})")
            if ae > 0:
                print(f"  eccentricity {e * 1e6:6.1f} um  "
                      f"(ratio {ae / gap:.3f})")
                h_min = gap - ae
                h_max = gap + ae
                print(f"  gap range  {h_min * 1e3:.4f} to {h_max * 1e3:.4f} mm")
                print(f"  stress range (narrow-gap est.) "
                      f"{-100 * ae / h_max:+.1f}% to "
                      f"{100 * ae / h_min:+.1f}%")
            print()

        if gui:
            gmsh.fltk.run()

        gmsh.write(out_path)
        if verbose:
            print(f"wrote {out_path}")

        return stats

    finally:
        gmsh.finalize()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="annulus.msh")
    ap.add_argument("--Ri", type=float, default=DEF_RI, help="inner radius, m")
    ap.add_argument("--Ro", type=float, default=DEF_RO, help="outer radius, m")
    ap.add_argument("-t", "--thickness", type=float, default=None,
                    help="axial thickness, m; default is the radial cell size")
    ap.add_argument("--nr", type=int, default=DEF_NR, help="radial cells")
    ap.add_argument("--na", type=int, default=DEF_NA,
                    help="azimuthal cells per quadrant")
    ap.add_argument("--bump", type=float, default=DEF_BUMP,
                    help="radial clustering; <1 clusters at both walls")
    ap.add_argument("-e", "--eccentricity", type=float, default=0.0, help="metres")
    ap.add_argument("--gui", action="store_true", help="open Gmsh before writing")
    args = ap.parse_args()

    build_annulus(
        Ri=args.Ri,
        Ro=args.Ro,
        t=args.thickness,
        nr=args.nr,
        na=args.na,
        bump=args.bump,
        eccentricity=args.eccentricity,
        out_path=args.output,
        gui=args.gui,
    )


if __name__ == "__main__":
    sys.exit(main())