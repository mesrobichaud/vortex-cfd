#!/usr/bin/env python3
"""
couette_figure.py -- velocity figure for an annulus-2D case.

    python post/couette_figure.py cases/smoke
    python post/couette_figure.py cases/smoke --time 0.1 -o figures/smoke.pdf

(a) Velocity magnitude on the mesh cells, with the annulus unrolled:
    angle about the cup axis against radius.
(b) Tangential velocity across the gap, u_theta / (Omega Ri) against the
    normalised gap position xi = (distance from bob) / (local gap width).
    Concentric case: cell-centre values against the analytical Couette
    solution. Eccentric case: profiles along the x axis at the narrowest and
    widest gap (linear interpolation between cells).

Default output: <case>/figures/couette.pdf and .png.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cfdtools import post


def couette_profile(xi, Ri, Ro):
    """Analytical u_theta / (Omega Ri) for concentric cylinders, inner rotating."""
    r = Ri + xi * (Ro - Ri)
    return (Ri / r) * (Ro**2 - r**2) / (Ro**2 - Ri**2)


def unroll(poly):
    """One cell's (x, y) vertices -> (angle in degrees, radius in mm)."""
    theta = np.degrees(np.unwrap(np.arctan2(poly[:, 1], poly[:, 0])))
    # Cells on the +-180 degree line: keep each cell on one side.
    if theta.mean() > 180.0:
        theta -= 360.0
    elif theta.mean() < -180.0:
        theta += 360.0
    return np.column_stack([theta, np.hypot(poly[:, 0], poly[:, 1]) * 1e3])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("case", help="case directory")
    ap.add_argument("--time", default="latest", help="saved time to plot (default: latest)")
    ap.add_argument("-o", "--output", help="output file (.pdf, .png, .svg); default <case>/figures/couette.pdf")
    args = ap.parse_args(argv)

    case = Path(args.case)
    p = post.case_parameters(case)
    Ri, Ro, ecc, omega = float(p["Ri"]), float(p["Ro"]), float(p["ecc"]), float(p["omega"])
    mesh, t = post.read_case(case, args.time)
    z = 0.5 * (mesh.bounds[4] + mesh.bounds[5])
    u_ref = omega * Ri

    post.use_style()
    fig, (ax_map, ax_prof) = plt.subplots(
        1, 2, figsize=(7.0, 3.0), layout="constrained",
        gridspec_kw={"width_ratios": [1.3, 1]})

    # (a) |U| on the cells, annulus unrolled about the cup axis.
    pc = post.plot_field(ax_map, mesh, "U", vmin=0.0, transform=unroll)
    ax_map.set_xlim(-180.0, 180.0)
    ax_map.set_ylim((Ri - ecc) * 1e3, Ro * 1e3)
    ax_map.set_xticks([-180, -90, 0, 90, 180])
    ax_map.set_xlabel(r"$\theta$ [deg]")
    ax_map.set_ylabel("r [mm]")
    cbar = fig.colorbar(pc, ax=ax_map, pad=0.02)
    cbar.set_label("|U| [m/s]")
    cbar.outline.set_linewidth(0.4)
    ax_map.set_title("(a)", loc="left")

    # (b) Profiles across the gap.
    if ecc == 0.0:
        c = post.cell_centres(mesh)
        U = mesh.cell_data["U"]
        u_theta = (-c[:, 1] * U[:, 0] + c[:, 0] * U[:, 1]) / np.hypot(c[:, 0], c[:, 1])
        # Radial position of each cell: midway between its inner and outer
        # nodes, which lie on the arcs. The centroid of a flat-sided cell sits
        # inside the arc and would shift every point towards the bob.
        r = np.array([0.5 * (rv.min() + rv.max()) for rv in
                      (np.hypot(v[:, 0], v[:, 1]) for v in post.cell_vertices(mesh))])
        # Every cell at the same radius has the same value; average per radius.
        r_key = np.round(r, 9)
        radii = np.unique(r_key)
        u_mean = np.array([u_theta[r_key == rk].mean() for rk in radii])
        xi = (radii - Ri) / (Ro - Ri)

        xi_line = np.linspace(0.0, 1.0, 200)
        ax_prof.plot(xi_line, couette_profile(xi_line, Ri, Ro),
                     color=post.INK, linewidth=1.0, linestyle="--", label="Analytical")
        ax_prof.plot(xi, u_mean / u_ref, linestyle="none", marker="o",
                     markersize=4, color=post.SERIES[0], label="OpenFOAM (cell centres)")
    else:
        # Narrow gap on +x, wide gap on -x (the bob is offset along +x).
        for sign, label, colour in [(+1, "Narrowest gap", post.SERIES[0]),
                                    (-1, "Widest gap", post.SERIES[1])]:
            start = (ecc + sign * Ri, 0.0, z)
            end = (sign * Ro, 0.0, z)
            gap = abs(end[0] - start[0])
            d, uy = post.sample_line(mesh, start, end, "U", which="y")
            ax_prof.plot(d / gap, sign * uy / u_ref, color=colour, label=label)

    ax_prof.set_xlim(0.0, 1.0)
    # Allow negative values: the wide gap of an eccentric annulus can have
    # reverse flow near the cup.
    ax_prof.set_ylim(min(0.0, ax_prof.get_ylim()[0]), 1.05)
    ax_prof.axhline(0.0, color=post.INK_MUTED, linewidth=0.6, zorder=1)
    ax_prof.set_xlabel(r"$\xi$ (0 = bob, 1 = cup)")
    ax_prof.set_ylabel(r"$u_\theta / (\Omega R_i)$")
    ax_prof.legend(loc="upper right")
    ax_prof.set_title("(b)", loc="left")

    fig.suptitle(f"{case.resolve().name}   t = {t:g} s   e = {ecc * 1e3:g} mm",
                 color=post.INK_SECONDARY, fontsize=8)

    out = Path(args.output) if args.output else case / "figures" / "couette.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    if out.suffix == ".pdf":
        fig.savefig(out.with_suffix(".png"))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
