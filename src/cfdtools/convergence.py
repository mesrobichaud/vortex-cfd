"""
convergence.py -- mesh-independence checks.

    from cfdtools.convergence import gci

    r = gci([T_coarse, T_medium, T_fine], [6600, 14850, 33500])
    r["p"], r["extrapolated"], r["gci_fine"]

gci follows the Grid Convergence Index procedure of Celik et al. (2008),
"Procedure for estimation and reporting of uncertainty due to discretization
in CFD applications", J. Fluids Eng. 130(7), 078001.
"""

import math


def gci(values, cells, dim=2, safety=1.25, tol=1e-10, max_iter=100):
    """Grid Convergence Index from one quantity on three meshes.

    Parameters
    ----------
    values : the quantity on each mesh, ordered coarse, medium, fine
    cells : cell count of each mesh, same order
    dim : 2 for a 2D mesh (one cell thick, e.g. the wedge), 3 for 3D.
        Sets the representative cell size h = (1 / cells) ** (1 / dim).
    safety : safety factor; 1.25 for three meshes (Celik et al.)

    Returns
    -------
    dict with
        p : observed order of accuracy
        extrapolated : estimate of the zero-cell-size value
        r21, r32 : refinement ratios, fine/medium and medium/coarse
        error_fine : relative change from medium to fine, |(f1 - f2) / f1|
        error_extrapolated : relative difference between fine and extrapolated
        gci_fine : relative uncertainty band on the fine-mesh value
        oscillatory : True if the values go up and down with refinement. The
            order and GCI are then unreliable; add a mesh or refine further.
    """
    f3, f2, f1 = (float(v) for v in values)
    n3, n2, n1 = (float(n) for n in cells)
    if not n3 < n2 < n1:
        raise ValueError("cells must increase from coarse to fine")

    h1, h2, h3 = ((1.0 / n) ** (1.0 / dim) for n in (n1, n2, n3))
    r21 = h2 / h1
    r32 = h3 / h2

    e21 = f2 - f1
    e32 = f3 - f2
    if e21 == 0.0:
        raise ValueError("medium and fine values are identical; the order can't be computed")
    ratio = e32 / e21
    s = math.copysign(1.0, ratio)

    # Observed order p, solved by fixed-point iteration (Celik et al., eq. 3).
    p = abs(math.log(abs(ratio))) / math.log(r21)
    for _ in range(max_iter):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if abs(p_new - p) < tol:
            p = p_new
            break
        p = p_new

    extrapolated = (r21 ** p * f1 - f2) / (r21 ** p - 1.0)
    error_fine = abs((f1 - f2) / f1)
    return {
        "p": p,
        "extrapolated": extrapolated,
        "r21": r21,
        "r32": r32,
        "error_fine": error_fine,
        "error_extrapolated": abs((extrapolated - f1) / extrapolated),
        "gci_fine": safety * error_fine / (r21 ** p - 1.0),
        "oscillatory": s < 0,
    }
