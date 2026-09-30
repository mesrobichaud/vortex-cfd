"""
post.py -- read OpenFOAM results and draw publication figures with matplotlib.

    from cfdtools import post

    mesh, t = post.read_case("cases/smoke")             # latest saved time
    params = post.case_parameters("cases/smoke")        # values the case was built with

    post.use_style()
    fig, ax = plt.subplots()
    pc = post.plot_field(ax, mesh, "U")                 # |U| on the cells of a 2D mesh
    fig.colorbar(pc, ax=ax, label="|U| [m/s]")

    d, values = post.sample_line(mesh, (0.010, 0, z), (0.0105, 0, z), "U")

plot_field draws each mesh cell as a polygon, so the figure shows the real
cells with no interpolation. Holes in the domain stay empty.
"""

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Colours. Sequential ramp: one hue, light to dark, for magnitudes. Series:
# fixed order, used for line plots. Ink: text, axes and reference lines.
# ---------------------------------------------------------------------------
SEQUENTIAL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
              "#256abf", "#184f95", "#104281", "#0d366b"]
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"


def sequential_cmap(name="cfd_blue"):
    """Single-hue sequential colormap (light = low, dark = high)."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, SEQUENTIAL)


def use_style(font_size=9):
    """matplotlib settings for print figures: white background, thin lines,
    muted axes, vector-friendly fonts."""
    import matplotlib as mpl
    from cycler import cycler

    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "savefig.dpi": 300,
        "pdf.fonttype": 42,              # embed TrueType, keeps text editable
        "font.size": font_size,
        "axes.titlesize": font_size,
        "axes.labelsize": font_size,
        "legend.fontsize": font_size - 1,
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK_MUTED,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "xtick.color": INK_SECONDARY,
        "ytick.color": INK_SECONDARY,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth": 1.5,
        "lines.markersize": 4,
        "legend.frameon": False,
        "axes.prop_cycle": cycler(color=SERIES),
    })


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def _foam_file(case):
    """The case's .foam file, created empty if missing (ParaView/VTK need one)."""
    case = Path(case)
    found = sorted(case.glob("*.foam"))
    if found:
        return found[0]
    path = case / f"{case.resolve().name}.foam"
    path.touch()
    return path


def read_case(case, time="latest"):
    """Read the internal mesh and cell fields at one saved time.

    Parameters
    ----------
    case : path to the case directory
    time : "latest", or a saved time value (the nearest saved time is used)

    Returns
    -------
    (mesh, time) : pyvista.UnstructuredGrid with fields as cell data, and the
    time actually read.
    """
    import pyvista as pv

    reader = pv.OpenFOAMReader(str(_foam_file(case)))
    reader.cell_to_point_creation = False       # keep raw cell values only
    times = reader.time_values
    if not times:
        raise ValueError(f"no saved times in {case}")
    t = times[-1] if time == "latest" else min(times, key=lambda v: abs(v - float(time)))
    reader.set_active_time_value(t)
    mesh = reader.read()["internalMesh"]
    return mesh, t


def case_parameters(case):
    """Parameter values the case was built with (PyFoamPrepareCaseParameters),
    as a dict. Includes the derived SI values from derivedParameters.py."""
    from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

    path = Path(case) / "PyFoamPrepareCaseParameters"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found; was the case built with pyFoamPrepareCase.py?")
    return dict(ParsedParameterFile(str(path), noHeader=True).getValueDict())


def component(values, which=None):
    """Reduce a cell array to one value per cell.

    which : None or "mag" (magnitude of a vector), 0/1/2 or "x"/"y"/"z".
    Scalars are returned unchanged.
    """
    values = np.asarray(values)
    if values.ndim == 1:
        return values
    if which in (None, "mag"):
        return np.linalg.norm(values, axis=1)
    index = {"x": 0, "y": 1, "z": 2}.get(which, which)
    return values[:, int(index)]


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_field(ax, mesh, field, which=None, normal="z", cmap=None,
               vmin=None, vmax=None, transform=None, rasterized=True):
    """Draw a field on the cells of a 2D (one-cell-thick) mesh.

    The mesh is cut through the middle along `normal`; each cut cell is drawn
    as a filled polygon coloured by its cell value.

    Parameters
    ----------
    ax : matplotlib Axes
    mesh : mesh from read_case
    field : cell field name, e.g. "U" or "p"
    which : component for vectors, see component()
    normal : "x", "y" or "z", the direction the 2D mesh is one cell thick in
    cmap : colormap; default is the single-hue sequential_cmap()
    transform : optional function mapping one cell's (n, 2) vertex array to
        new 2D coordinates, e.g. to unroll an annulus into (angle, radius).
        Without it the axes use equal scaling.
    rasterized : draw the cells as an image inside an otherwise vector PDF,
        which keeps large meshes small

    Returns
    -------
    matplotlib.collections.PolyCollection, for fig.colorbar(...)
    """
    from matplotlib.collections import PolyCollection

    axis = {"x": 0, "y": 1, "z": 2}[normal]
    lo, hi = mesh.bounds[2 * axis], mesh.bounds[2 * axis + 1]
    origin = list(mesh.center)
    origin[axis] = 0.5 * (lo + hi)
    cut = mesh.slice(normal=normal, origin=origin)

    keep = [i for i in range(3) if i != axis]
    xy = cut.points[:, keep]
    polygons = [xy[face] for face in cut.irregular_faces]
    if transform is not None:
        polygons = [transform(poly) for poly in polygons]
    values = component(cut.cell_data[field], which)

    pc = PolyCollection(polygons, array=values, cmap=cmap or sequential_cmap(),
                        edgecolors="face", linewidths=0.2, rasterized=rasterized)
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    if transform is None:
        ax.set_aspect("equal")
    ax.autoscale_view()
    ax.grid(False)
    return pc


def cell_centres(mesh):
    """Cell centre coordinates, shape (n_cells, 3)."""
    return mesh.cell_centers().points


def cell_vertices(mesh):
    """Vertex coordinates of each cell, as a list of (n, 3) arrays."""
    offsets = mesh.cell_offsets
    conn = mesh.cell_connectivity
    return [mesh.points[conn[a:b]] for a, b in zip(offsets[:-1], offsets[1:])]


def sample_line(mesh, start, end, field, which=None, n=200):
    """Values of a field along a straight line, linearly interpolated.

    Returns (distance from start, values). Points outside the mesh are dropped.
    """
    line = mesh.cell_data_to_point_data().sample_over_line(start, end, resolution=n)
    valid = line.point_data["vtkValidPointMask"].astype(bool)
    distance = np.linalg.norm(line.points - np.asarray(start, dtype=float), axis=1)
    values = component(line.point_data[field], which)
    return distance[valid], values[valid]
