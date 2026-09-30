"""
paraview_view.py -- ParaView view of a 2D case, run inside ParaView's Python.

Do not run this file directly; use the launcher, which finds ParaView and
passes the case:

    python -m cfdtools.pv cases/smoke            # ParaView window
    python -m cfdtools.pv cases/smoke --save     # image, no window

Sets a white background, parallel projection, no orientation axes, flat
colours, the Cool to Warm colour map and a colour bar, and looks at
the mesh along its thin direction.

Options come from the command line (pvbatch) or from the PV_CASE and
PV_ARGS environment variables (ParaView window, which cannot take script
arguments). This file cannot import cfdtools: ParaView's Python does not
have it installed.
"""

import argparse
import glob
import os
import shlex
import sys

from paraview.simple import (
    ColorBy, GetActiveViewOrCreate, GetAnimationScene, GetColorTransferFunction,
    GetOpacityTransferFunction, GetScalarBar, OpenFOAMReader, Render,
    ResetCamera, SaveScreenshot, Show,
)

# ParaView colour map preset.
COLORMAP = "Cool to Warm"
INK = [11 / 255, 11 / 255, 11 / 255]

# Units shown in the colour bar title.
UNITS = {"U": "m/s", "p": "m²/s²"}


def parse_args():
    ap = argparse.ArgumentParser(description="ParaView view of an OpenFOAM case.")
    ap.add_argument("case", nargs="?", default=os.environ.get("PV_CASE"),
                    help="case folder or .foam file (default: $PV_CASE)")
    ap.add_argument("--field", default="U", help="cell field to colour by (default: U)")
    ap.add_argument("--time", default="latest", help="saved time (default: latest)")
    ap.add_argument("--range", nargs=2, type=float, metavar=("MIN", "MAX"),
                    help="fixed colour range (default: data range at that time)")
    ap.add_argument("-o", "--output", help="image file to write (.png)")
    ap.add_argument("--size", nargs=2, type=int, default=[3000, 2000],
                    metavar=("W", "H"), help="image size in pixels (default: 3000 2000)")
    argv = sys.argv[1:] + shlex.split(os.environ.get("PV_ARGS", ""))
    args, _ = ap.parse_known_args(argv)
    if not args.case:
        ap.error("give a case folder, or set PV_CASE")
    return args


def foam_file(case):
    """The .foam file for a case folder, created empty if it is missing."""
    case = os.path.abspath(case)
    if case.endswith(".foam"):
        return case
    found = sorted(glob.glob(os.path.join(case, "*.foam")))
    if found:
        return found[0]
    path = os.path.join(case, os.path.basename(case) + ".foam")
    open(path, "a").close()
    return path


def is_batch():
    """True under pvbatch/pvpython (no window), False in the ParaView GUI."""
    try:
        from paraview.modules.vtkRemotingCore import vtkProcessModule
        return vtkProcessModule.GetProcessType() != vtkProcessModule.PROCESS_CLIENT
    except Exception:
        return True


def main():
    args = parse_args()
    foam = foam_file(args.case)
    case_dir = os.path.dirname(foam)

    reader = OpenFOAMReader(FileName=foam)
    reader.MeshRegions = ["internalMesh"]
    times = list(reader.TimestepValues)
    t = times[-1] if args.time == "latest" else min(times, key=lambda v: abs(v - float(args.time)))

    view = GetActiveViewOrCreate("RenderView")
    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    scene.AnimationTime = t
    view.ViewTime = t

    display = Show(reader, view)
    reader.UpdatePipeline(t)

    # Colour by the field (magnitude for vectors).
    ColorBy(display, ("CELLS", args.field, "Magnitude"))
    lut = GetColorTransferFunction(args.field)
    if args.range:
        lo, hi = args.range
    else:
        lo, hi = reader.CellData[args.field].GetRange(-1)
        if reader.CellData[args.field].GetNumberOfComponents() > 1:
            lo = 0.0        # a magnitude starts at zero
    lut.ApplyPreset(COLORMAP, True)
    lut.RescaleTransferFunction(lo, hi)
    lut.AutomaticRescaleRangeMode = "Never"
    GetOpacityTransferFunction(args.field).RescaleTransferFunction(lo, hi)

    # Flat colours: no shading on the 2D surface.
    display.Ambient = 1.0
    display.Diffuse = 0.0

    # Colour bar.
    display.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    vector = reader.CellData[args.field].GetNumberOfComponents() > 1
    # Plain text: ParaView's math text treats "|" as a table separator.
    name = f"{args.field} magnitude" if vector else args.field
    unit = UNITS.get(args.field)
    bar.Title = f"{name} [{unit}]" if unit else name
    bar.ComponentTitle = ""
    bar.TitleColor = INK
    bar.LabelColor = INK
    bar.TitleFontSize = 16
    bar.LabelFontSize = 16
    bar.AutomaticLabelFormat = 0
    bar.LabelFormat = "{:.3g}"
    bar.RangeLabelFormat = "{:.3g}"
    bar.WindowLocation = "Any Location"
    bar.Orientation = "Vertical"
    bar.Position = [0.74, 0.12]
    bar.HorizontalTitle = 1          # title above the bar, not rotated
    bar.ScalarBarLength = 0.6

    # View: white background, no orientation axes, parallel projection.
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1
    try:
        view.UseColorPaletteForBackground = 0
        view.BackgroundColorMode = "Single Color"
    except AttributeError:
        pass
    view.Background = [1.0, 1.0, 1.0]

    # Look along the thinnest direction of the mesh (the 2D normal).
    b = reader.GetDataInformation().GetBounds()
    extent = [b[1] - b[0], b[3] - b[2], b[5] - b[4]]
    normal = extent.index(min(extent))
    centre = [0.5 * (b[2 * i] + b[2 * i + 1]) for i in range(3)]
    position = list(centre)
    position[normal] += max(extent)
    view.CameraFocalPoint = centre
    view.CameraPosition = position
    view.CameraViewUp = [0, 0, 1] if normal != 2 else [0, 1, 0]
    if is_batch():
        # Lay out at 1000 px wide; SaveScreenshot scales it (text included)
        # up to --size, so text stays in proportion at any resolution.
        view.ViewSize = [1000, round(1000 * args.size[1] / args.size[0])]
    ResetCamera(view)

    # Shift the domain left so the colour bar has its own space on the right.
    dop = [view.CameraFocalPoint[i] - view.CameraPosition[i] for i in range(3)]
    up = view.CameraViewUp
    right = [dop[1] * up[2] - dop[2] * up[1],
             dop[2] * up[0] - dop[0] * up[2],
             dop[0] * up[1] - dop[1] * up[0]]
    norm = sum(v * v for v in right) ** 0.5
    width, height = view.ViewSize
    scale = view.CameraParallelScale * 1.05
    view.CameraParallelScale = scale
    shift = 0.15 * 2.0 * scale * width / height
    view.CameraFocalPoint = [view.CameraFocalPoint[i] + shift * right[i] / norm for i in range(3)]
    view.CameraPosition = [view.CameraPosition[i] + shift * right[i] / norm for i in range(3)]

    Render(view)

    if args.output or is_batch():
        out = args.output or os.path.join(case_dir, "figures", f"paraview_{args.field}.png")
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        SaveScreenshot(out, view, ImageResolution=args.size)
        print(f"wrote {out}")


main()
