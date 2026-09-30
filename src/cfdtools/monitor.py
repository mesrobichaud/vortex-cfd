"""
monitor.py -- live plots of an OpenFOAM case's postProcessing/*.dat files.

A stand-in for OpenFOAM's foamMonitor, which needs gnuplot's x11 terminal and
GNU coreutils and so does not run on macOS as shipped. Run it in a second
terminal while the solver runs:

    python -m cfdtools.monitor cases/smoke
    python -m cfdtools.monitor cases/smoke --refresh 2
    python -m cfdtools.monitor cases/smoke --save monitor.png   # one shot, no window

Every function object that writes a .dat file gets its own panel. solverInfo
files are shown as initial residuals on a log scale; everything else as each
column against time. Columns that are only numerical noise next to the largest
one in their file (e.g. the x/y torque of a 2D case) are left out. Only the
latest run is shown: files left by earlier runs of the same case are ignored.
"""

import argparse
import math
import sys
from pathlib import Path


# A column is dropped when its largest magnitude is below this fraction of the
# largest magnitude in the same file.
NOISE_RATIO = 1e-6


def _as_float(text):
    try:
        return float(text)
    except ValueError:
        return None


def _series_name(stem):
    """'moment_0' -> 'moment'. Re-running a case does not overwrite a
    function object's file; OpenFOAM writes <name>_<startTime>.dat beside it."""
    base, sep, tail = stem.rpartition("_")
    return base if sep and _as_float(tail) is not None else stem


def find_series(case):
    """Map 'object/file' -> the .dat paths of the current run, in time order.

    The current run is the most recently written file of each series. Earlier
    start-time directories are kept only if they started before it (a restart);
    anything else is left over from a previous run and ignored.
    """
    root = Path(case) / "postProcessing"
    if not root.is_dir():
        return {}

    found = {}
    for obj_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        for tdir in obj_dir.iterdir():
            start = _as_float(tdir.name)
            if start is None or not tdir.is_dir():
                continue
            for f in tdir.glob("*.dat"):
                key = f"{obj_dir.name}/{_series_name(f.stem)}"
                found.setdefault(key, []).append((start, f.stat().st_mtime, f))

    series = {}
    for key, files in found.items():
        current_start, _, current = max(files, key=lambda x: x[1])
        # Newest file in each earlier start directory, then the current one.
        earlier = {}
        for start, mtime, f in files:
            if start < current_start and mtime > earlier.get(start, (0.0, None))[0]:
                earlier[start] = (mtime, f)
        series[key] = [earlier[s][1] for s in sorted(earlier)] + [current]
    return series


def read_dat(paths):
    """Read OpenFOAM .dat files into {column name: [floats]}.

    Column names come from the last '#' line before the data. Files are joined
    in order; where a later file restarts at an earlier time, it replaces the
    overlap. Non-numeric columns (solver names, converged flags) are dropped.
    """
    names = None
    rows = []
    for path in paths:
        try:
            lines = path.read_text().splitlines()
        except OSError:
            continue
        segment = []
        for line in lines:
            if line.startswith("#"):
                tokens = line[1:].split()
                if tokens:
                    names = tokens
                continue
            tokens = line.replace("(", " ").replace(")", " ").split()
            # Skip a line the solver is still writing.
            if names and len(tokens) == len(names):
                segment.append(tokens)
        if segment:
            t0 = _as_float(segment[0][0])
            if t0 is not None:
                rows = [r for r in rows if float(r[0]) < t0]
            rows.extend(segment)

    if not names or not rows:
        return {}

    columns = {}
    for i, name in enumerate(names):
        try:
            columns[name] = [float(r[i]) for r in rows]
        except ValueError:
            continue
    return columns


def select(key, columns):
    """Pick the columns worth plotting. Returns (time, {name: values}, log_y)."""
    time = columns.pop("Time", None)
    if time is None:
        return None, {}, False

    is_residuals = any(n.endswith("_initial") for n in columns)
    if is_residuals:
        columns = {n: v for n, v in columns.items() if n.endswith("_initial")}

    peak = {n: max(abs(x) for x in v) for n, v in columns.items()}
    top = max(peak.values(), default=0.0)
    if top > 0.0:
        columns = {n: v for n, v in columns.items() if peak[n] >= NOISE_RATIO * top}
    return time, columns, is_residuals


def draw(fig, case, ncols=1):
    """Redraw every panel from the files on disk. Returns the panel count."""
    fig.clf()
    panels = []
    for key, paths in find_series(case).items():
        time, columns, log_y = select(key, read_dat(paths))
        if time and columns:
            panels.append((key, time, columns, log_y))

    if not panels:
        fig.text(0.5, 0.5, f"waiting for {Path(case) / 'postProcessing'} ...",
                 ha="center", va="center")
        return 0

    ncols = max(1, min(ncols, len(panels)))
    nrows = math.ceil(len(panels) / ncols)
    for i, (key, time, columns, log_y) in enumerate(panels, start=1):
        ax = fig.add_subplot(nrows, ncols, i)
        for name, values in columns.items():
            ax.plot(time, values, label=name, linewidth=1.2)
        if log_y:
            ax.set_yscale("log")
        ax.set_title(key)
        ax.set_xlabel("Time [s]")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize="small")
    fig.suptitle(str(Path(case).resolve().name))
    fig.tight_layout()
    return len(panels)


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m cfdtools.monitor",
        description="Live plots of an OpenFOAM case's postProcessing/*.dat files.",
    )
    ap.add_argument("case", nargs="?", default=".", help="case directory (default: .)")
    ap.add_argument("-r", "--refresh", type=float, default=1.0,
                    help="seconds between refreshes (default: 1)")
    ap.add_argument("-c", "--columns", type=int, default=1,
                    help="panels per row (default: 1)")
    ap.add_argument("--save", metavar="FILE",
                    help="draw once, save to FILE and exit (no window)")
    args = ap.parse_args(argv)

    if not Path(args.case).is_dir():
        ap.error(f"no such case directory: {args.case}")

    import matplotlib
    if args.save:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("dark_background")
    fig = plt.figure(figsize=(5, 7))

    if args.save:
        if draw(fig, args.case, args.columns) == 0:
            print(f"nothing to plot in {args.case}/postProcessing", file=sys.stderr)
            return 1
        fig.savefig(args.save, dpi=120)
        print(f"wrote {args.save}")
        return 0

    # Redrawing from scratch each time resets any zoom; pause the refresh by
    # closing and reopening if you need to inspect a region.
    draw(fig, args.case, args.columns)
    plt.show(block=False)
    while plt.fignum_exists(fig.number):
        # Not plt.pause(): it calls show() every time, which raises the window
        # above everything else on each refresh.
        fig.canvas.start_event_loop(args.refresh)
        if plt.fignum_exists(fig.number):
            draw(fig, args.case, args.columns)
            fig.canvas.draw_idle()
    return 0


if __name__ == "__main__":
    sys.exit(main())
