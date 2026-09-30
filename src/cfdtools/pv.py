"""
pv.py -- open a case in ParaView with the publication view set up, or
save that view as an image without opening a window.

    python -m cfdtools.pv                          # case = current folder
    python -m cfdtools.pv cases/smoke              # ParaView window
    python -m cfdtools.pv cases/smoke --save       # image, no window
    python -m cfdtools.pv cases/smoke --field p --time 0.1 --range 0 0.1 -o fig.png

The view itself is defined in paraview_view.py, which runs inside ParaView.
(This module is not called paraview.py: ParaView would import it in place of
its own paraview package when running paraview_view.py.)

ParaView is looked for in this order:
  1. the folder in $PARAVIEW_BIN (containing paraview / pvbatch),
  2. paraview / pvbatch on the PATH,
  3. macOS: /Applications/ParaView-*.app (newest version first).
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).with_name("paraview_view.py")


def _version_key(path):
    return [int(n) for n in re.findall(r"\d+", Path(path).name)]


def find_executable(name):
    """Full path of a ParaView executable ('paraview' or 'pvbatch')."""
    env = os.environ.get("PARAVIEW_BIN")
    if env:
        candidate = Path(env) / name
        if candidate.exists():
            return str(candidate)

    found = shutil.which(name)
    if found:
        return found

    if sys.platform == "darwin":
        apps = sorted(glob.glob("/Applications/ParaView*.app"), key=_version_key, reverse=True)
        for app in apps:
            for sub in ("Contents/MacOS", "Contents/bin"):
                candidate = Path(app) / sub / name
                if candidate.exists():
                    return str(candidate)

    sys.exit(f"cannot find '{name}'. Install ParaView, put it on the PATH, "
             f"or set PARAVIEW_BIN to the folder that contains {name}.")


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m cfdtools.pv",
        description="Open an OpenFOAM case in ParaView with a publication view, "
                    "or save that view as an image.",
    )
    ap.add_argument("case", nargs="?", default=".", help="case folder (default: current folder)")
    ap.add_argument("--save", action="store_true",
                    help="write an image without opening a window "
                         "(default file: <case>/figures/paraview_<field>.png)")
    ap.add_argument("--field", default="U", help="cell field to colour by (default: U)")
    ap.add_argument("--time", default="latest", help="saved time (default: latest)")
    ap.add_argument("--range", nargs=2, type=float, metavar=("MIN", "MAX"),
                    help="fixed colour range (default: data range)")
    ap.add_argument("--size", nargs=2, type=int, default=[3000, 2000], metavar=("W", "H"),
                    help="image size in pixels (default: 3000 2000)")
    ap.add_argument("-o", "--output", help="image file; with a window, also saves it")
    args = ap.parse_args(argv)

    case = Path(args.case).resolve()
    if not case.is_dir():
        ap.error(f"no such case folder: {args.case}")

    options = ["--field", args.field, "--time", str(args.time),
               "--size", str(args.size[0]), str(args.size[1])]
    if args.range:
        options += ["--range", str(args.range[0]), str(args.range[1])]
    if args.output:
        options += ["-o", str(Path(args.output).resolve())]

    if args.save:
        cmd = [find_executable("pvbatch"), str(SCRIPT), str(case)] + options
        result = subprocess.run(cmd, capture_output=True, text=True)
        # Show pvbatch's own messages only on failure; it prints harmless
        # library warnings on every run.
        for line in result.stdout.splitlines():
            if line.startswith("wrote "):
                print(line)
        if result.returncode != 0:
            sys.stderr.write(result.stdout + result.stderr)
        return result.returncode

    # ParaView window: the script cannot take arguments here, so the case and
    # options go through environment variables.
    env = dict(os.environ, PV_CASE=str(case), PV_ARGS=" ".join(_quote(o) for o in options))
    subprocess.Popen([find_executable("paraview"), f"--script={SCRIPT}"], env=env,
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                     start_new_session=True)
    print(f"opening {case} in ParaView")
    return 0


def _quote(text):
    import shlex
    return shlex.quote(text)


if __name__ == "__main__":
    sys.exit(main())
