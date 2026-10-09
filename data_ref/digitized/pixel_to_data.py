#!/usr/bin/env python3
"""
pixel_to_data.py -- convert a pixel position on a digitized figure to data values.

Each axis is calibrated from its tick marks: between two neighbouring ticks the
value is interpolated linearly (in log10 of the value for a logarithmic axis).
Beyond the first or last tick the nearest segment is extended. The tick
positions are those in the header of each figure's CSV, and refer to the
images in data_ref/figures.

Figures:
    cho-kensey   Cho-kensey-fig1.jpg   x: shear rate (1/s)    y: viscosity (mPa s)
    chien-fig1   Chien1966-fig1.png    x: hematocrit (%)      y: viscosity (cP)
    chien-fig4   Chien-fig4.jpg        x: shear rate (1/s)    y: viscosity (cP)

chien-fig4 is tilted slightly: the viscosity ticks are 6-7 px higher at the
right axis than at the left, so their positions are interpolated across the
width before converting y.

Use:
    python pixel_to_data.py cho-kensey 1050 750
    python pixel_to_data.py chien-fig1 600 700 650 720     # several x y pairs

In Python:
    from pixel_to_data import pixel_to_data
    x, y = pixel_to_data("cho-kensey", 1050, 750)
"""

import argparse
import sys

import numpy as np


# Each axis: (tick pixel positions, tick values, "log" or "linear"). Log values are log10.
FIGURES = {
    "cho-kensey": dict(
        x=([524.5, 742.5, 961.5, 1173.5, 1386.5, 1603.5], [-1, 0, 1, 2, 3, 4], "log"),
        y=([117.0, 657.0, 1207.0], [2, 1, 0], "log"),
        units=("1/s", "mPa s"),
    ),
    "chien-fig1": dict(
        x=([367.0, 468.5, 570.0, 672.0, 774.0, 875.5], [0, 20, 40, 60, 80, 100], "linear"),
        y=([27.5, 286.0, 544.5, 803.5, 1062.5, 1324.0], [4, 3, 2, 1, 0, -1], "log"),
        units=("%", "cP"),
    ),
    "chien-fig4": dict(
        x=([233.5, 490.5, 746.0, 999.5, 1257.0, 1508.0], [-3, -2, -1, 0, 1, 2], "log"),
        # y ticks at the left axis (measured at x = 211) and at the right axis (measured at x = 1488)
        y=([69.0, 319.5, 570.5, 823.5], [3, 2, 1, 0], "log"),
        y_right=[62.0, 312.5, 563.5, 817.0],
        y_tick_x=(211.0, 1488.0),
        units=("1/s", "cP"),
    ),
}


def axis_value(p, ticks, values, kind):
    """Value at pixel p on an axis with ticks at pixel positions `ticks` having `values`."""
    ticks = np.asarray(ticks, dtype=float)
    values = np.asarray(values, dtype=float)
    if ticks[0] > ticks[-1]:                       # np.interp needs increasing positions
        ticks, values = ticks[::-1], values[::-1]
    if p < ticks[0]:
        v = values[0] + (p - ticks[0]) * (values[1] - values[0]) / (ticks[1] - ticks[0])
    elif p > ticks[-1]:
        v = values[-1] + (p - ticks[-1]) * (values[-1] - values[-2]) / (ticks[-1] - ticks[-2])
    else:
        v = np.interp(p, ticks, values)
    return float(10 ** v) if kind == "log" else float(v)


def pixel_to_data(figure, px, py):
    """(x value, y value) at pixel (px, py) on `figure` (see FIGURES)."""
    if figure not in FIGURES:
        raise ValueError(f"unknown figure {figure!r}; choose from {sorted(FIGURES)}")
    cal = FIGURES[figure]
    x = axis_value(px, *cal["x"])
    y_ticks, y_values, y_kind = cal["y"]
    if "y_right" in cal:                           # tilted scan: tick positions at this x
        x_left, x_right = cal["y_tick_x"]
        f = (px - x_left) / (x_right - x_left)
        y_ticks = [l + f * (r - l) for l, r in zip(y_ticks, cal["y_right"])]
    y = axis_value(py, y_ticks, y_values, y_kind)
    return x, y


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("figure", choices=sorted(FIGURES))
    ap.add_argument("pixels", type=float, nargs="+", help="x y [x y ...] pixel positions")
    args = ap.parse_args()
    if len(args.pixels) % 2:
        ap.error("give pixel positions in x y pairs")
    ux, uy = FIGURES[args.figure]["units"]
    for px, py in zip(args.pixels[0::2], args.pixels[1::2]):
        x, y = pixel_to_data(args.figure, px, py)
        print(f"({px:g}, {py:g}) px  ->  x = {x:.4g} {ux}, y = {y:.4g} {uy}")


if __name__ == "__main__":
    sys.exit(main())
