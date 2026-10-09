#!/usr/bin/env python3
"""
make_check.py -- draw digitized points over the figure they were read from.

Reads a digitized CSV (with pixel_x and pixel_y columns, and optionally a
confidence column) and the source image named in its header, and writes a copy
of the image with every point marked by a cross and numbered by its row in the
CSV (0 = first data row). Cross colour shows the confidence: green high, orange
medium, red low. Points without a confidence column are drawn green.

The output has the same size as the source image and is not cropped or scaled,
so a pixel position read in an image editor (e.g. GIMP) is the same pixel
position as in the CSV and can be converted with pixel_to_data.py.

The image is found from the header line "# Source image: <file>", relative to
the CSV's folder. Pass --image to use another file.

Use (from the data_ref/digitized folder):
    python make_check.py chien-1966-fig1-digitized.csv
    python make_check.py my-points.csv --image ../figures/Some-figure.png -o out.png

By default the output goes to the digitization-checks folder next to this
script, as <csv name>-check.png.
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont


COLOURS = {"high": (0, 200, 0), "medium": (255, 140, 0), "low": (230, 0, 0)}
OTHER = (0, 160, 255)           # a confidence value not in COLOURS
LABEL = (220, 0, 220)
ARM = 8                         # half-length of the cross, pixels
WIDTH = 2                       # line width of the cross, pixels
FONT_SIZE = 14


def source_image(csv_path):
    """The file named on the '# Source image:' header line, relative to the CSV."""
    for line in Path(csv_path).read_text().splitlines():
        if not line.startswith("#"):
            break
        match = re.search(r"Source image:\s*([^\s(]+)", line)
        if match:
            return Path(csv_path).parent / match.group(1).rstrip(".,;")
    return None


def make_check(csv_path, image_path=None, out_path=None):
    csv_path = Path(csv_path)
    image_path = Path(image_path) if image_path else source_image(csv_path)
    if image_path is None:
        raise SystemExit(f"no '# Source image:' line in {csv_path}; pass --image")
    if not image_path.is_file():
        raise SystemExit(f"image not found: {image_path}")
    out_path = Path(out_path) if out_path else Path(__file__).parent / "digitization-checks" / f"{csv_path.stem}-check.png"

    data = pd.read_csv(csv_path, comment="#")
    for column in ("pixel_x", "pixel_y"):
        if column not in data.columns:
            raise SystemExit(f"{csv_path} has no {column} column")
    if "confidence" not in data.columns:
        data["confidence"] = "high"

    image = Image.open(image_path).convert("RGB")      # same size as the source: pixels match the CSV
    draw = ImageDraw.Draw(image)
    font = ImageFont.load_default(size=FONT_SIZE)
    for i, row in data.iterrows():
        x, y = row["pixel_x"], row["pixel_y"]
        if pd.isna(x) or pd.isna(y):
            continue
        colour = COLOURS.get(str(row["confidence"]), OTHER)
        draw.line([(x - ARM, y), (x + ARM, y)], fill=colour, width=WIDTH)
        draw.line([(x, y - ARM), (x, y + ARM)], fill=colour, width=WIDTH)
        draw.text((x + ARM + 2, y - ARM - FONT_SIZE), str(i), fill=LABEL, font=font)

    image.save(out_path)
    return out_path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="digitized CSV with pixel_x and pixel_y columns")
    ap.add_argument("--image", help="source image (default: the '# Source image:' header line)")
    ap.add_argument("-o", "--output", help="output PNG (default: digitization-checks/<csv name>-check.png next to this script)")
    args = ap.parse_args()
    print(f"wrote {make_check(args.csv, args.image, args.output)}")


if __name__ == "__main__":
    sys.exit(main())
