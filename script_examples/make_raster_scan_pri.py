#!/usr/bin/env python3
"""
Generate a raster scan (.pri) over a rectangular area, for simulated SEM images
(see sem_image.py). Pixel (px, py) is scanned with a fixed number of electrons,
optionally with Poisson shot noise and a Gaussian beam spot.

Example: 100 x 100 pixels over 1 x 1 um, 1000 electrons per pixel at 5 keV
    python make_raster_scan_pri.py --output scan_5kV.pri --energy 5000 \
        --x_range -500 500 --y_range -500 500 --pixels 100 100 --epx 1000
"""

import argparse
import numpy as np
from nebula_formats import PRI_DTYPE


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate a raster scan primary electron file.")
    parser.add_argument("--output", required=True, help="Output .pri file")
    parser.add_argument("--energy", type=float, default=5000, help="Beam energy in eV (default: 5000)")
    parser.add_argument("--x_range", type=float, nargs=2, default=[-500, 500], metavar=("XMIN", "XMAX"),
                        help="Scan range in x (nm), pixel centres (default: -500 500)")
    parser.add_argument("--y_range", type=float, nargs=2, default=[-500, 500], metavar=("YMIN", "YMAX"),
                        help="Scan range in y (nm), pixel centres (default: -500 500)")
    parser.add_argument("--pixels", type=int, nargs=2, default=[100, 100], metavar=("NX", "NY"),
                        help="Number of pixels in x and y (default: 100 100)")
    parser.add_argument("--epx", type=int, default=1000, help="Electrons per pixel (default: 1000)")
    parser.add_argument("--z", type=float, default=10,
                        help="Starting z in nm, between the sample surface and the top detector (default: 10)")
    parser.add_argument("--sigma", type=float, default=0,
                        help="Standard deviation (nm) of a Gaussian beam spot (default: 0)")
    parser.add_argument("--poisson", action="store_true",
                        help="Draw the number of electrons per pixel from a Poisson distribution")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    return parser.parse_args()


def main():
    args = parse_arguments()
    rng = np.random.default_rng(args.seed)
    nx, ny = args.pixels
    xs = np.linspace(*args.x_range, nx)
    ys = np.linspace(*args.y_range, ny)

    total = 0
    with open(args.output, "wb") as f:
        for iy, y in enumerate(ys):           # one scan line at a time keeps memory low
            counts = rng.poisson(args.epx, nx) if args.poisson else np.full(nx, args.epx)
            n = int(counts.sum())
            line = np.zeros(n, dtype=PRI_DTYPE)
            line['x'] = np.repeat(xs, counts)
            line['y'] = y
            if args.sigma > 0:
                line['x'] += rng.normal(0, args.sigma, n)
                line['y'] += rng.normal(0, args.sigma, n)
            line['z'] = args.z
            line['dz'] = -1
            line['E'] = args.energy
            line['origin_z'] = args.z
            line['px'] = np.repeat(np.arange(nx), counts)
            line['py'] = iy
            line.tofile(f)
            total += n

    step_x = (xs[1] - xs[0]) if nx > 1 else 0
    step_y = (ys[1] - ys[0]) if ny > 1 else 0
    print(f"Wrote {total} electrons ({nx} x {ny} pixels, pixel size {step_x:g} x {step_y:g} nm, "
          f"{args.energy:g} eV) to {args.output}")


if __name__ == "__main__":
    main()
