#!/usr/bin/env python3
"""
Generate a stationary point probe (.pri): N electrons entering at one position,
for interaction volumes, yields, escape depths and angular distributions.

Every electron gets a unique px index (py = 0). Secondaries inherit the pixel
indices of their primary, so px identifies the whole cascade of one primary,
which plot_trajectories.py uses to separate the paths.

Example:
    python make_point_probe_pri.py --output point_10kV.pri --energy 10000 --num_e 100000
"""

import argparse
import numpy as np
from nebula_formats import PRI_DTYPE


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate a point probe primary electron file.")
    parser.add_argument("--output", required=True, help="Output .pri file")
    parser.add_argument("--energy", type=float, default=10000, help="Beam energy in eV (default: 10000)")
    parser.add_argument("--num_e", type=int, default=100000, help="Number of electrons (default: 100000)")
    parser.add_argument("--x", type=float, default=0, help="Beam x position in nm (default: 0)")
    parser.add_argument("--y", type=float, default=0, help="Beam y position in nm (default: 0)")
    parser.add_argument("--z", type=float, default=10,
                        help="Starting z in nm; must be inside the geometry, between the sample "
                             "surface and the top detector (default: 10)")
    parser.add_argument("--sigma", type=float, default=0,
                        help="Standard deviation (nm) of a Gaussian beam spot (default: 0, ideal point)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for the beam spot")
    return parser.parse_args()


def main():
    args = parse_arguments()
    rng = np.random.default_rng(args.seed)
    N = args.num_e

    pri = np.zeros(N, dtype=PRI_DTYPE)
    pri['x'] = args.x + (rng.normal(0, args.sigma, N) if args.sigma > 0 else 0)
    pri['y'] = args.y + (rng.normal(0, args.sigma, N) if args.sigma > 0 else 0)
    pri['z'] = args.z
    pri['dz'] = -1               # straight down; dx = dy = 0
    pri['E'] = args.energy
    pri['origin_z'] = args.z     # origin of the primary = starting position
    pri['px'] = np.arange(N)     # unique electron ID
    pri['py'] = 0
    # lineage and event_id stay 0

    pri.tofile(args.output)
    print(f"Wrote {N} electrons of {args.energy:g} eV at ({args.x:g}, {args.y:g}, {args.z:g}) nm "
          f"to {args.output}")


if __name__ == "__main__":
    main()
