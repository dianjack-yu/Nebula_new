#!/usr/bin/env python3
"""
Plot the distribution of origin_z (escape depth statistics) of detected electrons.

origin_z is the z position where the current state of an electron originated:
  - SE (lineage 11, 12, 13, 21, 22, 23): where the secondary was generated
  - BSE (lineage 1): where the primary turned around (direction changed to dz > 0)
  - primaries (lineage 0): the origin_z given in the .pri file
With the surface at z = 0, -origin_z is the generation depth below the surface.

Examples:
    # SE generation depth of reflected SE, one curve per file
    python origin_z_distribution.py film_*.det --energy_max 50 --reflected --z_range -100 0

    # SE1 vs SE2 in one file, as electron yield per nm per incident electron, with CSV output
    python origin_z_distribution.py film20_10kV.det --pri point_10kV.pri --energy_max 50 \
        --reflected --split SE1,SE2 --normalize --save_csv
"""

import argparse
import os
import numpy as np
import matplotlib
import nebula_formats as nf

GROUPS = {
    "SE1": nf.SE1_LINEAGES,
    "SE2": nf.SE2_LINEAGES,
    "SE": nf.SE1_LINEAGES + nf.SE2_LINEAGES,
    "BSE": (1,),
    "primary": (0,),
}


def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot origin_z distributions of detected electrons.")
    parser.add_argument("det_files", nargs="+", help="One or more .det files")
    parser.add_argument("--energy_min", type=float, default=0.0, help="Minimum energy in eV")
    parser.add_argument("--energy_max", type=float, default=np.inf, help="Maximum energy in eV")
    parser.add_argument("--reflected", action="store_true",
                        help="Only electrons moving up (dz > 0) when detected")
    parser.add_argument("--lineage", type=int, nargs="+", help="Lineage tags to include (default: all)")
    parser.add_argument("--split", default=None,
                        help="Plot separate curves for lineage groups, comma separated, from: "
                             + ", ".join(GROUPS))
    parser.add_argument("--px", nargs="+", help="px values/ranges to include (default: all)")
    parser.add_argument("--py", nargs="+", help="py values/ranges to include (default: all)")
    parser.add_argument("--z_range", type=float, nargs=2, metavar=("ZMIN", "ZMAX"),
                        help="Histogram range in nm (default: data range)")
    parser.add_argument("--bins", type=int, default=100, help="Number of bins (default: 100)")
    parser.add_argument("--normalize", action="store_true",
                        help="Yield per nm per incident electron: counts / (N_incident x bin width). "
                             "Needs --pri or --n_incident")
    parser.add_argument("--pri", help="Primary file of the simulation (gives N_incident)")
    parser.add_argument("--n_incident", type=int, help="Number of incident electrons")
    parser.add_argument("--cumulative", action="store_true", help="Plot the cumulative percentage")
    parser.add_argument("--save_csv", action="store_true", help="Save each curve as CSV")
    parser.add_argument("--output", default="origin_z_distribution.png", help="Output figure")
    parser.add_argument("--show", action="store_true", help="Show the plot")
    return parser.parse_args()


def select(data, args):
    mask = (data['E'] >= args.energy_min) & (data['E'] <= args.energy_max)
    if args.reflected:
        mask &= data['dz'] > 0
    if args.lineage:
        mask &= np.isin(data['lineage'], args.lineage)
    px = nf.parse_index_list(args.px)
    py = nf.parse_index_list(args.py)
    if px is not None:
        mask &= np.isin(data['px'], px)
    if py is not None:
        mask &= np.isin(data['py'], py)
    return data[mask]


def main():
    args = parse_arguments()
    if not args.show:
        matplotlib.use("agg")
    import matplotlib.pyplot as plt

    n_incident = args.n_incident or (nf.count_primaries(args.pri) if args.pri else None)
    if args.normalize and not n_incident:
        raise SystemExit("--normalize needs --pri or --n_incident")
    groups = [g.strip() for g in args.split.split(",")] if args.split else [None]
    for g in groups:
        if g is not None and g not in GROUPS:
            raise SystemExit(f"Unknown group '{g}', choose from {', '.join(GROUPS)}")

    # Common bin edges for all curves
    selected = [(p, select(nf.read_det(p), args)) for p in args.det_files]
    if args.z_range:
        z_lo, z_hi = args.z_range
    else:
        all_z = np.concatenate([d['origin_z'] for _, d in selected])
        if len(all_z) == 0:
            raise SystemExit("No electrons left after filtering.")
        z_lo, z_hi = float(all_z.min()), float(all_z.max())
    edges = np.linspace(z_lo, z_hi, args.bins + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for path, data in selected:
        name = os.path.splitext(os.path.basename(path))[0]
        for g in groups:
            d = data if g is None else data[np.isin(data['lineage'], GROUPS[g])]
            counts, _ = np.histogram(d['origin_z'], bins=edges)
            if args.cumulative:
                y = 100 * np.cumsum(counts) / max(counts.sum(), 1)
                column = "cumulative_percent"
            elif args.normalize:
                y = counts / (n_incident * (edges[1] - edges[0]))
                column = "yield_per_nm"
            else:
                y = counts
                column = "counts"
            label = name if g is None else f"{name} {g}"
            ax.plot(centres, y, lw=1.5, label=f"{label} (n={len(d)})")
            if args.save_csv:
                csv = os.path.splitext(path)[0] + (f"_{g}" if g else "") + f"_origin_z_{column}.csv"
                np.savetxt(csv, np.column_stack([centres, y]), delimiter=",",
                           header=f"origin_z_nm,{column}", comments="")
                print(f"Saved {csv}")

    ax.set_xlabel("origin_z (nm)")
    ax.set_ylabel({"cumulative_percent": "Cumulative (%)",
                   "yield_per_nm": "Yield per nm per incident electron",
                   "counts": "Counts"}[column])
    ax.set_title(f"origin_z, E = {args.energy_min:g}-{args.energy_max:g} eV"
                 + (", reflected only" if args.reflected else ""))
    ax.grid(alpha=0.3)
    ax.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
