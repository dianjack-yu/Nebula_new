#!/usr/bin/env python3
"""
Inspect detected electron files (.det): print yields and counts per lineage,
plot histograms of every field and, optionally, the angular distribution.

Electrons moving up (dz > 0) when detected are counted as reflected, the others
as transmitted. SE and BSE are separated by energy (< 50 eV and >= 50 eV) and,
independently, by lineage tag.

Examples:
    python inspect_det.py point_10kV.det --pri point_10kV.pri
    python inspect_det.py *.det --pri point_10kV.pri --angular --energy_max 50
    python inspect_det.py scan.det --lineage 11 12 13 --px 0-49 --show
"""

import argparse
import os
import numpy as np
import matplotlib
import nebula_formats as nf


def parse_arguments():
    parser = argparse.ArgumentParser(description="Summarise and plot detected electron files.")
    parser.add_argument("det_files", nargs="+", help="One or more .det files")
    parser.add_argument("--pri", help="Primary file used for the simulation, to compute yields per incident electron")
    parser.add_argument("--n_incident", type=int, help="Number of incident electrons (alternative to --pri)")
    parser.add_argument("--energy_min", type=float, default=0.0, help="Minimum energy in eV")
    parser.add_argument("--energy_max", type=float, default=np.inf, help="Maximum energy in eV")
    parser.add_argument("--lineage", type=int, nargs="+", help="Lineage tags to include (default: all)")
    parser.add_argument("--px", nargs="+", help="px values/ranges to include, e.g. 5 20-40 (default: all)")
    parser.add_argument("--py", nargs="+", help="py values/ranges to include (default: all)")
    parser.add_argument("--bins", type=int, default=50, help="Histogram bins (default: 50)")
    parser.add_argument("--angular", action="store_true",
                        help="Also plot the polar angular distribution of reflected electrons "
                             "(per unit solid angle, with the Lambert cosine law)")
    parser.add_argument("--show", action="store_true", help="Show the plots instead of only saving PNGs")
    return parser.parse_args()


def filter_tag(args):
    """Suffix for output files describing the active filters, so runs don't overwrite each other."""
    parts = []
    if args.energy_min > 0 or np.isfinite(args.energy_max):
        hi = "inf" if np.isinf(args.energy_max) else f"{args.energy_max:g}"
        parts.append(f"E{args.energy_min:g}-{hi}")
    if args.lineage:
        parts.append("L" + "-".join(map(str, sorted(args.lineage))))
    if args.px:
        parts.append("px" + "_".join(args.px))
    if args.py:
        parts.append("py" + "_".join(args.py))
    return "".join("_" + p for p in parts)


def select(data, args):
    mask = (data['E'] >= args.energy_min) & (data['E'] <= args.energy_max)
    if args.lineage:
        mask &= np.isin(data['lineage'], args.lineage)
    px = nf.parse_index_list(args.px)
    py = nf.parse_index_list(args.py)
    if px is not None:
        mask &= np.isin(data['px'], px)
    if py is not None:
        mask &= np.isin(data['py'], py)
    return data[mask]


def print_summary(name, data, n_incident):
    def fmt(n):
        return f"{n:10d}" + (f"   yield {n / n_incident:.4f}" if n_incident else "")

    up = data['dz'] > 0
    se = data['E'] < nf.SE_ENERGY_LIMIT
    print(f"\n=== {name}: {len(data)} detected electrons"
          + (f", {n_incident} incident ===" if n_incident else " ==="))
    print("By energy:")
    print(f"  SE  (< 50 eV), reflected    {fmt(np.count_nonzero(up & se))}")
    print(f"  BSE (>= 50 eV), reflected   {fmt(np.count_nonzero(up & ~se))}")
    print(f"  < 50 eV, transmitted        {fmt(np.count_nonzero(~up & se))}")
    print(f"  >= 50 eV, transmitted       {fmt(np.count_nonzero(~up & ~se))}")
    print("By lineage (reflected / transmitted):")
    for tag in np.unique(data['lineage']):
        sel = data['lineage'] == tag
        print(f"  {nf.lineage_label(int(tag)):22s} {np.count_nonzero(sel & up):10d} {np.count_nonzero(sel & ~up):10d}")


def plot_fields(name, data, bins):
    import matplotlib.pyplot as plt
    fields = list(nf.DET_DTYPE.names) + ['polar angle (deg)']
    cols = 4
    rows = -(-len(fields) // cols)
    fig, axes = plt.subplots(rows, cols, figsize=(16, 3.2 * rows))
    angles = nf.polar_angle_deg(data)
    for ax, field in zip(axes.flat, fields):
        values = angles if field.startswith('polar') else data[field]
        if field == 'lineage':
            tags, counts = np.unique(values, return_counts=True)
            ax.bar([str(t) for t in tags], counts, color='steelblue')
        else:
            lo, hi = np.percentile(values, [0.5, 99.5])
            if lo == hi:
                lo, hi = lo - 0.5, hi + 0.5
            ax.hist(values, bins=bins, range=(lo, hi), color='steelblue')
        ax.set_title(field)
        ax.grid(axis='y', alpha=0.3)
    for ax in list(axes.flat)[len(fields):]:
        ax.axis('off')
    fig.suptitle(f"{name} (n = {len(data)})")
    fig.tight_layout()
    return fig


def plot_angular(all_data, bins):
    """Polar plot of reflected electrons per unit solid angle, normalised to the maximum."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': 'polar'})
    edges = np.linspace(0, 90, bins + 1)
    centres = np.radians(0.5 * (edges[1:] + edges[:-1]))
    # solid angle of each ring: 2 pi (cos(theta1) - cos(theta2))
    solid_angle = 2 * np.pi * (np.cos(np.radians(edges[:-1])) - np.cos(np.radians(edges[1:])))
    for name, data in all_data:
        angles = nf.polar_angle_deg(data[data['dz'] > 0])
        counts, _ = np.histogram(angles, bins=edges)
        density = counts / solid_angle
        if density.max() > 0:
            density /= density.max()
        theta = np.concatenate([-centres[::-1], centres])     # mirror to show both sides
        ax.plot(theta, np.concatenate([density[::-1], density]), label=name)
    theta = np.linspace(-np.pi / 2, np.pi / 2, 200)
    ax.plot(theta, np.cos(theta), 'k--', lw=1, label='Lambert (cos)')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    ax.set_title("Angular distribution of reflected electrons")
    ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1.1), fontsize='small')
    fig.tight_layout()
    return fig


def main():
    args = parse_arguments()
    if not args.show:
        matplotlib.use("agg")
    import matplotlib.pyplot as plt

    n_incident = args.n_incident or (nf.count_primaries(args.pri) if args.pri else None)
    all_data = []
    for path in args.det_files:
        name = os.path.splitext(os.path.basename(path))[0]
        data = select(nf.read_det(path), args)
        if len(data) == 0:
            print(f"\n{path}: no electrons left after filtering")
            continue
        print_summary(name, data, n_incident)
        fig = plot_fields(name, data, args.bins)
        out = os.path.splitext(path)[0] + filter_tag(args) + "_fields.png"
        fig.savefig(out, dpi=120)
        print(f"Saved {out}")
        all_data.append((name, data))

    if args.angular and all_data:
        fig = plot_angular(all_data, bins=18)
        out = os.path.splitext(args.det_files[0])[0] + filter_tag(args) + "_angular.png"
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print(f"Saved {out}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
