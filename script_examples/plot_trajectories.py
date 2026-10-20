#!/usr/bin/env python3
"""
Plot electron trajectories from an event file (.evnt) of nebula_cpu_all_events,
and print how often each event type occurred.

Each row of a .evnt file is one event: the position after the event, the energy
before it, the energy lost, origin_z, the pixel indices of the primary, the lineage
tag and the event ID. A path is the sequence of events with the same
(px, py, lineage, origin_z). The primary path (lineage 0) is joined to its BSE
continuation (lineage 1) at the turnaround point. Use make_point_probe_pri.py so
that every primary has a unique px.

Every particle starts with a zero-length row at its starting point. Secondaries
get two: one with their creation event ID (21, 22 or 41) when they are created,
and one with event ID 0 when their own simulation starts. In .evnt files event
ID 0 therefore marks the start of every particle, not only of primaries.

Examples:
    # first 20 primaries, coloured by lineage
    python plot_trajectories.py point_1kV.evnt

    # electrons 0-4, coloured by energy, with elastic (1) and detection events marked
    python plot_trajectories.py point_1kV.evnt --electrons 0-4 --color energy --events 1 991 992 993

    # only secondaries
    python plot_trajectories.py point_1kV.evnt --lineage 11 12 13 21 22 23
"""

import argparse
import os
import numpy as np
import matplotlib
import nebula_formats as nf

LINEAGE_COLORS = {
    0: 'tab:blue', 1: 'tab:orange',
    11: 'tab:green', 12: 'tab:olive', 13: 'tab:cyan',
    21: 'tab:red', 22: 'tab:purple', 23: 'tab:brown',
}
EVENT_MARKERS = ['o', 'v', '^', 's', 'P', '*', 'D', 'X', 'd', 'p', 'h', '8', '<', '>']


def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot trajectories from a .evnt file.")
    parser.add_argument("evnt_file", help="Event file written by nebula_cpu_all_events")
    parser.add_argument("--electrons", nargs="+", default=None,
                        help="px values/ranges of the primaries to plot, e.g. 0-9 15 "
                             "(default: the first 20 in the file; 'all' for all)")
    parser.add_argument("--lineage", type=int, nargs="+", help="Lineage tags to plot (default: all)")
    parser.add_argument("--color", choices=["lineage", "energy"], default="lineage",
                        help="Colour paths by lineage tag or by energy (default: lineage)")
    parser.add_argument("--events", type=int, nargs="+", default=[],
                        help="Event IDs to mark on the paths, e.g. 1 2 991")
    parser.add_argument("--xlim", type=float, nargs=2, help="Limits of the horizontal axes (nm)")
    parser.add_argument("--zlim", type=float, nargs=2, help="Limits of the z axis (nm)")
    parser.add_argument("--output", help="Output PNG (default: next to the input file)")
    parser.add_argument("--show", action="store_true", help="Show the plot")
    return parser.parse_args()


def order_rows(data):
    """Group rows by electron and path while keeping the time order within a path."""
    order = np.lexsort((np.arange(len(data)), data['origin_z'], data['lineage'], data['py'], data['px']))
    data = data[order]

    # Close each primary path at the point where it became a BSE: copy the first
    # lineage-1 row as the last lineage-0 row of the same electron.
    same = (data['px'][1:] == data['px'][:-1]) & (data['py'][1:] == data['py'][:-1])
    turn = np.flatnonzero(same & (data['lineage'][:-1] == 0) & (data['lineage'][1:] == 1)) + 1
    if len(turn):
        closing = data[turn].copy()
        closing['lineage'] = 0
        closing['origin_z'] = data['origin_z'][turn - 1]
        data = np.insert(data, turn, closing)
    return data


def split_paths(data):
    keys = ('px', 'py', 'lineage', 'origin_z')
    change = np.zeros(len(data) - 1, dtype=bool)
    for k in keys:
        change |= data[k][1:] != data[k][:-1]
    return np.split(data, np.flatnonzero(change) + 1)


def print_event_counts(data):
    ids, counts = np.unique(data['event_id'], return_counts=True)
    print("Event counts in the selection:")
    for i, n in zip(ids, counts):
        print(f"  {i:4d}  {nf.EVENT_NAMES.get(int(i), 'unknown'):38s} {n}")


def main():
    args = parse_arguments()
    if not args.show:
        matplotlib.use("agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    data = nf.read_evnt(args.evnt_file)
    if args.electrons is None:
        selected = np.unique(data['px'])[:20]
    else:
        selected = nf.parse_index_list(args.electrons)
    if selected is not None:
        data = data[np.isin(data['px'], selected)]
    if args.lineage:
        data = data[np.isin(data['lineage'], args.lineage)]
    if len(data) == 0:
        raise SystemExit("No events left after filtering.")
    print(f"{len(data)} events from {len(np.unique(data['px']))} primaries")
    print_event_counts(data)

    paths = [p for p in split_paths(order_rows(data)) if len(p) > 1]
    print(f"{len(paths)} paths")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    e_norm = plt.Normalize(data['E'].min(), data['E'].max())
    for ax, (h, title) in zip(axes, [('x', 'XZ'), ('y', 'YZ')]):
        segments, colors, energies = [], [], []
        for p in paths:
            pts = np.column_stack([p[h], p['z']])
            segments.extend(np.stack([pts[:-1], pts[1:]], axis=1))
            colors.extend([LINEAGE_COLORS.get(int(p['lineage'][0]), 'black')] * (len(p) - 1))
            energies.extend(0.5 * (p['E'][:-1] + p['E'][1:]))
        if args.color == 'energy':
            lc = LineCollection(segments, cmap='viridis', norm=e_norm, linewidths=0.8)
            lc.set_array(np.asarray(energies))
        else:
            lc = LineCollection(segments, colors=colors, linewidths=0.8)
        ax.add_collection(lc)

        for k, ev in enumerate(args.events):
            sel = data[data['event_id'] == ev]
            ax.plot(sel[h], sel['z'], linestyle='none', marker=EVENT_MARKERS[k % len(EVENT_MARKERS)],
                    ms=3, label=f"{ev}: {nf.EVENT_NAMES.get(ev, 'unknown')}")

        ax.autoscale()
        if args.xlim:
            ax.set_xlim(args.xlim)
        if args.zlim:
            ax.set_ylim(args.zlim)
        ax.set_aspect('equal', adjustable='datalim')
        ax.axhline(0, color='gray', lw=0.5)
        ax.set_xlabel(f"{h} (nm)")
        ax.set_ylabel("z (nm)")
        ax.set_title(f"{title} projection")

    if args.color == 'energy':
        fig.colorbar(lc, ax=axes, label='Energy (eV)', shrink=0.8)
    else:
        present = sorted(set(int(p['lineage'][0]) for p in paths))
        handles = [Line2D([], [], color=LINEAGE_COLORS.get(t, 'black'), label=nf.lineage_label(t))
                   for t in present]
        axes[1].legend(handles=handles, title="Lineage", loc='lower right', fontsize='small')
    if args.events:
        axes[0].legend(loc='lower left', fontsize='small', title="Events")

    out = args.output or os.path.splitext(args.evnt_file)[0] + f"_trajectories_{args.color}.png"
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print(f"Saved {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
