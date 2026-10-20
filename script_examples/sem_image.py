#!/usr/bin/env python3
"""
Form a simulated SEM image from a raster scan (.det from make_raster_scan_pri.py):
count the detected electrons in each (px, py) pixel after filtering by energy,
emission angle, origin z and lineage, like an energy/angle-selective detector.

Examples:
    # SE image (< 50 eV) and BSE image (>= 50 eV)
    python sem_image.py scan_5kV.det --energy_max 50
    python sem_image.py scan_5kV.det --energy_min 50

    # only SE1 electrons, saved as PNG and 32-bit TIFF
    python sem_image.py scan_5kV.det --energy_max 50 --lineage 11 12 13 --tiff

    # all .det files in a folder, divided by the electrons per pixel
    python sem_image.py ./results --energy_max 50 --epx 1000
"""

import argparse
import os
import numpy as np
import matplotlib
import nebula_formats as nf


def parse_arguments():
    parser = argparse.ArgumentParser(description="Form SEM images from raster scan .det files.")
    parser.add_argument("input", help="A .det file or a folder containing .det files")
    parser.add_argument("--output_folder", default=None, help="Output folder (default: next to the input)")
    parser.add_argument("--energy_min", type=float, default=0.0, help="Minimum energy in eV (default: 0)")
    parser.add_argument("--energy_max", type=float, default=np.inf, help="Maximum energy in eV (default: none)")
    parser.add_argument("--angle_min", type=float, default=0.0,
                        help="Minimum angle from the +z axis in degrees (default: 0)")
    parser.add_argument("--angle_max", type=float, default=90.0,
                        help="Maximum angle from the +z axis in degrees; 90 keeps reflected "
                             "electrons only, 180 also transmitted ones (default: 90)")
    parser.add_argument("--origin_z_min", type=float, default=-np.inf, help="Minimum origin_z in nm")
    parser.add_argument("--origin_z_max", type=float, default=np.inf, help="Maximum origin_z in nm")
    parser.add_argument("--lineage", type=int, nargs="+", help="Lineage tags to include (default: all)")
    parser.add_argument("--epx", type=float, default=None,
                        help="Electrons per pixel of the scan; if given, the image shows the yield "
                             "per incident electron instead of counts")
    parser.add_argument("--tiff", action="store_true", help="Also save a 32-bit float TIFF (needs tifffile)")
    parser.add_argument("--show", action="store_true", help="Show the image")
    return parser.parse_args()


def settings_tag(args):
    def num(v):
        return "inf" if np.isinf(v) else f"{v:g}"
    tag = f"E{num(args.energy_min)}-{num(args.energy_max)}_A{num(args.angle_min)}-{num(args.angle_max)}"
    if np.isfinite(args.origin_z_min) or np.isfinite(args.origin_z_max):
        tag += f"_OZ{num(args.origin_z_min)}-{num(args.origin_z_max)}"
    tag += "_Lall" if args.lineage is None else "_L" + "-".join(map(str, sorted(args.lineage)))
    return tag


def make_image(data, args):
    angles = nf.polar_angle_deg(data)
    mask = ((data['E'] >= args.energy_min) & (data['E'] <= args.energy_max)
            & (angles >= args.angle_min) & (angles <= args.angle_max)
            & (data['origin_z'] >= args.origin_z_min) & (data['origin_z'] <= args.origin_z_max))
    if args.lineage:
        mask &= np.isin(data['lineage'], args.lineage)

    # Image size from all pixels in the file, so filtered images keep the same shape
    nx = data['px'].max() + 1
    ny = data['py'].max() + 1
    image = np.zeros((ny, nx), dtype=np.float32)
    np.add.at(image, (data['py'][mask], data['px'][mask]), 1)
    if args.epx:
        image /= args.epx
    return image, np.count_nonzero(mask)


def process(path, args):
    import matplotlib.pyplot as plt
    data = nf.read_det(path)
    if len(data) == 0:
        print(f"{path}: empty file")
        return
    image, n = make_image(data, args)
    print(f"{path}: {n} of {len(data)} detected electrons pass the filters, "
          f"image {image.shape[1]} x {image.shape[0]} pixels")

    folder = args.output_folder or os.path.dirname(os.path.abspath(path))
    os.makedirs(folder, exist_ok=True)
    stem = os.path.join(folder, f"{os.path.splitext(os.path.basename(path))[0]}_{settings_tag(args)}")

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(image, cmap='gray', origin='lower')
    fig.colorbar(im, ax=ax, label='Yield per incident electron' if args.epx else 'Counts')
    ax.set_xlabel('px')
    ax.set_ylabel('py')
    ax.set_title(os.path.basename(stem), fontsize=9)
    fig.tight_layout()
    fig.savefig(stem + ".png", dpi=150)
    print(f"  Saved {stem}.png")

    if args.tiff:
        try:
            import tifffile
            tifffile.imwrite(stem + ".tif", image)
            print(f"  Saved {stem}.tif")
        except ImportError:
            print("  tifffile is not installed; skipped TIFF (pip install tifffile)")

    if not args.show:
        plt.close(fig)


def main():
    args = parse_arguments()
    if not args.show:
        matplotlib.use("agg")
    import matplotlib.pyplot as plt

    if os.path.isdir(args.input):
        paths = sorted(os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith(".det"))
        if not paths:
            print(f"No .det files in {args.input}")
    else:
        paths = [args.input]
    for path in paths:
        process(path, args)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
