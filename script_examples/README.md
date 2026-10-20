# Nebula 1.1.0 example scripts

Python and shell scripts to set up, run and analyse Nebula simulations with the
file formats of version 1.1.0 (origin_z, lineage tag and event ID).
They need Python 3 with numpy and matplotlib; `tifffile` is optional (TIFF output of `sem_image.py`).
Run every script with `--help` to see all options.

| Script | Purpose |
|---|---|
| `nebula_formats.py` | File formats (dtypes), event ID and lineage definitions, helpers. Imported by the other scripts; keep it in the same folder. |
| `make_film_tri.py` | Geometry: stack of flat films on a substrate (e.g. a membrane on a liquid) inside a detector box. |
| `make_point_probe_pri.py` | Primaries: stationary point probe, one unique `px` per electron. |
| `make_raster_scan_pri.py` | Primaries: raster scan over a pixel grid, for images. |
| `run_example.sh` | Complete workflow: geometry, point probes and a raster scan on the GPU, an all-events run on the CPU, and the analysis below. |
| `inspect_det.py` | Detected electrons: SE/BSE yields, counts per lineage, histograms of every field, angular distribution. |
| `origin_z_distribution.py` | Escape depth statistics: distribution of `origin_z`, per file or per lineage group (SE1, SE2, BSE, ...). |
| `sem_image.py` | SEM image from a raster scan, filtered by energy, emission angle, `origin_z` and lineage. |
| `plot_trajectories.py` | Electron trajectories and event statistics from `nebula_cpu_all_events` output (`.evnt`). |

## Quick start

```bash
# 1. geometry: 20 nm membrane (material 0) on water (material 1)
python make_film_tri.py --layers 20 --output film20_water.tri

# 2. primaries: 100k electrons of 3 keV entering at (0, 0, 10) nm
python make_point_probe_pri.py --output point_3kV.pri --energy 3000 --num_e 100000

# 3. simulate; materials in the order of the geometry, positional arguments first
nebula_gpu film20_water.tri point_3kV.pri siliconnitride.mat water.mat \
    --phonon-domain 200 200 > film20_water_3kV.det

# 4. analyse
python inspect_det.py film20_water_3kV.det --pri point_3kV.pri --angular
python origin_z_distribution.py film20_water_3kV.det --pri point_3kV.pri \
    --energy_max 50 --reflected --split SE1,SE2 --normalize --z_range -15 0
```

Or run the whole workflow (about two minutes on an RTX 4080):

```bash
MEMBRANE_MAT=/path/to/siliconnitride.mat LIQUID_MAT=/path/to/water.mat \
NEBULA_BIN=/path/to/nebula/build/bin ./run_example.sh
```

Material files (`.mat`) are generated with [cstool](https://github.com/nebula-simulator/cstool).
The flat film of the example gives a uniform SEM image; use a structured geometry
for meaningful image contrast.

## Things to know

- **Material order.** Material index *i* in the `.tri` file is the *i*-th `.mat` file on the
  command line. `--phonon-domain` takes one value (eV) per material in the same order and collects
  all values up to the next `--` flag, so put it after the positional arguments.
- **Starting position.** Primaries must start inside the detector box and, by default, in vacuum
  (`--initial-material -123`). With `make_film_tri.py` this means 0 < z < `--det_top`.
- **Pixel indices.** Secondaries inherit `px` and `py` from their primary. In a point probe a
  unique `px` per primary lets you follow each cascade; in a raster scan `(px, py)` is the image pixel.
- **Reflected and transmitted.** A detected electron with `dz > 0` was moving up (reflected), otherwise
  down (transmitted). SE and BSE are separated either by energy (< 50 eV / >= 50 eV) or by lineage tag.
- **origin_z.** Where the current state of the electron originated: the generation depth of an SE,
  the turnaround depth of a BSE, or the value in the `.pri` file for a primary.
- **Event files** of `nebula_cpu_all_events` contain one row per event and grow quickly; simulate
  only a few low-energy primaries with it. Event ID 0 marks the start of every particle,
  secondaries included.

See the README of the Nebula repository for the complete event ID and lineage tag tables.
