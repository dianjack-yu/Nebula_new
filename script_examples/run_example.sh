#!/bin/bash
# Example workflow: a thin membrane on water (liquid cell), simulated with Nebula 1.1.0.
#
#   1. geometry:    20 nm membrane (material 0) on water (material 1)
#   2. point probe: SE/BSE yields, escape depths and angular distributions at several energies (GPU)
#   3. raster scan: SE and BSE images (GPU)
#   4. all events:  trajectories of a few electrons (single-threaded CPU)
#
# Usage:
#   MEMBRANE_MAT=/path/to/siliconnitride.mat LIQUID_MAT=/path/to/water.mat ./run_example.sh
# Set NEBULA_BIN to the folder with the Nebula executables if they are not on your PATH,
# and PYTHON to a Python with numpy and matplotlib.
# Material files are made with cstool (https://github.com/nebula-simulator/cstool).

set -u
cd "$(dirname "$0")"

MEMBRANE_MAT=${MEMBRANE_MAT:-siliconnitride.mat}
LIQUID_MAT=${LIQUID_MAT:-water.mat}
NEBULA_BIN=${NEBULA_BIN:-}
PYTHON=${PYTHON:-python3}
OUT=${OUT:-example_output}

ENERGIES=(1000 3000 10000)     # point probe beam energies (eV)
N_POINT=100000                 # electrons per point probe
PHONON="200 200"               # phonon domain (eV) per material, same order as the .mat files
LOG="$OUT/error_log.txt"

nebula() { "${NEBULA_BIN:+$NEBULA_BIN/}$@"; }

mkdir -p "$OUT"
echo "Execution log - $(date)" > "$LOG"

# --- Check inputs before spending time on simulations ---
missing=0
for f in "$MEMBRANE_MAT" "$LIQUID_MAT"; do
    [[ -f "$f" ]] || { echo "Error: material file '$f' not found" | tee -a "$LOG"; missing=1; }
done
for exe in nebula_gpu nebula_cpu_all_events; do
    command -v "${NEBULA_BIN:+$NEBULA_BIN/}$exe" > /dev/null \
        || { echo "Error: $exe not found (set NEBULA_BIN)" | tee -a "$LOG"; missing=1; }
done
[[ $missing -eq 0 ]] || { echo "Aborted, see $LOG"; exit 1; }

# --- 1. Geometry ---
TRI="$OUT/film20_water.tri"
$PYTHON make_film_tri.py --layers 20 --output "$TRI"

# --- 2. Point probes ---
for E in "${ENERGIES[@]}"; do
    PRI="$OUT/point_${E}eV.pri"
    DET="$OUT/film20_water_point_${E}eV.det"
    $PYTHON make_point_probe_pri.py --output "$PRI" --energy "$E" --num_e $N_POINT

    # Low energies create fewer secondaries per primary, so a smaller batch factor is safer
    case $E in
        1000|3000) BATCH=0.7 ;;
        *)         BATCH=0.9 ;;
    esac

    echo "nebula_gpu: point probe, $E eV"
    # Positional arguments first: --phonon-domain takes all values up to the next -- flag
    if ! nebula nebula_gpu "$TRI" "$PRI" "$MEMBRANE_MAT" "$LIQUID_MAT" \
            --batch-factor $BATCH --phonon-domain $PHONON > "$DET"; then
        echo "nebula_gpu failed for $DET" | tee -a "$LOG"
    fi
done

# --- 3. Raster scan for images ---
SCAN_PRI="$OUT/scan_3000eV.pri"
SCAN_DET="$OUT/film20_water_scan_3000eV.det"
$PYTHON make_raster_scan_pri.py --output "$SCAN_PRI" --energy 3000 \
    --x_range -100 100 --y_range -100 100 --pixels 50 50 --epx 200
echo "nebula_gpu: raster scan"
nebula nebula_gpu "$TRI" "$SCAN_PRI" "$MEMBRANE_MAT" "$LIQUID_MAT" --phonon-domain $PHONON \
    > "$SCAN_DET" || echo "nebula_gpu failed for $SCAN_DET" | tee -a "$LOG"

# --- 4. All events for a few electrons (keep this small: every event is written) ---
EV_PRI="$OUT/point_1000eV_20e.pri"
$PYTHON make_point_probe_pri.py --output "$EV_PRI" --energy 1000 --num_e 20
echo "nebula_cpu_all_events: 20 electrons"
nebula nebula_cpu_all_events "$TRI" "$EV_PRI" "$MEMBRANE_MAT" "$LIQUID_MAT" \
    --detect-filename "$OUT/point_1000eV_20e.det" \
    --events-filename "$OUT/point_1000eV_20e.evnt" \
    --phonon-domain $PHONON || echo "nebula_cpu_all_events failed" | tee -a "$LOG"

# --- 5. Analysis ---
for E in "${ENERGIES[@]}"; do
    $PYTHON inspect_det.py "$OUT/film20_water_point_${E}eV.det" --pri "$OUT/point_${E}eV.pri"
done
$PYTHON inspect_det.py "$OUT"/film20_water_point_*eV.det --angular --energy_max 50
$PYTHON origin_z_distribution.py "$OUT/film20_water_point_3000eV.det" --pri "$OUT/point_3000eV.pri" \
    --energy_max 50 --reflected --split SE1,SE2 --normalize --z_range -15 0 --bins 60 \
    --output "$OUT/origin_z_3000eV.png"
$PYTHON sem_image.py "$SCAN_DET" --energy_max 50 --epx 200
$PYTHON sem_image.py "$SCAN_DET" --energy_min 50 --epx 200
$PYTHON plot_trajectories.py "$OUT/point_1000eV_20e.evnt"

echo "Done. Results are in $OUT/, errors (if any) in $LOG"
