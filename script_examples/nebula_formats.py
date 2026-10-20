"""
Shared file formats and tag definitions for Nebula 1.1.0.

All example scripts import this module, so keep it in the same folder.
The definitions match the README of the Nebula repository.
"""

import os
import numpy as np

# ---------------------------------------------------------------------------
# Binary file formats (native byte order, 32-bit floats and ints)
# ---------------------------------------------------------------------------

# Primary electrons (.pri), read by all Nebula programs.
# 'lineage' and 'event_id' are read but ignored: primaries always start at 0.
PRI_DTYPE = np.dtype([
    ('x',  '=f'), ('y',  '=f'), ('z',  '=f'),  # Position (nm)
    ('dx', '=f'), ('dy', '=f'), ('dz', '=f'),  # Direction (unnormalised)
    ('E',  '=f'), ('origin_z', '=f'),          # Energy (eV) and origin z position (nm)
    ('px', '=i'), ('py', '=i'),                # Pixel indices
    ('lineage', '=i'), ('event_id', '=i'),     # Lineage tag and event ID
])

# Detected electrons (.det): stdout of nebula_gpu / nebula_cpu_mt,
# --detect-filename of nebula_cpu_edep / nebula_cpu_all_events.
DET_DTYPE = np.dtype([
    ('x', '=f'), ('y', '=f'), ('z', '=f'),
    ('dx', '=f'), ('dy', '=f'), ('dz', '=f'),
    ('E', '=f'),
    ('origin_z', '=f'),
    ('px', '=i'), ('py', '=i'),
    ('lineage', '=i'),
])

# Event records (.evnt): --events-filename of nebula_cpu_all_events.
# One row per scattering or boundary event; E is the energy before the event,
# delE the energy lost in it, and x, y, z the position after it.
EVNT_DTYPE = np.dtype([
    ('x', '=f'), ('y', '=f'), ('z', '=f'),
    ('E', '=f'), ('delE', '=f'), ('origin_z', '=f'),
    ('px', '=i'), ('py', '=i'), ('lineage', '=i'), ('event_id', '=i'),
])

# Energy deposits: --deposit-filename of nebula_cpu_edep (EVNT_DTYPE without origin_z).
EDEP_DTYPE = np.dtype([
    ('x', '=f'), ('y', '=f'), ('z', '=f'),
    ('E', '=f'), ('delE', '=f'),
    ('px', '=i'), ('py', '=i'), ('lineage', '=i'), ('event_id', '=i'),
])

# ---------------------------------------------------------------------------
# Tag definitions
# ---------------------------------------------------------------------------

EVENT_NAMES = {
    0:   "Primary creation",
    1:   "Elastic scattering",
    11:  "Elastic, acoustic phonon loss",
    12:  "Elastic, atomic recoil loss",
    2:   "Fermi sea excitation",
    21:  "SE by plasmon decay",
    22:  "SE by single excitation",
    3:   "Sub-bandgap energy loss",
    31:  "Optical phonon loss (Kieft)",
    4:   "Inner-shell excitation",
    41:  "SE from inner-shell excitation",
    42:  "Excitation w/o momentum cons. (Kieft)",
    991: "Detection",
    992: "SE detection",
    993: "BSE/TE detection",
    994: "Termination",
    995: "Mirror",
    996: "Transmission without refraction",
    997: "Transmission with refraction",
    998: "Absorption at interface",
    999: "Total internal reflection",
}

DETECTION_EVENTS = (991, 992, 993)

LINEAGE_NAMES = {
    0:  "Primary",
    1:  "BSE",
    11: "SE1",
    21: "SE2",
    12: "SE1 2nd gen",
    22: "SE2 2nd gen",
    13: "SE1 3rd+ gen",
    23: "SE2 3rd+ gen",
}

# Lineage groups that are often analysed together
SE1_LINEAGES = (11, 12, 13)
SE2_LINEAGES = (21, 22, 23)

# Conventional energy limit between secondary and backscattered electrons (eV)
SE_ENERGY_LIMIT = 50.0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def read_binary(path, dtype):
    """Read a Nebula binary file, warning if its size does not fit the format."""
    size = os.path.getsize(path)
    if size % dtype.itemsize != 0:
        print(f"WARNING: size of '{path}' ({size} bytes) is not a multiple of "
              f"{dtype.itemsize} bytes. Was it written by an older Nebula version?")
    return np.fromfile(path, dtype=dtype, count=size // dtype.itemsize)


def read_det(path):
    return read_binary(path, DET_DTYPE)


def read_evnt(path):
    return read_binary(path, EVNT_DTYPE)


def count_primaries(pri_path):
    """Number of primary electrons in a .pri file."""
    return os.path.getsize(pri_path) // PRI_DTYPE.itemsize


def polar_angle_deg(data):
    """Angle (degrees) between the electron direction and the +z axis."""
    norm = np.sqrt(data['dx']**2 + data['dy']**2 + data['dz']**2)
    norm[norm == 0] = np.inf
    return np.degrees(np.arccos(np.clip(data['dz'] / norm, -1, 1)))


def parse_index_list(tokens):
    """
    Parse pixel selections such as ['5'], ['20-40'], ['0-9', '50'] or ['all'].
    Returns a sorted list of ints, or None meaning 'all'.
    """
    if not tokens or any(str(t).lower() == 'all' for t in tokens):
        return None
    values = set()
    for token in map(str, tokens):
        lo, sep, hi = token.partition('-')
        if sep and lo and hi:
            values.update(range(int(lo), int(hi) + 1))
        else:
            values.add(int(token))
    return sorted(values)


def lineage_label(tag):
    return f"{tag} ({LINEAGE_NAMES.get(tag, 'unexpected')})"
