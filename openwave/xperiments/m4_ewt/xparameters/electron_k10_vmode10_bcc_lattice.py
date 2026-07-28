# electron_k10_vmode10_bcc_lattice.py
from openwave.xperiments.m4_ewt.xparameters.utils.geometry import bcc_lattice_positions
from openwave.common import constants

# ================================================================
# V_MODE=10: Gaussian profile + quintic saturation
# K=10: BCC lattice positions (nearest nodes)
# ================================================================
UNIVERSE_EDGE = 2e-15
TARGET_VOXELS = 55_000_000
K = 10
PERTURBATION = 0.02

# BCC lattice — generic placement, NOT 1-3-6
# Spacing chosen to give similar cluster size as golden angle
LOCK_SPACING = constants.EWAVE_LENGTH / UNIVERSE_EDGE
SPACING = 0.35 * LOCK_SPACING
POSITIONS = bcc_lattice_positions(K, center=(0.5, 0.5, 0.5), spacing=SPACING)
PHASES = [180] * K

XPARAMETERS = {
    "meta": {
        "X_NAME": "Electron K=10 V_MODE=10 BCC LATTICE (NOT 1-3-6)",
        "DESCRIPTION": "K=10 BCC lattice positions, control test",
    },
    "camera": {
        "INITIAL_POSITION": [0.94, 0.91, 0.69],
    },
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],
        "TARGET_VOXELS": TARGET_VOXELS,
    },
    "wave_centers": {
        "COUNT": K,
        "POSITION": POSITIONS,
        "PHASE_OFFSETS_DEG": PHASES,
        "APPLY_MOTION": True,
    },
    "engine": {
        "SEED_MODE": 2,
        "SEED_BOOST": 0.01,
        "V_MODE": 10,
        "V_C1": -0.3,
        "V_C2": 0.05,
        "WC_INTERACT_MODE": 3,
        "WC_BOOST": 1.0,
        "WC_RADIUS": 2,
        "WC_SIGMA": 1.5,
        "R_WALL": 100.0,
        "WALL_HEIGHT": 1.2,
        "DEFICIT_DEPTH": 0.9,
        "R_SOLITON": 35.0,
        "SIGMA": 3.0,
        "PRESSURE_STRENGTH": 0.001,
        "CFL_SAFETY": 0.1,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 1,
        "WARP_MESH": 100,
        "SHOW_GRANULES": False,
        "PARTICLE_SHELL": True,
        "SIM_SPEED": 1.0,
        "PAUSED": False,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 1,
    },
    "analytics": {
        "INSTRUMENTATION": True,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
}
