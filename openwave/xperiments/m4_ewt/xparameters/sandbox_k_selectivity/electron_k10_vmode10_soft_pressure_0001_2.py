# electron_k10_vmode10_soft_pressure_0001.py
"""
K=10, 1-3-6 geometry, soft interaction, V_MODE=10, ultra-low boost.
Increased vacuum pressure (PRESSURE_STRENGTH=0.001) to test enhanced stability.
"""

from openwave.xperiments.m4_ewt.xparameters.utils.geometry import (
    generate_positions_by_EWT_geometry_locked,
)
from openwave.common import constants

UNIVERSE_EDGE = 4e-15
TARGET_VOXELS = 55_000_000
K = 10

PERTURBATION = 0.1
POSITIONS = generate_positions_by_EWT_geometry_locked(
    UNIVERSE_EDGE, K, center=(0.5, 0.5, 0.5), perturbation=PERTURBATION
)
PHASES = [180] * K

XPARAMETERS = {
    "meta": {
        "X_NAME": "K=10 1-3-6 soft boost=0.00001 V_MODE=10 pressure=0.001",
        "DESCRIPTION": "Increased pressure test: does stronger vacuum pressure improve soliton stability?",
    },
    "camera": {"INITIAL_POSITION": [1.20, 1.80, 1.20]},
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
        "SEED_BOOST": 0.001,
        "V_MODE": 10,
        "V_C1": -0.3,
        "V_C2": 0.5,
        "WC_INTERACT_MODE": 1,
        "WC_BOOST": 0.002,
        "WC_RADIUS": 30,
        "WC_SIGMA": 1.5,
        "R_WALL": 100.0,
        "WALL_HEIGHT": 1.2,
        "DEFICIT_DEPTH": 0.9,
        "R_SOLITON": 35.0,
        "SIGMA": 3.0,
        "PRESSURE_STRENGTH": 0.001,
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
        "PARTICLE_SHELL": False,
        "SIM_SPEED": 1.0,
        "PAUSED": False,
    },
    "color_defaults": {"COLOR_THEME": "OCEAN", "WAVE_MENU": 1},
    "analytics": {"INSTRUMENTATION": True, "EXPORT_VIDEO": False, "VIDEO_FRAMES": 24},
}
