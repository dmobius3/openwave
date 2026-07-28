# electron_k09_ewt_vmode10_soft_pressure_00001.py
"""
K=9, EWT geometry (tricapped trigonal prism), soft interaction, V_MODE=10, ultra-low boost.
Vacuum pressure force enabled (PRESSURE_STRENGTH=0.0001).
Tests whether the EWT-native geometry with K=9 behaves like 1-3-6 (stable) or golden-angle (metastable).
"""

from openwave.xperiments.m4_ewt.xparameters.geometry import generate_positions_by_EWT_geometry
from openwave.common import constants

UNIVERSE_EDGE = 2e-15
TARGET_VOXELS = 55_000_000
K = 9
PERTURBATION = 0.02

POSITIONS = generate_positions_by_EWT_geometry(
    UNIVERSE_EDGE, K, center=(0.5, 0.5, 0.5), perturbation=PERTURBATION
)
PHASES = [180] * K

XPARAMETERS = {
    "meta": {
        "X_NAME": "K=9 EWT prism soft boost=0.00001 V_MODE=10 pressure=0.0001",
        "DESCRIPTION": "K-selectivity pressure test: K=9 EWT geometry with vacuum pressure force",
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
        "SEED_BOOST": 0.01,
        "V_MODE": 10,
        "V_C1": -0.3,
        "V_C2": 0.1,
        "WC_INTERACT_MODE": 3,
        "WC_BOOST": 0.00001,
        "WC_RADIUS": 2,
        "WC_SIGMA": 2.0,
        "VELOCITY_DAMPING": 0.95,
        "R_WALL": 100.0,
        "WALL_HEIGHT": 1.2,
        "DEFICIT_DEPTH": 0.9,
        "R_SOLITON": 35.0,
        "SIGMA": 3.0,
        "PRESSURE_STRENGTH": 0.0001,
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
        "PARTICLE_SHELL": False,
        "SIM_SPEED": 1.0,
        "PAUSED": False,
    },
    "color_defaults": {"COLOR_THEME": "OCEAN", "WAVE_MENU": 1},
    "analytics": {"INSTRUMENTATION": True, "EXPORT_VIDEO": False, "VIDEO_FRAMES": 24},
}