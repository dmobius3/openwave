# electron_k10_vmode1_neumann_boost0015.py
"""
K=10, 1-3-6 geometry, neumann interaction, WC_BOOST=0.015.
Step 3 of K-selectivity sweep - near equilibrium.
"""

from openwave.xperiments.m4_ewt.xparameters.geometry import tetrahedron_10
from openwave.common import constants

UNIVERSE_EDGE = 2e-15
TARGET_VOXELS = 55_000_000
K = 10

POSITIONS = tetrahedron_10(UNIVERSE_EDGE, center=(0.5, 0.5, 0.5))
PHASES = [180] * K

XPARAMETERS = {
    "meta": {
        "X_NAME": "K=10 1-3-6 neumann boost=0.015",
        "DESCRIPTION": "Step 3: near-equilibrium boost",
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
        "V_MODE": 1,
        "V_C1": -0.1,
        "V_C2": 0.1,
        "WC_INTERACT_MODE": 2,
        "WC_BOOST": 0.015,
        "WC_RADIUS": 2,
        "WC_SIGMA": 1.5,
        "VELOCITY_DAMPING": 0.990,
        "R_WALL": 100.0,
        "WALL_HEIGHT": 1.2,
        "DEFICIT_DEPTH": 0.9,
        "R_SOLITON": 35.0,
        "SIGMA": 3.0,
        "PRESSURE_STRENGTH": 0.0,
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
