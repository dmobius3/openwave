# electron_k10_vmode5_c1_m5_boost_005.py
from openwave.xperiments.m4_ewt.xparameters.formation02 import generate_K_positions
from openwave.common import constants

UNIVERSE_EDGE = 1e-15
TARGET_VOXELS = 55_000_000
K = 10
PERTURBATION = 0.1

POSITIONS = generate_K_positions(UNIVERSE_EDGE, K, center=(0.5,0.5,0.5), perturbation=PERTURBATION)
PHASES = [180]*K

XPARAMETERS = {
    "meta": {"X_NAME": "Electron K=10 V_MODE=5 c1=-5 wall=2.0"},
    "camera": {"INITIAL_POSITION": [0.94, 0.91, 0.69]},
    "universe": {"SIZE": [UNIVERSE_EDGE]*3, "TARGET_VOXELS": TARGET_VOXELS},
    "wave_centers": {
        "COUNT": K,
        "POSITION": POSITIONS,
        "PHASE_OFFSETS_DEG": PHASES,
        "APPLY_MOTION": True,
    },
    "engine": {
        "SEED_MODE": 2,
        "SEED_BOOST": 0.05,
        "V_MODE": 5,
        "V_C1": -5.0,              # 2.5× silniejsze
        "V_C2": 0.0,
        "WC_INTERACT_MODE": 0,
        "WC_BOOST": 1.0,
        "WC_RADIUS": 2,
        "WC_SIGMA": 1.5,
        "R_WALL": 100.0,
        "WALL_HEIGHT": 2.0,        # œciana odbijaj¹ca
        "DEFICIT_DEPTH": 0.8,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "FLUX_MESH_PLANES": [0.5,0.5,0.5],
        "SHOW_FLUX_MESH": 1,
        "WARP_MESH": 100,
        "SHOW_GRANULES": False,
        "PARTICLE_SHELL": True,
        "SIM_SPEED": 1.0,
        "PAUSED": False,
    },
    "color_defaults": {"COLOR_THEME": "OCEAN", "WAVE_MENU": 1},
    "analytics": {"INSTRUMENTATION": True, "EXPORT_VIDEO": False, "VIDEO_FRAMES": 24},
}