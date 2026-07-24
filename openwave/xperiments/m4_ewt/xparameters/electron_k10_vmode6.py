"""
XPERIMENT PARAMETERS – K=10 with V_MODE=6 (flat density profile)

V_MODE=6: flat-bottom density deficit + sigmoid wall.
- Constant deficit inside soliton (r < R_SOLITON)
- Rapid rise to vacuum outside (sigmoid transition)
- Gives uniform nonlinearity inside, prevents asymmetric drift
"""

from openwave.xperiments.m4_ewt.xparameters.formation02 import generate_K_positions
from openwave.common import constants

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 55_000_000
K = 10
PERTURBATION = 0.1

# Generate 1-3-6 tetrahedron positions
POSITIONS = generate_K_positions(
    UNIVERSE_EDGE, K, center=(0.5, 0.5, 0.5), rotation=(0, 0, 0), perturbation=PERTURBATION
)
PHASES = [180] * K  # all in-phase (electron)

XPARAMETERS = {
    "meta": {
        "X_NAME": f"Electron K={K} V_MODE=6 flat deficit R=20",
        "DESCRIPTION": "Flat-bottom density profile with sigmoid wall",
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
    "SEED_BOOST": 0.02,          # mniejszy impuls
    "V_MODE": 6,
    "V_C1": -3.0,                # mniejsza nieliniowość
    "V_C2": 0.0,
    "WC_INTERACT_MODE": 0,
    "R_WALL": 100.0,
    "WALL_HEIGHT": 1.0,
    "DEFICIT_DEPTH": 0.9,
    "R_SOLITON": 30.0,           # większy promień – więcej miejsca
    "SIGMA": 10.0,               # łagodniejsze przejście
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