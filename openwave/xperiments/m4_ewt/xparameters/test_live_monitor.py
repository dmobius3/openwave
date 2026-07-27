"""
Minimal test for live monitoring - single Dirichlet-driven wave center.
No nonlinearity, no motion, just one strong source in the center.
"""

UNIVERSE_EDGE = 2e-15
TARGET_VOXELS = 55_000_000

XPARAMETERS = {
    "meta": {
        "X_NAME": "Test Live Monitor (K=1, Dirichlet)",
        "DESCRIPTION": "Single strong WC to verify live plotting works",
    },
    "camera": {"INITIAL_POSITION": [1.40, 1.40, 1.20]},
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],
        "TARGET_VOXELS": TARGET_VOXELS,
    },
    "wave_centers": {
        "COUNT": 1,
        "POSITION": [[0.50, 0.50, 0.50]],
        "PHASE_OFFSETS_DEG": [0],
        "APPLY_MOTION": False,
    },
    "engine": {
        "SEED_MODE": 2,
        "SEED_BOOST": 1.0,
        "V_MODE": 0,
        "V_C1": 0.0,
        "V_C2": 0.0,
        "WC_INTERACT_MODE": 1,
        "WC_BOOST": 1.0,
        "WC_RADIUS": 2,
        "WC_SIGMA": 1.5,
        "VELOCITY_DAMPING": 0.999,
        "R_WALL": 100.0,
        "WALL_HEIGHT": 1.2,
        "DEFICIT_DEPTH": 0.9,
        "R_SOLITON": 35.0,
        "SIGMA": 3.0,
        "PRESSURE_STRENGTH": 0.0,
        "CFL_SAFETY": 0.95,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 1,
        "WARP_MESH": 20,
        "SHOW_GRANULES": False,
        "PARTICLE_SHELL": False,
        "SIM_SPEED": 1.0,
        "PAUSED": False,
    },
    "color_defaults": {"COLOR_THEME": "OCEAN", "WAVE_MENU": 1},
    "analytics": {"INSTRUMENTATION": True, "EXPORT_VIDEO": False, "VIDEO_FRAMES": 24},
}