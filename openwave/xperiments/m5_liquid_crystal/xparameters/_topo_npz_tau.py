"""
XPERIMENT PARAMETERS — the TAU-CANDIDATE from disk (M5.23.2 arm (3)): the
M5.21.6 level-C census state loaded npz -> launcher seed.

Same route as `_topo_npz_mu` on `m5_21_6_end_f64_C.npz` (64^3 3x3, free
bc, h = 1.0, the C excited level with its M_it* history — the state whose
free decay produced the M5.21.10 symmetric ejection pair; the tracer's
first physics target reads those snapshots directly). See `_topo_npz_mu`'s
header for the fit/convention notes and the missing-p32 provenance.

No clock is set (no certified fixed-J rung for tau — statics honest).
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # 63^3 after odd-rounding -> CROP fit of the 64^3 data

TOPOLOGY_SEED = {
    "MODE": "npz_file",
    "PATH": "research/data/m5_21_6_end_f64_C.npz",
    "BIAXIAL_DELTA": 0.30,
    "AUTO_RELAX_STEPS": 0,
    "INTEGRATOR_4D": "canonical",
    "DT_ETA_CAP": 0.005,
    "ETA_SUBSTEPS": 64,
    "ETA_SPONGE_GAMMA": 0.5,
    "ETA_SPONGE_WIDTH": 10.0,
    "CANON_RELAX_ITERS": 0,
    "FIXEDJ_OMEGA": 0.0,
    "FIXEDJ_RENV": 10.0,
    "FIXEDJ_LOG_EVERY": 0,
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "Tau-candidate from disk (census C)",
        "DESCRIPTION": "The M5.21.6 level-C state loaded npz -> launcher (M5.23.2 arm 3)",
    },
    "camera": {
        "INITIAL_POSITION": [1.10, 1.46, 0.81],
    },
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],
        "TARGET_VOXELS": TARGET_VOXELS,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "VIZ_STRIDE": 1,
        "SHOW_GLYPHS": 0,
        "SHOW_ELLIPSOID": True,
        "ELLIPSOID_SHELL": True,
        "ELLIPSOID_RODS": True,
        "ELLIPSOID_RODRINGS": True,
        "ELLIPSOID_RADIUS": 0.1,
        "ELLIPSOID_SIZE": 0.025,
        "ELLIPSOID_COUNT": 299,
        "SHOW_ISO": False,
        "ISO_LEVEL": 0.30,
        "ISO_COVER": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 0,
        "WARP_MESH": 0,
        "SHOW_GRANULES": False,
        "SIM_SPEED": 1.0,
        "PAUSED": True,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 4,
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
