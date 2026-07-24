"""
XPERIMENT PARAMETERS — the MU-CANDIDATE from disk (M5.23.2 arm (3)): the
M5.21.6 level-B census state loaded npz -> launcher seed.

The author's 2026-07-20 ask ("would be great to visualize them, also for
muon and taon") on the on-disk data: `m5_21_6_end_f64_B.npz` (64^3 3x3,
free bc, h = 1.0, the B excited level of the biaxial-hedgehog census with
its M_it* convergence history). 63^3 launcher arena = CROP fit. The 3x3
spatial block takes the launcher's own 4D embed + covariant flip; ETA_DX
defaults to the npz h (= 1.0).

The stub's originally-named 32^3 pinned endpoints (m5_21_6_end_p32_B/C)
are MISSING locally (pre-2026-07-20 delete-rule casualty, m5_23_2 Issue
list); the f64 states are the on-disk μ/τ pair of record AND free-bc, the
launcher-compatible boundary. No clock is set: no certified fixed-J rung
exists for the mu state (honest display: statics + free evolution only).
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # 63^3 after odd-rounding -> CROP fit of the 64^3 data

TOPOLOGY_SEED = {
    "MODE": "npz_file",
    "PATH": "research/data/m5_21_6_end_f64_B.npz",
    "BIAXIAL_DELTA": 0.30,
    "AUTO_RELAX_STEPS": 0,
    "INTEGRATOR_4D": "canonical",
    "DT_ETA_CAP": 0.005,
    "ETA_SUBSTEPS": 64,
    "ETA_SPONGE_GAMMA": 0.5,
    "ETA_SPONGE_WIDTH": 10.0,
    "CANON_RELAX_ITERS": 0,
    "FIXEDJ_OMEGA": 0.0,  # no certified clock rung for mu — statics honest
    "FIXEDJ_RENV": 10.0,
    "FIXEDJ_LOG_EVERY": 0,
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "Mu-candidate from disk (census B)",
        "DESCRIPTION": "The M5.21.6 level-B state loaded npz -> launcher (M5.23.2 arm 3)",
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
        "SHOW_ELLIPSOID": True,  # the ellipsoid composition view
        "ELLIPSOID_SHELL": True,
        "ELLIPSOID_RODS": True,
        "ELLIPSOID_RODRINGS": True,
        "ELLIPSOID_RADIUS": 0.1,
        "ELLIPSOID_SIZE": 0.025,
        "ELLIPSOID_COUNT": 299,
        "SHOW_ISO": False,  # arm (4) toggle: energy isosurface
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
