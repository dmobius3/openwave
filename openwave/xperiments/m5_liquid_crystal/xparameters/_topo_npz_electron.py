"""
XPERIMENT PARAMETERS — the CERTIFIED ELECTRON from disk (M5.23.2 arm (3)):
the M5.21.9 fixed-J isorotation endpoint loaded npz -> launcher seed.

Loads `research/data/m5_21_9_fixedj_conj_om0.2_end.npz` (32^3 4x4, the
omega = 0.2 conjugation rung of record, kin = 0.1206) into a 31^3 launcher
arena (CROP fit: the last plane per axis is pin-shell content). The 4x4
data is flipped to the covariant vacuum at load (measured: V4 = 759/cell
raw vs ~0 flipped, the loader convention pin). FIXEDJ_OMEGA re-kicks the
carried clock at the recorded omega* so the state EVOLVES as the certified
rotating electron — the research object itself on screen, not an ansatz.

ETA_DX is intentionally omitted: the loader defaults it to the npz grid
unit (h = 1.5, research-twin geometry).

NOTE the loaded state is an INITIAL CONDITION here: the research arena was
pinned, this arena is free + sponge — drift beyond the certified hold
window is physics of THIS arena, not an inherited certificate.
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 32**3  # 31^3 after odd-rounding -> CROP fit of the 32^3 data

TOPOLOGY_SEED = {
    "MODE": "npz_file",
    "PATH": "research/data/m5_21_9_fixedj_conj_om0.2_end.npz",
    "BIAXIAL_DELTA": 0.30,
    "AUTO_RELAX_STEPS": 0,
    "INTEGRATOR_4D": "canonical",
    "DT_ETA_CAP": 0.005,
    "ETA_SUBSTEPS": 64,
    "ETA_SPONGE_GAMMA": 0.5,
    "ETA_SPONGE_WIDTH": 6.0,  # smaller arena than the 63^3 demos
    "CANON_RELAX_ITERS": 0,  # the state is already converged
    "FIXEDJ_OMEGA": 0.19923,  # the recorded omega* of this endpoint
    "FIXEDJ_RENV": 10.0,
    "FIXEDJ_LOG_EVERY": 30,
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "Electron from disk (fixed-J endpoint)",
        "DESCRIPTION": "The M5.21.9 certified rotating electron loaded npz -> launcher (M5.23.2 arm 3)",
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
        "SHOW_GLYPHS": 3,  # the delta clock-hand axis is the rotating observable
        "SHOW_ISO": True,  # boot straight into the energy rod-tube surface
        "ISO_LEVEL": 0.50,
        "ISO_COVER": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 2,
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
