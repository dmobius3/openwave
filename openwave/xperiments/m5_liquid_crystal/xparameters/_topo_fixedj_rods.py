"""
XPERIMENT PARAMETERS — the J/mu TWIST DEMO: disclination rods under the live
fixed-J clock (M5.23.2 arm (2), on the M5.23.1 production port).

The angular-momentum / magnetic-dipole demonstration from simulation: the
VIZ.5 Stage D rod render (rod-line + rod-ring ellipsoid samples, shapes read
LIVE from the M.u map) over the certified fixed-J isorotation flow
(RELAX -> SET-J -> EVOLVE). As the carried clock rotates the core's delta
frame, the rod/ring ellipsoid SHAPES twist on screen — simulated dynamics
only (the standing no-display-only-kinematics directive): the rotation you
see IS the evolved field, not an animation.

The visible axis rate is omega*/|a0_raw| (M5.23.1 measured ~0.016 rad/tau at
omega = 0.2 on this arena family), NOT omega* itself — the unit-global-norm
a0 convention spreads the rotation over the flow support. Raise ETA_SUBSTEPS
for a faster visual.

Sample positions are the Stage D analytic layout (centered rod along the
spin axis); the sampled SHAPES are live field reads. Tracer-driven rod
PLACEMENT is the staged follow-up (m5_23_2 task_details, deferred).

The energy isosurface (M5.23.2 arm (4)) is available here too: check
"Iso-Surface (energyH)" — the electron's surface renders as the rod TUBE
(measured m5_23_2: the density concentrates along the boundary-to-boundary
rod; the tube is open at the box, that is the physics).
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # 63^3 after odd-rounding

TOPOLOGY_SEED = {
    "MODE": "biaxial_hedgehog",
    "CENTER": [0.50, 0.50, 0.50],
    "R0_FRACTION": 0.06,
    "RHOC_VOXELS": 3.0,
    "BIAXIAL_DELTA": 0.30,
    "AUTO_RELAX_STEPS": 0,
    "INTEGRATOR_4D": "canonical",  # M5.24 verified-L stack
    "DT_ETA_CAP": 0.005,  # certified full-3D tau-step
    "ETA_DX": 1.5,  # research-twin geometry
    "ETA_SUBSTEPS": 64,
    "ETA_SPONGE_GAMMA": 0.5,
    "ETA_SPONGE_WIDTH": 10.0,
    "CANON_RELAX_ITERS": 300,  # condition the seed BEFORE the kick
    "FIXEDJ_OMEGA": 0.2,  # the certified research rung
    "FIXEDJ_RENV": 10.0,
    "FIXEDJ_LOG_EVERY": 30,
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "J/mu twist demo (rods + clock)",
        "DESCRIPTION": "Disclination rods twisting under the live fixed-J ZBW clock (M5.23.2 arm 2)",
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
        "SHOW_GLYPHS": 0,  # the ellipsoid views replace the plane glyphs
        "SHOW_ELLIPSOID": True,  # VIZ.5 on at boot: the twist is the show
        "ELLIPSOID_SHELL": False,  # rods + rings only (shell occludes the twist)
        "ELLIPSOID_RODS": True,
        "ELLIPSOID_RODRINGS": True,
        "ELLIPSOID_RADIUS": 0.1,
        "ELLIPSOID_SIZE": 0.025,
        "ELLIPSOID_COUNT": 299,
        "SHOW_ISO": False,  # arm (4): toggle for the energy rod-tube surface
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
        "WAVE_MENU": 4,  # Hamiltonian energy
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
