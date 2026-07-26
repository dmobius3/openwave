"""
XPERIMENT PARAMETERS — the FIXED-J isorotation clock on the CANONICAL stack
(M5.23.1: the production port of the M5.21.9 constraint-carried electron).

The two-stack consensus (research m5_21_9_note.md; conceded by the author
2026-07-19) is that free evolution NEVER spins up the ZBW clock: the clock
exists only as a FIXED-J isorotation state, ω* = J/(2·kin), the standard
rotating-soliton construction (Coleman → Radu-Volkov). This xperiment runs
the launcher flow that carries it live:

    seed (biaxial hedgehog) → flip to the covariant vacuum →
    RELAX (FIRE, CANON_RELAX_ITERS) → SET-J (Ṁ(0) = ω*·a0 on the
    conjugation-tangent clock flow) → EVOLVE (the certified η + V4 leapfrog)

so the rotation on screen is SIMULATED DYNAMICS ONLY (the standing
no-display-only-kinematics directive, m5_visualization.md) — the first
honest live ZBW clock. The δ clock-hand glyph (VIZ.3, middle eigenvector)
is the natural thing to watch.

CONVENTIONS (research-pinned, engine2_pde M5.23.1 section): kin/J printed
at SET-J are FLOW-PARAMETER values on the conjugation tangent (kin = 0.1206
on the certified research state); absolute J / ħ/2 / g claims need the
physical-rate convention (research m5_21_5_note § 5). The hold at THIS
arena (63³ launcher grid vs the 32³ research states) is a measurement, not
an inherited certificate — the m5_23_1 selftest S5 gates it headlessly.

Knobs beyond _topo_canonical4d: FIXEDJ_OMEGA (the target clock rate; 0
disables the seed-time kick, the SET J button still works), FIXEDJ_RENV
(the clock-flow envelope radius, research units; the research value 10),
FIXEDJ_LOG_EVERY (periodic console J readout, frames; 0 = off).
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # ~262k voxels (63³ after odd-rounding)

TOPOLOGY_SEED = {
    "MODE": "biaxial_hedgehog",
    "CENTER": [0.50, 0.50, 0.50],
    "R0_FRACTION": 0.06,
    "RHOC_VOXELS": 3.0,
    "BIAXIAL_DELTA": 0.30,
    "AUTO_RELAX_STEPS": 0,
    "INTEGRATOR_4D": "canonical",  # M5.24 — the verified-L era stack
    "DT_ETA_CAP": 0.005,  # certified full-3D τ-step (canonical § 3)
    "ETA_DX": 1.5,  # research grid unit: research-twin geometry
    "ETA_SUBSTEPS": 64,  # physics steps per rendered frame (viz speed)
    "ETA_SPONGE_GAMMA": 0.5,  # boundary damping (absorbs radiation)
    "ETA_SPONGE_WIDTH": 10.0,
    "CANON_RELAX_ITERS": 300,  # seed-time FIRE: condition the seed BEFORE the kick
    "FIXEDJ_OMEGA": 0.2,  # M5.23.1 — the target clock rate ω (research rung 1)
    "FIXEDJ_RENV": 10.0,  # clock-flow envelope radius (research units)
    "FIXEDJ_LOG_EVERY": 30,  # console J readout every N frames (0 = off)
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "Hedgehog FIXED-J clock",
        "DESCRIPTION": "The live ZBW isorotation clock: RELAX -> SET-J -> EVOLVE on the canonical stack (M5.23.1)",
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
        "SHOW_GLYPHS": 3,  # all glyph planes — the δ clock-hand axis is the rotating observable
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 2,
        "WARP_MESH": 0,
        "SHOW_GRANULES": False,
        "SIM_SPEED": 1.0,
        "PAUSED": True,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 4,  # Hamiltonian energy — the verified-H view
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
