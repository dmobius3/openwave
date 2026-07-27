"""
XPERIMENT PARAMETERS - Nonlinearity test with gamma (Variant A)

Goal: verify whether the cubic nonlinearity with gamma = 1/epsilon_M ~ 2.4148e4
      stabilises the K=10 electron configuration (1-3-6 tetrahedron).

Variant A: F(psi) = gamma * psi^3  (pure cubic self-trapping)
Coefficient gamma follows from BCC geometry: gamma = N_final * pi^3 ~ 2.4148e4

A perturbation of 0.1*lambda is applied on purpose - at perfect placement
all K values are stable.
"""

from openwave.common import constants
from openwave.xperiments.m4_ewt.xparameters.formation02 import generate_K_positions

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 55_000_000

K = 10
PERTURBATION = 0.1  # 10% lambda - essential for the selectivity test

POSITIONS = generate_K_positions(
    UNIVERSE_EDGE, K, center=(0.5, 0.5, 0.5), rotation=(0, 0, 0), perturbation=PERTURBATION
)
PHASES = [180] * K  # all in the same phase (electron)

XPARAMETERS = {
    "meta": {
        "X_NAME": f"  /Soliton gamma-test (K={K})",
        "DESCRIPTION": "Nonlinearity test with gamma - Variant A (pure cubic)",
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
        "SEED_MODE": 0,  # gaussian pulse
        "SEED_BOOST": 10.0,  #
        "V_MODE": 4,  # density-modulated cubic
        "V_C1": 300.0,  # start od 300
        "V_C2": 0.0,
        "WC_INTERACT_MODE": 0,
        "WC_BOOST": 1.0,
        "WC_RADIUS": 2,
        "WC_SIGMA": 1.5,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 1,
        "WARP_MESH": 30,
        "PARTICLE_SHELL": True,
        "TIMESTEP": 5.0,
        "PAUSED": False,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 1,  # energy - best reveals nonlinear effects
    },
    "analytics": {
        "INSTRUMENTATION": True,  # enable diagnostics
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
}
