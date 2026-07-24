"""
XPERIMENT PARAMETERS - Nonlinearity test with EMC wall (Variant B + wall)

Goal: test whether the density-modulated cubic nonlinearity with a
      Degraded EMC Wall (v_mode=5) stabilises the K=10 electron
      configuration (1-3-6 tetrahedron) better than neighbouring K.

Variant B + wall: F(psi) = gamma * psi^3 * (1 - rho(r)/rho_0)
where rho(r) is a radial EMC density profile with a core deficit
and a Gaussian wall at r = R_WALL.

Parameters:
  V_MODE = 5   (deficit + wall)
  R_WALL = 100 (wall radius in grid units)
  WALL_HEIGHT = 2.0 (peak density relative to background)
  DEFICIT_DEPTH = 0.9 (core density reduction)
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
        "X_NAME": f"  /Soliton EMC-wall (K={K})",
        "DESCRIPTION": "Nonlinearity test with EMC wall - K=10 (electron)",
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
        "SEED_MODE": 0,
        "SEED_BOOST": 0.01,                
        "V_MODE": 1,                      
        "V_C1": -5.0,                 
        "V_C2": 0.0,
        "WC_INTERACT_MODE": 0,
        "WC_BOOST": 1.0,
        "WC_RADIUS": 2,
        "WC_SIGMA": 1.5,
        "R_WALL": 100.0,
        "WALL_HEIGHT": 1.2,
        "DEFICIT_DEPTH": 0.8,
        "CFL_SAFETY": 1e-7,                # sub-CFL
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
        "WAVE_MENU": 1,          # displacement - best reveals wave structure
    },
    "analytics": {
        "INSTRUMENTATION": True, # enable diagnostics
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
}