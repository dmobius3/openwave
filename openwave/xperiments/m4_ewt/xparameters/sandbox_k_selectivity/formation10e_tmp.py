"""
XPERIMENT PARAMETERS

Progressive particle formation test.
Verifies prediction: K=2..9 are UNSTABLE (decay/fly apart),
K=10 is the first stable standalone particle (electron tetrahedron).

Switch K below. Positions come from xparameters/utils/geometry.py, which
dispatches on K:
  K=10:  the 1-3-6 tetrahedron (the electron)
  K=11:  a golden-angle (Fibonacci) sphere
  other: the golden-angle fallback, every point on a sphere of radius 0.35 lambda

Spacing caveat: only K=10 spans the lock-in wells at r = n*lambda, the energy
minima where same-phase wave centres hold. Under the fallback all pair
separations land between 0.33 and 0.70 lambda, inside the first well, so K=2..9
are not tested at the lock-in spacing. The named geometries this file used to
build (line, triangle, tetrahedron, bipyramid, octahedron, cube, tricapped
prism) are no longer generated.
"""

from openwave.xperiments.m4_ewt.xparameters.utils.geometry import (
    generate_positions_by_EWT_geometry,
)

UNIVERSE_EDGE = 2e-15  # m, universe edge length in meters
TARGET_VOXELS = 55_000_000  # Target voxel count (impacts performance)

# ════════════════════════════════════════════════════════════════════════════
# SELECT K VALUE HERE. K=10 is the 1-3-6 tetrahedron, K=11 a golden-angle
# sphere, every other K the golden-angle fallback.
# ════════════════════════════════════════════════════════════════════════════
K = 10

# Perturbation: shift each WC by random ±PERTURBATION fraction of λ.
# At 0.0: perfect lattice (all K stable). At 0.2+: real test.
PERTURBATION = 0.1  # fraction of λ (0.0 = perfect, 0.3 = 30% random displacement)

POSITIONS = generate_positions_by_EWT_geometry(
    UNIVERSE_EDGE, K, center=(0.5, 0.5, 0.5), perturbation=PERTURBATION
)
PHASES = [180] * K  # all same phase (electron-like)

XPARAMETERS = {
    "meta": {
        "X_NAME": f"  /Electron (K={K})",
        "DESCRIPTION": f"K={K} stability test — {'STABLE' if K == 10 else 'expect UNSTABLE'}",
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
        # Base wave seed (P1)
        "SEED_MODE": 2,  # 0 = gaussian pulse, 1 = radial cosine, 2 = full (domain-filling base wave)
        "SEED_BOOST": 0.05,  # seed amplitude multiplier
        # Non-linear potential V(ψ) (P2)
        "V_MODE": 7,  # 0 = linear/off, 1 = cubic ψ³, 2 = saturating, 3 = double-well
        "V_C1": 0,
        "V_C2": 1,
        # Wave-center interaction (P3)
        "WC_INTERACT_MODE": 3,  # 0 = free, 1 = dirichlet, 2 = neumann, 3 = soft
        "WC_BOOST": 0.05,  # WC drive amplitude multiplier
        "WC_RADIUS": 2,  # WC drive ball radius (voxels)
        "WC_SIGMA": 1.5,  # soft-mode Gaussian width (voxels)
        "PRESSURE_STRENGTH": 0.001,
        "CFL_SAFETY": 0.1,
        "R_SOLITON": 35,
        "DEFICIT_DEPTH": 0.9,
        "VELOCITY_DAMPING": 0.90,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 1,
        "WARP_MESH": 30,
        "SHOW_GRANULES": False,  # Toggle to show/hide granule particles (rendered as points)
        "TIMESTEP": 5.0,
        "PAUSED": False,
        "PARTICLE_SHELL": True,
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
