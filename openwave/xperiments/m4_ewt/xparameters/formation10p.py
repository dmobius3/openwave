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

Spacing note, MEASURED (all 45 pair separations, univ_edge = 1e-15, in units of
lambda). An earlier version of this note claimed only K=10 sits at the lock-in
wells r = n*lambda; that was written without measuring and is backwards:

  K=2..4    every pair exactly 1.000              fully at the first well
  K=5..9    0.76 to 1.73, some pairs on wells     golden-angle fallback
  K=10      0.701 to 3.336, NO pair on a well     the 1-3-6 tetrahedron
  K=11      0.326 to 0.698, NO pair on a well     golden-angle sphere

So K=10 is currently the case furthest from lock-in spacing, not the closest,
and the 0.33-0.70 lambda band the old note attributed to K=2..9 is really K=11's.

Two consequences worth knowing before reading a K=10 run as a lock-in result:
tetrahedron_10 in utils/geometry.py sets its radii as normalized constants
(r1 = 0.02, r2 = 0.04) instead of deriving them from LOCK_SPACING, so (a) none of
its separations land on a well, and (b) the geometry does not scale with
UNIVERSE_EDGE, so at 2e-15 every number above doubles in units of lambda. Which
radii the 1-3-6 electron should use is the model author's call and is being
explored across the electron_k*_vmode10_* xparameters; this note only records
what the shipped generator currently does.

Reproduce: generate_positions_by_EWT_geometry(1e-15, K, center=(0.5, 0.5, 0.5),
rotation=(0, 0, 0), perturbation=0.0), then take math.dist over every pair and
divide by constants.EWAVE_LENGTH / 1e-15.

The named geometries this file used to build (line, triangle, tetrahedron,
bipyramid, octahedron, cube, tricapped prism) are no longer generated.
"""

from openwave.xperiments.m4_ewt.xparameters.utils.geometry import (
    generate_positions_by_EWT_geometry,
)

UNIVERSE_EDGE = 1e-15  # m, universe edge length in meters
TARGET_VOXELS = 75_000_000  # Target voxel count (impacts performance)

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
PHASES = [0] * K  # all same phase (electron-like)

XPARAMETERS = {
    "meta": {
        "X_NAME": f"  /Positron (K={K})",
        "DESCRIPTION": f"K={K} stability test — {'STABLE' if K == 10 else 'expect UNSTABLE'}",
    },
    "camera": {
        "INITIAL_POSITION": [0.29, 1.28, 0.22],  # [x, y, z] in normalized coordinates
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
        "SEED_BOOST": 1.0,  # seed amplitude multiplier
        # Non-linear potential V(ψ) (P2)
        "V_MODE": 0,  # 0 = linear/off, 1 = cubic ψ³, 2 = saturating, 3 = double-well
        "V_C1": 0.0,  # primary coefficient (k for modes 1/2, a for mode 3); c1 < 0 = focusing
        "V_C2": 0.0,  # secondary coefficient (q for mode 2, b for mode 3)
        # Wave-center interaction (P3)
        "WC_INTERACT_MODE": 0,  # 0 = free, 1 = dirichlet, 2 = neumann, 3 = soft
        "WC_BOOST": 1.0,  # WC drive amplitude multiplier
        "WC_RADIUS": 2,  # WC drive ball radius (voxels)
        "WC_SIGMA": 1.5,  # soft-mode Gaussian width (voxels)
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
        "PARTICLE_SHELL": True,
        "TIMESTEP": 5.0,
        "PAUSED": False,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 1,
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
}
