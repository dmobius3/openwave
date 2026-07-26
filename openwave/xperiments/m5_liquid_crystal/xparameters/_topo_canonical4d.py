"""
XPERIMENT PARAMETERS — Biaxial Hedgehog on the CANONICAL 4D stack (M5.24)

The biaxial hedgehog seed evolved by the production port of the dynamics
formulation of record (research m5_theory_canonical.md § 2 row 6): the
canonical kinetic ½‖Ṁ‖²_F + the η-signed curvature (F_ij = [A_i, A_j]_η,
η = diag(−1,1,1,1)) + the universal spectral potential

    V4 = w · Σ_{p=1..4} (Tr_η(M^p) − C_p)² ,   C_p = g^p + 1 + δ^p ,
    w = 7.24e-4 (the locked WSCALE)

on the M5.21.2b symmetrized stencil. At seed time the launcher flips the
time axis to the covariant vacuum M[0,0] = −g (canonical § 1: diag(+g, …)
is NOT a vacuum of V4). The far field diag(1, δ, 0) of THIS seed family is
exactly the covariant vacuum spatial spectrum, so both the curvature and V4
are zero at the boundary: the energy view renders from a true-zero floor
(no v0 subtraction hack).

WHAT REPLACES WHAT (vs _topo_biaxial1_von): the Eq.13 LdG amplitude well
(LDG_STIFFNESS_K) is retired on this path — V4 is the potential of record
(it confines amplitude AND holds the biaxial spectrum, the author-era open
Q7 answered by the spectral form). The plain-commutator Eq.18 flux is
replaced by the η-bracket flux of the verified Lagrangian (M5.18, 15/15
machine checks).

HONEST LABEL (canonical § 2): this stack is the RUNNABLE REGULARIZATION,
"fine for statics-adjacent dynamics questions and films" — the launcher's
use class. It is NOT the true-L evolution (the free-EL IVP is ill-posed).

dt: compute_timestep caps dt_eff = c·dt at DT_ETA_CAP (default 0.005 τ, the
full-3D certified step, canonical § 3; the V4 stiff g-mode ω ≈ 78 at g = 8
bounds dt_eff < 0.026, and the axisym 0.02 proved margin-critical in 3D f32).

ROUND 2 KNOBS: ETA_SUBSTEPS = physics steps per rendered frame (the viz-speed
lever — dt itself is walled by the stiff mode, so speed comes from substeps;
raise it if the motion is still too slow to watch, at GPU cost per frame).
ETA_DX = the grid spacing the kernels use, in research units (1.5 = the
m5_21_2b/3 h): makes the 63³ box research-twin geometry instead of the
physical dx_am, whose ~15.6 value weakens curvature ~100× against V4.
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # ~262k voxels — small for fast smoke-test

TOPOLOGY_SEED = {
    "MODE": "biaxial_hedgehog",
    "CENTER": [0.50, 0.50, 0.50],
    "R0_FRACTION": 0.06,
    "RHOC_VOXELS": 3.0,
    "BIAXIAL_DELTA": 0.30,
    "AUTO_RELAX_STEPS": 0,
    "INTEGRATOR_4D": "canonical",  # M5.24 — the verified-L era stack
    "DT_ETA_CAP": 0.005,  # certified full-3D τ-step (canonical § 3)
    "ETA_DX": 1.5,  # research grid unit (m5_21_2b/3 h): 63³ box = research-twin geometry
    "ETA_SUBSTEPS": 64,  # physics steps per rendered frame (viz speed; dt itself is walled)
    "ETA_SPONGE_GAMMA": 0.5,  # boundary damping γmax (0 = off): absorbs radiation instead of reflecting
    "ETA_SPONGE_WIDTH": 10.0,  # sponge ramp width, cells
    "CANON_RELAX_ITERS": 0,  # seed-time FIRE iters (the RELAX button runs 500-iter chunks live)
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "Hedgehog CANONICAL 4D",
        "DESCRIPTION": "Biaxial hedgehog on the verified-L era stack: eta-curvature + spectral V4 (M5.24)",
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
        "SHOW_GLYPHS": 2,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 2,
        "WARP_MESH": 0,
        "SHOW_GRANULES": False,
        "SIM_SPEED": 1.0,
        "PAUSED": True,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 4,  # Hamiltonian energy — the verified-H view, true-zero vacuum floor
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
