"""
XPERIMENT PARAMETERS — Charged Disclination Ring (STATIC VIEW)
(M5.21.2 — the third defect type, seeded for 3D inspection)

THE HEDGEHOG'S SAME-CHARGE PARTNER. A closed half-disclination cord of
radius a in the equatorial plane, interior smoothly ESCAPED along +z,
far field IDENTICAL to the radial hedgehog: an enclosing sphere reads
unit hedgehog charge (q = 1) — the far sphere cannot tell ring from
point. Source: Alexander, Chen, Matsumoto & Kamien, Rev. Mod. Phys. 84,
497 (2012) § IV.B (disclination loops carry odd hedgehog charge; local
copy in theory/liquid_crystal_defects/, entry in theory/_CITATIONS.md).
Toulouse-Kleman placement: locally the s=1, d=1 vortex-line cell (a
closed cord, half winding around it), globally the s=2, d=0 erizo cell.

WHAT TO LOOK FOR (glyphs + flux planes): the escaped vertical column
through the ring interior (director = +z, no defect at the origin), the
cord circle where the eigenvalues melt isotropic (the ONLY singular
locus), and the radial hedgehog far field. Compare side by side with
the biaxial hedgehog xperiment: same far field, different core.

⚠️ VIEWING ONLY for now — do NOT expect a stable evolution: the
M5.21.2 statics census (research/findings/m5_21_2_census.md) is still
qualifying which core (point vs ring) the functional prefers, and the
task caught a stencil artifact that demanded a re-run. Boots PAUSED;
leave Evolve PDE alone until the census + successor validations land.

A_FRACTION 0.10 puts the cord at ~6.4 voxels of a 64-box (comfortably
resolved vs the 3-voxel melts). CORD_VOXELS / RHOC_VOXELS mirror the
hedgehog defaults.
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # ~262k voxels

TOPOLOGY_SEED = {
    "MODE": "charged_ring",
    "CENTER": [0.50, 0.50, 0.50],
    "A_FRACTION": 0.10,  # ring radius / box edge (~6.4 voxels at 64)
    "CORD_VOXELS": 3.0,  # half-disclination cord melt width
    "RHOC_VOXELS": 3.0,  # z-axis azimuth melt (house style)
    "BIAXIAL_DELTA": 0.30,
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "Charged Ring (static)",
        "DESCRIPTION": "Charged disclination ring — hedgehog-charge loop, static 3D view (M5.21.2)",
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
        "WAVE_MENU": 4,  # Hamiltonian energy — the cord circle lights up
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
