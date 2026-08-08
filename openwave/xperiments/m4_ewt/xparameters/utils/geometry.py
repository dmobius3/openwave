"""
GEOMETRY GENERATORS FOR WAVE CENTER POSITIONS

This module provides functions to generate initial positions for Wave Centers (WCs)
in various geometric configurations. Used by xparameters files to define particle
topologies for M4 simulations.

Supported geometries:
    - 1-3-6 tetrahedron (electron, K=10)
    - Golden angle phyllotaxis (minimal interference, K=11)
    - BCC lattice (cubic grid, for testing non-optimal configurations)
"""

import math
import random
from openwave.common import constants


def generate_K_positions(
    univ_edge, K, center=(0.5, 0.5, 0.5), rotation=(0, 0, 0), perturbation=0.0
):
    """Main dispatcher for K-position generation."""
    LOCK_SPACING = constants.EWAVE_LENGTH / univ_edge

    if K == 10:
        positions = tetrahedron_10(univ_edge, center=center, rotation=rotation)
    elif K == 11:
        radius = 0.35 * LOCK_SPACING
        positions = golden_angle_positions(K, radius, center)
    else:
        positions = _generic_positions(K, center, LOCK_SPACING)

    if perturbation > 0:
        rng = random.Random(42)
        positions = _apply_perturbation(positions, perturbation, LOCK_SPACING, rng)

    return positions


def tetrahedron_10(univ_edge, center=(0.5, 0.5, 0.5), rotation=(0, 0, 0), perturbation=0.0):
    """
    Generate proper 1-3-6 tetrahedral electron geometry.
    Returns 10 positions: 1 center, 3 inner, 6 outer.

    Scaling note (maintainer, measured): r1 and r2 below are normalized
    constants, not multiples of LOCK_SPACING, so this geometry does NOT scale
    with univ_edge even though it takes it as an argument. Measured at
    univ_edge = 1e-15, none of the 45 pair separations land on a lock-in well
    r = n*lambda (they run 0.701 to 3.336 lambda); at 2e-15 the same
    normalized coordinates come out doubled in units of lambda. Which radii the
    1-3-6 electron should use is the model author's call, being explored across
    the electron_k*_vmode10_* xparameters, so this is recorded rather than
    changed. LOCK_SPACING is currently computed and unused.
    """
    import math
    import random

    cx, cy, cz = center
    # Normalized wavelength for scaling
    LOCK_SPACING = constants.EWAVE_LENGTH / univ_edge
    # Radii: inner ~0.015-0.02, outer ~0.03-0.04 (in normalized units)
    # NOTE: normalized constants, NOT derived from LOCK_SPACING; see the
    # scaling note in this function's docstring before reading a K=10 run
    # as a lock-in result.
    r1 = 0.02  # inner radius
    r2 = 0.04  # outer radius
    h = r2 * math.sqrt(2 / 3)  # height offset for outer layer (approx 0.03266)

    positions = []
    # 1. Center
    positions.append([cx, cy, cz])

    # 2. Inner 3: in XY plane, Z=0, 120° apart
    angles_inner = [math.radians(90), math.radians(210), math.radians(330)]
    for a in angles_inner:
        positions.append([cx + r1 * math.cos(a), cy + r1 * math.sin(a), cz])

    # 3. Outer 6: two layers of 3, with Z offset, rotated 60° relative to inner
    angles_outer = [math.radians(30), math.radians(150), math.radians(270)]
    # Lower layer (Z = -h)
    for a in angles_outer:
        positions.append([cx + r2 * math.cos(a), cy + r2 * math.sin(a), cz - h])
    # Upper layer (Z = +h)
    for a in angles_outer:
        positions.append([cx + r2 * math.cos(a), cy + r2 * math.sin(a), cz + h])

    # Optional rotation (currently no rotation implemented, but could be added)
    # For simplicity, we skip rotation; rotation can be applied later.

    # Perturbation
    if perturbation > 0:
        rng = random.Random(42)
        positions = _apply_perturbation(positions, perturbation, LOCK_SPACING, rng)

    return positions


def tetrahedron_10_locked(
    univ_edge,
    center=(0.5, 0.5, 0.5),
    rotation=(0, 0, 0),
    perturbation=0.0,
):
    """
    Generate the 1-3-6 tetrahedral electron geometry locked on n·λ wells.

    Positions are scaled by LOCK_SPACING so that all separations are
    invariant with respect to the universe edge.

    Geometry:
        - The inner 3 WCs sit at exactly 1.000·λ from the centre.
        - The outer 6 WCs sit at exactly 2.000·λ from the centre (r2 is
          derived from √(r2² + h²) = 2λ, with h = r2·√(2/3), giving
          r2 = 2λ·√(3/5) ≈ 1.549λ).

    Parameters
    ----------
    univ_edge : float
        Universe edge length in metres.
    center : tuple[float, float, float]
        Normalised (cx, cy, cz) centre of the geometry.
    rotation : tuple[float, float, float]
        Euler angles in degrees (Z, Y, X) applied after construction.
    perturbation : float
        Random displacement fraction of LOCK_SPACING (0 = exact positions).
    """
    import math
    import random

    cx, cy, cz = center
    LOCK_SPACING = constants.EWAVE_LENGTH / univ_edge

    # Radii in normalised units, chosen as multiples of the lock-in spacing
    r1 = 1.0 * LOCK_SPACING  # inner shell: exactly 1·λ
    # outer shell: exactly 2·λ from centre (solved from √(r2² + h²) = 2λ with h = r2·√(2/3))
    r2 = 2.0 * math.sqrt(3.0 / 5.0) * LOCK_SPACING
    h = r2 * math.sqrt(2.0 / 3.0)  # vertical offset for the outer layers

    positions = []

    # 1. Centre
    positions.append([cx, cy, cz])

    # 2. Inner 3 – equilateral triangle in the XY plane
    angles_inner = [math.radians(90), math.radians(210), math.radians(330)]
    for a in angles_inner:
        positions.append(
            [
                cx + r1 * math.cos(a),
                cy + r1 * math.sin(a),
                cz,
            ]
        )

    # 3. Outer 6 – two layers of 3, rotated 60° relative to the inner triangle
    angles_outer = [math.radians(30), math.radians(150), math.radians(270)]
    # Lower layer (Z = -h)
    for a in angles_outer:
        positions.append(
            [
                cx + r2 * math.cos(a),
                cy + r2 * math.sin(a),
                cz - h,
            ]
        )
    # Upper layer (Z = +h)
    for a in angles_outer:
        positions.append(
            [
                cx + r2 * math.cos(a),
                cy + r2 * math.sin(a),
                cz + h,
            ]
        )

    # ---- Apply rotation ----
    if rotation != (0, 0, 0):
        rz, ry, rx = [math.radians(a) for a in rotation]
        cos_z, sin_z = math.cos(rz), math.sin(rz)
        cos_y, sin_y = math.cos(ry), math.sin(ry)
        cos_x, sin_x = math.cos(rx), math.sin(rx)

        rotated = []
        for x, y, z in positions:
            dx, dy, dz = x - cx, y - cy, z - cz
            # Z
            x1 = dx * cos_z - dy * sin_z
            y1 = dx * sin_z + dy * cos_z
            z1 = dz
            # Y
            x2 = x1 * cos_y + z1 * sin_y
            y2 = y1
            z2 = -x1 * sin_y + z1 * cos_y
            # X
            x3 = x2
            y3 = y2 * cos_x - z2 * sin_x
            z3 = y2 * sin_x + z2 * cos_x
            rotated.append([cx + x3, cy + y3, cz + z3])
        positions = rotated

    # ---- Apply perturbation ----
    if perturbation > 0:
        rng = random.Random(42)
        positions = _apply_perturbation(positions, perturbation, LOCK_SPACING, rng)

    return positions


def generate_positions_by_EWT_geometry_locked(
    univ_edge: float,
    K: int,
    center=(0.5, 0.5, 0.5),
    rotation=(0, 0, 0),
    perturbation: float = 0.0,
):
    """
    Generate K WC positions using EWT geometry locked on n·λ wells.

    For K = 2..9 this delegates to generate_positions_by_EWT_geometry,
    which returns bit-identical points.  For K = 10 it uses
    tetrahedron_10_locked instead of the unscaled tetrahedron_10.

    Parameters
    ----------
    univ_edge : float
        Universe edge length in metres.
    K : int
        Number of wave centres (2..10).
    center : tuple
        Normalised centre.
    rotation : tuple
        Euler angles in degrees (only used for K=10).
    perturbation : float
        Random displacement fraction of LOCK_SPACING.

    Returns
    -------
    list of [x, y, z] normalised positions.
    """
    if K == 10:
        return tetrahedron_10_locked(
            univ_edge, center=center, rotation=rotation, perturbation=perturbation
        )
    # For all other K, delegate to the existing, verified generator
    return generate_positions_by_EWT_geometry(
        univ_edge, K, center=center, rotation=rotation, perturbation=perturbation
    )


def golden_angle_positions(K, radius, center):
    """Generate K points on a sphere via Fibonacci spiral."""
    points = []
    phi = math.pi * (3.0 - math.sqrt(5.0))
    cx, cy, cz = center

    for i in range(K):
        y = 1.0 - (2.0 * i + 1.0) / K
        r_y = math.sqrt(1.0 - y * y)
        theta = phi * i
        points.append(
            [
                cx + math.cos(theta) * r_y * radius,
                cy + y * radius,
                cz + math.sin(theta) * r_y * radius,
            ]
        )
    return points


def bcc_lattice_positions(K, center=(0.5, 0.5, 0.5), spacing=0.04):
    """Generate K positions on a cubic/BCC lattice."""
    cx, cy, cz = center
    offsets = []
    max_radius = int(math.ceil((K ** (1 / 3)) / 2)) + 1

    for ix in range(-max_radius, max_radius + 1):
        for iy in range(-max_radius, max_radius + 1):
            for iz in range(-max_radius, max_radius + 1):
                if ix == 0 and iy == 0 and iz == 0:
                    continue
                dist = math.sqrt(ix * ix + iy * iy + iz * iz)
                offsets.append((ix, iy, iz, dist))

    offsets.sort(key=lambda x: x[3])

    positions = []
    for ix, iy, iz, _ in offsets[:K]:
        positions.append([cx + ix * spacing, cy + iy * spacing, cz + iz * spacing])
    return positions


def _generic_positions(K, center, lock_spacing):
    """Fallback for arbitrary K."""
    return golden_angle_positions(K, 0.35 * lock_spacing, center)


def _apply_perturbation(positions, perturbation, lock_spacing, rng):
    """Apply random perturbation."""
    return [
        [
            p[0] + rng.uniform(-perturbation, perturbation) * lock_spacing,
            p[1] + rng.uniform(-perturbation, perturbation) * lock_spacing,
            p[2] + rng.uniform(-perturbation, perturbation) * lock_spacing,
        ]
        for p in positions
    ]


def generate_shell_structure(
    K,
    univ_edge,
    core_rotation=(0, 0, 0),
    shell_mode="sphere",  # "sphere" or "torus"
    shell_radius_factor=1.5,  # relative to core radius
    torus_thickness_factor=0.3,  # for torus: thickness / radius
    num_shells=1,
    perturbation=0.0,
):
    """
    Generate a core + shell structure for higher K (muon, tau, etc.).

    Core is always the 1-3-6 tetrahedron (10 wave centers).
    Shells are placed on either a sphere (golden angle) or a torus.

    Args:
        K: Total number of wave centers (core + shells)
        univ_edge: Universe edge length (m) – used to compute λ spacing
        core_rotation: Rotation of the core (rx, ry, rz) in degrees
        shell_mode: "sphere" or "torus" – geometry of the shell
        shell_radius_factor: Radius of shell relative to core radius
        torus_thickness_factor: For torus: thickness = radius * factor
        num_shells: Number of shells (1 or 2)
        perturbation: Random displacement fraction of λ

    Returns:
        List of K positions [(x1,y1,z1), ...]
    """
    LOCK_SPACING = constants.EWAVE_LENGTH / univ_edge
    cx, cy, cz = 0.5, 0.5, 0.5

    # 1. Generate core (1-3-6 tetrahedron)
    core_positions = tetrahedron_10(univ_edge, center=(cx, cy, cz), rotation=core_rotation)
    core_size = len(core_positions)

    if K <= core_size:
        return core_positions[:K]

    # 2. Determine shell size
    remaining = K - core_size

    # Core radius (approximate distance from center to outer vertices)
    core_radius = LOCK_SPACING * 2.0 / math.sqrt(3)

    # 3. Generate shell positions
    shell_positions = []
    shell_radius = shell_radius_factor * core_radius

    if shell_mode == "sphere":
        # Distribute points on sphere using golden angle
        points = golden_angle_positions(remaining, shell_radius, (cx, cy, cz))
        shell_positions = points

    elif shell_mode == "torus":
        # Distribute points on a torus surface
        shell_positions = _torus_positions(
            remaining, shell_radius, torus_thickness_factor, (cx, cy, cz)
        )

    elif shell_mode == "flat_disk":
        # Flat disk (2D) – useful for testing planar configurations
        shell_positions = _disk_positions(remaining, shell_radius, (cx, cy, cz))

    else:
        raise ValueError(f"Unknown shell_mode: {shell_mode}")

    # 4. Apply perturbation (if any)
    if perturbation > 0:
        rng = random.Random(42)
        shell_positions = _apply_perturbation(shell_positions, perturbation, LOCK_SPACING, rng)

    # 5. Combine core + shells
    return core_positions + shell_positions


def _torus_positions(N, radius, thickness_factor, center):
    """
    Generate N points on a torus surface.

    The torus is oriented with its axis along Z.
    Points are distributed using golden angle on the torus surface.

    Args:
        N: Number of points
        radius: Major radius (distance from center to tube center)
        thickness_factor: Minor radius = radius * thickness_factor
        center: (cx, cy, cz) center of the torus

    Returns:
        List of N positions
    """
    cx, cy, cz = center
    minor_radius = radius * thickness_factor

    if N == 0:
        return []

    # Golden angle for even distribution on torus
    phi = math.pi * (3.0 - math.sqrt(5.0))

    positions = []
    for i in range(N):
        # Two angles: theta (around torus) and phi (around tube)
        theta = 2.0 * math.pi * i / N
        psi = phi * i

        # Torus surface parametrization
        x = (radius + minor_radius * math.cos(psi)) * math.cos(theta)
        y = (radius + minor_radius * math.cos(psi)) * math.sin(theta)
        z = minor_radius * math.sin(psi)

        positions.append([cx + x, cy + y, cz + z])

    return positions


def _disk_positions(N, radius, center):
    """
    Generate N points on a flat disk (2D) using golden angle spiral.

    This is useful for testing planar configurations.

    Args:
        N: Number of points
        radius: Disk radius
        center: (cx, cy, cz) center of the disk

    Returns:
        List of N positions (all with Z = center_z)
    """
    cx, cy, cz = center

    if N == 0:
        return []

    positions = []
    # Use Fermat spiral for disk distribution
    golden_ratio = (1.0 + math.sqrt(5.0)) / 2.0

    for i in range(N):
        # Fermat spiral: r = sqrt(i / N) * radius, theta = i * golden_angle
        r = math.sqrt(i / N) * radius
        theta = 2.0 * math.pi * i / golden_ratio

        x = cx + r * math.cos(theta)
        y = cy + r * math.sin(theta)
        z = cz

        positions.append([x, y, z])

    return positions

    # ================================================================


# EWT LEGACY GEOMETRIES (from formation02.py)
# ================================================================


def tricapped_trigonal_prism_positions(
    K: int, center=(0.5, 0.5, 0.5), lock_spacing: float = 0.05, perturbation: float = 0.0
):
    """
    Generate K=9 positions in a tricapped trigonal prism geometry.

    Structure: 6 points form a trigonal prism (3 top, 3 bottom),
    plus 3 caps at the side-face centers.

    Args:
        K: Number of wave centres (must be 9)
        center: (cx, cy, cz) center in normalized coords
        lock_spacing: Characteristic spacing (λ in normalized coords)
        perturbation: Random displacement fraction of lock_spacing

    Returns:
        List of K positions [[x1,y1,z1], ...]
    """
    import math
    import random

    if K != 9:
        raise ValueError(f"tricapped_trigonal_prism_positions only supports K=9, got K={K}")

    cx, cy, cz = center
    s = lock_spacing

    angles = [math.radians(90), math.radians(210), math.radians(330)]
    r = s / math.sqrt(3)  # circumradius for edge = s
    h = s / 2  # half-height of prism

    positions = []

    # Bottom and top triangles
    for z_off in [-h, h]:
        for a in angles:
            positions.append([cx + r * math.cos(a), cy + r * math.sin(a), cz + z_off])

    # Caps: 3 points at the side-face centers (rotated by 60°)
    cap_angles = [math.radians(30), math.radians(150), math.radians(270)]
    cap_r = r * 1.5
    for a in cap_angles:
        positions.append([cx + cap_r * math.cos(a), cy + cap_r * math.sin(a), cz])

    if perturbation > 0:
        rng = random.Random(42)
        positions = _apply_perturbation(positions, perturbation, lock_spacing, rng)

    return positions


def generate_positions_by_EWT_geometry(
    univ_edge: float,
    K: int,
    center=(0.5, 0.5, 0.5),
    rotation=(0, 0, 0),
    perturbation: float = 0.0,
):
    """
    Generate K WC positions using the original EWT geometry definitions.

    This is the "legacy" geometry generator from formation02.py.
    It uses specific polyhedral geometries for each K:
        K=2:  line
        K=3:  equilateral triangle
        K=4:  regular tetrahedron
        K=5:  trigonal bipyramid
        K=6:  octahedron
        K=7:  pentagonal bipyramid
        K=8:  cube (dual tetrahedra)
        K=9:  tricapped trigonal prism
        K=10: 1-3-6 tetrahedron (electron)

    For K outside 2..10, falls back to golden_angle_positions.

    Args:
        univ_edge: Universe edge length (m) – used to compute λ spacing
        K: Number of wave centres
        center: (cx, cy, cz) center in normalized coords
        rotation: (rx, ry, rz) rotation in degrees (only used for K=10)
        perturbation: Random displacement fraction of λ

    Returns:
        List of K positions [[x1,y1,z1], ...]
    """
    import math
    import random

    LOCK_SPACING = constants.EWAVE_LENGTH / univ_edge
    cx, cy, cz = center
    s = LOCK_SPACING

    if K == 2:
        positions = [[cx - s / 2, cy, cz], [cx + s / 2, cy, cz]]

    elif K == 3:
        angles = [math.radians(90), math.radians(210), math.radians(330)]
        r = s / math.sqrt(3)
        positions = [[cx + r * math.cos(a), cy + r * math.sin(a), cz] for a in angles]

    elif K == 4:
        h = s * math.sqrt(2 / 3)
        r = s / math.sqrt(3)
        positions = [
            [cx, cy + r, cz],
            [cx - s / 2, cy - r / 2, cz],
            [cx + s / 2, cy - r / 2, cz],
            [cx, cy, cz + h],
        ]

    elif K == 5:
        angles = [math.radians(90), math.radians(210), math.radians(330)]
        r = s / math.sqrt(3)
        positions = [[cx + r * math.cos(a), cy + r * math.sin(a), cz] for a in angles]
        positions.append([cx, cy, cz + s / 2])
        positions.append([cx, cy, cz - s / 2])

    elif K == 6:
        d = s / math.sqrt(2)
        positions = [
            [cx + d, cy, cz],
            [cx - d, cy, cz],
            [cx, cy + d, cz],
            [cx, cy - d, cz],
            [cx, cy, cz + d],
            [cx, cy, cz - d],
        ]

    elif K == 7:
        angles = [math.radians(i * 72) for i in range(5)]
        r = s / (2 * math.sin(math.pi / 5))
        positions = [[cx + r * math.cos(a), cy + r * math.sin(a), cz] for a in angles]
        positions.append([cx, cy, cz + s / 2])
        positions.append([cx, cy, cz - s / 2])

    elif K == 8:
        d = s / 2
        positions = [
            [cx + d, cy + d, cz + d],
            [cx + d, cy + d, cz - d],
            [cx + d, cy - d, cz + d],
            [cx + d, cy - d, cz - d],
            [cx - d, cy + d, cz + d],
            [cx - d, cy + d, cz - d],
            [cx - d, cy - d, cz + d],
            [cx - d, cy - d, cz - d],
        ]

    elif K == 9:
        positions = tricapped_trigonal_prism_positions(K, center, LOCK_SPACING, perturbation)
        # perturbation already applied inside helper
        return positions

    elif K == 10:
        positions = tetrahedron_10(univ_edge, center=center, rotation=rotation)

    else:
        # Fallback: golden angle on a small sphere
        positions = golden_angle_positions(K, 0.35 * LOCK_SPACING, center)

    # Apply perturbation for all cases except K=9 (already handled)
    if perturbation > 0 and K != 9:
        rng = random.Random(42)
        positions = _apply_perturbation(positions, perturbation, LOCK_SPACING, rng)

    return positions
