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


def generate_K_positions(univ_edge, K, center=(0.5, 0.5, 0.5), rotation=(0, 0, 0), perturbation=0.0):
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


def tetrahedron_10(univ_edge, center=(0.5, 0.5, 0.5), rotation=(0, 0, 0)):
    """Generate 10 positions in 1-3-6 tetrahedral arrangement."""
    LOCK_SPACING = constants.EWAVE_LENGTH / univ_edge
    layer_h = LOCK_SPACING * math.sqrt(2 / 3)
    cx, cy, cz = center

    Rb = LOCK_SPACING * 2 / math.sqrt(3)
    Rm = LOCK_SPACING / math.sqrt(3)

    angles_v = [math.radians(90), math.radians(210), math.radians(330)]
    angles_m = [math.radians(30), math.radians(150), math.radians(270)]

    local_positions = []
    for a in angles_v:
        local_positions.append([Rb * math.cos(a), Rb * math.sin(a), 0.0])
    for a in angles_m:
        local_positions.append([Rm * math.cos(a), Rm * math.sin(a), 0.0])
    for a in angles_v:
        local_positions.append([Rm * math.cos(a), Rm * math.sin(a), layer_h])
    local_positions.append([0.0, 0.0, 2 * layer_h])

    # Apply rotation if needed
    if rotation != (0, 0, 0):
        rx, ry, rz = math.radians(rotation[0]), math.radians(rotation[1]), math.radians(rotation[2])
        Rx = [[1, 0, 0], [0, math.cos(rx), -math.sin(rx)], [0, math.sin(rx), math.cos(rx)]]
        Ry = [[math.cos(ry), 0, math.sin(ry)], [0, 1, 0], [-math.sin(ry), 0, math.cos(ry)]]
        Rz = [[math.cos(rz), -math.sin(rz), 0], [math.sin(rz), math.cos(rz), 0], [0, 0, 1]]

        def mat_vec_mul(M, v):
            return [M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2],
                    M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2],
                    M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2]]

        rotated = []
        for p in local_positions:
            temp = mat_vec_mul(Rx, p)
            temp = mat_vec_mul(Ry, temp)
            temp = mat_vec_mul(Rz, temp)
            rotated.append(temp)
        local_positions = rotated

    return [[cx + p[0], cy + p[1], cz + p[2]] for p in local_positions]


def golden_angle_positions(K, radius, center):
    """Generate K points on a sphere via Fibonacci spiral."""
    points = []
    phi = math.pi * (3.0 - math.sqrt(5.0))
    cx, cy, cz = center

    for i in range(K):
        y = 1.0 - (2.0 * i + 1.0) / K
        r_y = math.sqrt(1.0 - y * y)
        theta = phi * i
        points.append([cx + math.cos(theta) * r_y * radius,
                       cy + y * radius,
                       cz + math.sin(theta) * r_y * radius])
    return points


def bcc_lattice_positions(K, center=(0.5, 0.5, 0.5), spacing=0.04):
    """Generate K positions on a cubic/BCC lattice."""
    cx, cy, cz = center
    offsets = []
    max_radius = int(math.ceil((K ** (1/3)) / 2)) + 1

    for ix in range(-max_radius, max_radius + 1):
        for iy in range(-max_radius, max_radius + 1):
            for iz in range(-max_radius, max_radius + 1):
                if ix == 0 and iy == 0 and iz == 0:
                    continue
                dist = math.sqrt(ix*ix + iy*iy + iz*iz)
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
    return [[
        p[0] + rng.uniform(-perturbation, perturbation) * lock_spacing,
        p[1] + rng.uniform(-perturbation, perturbation) * lock_spacing,
        p[2] + rng.uniform(-perturbation, perturbation) * lock_spacing
    ] for p in positions]

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
        shell_positions = _torus_positions(remaining, shell_radius, torus_thickness_factor, (cx, cy, cz))

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