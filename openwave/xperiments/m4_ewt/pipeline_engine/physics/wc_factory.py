"""
Wave-center geometry factory.

build_wc_state(K, geometry, nx, ny, nz, ...) -> WCState

Pure Python. No Taichi. Testable without any simulation.
Geometries provided:
    single                 K = 1, one center at domain center
    pair                   K = 2, along x-axis, separated by spacing
    line                   K >= 2, evenly spaced along x-axis
    golden                 K >= 1, golden-angle phyllotaxis on a sphere
    tetrahedron_10_locked  K = 10, 1-3-6 with radii = n * spacing
"""

from __future__ import annotations

import math

from .wc_types import WC, WCState


def build_wc_state(
    K: int,
    geometry: str,
    nx: int,
    ny: int,
    nz: int,
    wavelength: float,
    spacing: float | None = None,
    shell_radius: float | None = None,
    phases: list[float] | None = None,
    amplitude: float = 1.0,
) -> WCState:
    """
    Build a WCState for the given K and geometry.

    Args:
        K: number of wave centers
        geometry: one of the supported geometry names
        nx, ny, nz: grid dimensions
        wavelength: float
            Reference length in grid units. No default: the caller
            supplies it, typically from a unit system or an experiment
            configuration.
        shell_radius: for 'golden'; defaults to 0.35 * spacing
        phases: list of phases in radians, one per center;
                defaults to all zeros
        amplitude: per-center amplitude multiplier
    """
    cx = nx * 0.5
    cy = ny * 0.5
    cz = nz * 0.5

    if spacing is None:
        spacing = wavelength

    if phases is None:
        phases = [0.0] * K

    if len(phases) != K:
        raise ValueError(f"phases length {len(phases)} does not match K={K}")

    if geometry == "single":
        if K != 1:
            raise ValueError("geometry 'single' requires K = 1")
        centers = [WC(cx, cy, cz, phases[0], True, amplitude)]

    elif geometry == "pair":
        if K != 2:
            raise ValueError("geometry 'pair' requires K = 2")
        d = spacing * 0.5
        centers = [
            WC(cx - d, cy, cz, phases[0], True, amplitude),
            WC(cx + d, cy, cz, phases[1], True, amplitude),
        ]

    elif geometry == "line":
        if K < 2:
            raise ValueError("geometry 'line' requires K >= 2")
        total = (K - 1) * spacing
        start = -total * 0.5
        centers = []
        for i in range(K):
            x = cx + start + i * spacing
            centers.append(WC(x, cy, cz, phases[i], True, amplitude))

    elif geometry == "golden":
        if K < 1:
            raise ValueError("geometry 'golden' requires K >= 1")
        if shell_radius is None:
            shell_radius = 0.35 * spacing
        phi = math.pi * (3.0 - math.sqrt(5.0))
        centers = []
        for i in range(K):
            y = 1.0 - (2.0 * i + 1.0) / K
            r_y = math.sqrt(max(0.0, 1.0 - y * y))
            theta = phi * i
            px = math.cos(theta) * r_y * shell_radius
            pz = math.sin(theta) * r_y * shell_radius
            py = y * shell_radius
            centers.append(WC(cx + px, cy + py, cz + pz, phases[i], True, amplitude))

    elif geometry == "tetrahedron_10_locked":
        if K != 10:
            raise ValueError("geometry 'tetrahedron_10_locked' requires K = 10")
        centers = _tetrahedron_10_locked(cx, cy, cz, spacing, phases, amplitude)

    else:
        raise ValueError(f"unknown geometry: '{geometry}'")

    return WCState(centers=centers)


def _tetrahedron_10_locked(
    cx: float,
    cy: float,
    cz: float,
    spacing: float,
    phases: list[float],
    amplitude: float,
) -> list[WC]:
    """
    1-3-6 arrangement with radii = n * spacing.
    Inner 3 at 1.0 * spacing from center, outer 6 at 2.0 * spacing.
    Matches the 'locked' geometry from M4's earlier work.
    """
    r1 = 1.0 * spacing
    r2 = 2.0 * math.sqrt(3.0 / 5.0) * spacing
    h = r2 * math.sqrt(2.0 / 3.0)

    centers: list[WC] = []
    # 1. Center
    centers.append(WC(cx, cy, cz, phases[0], True, amplitude))

    # 2. Inner 3 (in XY plane, 120 deg apart)
    angles_inner = [math.radians(90.0), math.radians(210.0), math.radians(330.0)]
    for idx, a in enumerate(angles_inner):
        centers.append(
            WC(
                cx + r1 * math.cos(a),
                cy + r1 * math.sin(a),
                cz,
                phases[1 + idx],
                True,
                amplitude,
            )
        )

    # 3. Outer 6 (two layers of 3, rotated 60 deg from inner)
    angles_outer = [math.radians(30.0), math.radians(150.0), math.radians(270.0)]
    for idx, a in enumerate(angles_outer):
        centers.append(
            WC(
                cx + r2 * math.cos(a),
                cy + r2 * math.sin(a),
                cz - h,
                phases[4 + idx],
                True,
                amplitude,
            )
        )
    for idx, a in enumerate(angles_outer):
        centers.append(
            WC(
                cx + r2 * math.cos(a),
                cy + r2 * math.sin(a),
                cz + h,
                phases[7 + idx],
                True,
                amplitude,
            )
        )

    return centers
