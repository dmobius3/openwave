"""
Physics feature types. These are contracts between physics processors,
not part of the engine. The engine never imports this module.
"""

from __future__ import annotations

from dataclasses import dataclass

import taichi as ti


@dataclass
class WaveGrid:
    """Grid + wave parameters for a scalar/vector wave field."""

    nx: int
    ny: int
    nz: int
    dx: float
    c: float

    @property
    def max_size(self) -> int:
        return max(self.nx, self.ny, self.nz)

    @property
    def center(self) -> tuple[int, int, int]:
        return (self.nx // 2, self.ny // 2, self.nz // 2)


@dataclass
class PsiField:
    """Three time levels of the wave field. Leapfrog needs all three."""

    psi: ti.Vector.field  # psi(t)
    psi_prev: ti.Vector.field  # psi(t-dt)
    psi_new: ti.Vector.field  # scratch / psi(t+dt)


@dataclass
class WaveStats:
    """Mutable holder for measured quantities. Written by trackers."""

    amp_max: float = 0.0
    mass: float = 0.0
