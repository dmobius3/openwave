"""
Wave center types. Pure Python, no Taichi.

WC      -- a single wave center: position, phase, active flag, amplitude.
WCState -- all wave centers for the current run. Passed in via initial_features.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class WC:
    """A single wave center in grid coordinates."""

    x: float
    y: float
    z: float
    phase: float = 0.0  # radians; 0 = electron-like, pi = positron-like
    active: bool = True
    amplitude: float = 1.0


@dataclass
class WCState:
    """All wave centers in the current run. Mutable; processors may update positions."""

    centers: list[WC]

    @property
    def K(self) -> int:
        return len(self.centers)

    @property
    def active_count(self) -> int:
        return sum(1 for wc in self.centers if wc.active)
