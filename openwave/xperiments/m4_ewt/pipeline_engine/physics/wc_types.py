"""
Wave center types. Pure Python, no Taichi.

WC      -- a single wave center: position, phase, active flag, amplitude.
WCState -- all wave centers for the current run. Passed in via initial_features.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class WC:
    """A single wave center in grid coordinates.

    phase: radians. Used only as a sign convention for the seed envelope
           (cos(phase)). This is an imposed sign, NOT a charge model. The
           emergent-charge question is M4.2.

    Reflector attributes (plan Section 6, item 1.5):
      reflect_coeff_long  -- amplitude reflection coefficient for the
                             longitudinal mode. Default 1.0.
      reflect_coeff_trans -- amplitude conversion coefficient into the
                             transverse mode. Default 0.0.
      phase_shift         -- reflection phase shift in radians. Default 0.0.
    """

    x: float
    y: float
    z: float
    phase: float = 0.0
    active: bool = True
    amplitude: float = 1.0
    reflect_coeff_long: float = 1.0
    reflect_coeff_trans: float = 0.0
    phase_shift: float = 0.0


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
