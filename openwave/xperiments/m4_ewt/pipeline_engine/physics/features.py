"""
Physics feature types. These are contracts between physics processors,
not part of the engine. The engine never imports this module.

models in the platform may use these features; per CROSS_MODEL_TESTING.md,
Section 2, a borrowed field is native and untwisted unless its docstring
says otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass

import taichi as ti

# =============================================================================
# Grid and statistics
# =============================================================================


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
class WaveStats:
    """Mutable holder for measured quantities. Written by trackers."""

    amp_max: float = 0.0
    mass: float = 0.0


# =============================================================================
# Vector fields
# =============================================================================


@dataclass
class PsiTripleBuffer:
    """
    Three time levels of a vector field, for leapfrog integration.

    Base class. Subclasses name a specific mode (base wave, longitudinal,
    transverse). Processors should require the concrete subclass, not
    this base.

    Representation: a single ti.Vector.field(3, f32) per time level.
    Components are (x, y, z) in the engine's spatial axes. Not an
    internal triplet, not a complex scalar, not a director field.
    """

    psi: ti.Vector.field  # psi(t)
    psi_prev: ti.Vector.field  # psi(t - dt)
    psi_new: ti.Vector.field  # scratch / psi(t + dt)


@dataclass
class PsiBaseField(PsiTripleBuffer):
    """
    Background oscillation of the medium.

    In EWT this is the always-on base wave. It is not a particle
    disturbance; it is the ground state of the medium. Wave parameters
    (amplitude, wavelength, frequency) come from the UnitSystem feature
    and the experiment configuration.
    """

    pass


@dataclass
class PsiLongField(PsiTripleBuffer):
    """
    Longitudinal mode. Carries mass and charge in the EWT picture.

    This is the mode that responds to the EMC density gradient and
    participates in the push-out mechanism.
    """

    pass


@dataclass
class PsiTransField(PsiTripleBuffer):
    """
    Transverse mode. Carries spin and magnetism in the EWT picture.

    Coupled to the longitudinal mode at wave centres through the
    fine-structure constant. The conversion coefficient is loaded from
    GeometricConstants, not derived from the field.
    """

    pass


# =============================================================================
# Scalar fields
# =============================================================================


@dataclass
class EMCDensityField:
    """
    EMC packing density rho(r), single buffer.

    High away from matter (statutory background), low inside a soliton
    (push-out deficit). The gradient drives the pressure force on wave
    centres and modulates the local wave speed.

    Representation: scalar ti.field(f32), one value per voxel. Not a
    vector, not a tensor.
    """

    rho: ti.field


@dataclass
class EMCFluxField:
    """
    EMC flux through a surface, single buffer.

    Optional. Used by boundary processors and energy-budget trackers to
    account for EMC leaving or entering the simulated domain.

    Representation: scalar ti.field(f32), one value per voxel.
    """

    flux: ti.field
