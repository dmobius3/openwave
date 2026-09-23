"""Wave field allocator. Lifecycle-only processor."""

from __future__ import annotations

from ..pipeline import BaseProcessor

import taichi as ti

from .features import (
    EMCDensityField,
    EMCFluxField,
    PsiBaseField,
    PsiLongField,
    PsiTransField,
    WaveGrid,
    WaveStats,
)


def _triple_buffer(shape: tuple[int, int, int]):
    """Allocate one triple-buffered vector field."""
    return (
        ti.Vector.field(3, dtype=ti.f32, shape=shape),
        ti.Vector.field(3, dtype=ti.f32, shape=shape),
        ti.Vector.field(3, dtype=ti.f32, shape=shape),
    )


class AllocateWaveField(BaseProcessor):
    """
    Allocates the grid feature, the five field features, and the stats
    holder. Runs once at setup.

    All fields are allocated together so that a pipeline can use any
    subset without reallocation. Memory cost: three triple-buffered
    vector fields and two scalar fields, 29 f32 values per voxel against
    9 for a single triple buffer (about 3.2x; 232 MiB against 72 MiB on
    a 128^3 grid).
    """

    name = "AllocateWaveField"
    stage = None
    order = 10
    provides = (
        WaveGrid,
        PsiBaseField,
        PsiLongField,
        PsiTransField,
        EMCDensityField,
        EMCFluxField,
        WaveStats,
    )

    def __init__(self, nx: int, ny: int, nz: int, dx: float, c: float):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.c = c

    def setup(self, ctx) -> None:
        grid = WaveGrid(nx=self.nx, ny=self.ny, nz=self.nz, dx=self.dx, c=self.c)
        shape = (self.nx, self.ny, self.nz)
        ctx.data.set(grid)

        for cls in (PsiBaseField, PsiLongField, PsiTransField):
            psi, psi_prev, psi_new = _triple_buffer(shape)
            ctx.data.set(cls(psi=psi, psi_prev=psi_prev, psi_new=psi_new))

        ctx.data.set(EMCDensityField(rho=ti.field(dtype=ti.f32, shape=shape)))
        ctx.data.set(EMCFluxField(flux=ti.field(dtype=ti.f32, shape=shape)))
        ctx.data.set(WaveStats())