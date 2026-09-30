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
    TrackerFields,
    WaveGrid,
    WaveSpeedField,
    WaveStats,
)

class AllocateTrackers(BaseProcessor):
    """
    Allocates the TrackerFields feature. Must run AFTER AllocateWaveField,
    since the per-voxel fields are shaped by WaveGrid and rho_local
    aliases EMCDensityField.rho.

    All per-voxel trackers start at zero. The processor that populates
    them lives in the physics layer.

    rho_local shares the buffer with EMCDensityField.rho: same ti.field
    object, not a copy. Any write through one is visible through the
    other. This avoids duplicating ~40 MB on a 128^3 grid and makes the
    "tracker mirrors the source" contract explicit.
    """

    name = "AllocateTrackers"
    stage = None
    order = 20
    provides = (TrackerFields,)
    requires = (WaveGrid, EMCDensityField)

    def setup(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        emc = ctx.data.require(EMCDensityField)
        shape = (grid.nx, grid.ny, grid.nz)

        ctx.data.set(
            TrackerFields(
                amp_local=ti.field(dtype=ti.f32, shape=shape),
                freq_local=ti.field(dtype=ti.f32, shape=shape),
                energy_long_local=ti.field(dtype=ti.f32, shape=shape),
                energy_trans_local=ti.field(dtype=ti.f32, shape=shape),
                rho_local=emc.rho,  # shared buffer, not a copy
                last_crossing=ti.field(dtype=ti.f32, shape=shape),
                amp_global=ti.field(dtype=ti.f32, shape=()),
                freq_global=ti.field(dtype=ti.f32, shape=()),
                energy_global=ti.field(dtype=ti.f32, shape=()),
            )
        )

class AllocateWaveSpeed(BaseProcessor):
    """
    Allocates WaveSpeedField.c2_local. Must run AFTER AllocateWaveField,
    since the per-voxel field is shaped by WaveGrid.

    Separate from AllocateWaveField so pipelines that use the constant
    Laplacian (LaplacianProcessor) do not pay for the extra buffer.
    Pipelines that use LaplacianVariableCoeffProcessor must include
    this allocator.
    """

    name = "AllocateWaveSpeed"
    stage = None
    order = 25
    provides = (WaveSpeedField,)
    requires = (WaveGrid,)

    def setup(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        ctx.data.set(
            WaveSpeedField(
                c2_local=ti.field(dtype=ti.f32, shape=(grid.nx, grid.ny, grid.nz))
            )
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


