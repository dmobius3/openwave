"""
Three-plane sampling for scalar 3D fields.

Estimates the average of a scalar field over three centre-axis planes
instead of the whole volume. Samples ~3 N^2 voxels instead of N^3
(~3% at N = 100), a deliberate performance compromise documented in
the plan (item 1.2). Assumes the field is near-isotropic: valid for
solitons, waves, and slowly varying fields.

Engine layer: takes raw grid dimensions, not a physics WaveGrid.
Domain-agnostic and Taichi-only.

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects, not stringified annotations.
"""

import taichi as ti


class ThreePlaneSampler:
    """
    Cached 2D slice buffers plus kernels for three-plane sampling.

    One sampler serves one grid shape. Slice buffers are three 2D f32
    fields (~3 N^2 values). Allocate once in setup(), reuse in process().
    """

    def __init__(self, nx: int, ny: int, nz: int):
        if nx < 2 or ny < 2 or nz < 2:
            raise ValueError(f"grid dimensions must be >= 2, got ({nx}, {ny}, {nz})")
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self._slice_xy = ti.field(dtype=ti.f32, shape=(nx, ny))
        self._slice_xz = ti.field(dtype=ti.f32, shape=(nx, nz))
        self._slice_yz = ti.field(dtype=ti.f32, shape=(ny, nz))

    def average(self, field: ti.field) -> float:
        """
        Mean of the three centre-plane slices, edges excluded.

        Edges are skipped because the interior processors treat the
        outer shell as Dirichlet; including it would bias the average.

        Midpoint convention: mid = nx // 2 (integer division). On an
        even grid this places the plane just above the geometric
        centre; on an odd grid it is the exact centre voxel. The
        convention is fixed and tested
        (test_three_plane_sampler_linear_field_strict).
        """

        mid_x, mid_y, mid_z = self.nx // 2, self.ny // 2, self.nz // 2

        _copy_xy(field, self._slice_xy, mid_z)
        _copy_xz(field, self._slice_xz, mid_y)
        _copy_yz(field, self._slice_yz, mid_x)

        xy = self._slice_xy.to_numpy()[1:-1, 1:-1]
        xz = self._slice_xz.to_numpy()[1:-1, 1:-1]
        yz = self._slice_yz.to_numpy()[1:-1, 1:-1]

        n = xy.size + xz.size + yz.size
        return float((xy.sum() + xz.sum() + yz.sum()) / n)


@ti.kernel
def _copy_xy(src: ti.template(), dst: ti.template(), mid_z: ti.i32):
    for i, j in dst:
        dst[i, j] = src[i, j, mid_z]


@ti.kernel
def _copy_xz(src: ti.template(), dst: ti.template(), mid_y: ti.i32):
    for i, k in dst:
        dst[i, k] = src[i, mid_y, k]


@ti.kernel
def _copy_yz(src: ti.template(), dst: ti.template(), mid_x: ti.i32):
    for j, k in dst:
        dst[j, k] = src[mid_x, j, k]