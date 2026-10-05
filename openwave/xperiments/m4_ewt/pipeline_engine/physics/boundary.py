"""Outer boundary conditions (plan item 1.7): dirichlet, periodic, reflecting."""

from ..pipeline import BaseProcessor, Stage

import taichi as ti

from .features import BoundaryCondition, PsiLongField, WaveGrid



class DirichletBoundaryProcessor(BaseProcessor):
    """Zero the outer shell of the grid every step."""

    name = "DirichletBoundary"
    stage = Stage.POST_UPDATE
    order = 10
    requires = (WaveGrid, PsiLongField)

    def process(self, ctx) -> None:
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiLongField)
        _dirichlet(field.psi, field.psi_prev, grid.nx, grid.ny, grid.nz)


@ti.kernel
def _dirichlet(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    z = ti.Vector([0.0, 0.0, 0.0])
    for i, j, k in ti.ndrange(nx, ny, nz):
        on_edge = i == 0 or i == nx - 1 or j == 0 or j == ny - 1 or k == 0 or k == nz - 1
        if on_edge:
            psi[i, j, k] = z
            prev[i, j, k] = z

@ti.kernel
def _boundary_dirichlet(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    z = ti.Vector([0.0, 0.0, 0.0])
    for i, j, k in ti.ndrange(nx, ny, nz):
        on_edge = (
            i == 0 or i == nx - 1
            or j == 0 or j == ny - 1
            or k == 0 or k == nz - 1
        )
        if on_edge:
            psi[i, j, k] = z
            prev[i, j, k] = z


@ti.kernel
def _boundary_periodic(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    for j, k in ti.ndrange(ny, nz):
        psi[0, j, k] = psi[nx - 2, j, k]
        psi[nx - 1, j, k] = psi[1, j, k]
        prev[0, j, k] = prev[nx - 2, j, k]
        prev[nx - 1, j, k] = prev[1, j, k]
    for i, k in ti.ndrange(nx, nz):
        psi[i, 0, k] = psi[i, ny - 2, k]
        psi[i, ny - 1, k] = psi[i, 1, k]
        prev[i, 0, k] = prev[i, ny - 2, k]
        prev[i, ny - 1, k] = prev[i, 1, k]
    for i, j in ti.ndrange(nx, ny):
        psi[i, j, 0] = psi[i, j, nz - 2]
        psi[i, j, nz - 1] = psi[i, j, 1]
        prev[i, j, 0] = prev[i, j, nz - 2]
        prev[i, j, nz - 1] = prev[i, j, 1]


@ti.kernel
def _boundary_reflecting(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
):
    for j, k in ti.ndrange(ny, nz):
        psi[0, j, k] = psi[1, j, k]
        psi[nx - 1, j, k] = psi[nx - 2, j, k]
        prev[0, j, k] = prev[1, j, k]
        prev[nx - 1, j, k] = prev[nx - 2, j, k]
    for i, k in ti.ndrange(nx, nz):
        psi[i, 0, k] = psi[i, 1, k]
        psi[i, ny - 1, k] = psi[i, ny - 2, k]
        prev[i, 0, k] = prev[i, 1, k]
        prev[i, ny - 1, k] = prev[i, ny - 2, k]
    for i, j in ti.ndrange(nx, ny):
        psi[i, j, 0] = psi[i, j, 1]
        psi[i, j, nz - 1] = psi[i, j, nz - 2]
        prev[i, j, 0] = prev[i, j, 1]
        prev[i, j, nz - 1] = prev[i, j, nz - 2]


class BoundaryProcessor(BaseProcessor):
    """
    Apply the registered BoundaryCondition to one field's triple buffer.

    field_type selects which field the boundary applies to. Requires
    WaveGrid, BoundaryCondition, and the concrete field class.

    kind_override is for tests only. When set, it overrides the
    feature's kind. Production code should leave it None and let the
    feature decide.

    Flux accounting through the boundary is not implemented here. It
    belongs to plan Section 6, item 1.8 (EnergyBudget, Phase C). The
    tests in _test_boundary.py check the field's behaviour at the
    edge, not an energy budget.

    The three kernels are first-order: dirichlet zeroes the outer
    shell, a fixed-end wall that reflects with the sign inverted (not
    a Mur or PML absorber); periodic copies the opposite face;
    reflecting copies the interior neighbour.
    """

    name = "BoundaryProcessor"
    stage = Stage.POST_UPDATE
    order = 10
    provides = ()

    def __init__(self, field_type=PsiLongField, kind_override=None):
        self.field_type = field_type
        self.kind_override = kind_override
        self.requires = (WaveGrid, BoundaryCondition, field_type)

    def process(self, ctx):
        grid = ctx.data.require(WaveGrid)
        bc = ctx.data.require(BoundaryCondition)
        field = ctx.data.require(self.field_type)
        kind = self.kind_override or bc.kind
        if kind == "dirichlet":
            _boundary_dirichlet(field.psi, field.psi_prev,
                                grid.nx, grid.ny, grid.nz)
        elif kind == "periodic":
            _boundary_periodic(field.psi, field.psi_prev,
                               grid.nx, grid.ny, grid.nz)
        elif kind == "reflecting":
            _boundary_reflecting(field.psi, field.psi_prev,
                                 grid.nx, grid.ny, grid.nz)
        else:
            raise ValueError(f"BoundaryProcessor: unknown kind {kind!r}")
