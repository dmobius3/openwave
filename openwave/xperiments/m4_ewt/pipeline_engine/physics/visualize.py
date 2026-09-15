"""
Interactive visualization via Taichi's ui.Window.

Renders a 2D slice of the field as a scaled image. Sets ctx.sim.should_stop
when the window is closed, so the runner exits cleanly.
"""

from pathlib import Path

from ..pipeline import BaseProcessor, Stage

import taichi as ti

from .features import PsiField, WaveGrid


class TaichiWindowProcessor(BaseProcessor):
    """
    Opens a Taichi UI window and renders an x-y slice at z = nz//2.
    Runs at the end of MEASURE. Reads only; does not mutate the field.

    Resources (window, canvas, image buffer) are allocated in setup()
    and released when the process exits. process() reads them only.
    """

    name = "TaichiWindow"
    stage = Stage.MEASURE
    order = 2000
    requires = (WaveGrid, PsiField)

    def __init__(
        self,
        size=(512, 512),
        scale=8,
        amp_scale=3.0,
        title="M4 wave demo",
        screenshot_every=0,
        screenshot_dir="out_demo",
    ):
        self.size = size
        self.scale = scale
        self.amp_scale = amp_scale
        self.title = title
        self.screenshot_every = screenshot_every
        self.screenshot_dir = Path(screenshot_dir)

    def setup(self, ctx) -> None:
        w, h = self.size
        self._img = ti.Vector.field(3, dtype=ti.f32, shape=(w, h))
        self._window = ti.ui.Window(self.title, (w, h))
        self._canvas = self._window.get_canvas()
        if self.screenshot_every > 0:
            self.screenshot_dir.mkdir(parents=True, exist_ok=True)

    def process(self, ctx) -> None:
        if not self._window.running:
            ctx.sim.should_stop = True
            return
        if self._window.is_pressed(ti.ui.ESCAPE):
            ctx.sim.should_stop = True
            return

        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiField)
        w, h = self.size
        z = grid.nz // 2

        _render_slice(
            field.psi,
            self._img,
            grid.nx,
            grid.ny,
            z,
            w,
            h,
            self.scale,
            self.amp_scale,
        )
        self._canvas.set_image(self._img)
        self._window.show()

        if self.screenshot_every > 0 and ctx.sim.step % self.screenshot_every == 0:
            path = self.screenshot_dir / f"frame_{ctx.sim.step:06d}.png"
            self._window.save_image(str(path))

    def teardown(self, ctx) -> None:
        # Taichi handles window cleanup on process exit; nothing to do.
        pass


@ti.kernel
def _render_slice(
    psi: ti.template(),
    img: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    z: ti.i32,
    w: ti.i32,
    h: ti.i32,
    scale: ti.i32,
    amp: ti.f32,
):
    for px, py in ti.ndrange(w, h):
        i = px // scale
        j = py // scale
        if i < nx and j < ny:
            v = psi[i, j, z].norm() * amp
            v = ti.min(v, 1.0)
            # Simple blue -> red heatmap
            img[px, py] = ti.Vector([v, v * 0.4, 1.0 - v])
        else:
            img[px, py] = ti.Vector([0.0, 0.0, 0.0])
