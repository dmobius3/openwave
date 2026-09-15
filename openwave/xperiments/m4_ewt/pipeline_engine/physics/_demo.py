"""
Runnable demo: linear wave on a 3D grid, with live Taichi visualization.

Run as a module from the OpenWave root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._demo

Close the window (or press ESC) to stop the run.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import taichi as ti

from ..pipeline import Pipeline
from ..runner import Runner
from ..sinks import JsonSessionSink
from ..loggers import LogProcessor

from .allocator import AllocateWaveField
from .seed import SeedPulse
from .evolution import LaplacianProcessor, LeapfrogProcessor
from .boundary import DirichletBoundaryProcessor
from .measure import AmplitudeTracker
from .visualize import TaichiWindowProcessor
from .features import WaveStats


def _payload(ctx) -> Mapping[str, Any]:
    stats = ctx.data.require(WaveStats)
    return {"amp_max": stats.amp_max, "mass": stats.mass}


class WaveDemoPipeline(Pipeline):
    def __init__(self) -> None:
        super().__init__()
        self.add(AllocateWaveField(nx=64, ny=64, nz=64, dx=1.0, c=1.0))
        self.add(SeedPulse(amplitude=1.0, radius=4.0))
        self.add(LaplacianProcessor())
        self.add(LeapfrogProcessor())
        self.add(DirichletBoundaryProcessor())
        self.add(AmplitudeTracker(every=10))
        self.add(TaichiWindowProcessor(
            size=(512, 512),
            scale=8,
            amp_scale=3.0,
            title="M4 wave demo",
        ))
        self.add(LogProcessor(
            "SessionLog", "session", _payload,
            every=20, order=500,
        ))


def main() -> None:
    ti.init(arch=ti.cpu, log_level=ti.WARN)

    out = Path("out_demo")
    sinks = {"session": JsonSessionSink(out / "session.json", flush_every=100)}
    runner = Runner(sinks, check_stateless=True)

    ctx = runner.run(
        WaveDemoPipeline(),
        name="wave_demo",
        params={"grid": 64, "c": 1.0, "dt": 0.3},
        output_dir=out,
        dt=0.3,
        max_steps=100_000,
    )

    print("=" * 64)
    print("Wave demo finished")
    print("=" * 64)
    print(f"steps simulated : {ctx.sim.step}")
    stats = ctx.data.require(WaveStats)
    print(f"final amp_max   : {stats.amp_max:.6f}")
    print(f"final mass      : {stats.mass:.6f}")
    print(f"errors          : {ctx.diag.errors}")
    print(f"session log     : {sinks['session'].path}")


if __name__ == "__main__":
    main()