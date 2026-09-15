"""
Runnable demo: linear or cubic-nonlinear wave on a 3D grid, with optional
live Taichi visualization. Supports multi-center seeding via WCState.

Run as a module from the OpenWave root:

    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._demo
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._demo --k 2 --geometry pair
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._demo --k 10 --geometry tetrahedron_10_locked

Close the window (or press ESC) to stop the run.
"""

import argparse
from pathlib import Path
from typing import Any, Mapping

import taichi as ti

from ..pipeline import Pipeline
from ..runner import Runner
from ..sinks import JsonSessionSink
from ..loggers import LogProcessor

from .allocator import AllocateWaveField
from .seed import SeedMultiCenter
from .evolution import LaplacianProcessor, LeapfrogProcessor
from .nonlinearity import NonlinearCubic
from .boundary import DirichletBoundaryProcessor
from .measure import AmplitudeTracker
from .visualize import TaichiWindowProcessor
from .features import WaveStats
from .wc_types import WCState
from .wc_factory import build_wc_state


def _payload(ctx) -> Mapping[str, Any]:
    stats = ctx.data.require(WaveStats)
    wc_state = ctx.data.require(WCState)
    return {
        "amp_max": stats.amp_max,
        "mass": stats.mass,
        "K": wc_state.K,
        "active": wc_state.active_count,
    }


class WavePipeline(Pipeline):
    """
    Same pipeline for all K and geometries. The only difference between runs
    is the WCState passed in via initial_features.
    """

    def __init__(
        self,
        *,
        grid: int,
        gamma: float,
        with_window: bool = True,
        external_provides: tuple[type, ...] = (),
    ):
        super().__init__(external_provides=external_provides)
        self.add(AllocateWaveField(nx=grid, ny=grid, nz=grid, dx=1.0, c=1.0))
        self.add(SeedMultiCenter(radius=4.0))
        self.add(LaplacianProcessor())
        if gamma > 0.0:
            self.add(NonlinearCubic(gamma=gamma))
        self.add(LeapfrogProcessor())
        self.add(DirichletBoundaryProcessor())
        self.add(AmplitudeTracker(every=10))
        if with_window:
            self.add(
                TaichiWindowProcessor(
                    size=(512, 512),
                    scale=8,
                    amp_scale=3.0,
                    title=f"M4 wave (K, gamma={gamma})",
                )
            )
        self.add(
            LogProcessor(
                "SessionLog",
                "session",
                _payload,
                every=20,
                order=500,
            )
        )


def main() -> None:
    p = argparse.ArgumentParser(description="M4 wave demo")
    p.add_argument("--k", type=int, default=1, help="number of wave centers")
    p.add_argument(
        "--geometry",
        type=str,
        default="single",
        choices=["single", "pair", "line", "golden", "tetrahedron_10_locked"],
        help="geometry used to place the wave centers",
    )
    p.add_argument(
        "--spacing",
        type=float,
        default=None,
        help="characteristic spacing in grid units (default: wavelength)",
    )
    p.add_argument(
        "--wavelength", type=float, default=12.0, help="reference wavelength in grid units"
    )
    p.add_argument("--gamma", type=float, default=0.01, help="cubic coupling (0.0 = linear)")
    p.add_argument("--grid", type=int, default=64)
    p.add_argument("--steps", type=int, default=5000, help="max steps (headless mode)")
    p.add_argument("--no-window", action="store_true", help="run headless; exits after --steps")
    p.add_argument("--check-stateless", action="store_true", help="enable the stateless guard (dev/CI use)")
    args = p.parse_args()

    ti.init(arch=ti.cpu, log_level=ti.WARN)

    if not args.no_window:
        print("[info] window mode: --steps ignored; close the window to stop")

    # Build WCState from K and geometry. Same pipeline, different input.
    wc_state = build_wc_state(
        K=args.k,
        geometry=args.geometry,
        nx=args.grid,
        ny=args.grid,
        nz=args.grid,
        wavelength=args.wavelength,
        spacing=args.spacing,
    )

    out = Path(__file__).resolve().parent / "out_demo"
    sinks = {"session": JsonSessionSink(out / "session.json", flush_every=100)}
    runner = Runner(sinks, check_stateless=args.check_stateless)

    pipeline = WavePipeline(
        grid=args.grid,
        gamma=args.gamma,
        with_window=not args.no_window,
        external_provides=(WCState,),
    )

    max_steps = args.steps if args.no_window else 100_000

    ctx = runner.run(
        pipeline,
        name=f"wave_K{args.k}_{args.geometry}_gamma_{args.gamma}",
        params={
            "grid": args.grid,
            "gamma": args.gamma,
            "c": 1.0,
            "dt": 0.3,
            "K": args.k,
            "geometry": args.geometry,
        },
        output_dir=out,
        dt=0.3,
        max_steps=max_steps,
        initial_features=[wc_state],
    )

    print("=" * 64)
    print(f"Wave demo finished (K={args.k}, geometry={args.geometry}, gamma={args.gamma})")
    print("=" * 64)
    print(f"[info] stateless guard: {'ON' if args.check_stateless else 'OFF'}")
    print(f"steps simulated : {ctx.sim.step}")
    stats = ctx.data.require(WaveStats)
    print(f"final amp_max   : {stats.amp_max:.6f}")
    print(f"final mass      : {stats.mass:.6f}")
    print(f"errors          : {ctx.diag.errors}")
    print(f"session log     : {sinks['session'].path}")


if __name__ == "__main__":
    main()
