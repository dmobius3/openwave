"""
Tests for UpdateEMCDensityProcessor and UpdateWaveSpeedProcessor.

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_emc

Note: no `from __future__ import annotations` here. Taichi kernel
definitions need live type objects.
"""

import sys
import traceback

import numpy as np
import taichi as ti

from .features import (
    EMCDensityField,
    PsiLongField,
    WaveGrid,
    WaveSpeedField,
)
from ..pipeline import BaseProcessor, Pipeline, Stage
from .allocator import AllocateWaveField, AllocateWaveSpeed
from .emc import UpdateEMCDensityProcessor, UpdateWaveSpeedProcessor
from .units import NaturalUnitSystem, UnitSystem

_TI_INITIALIZED = False


def _ti_init():
    global _TI_INITIALIZED
    if not _TI_INITIALIZED:
        ti.init(arch=ti.cpu, log_level=ti.ERROR)
        _TI_INITIALIZED = True


class _SeedPsi(BaseProcessor):
    """Fills PsiLong.psi with a chosen function. For tests only."""

    name = "_SeedPsi"
    stage = Stage.PRE_UPDATE
    order = 10
    requires = (WaveGrid, PsiLongField)
    provides = ()

    def __init__(self, mode: str = "zero", amp: float = 1.0, sigma: float = 3.0):
        self.mode = mode
        self.amp = amp
        self.sigma = sigma

    def process(self, ctx) -> None:
        if ctx.sim.step > 0:
            return
        grid = ctx.data.require(WaveGrid)
        field = ctx.data.require(PsiLongField)
        _seed_psi(
            field.psi,
            field.psi_prev,
            grid.nx,
            grid.ny,
            grid.nz,
            self.amp,
            self.sigma,
            0 if self.mode == "zero" else 1 if self.mode == "const" else 2,
        )


@ti.kernel
def _seed_psi(
    psi: ti.template(),
    prev: ti.template(),
    nx: ti.i32,
    ny: ti.i32,
    nz: ti.i32,
    amp: ti.f32,
    sigma: ti.f32,
    mode: ti.i32,
):
    cx = nx * 0.5
    cy = ny * 0.5
    cz = nz * 0.5
    for i, j, k in ti.ndrange(nx, ny, nz):
        v = ti.Vector([0.0, 0.0, 0.0])
        if mode == 1:
            v = ti.Vector([amp, 0.0, 0.0])
        elif mode == 2:
            dx = ti.cast(i, ti.f32) - cx
            dy = ti.cast(j, ti.f32) - cy
            dz = ti.cast(k, ti.f32) - cz
            r2 = dx * dx + dy * dy + dz * dz
            v = ti.Vector([amp * ti.exp(-r2 / (2.0 * sigma * sigma)), 0.0, 0.0])
        psi[i, j, k] = v
        prev[i, j, k] = v


def _build_emc_pipeline(
    nx=16,
    ny=16,
    nz=16,
    seed_mode="zero",
    amp=1.0,
    sigma=3.0,
    beta_rho=1.0,
    with_wave_speed=True,
    with_emc=True,
):
    class P(Pipeline):
        def __init__(self):
            super().__init__(external_provides=(UnitSystem,))
            self.add(AllocateWaveField(nx=nx, ny=ny, nz=nz, dx=1.0))
            if with_wave_speed:
                self.add(AllocateWaveSpeed())
            self.add(_SeedPsi(mode=seed_mode, amp=amp, sigma=sigma))
            if with_emc:
                self.add(UpdateEMCDensityProcessor(beta_rho=beta_rho))
            if with_wave_speed:
                self.add(UpdateWaveSpeedProcessor())

    return P()


def _run(pipeline, max_steps=1):
    from ..runner import Runner
    from ..sinks import InMemorySink

    units = NaturalUnitSystem()
    return Runner({"session": InMemorySink()}).run(
        pipeline,
        name="emc_test",
        params={},
        dt=0.1,
        max_steps=max_steps,
        initial_features=[units],
    )


# ----------------------------------------------------------------------
# UpdateEMCDensityProcessor
# ----------------------------------------------------------------------


def test_emc_zero_psi_gives_statutory_rho():
    """
    psi = 0 -> rho = 1.0 everywhere.

    Mutation: formula changed to rho = beta*|psi|^2 -> rho = 0, check fails.
    """
    _ti_init()
    ctx = _run(_build_emc_pipeline(seed_mode="zero", beta_rho=1.0), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    rho = ctx.data.require(EMCDensityField).rho.to_numpy()
    # Everywhere in the grid, rho = 1.0 (interior AND boundary; the
    # processor writes the full grid).
    assert np.allclose(rho, 1.0, atol=1e-6), (rho.min(), rho.max())


def test_emc_const_psi_gives_uniform_deficit():
    """
    psi = [amp, 0, 0] everywhere -> rho = 1 - beta*amp^2 everywhere.

    Mutation: beta ignored -> rho = 1, check fails.
    """
    _ti_init()
    amp = 0.5
    beta = 0.8
    ctx = _run(_build_emc_pipeline(seed_mode="const", amp=amp, beta_rho=beta), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    rho = ctx.data.require(EMCDensityField).rho.to_numpy()
    expected = 1.0 - beta * amp * amp
    assert np.allclose(rho, expected, atol=1e-6), (rho.min(), rho.max(), expected)


def test_emc_gaussian_psi_gives_gaussian_deficit():
    """
    psi = amp * Gaussian -> rho = 1 - beta*amp^2 * Gaussian^2.
    Center voxel has the deepest deficit.

    Mutation: |psi|^2 replaced by |psi| -> deficit curve different.
    """
    _ti_init()
    amp = 0.6
    beta = 1.0
    sigma = 3.0
    N = 16
    ctx = _run(
        _build_emc_pipeline(
            nx=N, ny=N, nz=N, seed_mode="gauss", amp=amp, sigma=sigma, beta_rho=beta
        ),
        max_steps=1,
    )
    assert ctx.diag.errors == [], ctx.diag.errors

    rho = ctx.data.require(EMCDensityField).rho.to_numpy()
    c = N // 2
    # Center: r=0, |psi|^2 = amp^2
    expected_center = 1.0 - beta * amp * amp
    assert abs(rho[c, c, c] - expected_center) < 1e-4, (rho[c, c, c], expected_center)
    # Far corner: r large, |psi|^2 ~ 0, rho ~ 1
    assert abs(rho[0, 0, 0] - 1.0) < 1e-3, rho[0, 0, 0]
    # Monotonic: center <= edge (deficit deepest at center)
    assert rho[c, c, c] < rho[0, 0, 0]


def test_emc_beta_zero_is_noop():
    """
    beta_rho = 0 -> rho = 1 regardless of psi.

    Mutation: beta_rho ignored -> rho != 1 for non-zero psi, check fails.
    """
    _ti_init()
    ctx = _run(_build_emc_pipeline(seed_mode="const", amp=2.0, beta_rho=0.0), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors
    rho = ctx.data.require(EMCDensityField).rho.to_numpy()
    assert np.allclose(rho, 1.0, atol=1e-6)


# ----------------------------------------------------------------------
# UpdateWaveSpeedProcessor
# ----------------------------------------------------------------------


def test_wave_speed_statutory_rho_gives_c0_squared():
    """
    Smoke check: rho = 1 (no deficit) -> c^2 = c0^2 = units.c**2.

    This test cannot tell c0^2 * rho from c0^2 / rho, because at rho = 1
    both give c0^2. The discriminating test is
    test_wave_speed_deficit_halves_c_squared, where rho = 0.5. This test
    catches only gross errors: missing c0 factor, missing rho factor,
    wrong sign of the exponent.
    """
    _ti_init()
    ctx = _run(_build_emc_pipeline(seed_mode="zero", beta_rho=1.0), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    c2 = ctx.data.require(WaveSpeedField).c2_local.to_numpy()
    units = NaturalUnitSystem()
    expected = float(units.c) ** 2
    assert np.allclose(c2, expected, atol=1e-6), (c2.min(), c2.max(), expected)


def test_wave_speed_deficit_halves_c_squared():
    """
    rho = 0.5 -> c^2 = 0.5 * c0^2.

    Mutation: linear factor dropped -> c^2 = c0^2, check fails.
    """
    _ti_init()
    # psi = 0.7071 = sqrt(0.5) with beta=1 gives rho = 1 - 0.5 = 0.5
    amp = np.sqrt(0.5)
    ctx = _run(_build_emc_pipeline(seed_mode="const", amp=float(amp), beta_rho=1.0), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    c2 = ctx.data.require(WaveSpeedField).c2_local.to_numpy()
    units = NaturalUnitSystem()
    expected = 0.5 * float(units.c) ** 2
    assert np.allclose(c2, expected, atol=1e-5), (c2.min(), c2.max(), expected)


def test_wave_speed_full_deficit_gives_zero():
    """
    psi = 1 with beta = 1 -> rho = 0 -> c^2 = 0.

    Mutation: c^2 floored at some minimum -> c^2 != 0, check fails.
    """
    _ti_init()
    ctx = _run(_build_emc_pipeline(seed_mode="const", amp=1.0, beta_rho=1.0), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors
    c2 = ctx.data.require(WaveSpeedField).c2_local.to_numpy()
    assert np.allclose(c2, 0.0, atol=1e-6), (c2.min(), c2.max())


def test_wave_speed_spatial_variation_matches_rho():
    """
    Gaussian psi -> Gaussian rho -> Gaussian c^2. Monotonic ordering
    must be preserved (deeper deficit = lower c^2).

    Mutation: c^2 formula uses wrong field (e.g. psi instead of rho)
    -> ordering breaks.
    """
    _ti_init()
    N = 16
    ctx = _run(
        _build_emc_pipeline(nx=N, ny=N, nz=N, seed_mode="gauss", amp=0.8, sigma=3.0, beta_rho=1.0),
        max_steps=1,
    )
    assert ctx.diag.errors == [], ctx.diag.errors

    rho = ctx.data.require(EMCDensityField).rho.to_numpy()
    c2 = ctx.data.require(WaveSpeedField).c2_local.to_numpy()

    c = N // 2
    # Center has lower rho -> lower c2 than corner.
    assert rho[c, c, c] < rho[0, 0, 0]
    assert c2[c, c, c] < c2[0, 0, 0]
    # c2 proportional to rho, same ratio everywhere.
    ratio = c2[c, c, c] / rho[c, c, c]
    units = NaturalUnitSystem()
    assert abs(ratio - float(units.c) ** 2) < 1e-4, (ratio, float(units.c) ** 2)


def test_negative_c2_regime_is_reachable():
    """
    The plan documents beta_rho |psi|^2 > 1 as a regime outside the
    model's domain of validity: rho < 0, hence c2_local < 0. The plan
    does not floor rho; the negative value is the signal.

    This test pins that the regime is reachable and that the sign of
    the output matches the sign of the input. It does not guard the
    regime, and it does not test that the simulation stays stable
    there.

    Mutation caught: a floor on rho, e.g. max(0.0, 1 - beta*|psi|^2),
    hides the regime and the test fails on the sign check.
    """
    _ti_init()
    # psi = 1.5, beta_rho = 1.0 -> rho = 1 - 2.25 = -1.25
    ctx = _run(_build_emc_pipeline(seed_mode="const", amp=1.5, beta_rho=1.0), max_steps=1)
    assert ctx.diag.errors == [], ctx.diag.errors

    rho = ctx.data.require(EMCDensityField).rho.to_numpy()
    c2 = ctx.data.require(WaveSpeedField).c2_local.to_numpy()

    assert rho.min() < 0.0, f"rho must reach negative values in this regime, got min {rho.min()}"
    assert c2.min() < 0.0, f"c2 must track rho in sign, got min {c2.min()}"
    # Center voxel: rho = 1 - 1.5^2 = -1.25, c2 = -1.25 * c0^2.
    c = rho.shape[0] // 2
    assert abs(rho[c, c, c] - (-1.25)) < 1e-5, rho[c, c, c]


# ----------------------------------------------------------------------
# Runner
# ----------------------------------------------------------------------


def main() -> int:
    tests = [
        test_emc_zero_psi_gives_statutory_rho,
        test_emc_const_psi_gives_uniform_deficit,
        test_emc_gaussian_psi_gives_gaussian_deficit,
        test_emc_beta_zero_is_noop,
        test_wave_speed_statutory_rho_gives_c0_squared,
        test_wave_speed_deficit_halves_c_squared,
        test_wave_speed_full_deficit_gives_zero,
        test_wave_speed_spatial_variation_matches_rho,
        test_negative_c2_regime_is_reachable,
    ]
    passed = 0
    for t in tests:
        try:
            t()
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            traceback.print_exc()
        except Exception as e:
            print(f"ERROR: {t.__name__}: {type(e).__name__}: {e}")
            traceback.print_exc()
        else:
            print(f"PASS: {t.__name__}")
            passed += 1
    print(f"\n{passed}/{len(tests)} tests passed")
    return 0 if passed == len(tests) else 1


if __name__ == "__main__":
    sys.exit(main())
