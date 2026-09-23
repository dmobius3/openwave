"""
Tests for the physics feature types.

No Taichi required for the type-level checks. The allocation test uses
Taichi and is skipped when Taichi cannot initialise.
"""

from __future__ import annotations

import sys

from openwave.xperiments.m4_ewt.pipeline_engine.context import FeatureBag
from openwave.xperiments.m4_ewt.pipeline_engine.physics.features import (
    EMCDensityField,
    EMCFluxField,
    PsiBaseField,
    PsiLongField,
    PsiTransField,
    PsiTripleBuffer,
    WaveGrid,
    WaveStats,
)


# ---------------------------------------------------------------------------
# Type hierarchy
# ---------------------------------------------------------------------------


def test_field_classes_are_triple_buffer_subclasses():
    for cls in (PsiBaseField, PsiLongField, PsiTransField):
        assert issubclass(cls, PsiTripleBuffer), cls.__name__


def test_field_classes_are_distinct():
    # FeatureBag keys by type, so the three field classes must not be
    # aliases of each other.
    classes = {PsiBaseField, PsiLongField, PsiTransField}
    assert len(classes) == 3


def test_grid_and_stats_are_independent():
    grid = WaveGrid(nx=4, ny=4, nz=4, dx=0.1, c=1.0)
    stats = WaveStats()
    assert grid.max_size == 4
    assert grid.center == (2, 2, 2)
    assert stats.amp_max == 0.0
    assert stats.mass == 0.0


# ---------------------------------------------------------------------------
# FeatureBag keying
# ---------------------------------------------------------------------------


def _make_dummy_triple():
    """Return three placeholder objects standing in for the Taichi fields."""

    class _Sentinel:
        pass

    return _Sentinel(), _Sentinel(), _Sentinel()


def test_feature_bag_keys_long_and_trans_separately():
    bag = FeatureBag()

    l_psi, l_prev, l_new = _make_dummy_triple()
    long_field = PsiLongField(psi=l_psi, psi_prev=l_prev, psi_new=l_new)

    t_psi, t_prev, t_new = _make_dummy_triple()
    trans_field = PsiTransField(psi=t_psi, psi_prev=t_prev, psi_new=t_new)

    bag.set(long_field)
    bag.set(trans_field)

    assert bag.require(PsiLongField) is long_field
    assert bag.require(PsiTransField) is trans_field


def test_feature_bag_base_key_resolves_to_most_recent():
    # Both PsiLongField and PsiTransField register under PsiTripleBuffer
    # through their MRO. The most recent set wins for the base key.
    bag = FeatureBag()

    l_psi, l_prev, l_new = _make_dummy_triple()
    long_field = PsiLongField(psi=l_psi, psi_prev=l_prev, psi_new=l_new)

    t_psi, t_prev, t_new = _make_dummy_triple()
    trans_field = PsiTransField(psi=t_psi, psi_prev=t_prev, psi_new=t_new)

    bag.set(long_field)
    bag.set(trans_field)

    assert bag.require(PsiTripleBuffer) is trans_field

    bag.set(long_field)
    assert bag.require(PsiTripleBuffer) is long_field


def test_feature_bag_scalar_fields_are_independent():
    bag = FeatureBag()
    rho_field = EMCDensityField(rho=object())
    flux_field = EMCFluxField(flux=object())
    bag.set(rho_field)
    bag.set(flux_field)
    assert bag.require(EMCDensityField) is rho_field
    assert bag.require(EMCFluxField) is flux_field


# ---------------------------------------------------------------------------
# Allocation (Taichi)
# ---------------------------------------------------------------------------


def test_allocation_provides_all_features():
    import taichi as ti

    ti.init(arch=ti.cpu, log_level=ti.ERROR)

    from openwave.xperiments.m4_ewt.pipeline_engine.pipeline import Pipeline
    from openwave.xperiments.m4_ewt.pipeline_engine.runner import Runner
    from openwave.xperiments.m4_ewt.pipeline_engine.sinks import InMemorySink
    from openwave.xperiments.m4_ewt.pipeline_engine.physics.allocator import (
        AllocateWaveField,
    )

    class AllocPipeline(Pipeline):
        def __init__(self) -> None:
            super().__init__()
            self.add(AllocateWaveField(nx=4, ny=4, nz=4, dx=0.1, c=1.0))

    sinks = {"session": InMemorySink()}
    runner = Runner(sinks)
    ctx = runner.run(AllocPipeline(), name="alloc_test", params={}, max_steps=0)

    assert not ctx.diag.errors, ctx.diag.errors
    for key in (
        WaveGrid,
        PsiBaseField,
        PsiLongField,
        PsiTransField,
        EMCDensityField,
        EMCFluxField,
        WaveStats,
    ):
        assert key in ctx.data, key.__name__


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def main() -> int:
    tests = [
        test_field_classes_are_triple_buffer_subclasses,
        test_field_classes_are_distinct,
        test_grid_and_stats_are_independent,
        test_feature_bag_keys_long_and_trans_separately,
        test_feature_bag_base_key_resolves_to_most_recent,
        test_feature_bag_scalar_fields_are_independent,
        test_allocation_provides_all_features,
    ]
    passed = 0
    for t in tests:
        try:
            t()
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            print(f"ERROR: {t.__name__}: {type(e).__name__}: {e}")
        else:
            print(f"PASS: {t.__name__}")
            passed += 1
    print(f"\n{passed}/{len(tests)} tests passed")
    return 0 if passed == len(tests) else 1


if __name__ == "__main__":
    sys.exit(main())