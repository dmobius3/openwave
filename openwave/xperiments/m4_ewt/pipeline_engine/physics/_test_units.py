"""
Tests for the UnitSystem contract and its implementations.

No Taichi required. Runs in milliseconds.
"""

from __future__ import annotations

import math
import sys

from openwave.common import constants
from openwave.xperiments.m4_ewt.pipeline_engine.context import FeatureBag
from openwave.xperiments.m4_ewt.pipeline_engine.physics.units import (
    NaturalUnitSystem,
    OpenWaveUnitSystem,
    SIUnitSystem,
    UnitSystem,
    make_unit_system,
)


# ---------------------------------------------------------------------------
# Construction and basic invariants
# ---------------------------------------------------------------------------


def test_natural_default_construction():
    units = NaturalUnitSystem()
    assert units.grid_voxels_per_lambda == 20
    assert units.cfl_safety == 0.9


def test_openwave_default_construction():
    units = OpenWaveUnitSystem()
    assert units.grid_voxels_per_lambda == 12
    assert units.cfl_safety == 0.9


def test_si_default_construction():
    units = SIUnitSystem()
    assert units.grid_voxels_per_lambda == 20
    assert units.cfl_safety == 0.9


# ---------------------------------------------------------------------------
# Inheritance
# ---------------------------------------------------------------------------


def test_all_implementations_are_unit_systems():
    for factory in (NaturalUnitSystem, OpenWaveUnitSystem, SIUnitSystem):
        units = factory()
        assert isinstance(units, UnitSystem), factory.__name__


# ---------------------------------------------------------------------------
# Core constants
# ---------------------------------------------------------------------------


def test_natural_core_constants():
    units = NaturalUnitSystem()
    assert units.c == 1.0
    assert units.wavelength == 1.0


def test_openwave_core_constants():
    units = OpenWaveUnitSystem()
    expected_c = constants.WAVE_SPEED * constants.RONTOSECOND / constants.ATTOMETER
    expected_lambda = constants.EWAVE_LENGTH / constants.ATTOMETER
    assert math.isclose(units.c, expected_c, rel_tol=1e-12)
    assert math.isclose(units.wavelength, expected_lambda, rel_tol=1e-12)


def test_si_core_constants():
    units = SIUnitSystem()
    assert units.c == constants.WAVE_SPEED
    assert units.wavelength == constants.EWAVE_LENGTH


# ---------------------------------------------------------------------------
# CFL bound
# ---------------------------------------------------------------------------


def test_cfl_bound_all_systems():
    for factory in (NaturalUnitSystem, OpenWaveUnitSystem, SIUnitSystem):
        units = factory()
        bound = units.dx / (units.c * math.sqrt(3.0))
        assert units.dt < bound, factory.__name__
        assert math.isclose(units.dt, 0.9 * bound, rel_tol=1e-12), factory.__name__


# ---------------------------------------------------------------------------
# Geometric constants are unit-independent
# ---------------------------------------------------------------------------


def test_geometric_constants_identical_across_systems():
    n = NaturalUnitSystem()
    o = OpenWaveUnitSystem()
    s = SIUnitSystem()
    for attr in ("A_pi", "eps_M", "N_geom", "gamma", "X_eff", "N_nu_eff"):
        v_n = getattr(n, attr)
        v_o = getattr(o, attr)
        v_s = getattr(s, attr)
        assert v_n == v_o == v_s, attr


def test_geometric_constants_expected_values():
    units = NaturalUnitSystem()
    assert math.isclose(units.A_pi, 137.036303775878, rel_tol=1e-12)
    assert math.isclose(units.N_geom, 778.8025178842, rel_tol=1e-10)
    assert math.isclose(units.eps_M, 4.1411697693e-05, rel_tol=1e-10)
    assert math.isclose(units.gamma, 1.0 / units.eps_M, rel_tol=1e-12)


def test_rho_0_is_statutory_density():
    for factory in (NaturalUnitSystem, OpenWaveUnitSystem, SIUnitSystem):
        units = factory()
        assert 1e51 < units.rho_0 < 1e53, factory.__name__


# ---------------------------------------------------------------------------
# Conversion round-trips
# ---------------------------------------------------------------------------


def test_length_round_trip():
    for factory in (NaturalUnitSystem, OpenWaveUnitSystem, SIUnitSystem):
        units = factory()
        for x in (1.0, 10.0, 100.0):
            physical = units.to_physical_length(x)
            back = physical / units.to_physical_length(1.0)
            assert math.isclose(back, x, rel_tol=1e-12), factory.__name__


def test_time_round_trip():
    for factory in (NaturalUnitSystem, OpenWaveUnitSystem, SIUnitSystem):
        units = factory()
        for t in (1.0, 10.0, 100.0):
            physical = units.to_physical_time(t)
            back = physical / units.to_physical_time(1.0)
            assert math.isclose(back, t, rel_tol=1e-12), factory.__name__


def test_cross_system_length_consistency():
    length_m = 1e-15
    for factory in (NaturalUnitSystem, OpenWaveUnitSystem, SIUnitSystem):
        units = factory()
        engine = length_m / units.to_physical_length(1.0)
        back_m = units.to_physical_length(engine)
        assert math.isclose(back_m, length_m, rel_tol=1e-12), factory.__name__


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_grid_voxels_per_lambda_validation():
    for factory in (NaturalUnitSystem, OpenWaveUnitSystem, SIUnitSystem):
        try:
            factory(grid_voxels_per_lambda=1)
        except ValueError:
            pass
        else:
            raise AssertionError(f"{factory.__name__}: expected ValueError")


def test_cfl_safety_validation():
    for factory in (NaturalUnitSystem, OpenWaveUnitSystem, SIUnitSystem):
        for bad in (0.0, -0.1, 1.5):
            try:
                factory(cfl_safety=bad)
            except ValueError:
                pass
            else:
                raise AssertionError(
                    f"{factory.__name__}: expected ValueError for {bad}"
                )


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------


def test_factory_returns_correct_types():
    assert isinstance(make_unit_system("natural"), NaturalUnitSystem)
    assert isinstance(make_unit_system("openwave"), OpenWaveUnitSystem)
    assert isinstance(make_unit_system("si"), SIUnitSystem)


def test_factory_rejects_unknown_name():
    try:
        make_unit_system("furlongs_per_fortnight")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")


# ---------------------------------------------------------------------------
# FeatureBag MRO keying
# ---------------------------------------------------------------------------


def test_feature_bag_keys_by_mro():
    bag = FeatureBag()
    units = NaturalUnitSystem()
    bag.set(units)

    assert bag.require(NaturalUnitSystem) is units
    assert bag.require(UnitSystem) is units

    # object is excluded from the key set.
    try:
        bag.require(object)
    except KeyError:
        pass
    else:
        raise AssertionError("require(object) should fail")


def test_feature_bag_overwrites_same_type():
    bag = FeatureBag()
    a1 = NaturalUnitSystem(grid_voxels_per_lambda=10)
    a2 = NaturalUnitSystem(grid_voxels_per_lambda=40)
    bag.set(a1)
    bag.set(a2)
    assert bag.require(NaturalUnitSystem) is a2
    assert bag.require(UnitSystem) is a2


def test_feature_bag_distinct_concrete_types_coexist():
    bag = FeatureBag()
    a = NaturalUnitSystem()
    b = OpenWaveUnitSystem()
    bag.set(a)
    bag.set(b)
    # Concrete keys are independent.
    assert bag.require(NaturalUnitSystem) is a
    assert bag.require(OpenWaveUnitSystem) is b
    # The shared base points to the most recent write.
    assert bag.require(UnitSystem) is b


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def main() -> int:
    tests = [
        test_natural_default_construction,
        test_openwave_default_construction,
        test_si_default_construction,
        test_all_implementations_are_unit_systems,
        test_natural_core_constants,
        test_openwave_core_constants,
        test_si_core_constants,
        test_cfl_bound_all_systems,
        test_geometric_constants_identical_across_systems,
        test_geometric_constants_expected_values,
        test_rho_0_is_statutory_density,
        test_length_round_trip,
        test_time_round_trip,
        test_cross_system_length_consistency,
        test_grid_voxels_per_lambda_validation,
        test_cfl_safety_validation,
        test_factory_returns_correct_types,
        test_factory_rejects_unknown_name,
        test_feature_bag_keys_by_mro,
        test_feature_bag_overwrites_same_type,
        test_feature_bag_distinct_concrete_types_coexist,
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