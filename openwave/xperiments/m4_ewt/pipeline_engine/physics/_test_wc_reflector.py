"""
Tests for the WC reflector attributes (plan item 1.5).

The reflector interface adds three attributes to WC:
    reflect_coeff_long, reflect_coeff_trans, phase_shift.

They are attributes, not behaviours: nothing reads them yet. These
tests pin the defaults, the unitarity of the default pair, and the
constructibility of non-unitary pairs (which are legal attribute
values, since the interface does not enforce unitarity at runtime).

Run from the project root:
    python -m openwave.xperiments.m4_ewt.pipeline_engine.physics._test_wc_reflector
"""

import sys
import traceback

from .wc_types import WC, WCState


def test_default_reflector_attributes():
    """
    A WC constructed with positional args only carries the documented
    reflector defaults.

    Mutation caught: a default changed to something other than
    (1.0, 0.0, 0.0).
    """
    wc = WC(1.0, 2.0, 3.0)
    assert wc.reflect_coeff_long == 1.0, wc.reflect_coeff_long
    assert wc.reflect_coeff_trans == 0.0, wc.reflect_coeff_trans
    assert wc.phase_shift == 0.0, wc.phase_shift


def test_default_reflector_is_unitary():
    """
    The default pair (1.0, 0.0) satisfies unitarity trivially:
    reflect_coeff_long^2 + reflect_coeff_trans^2 = 1.0.

    Mutation caught: a default pair that violates unitarity, e.g.
    (1.0, 1.0) or (0.5, 0.5).
    """
    wc = WC(0.0, 0.0, 0.0)
    total = wc.reflect_coeff_long**2 + wc.reflect_coeff_trans**2
    assert abs(total - 1.0) < 1e-12, total


def test_custom_reflector_attributes_persist():
    """
    Custom reflector values survive construction.

    Mutation caught: a field that is silently dropped or overwritten
    by __post_init__ or a base class.
    """
    wc = WC(
        x=1.0,
        y=2.0,
        z=3.0,
        reflect_coeff_long=0.8,
        reflect_coeff_trans=0.6,
        phase_shift=0.25,
    )
    assert wc.reflect_coeff_long == 0.8
    assert wc.reflect_coeff_trans == 0.6
    assert wc.phase_shift == 0.25
    # This pair is unitary: 0.8^2 + 0.6^2 = 0.64 + 0.36 = 1.0.
    total = wc.reflect_coeff_long**2 + wc.reflect_coeff_trans**2
    assert abs(total - 1.0) < 1e-12, total


def test_non_unitary_pair_is_legal():
    """
    The interface does not enforce unitarity at runtime. A pair whose
    squares do not sum to 1 is a legal attribute value; the scattering
    operator that will read it (plan Section 7, item 2.2 B2d) is what
    must enforce unitarity.

    This is the arena note for the tests above: they check the default
    pair and the persistence, not a runtime rule.
    """
    wc = WC(0.0, 0.0, 0.0, reflect_coeff_long=0.5, reflect_coeff_trans=0.5)
    total = wc.reflect_coeff_long**2 + wc.reflect_coeff_trans**2
    # 0.25 + 0.25 = 0.5, not 1.0. Construction does not raise.
    assert abs(total - 0.5) < 1e-12, total


def test_wc_state_keeps_reflector_attributes_on_each_center():
    """
    WCState is a list of WC. Each center carries its own reflector
    attributes; the state does not shadow or overwrite them.

    Mutation caught: a WCState that ignores or replaces its centers'
    reflector fields.
    """
    wc_a = WC(0.0, 0.0, 0.0, reflect_coeff_long=1.0, reflect_coeff_trans=0.0)
    wc_b = WC(1.0, 0.0, 0.0, reflect_coeff_long=0.6, reflect_coeff_trans=0.8)
    state = WCState(centers=[wc_a, wc_b])
    assert state.centers[0].reflect_coeff_long == 1.0
    assert state.centers[1].reflect_coeff_long == 0.6
    assert state.centers[1].reflect_coeff_trans == 0.8


def main() -> int:
    tests = [
        test_default_reflector_attributes,
        test_default_reflector_is_unitary,
        test_custom_reflector_attributes_persist,
        test_non_unitary_pair_is_legal,
        test_wc_state_keeps_reflector_attributes_on_each_center,
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
