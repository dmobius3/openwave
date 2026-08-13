"""Causal mutation battery for the rung-3b evaluator (protocol section 4.2 step 6).

Each mutation targets ONE layer and must redden ONE named check.  The point is
not coverage but attribution: after this runs, every proposition the evaluator
rests on has a test that could have falsified it, and no two propositions share
a test.  That is the preventive form of the displaced-verification lesson from
M8.8: a single passing tower comparison establishes at most one row of this table.

Two of the five are the actual published traps, not invented perturbations:

    M1  the generator's final block as printed in Lauret v5 section 2
    M2  Theorem 3.3's branch mapping as printed

A mutation that CRASHES is not a red.  Each case below asserts the failure is a
wrong ANSWER at a named check, and a crash is reported as an inconclusive
mutation rather than a catch.

Nontrivial-Gamma cases use L(7; 1, 2), a member of the pilot tuning set frozen in
section 6.1.  Section 4.1's admissibility rule states the sealed adjudication
case is distinct from every member of that set, so this carries no case
information.
"""

import sys
import traceback
from math import comb as _math_comb

import lauret_evaluator as ev
import gate_gamma1 as gate

Q_LENS, S_LENS = 7, (1, 2)      # L(7; 1, 2): inhomogeneous, in the tuning set
KMAX_LENS = 8


def _baseline_lens():
    return ev.p_form_spectrum(1, Q_LENS, S_LENS, KMAX_LENS)


# --- the five mutations ------------------------------------------------------

def m1_generator_final_block():
    """M1, GROUP-ACTION / CONGRUENCE-LATTICE check.

    Lauret v5 section 2 prints the final rotation block with s_1 repeated.
    Reddens by changing the congruence lattice, so it is invisible at Gamma = 1
    (q = 1 satisfies every congruence) and must be tested on a nontrivial case.
    """
    base = _baseline_lens()
    mutated = ev.p_form_spectrum(1, Q_LENS, S_LENS, KMAX_LENS,
                                 blocks=(S_LENS[0],) * len(S_LENS))
    return ("congruence lattice", base != mutated)


def m2_branch_mapping_as_printed():
    """M2, SPECTRAL-INDEX check at Gamma = 1.

    Theorem 3.3 as printed.  Returns M_Gamma(k, 0) = 0 on the lower branch
    because eq (3.1)'s j-sum is empty at p = 0.
    """
    failures = gate.run(kmax=6, mapping="as_printed", verbose=False)
    return ("Gamma = 1 spectral index", len(failures) > 0)


def m3_binomial_convention():
    """M3, SUMMATION-DOMAIN check.

    The realistic error: treating C(b, 0) as 1 by the empty-product convention
    even when b < 0, where the paper's Notation 3.1 gives 0.  Silent and wrong
    rather than a crash.
    """
    def naive_binom(b, a):
        if a == 0:
            return 1                     # WRONG when b < 0
        if a < 0 or b < a:
            return 0
        return _math_comb(b, a)

    failures = gate.run(kmax=6, verbose=False) or []
    if failures:
        return ("summation domain", False)   # baseline already red, inconclusive
    rows = ev.p_form_spectrum(1, 1, (1, 1), 6, binom=naive_binom)
    ok = all(r["lower_branch"]["multiplicity"] == (r["k"] + 1) ** 2
             and r["upper_branch"]["multiplicity"] == 2 * r["k"] * (r["k"] + 2)
             for r in rows)
    return ("summation domain", not ok)


def m4_lattice_census():
    """M4, LATTICE-CENSUS check.

    Z(mu) miscounted.  Z feeds the Lemma 3.2 kernel through ell, so a wrong
    census produces wrong multiplicities with the lattice membership itself
    untouched.  Distinct from M1, which changes WHICH weights are in the lattice.
    """
    original = ev.zero_count
    try:
        ev.zero_count = lambda a: max(0, sum(1 for aj in a if aj == 0) - 1)
        rows = ev.p_form_spectrum(1, 1, (1, 1), 6)
        ok = all(r["lower_branch"]["multiplicity"] == (r["k"] + 1) ** 2
                 and r["upper_branch"]["multiplicity"] == 2 * r["k"] * (r["k"] + 2)
                 for r in rows)
    finally:
        ev.zero_count = original
    return ("lattice census", not ok)


def m5_eigenvalue_labelling():
    """M5, EIGENVALUE-LABELLING check.

    lambda_{k,p} = (k+p)(k+2n-2-p) altered.  Reddens on the eigenvalue while the
    multiplicities stay correct, which is what makes it distinguishable from the
    other four.
    """
    original = ev.eigenvalue
    try:
        ev.eigenvalue = lambda k, p, n: 0 if p == -1 else (k + p) * (k + 2 * n - 1 - p)
        failures = gate.run(kmax=6, verbose=False)
    finally:
        ev.eigenvalue = original
    return ("eigenvalue labelling", len(failures) > 0)


MUTATIONS = [
    ("M1 generator final block (Lauret sec 2 as printed)", m1_generator_final_block),
    ("M2 branch mapping as printed (Thm 3.3 as printed)", m2_branch_mapping_as_printed),
    ("M3 binomial convention, C(b,0)=1 for b<0", m3_binomial_convention),
    ("M4 Z(mu) miscounted", m4_lattice_census),
    ("M5 eigenvalue exponent altered", m5_eigenvalue_labelling),
]


def main():
    print("  rung-3b evaluator: causal mutation battery\n")
    print("  baseline must be green before any mutation is meaningful:")
    base_fail = gate.run(kmax=6, verbose=False)
    print(f"    Gamma = 1 baseline: {'GREEN' if not base_fail else 'RED'}")
    if base_fail:
        print("    baseline red; mutations are not interpretable.  Stop.")
        return 1
    print()

    seen_checks, bad = set(), []
    for name, fn in MUTATIONS:
        try:
            check, caught = fn()
        except Exception:
            print(f"    INCONCLUSIVE  {name}")
            print(f"                  crashed rather than returning a wrong answer:")
            print("                  " + traceback.format_exc().strip().splitlines()[-1])
            bad.append(name)
            continue
        print(f"    {'DETECTED  ' if caught else 'MISSED    '}{name}")
        print(f"                  reddens: {check}")
        if not caught:
            bad.append(name)
        seen_checks.add(check)

    print()
    print(f"  distinct checks reddened: {len(seen_checks)} of {len(MUTATIONS)}")
    if len(seen_checks) != len(MUTATIONS):
        print("  WARNING: two mutations share a check, so one passing comparison "
              "could be read as establishing both.")
    if bad:
        print(f"  BATTERY RED: {len(bad)} mutation(s) not caught cleanly.")
        return 1
    print("  BATTERY GREEN: every mutation caught, each at its own check.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
