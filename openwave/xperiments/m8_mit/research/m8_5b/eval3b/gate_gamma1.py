"""The Gamma = 1 validation gate for the rung-3b evaluator.

Protocol section 4.2 step 4: the evaluator "is first checked at Gamma = 1, where
it must recover the unit-S^3 one-form tower of rung 2.  Failing that check voids
the rung before it is applied."

The reference tower is pinned in `../ref_pin_sheet.md` section C.3 from three
independent sources:

    exact 1-forms    lambda_{k,0} = k(k+2)    multiplicity (k+1)^2
    coexact 1-forms  lambda_{k,1} = (k+1)^2   multiplicity 2k(k+2)

This gate checks a RANGE of k rather than the single famous point
(lambda = 4, multiplicity 6).  One point can be hit by coincidence; several
consecutive entries establish the indexing PATTERN, which is what correction B
is about.  Gamma = 1 carries no case information: q = 1 makes every congruence
trivially satisfied, and the rung-2 tower is in neither packet.

Exit 0 if every checked entry agrees, 1 otherwise.
"""

import sys

from lauret_evaluator import p_form_spectrum, FORMULA_PROVENANCE

KMAX = 12
Q_TRIVIAL = 1
S_TRIVIAL = (1, 1)          # n = 2; with q = 1 the parameters are immaterial
P_ONE_FORMS = 1


def reference_tower(k):
    """The pinned unit-S^3 one-form tower of rung 2."""
    return {
        "exact": {"eigenvalue": k * (k + 2), "multiplicity": (k + 1) ** 2},
        "coexact": {"eigenvalue": (k + 1) ** 2, "multiplicity": 2 * k * (k + 2)},
    }


def run(kmax=KMAX, mapping="corrected", verbose=True):
    rows = p_form_spectrum(P_ONE_FORMS, Q_TRIVIAL, S_TRIVIAL, kmax,
                           mapping=mapping)
    failures = []
    if verbose:
        print(f"  Gamma = 1, n = 2, p = 1, branch mapping = {mapping!r}\n")
        print("    k | exact lam  mult   ref | coexact lam  mult   ref | ok")
        print("    --+-----------------------+------------------------+---")
    for row in rows:
        k = row["k"]
        ref = reference_tower(k)
        lo, hi = row["lower_branch"], row["upper_branch"]
        ok_lo = (lo["eigenvalue"] == ref["exact"]["eigenvalue"]
                 and lo["multiplicity"] == ref["exact"]["multiplicity"])
        ok_hi = (hi["eigenvalue"] == ref["coexact"]["eigenvalue"]
                 and hi["multiplicity"] == ref["coexact"]["multiplicity"])
        if not (ok_lo and ok_hi):
            failures.append((k, lo, hi, ref))
        if verbose:
            print(f"    {k:2d} | {lo['eigenvalue']:9d} {lo['multiplicity']:5d} "
                  f"{ref['exact']['multiplicity']:5d} |"
                  f" {hi['eigenvalue']:11d} {hi['multiplicity']:5d} "
                  f"{ref['coexact']['multiplicity']:5d} |"
                  f" {'ok' if ok_lo and ok_hi else 'FAIL'}")
    return failures


def main():
    print("  rung-3b evaluator, Gamma = 1 validation gate")
    print(f"  provenance: {FORMULA_PROVENANCE}\n")
    failures = run()
    print()
    if failures:
        print(f"  GATE RED: {len(failures)} of {KMAX} entries disagree with the "
              "pinned rung-2 tower.")
        print("  Section 4.2 step 4: the rung is VOID before it is applied.")
        return 1
    print(f"  GATE GREEN: {KMAX} consecutive entries agree with the pinned "
          "rung-2 tower on both branches.")
    print("  Eigenvalues AND multiplicities checked; the indexing pattern is "
          "established over a range, not at one point.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
