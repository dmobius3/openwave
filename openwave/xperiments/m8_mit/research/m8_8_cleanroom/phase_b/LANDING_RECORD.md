# M8.8 Phase B: G-D05 routed through the frozen path, landed unread against the answer packet

> **Claims no result.** The answer packet has not been opened, no comparison has run, no
> verdict exists. Lands before § 8 step 6. The three prior Phase B landings, `f441e0ec`,
> `a923c554` and `bfdca08f`, are preserved in history.

## What landed

| File | SHA-256 | Role |
| --- | --- | --- |
| `run_phase_b.py` | `bcbf6543a98a81f765e28adccbd1841ec7e0ea3783184c6e251a4372a1cdec52` | the qualifier |
| `MUTATION_RESULTS.json` | `f0a8b91038a3f7eebb7ac6a3accf85452c4eba492d87622b6ad51df895411a67` | the execution record |
| `QUALIFICATION_RECORD.md` | `a11954b9d8b04d6d2c14fe32064565b0df6601c96fb81ee39f1b4a393ea4e52a` | the qualifier's summary |

## The blocking finding, and why the previous arm could not fail

The maintainer's review established that G-D05 as landed at `bfdca08f` computed
`det(I)/(det(I)·det(I))` in the qualifier's own code and compared the result with the Phase A
value. That comparison never entered the frozen `compute_torsion_sq`, which builds its three
minors internally and exposes no parameter for them, so the "input boundary" the record named
existed only in the qualifier's re-typed copy of the formula.

The consequence is worse than a fidelity defect. Simulated here before the repair: a Phase A
that returns a stored value and consumes none of its determinants is reddened by the old arm
and detected by the new one. The gate declares that derivation artifacts are consumed rather
than bypassed, and the old arm passed the bypass it exists to catch.

The commissioner's earlier verification checked the mutation's description against the
manifest declaration and did not check that the declared boundary existed in the frozen
object. That is what let it land twice.

## What the repair does

The qualifier substitutes `det_gc` on the loaded frozen module for the duration of one call to
the frozen `compute_torsion_sq`, returning the determinant of the same-dimension identity for
every minor the frozen code hands it, and restores the original in a `finally` block. Both
compared values are returned by the frozen function; the qualifier computes neither. The
intercepted count is asserted to equal exactly three, which is structural rather than
empirical: under identity substitution no determinant is zero, so the `M3` and `M1` minor
searches each break on their first iteration and the central `M2` determinant is unconditional.

The count is deliberately not compared against the unmutated path. That comparison was
considered and rejected: the unmutated search length is data-dependent, running three, five or
seven determinants depending on the irrep, so equality holds for V1 by coincidence and fails
for four of the eight acyclic irreps.

## Also in this landing

The four convention gates now execute through the frozen `validate_fixture.compute_torsion`
rather than a qualifier-local copy. That module builds its fixture inside `main()` and exposes
no importable builder, so the qualifier reconstructs it; the reconstruction is verified faithful
by reproducing the frozen module's own fixture baseline exactly. The qualification record gains
a reproduction route, absent before. Four metadata fields that were written as literals now
carry values returned by the checks.

## Commissioner enumeration, pre-registered before the run

| Check | Result |
| --- | --- |
| Phase A | 13 of 13 exact against the Addendum 1 table; no generated files beneath the frozen directory |
| G-D05 reaches the frozen path | instrumented the real run from above: five `compute_torsion_sq` invocations, exactly one with the substitution active, consuming three determinants |
| That check can fail | reverting G-D05 to the `bfdca08f` arm drops the substituted-call count to zero and the check goes red |
| The exactly-three assert is load-bearing | forcing a fourth determinant inside the substitution window fires it and the run dies nonzero |
| Exception safety | a fault induced inside the mutated call propagates as a failure rather than a silent pass |
| G-T03a to G-T03d | routed through the frozen fixture function; the reconstructed fixture reproduces the frozen module's own baseline |
| `declared_mutation` fidelity | 19 of 19 verbatim-identical to the frozen manifest's § 4 table, by diff |
| Outcomes | all nineteen red; coverage proven before and after execution |
| Reproduction | the documented route run verbatim from a clean directory exits 0 and regenerates `MUTATION_RESULTS.json` byte-identical |

## Two things this record does not claim

The literal-field repair is partial. `pre_execution_set_equality` and
`manifest_parsed_at_runtime` are computed comparisons, but `phase_a_hashes_verified_pre` and
`parser_self_test_passed` now come from functions whose terminal statement is still
`return True`, guarded by an earlier exit. The fields cannot read true while their checks fail,
so the record is not misleading, but the shape is displaced one level rather than removed.

Per-gate evidence is not uniform, and it decomposes by what changed. G-D05 changed and was
independently instrumented and mutation-attacked. The four convention handlers changed and were
reverified against the frozen fixture path. The remaining fourteen handlers are unchanged since
`bfdca08f` and inherit the maintainer's pair-by-pair read of eighteen of nineteen declared
against implemented mutations at that commit, where G-D05 was the sole blocking item; here they
are reconfirmed only as red with their declared mutations diffed verbatim against the manifest.
That inheritance is stated rather than restated as fresh commissioner instrumentation, because
nineteen independent instrumentations were not performed.
