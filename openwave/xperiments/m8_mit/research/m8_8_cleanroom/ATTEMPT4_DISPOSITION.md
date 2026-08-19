# M8.8 clean-room attempt 4: RETIRED PRE-REVEAL

> **Disposition, 2026-08-19.** Retired before § 8 step 6 for multiple gate-contract
> nonconformances. Values unadjudicated, answer packet sealed, case unconsumed. Outputs
> preserved uncompared at `9aa740ea`.

## Why this retirement is different from the first three

Attempts 1, 2 and 3 each exposed a genuine ambiguity or a missing piece of architecture, and
in each case repairing the authored instruction was the right response. Attempt 4 mostly did
not expose anything new. It contradicted requirements the instruction and the frozen protocol
already state explicitly:

- a coverage row marked not-checkable must be genuinely not pre-checkable, not merely
  inconvenient to check;
- the production expected gate set comes from the frozen manifest, not a duplicated literal;
- mutations must actually run;
- manual pass and red attestations are refused;
- an exact gate may not be replaced by a weaker proxy because the exact route is expensive;
- when a frozen requirement cannot be implemented, the authorized outcome is to stop.

## The nonconformances

| # | Finding | Where |
| --- | --- | --- |
| 1 | The pre-implementation exemptions are not defensible. The fundamental representation, the symmetric-power construction, the decomposition, the quaternion-to-`SU(2)` convention and the theorem-side representation checks are all marked not-checkable because they require representation construction. The instruction expressly permits validation scripts to construct and check internal objects provided they do not compute the target or populate the raw output, and these facts are target-blind and checkable from the packets. This is precisely the class the previous repair existed to move left of the freeze | manifest § 7 |
| 2 | The production coverage checker hard-codes its expected gate set instead of reading the frozen manifest, which is the duplicated-list implementation the instruction prohibits by name, and for the stated reason: a manual copy can omit the same gate the runtime omits | `torsion.py` |
| 3 | A mutation verdict is assigned as a literal with a comment reading "by design"; nothing is mutated and nothing is observed to fail | `torsion.py` GATE-M05-mut |
| 4 | A gate verdict is an attestation: a comment saying the path is verified by construction, followed by an unconditional pass, where the protocol requires the harness to inspect the route-native intermediates and mutation-test the dependencies linking them to the reported value | `torsion.py` GATE-R03 |
| 5 | The exact saturation requirement is replaced by nontrivial-irrep acyclicity plus Euler-characteristic reasoning, on the stated grounds that the exact route is expensive. The protocol forbids exactly this substitution and supplies the counterexample showing per-irrep acyclicity and ranks all surviving while degree-1 integral homology is wrong. The frozen manifest had itself promised exact certificates in production, so a frozen commitment was weakened after the manifest became immutable | `torsion.py` GATE-M05 against protocol § 9 |
| 6 | GATE-M04 does not test its own claim. Its pass predicate is the same condition as GATE-M06, and its mutation declares red because an entry changed rather than because augmented homology was recomputed and failed, though the manifest promised recomputation | `torsion.py` |
| 7 | The environment record states that no third-party package is imported, while a validation artifact imports NumPy and uses floating-point rank and determinant work during validation. Not the retiring defect, but the environment record describes the reproduction as a whole rather than the production script alone | `ENVIRONMENT.md` |

Further mutation-semantic seams were noted and are recorded as supporting evidence rather
than litigated: the unitarity mutation tests departure from the identity form while the
declared gate permits any invariant positive-definite form, and the row-identity mutation
resolves red because two baseline signatures differ rather than by running the mutated
assignment through the declared gate.

## The structural finding, which is the durable one

The coverage enforcement added for this attempt verifies that each registered gate has a
result key, a mutation key, and the value `RED`. A hard-coded literal satisfies all three. So
the mechanism built to guarantee that mutations exist cannot distinguish an executed mutation
from an asserted one: the same defect class it was created to close, one level up. This
project has recorded the identical failure once before, when a certificate wrote its premises
as literals rather than as outcomes of the checks above it.

That finding is recorded rather than repaired, for the reason in the next section.

## The process ruling: the instruction is not amended again

`TASK.md` is unchanged for attempt 5, at the identical cleared bytes
`e3d9b90861bb81862843988e8bd5da925b4d48bc48c0d2335becd3137df9cb17`, already committed and
already semantically cleared. Adding another paragraph restating requirements the document
already makes explicit would overfit the instruction to each failed implementer rather than
close a real seam.

That makes attempt 5 a different and more informative test than its predecessors: **can a
fresh implementer comply with the mature protocol as written**, rather than can the protocol
be amended until some implementation passes.

**Stopping rule, set in advance.** If attempt 5 fails on another direct violation of
requirements the instruction already states clearly, commissioning stops and the question
becomes whether this clean-room workflow is practical for this task at all, rather than
whether one more instruction revision would help.

## Preserved

| Object | State |
| --- | --- |
| Attempt-4 outputs | `9aa740ea`, pushed, no pull request, no comparison ever run |
| Attempt-4 room | `~/Desktop/OpenWave/M8.8_ATTEMPT4_RETIRED/`, intact |
| Attempt-5 room | rebuilt fresh; the three technical inputs re-extracted from the frozen commits and the task file taken from the committed cleared bytes rather than from any room |
| Canonical answer packet | SEALED, never opened |
