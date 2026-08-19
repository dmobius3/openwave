# M8.8 clean-room attempt 5: the output, landed unread against the answer packet

> **Claims no result.** The answer packet has not been opened, no comparison has run, no
> verdict exists. Ordering is carried by ancestry. Findings are enumerated below rather than
> deferred.

## What this attempt fixed, and it is substantial

Attempt 5 ran on the identical instruction as attempt 4, with no amendment, and independently
corrected most of what retired that run.

| Previously retired for | Attempt 5 |
| --- | --- |
| Representation constructions exempted from pre-implementation validation | validated, in `validate_representations.py`: homomorphism, Weyl characters, orthogonality, idempotent projector, distinct row signatures, nine irreps at the right dimensions. The Galois row explicitly records `Sym³` verified self-conjugate, which is the exact fact whose absence retired an earlier attempt |
| Exact saturation replaced by a weaker proxy | done exactly, in `validate_saturation.py`: unimodular minors of determinant `±1` for `im ∂₁`, `im ∂₂` and `im ∂₃`, which is the certificate the protocol demands in every degree |
| Verdicts written as hard-coded literals | none. A sweep for literal `PASS`/`RED` assignments across all nine scripts returns nothing |
| Dishonest not-checkable exemptions | two rows only, both because the torsion output is itself the target and cannot be checked target-blind. That reason is sound |
| Set equality unproven | `validate_manifest.py` reads the manifest and proves registry against coverage-table equality |

Ordering is clean: seven validation artifacts from 16:33 to 17:08, manifest final 17:12,
production 17:13, output 17:13. The firewall held, all four seeded inputs byte-unchanged.
Production reruns from the committed bytes at exit 0 and `RAW_OUTPUT.json` regenerates
byte-identical. A9 clean, standard library only.

## The finding that retires it

**Fifteen of nineteen declared gate mutations were never executed.**

Each `gate_results` entry carries `outcome`, `artifact`, and a `mutation` field. That field is
a *description* of the mutation, not a record that one ran and was observed to fail. Running
all seven validation artifacts, only `validate_fixture.py` executes mutations, four of them,
the convention mutations. The other six execute none, and `compute_torsion.py` contains no
gate identifiers, no coverage logic, and no mutation execution at all.

So the run reports every gate passing while fifteen of the declared mutation arms exist only
as prose. The protocol's § 9 states the requirement in bold and refuses this explicitly:
every gate carries a runnable mutation that must redden it, with enforced coverage and a
nonzero exit, and manual attestation is not accepted. The instruction repeats it: the
mutation's red outcome must be recorded, and incomplete coverage must force a nonzero exit.

**There is no coverage checker.** The requirement that the run be unable to report success
while coverage is incomplete is not implemented at any level, and the run exits 0.

This is materially the defect that retired the previous attempt, expressed as a description
rather than as a literal. A string standing in for an execution is the same object either
way.

## What is deliberately absent

No adjudication record, no comparison, no orientation selection, no answer packet.
