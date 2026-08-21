# M8.8 Phase B: mutation-semantic fidelity repaired, landed unread against the answer packet

> **Claims no result.** The answer packet has not been opened, no comparison has run, no
> verdict exists. Lands before § 8 step 6. Both prior Phase B landings, `f441e0ec` and
> `a923c554`, are preserved in history.

## What landed

| File | SHA-256 | Role |
| --- | --- | --- |
| `run_phase_b.py` | `2a741b5403bf7f30e4cf96ff4c4e9016038f520b41c825dfe11a0a0e568c24d3` | the qualifier: parses the § 4 registry at runtime and executes all nineteen declared mutations |
| `MUTATION_RESULTS.json` | `490938f89b0dbcfd9b6b9abd7c01d310624a415ae9345f983b929596d02fb0ae` | the execution record, carrying declared and implemented mutation text per gate |
| `QUALIFICATION_RECORD.md` | `2133a653e3dc76d128e757da58228824891d5d4fce4242f3cfb82e190211549a` | the qualifier's summary |

## The two fidelity defects are cured, and one of them had a real causal path

The prior landing executed, for two gates, a mutation related to the declared one rather than
the declared one itself. Both now execute what the manifest declares, at the boundary it
declares.

**G-D04**, declared as a swap of `φ ↦ 1−φ` on the characters of a Galois pair. Previously a
substitution of one already-computed value for another at the comparison, which tested the
equality checker rather than the Galois action. Now the Galois map is applied to the
representation matrices that carry those characters, the torsion is recomputed through the
same consumed path, and the resulting pair no longer satisfies the relation.

The instruction offered an honest stop here, because the frozen production path pairs Galois
partners from a literal list and computes torsion before that check runs, so it was genuinely
open whether the declared mutation could reach the relation at all. It can. The escape was
not needed, and the fact that it was available is why the answer means something.

**G-D05**, declared as replacing the torsion formula's input with identity matrices and
verifying the output changes. Previously the representation matrices were replaced and the
computation aborted as non-acyclic, which is not a changed output. Now the three determinant
sub-matrices are replaced at the formula's own input boundary and the output changes to a
different value. No record in the set now reports an abort in place of an outcome.

## Commissioner enumeration

| Check | Result |
| --- | --- |
| Phase A | 13 of 13 exact; no generated files beneath it, before or after my testing |
| Declared against implemented, all nineteen | read pair by pair; each implemented mutation performs what its declaration states, at the declared boundary |
| `declared_mutation` fidelity | verbatim from the frozen manifest for all nineteen |
| Coverage | parsed == implemented handlers == executed records, proven before and after execution |
| Outcomes | all nineteen red; none records an abort as a changed outcome |
| Reproduction | clean run exits 0 and regenerates `MUTATION_RESULTS.json` byte-identical |
| Linkage, retested after the repair | dropping one identifier from the parsed set below the integrity layer drives a coverage failure, not a hash failure |

## One sentence the record no longer overstates

The prior summary asserted that no mutation was skipped, substituted, or narrowed, which was
false for the two gates above. It now claims only that no mutation was skipped, and points at
the paired `declared_mutation` and `implemented_mutation` fields so any future divergence is
visible in the machine-readable record rather than resting on prose.
