# M8.8 clean-room attempt 4: the output, landed unread against the answer packet

> **Claims no result.** The answer packet has not been opened, no comparison has run, and no
> verdict exists. This commit discharges § 3's commit-before-unseal ordering, carried by
> ancestry.
>
> **This record states findings rather than deferring them.** The previous landing record
> flagged gate sufficiency as a review question when it was machine-checkable; the
> enumeration below is the correction to that.

## What landed

| File | SHA-256 | Role |
| --- | --- | --- |
| `TASK.md` | `e3d9b90861bb81862843988e8bd5da925b4d48bc48c0d2335becd3137df9cb17` | the operational instruction, committed before this run started |
| `METHOD_AND_GATE_MANIFEST.md` | `0160305315cfb4f50d0650cdf2cd27773329197bb0fd24d718d9b3cecef2af8d` | § 8 step 4, with registries, coverage table and the freeze line |
| `validate_pre_impl.py` | `a86e596bf82e20714558023134f2cf3854e1e00c04b181875a19c40e3785b968` | pre-implementation validation artifact |
| `validate_manifest.py` | `bdb408b30d4aeaf497021edda2b4a5deb91b0a825d16272e9caced5d86283000` | the registry/coverage set-equality validator |
| `torsion.py` | `866cd003d072b1583a0f79a1c55a8b5b74159f5cb5f12dc32eb3688ed1d4bf25` | production implementation |
| `ENVIRONMENT.md` | `d31867c3c4f81d2aa105e8efe0c7ac10cf5a639ae226cd7a04b3068bc48f813e` | environment record |
| `RAW_OUTPUT.json` | `cced827a8940e680c0418ca22a2c273575d8eb77626ef77d6a86d4f4575be272` | the § 5.5 raw output, carrying derivation_artifacts and gate_results |
| `CONSULTED_FILES.md` | `52ce859175c1d6c09831ba766f2ba1e3f60ac15046aa6c6989c3f9f1f8c87fb3` | consulted-files manifest |

## What this attempt got right

The architecture the last three rounds built is present and working. Ordering is
`validate_pre_impl.py` 15:24, `validate_manifest.py` 15:26, manifest final 15:27, production
15:49, outputs 15:50 to 15:51. The manifest carries a construction registry (9), a convention
registry (10) and a pre-reveal gate registry (17), and a coverage table with all 36 rows.
`validate_manifest.py` genuinely reads the manifest and proves exact set equality; run
independently it reports 36 against 36 and exits 0. The production coverage check blocks with
`sys.exit(1)` rather than warning. All four seeded inputs were byte-unchanged. Both
validators and the production run rerun from the committed bytes at exit 0, and
`RAW_OUTPUT.json` regenerates byte-identical.

## Findings, enumerated rather than deferred

| # | Finding | Evidence |
| --- | --- | --- |
| 1 | **A mutation verdict is a hard-coded literal.** `gate_results['GATE-M05-mut'] = 'RED'` with the comment "by design"; no mutation is executed | `torsion.py:1330` |
| 2 | **A gate verdict is pure attestation.** `gate_results['GATE-R03'] = 'PASS'` preceded only by a comment reading "Verified by construction"; no check runs | `torsion.py:1186` |
| 3 | **GATE-M05 substitutes a check the protocol names as insufficient.** Production replaces the exact saturation certificate with per-irrep acyclicity plus the Euler characteristic. The protocol's own § 9 states that the `6 − 5z` construction leaves per-irrep acyclicity, exact ranks and an `im ∂₃` certificate all passing while the damage lands in degree 1 | `torsion.py:1313-1327` against the protocol's § 9 row |
| 4 | **The manifest promised what production did not deliver.** The coverage table records GATE-M05 as "approx; exact in production"; production is explicitly approximate. Under the immutability clause this mismatch called for conform-or-STOP | manifest § 7 against `torsion.py:1316` |
| 5 | **Two gates evaluate one predicate.** `m04_ok = aug_d1_all_zero` is the same condition as GATE-M06, so GATE-M04 does not establish the augmented homology it claims | `torsion.py:1301-1334` |
| 6 | **A mutation arm is vacuous.** GATE-M04's mutation asserts that `x + 1 != x`, which is arithmetic rather than a test of the gate | `torsion.py:1309` |
| 7 | **The environment record overstates.** It says "No `numpy` ... is imported"; `validate_pre_impl.py` imports numpy. Production is genuinely standard-library only, so the claim is true of the implementation and false as written | `validate_pre_impl.py` |

## The structural finding, which outlives this attempt

The coverage enforcement added for this run verifies that each registered gate has a result
key, a mutation key, and the mutation value `RED`. **A hard-coded literal satisfies every one
of those conditions.** So the mechanism built to guarantee that mutations exist cannot
distinguish an executed mutation from an asserted one, which is the same defect class it was
created to close, one level further up.

This is also the failure this project has already recorded once, in its own history: a
certificate that recorded its premises as satisfied because they were Python literals in the
write rather than outcomes of the checks above them. Any repair should make the mutation
evidence an execution record, so that a verdict cannot be written without the mutation having
actually run and been observed to fail.

## What is deliberately absent

No adjudication record, no comparison, no orientation selection, no answer packet.
