# M8.9 S1b: implementation qualification of the frozen analytic-subspace control

> **You are qualifying an INSTRUMENT, not producing a scientific result.** Reproduce every gate and
> every mutation in the frozen contract, attack them for vacuity, and STOP on any contradiction.
> You may not redesign a threshold, relax a gate, or adjust anything in response to what the live
> target produces. If the contract cannot be implemented as written, that is your finding and you
> report it rather than working around it.

## 1. The contract governs, and it is append-only

[`contract/S1B_DECISION_RULE.md`](contract/S1B_DECISION_RULE.md), frozen region above the boundary
comment, SHA-256

    c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297

**Verify that hash before you write any code**, with the command the document carries. If it does
not match, stop and report; do not proceed on a document you cannot authenticate. Where this brief
and the contract appear to differ, the CONTRACT wins and the difference is reported rather than
resolved locally.

## 2. Why you and not the contract's author

The author wrote the rule across seven adversarial design passes. Several defects were caught in
the author's own gates: a conditioning test that the step before it guaranteed, a mutation arm whose
mutation could not change its object, a cover-space matrix applied to a quotient operator, an
adjudicating witness in one norm compared against uncertainty terms in another. Two more surfaced in
the freeze machinery itself. The author is therefore not the right party to certify that the
document can be implemented without judgment calls. That is your job.

## 3. Deliverables, in this order

**Q0. Room and contract provenance, before any code.** Run `room_import_gate.py` first; it is a
launch precondition and a failure blocks everything. It checks three independent things, and none
of them is the process exit code:

- **sentinel**, every REQUIRED module must emit `IMPORT_COMPLETE:<name>` AFTER `import_module()`
  returns. A module that calls `sys.exit(0)` during import has NOT imported, and exiting 0 without
  the sentinel is RED. This gate was itself false-green once for exactly that reason;
- **origin**, every REQUIRED module's resolved `__file__` must live under this room, so a missing
  dependency cannot silently satisfy itself from the live checkout or another `PYTHONPATH` entry;
- **manifest**, `ROOM_MANIFEST.json` pins a SHA-256 for EVERY file in this room, this brief and
  the gate script included, and excludes only itself. Its own SHA-256 is supplied to you OUT OF
  BAND in the handoff prompt and appears in no file here, so nothing in the room authenticates
  itself.

REQUIRED means: inside the transitive import closure of the S1b entry points, currently 15 modules
of the 50 supplied. Three modules outside that closure run a self-test battery and exit at import;
they are classified as SCRIPTS, not failures. **If a module INSIDE the closure ever does this, it
is an import defect and the launch is blocked.**

Verify BOTH the contract hash and the manifest against the value in your handoff prompt before you
write code, and re-verify every pinned file at the end so a mid-run edit cannot pass unnoticed.

Record the interpreter and the third-party package environment: `python`, `numpy`, `scipy`,
`mpmath` and anything else you use, with versions. Then affirm in writing that **no OpenWave or
project code, and no module in the required closure, was imported from outside this room.** The
standard library and third-party packages are necessarily external and are covered by the
environment record, not by this affirmation.

At handoff the gate reports: 10 entry points, a 15-module required closure, 50 modules supplied,
0 required-module failures, 35 modules outside the closure listed but not imported, contract hash
OK, and a manifest over every file in the room.
**The gate has been mutation-armed four ways and one arm found the gate itself vacuous.** Removing
a declared entry point from the room while the live checkout sat on `PYTHONPATH` originally
reported `ROOM READY`, because the module dropped out of the room scan and was therefore never
probed. Declared entry points are now probed BY NAME whether or not the scan finds them, and that
case reports `NOT SUPPLIED BY THE ROOM, resolved from <path>`. Treat the gate as fallible and say
so if you find another way past it.

**Q1. Implement every gate exactly as written.** G-REAL, G-RANK, G-SUBSPACE, G-ALIGN, G-SAMPLE,
G-BASIS, G-DISCRIM, G-WIRE. Each must be runnable, each must print its criterion and its measured
value, and each must be attached to the object the contract names, not a downstream proxy.

**Q2. Implement every mutation the contract specifies, and prove each has power.** For every arm:
show the parent is GREEN first, then apply the mutation, then show it goes RED, and record both. An
arm whose parent was never shown green proves nothing. This program has produced six vacuous arms;
assume yours is the seventh until you have shown otherwise.

**Q3. Adversarial vacuity audit, and this is the part that matters most.** For each gate, ask: can
this pass for a reason other than the property it claims to test? Construct the case that would
expose it. Specific hazards the contract already records, and you should look for more:
`G-SAMPLE` tested after orthonormalization is guaranteed by construction; `K_0` on a real 1x1 block
is zero by algebra; a relative rank cutoff on a numerically-zero matrix reports full rank; a
same-boundary swap preserves counts; the `n = 2` level cannot separate the correct realization from
the no-transpose one.

**Q4. Contract contradictions.** Report any place where the document is ambiguous, self-inconsistent,
or cannot be implemented without a choice it does not make. Do NOT choose. Stop and report.

**Q5. A dry run on synthetic inputs only, and it must test PRECEDENCE, not just reachability.**
Exercise the full path with manufactured matrices whose answers are known. One specimen per branch
proves only that all seven labels are reachable; it does not prove the adjudicator was assembled in
the right ORDER. A perfectly coded set of predicates can still be evaluated in the wrong sequence,
and this is the last place that can hide. Construct COLLISION cases where two rules both have a
true condition and show the earlier one wins:

| collision | must resolve to |
| --- | --- |
| a gate fails while the synthetic target otherwise looks SPECTRAL | `S1b-DEFECT` |
| one target AMBIGUOUS while another is qualified PERSISTENT | `S1b-NO_LABEL` |
| the two constructions disagree on `J` while one otherwise looks SPECTRAL | `S1b-NO_LABEL` |
| the constructions disagree on resolved `K` while another condition looks ADJOINT | `S1b-NO_LABEL` |

Also cover each of the seven outcomes at least once, so reachability and ordering are both shown.
**Do not run the live target.**

## 4. What you must not do

No live `R_0` spectrum at `n = 12` or `n = 20`. No S1b result. No threshold changes, no new
constants, no "improvements". No S2. Nothing here reopens M8.4 P1A, whose outcome is final.

The seed ladder `60, 120, 180` is preregistered; you do not enter a rung. `theta_Q` has no gate by
deliberate design, recorded in the contract as such: do not add one.

## 5. Reference material, and what its status is

`reference/` holds the merged M8.9 S1 record (`s1_note.md`, `s1_run.py`, `m8_9_s1_decision_rule.md`)
and the merged M8.4 P1A closeout. S1's note is where the contract's inherited G-WIRE criteria come
from, `‖L Q_0‖ <= 1e-8` against a measured 5.7e-10. `regenerate_estimator_table.py` shows the
established idiom for reconstructing a result from shipped packages. These are CONTEXT. The frozen
contract is the only thing that governs.

## 6. Deliverable form

A qualification note carrying: the contract hash verification, the Q0 manifest and environment, each
gate with its criterion and measured value, each mutation with its green parent and red child, the
vacuity audit with the cases you constructed, any contract contradictions, and the synthetic dry-run
routing table. Plus the code. Nothing is committed to any repository by you; the author lands it.

**A stop with evidence is a better outcome than a worked-around gate.** If you find yourself
reimplementing something because the supplied version will not import or will not do what you need,
that is a room defect: stop and say so.
