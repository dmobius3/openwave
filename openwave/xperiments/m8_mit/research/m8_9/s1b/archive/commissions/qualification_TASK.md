# M8.9 S1b: full implementation qualification

> **This supersedes both earlier commissions.** Round 1's is retained as
> `prior/TASK_ROUND1.md`. Round 1 finished without writing a note. Round 2 correctly STOPPED on a
> room defect and then lost its connection. Neither left usable evidence, so this is a FULL
> qualification, not a narrow rerun.
>
> The round-2 narrow commission's bytes were NOT preserved: it was overwritten in place by this
> file. That is the same failure addendum 2 § A2.4 exists to prevent, committed one more time
> while staging the fix for it, and it is recorded here rather than hidden. Nothing in this round
> depends on those bytes.
>
> Read the contract documents in this order, later governing earlier:
> `contract/S1B_ADDENDUM_3.md`, `contract/S1B_ADDENDUM_2.md`, `contract/S1B_ADDENDUM_1.md`, then
> `contract/S1B_DECISION_RULE.md`.
>
> Instrument qualification only. No live `R_0` spectrum at `n = 12` or `n = 20`. No threshold
> changes. If the contract cannot be implemented as written, STOP and report rather than choosing.

## 1. Provenance, before any other work

| document | frozen-region SHA-256 | boundary |
| --- | --- | --- |
| `contract/S1B_DECISION_RULE.md` | `c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297` | `<!-- FREEZE-BOUNDARY -->` |
| `contract/S1B_ADDENDUM_1.md` | `6da36a1c672772e8c731d4d66f20a77e81779d19d5cb883e617c028d43f46746` | `<!-- ADDENDUM-BOUNDARY -->` |
| `contract/S1B_ADDENDUM_2.md` | `14011c338ce331d9bb3424f0ed5aaff31a1c89a9ceafbc2bae41a91ebb41a222` | `<!-- ADDENDUM2-BOUNDARY -->` |
| `contract/S1B_ADDENDUM_3.md` | `e3304fe993f7a2ea523cfcbc4110c0f37174caea95383b17e465c8fa85b6e28c` | `<!-- ADDENDUM3-BOUNDARY -->` |

Run `room_import_gate.py` first; it is a launch precondition and checks all three plus the room
manifest. Verify the manifest against the value in your handoff prompt, which appears in no file
here. Record the interpreter and third-party versions. Affirm that no OpenWave or project code, and
no required-closure module, was imported from outside this room; the standard library and
third-party packages are external by necessity and are covered by the environment record.

## 2. THE FIRST THING YOU BUILD IS THE OUTPUT LEDGER

Addendum 2 § A2.4 requires it and addendum 3 § A3.2 fixes its paths, which are keyed by ATTEMPT and
not by round. Your `<attempt_id>` is supplied in your handoff prompt. Do not choose, derive or alter
it.

**Collision guard, first.** If `qualification_runs/<attempt_id>/` already exists, STOP and report,
without reading, writing, appending to or deleting anything inside it. Then, the moment Q0
provenance succeeds and before any other work, create

    qualification_runs/<attempt_id>/QUALIFICATION_NOTE.md   write Q0 into it NOW, append as you go
    qualification_runs/<attempt_id>/results/<gate>.json     one record per gate, on completion
    qualification_runs/<attempt_id>/OUTPUT_MANIFEST.json    written last

**Append after every completed gate, control and mutation. Do not batch.** Two runs have already
been lost because their evidence lived only in a terminal. A run that dies partway must leave
everything it established on disk and hashable.

You write ONLY inside your own attempt directory. Earlier attempts' ledgers are never read, written
or modified by you. The room's INPUT manifest is never modified by you; the OUTPUT manifest is
yours.

## 3. `prior/s1b_qualification.py` is ARCHIVAL ONLY

Verify its pinned hash `5a9e04845375c4d12c3a475607f20a1d5f13cc82829d64caa022c82d5e784802` and then
leave it alone. Per addendum 2 § A2.3 you may NOT execute it, relocate it, patch it, run it with an
altered `__file__`, bypass its manifest assertion, or use any output of it for qualification credit.
It cannot run in this room by construction and that is a recorded fact, not an obstacle to route
around. **No round-1 result receives evidentiary credit**, including results you may see quoted
anywhere.

## 4. What to qualify, in full

**Q1. All eight gates**, each runnable, each printing its criterion and measured value, each
attached to the object the contract names rather than a downstream proxy: G-REAL, G-RANK,
G-SUBSPACE, G-ALIGN, G-SAMPLE, G-BASIS, G-DISCRIM, G-WIRE.

**Q2. Every mutation, each from a demonstrated GREEN parent.** Show the parent green FIRST, apply
the mutation, show it red, record both. Includes the bridge arms of addendum 1 and the G-SUBSPACE
tilt arms under the sine form.

**Q3. All four ladder controls** of addendum 1 § A1.1, each reported with the RULE it reaches:
`k = 2` via rule 1, `k = 3` via rule 2, the synthetic triplet `(1e-8, 1e-12, 1e-12)` via rule 4,
and `(1e-8, 2e-11, 2e-11)` via rule 5, plus the two rule-4 tail edges.

**Q4. Sine-form angle machinery throughout.** `s = ‖(I - P_a) Q_b‖_2`, half-angle
`d = 2 sin(0.5 arcsin(clip(s, 0, 1)))`. The arccos route and the
`sqrt(2(1 - sqrt(1 - s^2)))` identity are BOTH PROHIBITED; addendum 1 § A1.4 says why each fails and
in which direction. Measure `theta_C` at both targets and `theta_Q`, which has never been measured
under the sine form.

**Q5. The shipped SVD feeds the qualifying path directly.** Run `invariant_dim_and_basis` ONCE per
target, serialize its returned dimension, basis and diagnostics immediately as qualification output,
hash that artifact, and have every downstream gate consume that exact serialization. Expect roughly
805 s at `n = 12` and 3919 s at `n = 20`; the `n = 20` call allocates about 44.8 GB nominally
against this machine's 25.8 GB and has been observed to complete anyway, so give it room and do not
kill it early. **Serialize before running anything downstream**, so the expensive call is never at
risk of being lost.

Record the bitwise `s` and `Vh` comparison against the economy variant as an equivalence
DIAGNOSTIC. Per addendum 2 it licenses nothing.

**Q6. `K_floor`** recomputed with measured `max_r kappa(W_n^r)`, never the `1e6` ceiling, and the
discrepancy term from the sine-form `theta_Q`. **Report which of the three terms dominates. The
addendum deliberately refuses to prejudge this**, so it is a result of your run, not a check against
an expected answer.

**Q7. All seven adjudicator outcomes reachable, AND the four precedence collisions**, which prove
ordering rather than reachability: a gate failing while the target otherwise looks SPECTRAL must
give DEFECT; one target AMBIGUOUS while another is qualified PERSISTENT must give NO_LABEL; the two
constructions disagreeing on `J` while one otherwise looks SPECTRAL must give NO_LABEL; the
constructions disagreeing on resolved `K` while another condition looks ADJOINT must give NO_LABEL.

**Q8. Adversarial vacuity audit.** For each gate ask whether it can pass for a reason other than the
property it claims, and construct the case that would expose it.

**Q9. Contract contradictions.** Report any place the documents are ambiguous, self-inconsistent, or
unimplementable without a choice they do not make. Do NOT choose. Stop and report.

## 5. What NOT to do

No live target. No threshold changes, no new constants, no improvements. No S2. No seed rung:
`60, 120, 180` is preregistered and entered only on a G-SAMPLE failure. No gate on `theta_Q`;
addendum 1 records its absence as deliberate. Nothing reopens M8.4 P1A, whose outcome is final.

## 6. Standing hazards, all already paid for

This program has produced eight vacuous mutation arms and three numerical-formula defects, two of
the latter found only by RUNNING a published check rather than reading it. Assume your arm is the
ninth until you have shown a green parent and a red child, with both recorded.

Named traps: a control green by algebra rather than by the property it claims, such as `K_0` on a
real 1x1 block or `k = 2` reaching rule 1 rather than rule 4; a relative cutoff applied to a
numerically zero matrix; a specimen sitting on a comparison boundary so its verdict turns on
rounding; an angle formula that manufactures 1e-08 structure from nothing; a half-angle identity
that annihilates 1e-10 structure into exact zero; a conditioning test applied after the step that
guarantees it; a gate that reports success because a module called `sys.exit(0)` during import.

## 7. End-of-run input re-verification, then the deliverable

**Immediately before writing `OUTPUT_MANIFEST.json`**, re-hash every file listed in
`ROOM_MANIFEST.json` and re-verify all four frozen-document hashes. **Any change to a
manifest-pinned input since Q0 is a STOP**, and you leave the partial ledger in place as evidence.
This is an INPUTS check, not a rerun of the import gate: by then the room legitimately contains
generated output the input manifest does not list, so the gate's accounting is no longer the right
instrument. Several failures in this program have involved a room changing underneath its own
report.

The output ledger of § 2, complete, plus your code. A stop with evidence beats a worked-around gate.
If a supplied component will not do what the contract needs, that is a room defect: stop and say so
rather than substituting. Round 1 substituted and disclosed; round 2 stopped correctly. This round
stops too.
