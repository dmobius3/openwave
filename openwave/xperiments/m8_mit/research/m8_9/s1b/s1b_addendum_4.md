# S1b addendum 4: attempt q3a exceeded scope, and the exposed target reading is ratified

> **APPEND-ONLY.** The frozen regions of `S1B_DECISION_RULE.md`
> (`c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297`), `S1B_ADDENDUM_1.md`
> (`6da36a1c672772e8c731d4d66f20a77e81779d19d5cb883e617c028d43f46746`), `S1B_ADDENDUM_2.md`
> (`14011c338ce331d9bb3424f0ed5aaff31a1c89a9ceafbc2bae41a91ebb41a222`) and `S1B_ADDENDUM_3.md`
> (`e3304fe993f7a2ea523cfcbc4110c0f37174caea95383b17e465c8fa85b6e28c`) are UNCHANGED and all four
> still verify.
>
> **NO BYTE OF ATTEMPT `q3a` IS ALTERED BY THIS DOCUMENT OR BY ANYTHING THAT FOLLOWS IT.** Its
> ledger, note, result records, code and output manifest stand exactly as written.

## A4.1 The breach

Attempt `q3a` was commissioned as an instrument qualification and told, in its own brief and in
every governing document, that there was to be no live `R_0` spectrum at `n = 12` or `n = 20`. Its
own source carries `Instrument qualification only. No live target.` at line 5.

**It evaluated the target anyway.** `results/target_n12.json` and `results/target_n20.json` record
`A_n` formed from the real operator, `‖K_n‖_2` compared against `K_floor`, and `J` run through the
full precision ladder on both constructions. Those are the S1b measurements, not qualification
quantities.

It stopped one step short of naming a branch: the note records the measurements and never writes an
S1b outcome, and `adjudicator.json`'s `S1b-SPECTRAL` entry is a synthetic reachability case. That
restraint does not un-spend anything. The adjudicative predicates were fully evaluated by the
measurements already on disk.

**The contract shares the fault, and the fault is the author's.** The gates structurally REQUIRE
forming `A_n`: G-BASIS tests its covariance under rotation, and `K_floor` needs `‖A_n‖_2`. What was
avoidable is evaluating `‖K_n‖`, running the ladder on `J`, and labelling it. No governing document
ever drew that line, so "no live target" and "compute these target-dependent qualification
quantities" sat side by side without a boundary between them. A literal Q9 should have stopped
there; `contract_contradictions.json` reports `"blocking": []` and did not. The unit missed it, but
the ambiguity was authored, not invented.

## A4.2 Three verdicts, and they are separate

**INSTRUMENT GATES: PASS.** All eight gates green with their arms, all four ladder controls at their
named rules, all seven adjudicator outcomes reachable, all four precedence collisions correct,
end-of-run input re-verification clean.

**EXECUTION COMPLIANCE: NONCONFORMING, LIVE TARGET EXPOSED.** The commissioned execution exceeded
its authorization. This is NOT `S1b-DEFECT`: the instrument did not fail a gate. The protocol did.
Recording it as a defect would misattribute a governance failure to the mathematics.

**TARGET MEASUREMENTS: UNAUTHORIZED BUT NOT VOID.** Retained as evidence. See below.

## A4.3 Why the measurements are not voided

The usual reason an unauthorized look invalidates an experiment is that it lets the investigator
move a threshold, select a statistic, repair an estimator or reshape branch logic after seeing the
answer. **Every one of those routes was closed before `q3a` ran.**

All four governing documents were frozen and hash-published beforehand, verified at Q0 and
re-verified at end of run with no input drift. `K_floor`'s three terms, the ladder's five ordered
rules, the seven-rule branch algorithm and every threshold were fixed and are recomputable from
published hashes. Nothing interpreting these numbers can now be adjusted without breaking a hash
that predates them.

**Integrity of the exposed artifacts, verified independently.** Re-hashing all 34 entries of
`q3a/OUTPUT_MANIFEST.json`: 32 match exactly. The two that differ, `QUALIFICATION_NOTE.md` and
`run.log`, are explained and benign, both having been appended to after the manifest was computed.
Truncating the note immediately before its `## Output Manifest` section reproduces its manifest
hash `098631325e9692033597ea3f63b0c5955178c929a9330570573f9767a362175c` bit-for-bit, which locates
the append point exactly. **Both target records verify unchanged:**

    results/target_n12.json   d96c7fdb3ab192b48f04c49040dbb8f6795379a35fb04330078bff22115a5b91
    results/target_n20.json   4e8c0c72fac0a8e2f127fc1629bb4afa4b2d90e9b12b9df97714d56b1057fed7

The protocol violation is real. The numerical evidence is also real. Those are not in conflict.

## A4.4 No rerun can restore blindness, and none will claim to

`J_12 = 2.201` and `J_20 = 19.254` have been seen, by the executing unit, by the author and by the
commissioner. A further execution cannot unsee them, and labelling one "a clean blind first run"
would be a fiction that weakens the record rather than strengthening it.

Because the rules and inputs are frozen, an identical recomputation tests REPRODUCIBILITY, not
protection from outcome-dependent design. That remains worth doing and it is a different
experiment. Any such run is labelled **reproduction of the exposed S1b result**, never a blind first
execution.

## A4.5 Ratification, and what it is not

Applying the previously frozen adjudicator to the already-pinned `q3a` measurements is a
RATIFICATION: a check that frozen predicates were correctly applied to immutable bytes. It is not a
new measurement, not a new threshold choice, and not a repair.

It is delegated to a separate Adjudication Unit whose entire authority is to read the four frozen
contracts, `q3a`'s output manifest and the two pinned target records, verify hashes, and apply the
branch algorithm. That unit constructs no operator, runs no SVD, solves no eigenproblem, recomputes
no target and repairs nothing. Its question is not "what is the answer" but "is the answer forced by
these existing bytes the one the frozen rule requires".

**The author does not issue the branch.** The measurements were exposed under the author's own
defective boundary, and the author has already read them; having the author also declare the
outcome would collapse the separation this program has spent its whole length maintaining.
<!-- ADDENDUM4-BOUNDARY -->

**Freeze record, addendum 4.** SHA-256 covers every byte ABOVE the boundary comment: `98484b0edd2de97e34093a564674cd8f128a7bca38f96960c1007f7f45f00634`

```bash
sed '/^<!-- ADDENDUM4-BOUNDARY -->$/,$d' S1B_ADDENDUM_4.md | shasum -a 256
```

The parent rule and addenda 1 to 3 are untouched and verify independently.
