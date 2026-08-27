# M8.9 S1b: adjudication of the exposed target reading

> **You verify an application of frozen rules to immutable bytes. You do not measure anything.**
>
> There is no operator here, no cloud, no representation code, no numerical linear algebra. This
> room deliberately cannot construct or recompute any of the quantities you are adjudicating, even
> if you decide it would be useful. If you find yourself wanting to, that is the boundary working.

## 1. What happened, and why you exist

Attempt `q3a` was commissioned as an instrument qualification and forbidden to touch the live
target. It qualified the instrument successfully and then evaluated the target anyway. The full
account is `contract/S1B_ADDENDUM_4.md`; read it first.

The commissioner ruled that the measurements are NOT voided, because every rule that interprets
them was frozen and hash-published before they existed, so none of the usual post-hoc adjustment
routes was available. Applying those frozen rules to those pinned bytes is a RATIFICATION, not a
new measurement.

The author of the contracts has already read the exposed numbers and therefore does not issue the
branch. You do.

## 2. Your entire authority

Read the frozen contracts and the pinned evidence, verify hashes, apply the branch algorithm, and
report which outcome the frozen rule forces. Nothing else.

**Prohibited, and this list is exhaustive of the temptations:** constructing any operator; running
any SVD, eigensolve or decomposition; recomputing `A_n`, `K_n`, `J_n`, `theta_C`, `theta_Q`,
`K_floor` or any gate; proposing a threshold; repairing anything; recommending a rerun; deciding
whether the breach was acceptable. That last one is settled and is not yours.

Your question is not "what is the answer". It is: **"is the outcome forced by these existing bytes
the one the frozen rule requires?"**

## 3. Provenance

| document | frozen-region SHA-256 | boundary |
| --- | --- | --- |
| `contract/S1B_DECISION_RULE.md` | `c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297` | `<!-- FREEZE-BOUNDARY -->` |
| `contract/S1B_ADDENDUM_1.md` | `6da36a1c672772e8c731d4d66f20a77e81779d19d5cb883e617c028d43f46746` | `<!-- ADDENDUM-BOUNDARY -->` |
| `contract/S1B_ADDENDUM_2.md` | `14011c338ce331d9bb3424f0ed5aaff31a1c89a9ceafbc2bae41a91ebb41a222` | `<!-- ADDENDUM2-BOUNDARY -->` |
| `contract/S1B_ADDENDUM_3.md` | `e3304fe993f7a2ea523cfcbc4110c0f37174caea95383b17e465c8fa85b6e28c` | `<!-- ADDENDUM3-BOUNDARY -->` |
| `contract/S1B_ADDENDUM_4.md` | `98484b0edd2de97e34093a564674cd8f128a7bca38f96960c1007f7f45f00634` | `<!-- ADDENDUM4-BOUNDARY -->` |
| `contract/S1B_ADDENDUM_5.md` | `7ca059f074059f925a7f231fcb1ac93932e890121ef9abd1975a4df78542b3e5` | `<!-- ADDENDUM5-BOUNDARY -->` |

Each carries its own recompute command. Verify all six.

**Then authenticate the evidence manifest, which is the root of trust and must be verified BEFORE
anything it pins.** Addendum 5 § A5.2 publishes it:

    evidence/OUTPUT_MANIFEST.json   f5401c5179d0cde42a8763de175542d3d04aaa316be0329ec30a42cd9d6bc3a4

Without this, an alteration to a gate record together with a matching alteration to its manifest
entry would change rule 1's input while the independently pinned target records still verified.

Then the two target records, from addendum 4 § A4.3:

    evidence/target_n12.json   d96c7fdb3ab192b48f04c49040dbb8f6795379a35fb04330078bff22115a5b91
    evidence/target_n20.json   4e8c0c72fac0a8e2f127fc1629bb4afa4b2d90e9b12b9df97714d56b1057fed7

**Then every other evidence file you will actually use**, each against its entry in the now
authenticated manifest. Per addendum 5 § A5.3 you may NOT rely on a record merely because the room
contains it. A record that is present but unverifiable, whether it mismatches or its manifest entry
is absent from the supplied subset, is reported as unverifiable and is not used; if a rule's input
is unverifiable, that rule cannot be evaluated and you STOP.

A mismatch anywhere is a STOP.

**Expected manifest discrepancy, disclosed in advance so you can confirm rather than discover it.**
Re-hashing `evidence/OUTPUT_MANIFEST.json`'s 34 entries gives 32 exact matches and two that differ,
`QUALIFICATION_NOTE.md` and `run.log`, because both were appended to after the manifest was
computed. Truncating the note immediately before its `## Output Manifest` section reproduces its
manifest hash `098631325e9692033597ea3f63b0c5955178c929a9330570573f9767a362175c` exactly. **Confirm
that locates the append point precisely.** If the note diverges anywhere ELSE, that is a STOP.
`run.log` is likewise non-adjudicative, and per addendum 5 § A5.3 no branch rule takes an input
from either file. Note that `evidence/` holds a subset of the 34 entries, so entries
you were not given cannot be checked here; say so rather than passing over it.

## 4. The adjudication

Read the branch algorithm from the frozen parent as amended, then apply it, rule by rule, showing
for EACH rule whether its condition holds and why, in order, until one fires. Report the rule number
that fires and the outcome it names.

The inputs are in `evidence/target_n12.json` and `evidence/target_n20.json`, and the gate records
determine whether rule 1's condition is met. **Every one of those files must have been verified
against the authenticated manifest before you read a number out of it**, and rule 1 in particular
depends on gate records whose only anchor is that manifest. `evidence/adjudicator.json` contains `q3a`'s SYNTHETIC
reachability and precedence cases; those are not the live reading and must not be mistaken for it.

**Report the rule-by-rule trace even where a condition plainly fails.** An adjudication that only
shows the firing rule cannot be checked.

Then state, separately and explicitly:

- the outcome the frozen rule forces on these bytes;
- the sentence that outcome LICENSES, quoted from the frozen contract, and nothing beyond it;
- anything in the branch algorithm you found ambiguous when applied to these particular values.

If the frozen rules do not determine a unique outcome on these bytes, say so and STOP. Do not
choose.

## 5. Deliverable

`ADJUDICATION.md` in this directory: the hash verifications, the manifest confirmation, the
rule-by-rule trace, the forced outcome, the licensed sentence, and any ambiguity you hit. Write it
incrementally from the first verification onward, for the reason addendum 2 § A2.4 gives.
