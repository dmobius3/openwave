# T7: Review-cycle overhead, one document one merge

> Platform task **T7**. Status: ✅ **CLOSED 2026-08-01** · Roadmap:
> [`../platform_roadmap.md`](../platform_roadmap.md) · Rule landed in
> [`CONTRIBUTING.md § How review works`](../../CONTRIBUTING.md#how-review-works).

Filed and closed the same day: the deliverable is a rule, and it landed with this record.

## PLANNING

The question was whether the M8 column's process weight had become overhead on the model
author and the maintainer, given that it visibly carries more review activity than M5. The
task measures the two columns rather than arguing from impression, and changes only what the
measurement implicates.

[T5](t5_task_details.md) was filed on this same pain, "the coordination edits alone took
three reviewed merges (#382, #383, #384); for a fork-based author each such edit is a full
fork → PR → DCO cycle, and the cost is linear in edits × authors". T5 evaluated moving
coordination state to GitHub Issues and **declined**, so the pain it named was left standing.
T7 is the answer to it.

## FINDINGS

### The premise is wrong: documentation volume is not the difference

| Column | Markdown | Python | Ratio |
| --- | --- | --- | --- |
| M5 | 31,910 lines | 97,578 lines | 1 : 3.06 |
| M8 | 4,686 lines | 15,205 lines | 1 : 3.24 |

M8 is marginally more code-heavy per documentation line than M5. Whatever is heavier about
M8, it is not the amount written.

### The difference is round trips, and it comes from having two parties

| Pull requests, the two weeks to 2026-08-01 | M5 | M8 |
| --- | --- | --- |
| Result drops (a run closed: data, plots, findings; +3k to +76k lines) | ~20 | 5 |
| Sub-200-line document sync, clarify or amend PRs | ~4 | **12 of 26** |
| Authors involved | 1 | 2 |

Two chains show the shape:

```text
#374 author amends condition 3  →  #375 author syncs the docs  →
#376 maintainer clarifies the gate  →  #378 author "recovering orphaned fixes from #374 and #375"

#382 freeze  →  #383 clarify  →  #384 adopt the author's three directions from #382  →
#385 author addendum  →  #386 record the addendum as landed
```

Four merges for one amendment, the last existing only because fixes were orphaned by the
first two; five merges to agree on one document. M5 does not do this because M5 has one
author: alignment there costs an edit, not a merge cycle.

### Where the load actually falls

Nineteen of the 26 M8 pull requests are the maintainer's, and the small ones are the ones
nobody chose. The author's seven are large, self-contained, and filed on the author's own
schedule; every additional constraint in them (the context firewall, the claim ceiling, the
two-packet split) was volunteered rather than requested. The overhead risk is on the
maintainer side, not the author side, which is the opposite of what the question assumed.

### The habit, and why it is cheap to change

Making the agreement in the merge history instead of in the PR thread: seeing a needed change
during review and opening a PR for it afterwards, rather than settling it in the thread or
applying it at merge. Changing it moves **when** an edit lands, never **whether** it is
checked, so no standard has to move with it.

The precedent already existed and was used once without being written down, at the M8.2
close-out: "maintainer edits were applied at merge and announced in the thread rather than
sent back as a review round."

### What changed

The rule has two audiences and therefore two homes, one canonical for each side.

| File | Change |
| --- | --- |
| [`PR_REVIEW_STANDARDS.md § 10.3`](../PR_REVIEW_STANDARDS.md#103-one-document-one-merge) | **Canonical for the reviewer.** New subsection under Maintainer edits: the two PR chains as evidence, rules R1-R4, the maintainer's own coordination PRs named as the same failure, and an explicit statement that no gate is traded away to save a cycle. Plus a lessons-log row keyed to [#378](https://github.com/openwave-labs/openwave/pull/378) |
| [`CONTRIBUTING.md § How review works`](../../CONTRIBUTING.md#how-review-works) | **Canonical for the contributor.** What to expect: one pass of findings, maintainer fixes applied at merge and announced, frozen documents by batched addendum. Asks contributors to enable "Allow edits by maintainers", which is what makes the merge-time fix possible at all |
| [`../platform_roadmap.md`](../platform_roadmap.md) | `HOW REVIEW WORKS` convention pointing at both; T7 row; change-log entry |
| [`../../ONBOARDING_MODELS.md`](../../ONBOARDING_MODELS.md) | Pointer from the light-PR-review paragraph in STEP 3, the last thing a new model author reads before a first PR |

**Why the reviewer's copy is the load-bearing one.** The habit was executed by AI coding agents acting as reviewers, not by a human skipping a step. An agent's review procedure is [`PR_REVIEW_STANDARDS.md`](../PR_REVIEW_STANDARDS.md), reachable as entry 12 of the [`CLAUDE.md`](../../CLAUDE.md) document map, and a rule that lives only in the contributor-facing document would never be read by the party performing it.

### The rule was discoverable and still went unread

Raised by the user during this task: the review standards were not in the working context when [#402](https://github.com/openwave-labs/openwave/pull/402) was fetched and characterized, and the habit being fixed was performed by reviewing agents in the first place.

The document was never hidden. It sits at entry 12 of the [`CLAUDE.md`](../../CLAUDE.md) map, and [`CONTRIBUTING.md`](../../CONTRIBUTING.md) points at it. That is the defect: **a map states what exists, while a mandatory rule fires when a condition is met**, and only the second survives an agent that loads the map, takes the rows its task seems to need, and starts reading a diff. The same file already carried three rules in the second form (AI hygiene, the adversarial audit, the method note), so the fix was to move pull-request review into that block rather than to invent a mechanism.

| Change to [`CLAUDE.md`](../../CLAUDE.md) | Effect |
| --- | --- |
| New `### Pull-request review (any PR in play): MANDATORY`, beside the other three standing rules | Fires on a condition instead of waiting to be selected from a list |
| The trigger stated wider than the word "review": summarizing or characterizing a pull request counts, as does any read that could inform a merge, approval or comment | Closes the exact gap this task hit. A characterization without the gates reaches the same merge decision as a review |
| The section names what the document holds (intake, tiers, gates A-G, the commitment sweep, § 10, the verdict-is-a-review rule) | An agent that skips it can at least see what it is skipping |
| Doc-map row 12 and the Code Style table both marked **MANDATORY**, read before the diff | Three surfaces agree, so no single one has to be the lucky hit |

**Already present, and left alone.** § 10 and fairness rule G3 already required mechanical fixes to be made at merge rather than sent back, written after [#340](https://github.com/openwave-labs/openwave/pull/340) taught the same lesson inside a single review. What was missing was the cross-PR case: nothing said that a point visible during review must not become its own pull request afterwards, nothing addressed batching addenda to a frozen document, and nothing covered the maintainer's own coordination PRs. § 10.3 adds those three and changes none of the existing rules.

### Two stale references found and fixed

[T5](t5_task_details.md) closed two days earlier and made a task a roadmap row rather than a GitHub issue, but the review standards were not swept for it. Both surviving references told a reviewer to open issues: the 🔶 `Approve with follow-ups` verdict row, and the author-gated fallback in § 8.1. Corrected to roadmap rows, citing T5. Found only because this task read the document end to end, which is its own small argument for sweeping a rule's consumers when it lands.

### What deliberately did NOT change

Nothing in the verification layer. Blind protocols, adversarial audits, mutation tests,
pre-registration, freezes, claim ceilings and evidence classes all stand exactly as written,
because their measured catch rate is what justifies them: four checks that could not fail and
one mutation running green with the wrong Dynkin diagram installed (M8.1.1), an algebraically
forced identity mistaken for an empirical law (M8.1.1, a platform-side overclaim), a defect in
a published page (M8.3), a circular target that would have been pre-registered as cross-model
evidence (M8.6), and an `Ad(v⁻¹)` inverse no character check can detect
([PR #402](https://github.com/openwave-labs/openwave/pull/402) § 2). Three of those were
platform-side errors, which is the tell that the apparatus is load-bearing rather than
defensive.

Pre-registration size is not treated as overhead either. A frozen contract exists to end
renegotiation, so it removes round trips rather than adding them.

**Not propagated to other repositories.** The rule addresses two-party review over a fork.
Where a repository has a single author and no external pull requests, it has nothing to
batch.

## Cross-links

| Direction | Doc |
| --- | --- |
| The rule | [`CONTRIBUTING.md § How review works`](../../CONTRIBUTING.md#how-review-works) |
| The pain, named and left standing | [T5](t5_task_details.md) |
| The unwritten precedent | [`m8_roadmap.md § CHANGE-LOG`](../../openwave/xperiments/m8_mit/research/m8_roadmap.md), M8.2 close-out |
| First PR reviewed under the rule | [#402](https://github.com/openwave-labs/openwave/pull/402), M8.5-B pre-registration |
