# T9: The scope line: what a result licenses and what it does not

> Platform task **T9**. Status: ✅ **CLOSED 2026-10-27** (go 11:29 EDT, batched with T10) · Roadmap:
> [`../platform_roadmap.md`](../platform_roadmap.md) · Deliverable lands in
> [`METHOD_NOTE.md`](../METHOD_NOTE.md) (the model-owner-facing standard) and
> [`PR_REVIEW_STANDARDS.md`](../PR_REVIEW_STANDARDS.md) (the review checklist).

## PLANNING

**Scope.** One mandatory line in every results post, task review and method note: **Scope**, stating what the number or verdict licenses and what it does not (the object it was computed on, the box and spacing or the sample, the claims it does not support: no particle, no force, no physical unit unless anchored). The same line is checked at PR review. The idea comes from the claim registry of the repository that hosts the M5 discussion, where every accepted claim carries a `does_not_license` field and a `maximum_verdict` fixed before the run; the M5.32 rungs already state this in prose (the "what the poll does not do" and "not in R27" paragraphs), and the line makes it a fixed slot a reader can find.

**Definition of done.** The `METHOD_NOTE.md` template gains the slot (one line, with two filled examples from existing records: a reproduction verdict and a lattice gap); `PR_REVIEW_STANDARDS.md` gains the check ("the scope line is present and matches the result"); the M5 roadmap's task-review template points at it; one existing record (the M5.35 review) gets the line retrofitted as the worked example. No change to any result.

**Gating.** None beyond the user's go; a doc change in two standards and one worked example.

**Blindspot pass.** (i) A scope line that repeats the title is noise: the template says the line names the exclusions, not the inclusions. (ii) The line must not become a disclaimer wall: one sentence, at most two. (iii) It is not a status label: the ✅ / 🔶 / ⚠️ labels stay as they are.

## FINDINGS (2026-10-27)

| Deliverable | Where | State |
| --- | --- | --- |
| The slot | [`METHOD_NOTE.md`](../METHOD_NOTE.md): row 7 of the required shape, a pre-send check, and the section "The scope line" with the two filled examples (the M5.35 reproduction verdict; the R29 lattice-gap read) | ✅ |
| The review check | [`PR_REVIEW_STANDARDS.md`](../PR_REVIEW_STANDARDS.md) Gate C, check C9: a missing line is a requested item, a line that licenses more than the artifact earns is a C7 finding | ✅ |
| The task-review pointer | the M5 roadmap's CONVENTIONS, one entry pointing at the section (the review template itself lives outside this repo, so the convention is the pointer) | ✅ |
| The worked example | the M5.35 record's TASK REVIEW carries the line, retrofitted and marked as such | ✅ |

**Scope**: a doc standard and one retrofit; it licenses the slot's shape and where it is checked, and does not change any result or any cell.
