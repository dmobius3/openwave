# Pull Request Review Standards

> The procedure a maintainer (or an AI coding agent acting for one) follows when reviewing a pull request, especially one from an external contributor. It exists so that review quality does not depend on who is on duty, and so that the repository's DNA survives contact with contributions it did not plan for.
>
> **This is a live file.** Every PR that teaches us something new gets folded back into it: a new gate, a sharper command, or a row in the [lessons log](#14-lessons-log). If a review turns up a failure mode not covered here, add it before closing the PR.

## What review is, and is not

[`MODELS.md`](../MODELS.md) promises contributors a **light PR review, focused on reproducibility and honest documentation, not ideological gatekeeping**. That promise is binding. This document does not raise the bar on *which physics is allowed*; it makes concrete what "reproducible and honestly documented" means, and it adds the checks that protect the shared surfaces every other contributor depends on.

| Review IS | Review is NOT |
| --- | --- |
| Does the claim trace to something runnable, and does the artifact actually say what the text says | Does the maintainer agree with the framework |
| Does the change stay inside its own model's blast radius | Does the code match the maintainer's personal style |
| Does a shared file, another author's column, or a front-door doc change safely | A demand that the contributor fix pre-existing repository debt |
| Is the claim's strength proportional to the evidence | A demand that the result be positive |

A documented negative is a merge-worthy contribution. An overstated positive is not.

## Contents

| Section | What it covers |
| --- | --- |
| [1. Intake](#1-intake) | the five things to establish before reading a line of code |
| [2. Blast-radius map](#2-blast-radius-map) | which paths get which level of scrutiny |
| [3. Gate A: safety and hygiene](#3-gate-a-safety-and-hygiene) | large files, encoding, copyright, secrets, dangling references |
| [4. Gate B: scope containment](#4-gate-b-scope-containment) | model-folder discipline, shared files, root documents |
| [5. Gate C: claim to artifact](#5-gate-c-claim-to-artifact) | recompute the headline number from the shipped data yourself |
| [6. Gate D: the adversarial pass](#6-gate-d-the-adversarial-pass) | is the mechanism wired, is the signal above the noise, do the knobs manufacture the result |
| [7. Gate E: MODELS.md cell changes](#7-gate-e-modelsmd-cell-changes) | the evidence bar for moving a cell, and the linter a maintainer runs (§ 7.1) |
| [8. Gate F: other authors' work](#8-gate-f-other-authors-work) | not damaging a column you do not own |
| [9. Gate G: policy sweep](#9-gate-g-policy-sweep) | AI hygiene, conduct, contributing, reproduce, onboarding, style |
| [10. Maintainer edits](#10-maintainer-edits) | when to fix it yourself instead of asking, and how to push to a fork |
| [11. Fairness rules](#11-fairness-rules) | measure against the enforced baseline, not the aspirational one |
| [12. Verdict and how to write it](#12-verdict-and-how-to-write-it) | the ladder, and open-source review etiquette |
| [13. Command appendix](#13-command-appendix) | copy-paste checks |
| [14. Lessons log](#14-lessons-log) | what past PRs taught us |

## 1. Intake

Establish these five before reading any code. They set how heavy the rest of the review needs to be.

| # | Establish | How | Why it matters |
| --- | --- | --- | --- |
| 1 | **Who** and **prior history** | `gh pr view <N> --json author` then `gh pr list --state all --author <login>` | A returning contributor's past review threads tell you what they already know and what they had to be asked twice |
| 2 | **DCO sign-off on every commit** | `gh pr checks <N>`, plus `git log main..pr-<N> --format='%(trailers:key=Signed-off-by,valueonly)'` | Apache 2.0 provenance. Non-negotiable, and it is a one-time config fix, not a rejection (see [lessons log](#14-lessons-log)) |
| 3 | **Blast radius** | `gh pr view <N> --json files` | Drives everything below. A PR entirely inside one model folder is a different review from one touching `common/` or a root document |
| 4 | **What the PR claims** | the PR body plus any added research note | Write the headline claim down in one sentence. Gates C and D are tested against that sentence |
| 5 | **Source branch** | `headRefName` | A PR from a fork's `main` is workable but fragile: the contributor cannot start a second PR without entangling it. Worth a friendly note, never a blocker |

## 2. Blast-radius map

Scrutiny scales with how many people a file can break. Classify every changed path.

| Tier | Paths | Standard |
| --- | --- | --- |
| **T1 model-local** | `openwave/xperiments/<model>/` (except the briefing) | Light review. The column's author owns the physics. Check reproducibility, honesty, and that it does not reach outside |
| **T2 model front door** | `__M<x>_model_briefing.md`, that model's `research/` roadmap and trackers | The column's public face. Check that status language matches what the artifacts support |
| **T3 shared code** | `openwave/common/`, `openwave/gui/`, anything imported by more than one model | Heavy review. A defect here is every column's defect. Trace every caller before approving |
| **T4 shared surfaces** | [`MODELS.md`](../MODELS.md), [`README.md`](../README.md), [`CLAUDE.md`](../CLAUDE.md), [`AI_HYGIENE.md`](../AI_HYGIENE.md), [`CONTRIBUTING.md`](../CONTRIBUTING.md), [`REPRODUCE.md`](../REPRODUCE.md), [`ONBOARDING_MODELS.md`](../ONBOARDING_MODELS.md), [`CODE_OF_CONDUCT.md`](../CODE_OF_CONDUCT.md), `SECURITY.md`, `TRADEMARK.md`, `LICENSE`, `NOTICE`, `.github/` | Maintainer-owned. Default answer to an external edit here is "let us handle it in a separate maintainer PR", unless the change is a plain factual correction |
| **T5 another model's folder** | any `xperiments/<other model>/` | See [Gate F](#8-gate-f-other-authors-work). Requires the other column's author in the loop |

**Print the tier split explicitly in the review.** A contributor should never discover after the fact that one line in their diff was the part that needed care.

## 3. Gate A: safety and hygiene

Mechanical checks. Run them all; they take under a minute and they catch the things a physics read will not.

| # | Check | Bar | Command |
| --- | --- | --- | --- |
| A1 | **No large data files** | Tracked data is summary-scale: JSON, CSV, plots, manifests. Heavy arrays (`.npz`, `.npy`, raw fields, videos) stay local and gitignored, with a `_DATASETS.md` manifest recording the regeneration script instead | `git ls-tree -r -l pr-<N>` sorted by size, and diff the total against `main` |
| A2 | **Valid UTF-8, no BOM** | Every text file decodes as UTF-8 and starts without a byte-order mark | [see appendix](#13-command-appendix). A single CP1252 byte silently breaks `grep`, `black`, and some editors |
| A3 | **No third-party copyrighted material** | Citations, DOIs, and links only. Papers and author documents live in the gitignored `theory/` folder and are recorded in that model's `theory/_CITATIONS.md`. Contributor's own work, DCO-signed, is fine and does not need to be hidden | Read every added `.md` and `.pdf`. Check the reference list resolves to citations, not reproduced text |
| A4 | **Citations registered** | New DOIs or author documents referenced by the PR appear in the model's `theory/_CITATIONS.md`, and stale "unpublished / DOI n/a" entries get updated | `git diff main...pr-<N> -- '*_CITATIONS.md'` |
| A5 | **No secrets or personal data** | No keys, tokens, absolute home paths, emails beyond the commit trailer | grep the diff for `key`, `token`, `secret`, `/Users/`, `C:\Users` |
| A6 | **Deleted files leave no dangling references, and no capability** | Every deletion is either unreferenced or every referencing import and link is repointed in the same PR. Then ask the second question: what did the deletion *remove*? A reference check is clean by construction when a file is self-contained, which is exactly the case for an xperiment configuration, a validation script, or a research note. Deleting one takes an entry off the menu and nothing complains | grep the deleted basenames across the PR branch, `.py` and `.md` both, then open each deleted file and ask what it did |
| A7 | **New files follow the folder convention** | A new subfolder inside a model folder must match what the other columns use (`research/`, `theory/`, `data/`, `plots/`, `xparameters/`, `utils/`). A novel folder name is a convention fork. **The `utils/` rule holds at two levels.** At the model root: the launcher, the medium and the engines stay at the top, and supporting scripts (instrumentation, plotting, sampling, monitoring) go in `<model>/utils/`. Inside `xparameters/`: only launchable xperiment parameter modules, each defining `XPARAMETERS`, sit at the top level, and supporting scripts go in `xparameters/utils/`, because the launcher offers every top-level `.py` there as a selectable xperiment | `ls -d openwave/xperiments/*/*/` and compare, then `git check-ignore -v <new-folder>/<file>` on any folder name introduced by the PR |
| A8 | **Repository language is English** | Code comments, docstrings, and documentation in English, so every author and cold reader can read every column | scan added comment lines for non-English text |

## 4. Gate B: scope containment

The default expectation is that a model contribution touches only its own model folder.

| Question | If the answer is yes |
| --- | --- |
| Does the PR change a **T3 shared file**? | Enumerate every consumer (`grep -rn "<module>"`), and confirm the change is backward-compatible for all of them. Behaviour changes to shared code need their own justification in the PR body, not a line buried in a feature diff |
| Does the PR change a **T4 shared surface**? | Ask for it to be split out, unless it is a plain factual correction. Front-door documents carry the project's voice and are maintainer-owned |
| Does the PR change **defaults** that other experiments read? | Defaults are shared state. A default that only the new work needs belongs in the new work's own configuration file, not in the shared default block |
| Does the PR **remove an existing gate or option**? | Removing a flag (for example an instrumentation on/off switch) changes behaviour for everyone using that model. Ask for it to be preserved or for the removal to be stated as intentional |
| Does the PR **silently change a tuning constant** used by existing results? | Any constant that existing published cells were earned under is load-bearing. Changing it invalidates those cells until re-run. Ask for it to be moved into the per-experiment configuration instead |
| Does the PR **relocate, extract or replace a shared function**? | "Moved to a shared module" is a claim about behaviour, not a description of a diff. Reading the call sites cannot test it, because they look identical either way. Run the old implementation and the new one on the same inputs and compare the outputs numerically. A move that also rewrites the body is two changes wearing one commit message, and the second one is invisible |

## 5. Gate C: claim to artifact

This is the gate [`REPRODUCE.md`](../REPRODUCE.md) and [`AI_HYGIENE.md`](../AI_HYGIENE.md) exist to enforce, and it is where most real problems surface.

**The rule: do not read the numbers, recompute them.** If the PR ships the data, write your own short script, from the raw artifact, using your own definition of the metric, and compare. This is the [adversarial audit](../AI_HYGIENE.md#1-the-stance) applied to review.

| # | Check | Bar |
| --- | --- | --- |
| C1 | **Every claim links something runnable** | A number in a research note traces to a script in the PR, or to a configuration file plus a command. Prose-only numbers do not merge |
| C2 | **The analysis step is shipped too** | Shipping raw logs is not enough if the note's tables were produced by an unshared script. The script that turns logs into the reported table is part of the claim |
| C3 | **The reported numbers match the shipped data** | Recompute at least the headline metric and one table. Report any disagreement with both numbers |
| C4 | **The metric is defined unambiguously** | "Drift" must say whether it is a Euclidean norm, a per-axis component, or a maximum over axes. Mixed definitions inside one document are a defect even when each number is individually defensible |
| C5 | **Internal consistency** | The same quantity carries the same value in the abstract, the summary table, and the detail table. Cross-check them against each other before checking either against the data |
| C6 | **The reproduction route is written down** | A reader with a clean clone can get from the claim to the command. Per [`REPRODUCE.md`](../REPRODUCE.md) that lives in the research note, once |
| C7 | **Claim strength matches evidence strength** | "Within the explored parameter range, X" is a result. "Proof of X" from an empirical sweep is not. Words like *proof*, *confirmed*, *uniquely*, and *demonstrates* each need the artifact that earns them |

## 6. Gate D: the adversarial pass

Gate C asks whether the numbers are real. Gate D asks whether they mean what the author says they mean. Run this on any PR that claims a physics result.

| # | Question | Why it catches things |
| --- | --- | --- |
| D1 | **Is the proposed mechanism actually wired in the code path that produced the result?** | Follow the named mechanism from the configuration file, through the launcher, into the kernel, and confirm it is non-zero for the mode actually used. A new term added to an `if / elif` chain with no branch for the selected mode contributes exactly zero |
| D2 | **Is the independent variable coupled to the dynamics at all?** | If the claim is "X selects the outcome", find the term through which X enters. If diagnostics are identical across every value of X, the likeliest explanation is that X is not in the loop, not that the diagnostic is insensitive |
| D3 | **Is the discriminator larger than the known systematic error?** | If the note itself lists integrator drift, rounding, or a truncation asymmetry as open bugs, compare their size to the gap between the "pass" and "fail" cases. A signal of the same order as the acknowledged error is not a result yet |
| D4 | **Do the knobs manufacture the outcome?** | Damping, clamping, a reduced timestep, and a shortened run all suppress the very motion being measured. Check what the tuning constants do over the full run length before accepting "it stayed put" as physics |
| D5 | **Has anything converged?** | Plot the discriminating metric against time. If every configuration is still rising monotonically at the last sample, the finding is a rate difference, not a stability difference. Say so |
| D6 | **Is the run long enough in physical time?** | Step counts are not time. A run at a tenth of the usual CFL safety factor covers a tenth of the physical duration for the same step count |
| D7 | **Do any two runs agree suspiciously well?** | Bit-identical trajectories from configurations that should differ are a wiring diagnostic, not a coincidence |
| D8 | **Are the controls controls?** | A control must differ from the test case in exactly the intended variable. Check the configuration files, not the file names |
| D9 | **What would falsify this?** | If the note cannot say, ask. A model author who can name the falsifier is describing a result; one who cannot is describing a hope |
| D10 | **Can the shipped self-checks fail?** | Mutation-test every line a script prints as PASS: change the thing it checks to something wrong and confirm it goes red. A check whose two sides evaluate the same expression always passes, and to a later reader it is indistinguishable from a verified result. Where a quantity has no independent target to compare against, the honest label is *asserted*, not a self-check that cannot discriminate |

Findings from this gate are **questions to the author**, not verdicts. The author owns the physics; the reviewer owns the demand that the claim and the artifact agree.

## 7. Gate E: MODELS.md cell changes

[`MODELS.md`](../MODELS.md) is the platform's product. Cells move in both directions, and they are earned one at a time.

| Situation | What the maintainer does |
| --- | --- |
| The PR **proposes a cell change** | Open the linked artifact and re-derive the icon from it, using the status semantics in [`MODELS.md`](../MODELS.md). A cell claim is not accepted on the PR body's description of the artifact |
| The PR **produces a result but proposes no cell change** | Correct default for an external contributor. The maintainer decides separately whether the result now merits a cell move, in a maintainer PR, after Gates C and D pass |
| The result is real but **weaker than the cell it would claim** | ⚠️ partial with the caveat inline, not ✅. The caveat belongs in the cell, not in an appendix |
| The result is a **documented negative** | ❌ is a result. It merges on the same evidence bar as a positive |
| The result **contradicts an existing cell in the same column** | Both cells are the author's. Route through the column's author before either changes |
| The PR is **mid-flight work** | 🔶 lives on per-model pages, not in the shared matrix. 🚧 stays until something runnable exists |

### 7.1 The MODELS.md linter

**Run `python3 dev_docs/utils/check_models_md.py` before merging anything that touches [`MODELS.md`](../MODELS.md), and read its output rather than only its exit code.** Nothing else enforces it: the repository has no CI, so this script runs when a reviewer runs it and at no other time. That is a deliberate choice (the checks are instant and a reviewer is already at a terminal), and it has one failure mode, which has already happened once: the `regime` column was added, the script's positional table detection stopped finding the summary-status table, and it sat reporting 131 violations that nobody saw because nobody invoked it. A check that is not part of a procedure is not a check.

What it covers, so a reviewer knows what it does not:

| # | Check | Catches |
| --- | --- | --- |
| 1 | Cell budget | A per-model summary cell over 65 words of prose (links, status tag and `<br>→` pointer tails excluded) |
| 2 | Icon sync | The at-a-glance matrix disagreeing with the same criterion's status tag in the model's own table, in either direction, including rows missing from one side |
| 3 | Score-board | A count that does not equal the tally of that icon over that model's rows, a total that does not equal the criteria count, or an icon used in rows with no score-board row |
| 4 | Regime | A criterion whose `regime` is not `static`, `dynamic` or `both` |
| 5 | Simplest test | A criterion with an empty test, or a criteria set that does not match the matrix in either direction |
| 6 | Row shape | A data row whose cell count differs from its header, which is what an unescaped `\|` inside a cell looks like from the parser's side |

It does **not** check prose accuracy, link targets, or whether a cell's claim is supported by the artifact it links. Those are [Gate C](#5-gate-c-claim-to-artifact) and [Gate E](#7-gate-e-modelsmd-cell-changes), and they are yours.

⚠️ **Check 1's 65 and the roadmap linter's 65 ([§ 7.2](#72-the-roadmap-linter)) are the same number in different units.** This linter counts prose only; the roadmap linter counts link labels too and has no status tag to strip, so a cell at 65 here renders about a third larger than a row at 65 there. That is derived, not an accident: the two-column per-model table gives its summary column 1.36× the width of a four-column roadmap row, which absorbs the difference at equal rendered lines. Do not "fix" one to match the other. Measure: `python3 dev_docs/utils/models_cell_stats.py` (derivation of record: [T3](tasks/t3_task_details.md)).

Two operational notes. A clean run prints `clean` and exits 0; anything else lists line-numbered violations. And if a PR legitimately introduces a new criterion-level column (the way `regime` was), the script will refuse it by name until the column is registered in `REGIMES`-style fashion beside it, which is intentional: adding a column must not be able to silently switch a check off.

### 7.2 The roadmap linter

**Run `python3 dev_docs/utils/check_roadmaps.py` before merging anything that touches a roadmap.** Same no-CI caveat as 7.1: it runs when a reviewer runs it. It enforces the word budgets in [`ROADMAP_STANDARDS.md`](ROADMAP_STANDARDS.md), whose premise is that a roadmap row is a preview and the task document is the record. A row that needs more than its budget is a row whose content belongs one link deeper, so the fix is almost never a bigger budget.

The budgets it checks: description 65 words, title 15, every other cell in the row 35, section blockquote 50, intro blockquote 80, change-log entry 200. It also requires a column named exactly `Description` in each task table, and reports rows whose cell count does not match their header, which is what an unescaped `|` looks like from the parser's side. `ARCHIVE` and `LEGACY` sections are skipped as frozen history.

## 8. Gate F: other authors' work

Every column carries a person's name. Protecting that is a maintainer duty, and it is the part of review an author cannot do for themselves.

| # | Rule |
| --- | --- |
| F1 | A contributor extending a column they do not own contributes **under their own name**, and the research note says whose extension it is. The scoring rules for borrowed structure are in [`CROSS_MODEL_TESTING.md`](CROSS_MODEL_TESTING.md) |
| F2 | A change that **invalidates results another author earned** (a tuning constant, a shared kernel, a removed diagnostic) is a blocking finding until either the results are re-run or the change is scoped so they are untouched |
| F3 | Never let a review thread turn into a defence of someone else's physics. The maintainer checks reproducibility and containment; the author answers challenges to the model, per [`ONBOARDING_MODELS.md`](../ONBOARDING_MODELS.md) |
| F4 | A **regression in shared infrastructure** is the maintainer's to catch, because no single author will see it. Trace the callers |

### 8.1 When the contributor is not the model author

Check this on every PR that touches `openwave/xperiments/<model>/`: is the contributor the author of that column? The Identity table in the model briefing is the record (`Author` and `Author contact` rows; an `Extension` row names contributors who extend the column without owning it).

**If the contributor is not the column's author, the review carries a model-author note**: a short block naming the author, **@-mentioning the author's GitHub handle**, and stating which findings are the author's call and which are not. The contributor learns from the review itself who holds authority over what, rather than discovering it after a merge.

The @-mention is the mechanism, not a courtesy. GitHub notifies on `@handle` and on nothing else: writing the author's name in prose delivers no notification, and neither does a link to the briefing. Rules for it:

| # | Rule |
| --- | --- |
| M1 | Take the handle from the model briefing's **`Author contact`** row, which is where it is maintained. Do not guess it from the author's name |
| M2 | Put the `@handle` in the **first sentence** of the model-author note, not buried in a table cell or a closing line. It is what makes the routing real |
| M3 | Mention **once**, in the review body. Repeating it in every subsequent comment turns the channel into noise, which is the thing [`AI_HYGIENE.md`](../AI_HYGIENE.md) § 3 warns about |
| M4 | Say in the same sentence **what is being asked**. A bare mention reads as a summons; "tagging @handle on the two rows below, the rest is mine" is actionable |
| M5 | If the briefing carries **no handle**, that is a gap in the briefing. Route by the contact that does exist, say so in the thread, and open a follow-up to fill the Identity table |
| M6 | Mention the author, not the whole org or a team alias. One named person owns the column |

What that note routes is decided per finding, not per PR:

| Finding type | Author-gated? | Who decides |
| --- | --- | --- |
| Engine, kernel, or solver change to the column | ✅ yes | The author. It changes what the column *is* |
| Change to defaults or tuning constants existing cells were earned under | ✅ yes | The author, since it is the author's results that are revalued |
| A headline physics claim about the column, or a `MODELS.md` cell move | ✅ yes | The author, on the science. The maintainer still holds the evidence bar in [Gate E](#7-gate-e-modelsmd-cell-changes) |
| Reinterpreting what an existing cell or a past result means | ✅ yes | The author. This is intent and provenance, and [`AI_HYGIENE.md`](../AI_HYGIENE.md) makes it structurally unanswerable by anyone else |
| A new experiment, configuration, or geometry added alongside the existing ones | ❌ no | The maintainer. Additive, changes nothing the author already owns |
| Code regression, dangling import, encoding defect, formatting | ❌ no | The maintainer. Mechanical correctness is not a physics question |
| Documentation, cross-links, citation registration, folder conventions | ❌ no | The maintainer |
| Numbers in a note disagreeing with the shipped data | ❌ no | The maintainer, against the artifact. Only the *interpretation* of the corrected numbers is author-gated |

**A PR with no author-gated findings does not need the author in the loop.** Merge it. Pulling an author into every additive contribution burns the one channel that matters for the questions only the author can answer, which is the failure mode [`AI_HYGIENE.md`](../AI_HYGIENE.md) § 3 warns about ("keep such asks as small as possible, one item per ask").

**Where there is at least one author-gated finding, the author has the final say on those findings**, and the PR does not merge past them on the maintainer's opinion alone. The maintainer's own blocking findings (the ❌ rows above) stay the maintainer's and are fixed regardless.

Author-gated does not mean author-blocked forever. If the author does not respond, the maintainer may merge the non-gated part and keep the gated part open as an issue, saying so in the thread.

## 9. Gate G: policy sweep

Fast pass over the standing policies. Most PRs clear this in a minute.

| Policy | What to check |
| --- | --- |
| [`AI_HYGIENE.md`](../AI_HYGIENE.md) | Claims script-backed, not model-fluent. Status attached to AI-assisted findings. No aggregate self-ranking or cross-program scoreboards. Substantive claims carry an adversarial check, and the review itself is one |
| [`CODE_OF_CONDUCT.md`](../CODE_OF_CONDUCT.md) | Applies to the review thread as much as the contribution. Critique the artifact, never the contributor |
| [`CONTRIBUTING.md`](../CONTRIBUTING.md) | DCO, Apache 2.0, fork and branch flow, PEP 8 and Black as the style target (see [fairness](#11-fairness-rules) on how hard to press this) |
| [`REPRODUCE.md`](../REPRODUCE.md) | Reproduction route recorded once, in the research note. Nothing result-shaped lands in `REPRODUCE.md` itself |
| [`ONBOARDING_MODELS.md`](../ONBOARDING_MODELS.md) | For a new column: scaffold set complete, briefing filled, author responsibilities accepted, collaborator invitation sent |
| [`MARKDOWN_STYLE_GUIDE.md`](MARKDOWN_STYLE_GUIDE.md) | Blank lines around headings and lists, single trailing newline, fenced blocks carry a language |
| [`METHOD_NOTE.md`](METHOD_NOTE.md) | Required shape for anything that will be reported to a theory owner or an external physicist: equations first, equation-to-code map, the not-computed list, the audit recorded |
| [`CODING_STANDARDS.md`](CODING_STANDARDS.md) | Naming, docstrings, type hints, Taichi kernel guidance |

## 10. Maintainer edits

When to fix something yourself instead of asking for it.

A review round trip costs both sides more than most of the fixes it asks for. "Line 10 is CP1252, please re-save as UTF-8", a day of waiting, then a re-review, is more expensive than saving the file. So the default is: **if a fix is mechanical and unambiguous, make it, and say so in the review.**

This is not a licence to rewrite a contribution. The test is whether the fix requires knowing *why* the contributor did something.

| Fix it yourself | Leave it to the contributor |
| --- | --- |
| Encoding, BOM, trailing newline, formatter output | Anything whose right answer depends on their intent |
| A duplicated definition where the later one silently wins | Which of two implementations to keep, when both are plausible |
| A label that contradicts the value beside it, when only one reading is possible | The same contradiction, when either half could be the error |
| Restoring comments or documentation the diff removed in passing | Removing or rewriting their prose, in code or in documents |
| A missing gitignore entry, a stale citation, a dangling reference | Tuning constants, thresholds, and anything the physics rests on |
| The mechanical half of a finding, so only the real question is left | The finding itself, when it is a design or physics decision |

### 10.1 Prerequisites, in order

1. **Confirm you can push.** A PR from a fork is editable only if the contributor left the box ticked:

```bash
gh api repos/<owner>/<repo>/pulls/<N> -q '{editable: .maintainer_can_modify, head: .head.label}'
```

If that returns `false`, this section is moot: post the findings with commands instead.

2. **Push to their branch**, not to one of your own:

```bash
git remote add <contributor> https://github.com/<contributor>/<repo>.git
git push <contributor> <local-branch>:<their-head-ref>
```

3. **Check what their head ref is.** If it is their fork's default branch (`head: <user>:main`), say in the review that they must `git pull` before their next commit, or their default branch diverges from the PR. It is also the natural moment to suggest topic branches for future PRs, as a kindness rather than a rule.

### 10.2 Rules that keep this safe

| Rule | Why |
| --- | --- |
| **Announce every edit in the review**, what changed and why | The point of review is that the standard becomes visible. Silently fixing and merging teaches nothing, which is exactly what the old "ask, do not fix" rule protected. Writing it down protects the same thing at a fraction of the cost |
| **Formatter runs go last, in their own commit** | Reformatting a dozen files buries every semantic change before it in the diff, and it conflicts with whatever the contributor has locally |
| **Never put a repository-wide change inside a contributor's PR** | It inflates their blast radius, entangles the review, and stalls the platform fix if the PR stalls. Conventions are maintainer decisions and the history should say so |
| **When one change spans both sides, land the platform half on `main` first** | A convention that must hold across every model, applied to a file that exists only in the PR, is two changes. Enforce it on `main`, then apply it to the PR branch. The PR picks up the rest on merge |
| **Relocating a module can change what it computes** | Anything resolving a path from `Path(__file__).parent` silently means something else once the file moves a level down. Before pushing a move, grep the moved files for `__file__` and re-verify every path they derive, including the ones a second module derives independently and expects to match |
| **A fix that needs a guess is not mechanical** | If you cannot state why they wrote it that way, you are not fixing it, you are overwriting it |

## 11. Fairness rules

These exist so review stays honest in both directions.

| # | Rule |
| --- | --- |
| G1 | **Measure the PR against the enforced baseline, not the aspirational one.** Before citing a standard, check whether `main` itself passes it. If the repository has drift, that is the maintainer's debt, not the contributor's blocker. Say so out loud in the review |
| G2 | **Separate blocking from optional, explicitly.** Every finding is labelled: blocking, requested, or note. A contributor should be able to count the blockers on one hand and know exactly what merges the PR |
| G3 | **Fix it yourself when the fix is mechanical, and say so.** A review round trip costs both sides more than most of what it asks for. Make the unambiguous fixes on the contributor's branch, list every one of them in the review, and leave the findings that need their intent. Rules in [§ 10](#10-maintainer-edits) |
| G4 | **Give the command, not just the complaint.** The one-time DCO config, the exact reformat invocation, the recompute script. Reviews that ship commands get resolved in one round |
| G5 | **A negative result is not a weak PR.** Neither is a small one |
| G6 | **State what you verified and what you did not.** A review that implies more checking than happened is the same failure mode as an overstated claim |
| G7 | **Findings from an AI-assisted review are findings, not verdicts,** until they carry the command that reproduces them. Same contract as everything else here |

## 12. Verdict and how to write it

| Verdict | When | Action |
| --- | --- | --- |
| ✅ **Approve and merge** | All gates clear, or only notes remain | Merge. Record anything learned in the [lessons log](#14-lessons-log) |
| 🔶 **Approve with follow-ups** | Gates clear; requested items are real but not load-bearing | Merge, and open issues for the follow-ups so they do not evaporate |
| ⚠️ **Changes requested** | Blocking findings exist, all of them fixable by the contributor | Post the findings with commands. Re-review only the deltas |
| 🚧 **Split requested** | The contribution is sound but mixes tiers (model work plus shared-surface edits) | Ask for the shared-surface part to come out; merge the rest |
| ❌ **Decline** | Provenance cannot be established, or the contribution cannot be made reproducible | Rare. Explain which gate failed and what would change the answer |

**Review structure that works:**

1. One sentence on what the PR does and what it claims.
2. The tier split, so the blast radius is visible.
3. ✅ what passed, briefly. Contributors deserve to see the checks that cleared.
4. Blocking findings, each with: the file and line, what you ran, what you got, and what would fix it.
5. Requested and note items, clearly separated.
6. The physics questions from Gate D, framed as questions to the author.
7. The **model-author note**, when the contributor is not the column's author ([§ 8.1](#81-when-the-contributor-is-not-the-model-author)): the author's `@handle` in the opening sentence so the notification actually fires, then which findings are the author's call.
8. What you did not check.

Findings are about artifacts. "This table disagrees with the shipped data" is a finding; "you were careless" is not. Where a finding could read as a challenge to the author's model rather than to the artifact, say which one you mean.

## 13. Command appendix

Fetch and isolate the PR:

```bash
gh pr view <N> --json number,title,author,body,files,changedFiles,additions,deletions,headRefName,mergeable
gh pr checks <N>
git fetch origin pull/<N>/head:pr-<N>
git log main..pr-<N> --format='%H%n  %an <%ae>%n  signoff: %(trailers:key=Signed-off-by,valueonly)%n  %s'
git worktree add --detach /tmp/pr-<N> pr-<N>     # review without disturbing your tree
```

Blast radius and size:

```bash
git diff main...pr-<N> --stat
git diff main...pr-<N> --name-only --diff-filter=D          # deletions
git ls-tree -r -l pr-<N> | awk '{print $4, $5}' | sort -rn | head -20
```

Encoding audit (A2), the check that catches what `grep` cannot report:

```bash
python3 - <<'PY'
import subprocess, os
files = subprocess.run(["git","diff","main...pr-<N>","--name-only","--diff-filter=d"],
                       capture_output=True, text=True).stdout.split()
for f in files:
    if not os.path.exists(f):
        continue
    d = open(f, "rb").read()
    if d.startswith(b"\xef\xbb\xbf"):
        print("BOM:", f)
    try:
        d.decode("utf-8")
    except UnicodeDecodeError as e:
        print("NOT UTF-8: %s byte %s line %d" % (f, hex(d[e.start]), d[:e.start].count(b"\n") + 1))
PY
```

Dangling references after a deletion (A6):

```bash
for f in $(git diff main...pr-<N> --name-only --diff-filter=D); do
  base=$(basename "$f" .py)
  echo "--- $base"; git grep -n "$base" pr-<N> -- '*.py' '*.md'
done
```

Style and import sanity, on the PR worktree:

```bash
cat filelist.txt | tr '\n' '\0' | xargs -0 python3 -m black --check
python3 -m py_compile <changed files>
python3 dev_docs/utils/check_models_md.py    # mandatory when MODELS.md is touched, see 7.1
python3 dev_docs/utils/check_roadmaps.py     # mandatory when a roadmap is touched, see 7.2
```

### Recording the verdict: submit it as a review, not as a comment

The [verdict](#12-verdict-and-how-to-write-it) has to be submitted as a **review**, which carries a state, and not as a conversation comment, which does not. GitHub makes this easy to get wrong: both render identically and both notify, but only a review sets `reviewDecision`, shows the status badge on the PR list, and satisfies or blocks the CODEOWNERS gate. A review body that says "changes requested" while the PR state says `REVIEW_REQUIRED` leaves the contributor without the signal they watch for.

| Surface | Where in the UI | What it produces |
| --- | --- | --- |
| **Review** (correct) | **Files changed** tab → **Review changes** → choose Comment / Approve / **Request changes** → **Submit review** | A review with a state. Sets `reviewDecision`, badges the PR, drives the merge gate, and auto-dismisses on a force-push |
| Conversation comment | The comment box at the bottom of the **Conversation** tab | An issue comment. Notifies, but carries no state and does not touch the merge gate |

The radio buttons live only behind **Review changes** on the Files changed tab. There is no path to them from the Conversation tab, which is why a full review can be written and posted without ever being offered the choice. Inline comments anchored to a specific file and line are also review-only, so any finding worth pinning to the code has to go through the same button.

Terminal equivalents:

```bash
gh pr review <N> --request-changes --body "..."   # blocking findings exist
gh pr review <N> --approve       --body "..."     # all gates clear
gh pr review <N> --comment       --body "..."     # review-shaped, deliberately no state
gh pr comment <N> --body "..."                    # plain conversation comment, no state
gh pr view <N> --json reviewDecision,reviews      # verify the state actually landed
```

A long review body can be posted as a comment first and the state submitted separately with a one-line review pointing at it. Verify with the last command either way: the state is the part that is easy to lose.

**Reconcile the state at merge, not only at review time.** A `CHANGES_REQUESTED` review is not retired by the changes being made, by a follow-up comment, or by the merge itself: it stands until a later review supersedes it. So a PR whose findings were resolved as [maintainer edits](#10-maintainer-edits) merges with its record still reading changes-requested, and a contributor looking back at their own merged work sees a rejection badge on it. The last thing before or after clicking merge is therefore an approving review, and `gh pr review <N> --approve` works on an already-merged PR, so this is fixable after the fact:

```bash
gh pr merge <N> --merge                                # or the UI
gh pr review <N> --approve --body "Approving for the record: merged in <sha>. ..."
gh pr view <N> --json reviewDecision                   # must read APPROVED
```

## 14. Lessons log

One row per PR that taught us something. Newest at the bottom.

| PR | Lesson | Where it landed |
| --- | --- | --- |
| [#196](https://github.com/openwave-labs/openwave/pull/196) | A DCO block is a configuration problem, not a rejection. Shipping the exact four commands, customized to the contributor's branch, resolved it in one round. Warn that a force-push dismisses the standing approval | Intake row 2, fairness rule G4 |
| [#297](https://github.com/openwave-labs/openwave/pull/297) | Two things a diff read catches that a claim read does not: a backend left on CPU, and inline comments stripped by an auto-formatter. Those comments were there for cold readers, so their removal is a real loss even though no behaviour changed | Gate B (removing an existing option), Gate G |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | Recomputing the headline table from the shipped data disagreed with the note, and following the named mechanism through the code showed it contributing zero in the mode actually used. Both were invisible from the PR body. Also: one CP1252 byte in a shared module made `grep` silently skip the file during review | Gate C (recompute, do not read), Gate D rows 1 and 2, Gate A row A2 |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | The contributor was extending a column owned by someone else, and the PR changed that column's engine defaults while claiming its headline open problem. Nothing in the process said when the column's author has to be in the loop and when that would just be noise, so the per-finding routing was written down | Gate F § 8.1, review structure item 7 |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | The launcher listed every top-level `.py` in `xparameters/` as a selectable xperiment, so a helper module the PR added appeared as a menu entry that failed on selection. Found by running the GGUI, not by reading the diff. The rule now lives in code in all five launchers, and helpers live in `xparameters/utils/` | Gate A row A7, `_discover_xperiments` in every launcher |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | The first review asked for seven mechanical fixes and scheduled a second round to verify them. Most were an encoding byte, a BOM and a formatter run, so the round trip cost more than the edits would have. Fairness rule G3 was inverted and the maintainer-edit rules written down | § 10, fairness rule G3 |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | Two deleted files came back clean on a dangling-reference check, because both were self-contained xperiment configurations. Each defined an `XPARAMETERS` dict, so deleting them removed two entries from the launcher menu and nothing pointed at the loss. Reference checks answer "what breaks", never "what is gone" | Gate A row A6 |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | A generator extracted from one configuration into a shared module read as a clean refactor at every call site, and was not the same function. Running both on the same inputs moved every K from 2 to 10, up to 2.87 lambda, and left the K=2 to 9 cases sitting inside the first lock-in well instead of at it. Three existing experiments changed silently. Diff the outputs, never the call sites | Gate B, new row on relocated functions |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | A duplicated kernel definition was masking a compile error: Python kept the later one, so the earlier guarded version had never been compiled, and deleting the shadow surfaced a `TaichiSyntaxError` that would have failed every run. Treat pyflakes' `redefinition of unused` as blocking rather than cosmetic, and compile what the deletion exposes | Gate A, Gate C |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | First run under § 10. Maintainer edits collapsed a three-round review into one pass, and the two defects that mattered most were found by running the thing: a `os.execv` xperiment-switch path that skipped every teardown and orphaned a GUI child process, and the compile error above. Neither is visible in a diff. Also, relocating modules broke `Path(__file__).parent` resolution in three of them, the reviewer's own breakage, caught before the commit | § 10.2, Gate D |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | The first name chosen for the support-module folder, `lib/`, collided with `.gitignore` line 20 inherited from the Python template for build output. Committing it would have deleted `instrumentation.py` from four models and silently ignored the replacement. Caught by `git check-ignore` before the commit and resolved by renaming to `utils/` rather than adding ignore negations: a negation covers only the exact depths it names, and the failure it hides is silent | Gate A row A7 command, the `utils/` convention |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | The verdict went out as a conversation comment rather than a submitted review, so the PR sat at `REVIEW_REQUIRED` while its body said changes requested. The merge gate held on CODEOWNERS regardless, but the contributor never got the signal. Fixed with `gh pr review --request-changes` and written into the appendix | § 13 "Recording the verdict" |
| [#350](https://github.com/openwave-labs/openwave/pull/350) | A shipped script printed a PASS line that could not fail: its two sides evaluated the same expression, so replacing the rule under test with deliberate nonsense still reported PASS while the table filled with wrong values. Nothing downstream depended on it, but under a verification banner it reads exactly like a certified result. Mutation-testing every PASS line is now part of the adversarial pass | Gate D row D10 |
| [#350](https://github.com/openwave-labs/openwave/pull/350) | What made this review conclusive was recomputing the headline table by a genuinely different method rather than re-running the contributor's script: the group rebuilt as explicit quaternions with characters from Burnside class-sums, against the PR's McKay recursion. Agreement on 9/9 rows then meant something. "Recompute, do not read" is only as strong as the independence of the second route | Gate C, the recompute rule |
| [#350](https://github.com/openwave-labs/openwave/pull/350) | Also the good case worth naming: the contributor raised a cross-model question as a platform issue *before* the work depended on the answer, and took the two family questions to the column authors directly. That is what made the author-gated findings empty and the review light. Sequencing, not effort, is what keeps [Gate F](#8-gate-f-other-authors-work) cheap | Gate F § 8.1, as the worked example |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | The mirror image of the row above, on the same PR: this time the state was submitted correctly and then outlived the verdict. The findings were resolved as maintainer edits and the PR was merged, but no later review superseded the standing `CHANGES_REQUESTED`, so a first-time contributor's merged work carried a rejection badge. Nothing about making the changes, commenting, or merging retires a review state; only another review does | § 13 "Reconcile the state at merge" |
| [#340](https://github.com/openwave-labs/openwave/pull/340) | A claim we had written ourselves, in a maintainer commit, was the thing that blocked the merge. The docstring asserted only K=10 sat at the lock-in wells; measuring all 45 pair separations showed the opposite (K=2..4 entirely on the well, K=10 at none of 45) and that the band we attributed to K=2..9 was really K=11's. `git blame` on the contradiction before routing it saved a review round, because the answer was that the wrong half was ours | Gate C, [§ 10](#10-maintainer-edits) |

---

## See also

| Doc | Why |
| --- | --- |
| [`../CONTRIBUTING.md`](../CONTRIBUTING.md) | The contributor's side of this process: setup, fork and branch flow, DCO |
| [`../AI_HYGIENE.md`](../AI_HYGIENE.md) | The adversarial-audit cardinal rule that Gates C and D implement |
| [`../MODELS.md`](../MODELS.md) | The coverage matrix, its status semantics, and the light-review promise |
| [`../REPRODUCE.md`](../REPRODUCE.md) | What "reproducible" means here, and where reproduction commands live |
| [`../ONBOARDING_MODELS.md`](../ONBOARDING_MODELS.md) | Model-author responsibilities and the scaffolding sequence for a new column |
| [`CROSS_MODEL_TESTING.md`](CROSS_MODEL_TESTING.md) | Scoring rules when a contribution reaches into another column |
| [`METHOD_NOTE.md`](METHOD_NOTE.md) | The reporting shape required for model-owner-facing results |

---

**Deep readers and AI agents**: the full map of OpenWave's key documents, and the order to read them in, is in [`../CLAUDE.md`](../CLAUDE.md). The AI-collaboration contract is [`../AI_HYGIENE.md`](../AI_HYGIENE.md).
