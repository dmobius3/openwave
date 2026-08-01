# The context-isolated clean room

> The platform's procedure for running a reproduction in a context that provably could not
> have been handed the answers. Written from an exercise, not in the abstract: the M8.5-A
> reproduction ran this procedure end to end (packet [#389](https://github.com/openwave-labs/openwave/pull/389),
> commitment [#391](https://github.com/openwave-labs/openwave/pull/391), adjudication
> [#392](https://github.com/openwave-labs/openwave/pull/392), audit
> [#394](https://github.com/openwave-labs/openwave/pull/394)), and every rule below either
> held there or was added because something slipped. Planning record and deviations log:
> [T4](tasks/t4_task_details.md); the exercise's own record:
> [`m8_5_task_details.md`](../openwave/xperiments/m8_mit/research/tasks/m8_5_task_details.md).

## 1. The claim ladder

Four terms, a ladder rather than synonyms. Using the wrong rung oversells a result.

| Term | What it names |
| --- | --- |
| clean room | the ENVIRONMENT and procedure: the isolated directory, the audited packet, the ordering record. What this standard specifies |
| context firewall | the RULE SET the room enforces: fresh context, no quarantined access, no web tools, the consulted-files manifest |
| context-isolated independent-method reproduction | the CLAIM a successful run earns, and the ceiling |
| blind | a STRONGER claim the room alone never grants. Blind is a property of what the verifying context has seen, not of the environment; a clean room makes isolation auditable, it cannot erase prior exposure. For an AI implementer the opaque training corpus caps the label structurally: prior exposure to published values cannot be excluded, so isolation buys provenance, never the label |

## 2. When to run one

When a result must be reproduced by a context that the record can show was never handed the
targets: an author-supplied computation reproduced in-platform, a single-implementation
result re-derived before a claim rests on it, any case where the implementing context would
otherwise share state with a context that knows the answers. The cost is real (a sealed
packet, ordered commits, an operator at the keyboard), so the trigger is the claim, not
habit: if nobody will rely on the isolation, a plain independent implementation is cheaper.

## 3. Building the room

### 3.1 Location, and the ancestor-path check

The room is a directory outside every repository working tree with no agent-instruction file
anywhere on its ancestor path. The check is a procedure, not advice:

```bash
d=/path/to/room
while [ "$d" != "/" ]; do
  for f in "$d/CLAUDE.md" "$d/.claude/CLAUDE.md"; do [ -e "$f" ] && echo "FOUND: $f"; done
  d=$(dirname "$d")
done
```

Read every file it finds, not merely whether one exists. The M8.5-A near-miss: the obvious
sibling-of-the-repository folder sat under a 7-line instruction file whose entire content
told the agent to go read the quarantined tree, on the first turn, before opening the
packet. Nothing inside the room would have revealed the breach. A future room is placed by
running this check, never by copying a path that worked before (the concrete path is a local
choice and stays out of the repository).

### 3.2 The interpreter is a second, independent route

The ancestor-path check covers what loads by instruction; the Python environment is a
separate route into the quarantined tree, and checking one says nothing about the other. An
editable install makes the tree reachable through a bare `import` regardless of where the
room sits, and an environment named after the repository puts that name in the shell prompt.
Verify from the room's shell, before launch:

```bash
python3 -c "import <platform_package>" 2>&1 | tail -1   # must FAIL
pip show <platform_package> | head -3                    # must return nothing
```

Record the interpreter, versions, and available numerics libraries for the environment
record. The M8.5-A launch shell arrived with the platform installed and switched
environments on the spot; the rule exists because the miss was live, not hypothetical.

### 3.3 What loads unavoidably

The user-global instruction file loads into any session. Grep it for the subject vocabulary
before the run and disclose it in the consulted-files manifest as a named load. Read the
grep's summary with suspicion: "nothing relevant" is the conclusion to be most careful with,
since the file is written for a general working context and nobody edits it with a clean
room in mind. At M8.5-A it carried nothing answer-bearing and still named the platform
repository several times. Project memory is keyed to the working directory, so a fresh
folder starts with none.

### 3.4 Launch

```bash
claude --disallowedTools "WebSearch,WebFetch"
```

In DEFAULT permission mode, never a bypass mode. The approval prompt is what makes a tool
call reaching outside the room visible to the operator, and it is the one containment guard
that does not rely on the implementer's cooperation. The room is isolated by withheld tools
and by instruction, not by a sandbox: the filesystem outside it is reachable in principle,
and the prompts are what stand in the way. The implementer writes only inside the room; a
maintainer copies artifacts out and performs every commit.

## 4. The packet

The packet is everything that enters the room. A maintainer who knows the answers assembles
it, which is exactly what the audit anticipates: the containment is that no part is capable
of carrying a target, and every part gets a MECHANICAL audit rather than a judgment call.

| Part | Content rule | Audit |
| --- | --- | --- |
| the frozen protocol | a verbatim copy of the repository file | `diff` against the repository; identity is the whole audit |
| the group/data input | raw inputs only, no derived annotations, no labels, dimensions, or precomputed values | run the checks from the packet alone (closure counts, exact norms, a leakage scan for the derived-object vocabulary), each mutation-tested so no green line is a line that cannot fail |
| the task file | operational instructions only; its mathematical content is generic identities, never values | read against the leakage vocabulary; committed so the record shows what the room was told |

**Generators, not closures.** Hand over the seed, not the built object: an input the
implementer must construct from makes the construction gate a check that can fail, where
handing the full object makes it true by construction. Where the protocol permits supplying
a derivable quantity (M8.5-A: the group order), withholding it is the stricter choice, and
the withholding is STATED in the packet so silence does not read as an incomplete delivery.

**Author-side evidence is not the audit.** An author may ship a verification script and
report alongside the packet; they stay outside the room, and the maintainer audit runs from
the packet alone. An audit that reads the supplier's report and concurs is a review of the
report. Set the two results side by side afterwards as a cross-check.

**Canonicalize, then hash.** Record the INCOMING hash of the delivered bytes, declare a
canonical serialization, and issue the AUTHORITATIVE hash of the canonical form, stating
whether canonicalization changed the bytes (at M8.5-A it was a no-op because the author
supplied exact arithmetic; the chain is recorded either way). A hash over an uncanonical
rendering pins a transcription, not the object.

**The task file's small print, learned the hard way:**

| Instruction to include | Why |
| --- | --- |
| the room boundary stated as an obligation ("this directory is your whole world, and that is a rule, not a description of the machine"), with an out-of-room read named recoverable and its concealment named not | an unstated rule can be neither honored nor breached; the filesystem is not sandboxed |
| protocol links point into a repository the implementer does not have and must not seek | a verbatim protocol copy carries relative links; say they are dead ends |
| print repo-safe relative paths in all output | a script echoing an absolute path forces a disclosed redaction at copy-out, since the room path stays out of the repository |
| fail loud, report the honest outcome, tolerances fixed in source | a disagreement is a finding, not a failure of the task |
| a consulted-files manifest including anything that loaded WITHOUT being asked for | the manifest is corroborated against the transcript at copy-out |

**Commit before the room opens.** The packet, the task file, the audit script, and the audit
output land in the repository first, so the timestamps predate the run and nothing reads as
retrofitted to the outcome. Prepared artifacts live in the repository, never only in a
scratch directory or a conversation.

## 5. The opening prompt

The message that starts the room session is deliberately near-empty:

> Read `PROTOCOL.md`, `TASK.md`, and `GROUP_INPUT.md` in this directory, in that order, then
> do what `TASK.md` says. Everything you need is here; nothing outside this directory is
> available or needed.

Anything beyond this belongs in the task file, where it enters the audited packet and the
commit record. An instruction delivered as chat text is invisible to the packet hash, absent
from the manifest, and composed by a context holding the targets, which is the one channel
the firewall forbids. The near-empty prompt is a firewall property, not a style preference.

## 6. The operator runbook

The human at the room's keyboard is a firewall component. The rules, applied per event:

| Event | Operator action |
| --- | --- |
| at launch, before the opening prompt | verify the room's file hashes against the sealed record; verify the session banner's working directory is the room; record the banner's model for the environment record |
| prompt for a path inside the room, or for running the implementation there | approve |
| prompt for any path outside the room, or any network call | decline, and log the prompt text in the deviations record |
| no blanket grant, ever | never pre-approve tool classes, never "don't ask again", never a bypass mode; each prompt is read before it is answered, and the reading is the guard |
| the implementer asks a question | never answered from the operator's knowledge, however innocuous; relayed to the maintainer session, answered from the protocol text or not at all, logged as a deviation |
| an apparent stall | check for a PENDING PERMISSION PROMPT before interrupting: a frozen token counter with no visible prompt was, at M8.5-A, a prompt that had not surfaced in the terminal. If truly stalled, one interrupt and the single word `continue`, nothing more |
| the implementer finishes | the final message is relayed verbatim; nothing in the room is touched until copy-out |

Remote approval (answering prompts from a phone) is permitted: it relocates the operator's
channel and grants the session nothing. The two rules that erode on a small screen are the
ones to hold consciously: every prompt read in full, and no question ever answered from the
phone.

Information flow between the room and the maintainer session is one-directional: the
maintainer session may see anything from the room (it is already disqualified); the room
receives only the audited packet, fixed before launch and never amended mid-run. A declined
prompt is evidence, not only prevention: a room session that asked for an outside path has
stated something about its behavior under isolation, and the record wants that even when the
answer was no.

## 7. Copy-out and the commitment

Nothing is unsealed until the commitment is committed and merged. The commitment carries:

| Content | Rule |
| --- | --- |
| hashes of every deliverable | computed independently by the maintainer at copy-out and matched against the implementer's own record |
| the environment record and consulted-files manifest | folded in, with the implementer's prior-knowledge disclosure |
| the transcript check | extract every tool call from the session transcript and classify its paths against the manifest; the manifest is corroborated, never merely attested. Session-temp writes of the run's own output are disclosed, not hidden |
| any redaction | disclosed precisely, with BOTH hashes recorded (the original stays the commitment; the transcript retains the original bytes) |
| the pre-declared choices | anything the protocol requires fixed before unsealing (M8.5-A: whether the optional module ran) is declared here and unavailable afterwards |
| the operator log | launch checks, prompts, stalls, deviations |

## 8. Unseal and adjudicate, in a second commit

The comparison lands as its own separately merged change that NAMES the commitment it
postdates, so the ordering is auditable from the record alone: the result was public before
the comparison, or its target transcription, existed anywhere.

| Rule | Why |
| --- | --- |
| the adjudicator makes its OWN transcription of the pinned target | independent of any earlier reconstruction's fixture, so a shared typo cannot be silently self-consistent |
| reference data (a label-to-signature map) may come from the unsealed objects, declared | that is what the protocol opens them to after commitment; check it for internal consistency before any cell comparison |
| the comparison itself carries mutations on BOTH sides | perturb one transcribed target cell and one candidate cell; the comparison must go red both times, so it is a check that can fail |
| exact comparison where the object is exact | a tolerance where an equality belongs hides a class of defect |

## 9. Close-out

The method note keeps the clean-room draft byte-frozen (its hash pinned in the commitment;
the status block says so) and appends the maintainer layer: the adjudication outcome, and
the recorded adversarial audit by an independent second agent briefed to REFUTE with its own
methods (at M8.5-A: its own group construction, its own character-extraction route, its own
tables; eight claims attacked, none refuted, two strengthened, six unflagged weaknesses
recorded).

A communications rule that bit during the exercise: **an edit is not a notification, and a
notification is not a record.** Editing a merged PR body notifies nobody, and a mention
added by edit does not notify either. The durable record goes in the body; anything that
needs a human's answer goes in a NEW comment.

## 10. Teardown, last

Delete the room only after this standard's per-run write-up (or its refutation) is done: the
room is the last place to check a detail the transcript renders ambiguously. A room left on
disk is a room that gets reused, and its second occupant is not context-isolated: by then
the folder holds the packet, the implementation, and the unsealed comparison. Keep the
session transcript; it lives under the agent's project-history directory, outside the room,
and it is what makes the manifest checkable rather than an attestation.

## 11. What is inherited, and what was M8.5-A's

| Any column inherits | M8.5-A-specific |
| --- | --- |
| the claim ladder, the room build (§ 3), the packet discipline (§ 4), the opening prompt, the operator runbook, the commitment-before-unseal ordering, the two-commit adjudication, the teardown rule | the particular protocol, gates G1-G12, the group input and its `Q(φ)` arithmetic, the § 6 coexact module and its verdict categories, the pinned pre-registration target |
| the mutation-test posture: every green line shown able to go red, including the audit's own checks and the adjudication comparison | the specific mutation lists |
| the reproduce-vs-derive ceiling (§ 12) | the worked instance of it |

## 12. The ceiling: reproduce, not derive

Per [`PR_REVIEW_STANDARDS.md`](PR_REVIEW_STANDARDS.md) Gate B (obligations row): the
platform reproduces computations to a frozen spec and does not derive on request. A part of
a protocol requiring ORIGINAL derivation is optional and declinable per run, and declining
is a recorded outcome, never a breach. The worked instance: the M8.5-A § 6 coexact module
pre-declared `not resolved` as a permitted verdict; the run chose to derive and its
derivation survived the audit, but the permission to decline is what kept the obligation
priced. The claim ceiling for any AI-implemented run stays context-isolated
independent-method reproduction; blind is never granted by machinery (§ 1).

## 13. Failure modes actually hit

Every entry below happened during the first exercise; the standard exists because they did.

| # | Failure mode | Where caught | Guard now |
| --- | --- | --- | --- |
| 1 | the sibling-folder ancestor path carrying an instruction to read the quarantined tree | before launch, by the § 3.1 check | § 3.1, run the check, read the content |
| 2 | the launch shell arriving with the platform importable and the repo named in the prompt | live at launch | § 3.2 verify commands |
| 3 | the user-global instruction file summarized as "nothing relevant" while naming the repository | block A re-read | § 3.3, disclose and read the summary skeptically |
| 4 | prepared artifacts (task file, audit script, prompt) living only in session state | planning review | § 4, commit before the room opens |
| 5 | the audit's leakage check silently skipped when an earlier format check aborted | the audit's own mutation suite | § 4, leakage scan computed first, appended unconditionally |
| 6 | a 15-minute "stall" that was an unsurfaced permission prompt | during the run | § 6 stall row |
| 7 | the implementation printing the room's absolute path into committed output | copy-out | § 4 task-file small print; § 7 redaction rule |
| 8 | the byte-identity double run writing to session temp space against an absolute "writes only inside the room" statement | transcript check | § 7, disclosed session-temp writes |
| 9 | the audit reply landed as a merged-PR body edit, notifying nobody | close-out | § 9 communications rule |

None of the nine reached the result; five were caught by machinery (checks, mutation
suites, the transcript audit), four by a person reading carefully. That ratio is the
argument for keeping both.

## 14. Checklist

| Phase | Gate |
| --- | --- |
| build | ancestor-path check run and contents read; interpreter check green; unavoidable loads grepped and disclosed; room outside every working tree |
| packet | per-part mechanical audits green and mutation-tested; author evidence outside; incoming + authoritative hashes recorded; withholdings stated; task file carries the § 4 small print |
| commit 1 | packet + task file + audit script + output merged BEFORE launch |
| launch | hashes re-verified; banner cwd + model recorded; default permission mode; web tools withheld; near-empty opening prompt, nothing else |
| run | operator runbook applied per event; deviations logged as they happen; packet never amended |
| commit 2 | the commitment: hashes, environment, manifest, transcript check, redactions, pre-declarations, operator log; merged BEFORE unsealing |
| commit 3 | adjudication naming commit 2; own transcription; two-sided comparison mutations; exact comparison |
| close | method note draft frozen + maintainer layer; independent adversarial audit recorded; questions to humans in NEW comments |
| teardown | per-run write-up done, then delete the room; keep the transcript |
