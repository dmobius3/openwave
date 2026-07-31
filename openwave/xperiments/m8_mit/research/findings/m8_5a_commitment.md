# M8.5-A: the § 9 commitment

> The commitment record required by
> [`m8_5a_reproduction_protocol.md`](m8_5a_reproduction_protocol.md) § 9: the clean-room
> run's environment record, consulted-files manifest, commitment hashes, schema version, and
> the pre-declared § 6 outcome, committed BEFORE any quarantined object is unsealed for
> comparison. Adjudication against the pinned § 6.1 table is deliberately absent here: it is
> the § 7 step, in a second, separately dated commit that names this one. Task record:
> [`../tasks/m8_5_task_details.md`](../tasks/m8_5_task_details.md); procedure:
> [`dev_docs/tasks/t4_task_details.md`](../../../../../dev_docs/tasks/t4_task_details.md).

## 1. Commitment hashes (SHA-256)

Computed independently by the maintainer at copy-out and matching the implementer's own
`environment.md` record exactly, except where a redaction is disclosed in § 5.

| Artifact | Repo file | SHA-256 |
| --- | --- | --- |
| the implementation | [`../scripts/m8_5a_reproduction.py`](../scripts/m8_5a_reproduction.py) | `3433254983fb0d70ffbaafbdf8f669f2752521f9e273a22e0a0e1ec0ac130ee2` |
| raw output, as produced in the room | (room file; superseded in-repo by the redacted copy below) | `feef1e56604d0424907e7b7bf1fc1e03d0dec02c631fe573ae887806a041b3ff` |
| raw output, committed copy (one line redacted, § 5) | [`../data/m8_5a_raw_output.txt`](../data/m8_5a_raw_output.txt) | `b9a6323b0b9e22c35f6d4308ec973431e16df83671893b92c432a9ef486ca51f` |
| the result | [`../data/m8_5a_result.json`](../data/m8_5a_result.json) | `47c74573cdacecea868d74e21c8531014296c881fd5034ae5af9c15a5284fd02` |
| the method-note draft | [`m8_5a_method_note.md`](m8_5a_method_note.md) | `c8ac98a218c0449ff0dfa28eb280b879d9ce44c742c94cf13e7e635ea60ce6a8` |
| the input packet (unchanged through the run) | [`../data/m8_5a_packet.json`](../data/m8_5a_packet.json) | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` |

Schema version: `m8_5a-v1`.

**The § 6 declaration, fixed here and unavailable after unsealing: the coexact module RAN.**
Implementer-reported verdict: `structurally derived and reproduced`, with the derivation and
its self-flagged weakest links in the method-note draft, for the adversarial audit to attack.
Per the standing rule, the numerical match itself upgrades nothing.

## 2. Environment record (implementer's, verbatim content)

| Field | Value |
| --- | --- |
| Interpreter | CPython 3.13.9 (`python3`), conda base environment |
| Libraries | Python standard library only: `fractions`, `json`, `hashlib`, `argparse`, `platform`, `sys`, `os`, `copy`. No numpy, no third-party packages |
| OS / hardware | macOS, Darwin kernel 25.5.0, arm64 (Apple silicon) |
| Seeds | none: the pipeline is deterministic and fully exact over `Q(φ)`; a double-run byte-identity check confirmed determinism |
| Implementer | Claude Code session, model Fable 5, high effort, recorded from the launch banner by the operator |
| Date of run | 2026-07-31 |

## 3. Consulted-files manifest (implementer's, folded)

Deliberately consulted, all inside the room: `PROTOCOL.md`, `TASK.md`, `GROUP_INPUT.md` (read
in full), and `m8_5a_packet.json` (parsed, hash-verified against the `GROUP_INPUT.md` pin).
The implementer reports: no other file inside or outside the room read, listed, or searched;
no web access (disabled by design); the relative links inside `PROTOCOL.md` not followed.

Loaded automatically without being asked for, disclosed by the implementer:

| Item | Implementer's assessment |
| --- | --- |
| the user-global Claude instruction file | writing-style and workflow preferences; mentions the platform and correspondent names; no `2I` group data, characters, dimensions, distances, first occurrences, or M8 target values |
| the harness system prompt | tool definitions, environment description, memory-directory path; no memory files loaded, none carrying M8 facts exist for a fresh directory |

**Prior-knowledge disclosure (the § 2 ceiling), implementer's statement:** the implementer is
an AI model with an opaque training corpus, so prior exposure to published `2I` data cannot be
excluded; general facts about binary polyhedral groups and McKay graphs are part of general
mathematical training; no reference table was looked up or reconstructed from memory, and
every reported number is derived by the committed source from the packet generators, bound by
the gates and mutation harness. The claim label is capped at context-isolated
independent-method reproduction accordingly.

## 4. The maintainer's transcript check (the manifest is corroborated, not attested)

The session transcript was audited against the manifest before this commitment, per the T4
block C rule. Method: every `tool_use` event in the transcript JSONL was extracted and its
paths classified.

| Checked | Found |
| --- | --- |
| tool calls, total | 19: 3 reads, 4 writes, 2 edits, 10 shell commands |
| web tool calls | 0 |
| file-tool paths outside the room | none |
| shell commands touching paths outside the room | none reading; two writing, disclosed next row |
| the one deviation from "writes only inside the room" | the double-run byte-identity check redirected its two run outputs into the session's own harness-provided scratchpad directory and diffed them there. Session-private temp space, written not read, containing only the run's own stdout. Recorded because the write-scope rule was stated absolutely |
| permission prompts declined by the operator | none; every prompt was an inside-room read, write, or run |
| questions from the implementer mid-run | none |

## 5. Redaction disclosure

Line 62 of the room's `raw_output.txt` printed the clean-room directory's absolute path (the
implementation echoing where it wrote `result.json`). The T4 rule keeps that path out of the
repository, so the committed copy replaces the path with `<clean-room-directory>`; nothing
else differs. Both hashes are in § 1: the original hash stays the commitment (it is also
recorded inside the implementer's own `environment.md`), and the transcript retains the
original line for any later verification. Lesson for the T4 standard: instruct the
implementation to print repo-safe relative paths, so the raw output commits byte-identical.

## 6. Operator log (block B)

| Event | Note |
| --- | --- |
| launch | new terminal, room hash-verified (4 files), conda base environment verified free of the platform (`import` fails, `pip show` empty) after a first shell arrived with an environment carrying it (T4 deviations log) |
| opening prompt | the near-empty prompt, verbatim from the T4 record, nothing else typed |
| stall, ~15 min | a pending permission prompt had not surfaced in the terminal; the token counter sat frozen until the prompt appeared. Lesson: check for a pending prompt before interrupting a "stalled" room |
| prompts | all inside-room; each approved individually, no blanket grants used |
| finish | final message relayed to the maintainer session verbatim; room untouched until this copy-out |
