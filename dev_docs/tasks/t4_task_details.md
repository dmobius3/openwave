# T4: The context-isolated clean room

> Roadmap row: [`../platform_roadmap.md`](../platform_roadmap.md). Status: 🚧 PLANNED
> (filed 2026-07-30). Owner: maintainers.

## PLANNING

### Why this is platform work

The [M8.5-A protocol](../../openwave/xperiments/m8_mit/research/findings/m8_5a_reproduction_protocol.md)
§ 4 specifies a context firewall but not how to build one on a real machine. Answering that
produced a set of decisions that are not about `2I`, the McKay ladder, or anything else in the
M8 column: where the room lives, which tools are withheld, what enters the packet, who audits
it, and what order the commits go in. Every one of those recurs the next time any column owes a
context-isolated reproduction.

The [roadmap's ownership test](../platform_roadmap.md#conventions) splits it cleanly. The
procedure is shared machinery and files here. The M8.5-A result is evidence about one model's
physics and files on the [M8 roadmap](../../openwave/xperiments/m8_mit/research/m8_roadmap.md)
under M8.5-A, where its planning already sits in
[`m8_5_task_details.md`](../../openwave/xperiments/m8_mit/research/tasks/m8_5_task_details.md).

**Deliverable**: a `dev_docs/` standard, written from an exercise rather than in the abstract,
in the shape of [`METHOD_NOTE.md`](../METHOD_NOTE.md) and
[`CROSS_MODEL_TESTING.md`](../CROSS_MODEL_TESTING.md). It is written after the first run, not
before, so it records what actually held.

### The room

| Setting | Decision |
| --- | --- |
| Location | any directory with no `CLAUDE.md` anywhere on its ancestor path, and outside every repository working tree. The concrete path is a local choice and stays out of the repository |
| Launch | `claude --disallowedTools "WebSearch,WebFetch"` |
| Write scope | the implementer writes only inside the room and never reaches the repository; a maintainer copies artifacts out and performs every commit |
| Teardown | the room is deleted at the end of the run, block F below. Not before: a rerun under § 3 needs it, and it is the only copy until block C has copied the artifacts out |

**The ancestor-path rule is the load-bearing one, and it is easy to get wrong.** A folder
alongside the repository looks isolated and is not: any `CLAUDE.md` above it loads on the first
turn, before the agent reads a word of the packet. The obvious sibling location for this run
sits under a 7-line `CLAUDE.md` whose entire content instructs the agent to go read the
`openwave` directory, which is the quarantined tree. Nothing in the room's own contents would
have revealed the breach. Check the whole ancestor path, and check what the file says rather
than only whether one exists.

**What loads unavoidably**: `~/.claude/CLAUDE.md`, the user-level instructions. Grepped for the
subject terms before the run and found nothing relevant, so it enters the consulted-files
manifest as a disclosed load rather than an unmentioned one. Project memory is keyed to the
working directory, so a new folder starts with none.

### The packet, and who may write it

§ 4 defines the packet as the protocol plus the group input. A maintainer who already knows the
answers is the one assembling it, which is exactly the situation the § 4 audit anticipates. The
containment is that no part of the packet is capable of carrying a target:

| Part | Audit |
| --- | --- |
| `PROTOCOL.md` | a verbatim copy of the frozen protocol, so `diff` against the repository file is the whole audit and no authorship enters |
| `GROUP_INPUT.md` | two generators as unit quaternions, supplied by the model author. The audit is mechanical: closure gives exactly the expected order, and the file carries no labels, dimensions, distances, or character values |
| `TASK.md` | operational instructions only: the added gates, the § 6 declaration, output paths. Its mathematical content is generic identities, not values |

**Generators, not the full element set.** G1 gates group order and closure. Handing over the
complete element set makes G1 true by construction, which is precisely the class of check the
column banned after `m7_trivial_ok`. Handing two generators makes the agent build the closure
and G1 becomes a check that can fail.

**Canonicalize before hashing.** The packet hash is what preserves provenance across the run,
and it is only meaningful if the same group always produces the same bytes. Generator
components are generally irrational, so a supplied decimal rendering fixes a hash to a
transcription rather than to a group. Canonicalization (a declared exact form, or a declared
precision and ordering) happens in block A before the hash is taken, and the canonical form
is what enters the room.

**Implementer eligibility.** Any context that has read the M8.2 generator, the maintainer
reconstruction, or the pre-registration § 6.1 is permanently disqualified, whatever it later
forgets. That includes the session that writes the packet.

### Running it beside a live maintainer session

The two sessions are independent processes and may run side by side. The information rule is
one-directional:

| Direction | Rule |
| --- | --- |
| clean room → maintainer session | unrestricted. The maintainer session is already disqualified, so nothing is lost by showing it anything |
| maintainer session → clean room | ONLY the audited packet, fixed before launch and not amended mid-run |

If the implementer asks a question mid-run, it is answered from the protocol text or not at
all, and the exchange is recorded as a deviation. Relaying a sentence composed by a context
holding the targets defeats the firewall no matter how innocuous the sentence looks.

### Run plan

| Block | Work | Commit boundary |
| --- | --- | --- |
| A | receive the author's generators; audit the packet (closure to the expected order, no derived annotations, `diff` on the protocol copy), canonicalize it, record its SHA-256 | the packet and its hash, committed first so the timestamp predates the run |
| B | the clean-room session: implementation, gates, runnable mutation harness, the § 6 module, raw output and JSON | none; the agent writes only inside the room |
| C | copy artifacts out, generate SHA-256, write the environment record and method-note draft | **the § 9 commitment.** Nothing is unsealed before it lands |
| D | unseal; build the § 7 comparison harness with its transcription mutation; adjudicate | second, separately dated commit |
| E | method note per [`METHOD_NOTE.md`](../METHOD_NOTE.md), then the adversarial audit by a second agent | final |
| F | **delete the room.** Keep its session transcript | none |

Blocks A through D are one sitting. Block E is not promised in the same sitting: an adversarial
audit that is rushed to fit a day is not an adversarial audit.

**Block F is a required step, not housekeeping.** A room left on disk is a room that gets
reused, and its second occupant is not context-isolated: the folder holds the packet, the
implementation, and by then the unsealed comparison. Delete it once block E closes.

Keep the session transcript, which lives under the agent's own project-history directory rather
than inside the room, so deleting the room does not touch it. It is worth keeping because it
corroborates the § 4 manifest: the manifest states what the implementer consulted, and the
transcript is the independent record of what it actually opened. A manifest with no way to
check it is an attestation, and this column already ruled against those.

### Suggested definition of done

| # | Item |
| --- | --- |
| 1 | The `dev_docs/` standard exists and is written from the run, not from the plan |
| 2 | It names the failure modes the exercise actually hit, including any this doc did not predict |
| 3 | The ancestor-path check is stated as a procedure a reader can execute, not as advice |
| 4 | It states which parts are protocol-specific to M8.5-A and which any column inherits |

### Blindspots

| Risk | Guard |
| --- | --- |
| The standard is written from this plan rather than from the run, so it records intentions | block E writes it, after adjudication, with a deviations log kept throughout |
| A future room is placed by copying the path rather than applying the rule | the standard leads with the ancestor-path check and the sibling-folder near-miss that motivated it |
| The packet audit becomes a rubber stamp because the same person assembles and audits | each part has a mechanical audit (`diff`, closure count) rather than a judgment call |
| Withholding web tools is mistaken for an offline requirement | the standard states plainly that the network stays up and only the search tools are withheld |

### Ownership + gating

**Gated by**: the M8.5-A packet being ready. No author reply is required: the packet's group
input is forced by protocol § 5.3, and the fallback for unanswered provenance is recorded in
[`m8_5_task_details.md`](../../openwave/xperiments/m8_mit/research/tasks/m8_5_task_details.md).

## DEVIATIONS LOG

(none)

## FINDINGS

(pending)
