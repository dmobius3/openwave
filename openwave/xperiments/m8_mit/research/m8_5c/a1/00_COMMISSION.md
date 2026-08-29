# M8.5-C Build Unit commission

You are the fresh-context **Build Unit** for M8.5-C, the ONE preregistered qualification
attempt for the spectral/Galerkin chassis. Target-free: you determine whether the chassis
is numerically trustworthy, without possessing any target. There is exactly one attempt.

**What governs.** The frozen region of
`openwave/xperiments/m8_mit/research/findings/m8_5c_protocol.md` governs every gate,
threshold, ladder, definition, arena, stop condition, and outcome sentence. The symmetry
derivation note beside it is the theorem record § 2 names. THIS COMMISSION IS
NON-GOVERNING (protocol § 12): it carries room mechanics only. If this commission and the
protocol ever disagree, the protocol wins; if the disagreement is material, STOP and
report it instead of choosing. If the protocol itself cannot be implemented AS WRITTEN at
any point BEFORE your first GATE record, STOP and report: a pre-execution defect
supersedes the protocol as a whole (§ 0, § 11), and that supersession decision belongs to
the author, never to you.

## 1. Provenance, before any other work

Your handoff prompt supplies two values that appear in no file in this room: the
ROOM_MANIFEST SHA-256 and your attempt id. First:

1. `shasum -a 256 ROOM_MANIFEST.json` must equal the handoff value. A mismatch means the
   room changed after commissioning: STOP.
2. Verify the protocol's frozen region with its own recorded command, run from
   `openwave/xperiments/m8_mit/research/findings/`:

       sed '/^<!-- M85C-FREEZE-BOUNDARY -->$/,$d' m8_5c_protocol.md | shasum -a 256

   It must equal the value in the protocol's own freeze record below its boundary, and
   the COMPLETE marker line must occur exactly once (`grep -c '^<!-- M85C-FREEZE-BOUNDARY -->$'`).
   Hash checks match the full marker line only; inline mentions are harmless (§ 12).
3. Verify every file listed in `ROOM_MANIFEST.json` against its recorded hash. The § 15
   freeze-time pins for these same files must agree; read § 15 and check.
4. Collision guard: `ledger/` and `build/` must be empty except for what this commission
   names. If a prior attempt's records exist, STOP without reading them.

## 2. Launch gate

    python3 room_import_gate.py --selftest      # must end SELFTEST GREEN
    python3 room_import_gate.py                 # must end GATE GREEN, exit 0

Both outputs go into the COMMISSIONING record verbatim, with your interpreter path and
numpy/scipy versions. The design-inputs README records the versions its evidence was
re-run under; record any difference. A RED gate means the room does not launch: STOP.

## 3. The ledger, before any build work

Create `ledger/COMMISSIONING.md` NOW and append as you go; never batch. Rooms have died
with their evidence living only in a terminal. Everything you establish must be on disk
and hashable at the moment you establish it.

    ledger/COMMISSIONING.md        preflight: provenance checks, gate outputs, environment,
                                   room-setup notes, and the INPUT_MANIFEST closing hash
    ledger/INPUT_MANIFEST.json     the § 12 input manifest, built by you (below)
    ledger/OUTPUT_LEDGER.jsonl     append-only, one JSON record per event, § 13's REQUIRED
                                   record types with their named fields
    build/                         every line of code you write lives here

You write ONLY under `ledger/` and `build/`. The supplied inputs are immutable: any byte
change to a whitelisted file voids the room. Never write outside the room directory,
never read outside it (§ 12 read closure: your executable must enforce this on itself,
mutation-armed), no network, no repository access. When a supplied document's relative
link points at a file that is not here (`m8_4_closeout.md`, briefing, roadmap), that is
the whitelist working as designed, not an error: those documents are named OUT (§ 12) and
you must not seek them elsewhere. The shipped M8.2 tables are the ADJUDICATOR's input,
never yours: your § 3.4(e) first-occurrence structure is RE-DERIVED in-room, and the
comparison happens outside this room.

## 4. Sequence to the execution boundary

Execution begins at your FIRST GATE RECORD in `ledger/OUTPUT_LEDGER.jsonl` (§ 11). Before
it, in order:

1. Preflight (§§ 1 and 2 above) into COMMISSIONING.md.
2. Build the § 3 basis objects per the deterministic construction (stacked-constraint SVD
   in the frozen element order, Lowdin, sign rule), the § 10 arena-constructor registry,
   and the § 4.3 packet; write the packet's canonical-serialization SHA-256 digest.
3. Close `ledger/INPUT_MANIFEST.json`: the whitelist table from ROOM_MANIFEST.json
   verbatim, the arena-constructor registry, the § 4.3 digest, the § 3 basis-object hash
   per sector and rung, the Control A reference and lattice-table hashes, and this
   commission's own hash. Hash the closed manifest and write that hash into
   COMMISSIONING.md. After this point the input manifest never changes.
4. Only then: the first GATE record. From it, one attempt means one attempt (§ 11).
   Gates, ladders, arms, dispositions, and stop conditions are §§ 5 through 11 of the
   protocol; report wall-clock per gate in RESOURCE records against the § 11 ceiling.

## 5. What ships

Per § 13: the output ledger, the logs, the input manifest, and every hash. You write no
findings prose, no verdict, and no § 1 sentence: adjudication is a separate unit that
receives the protocol, your ledger, and the locked M8.2 reference, nothing else. Your
COMMISSIONING record and this commission are archived as non-governing provenance.

## 6. Discipline reminders, earned the hard way in this column

- Append to the ledger after every completed check, control, and mutation. Do not batch.
- A verdict your code prints must be COMPUTED, and your process exit codes must carry it
  (the platform's CONTRIBUTING records this as the fourth unfalsifiable-check shape).
- Never redirect a rerun over an existing record's filename: a v2 run gets a v2 filename;
  the v1 record is evidence and is immutable.
- Every PASS needs an arm that can fail, run in the same session, with its green parent
  named in the same record (§ 8's table declares them; § 13's GATE record carries both).
- If a check cannot be implemented as written, STOP and report. Never substitute a weaker
  check silently. Attribution classes (instrument vs bases vs run) are the protocol's,
  and the ledger records which one fired.
