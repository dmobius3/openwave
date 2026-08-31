# M8.5-C2 Build Unit commission

You are the fresh-context **Build Unit** for M8.5-C2, the ONE preregistered qualification
attempt for the spectral chassis under the successor protocol. Target-free: you determine
whether the chassis is numerically trustworthy, without possessing any target.

**COMMISSIONING DOES NOT SPEND THE ATTEMPT; YOUR FIRST GATE LEDGER RECORD DOES.** Every
preflight step, the § 3 build, the § 4.3 packet digest, and the input-manifest close all
happen BEFORE that boundary, and a protocol defect found before it supersedes the protocol
as a whole at no cost to the attempt (§§ 0, 11). After the first GATE record, one attempt
means one attempt, and § 11's PROTOCOL-INVALID class is the only administrative exit.

**What governs.** The frozen region of
`openwave/xperiments/m8_mit/research/findings/m8_5c2_protocol.md` governs every gate,
threshold, ladder, definition, arena, stop condition, and outcome sentence. The symmetry
derivation note beside it is the theorem record § 2 names. THIS COMMISSION IS
NON-GOVERNING (§ 12): room mechanics only; on any material conflict with the protocol,
STOP and report rather than choosing.

## 1. The bootstrap sequence, frozen so it is executable exactly as written

Your handoff prompt supplies two values that appear in no file in this room: the
ROOM_MANIFEST SHA-256 and your attempt id (§ 13 namespaces it; the bare "A1" in the
protocol always means the SUPERSEDED protocol's archived attempt, never yours).

The collision guard requires `ledger/` empty, and the append-as-you-go rule requires the
commissioning record to exist from the start; those two cannot both be literal, so the
order below IS the resolution and is not yours to rearrange:

1. Verify `shasum -a 256 ROOM_MANIFEST.json` equals the handoff value; a mismatch means
   the room changed after commissioning: STOP, writing nothing.
2. Collision guard: `ledger/` and `build/` must be EMPTY. Prior records anywhere: STOP
   without reading them, writing nothing.
3. Immediately create `ledger/COMMISSIONING.md` and write steps 1 and 2's results as its
   first entry. These two facts, and only these two, are necessarily recorded together,
   because the guard itself requires the record not yet to exist; everything after this
   line appends the moment it is established, never batched.
4. Verify the protocol's frozen region with its own recorded command, run from
   `openwave/xperiments/m8_mit/research/findings/`:

       sed '/^<!-- M85C2-FREEZE-BOUNDARY -->$/,$d' m8_5c2_protocol.md | shasum -a 256

   It must equal the digest in the freeze record below the boundary, and the COMPLETE
   marker line must occur exactly once (full-line matching; inline mentions are
   harmless). Append the result.
5. Verify every file listed in `ROOM_MANIFEST.json` against its recorded hash, AND the
   § 15 freeze-time rows for those same files against the same values; append each
   result as computed.
6. Run `python3 room_import_gate.py --selftest`; it must end SELFTEST GREEN with all
   three arms firing. Append the output verbatim.
7. Run `python3 room_import_gate.py`; it must end GATE GREEN, exit 0. Append the output
   verbatim, with interpreter and numpy/scipy versions and any difference from the
   design-inputs README's recorded versions. RED at step 6 or 7 means the room does not
   launch: STOP.

The gate's selftest writes and removes scratch under `build/selftest_gate/`, so `build/`
being non-empty AFTER step 6 begins is gate residue, not a prior record; the guard at
step 2 is the only point that reads emptiness.

## 2. The ledger, from here forward

`ledger/COMMISSIONING.md` exists from step 3; append as you go, never batch. Everything
you establish must be on disk and hashable the moment you establish it.

    ledger/COMMISSIONING.md        preflight record, and the INPUT_MANIFEST closing hash
    ledger/INPUT_MANIFEST.json     the § 12 input manifest, built by you (below)
    ledger/OUTPUT_LEDGER.jsonl     append-only, one JSON record per event, § 13's
                                   REQUIRED record types with their named fields
    build/                         every line of code you write

**§ 13's ledger schema binds you and the adjudicator checks it mechanically:** every
RESOURCE record carries BOTH `wall_clock_seconds` and a strictly nondecreasing
`cumulative_seconds` whose sum identity the adjudicator recomputes; no record carries a
bare "completed" list, only per-gate states from {COMPLETE, PARENT-ONLY, NOT-REACHED}
where COMPLETE means parent AND every declared arm; and EVERY hash you record anywhere
NAMES THE BYTES IT COVERS (path plus the event after which it was computed, or the
recipe), a bare hash being a missing-record defect under § 13's own clause.

You write ONLY under `ledger/` and `build/`; the supplied inputs are immutable and any
byte change voids the room. No reading outside the room, no network, no repository. The
SUPERSEDED M8.5-C protocol, its addenda, and the A1 archive `research/m8_5c/a1/` are
named OUT by § 12: the room re-earns every green from gate 1, and a dangling relative
link to an OUT document is the whitelist working, not an error. The shipped M8.2 tables
are the ADJUDICATOR's input, never yours.

## 3. Sequence to the execution boundary

1. The § 1 bootstrap, already complete by this point, lives in COMMISSIONING.md.
2. Build the § 3 basis objects per the deterministic construction, the § 10
   arena-constructor registry, and the § 4.3 packet; write the packet's
   canonical-serialization digest. § 8's arena column names each arm's draw space, and
   read space where distinct; your implementations must match those cells, and gate 6's
   leakage read is RUNG-RELATIVE, `n ∈ {N+1, N+3, N+5}`, never literals.
3. Close `ledger/INPUT_MANIFEST.json`: the whitelist table from ROOM_MANIFEST.json
   verbatim, the arena registry, the § 4.3 digest, the § 3 basis-object hash per sector
   and rung, the Control A reference and lattice-table hashes, this commission's hash,
   every entry carrying its § 13 coverage statement. Hash the closed manifest into
   COMMISSIONING.md. It never changes afterward.
4. Only then: the first GATE record. From it, the attempt is live and §§ 5 through 11
   own everything, including STOP-QUAL and the ceiling's RESOURCE accounting.

## 4. What ships

Per § 13: the output ledger, the logs, the input manifest, and every hash. You write no
findings prose and no verdict; adjudication is a separate unit receiving the protocol,
your ledger, and the locked M8.2 reference, nothing else.

## 5. Discipline, earned the hard way in this lineage

- Append after every completed check, control, and mutation.
- A verdict your code prints must be COMPUTED and carried into the exit code.
- Never redirect a rerun over an existing record's filename; v2 runs get v2 filenames.
- Every arm fires against a named green parent in the same session; regression CONTROLS
  are named as controls, their green condition machine zero, never counted as arms.
- Prefer repo-relative paths from the room root in everything you write; no absolute
  paths in code (the captured stdout of a run is a record and keeps whatever it printed).
- If anything cannot be implemented AS WRITTEN, STOP and report. The gate-6 leakage
  finding that created this protocol was caught exactly that way, pre-freeze, by its own
  obligation; in-room, § 11's PROTOCOL-INVALID path exists for a demonstrated
  contradiction, and its bar is a maintainer reproduction with independent code.
