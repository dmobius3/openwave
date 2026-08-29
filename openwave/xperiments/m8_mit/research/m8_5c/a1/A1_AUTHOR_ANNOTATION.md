# Author annotation on attempt M85C-A1 (dated 2026-08-29, post-stop)

Author-side, written AFTER the Build Unit stopped and disconnected. The unit's files are
untouched: this annotation sits beside the ledger, never in it. The room is CLOSED; the
attempt does not resume under any disposition.

## Archive hashes, as the unit left them

| file | SHA-256 |
| --- | --- |
| `ledger/OUTPUT_LEDGER.jsonl` | `794e2969f2ca3afe35a897ed22bd5102fe7e48e5951eecbb563f7b0b846491b0` |
| `ledger/COMMISSIONING.md` | `498d3e000f1e3c826987b2eac5ed945d38c88032ffc5d77695e2d0ae0767f1c3` |
| `ledger/INPUT_MANIFEST.json` | `f4a3ddf8415c828a582fb5fcf47c3ca424fa53f29c4a635d7789bb0630f57d9a` |

## Two ledger-semantics corrections a later reader needs

1. **`gates_completed` in the STOP record lists `G4-PROJECTOR`. Read it as PARENT-ONLY.**
   Gate 4's parent check executed and measured green (4N rule exact to 7.7e-15); its arm
   set did NOT complete: the node-drop arm was found unexecutable (the STOP's own
   `dead_mutations` field) and the second arm, the section 4.3 dual-route agreement, never
   ran. This is the one field a mechanical adjudicator would key on, and the disposition
   question turns exactly on complete-and-failed versus incomplete, so the correction is
   recorded here rather than left to inference.
2. **`cumulative_seconds` in the RESOURCE records carries the PER-GATE wall clock**, not a
   running total. The STOP record's `cumulative_wall_seconds` of 333.4 is the true
   cumulative and equals the per-gate sum.

## The denominator, stated plainly

Of the protocol's ten checks: three completed with firing mutation arms (gates 1 to 3);
gate 4 is PARENT-ONLY, its parent green at 7.7e-15 and its mutation arm unexecutable, the
second arm never run; six (gates 5 to 10) were never reached. No chassis prediction was
falsified, and this was NOT nearly a qualification.

## Provenance timing of the post-stop measurements

The replacement-arm measurements (angular drop to even K fires at 1.8e-01; halved
Gauss-Legendre u-nodes fire at 7.1e-02; full rule exact as control) and the three-way
parity discrimination (even-only exact, odd-only exact, mixed fires at 9.8e-02) were taken
by the author AFTER the stop, while the disposition was open. They are successor design
inputs and A1 evidence about the mechanism; they are not part of the unit's run and carry
no execution credit.

## The parity-purity strengthening (redline unit 2; verified by the author against the unit's built fibres AFTER the stop, same post-stop provenance class as the replacement-arm measurements)

Since `-1` is central and acts on `V_n` as `(-1)^n`, a sector occurs only at levels with
`(-1)^n` equal to its central sign: every sector's level set is parity-pure at
`n = d_rho mod 2`, confirmed for all nine sectors against the room's built fibre
characters over n = 0..24. The node-drop arm is therefore dead on ALL NINE arenas at every
rung, and the design record's mixed-parity field set (levels 0 through 3) is a
configuration that cannot arise in any legal arena of this protocol.
