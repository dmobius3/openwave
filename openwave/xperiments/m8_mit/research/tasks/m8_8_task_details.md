# M8.8: Independent-method reproduction of the M8.3 torsion closed forms

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: 🚧 PLANNED STUB
> (2026-07-30). Protocol author: the model author (2026-08-02); implementer deliberately
> unassigned (fresh context). Full PLAN at go.

## Why this task exists

[M8.3](m8_3_task_details.md) ✅ landed the mass-formula reproducer and corrected a defect in
the published page. Its method note lists one claim it does not verify: the corrected
Reidemeister-torsion closed forms rest on a **single implementation**, since the script and the
mode-identity-theory artifact use closely related rather than disjoint methods
([`m8_3_method_note.md § 4`](../findings/m8_3_method_note.md)).

That obligation was recorded as "queued (M8.5)". The
[M8.5-A protocol § 0](../findings/m8_5a_reproduction_protocol.md) then placed it outside both
M8.5 sub-deliverables, which left the two documents pointing at each other with no task owning
the work. This row closes that loop; the scope statement below is the one both documents now
point at.

## Scope (stub level)

| Piece | Content |
| --- | --- |
| Target | the corrected torsion closed forms: the 9 irrep forms, both Galois ratios, both sector products |
| Method | disjoint from [`m8_3_mass_reproducer.py`](../scripts/m8_3_mass_reproducer.py) and from the mode-identity-theory artifact; the overlap is disclosed, not assumed absent |
| Not in scope | the mass comparison itself (M8.3 verified its arithmetic), and the dead-zone entries' physical status (open on the source page) |
| Precedent to follow | the M8.5-A shape: gates that can go red, a coverage-enforced mutation harness, an explicit list of what is not verified |

**The implementer must be a fresh context** (author direction, 2026-07-30). Whoever runs this
may not be the context that produced the M8.3 computation, on the same reasoning that shaped
M8.5-A: a context holding the closed forms and their derived fixtures cannot serve as its own
reproducer however separately the second implementation is written. The owner therefore stays
deliberately unassigned rather than defaulting to the nearest available context. Protocol
authorship can be assigned later without compromising the reproduction, and nothing on the
dynamics path waits on either.

## Protocol authorship and group input (2026-08-02)

The author took the protocol half of this task
([#405 comment](https://github.com/openwave-labs/openwave/pull/405#issuecomment-5159525049),
2026-08-02), on the M8.5-A terms: the answer-holding author writes the reproduction protocol,
the maintainer reviews and freezes it, and a fresh context implements under the resulting
firewall. Three maintainer calls were made with the go:

| Call | Decision |
| --- | --- |
| Protocol authorship | Go. The stub already provided for assigning it without compromising the reproduction |
| Group input | The audited M8.5-A raw generator packet ([`../data/m8_5a_packet.json`](../data/m8_5a_packet.json)) may serve, confirmed to contain no derived data relevant to the torsion target (raw `(a + b*phi)/2` quaternion generators, minimal polynomial, basis, format tag only). Condition: the M8.8 protocol names the packet by its published hash and declares its own forbidden-inputs list, and the packet audit is RE-RUN against that list at freeze, not inherited from the [M8.5-A audit record](../data/m8_5a_packet_audit.json) |
| Timing | Protocol drafting begins on the author's schedule; if its review ever lands together with the M8.5-B implementation review, B takes precedence |

The M8.5-A structural claim ceiling applies unchanged: the torsion closed forms are published
on the source page, so prior corpus exposure cannot be excluded for an AI implementer, and
isolation buys provenance rather than the label.

**Gated by**: M8.3 ✅ (2026-07-28) and an owner. The sign convention is fixed once on the `R7`
closed form and declared, so the remaining forms are genuine checks rather than circular
([`m8_3_method_note.md § 4`](../findings/m8_3_method_note.md)).

## DEVIATIONS LOG

(none)

## FINDINGS

(pending)
