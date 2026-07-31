# M8.8: Independent-method reproduction of the M8.3 torsion closed forms

> Roadmap row: [`../m8_roadmap.md`](../m8_roadmap.md). Status: 🚧 PLANNED STUB
> (2026-07-30). Owner unassigned. Full PLAN at go.

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

**Gated by**: M8.3 ✅ (2026-07-28) and an owner. The sign convention is fixed once on the `R7`
closed form and declared, so the remaining forms are genuine checks rather than circular
([`m8_3_method_note.md § 4`](../findings/m8_3_method_note.md)).

## DEVIATIONS LOG

(none)

## FINDINGS

(pending)
