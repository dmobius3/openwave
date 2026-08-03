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
| Target | the corrected torsion closed forms: the irrep forms, both Galois ratios, both sector products. `T^2(R0) = 1` is carried as a declared convention rather than a gated form (see the 2026-08-03 section), matching [`m8_3_method_note.md § 4`](../findings/m8_3_method_note.md) |
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

## Construction input: the second packet (2026-08-03)

Pre-draft recon raised a second input-boundary question
([#402 comment](https://github.com/openwave-labs/openwave/pull/402#issuecomment-5166892675),
2026-08-03). The M8.5-A packet pins the group, but a combinatorial torsion calculation also needs
the 3-cell: a balanced presentation gives a 2-complex with `chi = 1` while the closed orientable
quotient has `chi = 0`, so Fox calculus supplies `d1` and `d2` and the torsion product still
consumes `d3`, which requires the identity among relations or an equivalent resolution.

| Call | Decision |
| --- | --- |
| Construction packet | Go. A second public packet carrying a based cellular chain complex or a truncated periodic `Z[2I]`-resolution, boundary maps in abstract generators explicitly matched to the M8.5-A quaternion generators, carrying no evaluated irreducible matrices, determinants, torsion values, ratios, products, or target forms |
| Reference-allowance alternative | Declined, on the precedent that already declined a web allowlist for M8.5-A: the literature carrying the period-4 resolution for this group carries the torsion formulas with it, so an allowance re-opens the leak the firewall exists to close |
| Packet audit | Maintainer-side, independent, and mechanical, run from the packet alone, per the [M8.5-A § 4 rule](../findings/m8_5a_reproduction_protocol.md) that the author knows the answers and the audit is the guard against author-side leakage. Any author-side verification artifacts ship in a separate archive that stays outside the room and is set side by side afterwards, as in [M8.5-A](m8_5_task_details.md) |
| Packet roles | The construction packet is public before the run because it is permitted input; the canonical answer packet stays quarantined until commitment because it is only an adjudication reference |
| Claim | Narrows as the author states: an independent reproduction of the torsion calculation, not an independent derivation of the supplied model. The supplied complex is verified rather than trusted, since `d d = 0`, the rank and `chi` census, the integral homology of a `Z`-homology 3-sphere, and per-irrep acyclicity are all gates that can go red |

**Why the disjointness is stronger than it looks.** M8.3 computed the torsions analytically, from
the spectral-zeta definition (`log T^2 = zeta'_coexact(0) - 2 zeta'_scalar(0)`, the Ray-Singer
combination on `S^3/2I`), not from any chain complex. Nothing in a construction packet was an
input to M8.3, and the equality of the analytic and combinatorial routes is the Cheeger-Müller
theorem rather than a shared construction, so a disagreement localizes an implementation error in
one route instead of being unresolvable.

Two protocol points settled before drafting, to keep them out of the freeze review:

| Point | Content |
| --- | --- |
| Normalization anchor | Combinatorial torsion is pinned only once the basing is pinned. The packet supplies that, but the correspondence to M8.3's analytic normalization still needs one declared anchor. M8.3 fixed the overall sign once on the `R7` closed form and declared it, leaving the rest genuine checks; the protocol declares its anchor the same way and before the run |
| The trivial representation | `T^2(R0) = 1` is a declared convention in M8.3, not a gated result, because the twisted complex is not acyclic for the trivial representation. The same obstruction appears in the combinatorial route, so `R0` is carried as convention there too rather than counted as a reproduced form |

**Gated by**: M8.3 ✅ (2026-07-28) and an owner. The sign convention is fixed once on the `R7`
closed form and declared, so the remaining forms are genuine checks rather than circular
([`m8_3_method_note.md § 4`](../findings/m8_3_method_note.md)).

## DEVIATIONS LOG

(none)

## FINDINGS

(pending)
