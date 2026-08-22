# M8.8 Independent-Method Reproduction Protocol (torsion closed forms)

> **Status: DRAFT FOR FREEZE (author-written protocol; NO implementation code and NO
> reproduction run in this document).** Freezes WHAT must be reproduced, by WHAT class of
> method, under WHAT firewall, and what may be claimed at the end. It does not perform the
> reproduction. Pairs with the task record
> [`m8_8_task_details.md`](../tasks/m8_8_task_details.md) and inherits the procedure of
> [`dev_docs/CLEAN_ROOM_STANDARDS.md`](../../../../../dev_docs/CLEAN_ROOM_STANDARDS.md).
> Owner: protocol author, Blake Shatto; implementer deliberately unassigned (fresh context);
> reviewer: maintainers.

## 0. Scope

M8.8 discharges one obligation recorded in
[`m8_3_method_note.md § 4`](m8_3_method_note.md): the corrected Reidemeister-torsion closed
forms rest on a **single implementation**, because the M8.3 script and the
mode-identity-theory artifact use closely related rather than disjoint methods.

**The reproduction target** is the corrected torsion table on `S³/2I`: nine irrep closed
forms, both Galois ratios, both sector products. § 5 fixes it exactly and separates the
output surface from the evidentiary content, which are not the same thing.

**Out of scope, stated explicitly:**

1. **The mass comparison.** M8.3 verified its arithmetic. Nothing here re-derives the
   `(ρ, σ) → fermion` assignment, the structural sectors, or the PDG comparison.
2. **The dead-zone entries' physical status.** Open on the source page and outside any
   reproduction.
3. **Everything M8.5 owns.** The M8.5-A protocol § 0 placed this task outside both M8.5
   sub-deliverables; that boundary is unchanged.
4. **Independent derivation of the topological model.** The based complex is supplied and
   verified, not discovered here. § 6 and § 9 state this in those words.

## 1. The principle, and the real risk

M8.5-A's risk was context leakage: a verifying context that had been handed the answers.
M8.8 carries that risk and two more that matter as much.

**Risk 1, method non-disjointness.** M8.3 computed the torsions ANALYTICALLY, from the
spectral-zeta definition

```text
log T^2 = zeta'_coexact(0) - 2 zeta'_scalar(0)      (Ray-Singer, S^3/2I)
```

over twisted 0-form and coexact 1-form spectra, with a quasi-quadratic fit and a Hurwitz-zeta
reduction. The mode-identity-theory artifact is closely related. Both existing
implementations therefore sit on the same side, and a second implementation that reaches for
the same definition reproduces the arithmetic without testing the result.

**Risk 2, recall.** Nine closed forms in `φ` are a compact, published object. An implementer
can produce the right table without deriving it, and a protocol that gates only the final
answer cannot tell the two apart. § 7 is the response.

**What makes the reproduction worth doing.** The two routes are joined by the
Cheeger-Müller theorem, not by shared computational machinery. M8.3 supplied no
chain-complex fixture, and the construction packet inherits nothing from its analytic
computation: that object never existed on the analytic side. The target equality is
therefore a theorem, and a disagreement is assessed against a STRUCTURED DIAGNOSTIC
PARTITION rather than being unresolvable: every preregistered branch not excluded by its
evidence is reported, and no unique attribution is forced. Because
the packet is nevertheless AUTHOR-SUPPLIED, a model or provenance error remains a distinct
branch below, bounded by the independent packet audit and the § 9 model gates rather than
excluded by them.

**The disagreement taxonomy, and why all four branches stay open.** A mismatch is one of:

| Branch | What it means | What is meant to exclude it |
| --- | --- | --- |
| implementation error in one route | a bug on the analytic or the combinatorial side | the existing M8.3 verification record (pinned in § 11) for the analytic side; the § 9 gates for this run |
| convention-bridge mismatch | the outputs differ by the SOLE admitted bridge: the global sign of `log T²`, equivalently the inversion `T² ↔ (T²)⁻¹` | the anchor rule preregistered in § 5.4 |
| failed hypothesis of the theorem | Cheeger-Müller's classical statement requires an ACYCLIC representation | the per-irrep acyclicity gate of § 9 |
| supplied-model or provenance error | the public construction packet does not represent the intended based model, despite internal consistency | the § 4.2 frozen-scope audit, the maintainer-held provenance record, the generator correspondence, and the § 9 model gates |

The third branch is not hypothetical. The acyclicity hypothesis demonstrably fails for the
trivial representation, which is exactly why `T²(R0) = 1` is a declared convention rather
than a computed value (§ 5.2). The acyclicity gate is therefore the theorem's own hypothesis
check and not merely a model check, and no branch may be dismissed until its gate has passed.
The fourth branch is bounded, never excluded: the model gates establish that the supplied
complex behaves like the based complex of an integral homology 3-sphere with the right
homology; they do not by themselves identify the model, which is why the audit and the
maintainer-held provenance record sit beside them.

## 2. Claim ceiling and label

**The ceiling is context-isolated independent-method reproduction, and never blind.** The
corrected closed forms are published on the source page. Prior corpus exposure cannot be
excluded for an AI implementer, so isolation buys auditable provenance rather than ignorance.
This is the M8.5-A ceiling, inherited unchanged, and it is structural: no procedure available
here can raise it.

Frozen phrasings, so no write-up improvises:

- A successful run, which is § 8 outcome `reproduced`, or `convention difference` with the
  orientation recorded:
  > **The M8.3 torsion closed forms were reproduced by a context-isolated
  > independent-method run, from a based chain complex rather than the spectral-zeta
  > definition. The supplied topological model was verified, not independently derived.**
- Every other outcome keeps its own category (§ 8) and is never absorbed by the label above.

**`blind` is never used as a RESULT LABEL in any M8.8 deliverable.** It names a rung this
procedure cannot reach, and it appears in this protocol only in order to be excluded.

## 3. The three objects, the quarantine, and the ordering record

Inherited from M8.5-A § 3 unchanged in shape.

| Object | Status |
| --- | --- |
| 1, the M8.3 implementation and its outputs | QUARANTINED from the implementing context |
| 2, the mode-identity-theory artifact | QUARANTINED |
| 3, the canonical answer packet (§ 4.4) | QUARANTINED until commitment |

**The ordering record is the evidence the firewall held.** The implementing context commits
its source, its environment record, its route-native derivation artifacts (§ 7) and its raw
output BEFORE any quarantined object is unsealed to it. Commit timestamps and output hashes
are the durable record; an attestation is not.

## 4. Context firewall and construction inputs

### 4.1 Three packets, with deliberately different statuses

The distinction is the point: what the implementation may BUILD FROM is public and opened
before the run; what it may only COMPARE AGAINST stays quarantined until commitment.

| Packet | Role | Status |
| --- | --- | --- |
| **group packet**, the audited [`m8_5a_packet.json`](../data/m8_5a_packet.json) | exact `2I` generators in `Q(φ)`, coefficient-field conventions | PUBLIC, already audited, re-audited against § 4.3 |
| **construction packet** (§ 4.2) | a finite based 3-dimensional chain complex representing `S³/2I` | PUBLIC, audited before the room opens |
| **answer packet** (§ 4.4) | the canonical adjudication reference | QUARANTINED until commitment |

**Why a bare presentation is not enough.** `2I` has the balanced presentation
`⟨s, t | s³ = t⁵ = (st)²⟩`, whose presentation 2-complex has `χ = 1 − 2 + 2 = 1`, while a
closed orientable 3-manifold requires `χ = 0`. The gap is exactly the 3-cell. Fox calculus
on a presentation yields `∂₁` and `∂₂`; `∂₃` needs the identity among relations or equivalent
attaching data, and `∂₃` is precisely the boundary map the torsion product consumes.

**Why a bounded reference allowance was declined.** On the precedent that already declined a
web allowlist for M8.5-A: the literature carrying a period-4 resolution for this group
carries the torsion formulas alongside it, frequently on the same page. A reference allowance
broad enough to be useful re-opens the leak the firewall exists to close, and it is not
policeable.

### 4.2 The construction packet: format and key set

Declared here, not merely a contents policy, so the maintainer audit is scopeable at freeze
rather than after it.

**The packet's object is a FINITE BASED 3-DIMENSIONAL CHAIN COMPLEX representing `S³/2I`,
in degrees 0 through 3.** A periodic resolution may be its SOURCE, since `2I` acts freely on
`S³` and so has period-4 cohomology, but a resolution is exact after augmentation while the
finite complex retains its top homology, and the torsion consumes the BASES: the word
"truncated" pins neither where the cut falls nor the basis and orientation the top degree
carries. A resolution-derived packet must therefore declare its exact truncation and capping.

One JSON file, canonical form as § 10.2, exactly these keys:

| Key | Content |
| --- | --- |
| `format_version` | schema string |
| `group_packet_sha256` | the group packet this is matched to |
| `abstract_generators` | the presentation generators by name, each mapped to its canonical element ID in the group packet's closure. **These IDs name elements satisfying the DECLARED RELATORS. Do NOT read them as the group packet's listed generators, and do not treat a mismatch with those as a discrepancy.** A presentation generator need not BE a packet generator: the packet generators have orders 6 and 4, while a standard balanced presentation of `2I` uses an order-10 generator, so at least one ID here will differ from the packet's own list by construction. The IDs are also NOT uniquely determined by the relators: several elements typically satisfy them and generate, so this field records the choice made rather than a forced value, and the enumeration it addresses is fixed by the group packet's canonical ordering, not by this field |
| `model_kind` | `finite_cellular` or `resolution_derived` |
| `degree_range` | `[0, 3]`, frozen |
| `free_ranks` | the based free `Z[2I]`-module rank in each degree |
| `boundary_maps` | `∂ₙ` as matrices over `Z[2I]`, entries as canonical group-ring terms per the encoding below |
| `top_closure` | the top-degree basis and orientation data: the fundamental-class convention the torsion consumes |
| `truncation_rule` | present in EVERY packet: `null` when `model_kind` is `finite_cellular`, otherwise the exact cut and capping that produce the finite complex from the periodic resolution |
| `basing` | basis order, orientation, the augmentation convention, the module side (left or right), the vector convention (row or column), the evaluation convention for group-ring entries (`g ↦ ρ(g)`; any inverse, transpose or dual variant must be declared here or is forbidden), and the direction the boundary matrices act |
| `provenance_id` | an OPAQUE identifier and the source-content SHA-256; where no external construction source is used, the hash covers the frozen derivation instead, per the provenance classes below; no bibliographic pointer |

**Excluded, and the audit checks for each:** evaluated irreducible matrices, characters,
determinants, torsion values, ratios, sector products, any target form, and any decimal
rendering. The packet carries integers and canonical element IDs, nothing evaluated.

**Provenance is split across the firewall.** The packet carries only the opaque
`provenance_id` and the source-content hash. The full citation is held in the
maintainer-side construction-audit artifact until commitment and is published with the
adjudication evidence (§ 8 step 9). The reason is the one that declined the reference allowance: for a
fresh AI context, a title, author or theorem name is itself a pointer capable of activating
recalled target material from exactly the source class that prints the resolution beside
the torsion formulas. The implementer never needs the citation; the model gates verify the
complex, and the audit verifies where it came from.

**The construction-packet audit, scope frozen.** The audit is load-bearing for § 1 branch
four, because the internal model gates cannot identify the intended manifold complex
uniquely. It verifies, from the packet and the maintainer-held provenance record, each
check mutation-tested per the standing rule:

1. the recorded source-content hash is reproduced: for an `extracted` packet, by the
   retained source bytes; for a `derived` packet, by the derivation source and environment
   bytes, together with a check that the derivation consumes only permitted public inputs;
2. the boundary maps, bases, orientation, and top and bottom data are checked: for an
   `extracted` packet, independently against that source; for a `derived` packet, by
   rerunning the deterministic construction under its frozen seed and structural acceptance
   predicate, reproducing the canonical packet byte for byte, and executing the full model
   and mutation gates;
3. the abstract-generator correspondence is verified against the exact group packet;
4. the canonical element encoding and every matrix dimension are verified;
5. every forbidden answer-bearing category of this section is absent;
6. the maintainer-side construction-audit artifact, including its provenance class, its
   provenance record, and for an `extracted` packet its full citation, stays outside the
   room until commitment and is published with the
   adjudication evidence (§ 8 step 9), its bytes verified against the hash frozen at
   § 8 step 1.

**Provenance classes, frozen.** Every packet is `extracted` or `derived`. The class, like
the citation, is recorded in the maintainer-side construction-audit artifact and never in
the packet.

`extracted`: `provenance_id` identifies retained source bytes. The audit establishes
FAITHFUL EXTRACTION from the cited source, not independent discovery of the resolution.

`derived`: no external construction source is used, and `provenance_id` identifies a frozen
derivation. The audit establishes REPRODUCIBLE DERIVATION AND STRUCTURAL VERIFICATION. It
does NOT establish faithful extraction from literature, independent discovery of the
topological model, or uniqueness of the derived complex. Where the derivation selects among
structurally acceptable objects, it states the selection rule and the acceptance predicate.
The predicate must be answer-independent: it may read only structural properties, never a
reproduced quantity or anything computed from one. It must ALSO be exact for the property it
claims, not a finite-sample or finite-prime proxy for it. Answer-independence and
sufficiency are separate requirements, and a rejected packet has already been admitted by a
predicate that had the first and lacked the second: rank at one large prime certifies
rational generation and is structurally blind to a finite index. Where the exact predicate
is one-sided, the direction that can only reject must be reported as inconclusive rather
than as a verdict. The audit reproduces the selection; it makes no claim that alternative
acceptable selections are equivalent to it.

Basis invariance is argued exactly, not numerically, and it covers the TRIVIAL units ONLY.
A change of the top-cell basis by a unit `±g` multiplies the torsion by `det ρ(±g)^{±1}`.
Every element of a finite group has finite order, so `ρ(±g)` has finite order in `GL_d(ℂ)`,
so `det ρ(±g)` is a root of unity and its modulus is exactly 1. Hence `T²_target(ρ)` is
unchanged under the DECLARED top-cell basis changes. This holds in any realization, unitary
or not, and does not depend on the chosen basis. Numerical perturbation evidence is retained
as an implementation check on the code, never as the basis of the claim.

Read the scope narrowly. This argument does NOT extend to a different choice of generator,
because the units of `Z[2I]/(N)` are not exhausted by `±g`. That residue is item 6 of § 9's
explicitly-not-verified list, and it is the sharper statement of the disclaimer above: the
acceptance predicate pins the complex, not the reproduced quantity.

Neither class discharges § 1 branch four on its own. Without the audit, branch four would be
controlled only by internal consistency checks this protocol already concedes cannot
identify the model.

**Canonical nested encodings (frozen).** Free word reduction alone is NOT canonical: `s³`,
`t⁵` and `(st)²` are three distinct freely reduced words for the same central element, so a
word-level encoding lets semantically identical objects differ in bytes. A group-ring entry
is therefore an ordered list of terms `(coefficient, element_id)`, where element IDs are
assigned by the rank of each element's canonical exact coordinate tuple in lexicographic
order over the group packet's 120-element closure; both sides derive the SAME enumeration
from the group packet alone.

**THE ID SORT KEY, stated here because the nearest canonical form in this paragraph is a
DIFFERENT one.** Take the four quaternion components in the group packet's
`quaternion_basis` order. Each is `(A + B·φ)/2` with the denominator FIXED at 2 and NOT
reduced, per the group packet's own `coefficient_form`. The key is the eight signed integers
`A₁, B₁, A_i, B_i, A_j, B_j, A_k, B_k`, the numerator pair per component with the fixed
denominator dropped, compared entrywise as SIGNED integers with the first entry most
significant. Rank is the 0-based position and is the element ID.

**This is NOT the normalized `Q(φ)` triple defined at the end of this paragraph.** That
normalization governs `Q(φ)` VALUES and reduces by `gcd(a, b, c)`; applying it as the ID sort
key instead renames **26 of the 120 IDs**, ranks 70 to 74 and 99 to 119, a block containing
the identity and typically at least one declared generator. The misread is not
self-announcing: element orders and `⟨s,t⟩ = 2I` are blind to it, and `ε(∂₁) = 0` is blind
because `ε` sends every element to 1. What catches it is the relator check, where `order(st)`
reads 6 rather than 4, and the identity failing to sit at rank 119.

**The enumeration is checkable BEFORE anything is built.** SHA-256 over one JSON array of the
120 rank-ordered eight-integer arrays, no whitespace (separators `,` and `:`), integers as
bare decimal, ASCII, no trailing newline:

```text
27ff780d28d5d854d464ead87e8fc20244fac8334bda9f0600c6ee1b3c89561e   2389 bytes
rank 0 [-2,0,0,0,0,0,0,0]   rank 118 [1,0,1,0,1,0,1,0]   rank 119 [2,0,0,0,0,0,0,0]
```

It is a function of the group packet alone, so it is derivable from a § 4.3 permitted input
and quarantines nothing. Check it FIRST: a room that builds under the wrong key otherwise
meets a relator failure against a correct packet with no way to tell whether the packet or
its own enumeration is at fault, and the honest response to that dead end is a defect report
against bytes that are right.

Coefficients on identical IDs are combined, zero terms
removed, terms sorted by ID. Source WORDS, where the construction's source gives them, live
in the maintainer-side construction-audit artifact for provenance, never in the clean-room packet.
Matrices declare their dimensions and are entry-ordered row-major. Every `Q(φ)` value is
the normalized triple `(a, b, c)` for `(a + b·φ)/c` with `c > 0` and `gcd(a, b, c) = 1`.
Under these rules, and only under them, semantically identical objects have identical
bytes, and the audit validates the representation the implementation actually parses.

### 4.3 Permitted and forbidden construction inputs

**Permitted:** the two public packets above; generic algebra over group rings and number
fields; generic representation theory of finite groups; generic algebraic topology of chain
complexes and torsion; standard numerics and exact-arithmetic libraries.

**Forbidden as construction inputs, in any form (call, import, copy, read, execute):**
[`m8_3_mass_reproducer.py`](../scripts/m8_3_mass_reproducer.py); the mode-identity-theory
artifact; the M8.3 method note; the source page; the answer packet before its reveal step;
any spectral fixture, twisted spectrum, zeta value, or heat-kernel datum; and any published
`2I` torsion table, ratio or sector product.

**The M8.5-A packet is re-audited against this list**, not inherited on its M8.5-A audit.
The two tasks have different forbidden sets, and a packet clean for one is not automatically
clean for the other.

### 4.4 The answer packet: schema, controls, and what it does not claim

**It does not establish blindness and does not exclude prior corpus exposure.** It fixes the
adjudication reference and enforces task-time ordering: the fresh implementation commits its
code, its derivation artifacts and its raw output before the canonical reference is opened to
that context. It also gives the comparison harness a canonical machine-readable input,
removing manual transcription and post-hoc target selection. Those are its three jobs and it
claims nothing beyond them.

**Schema, declared at freeze exactly as the construction packet's is.** One JSON file,
canonical form as § 10.2, exactly these keys:

| Key | Content |
| --- | --- |
| `format_version` | schema string, versioned independently of the other packets |
| `target_id` | opaque identifier for this adjudication |
| `adjudicates` | the SHA-256 of the group packet and of the construction packet this reference is valid against |
| `rows` | one entry per irrep: label; label-free identity: the § 5.5 row signature, public and frozen before the run; evidentiary class (`declared_convention`, `free`, `free_orientation_selector`, the `R7` orientation being consumed by the § 5.4 selection rather than declared by the row); and the exact `Q(φ)` value encoded as integer triples `(a, b, c)` meaning `(a + b·φ)/c`, never a decimal |
| `identities` | one entry per slot: the structured row-signature factors and integer exponents defining the formula, plus the expected value in the same exact encoding |
| `indexing_map` | structured and machine-readable; the harness APPLIES it |
| `convention_map` | structured: the § 5.4 bridge and orientation-anchor rule, the § 4.2 basing reference, and the evaluation convention |

**Controls, inherited from the M8.5-B adjudication and binding here.** The harness loads the
canonical bytes directly, hash-verified BEFORE parsing, with no manual transcription anywhere
in the path. `adjudicates` must match the run's actual packets. A mutation of one loaded
reference cell, on an in-memory copy downstream of the completed hash check, must redden the
comparison, so the mutation tests the comparison layer rather than SHA-256. A synthetic
nonidentity indexing fixture proves the map-processing path is operative, since a live map
may be too simple to discriminate applying from ignoring it. And the canonical answer packet
is PUBLISHED once the adjudication record is committed, its bytes verified against the frozen
hash, so the commitment is publicly checkable rather than permanently private.

## 5. The frozen reproduction target

### 5.1 Output surface

All nine irrep closed forms, both Galois ratios, both sector products, in exact `Q(φ)`
arithmetic with `φ = (1 + √5)/2`. No floating-point acceptance anywhere: a form is reproduced
when it is exactly equal in `Q(φ)`, not when it agrees to a tolerance.

### 5.2 Evidentiary content, which is not the same as the output surface

| Class | Members | Standing |
| --- | --- | --- |
| declared ledger convention, BOTH routes | `T²(R0) = 1` | not reproduced evidence. The twisted complex is non-acyclic for the trivial representation, so the same obstruction appears combinatorially |
| global orientation | one bit, selected at adjudication by the § 5.4 rule | not evidence; one bit covers every row at once |
| fully free forms | `R1`, `R2`, `R3`, `R4`, `R5`, `R6`, `R8` | seven independent checks |
| orientation-free value, and orientation selector | `R7` | its RECIPROCAL CLASS `{x, x⁻¹}` is an independent check and is reported as such; which member of the class is correct is consumed by the § 5.4 post-reveal selection |

This decomposition resolves a contradiction in the source obligation, which calls `R0`
asserted and separately calls "the remaining 8 irreps" genuine checks. Both cannot hold. The
output surface is NOT shrunk to seven forms; only the evidentiary weight is apportioned.

### 5.3 The two identity families, at two levels

**Mandatory consistency layer, computed entirely after reveal.** The four identities are
part of the output contract and carry NO independent evidentiary weight, since they follow
by arithmetic from the table.

**They cannot be computed before reveal, and the raw output does not carry them.** Which
rows occupy each formula is not derivable from any permitted pre-reveal input: the Galois
pairs are, being conjugate under the Galois action on `Q(φ)`, but the SECTOR PARTITION is
an M8.3 structural assignment that § 0 places out of scope. Requiring the implementer to
emit four identities whose membership is defined only in the quarantined packet would be
the same defect as requiring a label it cannot know.

**The frozen rule.** The answer packet's `identities` entries carry, per slot, the
structured row-signature factors and integer exponents that define the formula, alongside
the expected value. After the § 5.4 selection, the harness applies the identity or the
global inversion to the COMMITTED ROWS, recomputes all four identities from those selected
rows under the packet's slot definitions, and compares POSITION-WISE.

The substantive control survives intact and is in fact strengthened: the identities are
derived only from committed rows, and with no identity field in the raw output there is
nothing an implementation could hardcode. Recomputation from selected rows is also what lets
`convention difference` be a genuine success rather than an unfalsifiable escape.

**Comparison is POSITION-WISE, never as an unordered set.** Under global inversion the two
sector products exchange values, so a set-wise comparison of that pair passes under inversion
with no orientation handling whatever: a false pass that would mask both a missing selection
step and a swapped sector assignment. Each identity is compared to the reference identity
occupying the same declared slot.

**Optional structural layer.** A ratio or product may additionally be reported as
STRUCTURALLY reproduced when derived without first inserting the individual forms: from the
Galois action on the character field, a pairing of chain complexes, determinant norms, or a
sector decomposition. This is not required for the verdict. A difficult auxiliary derivation
must not be able to block an otherwise valid reproduction of the forms themselves.

### 5.4 The bridge and the anchor, declared before the run

Combinatorial torsion is pinned only once the basing is pinned. The construction packet pins
the basing; the correspondence to M8.3's ANALYTIC normalization then needs a declared BRIDGE
and a declared ANCHOR RULE. The bridge and the rule are fixed here; the rule's one-bit
selection happens only after the answer packet opens, and by the adjudicator, never by the
implementation.

**The bridge (frozen).** The combinatorial routine computes `τ_ρ`, the Reidemeister torsion
of the based acyclic complex `C_* ⊗_{Z[2I]} V_ρ` in the declared bases: a complex number
whose residual indeterminacy under the declared conventions is a unit-modulus factor, so its
squared modulus is well defined. The target identification is

```text
T²_target(ρ) := |τ_ρ|²
```

exactly: the squared modulus, with determinants over `C` in the evaluated representation; no
field norm to a subfield; no absolute value anywhere else; no squaring beyond the displayed
one; and the involution below acting on `T²_target` AFTER this identification. For the
real-type irreps `|τ_ρ|² = τ_ρ²`, a consequence, not a second definition. The `R7` anchor
cannot substitute for this equation: it selects between a quantity and its inverse and can
distinguish neither torsion from torsion squared nor a modulus from a norm.

**Why this `F` and not another.** On a closed 3-manifold, Poincaré duality and the Hodge
split reduce the Ray-Singer combination to `2 log T_an = ζ'_coexact(0) − 2 ζ'_scalar(0)`,
whose right side is the analytic route's own defining formula for `log T²`; Cheeger-Müller
then gives `T_an = |τ|` per acyclic unitary sector, so `T² = |τ|²`. The reduction uses only
generic 3-manifold facts. The one convention that derivation does not fix, the overall sign
of the ζ'-combination, is exactly the involution below, which is why that involution and
nothing else is the admitted bridge.

**What `R7` checks freely, stated exactly.** Not a magnitude: `|log T²|` is not an element
of `Q(φ)` and the raw output carries no logarithm. The exact orientation-independent object
is the RECIPROCAL CLASS. Writing `x` for the committed native `R7` value and `r` for the
revealed reference, the free check is that the unordered class `{x, x⁻¹}` equals `{r, r⁻¹}`,
which is exactly: at least one of `x = r` or `x⁻¹ = r` holds. All comparisons stay in exact
`Q(φ)`.

**The anchor rule, preregistered; the selection, post-reveal.** The implementation declares
its NATIVE orientation convention in the § 8 method-and-gate manifest and commits every
`T²_target` value under that native convention, matching nothing. Only after the answer
packet opens may the adjudicator use `R7` to select between the committed table and its
global inverse: whichever of the two orientations agrees with the analytic convention at
`R7` is the selected one, and the full comparison proceeds under it. **The selection must be
unique, and that is a gate rather than an assumption.** Exactly one of `x = r` and
`x⁻¹ = r` may hold. Both holding means `r = r⁻¹`, a self-inverse reference that cannot
discriminate orientation at all; the adjudicator then records an INVALID ANCHOR and no
orientation is selected, since the protocol would otherwise silently pick one. Neither
holding is § 8 `disagreement`. The gate runs at reveal because the
reference value is quarantined until then. That selection
classifies the outcome and never modifies the implementation or the committed raw output.
The anchor value itself is never published beforehand: pre-publishing it would hand the
implementer the very value `R7` is supposed to reproduce freely.

**The bridge is one involution, and nothing else.** The sole admitted convention bridge
between the routes is the preregistered global sign of `log T²`, equivalently the inversion
`T² ↔ (T²)⁻¹` applied to every row at once, anchored once at `R7`. Outputs not related by
exactly that transformation are a DISAGREEMENT, never a convention difference. The rule is
declared in advance and the selection is mechanical, so nothing about the convention is
chosen after seeing the output beyond that single preregistered bit, and that bit is chosen
by the adjudicator against `R7`, not by the implementation against anything.

### 5.5 The raw-output schema and the public row signature

**The comparison consumes only the committed raw output.** Whatever the adjudication needs,
the raw output already carries; no new mathematical extraction from live implementation
state occurs after the answer packet opens.

**The row signature, public and frozen.** Each irrep is identified label-free by its
dimension together with its exact characters on the classes of the presentation generators
`s`, `t` and of the word `st`, pulled back through the construction packet's declared
generator correspondence, each value a normalized `Q(φ)` triple. The named elements live in
the PUBLIC packets, so the implementing context knows before commitment exactly which
identity fields its output must carry, and the quarantined packet contributes no selector of
its own. **Separation is verified at freeze rather than assumed**: the freeze review
confirms, author-side and outside the room, that this signature separates all nine irreps,
and extends the word list before the room opens if it does not.

| Field | Content |
| --- | --- |
| `schema_version` | frozen string |
| `group_packet_sha256`, `construction_packet_sha256` | the public inputs this run consumed. These are the pre-reveal binding to the adjudication; the opaque `target_id` lives in the quarantined packet and cannot be echoed by output committed before reveal |
| `manifest_sha256` | the § 8 method-and-gate manifest this run executed |
| `rows` | one entry per irrep: the § 5.5 public row signature and, for an ACYCLIC row, the exact `T²_target` triple under the manifest's declared native orientation; a non-acyclic row carries no computed value. Evidentiary classes are ANSWER-SIDE metadata, applied by the adjudicator after reveal and signature matching: the public signature does not say which row bears which label, so an implementer cannot assign them, and requiring it would demand quarantined information before reveal |
| `derivation_artifacts` | the § 7 route-native intermediates, by SHA-256 |
| `gate_results` | every PRE-REVEAL § 9 gate, by identifier, with outcome. The post-reveal gates (anchor uniqueness, selected-orientation equality, stage-two identity comparison, and the § 4.4 answer-packet controls) cannot run before the packet opens and are recorded in the § 8 step 8 comparison output |

## 6. Method disjointness

**The route class is frozen; the construction within it is not.** The implementation must
compute torsion from the supplied finite based chain complex (§ 4.2) by determinant data: a
topological or algebraic Reidemeister-torsion route.

**It may not use** spectral multiplicities, twisted spectra, zeta functions or their
derivatives, heat-kernel data, or any closed-form torsion fixture. Those constitute the
occupied route.

**It may choose** its own realization within the class, given the supplied complex: the
evaluation order, the determinant organization, the acyclicity argument, the exact-arithmetic
representation. Naming one construction would freeze a route that may prove awkward, and the
requirement that carries the evidence is disjointness from the analytic route, not a
particular chain construction.

**Overlap is disclosed, not asserted away.** The implementer records which mathematical
facts, libraries and conventions its route shares with the analytic one. Disclosed overlap is
a finding; undisclosed overlap is a structural failure.

**The theorem-side contract (frozen).** The equality invoked is the classical
Cheeger-Müller statement for ACYCLIC UNITARY representations, applied per irrep only where
its hypotheses hold, never as a forced comparison. The two hypotheses fail differently and
are categorized differently: unexpected NONTRIVIAL non-acyclicity is § 8 `hypothesis
failure`, since the theorem genuinely does not apply there; failure to construct or verify
the unitarity of a finite-group irrep is § 8 `structural failure`, since every such irrep is
unitarizable and a failure is therefore an implementation defect rather than a limit of the
theorem. Every irrep of a finite group is unitarizable, so this is a gate rather than an
obstruction: the implementation must construct an invariant positive-definite Hermitian form
(group averaging suffices) or verify exact unitarity, on the matrices it ACTUALLY CONSUMES.
The flat bundle must be the same one per label as the analytic route's, anchored by the
§ 5.5 PUBLIC row signature, which splits each Galois pair through the group packet's
embedding and is available before commitment;
dimension alone is insufficient, since 2, 3 and 4 each occur twice. The § 4.2 module and
evaluation conventions are part of this contract, and the gate on them is STRUCTURAL rather
than character-based: every `2I` character is real, so a `ρ(g)` versus `ρ(g⁻¹)ᵀ` convention
error is character-invisible; the gate is a declared convention plus a synthetic fixture on
which each wrong convention reddens at least one gate. For unitary representations the
contragredient swap is entrywise conjugation and leaves the squared-modulus torsion, ranks
and `∂∂ = 0` unchanged, so the fixture is chosen non-unitary where needed. **The fixture is
an exact representation and chain-complex instance**, expressed in a deliberately
non-unitary basis where required, and is processed through the SAME parser, group-ring
evaluator, module-action code and boundary-action code as the target run. It passes every
applicable fixture gate under the declared conventions. Each convention is then mutated
SEPARATELY, and each mutation must redden at least one preregistered gate while no other
fixture input changes. This is what makes the mutation causal: it must break a path that was
otherwise mathematically valid, rather than merely show that malformed matrices fail.

**What this does not establish**, stated here and repeated in § 9: the fresh context does not
independently derive the supplied CW model or periodic resolution. It reproduces the torsion
calculation that runs on that model.

## 7. Derivation-path requirement

The response to the recall risk of § 1. Route-independent and frozen now; route-native gates
instantiate after the implementer selects a construction.

**The requirement.** Every reported closed form must be produced from committed route-native
intermediate objects, through the SAME code path that produces the final calculation. The
harness gates those objects and mutation-tests the dependencies connecting them to the
reported result. **A correct final table with absent, inconsistent, or unused derivation
artifacts is a structural failure**, not a pass with a caveat.

**Why the same code path matters.** A parallel reconstruction of what the intermediates
should have been proves nothing about what the answer actually consumed. This is the M8.5-B
lesson: gates must inspect the objects the result was built from.

**Instantiation, once the construction is chosen**, is expected to gate the based twisted
chain groups, the evaluated boundary matrices, `∂∂ = 0` in each representation, ranks and
acyclicity, the determinant factors actually consumed by the torsion product, the basing and
sign normalization, and the Galois action on the consumed matrices or determinants. Mutations
alter those consumed objects and must redden the corresponding gates.

**"After method selection" is a point in the frozen sequence, not a license.** The § 8
manifest step commits the selected route and its instantiated gates BEFORE the target run,
and the sequence makes them irrevocable before the answer packet opens. The minimum
evidentiary machinery does not stay fluid during the run.

## 8. Adjudication

**Sequence, frozen.**

1. The maintainer records the CONTENT COMMIT: this protocol, the group-packet re-audit,
   the construction packet, the frozen SHA-256 of the maintainer-side construction-audit
   artifact, and the answer packet's hash. The audit ARTIFACT and the full citation stay
   outside the room per § 4.2 and enter the pull request at step 9; committing them here
   would place them in the base the room opens from.
2. The maintainer records the LOCK COMMIT: a separate manifest naming the content commit,
   every packet and audit hash, and the clean-room base. **A commit cannot contain its own
   hash**, so the lock manifest, never this file, carries the commit identifiers. From the
   lock commit the addenda-only rule of § 12 binds, whether or not the pull request has
   merged. Any rebase, force-push or other history rewrite after the lock VOIDS the freeze
   and requires a new lock record and a clean-room restart; if the merge strategy could
   rewrite history, both commits are preserved through a durable tag or release object.
3. The clean room opens from exactly the lock commit; the protocol text frozen by the
   process is the version at the content commit.
4. The implementer commits a METHOD-AND-GATE MANIFEST: the selected route within § 6's
   class, the instantiated route-native gates of § 7, the conventions consumed, and the
   declared native orientation of § 5.4.
5. The implementation, its environment record, its derivation artifacts and its raw output,
   under the § 5.5 schema, are committed.
6. The answer packet is opened and its hash verified.
7. The comparison harness LOADS the canonical answer packet directly, with no manual
   transcription anywhere in the path, APPLIES its structured indexing map, LOADS AND
   VALIDATES its convention map without yet applying the orientation, and consumes ONLY the
   § 5.5 committed raw output.
8. The adjudicator performs the § 5.4 orientation selection at `R7`. The harness then
   applies the selected identity or global inversion to the committed rows and recomputes
   all four identities from those selected rows (§ 5.3 stage two). The outcome is recorded
   by category below; the ordering record and the comparison output are committed.
9. The canonical answer packet and the MAINTAINER-SIDE CONSTRUCTION-AUDIT ARTIFACT,
   including its provenance class and record, and for an `extracted` packet its full
   citation (§ 4.2), are PUBLISHED, and each byte stream is verified against its frozen hash.

**Outcome categories**, each a result rather than a failure to report:

| Category | Meaning |
| --- | --- |
| `reproduced` | under the § 5.4 selection resolving to the IDENTITY orientation: the seven fully free forms and the `R7` value exactly equal in `Q(φ)`; both Galois ratios and both sector products consistent; every model, theorem-hypothesis, derivation-path, packet and comparison gate green |
| `partial disagreement` | some forms agree exactly, others do not; the disagreement is diagnosed against the four § 1 branches, and every branch not excluded by its evidence is reported, with multiple surviving branches retained rather than forced to one |
| `convention difference` | the same full agreement, under the § 5.4 selection resolving to the GLOBAL INVERSE: a reproduction under the opposite native orientation, carrying the same claim as `reproduced` with the orientation recorded. The two categories are disjoint and both are successes |
| `hypothesis failure` | a NONTRIVIAL target representation unexpectedly fails the acyclicity gate, so the theorem's hypothesis does not hold there and no comparison is claimed for it. Expected `R0` non-acyclicity is a PASS of the model gates and never this category |
| `structural failure` | ANY failed gate other than an exact row-value comparison, nontrivial acyclicity, or anchor uniqueness, which have their own categories elsewhere in this table; this includes every model, theorem-side, packet and derivation-path gate. Also: derivation artifacts absent, inconsistent, or unused; or undisclosed method overlap |
| `invalid anchor` | BOTH `x = r` and `x⁻¹ = r` hold at `R7`, so the revealed reference is self-inverse and cannot discriminate orientation. No orientation is selected and no comparison verdict is issued; the run is re-adjudicated only under a dated § 12 addendum naming a different anchor row |
| `disagreement` | NEITHER `x = r` nor `x⁻¹ = r` holds at `R7`. No orientation can be selected, so no row comparison is issued, and the mismatch is diagnosed against the four § 1 branches with every surviving branch reported |
| `not completed` | the run did not finish; recorded as such and never as a negative result |

## 9. Gates and the mutation-test requirement

**Every gate carries a runnable mutation that must redden it, with enforced coverage and a
nonzero exit.** Manual attestation is not accepted. A check that cannot fail is not a check.

**Model gates, on the supplied complex.** These are what make the model verified rather than
trusted, and each can go red:

| Gate | Establishes |
| --- | --- |
| `∂ₙ ∂ₙ₊₁ = 0` over `Z[2I]` | the complex is a complex |
| declared free ranks, and `χ = 0` | a closed 3-manifold's complex, not a presentation 2-complex |
| `H_*(Z ⊗_{Z[2I]} C_*) ≅ (Z, 0, 0, Z)`, the trivial module applied through the declared augmentation, tensor side per the declared module convention | the homology of `S³/2I`; a necessary identity check that provenance and the audit complete (§ 1, branch four). This gate and the universal-cover gate below are DIFFERENT computations on different objects, and neither implies the other. An earlier draft of this row asserted that the free complex's plain `Z`-homology returns the same four groups; that is false, and a rejected construction packet demonstrated it by passing this gate while failing the one below |
| `H_*(C_*) ≅ (Z, 0, 0, Z)` for `C_*` read as a complex of free `Z`-modules, that is the universal cover `S³`. Established INTEGRALLY: `im ∂₃ = ker ∂₂` as lattices, from containment (`∂₃∂₂ = 0`), the exact ranks, and an exact SATURATION certificate for `im ∂₃`. Admission rule, stated as a property rather than a list: the certificate must determine the elementary divisors exactly, by integer arithmetic, with no prime and no sample on the accept side. Examples under that rule: a maximal minor of determinant `±1`, Smith or Hermite data, or a unimodular-elimination certificate. Note the asymmetry: a determinant route always returns a number, while an elimination route may stall, and a stalled elimination is INCONCLUSIVE rather than a pass. A mod-`p` battery is a reject screen and is NOT sufficient on the accept side, since no finite prime set excludes an unseen index prime | that `C_*` is the chain complex of `S³` and not merely a complex with the right RATIONAL invariants. A complex can pass every other gate in this table while `im ∂₃` sits at finite index inside `ker ∂₂`: `∂₃ → k·∂₃` leaves `∂∂ = 0`, `ε(∂₃) = 0`, the augmented homology and per-irrep acyclicity all untouched, because each is a rank statement or is computed from `ε(∂₂)` alone, while multiplying every reproduced torsion by `k` to a power. This gate is the only one in the table that sees the integral lattice |
| the terminal map `C₀ → Z` is the declared augmentation `ε` | the declared augmentation, which is `ε`, not `∂₁` |
| `∂₁` matches the frozen 1-cell correspondence, and `ε ∂₁ = 0` | the generator correspondence, and exactness into the augmentation |
| per-irrep acyclicity, with the EXPECTED results frozen: `R0` non-acyclic is PASS; every nontrivial irrep acyclic is PASS; a nontrivial irrep non-acyclic is `hypothesis failure` | the theorem's hypothesis, per § 1. An expected failure is a passing result, so a valid run never triggers the failure category |

**Theorem-side gates**, the § 6 contract made runnable:

| Gate | Establishes |
| --- | --- |
| an invariant positive-definite Hermitian form constructed by group averaging, or exact unitarity verified, on the matrices actually consumed | the unitary hypothesis of the invoked equality |
| row identity resolved per the § 5.5 public row signature, splitting each Galois pair through the group packet's embedding | the same flat bundle per label as the analytic route, using only pre-reveal public data |
| every § 4.2 declared convention (module side, vector convention, evaluation map, boundary direction) exercised on the § 6 synthetic fixture: an EXACT representation and chain-complex instance, non-unitary in basis where required, run through the same parser and evaluation path as the target, GREEN under the declared conventions, then mutated one convention at a time with each mutation reddening at least one preregistered gate and no other input changed. The fixture cannot be only-unitary for the evaluation-map mutation: there the contragredient swap `ρ(g) ↦ ρ(g⁻¹)ᵀ` is entrywise conjugation, leaving `\|τ\|²`, every rank and `∂∂ = 0` unchanged, so that mutation could redden NOTHING | character-invisible convention errors, caught structurally since every `2I` character is real; causally, because the mutation must break an otherwise valid path rather than merely reject malformed input; and the suite's own health, since a mutation no gate can catch is a requirement that cannot be satisfied |

**Reproduction gates.** Exact `Q(φ)` equality for the seven free forms and for `R7` under
the selected orientation; the § 5.4 anchor-uniqueness gate; the § 5.3 two-stage identity
rule, compared POSITION-WISE; the § 7
derivation-path dependencies; the § 4.4 answer-packet controls; the § 6 overlap disclosure
present and non-empty.

**Explicitly NOT verified, and reported as such:**

1. The supplied CW model or periodic resolution is not independently derived. It is supplied
   and VERIFIED by the model gates above, which is more than "trusted" and less than
   "discovered".
2. `T²(R0) = 1` is a declared convention in both routes and is not reproduced.
3. The global orientation is selected at adjudication by the § 5.4 rule, not reproduced.
4. The mass comparison, the fermion assignment, and the dead-zone entries' physical status.
5. Absence of prior corpus exposure, which no procedure here can establish.
6. **That the reproduced quantity is independent of WHICH generator was installed as `∂₃`.**
   `C₃` has rank 1, so any two generators of `ker ∂₂` differ by `∂₃' = u·∂₃` for a unit `u`
   of `Z[2I]/(N)`, and on every nontrivial irrep, which kills `N`, the torsion moves by
   `det ρ(u)`. The basis-invariance argument covers only the TRIVIAL units `u = ±g`, where
   finite order forces `|det ρ(±g)| = 1` exactly. It does not cover the rest, and the rest
   are not empty: `2I` has 9 complex irreps in 7 Galois orbits with every character real, so
   `rank Wh(2I) = 9 − 7 = 2` and nontrivial units exist, with no reason for `|det ρ(u)|` to
   be 1. Under a `derived` packet the acceptance predicate therefore pins the COMPLEX and
   not `T²_target`. Bounded searches have turned up no second accepted generator, which is
   evidence about small support and not a theorem.

## 10. Provenance and the commitment

### 10.1 The commitment record

Committed BEFORE the answer packet opens: the § 8 method-and-gate manifest and its
SHA-256; implementation source and its SHA-256; the environment record; the consulted-files
manifest, including anything that loaded without being asked for; the raw output and its
hash; the derivation artifacts of § 7; and the schema version. The transcript's SHA-256 is recorded at teardown per the clean-room
standard § 10, with the file's location and its backup named, and never only in a
gitignored checkpoint.

### 10.2 Canonical form and hashing

All packets: JSON, keys sorted, two-space indent, ASCII, LF, single trailing newline. The
incoming hash of delivered bytes is recorded, the canonical serialization declared, and the
authoritative hash of the canonical form issued, stating whether canonicalization changed
the bytes. A hash over an uncanonical rendering pins a transcription, not an object.

### 10.3 Packaging

The protocol, the construction packet, its audit, and whatever else the run requires land as
a SINGLE pull request. One document, one merge. **The freeze is the § 8 lock commit, not
the merge event**: implementation and adjudication commits may accumulate in the same pull
request after it, and any later protocol change is a dated § 12 addendum even before the
merge.

## 11. Pins

| Pin | Value |
| --- | --- |
| M8.3 method note and verification record (the obligation discharged; the § 1 analytic-side control) | [PIN at landing] |
| group packet SHA-256 | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` |
| construction packet SHA-256 | [PIN at landing] |
| answer packet SHA-256 | [PIN at landing] |
| clean-room standard | [PIN at landing] |
| repository commit this protocol was drafted against | [PIN at landing] |
| the § 8 content and lock commits | carried by the LOCK MANIFEST, never by this file: a commit cannot contain its own hash |

## 12. Addenda (post-freeze only)

From the § 8 lock commit this document is FROZEN; the merge records the freeze, it does
not create it. Changes enter only as dated addenda in this section,
never in place. An in-place edit is a breach of the freeze, and the pinned paths of § 11 may
not move.

(none)
