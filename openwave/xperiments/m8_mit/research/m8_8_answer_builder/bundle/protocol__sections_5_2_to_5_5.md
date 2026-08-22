# Protocol § 5.2, § 5.3, § 5.4, § 5.5, verbatim

Sliced from `protocol__m8_8_reproduction_protocol.md`, sha256 `8bf24b9fe23e2c2f182c9fabb98a84680fc12a6ab7dec47ccc9a8b3d8bbe306d`.
Reproduce with: `python3 p1_snapshot.py <out> protocol__m8_8_reproduction_protocol.md ...`

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
