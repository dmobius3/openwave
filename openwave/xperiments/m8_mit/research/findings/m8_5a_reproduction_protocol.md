# M8.5-A Independent-Reproduction Protocol (context-isolated; scalar first-occurrence table)

> **Status: LOCKED 2026-07-30** ([PR #380](https://github.com/openwave-labs/openwave/pull/380)).
> Author-written protocol, filed before any implementation existed; NO implementation here.
> Freezes HOW the M8.2 lock's § 3 reproduction obligation is discharged; does not
> discharge it. FROZEN as of the lock date: later changes enter only as dated addenda
> in § 11, never in place. Pairs with the maintainers' task spec
> [`m8_5_task_details.md`](../tasks/m8_5_task_details.md) and the locked
> [`m8_2_preregistration.md`](m8_2_preregistration.md).
> Owner: Blake Shatto. Reviewer and step-2 implementer: maintainers.

## 0. Scope

This protocol governs **M8.5-A only**: the context-isolated, independent-method reproduction of
the frozen scalar (0-form) first-occurrence table, discharging the M8.2 lock's § 3 clause:

> "M8.5 must reproduce every certification table through an INDEPENDENTLY implemented
> decomposition. It may compare against `m8_2_first_occurrence.py` but may not call it, import
> its tables, or share its derived fixtures."

The scalar table is the only computed, frozen certification table (the M5 spin-2 table is
unbuilt, pending V1b; the M7 twisted coexact table was removed as a native target), so it is
the sole reproduction target here.

**Out of scope, stated explicitly:**

1. **M8.5-B** (the quotient simulation backend: prototypes, spectral benchmarks, resolution
   studies, route decision, backend certification) is specified separately when those
   implementations exist. Nothing in this protocol certifies a backend, and per the task spec,
   A gates any claim that B is certified.
2. The **coexact one-form entry rule** is not a certification target. It appears only as a
   separately labeled ASSERTED adjudication module (§ 6), whose numerical agreement cannot
   raise its evidential standing.
3. **Independent-method reproduction of the corrected M8.3 Reidemeister-torsion closed forms
   is not part of M8.5-A or M8.5-B.** It remains a separately tracked future verification
   obligation (queued in [`m8_3_method_note.md`](m8_3_method_note.md)).

## 1. The principle, and the real risk

**Protocol authorship does not compromise independence; target-aware implementation does.**
The answer-holding party writes the rules; the maintainer reviews and freezes them before a
fresh party runs the check.

The dominant risk is not importing `m8_2_first_occurrence.py`. It is **semantic and
mathematical circularity** through routes a file firewall does not close:

1. **Visible target structure.** Parts of the M8.2 record are public (paper and repository),
   so task-time isolation, not universal ignorance, is what the firewall can actually provide.
   The claim ceiling (§ 2) prices this in.
2. **Shared representation machinery.** If one module both builds the decomposition and
   supplies the expected answers, "the numbers match" is a tautology: the `m7_trivial_ok`
   failure one level up (a check whose two sides are the same expression). The gates (§ 8) and
   the derive-everything input rule (§ 4) close this.

## 2. Claim ceiling and label

A successful exact reproduction under M8.5-A is reported as **context-isolated
independent-method reproduction**, and nothing stronger. Other outcomes retain their § 7
category (partial disagreement, structural failure, or not completed): a disagreement or an
incomplete run is not a "reproduction," and the claim label never swallows the
negative-result taxonomy. No outcome, in either direction, earns a stronger label.

No stronger label is available for this table: the author's context generated it, and the
maintainer already reproduced every row during the PR #350 review. The ceiling is set by what
the verifying contexts have seen, not by how separate the implementation is (roadmap
§ CONVENTIONS). For an AI implementer the cap is structural: the training corpus is opaque, so
prior exposure to the published values cannot be excluded even under perfect task-time
isolation. Isolation still earns its keep (it maximizes provenance and semantic independence);
it never upgrades the label.

In the successful-match case, agreement of the fresh implementation with both existing
artifacts (§ 3) is reported as **three-way agreement**: a provenance statement, not a stronger
verification label, and never "three mutually independent methods" (§ 3's overlap disclosure
exists precisely because the methods may overlap).

## 3. The three objects, the quarantine, and the ordering record

| # | Object | Identity | Status for the implementer |
| --- | --- | --- | --- |
| 1 | the M8.2 generator | [`m8_2_first_occurrence.py`](../scripts/m8_2_first_occurrence.py), last modified at commit `269456b7` | FORBIDDEN until the § 9 commitment is complete; then an adjudication reference only |
| 2 | the maintainer reconstruction | [`m8_2_indep_reconstruction.py`](../scripts/m8_2_indep_reconstruction.py) + [`m8_2_indep_reconstruction_note.md`](m8_2_indep_reconstruction_note.md) + its raw/JSON outputs, landed at commit `12f4a94a`, self-verified against `d6790eef` (python 3.13.9, numpy 2.3.5, seed 11) | FORBIDDEN until the § 9 commitment is complete; then an adjudication reference only |
| 3 | the fresh M8.5-A implementation | the deliverable | built under §§ 4-5, committed under § 9 |

The quarantine extends beyond the two scripts: the locked pre-registration (its § 6.1 holds
the target values), the pinned MIT source papers, and every M8.2-derived table or JSON are
equally FORBIDDEN to the implementer until the commitment is complete.

**Ordering record (the enforcement).** The implementer commits source, environment record, raw
output, consulted-files manifest, and method-note draft FIRST (§ 9). Commit timestamps and the
SHA-256 hashes of the raw output are the durable record that the quarantine held. Only then is
anything quarantined unsealed, and from that point objects 1 and 2 serve as adjudication
references only. Any post-commitment change to code or output happens only as a separately
recorded, dated rerun; the original commitment stands as filed.

**Method-overlap disclosure (required; novelty is not).** The two obvious derivation routes
are already occupied: object 1 proceeds by the McKay recursion with its group-theory inputs as
literals; object 2 proceeds by explicit quaternions, brute-force conjugacy classes, and
Burnside class-sum diagonalization. The fresh implementation is NOT required to invent a third
fundamentally different method. It IS required to be context-isolated, implementation-independent
(no calling, importing, copying, or sharing derived fixtures with either object), and to carry
in its method note an explicit statement of where its mathematics overlaps either existing
route. Overlap is disclosed, never penalized.

## 4. Context firewall and construction inputs

**The firewall (mandatory step-2 environment):**

- a fresh context: no M8.2 internals loaded, no prior M8 session state;
- during implementation, no access to the openwave M8 tree beyond the audited packet, and no
  access to the mode-identity-theory repository;
- no web access, or a maintainer-approved allowlist that excludes the MIT papers and
  repositories and any `2I`-specific derived data;
- no persistent memory carrying M8 facts;
- a manifest of every file and reference consulted, committed with § 9.

**Construction inputs (permitted):** the group input in one of two forms: raw explicit
generators given as `SU(2)` matrices or unit quaternions, or an abstract presentation together
with a frozen explicit embedding into `SU(2)`. An abstract presentation without embedding data
is insufficient: it does not canonically distinguish the defining representation `Q` from its
Galois partner `Q′` (the two are exchanged by an outer automorphism), and the `standard` and
`galois` columns are separate outputs, so the embedding is what anchors them. Also permitted:
the group order, 120; generic finite-group algorithms (closure, conjugacy classes, character
extraction, tensor decomposition); generic `SU(2)` representation-theory and harmonic-analysis
references; standard numerics libraries.

**Forbidden as construction inputs:** objects 1 and 2 in any form (call, import, copy, read,
execute); any M8.2-derived JSON or table; the locked pre-registration; the MIT source papers;
and published `2I`-SPECIFIC derived data: exact character tables, McKay graphs or adjacency
data, tensor-product or multiplicity tables, distance tables, first-occurrence tables. Those
are answer keys; the implementer derives them.

**Packet audit (mandatory).** The packet that enters the clean room is this protocol plus the
group input (raw `SU(2)`-matrix or unit-quaternion generators, or a presentation with its
frozen embedding). The maintainer audits it before the clean room opens: no derived
annotations (no irrep labels, dimensions, distances, or character values), no worked examples
that encode targets, no target values. The group input must be raw elements or a bare
presentation-plus-embedding only; that resolves its dual status (a permitted construction
input, and the one input the author could accidentally decorate with answer-bearing
structure). The author knows the answers; the audit is the guard against author-side leakage.

## 5. The frozen mathematical object (scalar reproduction target)

Frozen so that two correct implementations are comparable:

1. **Arena.** `S³ ≃ SU(2)`; deck group `2I ⊂ SU(2)` acting by left multiplication; quotient
   `S³/2I`.
2. **Scalar tower.** `V_n` is the `(n+1)`-dimensional irreducible `SU(2)` representation,
   `n ≥ 0`, with eigenvalue `λ(n) = n(n+2)/R²` and `j = n/2`. (`L²(S³)` carries `V_n` with an
   extra multiplicity factor from the right action; a nonzero factor cannot change a first
   occurrence, so it is not reported.)
3. **Connections and coefficients.** `σ` ranges over the three declared flat-connection
   classes: trivial; the defining 2-dimensional representation `Q`; its Galois partner `Q′`.
   `Q` is IDENTIFIED, never declared: it is the 2-dimensional irreducible whose character
   equals `2cos θ(g)` on every conjugacy class, where `θ(g)` is the `SU(2)` rotation angle of
   `g` under the packet's frozen embedding (§ 4); `Q′` is the unique other 2-dimensional
   irreducible. Contract input, declared from the lock's § 2 and not re-derived: the
   coefficient representation is `τ_σ = Sym²(σ)`. The three columns are named `trivial`,
   `standard` (for `Q`), `galois` (for `Q′`).
4. **The reproduced quantity, exactly:**

   ```text
   n_first(ρ, σ) = min{ n ≥ 0 : ⟨ χ_{V_n}|_2I · χ_{τ_σ}, χ_ρ ⟩_2I > 0 }
   ```

   with `⟨a, b⟩_2I = (1/120) Σ_g a(g) · conj(b(g))`, multiplicities counted over ℂ. All `2I`
   characters are real, so conjugation and dual conventions cannot flip a result; the
   implementation states its convention anyway. The quantity inside the bracket is a
   nonnegative integer in exact arithmetic. **Integer-nearness rule (frozen):** every
   character inner product used as a representation multiplicity anywhere in the
   implementation (tensor, restriction, adjacency, and first-occurrence multiplicities alike)
   must lie within its stated tolerance of a nonnegative integer; any violation is a
   structural failure of the run. The floating-point occurrence test `> 1/2` applies only
   after that gate passes, so a numerically untrustworthy value can never silently become
   "present."
5. **Termination.** The search bound is `n ≤ 24`. Any `(ρ, σ)` with no occurrence within the
   bound is reported NOT-FOUND and the run FAILS LOUD; silent omission is not an option.
6. **Exactness.** First-occurrence outputs are exact integers; adjudication comparisons are
   exact equality, no tolerances. Floating point may appear only in the numerical
   representation-theory pipeline; every quantity interpreted as an integer is governed by the
   § 5.4 integer-nearness rule, and every float-dependent gate has a preregistered tolerance
   under § 8.
7. **Row identity is label-free.** Each irreducible is reported by its
   `(dimension, McKay-distance)` signature: the McKay tensor-multiplicity matrix is derived
   in-implementation as `A_ρτ = ⟨χ_ρ · χ_Q, χ_τ⟩_2I`, its entries subject to the § 5.4
   integer-nearness rule and rounded only after that gate passes; the graph is the support
   `A_ρτ > 0`, with the integer multiplicities derived, never assumed binary by thresholding;
   and the distance is graph distance from the trivial node by BFS. Irrep names are never
   imported and never matched on.

**Output schema (fixed):** `schema_version: "m8_5a-v1"`; one row per irreducible:
`{dim, mckay_distance, n_first: {trivial, standard, galois}}`; plus the environment record,
the § 8 gate results, the consulted-files manifest, and, if run, the § 6 module block.

## 6. The coexact module (ASSERTED rule; adjudication only)

**The asserted rule** (public in the record; task-time visibility is conceded and priced into
the § 2 ceiling): the coexact 1-form entry rule assigns level `d` for `d ≥ 2`, level `2` for
`d = 0`, and level `3` for `d = 1`, where `d` is the McKay distance; the source convention
prints coexact level `m` with eigenvalue `m²/R²`.

**Status: ASSERTED.** No independently derived coexact target exists anywhere in the record.
Object 2 reconstructed the scalar table only, and the pre-registration § 6.1 table agreeing
with mass-spectrum § 4 is the generator's own internal cross-check, not a second independent
result. Nothing in this module can convert a match into independent verification.

**The module, if run:** the implementer derives the coexact one-form tower on `S³` and its
`2I`-decomposition by its own harmonic analysis (generic references permitted), computes
`m_first(ρ, σ)` under the same definition pattern as § 5.4 applied to that tower, with the
same `≤ 24` bound and fail-loud rule, and supplies an explicit map from its level convention
to the source convention before adjudication.

**Pre-declaration:** whether the module runs is declared in the § 9 commitment, before
anything is unsealed. Adding it afterward is not permitted.

**Verdict categories (exactly four, pre-declared):**

| Verdict | Meaning |
| --- | --- |
| structurally derived and reproduced | the implementation supplies its own general argument yielding the rule, and its computed first occurrences agree |
| numerically consistent, not derived | computed first occurrences agree across the declared range; no general derivation supplied |
| contradicted | at least one cell disagrees under exact arithmetic, a verified convention map, and passing § 8 gates |
| not resolved | the implementation cannot adjudicate; the reason is recorded (a convention mismatch surviving the explicit mapping step, NOT-FOUND cells, or a gate failure confined to this module) |

**What counts as a general argument (frozen):** for "structurally derived and reproduced," the
argument must establish the rule for the full representation family from the operator and
representation structure; a pattern inferred or interpolated from the finite computed
first-occurrence range does not qualify, however clean the fit. The argument is separately
examined in the adversarial audit (§ 9).

**Standing rule (hard):** numerical agreement in any amount never upgrades the rule's standing
and is never reported as independent verification. Only "structurally derived and reproduced"
changes the rule's standing, and it is the derivation that does it, not the match.

## 7. Adjudication

**Target, pinned.** The reproduction target is the scalar first-occurrence table in the locked
pre-registration's § 6.1, exactly as at commit `ec877ee0` (upstream main at filing; § 6.1 is
unchanged since the close-out commit `269456b7`, the file's only later change being the
appended 2026-07-29 addendum, which touches no table). The adjudication source is § 6.1 and
only § 6.1: mass-spectrum § 4 is not a second target, since its agreement with § 6.1 is the
generator's own cross-check.

**Transcription.** The adjudicator transcribes § 6.1 at the pinned commit into the comparison
harness. The harness carries a transcription mutation (perturb one transcribed cell and the
comparison must go red): the `doc_typo` pattern from object 2, so the comparison itself is a
check that can fail.

**Matching.** Rows are matched label-free by the `(dim, mckay_distance)` signature. Before any
cell comparison, the harness verifies that the signatures on each side are pairwise distinct
(matching is otherwise ill-posed) and that the row counts agree.

**Comparison.** Exact integer equality, cell by cell, across all three columns.

**Scalar result categories** (any outcome is reportable; a disagreement is a finding, not a
failure of the task):

| Category | Meaning |
| --- | --- |
| reproduced | every cell matches |
| partial disagreement | each mismatching cell listed with both values |
| structural failure | a § 8 gate fails, signatures are not distinct, or NOT-FOUND cells exist |
| not completed | the run did not reach adjudication; the reason is recorded |

**Reporting.** Results are reported under this protocol as first frozen: no post-hoc
reclassification, no scope additions, no label movement (§ 2). Agreement of object 3 with
objects 1 and 2 on the scalar table is reported as three-way agreement.

## 8. Gates and the mutation-test requirement

Every line the implementation prints as PASS must be demonstrated to fail under a deliberate
defect BEFORE the deliverable ships, via a RUNNABLE mutation harness (a `--mutation-tests`
mode or equivalent), not manual attestation. Coverage is enforced: every gate goes red under
at least one mutation, and the harness exits nonzero if any mutation fails to redden its
target or any gate is uncovered. Precedent: the `m7_trivial_ok` defect (a check whose two
sides evaluate the same expression) and its replacement by falsifiable gates.

**Minimum gate set** (generic and theorem-anchored; no target values appear in any gate):

| # | Gate |
| --- | --- |
| G1 | group order 120 and closure of the constructed element set |
| G2 | conjugacy-class partition consistency (sizes divide the order and sum to it) |
| G3 | character-table orthonormality, with its stated float tolerance |
| G4 | sum of squared irreducible dimensions equals the group order |
| G5 | `Q` uniquely identified by `χ_Q(g) = 2cos θ(g)` on every class under the packet's frozen embedding; `Q′` the unique other 2-dimensional irreducible |
| G6 | McKay adjacency mark condition `A·dims = 2·dims` holds as a consequence |
| G7 | distances defined for every irreducible; `(dim, distance)` signatures pairwise distinct |
| G8 | the § 5.4 integer-nearness rule holds for every multiplicity computed (tensor, restriction, adjacency, first-occurrence); any violation is a structural failure |
| G9 | the McKay tensor-multiplicity matrix `A_ρτ` is symmetric |
| G10 | the § 7 comparison harness reddens under the transcription mutation |

The implementer may add gates; every added PASS line joins the mutation requirement. ASSERTED
inputs stay labeled ASSERTED and receive no self-check that cannot discriminate.

**Tolerance commitment (frozen).** Every floating-point tolerance used by G3, G5, G8, or any
added gate must be fixed in the source or configuration and justified in the pre-commitment
method-note draft (§ 9), before adjudication. Tolerances may not be widened after quarantined
targets are opened. A tolerance change requires a separately recorded, dated rerun that
preserves the original source, output, and verdict. The deliverable reports stability of every
gate under increased working precision or an equivalent numerical-conditioning check.

## 9. Provenance and the commitment

Committed, dated, and immutable, BEFORE anything quarantined is unsealed:

- the source;
- the environment record (interpreter and library versions, OS, hardware, seeds if any);
- the raw output and its SHA-256;
- the schema version;
- the consulted-files manifest;
- the declaration of whether the § 6 module ran;
- the method-note draft.

The final deliverable then follows [`dev_docs/METHOD_NOTE.md`](../../../../../dev_docs/METHOD_NOTE.md)
in full: equations first; an equation-to-code map with permalinks; each result beside its gate;
at most ~4 inspection artifacts; the recorded adversarial audit before it crosses to the
author; and an explicit list of claims NOT verified, which must include at minimum: the truth
of the coexact entry rule (unless the § 6 verdict is "structurally derived and reproduced");
the contract choice `τ_σ = Sym²(σ)` (its consequences are checked, the choice is not); the
fitness of any simulation backend (that is M8.5-B); and the M8.3 torsion closed forms (out of
scope, § 0).

## 10. Pins

| What | Pin |
| --- | --- |
| protocol filed against | upstream main `ec877ee0` (2026-07-30) |
| M8.2 lock | landed `f18daf27`; close-out `269456b7`; post-lock addendum `0776cc19` (append-only) |
| target table | pre-registration § 6.1 at `ec877ee0` (content unchanged since `269456b7`) |
| object 1 | `scripts/m8_2_first_occurrence.py`, last modified at `269456b7` |
| object 2 | `scripts/m8_2_indep_reconstruction.py` + note + raw/JSON, landed `12f4a94a`, self-verified against `d6790eef` (python 3.13.9, numpy 2.3.5, seed 11) |
| task spec | [`m8_5_task_details.md`](../tasks/m8_5_task_details.md), § "Independent reproduction, M8.5-A" (added 2026-07-28) |
| procedure record | the PR #350 close-out thread (2026-07-28) |

## 11. Addenda (post-freeze only)

(none yet)
