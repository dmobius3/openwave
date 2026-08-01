# M8.5-B Numerical Pre-registration (quotient backend: scalar and one-form infrastructure)

> **Status: DRAFT FOR FREEZE (filed 2026-07-30; author-written pre-registration; NO
> implementation code, and NO target or certification runs, in this document; the
> non-evidentiary engineering pilot is reported in § 6).** Freezes WHAT the M8.5-B backend must demonstrate, HOW its
> evidence is assembled, and what may be claimed at each stage; it does not build the backend.
> Pilot-selected values are recorded in § 6, and published reference pins with their
> convention maps are recorded in § 11. All slots are filled; after landing, this document is FROZEN and later changes enter only as
> dated addenda in § 12, never in place. Pairs with the maintainers' task spec
> [`m8_5_task_details.md`](../tasks/m8_5_task_details.md) and the frozen
> [`m8_5a_reproduction_protocol.md`](m8_5a_reproduction_protocol.md).
> Owner: Blake Shatto (author-driven, platform support per the task spec). Reviewer: maintainers.

## 0. Scope

M8.5-B builds and certifies the simulation infrastructure for the arena: fields on the
compact quotient `S³/2I`. Two candidate routes are prototyped far enough to choose one on
measured evidence (§ 2). No family-specific target dynamics are used as certification
benchmarks, and nothing certifies on the target problem.

**The four evidence classes (frozen).** Certification is assembled from evidence of four
distinct strengths, named and reported separately, never as equivalent (§ 9):

| # | Class | Where it applies |
| --- | --- | --- |
| 1 | published-value lookup | rungs 1-2 (untwisted `S³`); and rung 3a as a PRECOMMITTED POST-FREEZE published-value adjudication |
| 2 | source-derived theorem evaluation | rung 3b, as a PRECOMMITTED POST-FREEZE theorem evaluation: a pinned published RULE evaluated, because the verified source audit located no published one-form multiplicity lookup for the selected adjudication case |
| 3 | cross-route agreement | rungs 4-6 (`S³/2I`), where the two routes derive independently and are compared |
| 4 | post-run author-context reference agreement | the § 2 post-run references, consulted only after the frozen outputs are committed |

Class 2 is weaker than class 1 and is never described as a known-answer lookup: what is
published there is a theorem, not a number.

**Claim ceiling on rungs 3a and 3b (frozen).** Because the author-side source audit had
access to the adjudication case and its reference values before implementation, **these gates
carry no unseen-target claim**: the author has seen the answers, and no report may describe
them as held-out-from-the-author, previously unseen, or verified by a party ignorant of the
values. Their evidentiary force comes from four things and is stated in exactly those terms:
exclusion from all tuning; precommitment of the numerical contract; preservation of the
sealed reference packet; and comparison only after both routes' outputs are committed.

**The certification cap (frozen):** M8.5-B certifies the frozen scalar and one-form quotient
infrastructure only. It does not certify the M4 transport rule, the M5 tensor operator, the
M7 twisted operator, or any M8.4 dynamics. Those require their M8.2 § 7 execution appendices
and family-specific tests.

**Out of scope, stated explicitly:**

1. **The conic Möbius band**: its eigensolve belongs to M8.1 (closed) and its placement
   relative to `S³/2I` is deferred by the M8.2 lock behind a frozen edge operator plus
   bulk-to-edge correspondence, scored only in M8.4. Nothing here touches, places, or
   simulates the band, and no route design may pretend to locate it on a grid.
2. **Family-specific operators**: the M4 tangent-bundle transport choice, the M5
   degenerate-kinetic fluctuation-operator selection, and M7's family twisted operator are
   per-family execution-appendix content. B's architecture stays operator-agnostic; it must
   be able to host whichever choice each signed appendix freezes, and it selects none of them.
3. **All M8.4 dynamics** (Lagrangian families, vacua, defects, evolution).
4. **The M8.5-A reproduction** (its own frozen protocol governs it) and **the M8.3
   Reidemeister-torsion closed forms** (excluded from M8.5-A and M8.5-B alike; separately
   tracked).

## 1. Claim ceiling and frozen claim language

The task spec's rule is that A gates any claim that B is certified, while B may be built and
run first.

**The M8.5-A gate condition (frozen).** B may claim certification only if M8.5-A's CORE
SCALAR outcome is `reproduced`. An A outcome of `partial disagreement`, `structural failure`,
or `not completed` blocks B certification pending explicit adjudication: A closing is not the
same event as A succeeding. A's optional coexact-module outcome is separate: it affects the
combined coexact record under § 7 and never substitutes for the core scalar gate.

The following phrasings are frozen now so no write-up improvises:

- Only the pre-`2I` rungs complete (rungs 1, 2, 3a, 3b; the `2I` suite has not run):
  > **M8.5-B passed the preregistered pre-`2I` benchmark ladder (untwisted `S³` plus the
  > precommitted scalar adjudication and one-form theorem evaluation); the `S³/2I`
  > benchmark suite has not yet run, and no backend-certification claim is made.**
- A prototype that passes every mandatory rung before A's core outcome is known:
  > **M8.5-B prototype passed the preregistered numerical benchmark suite; backend
  > certification remains pending completion of M8.5-A and cross-artifact adjudication.**
- After A's core scalar outcome is `reproduced`, but before B's full definition of done:
  > **Backend candidate benchmarked under the frozen scalar/one-form scope; not yet
  > certified.**
- Only after A's core scalar outcome is `reproduced` and B's definition of done is met:
  > **M8.5-B certified for the frozen scalar and one-form quotient-infrastructure scope.**

The unqualified phrase "the `S³/2I` backend is certified" is never used: it would read as
covering the M4/M5/M7 operators and M8.4 physics, which § 0 excludes. The claim-label
discipline of the M8.5 task also binds every B write-up: results are described by the frozen
phrases above and by their § 4 rung and § 8 gate outcomes, never by any stronger verification
label.

## 2. The two routes, and what they may share

Both routes are defined for a FROZEN FINITE DECK GROUP, not for `2I` alone: the rung-3a/3b
adjudication case is a different group, and a route hard-wired to the 120-element group could
not run it. `Γ = 2I` is the specialization used for rungs 4-6.

**Route domain (frozen).** The routes are defined for a FINITE EFFECTIVE ACTION IN
`SU(2)_L × SU(2)_R`, not for a subgroup of `SU(2)`. A deck element is the isometry
`x ↦ u x v` given by the ACTION PAIR `[u, v]` of unit quaternions, and the routes consume
the raw action pairs. `Γ = 2I` is the specialization `[γ, 1]`, and any left action is
`v = 1`. The pair representation is two-to-one onto the effective action, since `[u, v]` and
`[−u, −v]` induce the same isometry, so the EFFECTIVE order is what the gates check.
**This is not a generalization for its own sake.** An INHOMOGENEOUS lens space is not
a left action: writing the action on `S³ ⊂ C²` as `(z₁, z₂) ↦ (ζ^{s₁} z₁, ζ^{s₂} z₂)`, left
multiplication alone forces `s₁ = s₂` and right multiplication alone forces `s₂ = −s₁`, so any
case with `s₂ ≢ ±s₁ (mod q)` requires both factors.

**Two consequences for the route architecture, both frozen in § 6.1b.** Route (a)'s
one-form reduction is `Ad`-equivariant rather than componentwise, because a right factor
carries the left-invariant coframe. Route (b) averages characters over both factors, since a
two-sided element acts on `V_n ⊗ V_n` through each of them.

**The coefficient law, derived and verified pointwise (frozen).** For `ω = f_j σ^j` in the
left-invariant coframe and `F(x) = u x v`, the pullback carries components by
`(F*ω)(x) = Ad(v) f(F(x))`, so `F*ω = ω` is equivalent to

```text
f(gamma . x)  =  Ad(v_gamma^-1) f(x)
```

**The inverse is load-bearing and no character check can find it.** `Ad(v)` and
`Ad(v⁻¹) = Ad(v)ᵀ` are orthogonal with the same `SU(2)` character, since `χ_n(v⁻¹) = χ_n(v)`, so
route (b) agreement and every multiplicity comparison are insensitive to substituting one
for another. The law above is therefore fixed by a manufactured pointwise pullback test that
separates the candidates at 1e-14 against 1e+2, and that test is a required gate (§ 8) rather
than a development convenience.

The scalar rung is the lighter of the two: it needs the isometry applied as a `4 × 4`
rotation and the two-sided character average.

| Route | Sketch |
| --- | --- |
| (a) `Γ`-equivariant numerical route | an `S³` grid (embedding or intrinsic charts) with the `\|Γ\|`-element identification imposed as an equivariance constraint or ghost-cell map; multiplicities MEASURED character-free: the numerical group action applied to computed eigenspaces, with invariant and covariant dimensions obtained by direct constraint solving, commutant decomposition, or transformation-matrix ranks |
| (b) spectral-harmonic route | fields expanded in `S³` harmonics restricted to `Γ`-invariant or covariant subspaces; coefficients represented and transformed (evolution belongs to M8.4 and is outside B); multiplicities PREDICTED by character averaging from the route's own derived character data |

**Method disjointness (frozen).** Route (b) owns the character-averaging prediction. Route
(a)'s primary multiplicity measurement is CHARACTER-FREE; character-based central projectors
may appear in route (a) only as a separately labeled post-run cross-check, never as the
primary measurement. This is what makes the two derivations disjoint in method, not merely
in code: without it, both routes would be character theory twice over, satisfying data
isolation while gutting the two-derivations claim.

**Route-specific construction rule (frozen; the `Q`/`Q′`/McKay clauses apply to the `Γ = 2I`
rungs 4-6, where those objects are defined).** Each route independently constructs its
irreducible sectors, `Q′`, and McKay graph; only the raw embedded group elements are shared.
Route (a) may not import character labels, character-derived `Q′` matrices, or a
character-derived adjacency graph. It identifies irreducible blocks from the numerical group
action by commutant or intertwiner decomposition; identifies `Q` by direct comparison with
the frozen `SU(2)` embedding and `Q′` as the unique other 2-dimensional block; builds
`Sym²(Q′)` from its own `Q′` block matrices; and computes
`A_ρτ = dim Hom_2I(τ, ρ ⊗ Q)` through character-free intertwiner-rank calculations. Route
(b) derives the corresponding objects through character averaging. The `Q` matrices are raw
inputs (they are the embedding); the `Q′` matrices are NOT, and each route constructs its
own, or the twisted-sector independence is weakened.

**Two-route sharing rule (frozen; B's core evidentiary mechanism).** The routes may share
ONLY the raw action pairs `[u_γ, v_γ]`, or a presentation together with the frozen
homomorphism into `SU(2)_L × SU(2)_R` that generates those pairs, plus coordinate conventions
and generic linear-algebra utilities. They may NOT share any characters, invariant dimensions,
multiplicity tables, projectors or projector ranks, reduced operators, symmetry-adapted bases,
spectral fixtures, or expected quotient answers.
Every quotient multiplicity match between the routes is then two disjoint derivations
agreeing, not one module reading its own output twice.

**Existing author-context code is post-run material only.** The `2I` character-theory
machinery already in the repository (the M8.2 scripts, quarantined for A's implementer, and
`m8_3_mass_reproducer.py`'s multiplicity formulas) is legal in B but only AFTER the frozen
outputs of both routes exist and are committed: as declared regression or adjudication
references. It may not supply expected answers, fixtures, or check values to either primary
route, and every such post-run use is declared in the method note.

**No `S³/2I` answer values in this document.** The untwisted `S³` rungs carry pinned
published values (§ 11), since those anchor gates against an external answer. The rung-3a/3b
adjudication case carries none here at all: its values are sealed (§ 4.1). The `S³/2I` multiplicity patterns appear nowhere in this
pre-registration or in any artifact shared between the routes, so that neither route can be
built toward a posted table.

**Row identity, scoped by rung.** For rungs 4-6 (`Γ = 2I`), sectors use the label-free
`(dimension, McKay-distance)` signature inherited from the A protocol, derived in-route, with
pairwise distinctness checked before any cross-route or reference comparison. For rungs 1
through 3b, records use that rung's frozen canonical sector identifier, or `null` where no
representation-sector decomposition is reported. The McKay signature is NOT assumed to
transfer outside `Γ = 2I`: a cyclic `Γ` can have irreducibles sharing a dimension, and its
McKay distance need not give a pairwise-unique signature.

## 3. Sequencing relative to M8.5-A (frozen ordering)

B has no context firewall of its own; A still does, and B's outputs are new `2I`-specific
answer-bearing artifacts. The durable ordering, enforced by commit timestamps:

1. The A protocol lands.
2. A's fresh implementation reaches its sealed commitment (the A protocol's § 9: source,
   environment, raw output, hashes, consulted-files manifest, method-note draft).
3. Only then may B generate and commit `2I`-specific prototype outputs (the § 4 rungs 4-7).
4. A unseals and adjudicates.
5. B's certification claim waits for A's core scalar outcome under the § 1 gate condition,
   not merely for A to close.

**Before step 2 completes, B work is restricted to:** untwisted `S³` machinery (rungs 1-2);
the rung-3a/3b adjudication gates, which are deliberately non-`2I` and therefore carry
no M8.2 target information; generic manifold numerics; synthetic or non-`2I` test cases;
performance and memory engineering carrying no M8.2 target information; and the § 6
engineering pilot. This pre-registration itself carries no answer values and may land before
step 2. If step 2 stalls, rungs 1, 2, 3a and 3b may still complete and be reported only under
§ 1's pre-`2I` ladder phrase (never the full-suite phrase); the `2I` rungs wait.

## 4. The certification ladder (frozen core)

Each rung is a pre-registered gate with a PASS line under § 8 discipline. **Rungs 1 through 6
(including both 3a and 3b) are MANDATORY for both routes; rung 7 is a joint comparison
deliverable. No rung may be declared inapplicable after the pilot.**

All gates are evaluated at UNIT RADIUS, in the exact form the sources state them (§ 5's
radius convention).

| # | Rung | Reference and evidence class |
| --- | --- | --- |
| 1 | untwisted `S³` scalar spectrum and multiplicities | published-value lookup; § 11 pins |
| 2 | untwisted `S³` one-forms: exact and coexact spectra and multiplicities, Hodge decomposition orthogonality, `H¹(S³) = 0`, and the implementation's chosen Weitzenböck convention verified as an identity | published-value lookup; § 11 pins. The Weitzenböck convention is stated and gated, not prescribed here |
| 3a | **precommitted scalar adjudication.** The complete scalar spectrum AND multiplicity sequence of the sealed non-`2I` adjudication case, through the certified band, reproduced by both routes | published-value adjudication (class 1), carrying no unseen-target claim (§ 0 ceiling). Case identity and action live in Packet I, and source, table, row, indexing map and values in Packet II (§ 4.1); this document carries only the opaque identifier and the two hashes |
| 3b | **precommitted one-form theorem evaluation.** Exact and coexact one-form multiplicities of the same sealed case, through the certified band, reproduced by both routes | **source-derived theorem evaluation (class 2), NOT a published-value lookup**: the verified source audit located no published one-form multiplicity lookup for the selected adjudication case. Procedure and its independence asymmetry frozen in § 4.2 |
| 4 | untwisted quotient scalar multiplicities on `S³/2I` | cross-route agreement; no values posted here |
| 5 | untwisted quotient one-form multiplicities and exact/coexact separation on `S³/2I` | cross-route agreement; no values posted here |
| 6 | BOTH nontrivial frozen coefficient sectors, `Q` and `Q′`, implemented and compared (coefficient `τ_σ = Sym²(σ)`, `Q` identified via `χ = 2cos θ` under the frozen embedding); failure or non-resolution in either sector is reported, and one sector cannot substitute for the other | cross-route agreement plus § 8 gates, never an author-context table |
| 7 | matched outputs from both routes sufficient for the measured trade-off table: accuracy vs cost vs implementation complexity at matched configurations (the § 6 common-band definition, never nominal resolution) | measured, not argued; route decision recorded with rationale; the losing prototype is kept as the M8.4 cross-check tool |

### 4.1 The sealed reference packet, and rung 3a

**Why sealing, not omission.** Omitting the integers from this document does not preserve the
adjudication case: a citation naming the source, table and row IS an answer-bearing pointer,
and either route could simply open the table before committing. So the identifying
information is sealed, not merely un-quoted.

**What this public document carries:** an opaque case identifier, the adjudication procedure,
the admissibility requirements below, and the **SHA-256 hashes of the two sealed reference packets**
(§ 11).

**TWO packets, because one cannot work.** A single packet would have to be opened to supply
the case input the prototypes need, while also staying sealed until those prototypes' outputs
are committed. That is circular. The seal is therefore split, and both hashes are published
before landing:

| Packet | Contents | Opens at |
| --- | --- | --- |
| **I, case-input** | ONE JSON file, exactly these fields: `case_id`, `family`, `parameters`, `generators`, `action_convention`, `format_version`. `generators` is a list of PAIRS `[u, v]` of unit quaternions denoting the isometry `q ↦ u q v`, and `action_convention` is `two_sided`. A left action is written `v = 1`; see the deck-group note in § 2 | step 2 below |
| **II, answer** | ONE JSON file, exactly these fields: `case_id`, `citation` (authors, title, venue, year, DOI, table, row), `indexing_map`, `reference_values`, `format_version` | step 4 below |

**Identifier scheme for the packet checks (frozen), so nothing drifts between documents.**
The packet build runs `S1` through `S6` (structural: free action, pair form, parameters
reproduce generators, closure order, freeness and central canonicalization, element-order
census) and `L1` through `L3` (leanness: unknown fields, prose, string length). These are
BUILD-TIME checks on the packet before sealing. `G1` and `G1b` in § 8 are RUN-TIME gates that
each route applies to the action it actually executes. They test overlapping properties at
different moments and are deliberately kept as separate identifiers; the method note, the
packet builder and its mutation output all use these names and no others.

**Packet contents are data and citations only (frozen).** Each packet is a single JSON file
with the fields listed above and no others, sealed with its SHA-256. Packets carry no working
notes, no rationale, no drafting or review history, no audit commentary, and **no free-text
field capable of holding any of those**. Every explanatory statement about the adjudication
lives in this pre-registration; every working record stays in local staging and is never
sealed, hashed, or published. The build is gated on a schema check that rejects unknown
fields, so leanness is enforced mechanically rather than by inspection, and the check is
mutation-tested against a packet carrying an extra field before either packet is sealed.

**Relation to the M8.5 reproduction standard.** M8.5-B makes no fresh-context or
unseen-target claim and therefore does not instantiate the standard's room procedure. It
inherits the claim ladder, commitment-before-reveal ordering, generator-based construction,
the canonicalize-then-hash rule, the two-commit adjudication structure, and the mutation-test
posture. Two adaptations are explicit. First, Packet I supplies the case parameter `q` for the
standalone theorem evaluator, while both primary routes construct the effective action and
derive its order from the supplied generator pair without using `q` as an expected-order
fixture; comparison with `q` is therefore a metadata-consistency check, while closure, the
effective census, central-kernel deduplication and freeness remain derived gates. Second, the
canonical Packet I and Packet II bytes are committed by publishing their SHA-256 hashes in
this pre-registration before either packet opens; Packet II itself remains outside the public
repository because committing its answer-bearing bytes would disclose the sealed case. Both
packets are reproducible from the retained author-side source, verified byte-for-byte against
the published hashes, so the commitment does not rest on a stored copy surviving.

**Frozen sequence.**

1. Pilot parameters, both routes' code, the standalone 3b evaluator's code, the generic
   tests, and the configuration are all committed.
2. Packet I is opened and its hash verified.
3. Both routes run on the now-disclosed case and commit their raw outputs.
4. Packet II is opened and its hash verified.
5. The adjudicator performs rung 3a against the revealed reference values.
6. The adjudicator runs the frozen 3b evaluator on the Packet I case input and compares its
   output against the ALREADY COMMITTED route outputs.

The 3b evaluator's code and its `Γ = 1` validation are pre-reveal (step 1); only its
case-specific execution happens at step 6, after route outputs are fixed. At step 5 the
adjudicator transcribes from Packet II into the comparison harness, which carries a
transcription mutation (perturb one transcribed cell; the comparison must go red), on the
pattern the A protocol uses, so the comparison is itself a check that can fail.

**Admissibility:** the case is non-`2I` and distinct from every member of the pilot tuning
set; it informs no pilot choice; and both routes must reproduce the COMPLETE scalar spectrum
and multiplicity sequence through the certified band, not a selected subset (§ 4's full-band
completeness rule applies here too).

**What the seal does and does not hide (deliberate boundary).** The case FAMILY is disclosed
(it is a cyclic quotient), because that fact is what justifies rung 3b's evidence class: the
verified source audit located no published one-form multiplicity lookup for a member of that
family, which is why 3b evaluates a theorem rather than looking a value up. The specific case
within the family, and every number attached to it, stays in the packets. The family is
infinite, so disclosing it does not by itself identify the specific adjudication case;
suppressing it would leave the class-2 label unexplained in the contract that has to justify
it.

**Honest label:** precommitted post-freeze published-value adjudication. Not held out from
the author, and carrying no unseen-target claim (§ 0 ceiling).

### 4.2 Rung 3b: precommitted theorem evaluation, and its independence asymmetry

The verified source audit located no published one-form multiplicity lookup for the selected
adjudication case; the pinned literature instead supplies a general multiplicity RULE that
must be evaluated. (That is a finding about the verified source set, not a proof of global
nonexistence.) This rung is therefore labeled **source-derived theorem evaluation** wherever
it is reported, never "known-answer". Its frozen procedure:

1. The published one-form multiplicity theorem is quoted precisely, with all its hypotheses,
   and the hypotheses are checked to hold for the adjudication case. **The theorem statement
   is GENERIC and may live in this public document (§ 11); the case-specific inputs it is
   evaluated on live in Packet I** (§ 4.1), since those inputs identify the case.
2. The source-to-protocol convention map is recorded (§ 11), covering at minimum the group
   action and its notation, the index conventions, the sector labeling, and real versus
   complex multiplicity counting.
3. A STANDALONE reference evaluator implements the theorem in exact or high-precision
   arithmetic. **It is separate from both prototype routes and never becomes a shared
   library**; no route may import it, and it may not import from either route.
4. The evaluator is first checked at `Γ = 1`, where it must recover the unit-`S³` one-form
   tower of rung 2. Failing that check voids the rung before it is applied. **This generic
   evaluator and its `Γ = 1` validation MAY be built during the pilot**, because that check
   uses only the rung-2 tower, which is in neither packet.
5. The case-specific input is drawn from Packet I and the evaluation performed only at
   step 6 of § 4.1's sequence, after both prototypes' outputs are committed.
6. Its transcription, group-action, and summation gates are each mutation-tested.

**Honest label:** precommitted post-freeze theorem evaluation, carrying no unseen-target
claim (§ 0 ceiling).

**Independence asymmetry (frozen; must be stated wherever 3b is reported).** Against route
(a) this is substantially method-independent, because route (a) is character-free by § 2's
construction rule while the published one-form rules are representation-theoretic. Against
route (b) it is a source-pinned theorem-CONSISTENCY check that likely overlaps route (b)'s
own character-averaging machinery. **3b is therefore not counted as a third independent
derivation of route (b)'s result**, and no report may present it as one.

**Full-band completeness (frozen; applies to rungs 4-6).** Cross-route comparison covers the
COMPLETE level-by-sector multiplicity array through the common certification band, INCLUDING
zero multiplicities. For every level, both routes report all irreducible signatures,
form/Hodge sectors, and frozen connection classes. Route (a)'s numerical eigenvalue clusters
must map to the allowed `S³` levels within the frozen tolerance. Missing expected clusters,
unmatched numerical clusters, omitted zero cells, and spurious levels below the common-band
ceiling are each STRUCTURAL FAILURES. Agreement on a selected subset is not agreement.

Notation is harmonized to the A protocol: the scalar tower is `V_n` with `j = n/2`, and every
spectrum is reported as the dimensionless `λ·R²`. Every reference formula, verbatim source
statement, theorem hypothesis, convention map, and use-warning lives in § 11; the runnable
gate uses the normalized protocol formula, never an unexplained verbatim expression. The
`S³/2I` rungs (4-6) stay value-free (§ 2), and rung 3a's integers stay out of this document
(§ 4.1). Route (a) must additionally measure and report the certified spectrum's degeneracy
splitting at every frozen resolution; route (b) must report in-band cutoff stability for every
observable. Both are deliverables, not spot checks.

## 5. Frozen conventions

Inherited from the A protocol, not re-opened: the deck action (`2I ⊂ SU(2)` by left
multiplication); the scalar tower and eigenvalue map; the three flat-connection classes with
`τ_σ = Sym²(σ)` as declared contract input; `Q` identified, never declared, via
`χ_Q(g) = 2cos θ(g)` under the frozen embedding, `Q′` the unique other 2-dimensional
irreducible; label-free `(dim, McKay-distance)` row identity; complex-representation
multiplicity counting with the convention stated.

Frozen here for B (the one-form side A never needed): the induced one-form action is
pullback on the geometric index plus the coefficient action on any twist index; restriction
multiplicity and quotient multiplicity (`dim Hom_2I`) are BOTH computed and reported,
labeled, wherever they differ; real vs complex dimension reporting is fixed per artifact and
stated in the schema; measured multiplicities (projector ranks, degeneracy counts) obey an
integer-nearness rule, ported from the A protocol § 5.4: every quantity interpreted as a
multiplicity must lie within its stated tolerance of a nonnegative integer, any violation is
a structural failure, and rounding happens only after that gate passes.

**Radius convention (frozen; a protocol derivation, not a cited claim).** All benchmark gates
are evaluated on the UNIT sphere, in the exact form the sources verify at `R = 1`. For
reporting on `S³(R)`, write `g_R = R² g_1`. Under this constant metric rescaling the scalar
and Hodge-Laplacian eigenvalues scale by `R^-2`, so the unit-radius eigenvalues equal the
reported dimensionless quantities `λ·R²`. **This is a protocol-level homothety derivation and
is labeled as such; it is not a verbatim claim attributed to the spectrum references.** (One
published statement of the rescaling does exist and is pinned in § 11, but its form-degree
range excludes the scalar case, so the protocol carries the derivation for both rather than
citing it unevenly.)

**Weitzenböck convention (frozen; quoted identity plus a labeled substitution).** Under the
frozen convention

```text
Delta_Hodge  =  d delta + delta d  =  nabla* nabla  +  Ric,      nabla* nabla = -trace(nabla^2)
```

**the identity above is QUOTED** (Petersen Thm 9.4.1, § 11.5, with four independent
corroborators agreeing on every sign). **The curvature substitution is ASSEMBLED, not
quoted**: the source audit located no single page stating `Ric = 2 g` for the unit `S³`, so
both it and the radius form `Ric = 2 R^-2 g` are DERIVED at protocol level from constant
sectional curvature `R^-2` together with `Ric = (n−1) K g` in dimension three. The
distinction between the quoted identity and the assembled substitution is recorded wherever
the gate is reported, and § 11.5 carries it in the same terms.

**Core output schema (frozen now, not pilot-selected).** Both routes emit these fields with
these meanings, so the two result surfaces are comparable by construction rather than by
later reconciliation. The pilot may add route-specific diagnostic fields; it may not change,
rename, or drop these, and their meanings are fixed at freeze.

```text
schema_version
route
run_id
configuration_id
arena_case_id
group_order
rung
sector_id
form_degree
hodge_sector
connection_class
harmonic_level
eigenvalue_R2
sector_signature {dimension, mckay_distance}
restriction_multiplicity
quotient_multiplicity
measured_rank
integer_nearness_margin
eigenvalue_residual
subspace_residual
degeneracy_splitting
convergence_statistic
gate_results
wall_time
peak_memory
```

**Null semantics (frozen).** `null` means NOT APPLICABLE under that rung's definition;
numeric zero means a computed mathematical zero; omission of a required record is a
STRUCTURAL FAILURE. The distinction is load-bearing because `sector_signature`,
`restriction_multiplicity`, and several residual fields are not naturally applicable to every
rung: outside `Γ = 2I` the McKay signature does not apply (§ 2), and those fields carry
`null` rather than a fabricated value.

Every level-by-sector cell through the common certification band exists as a record,
including zero-multiplicity cells (§ 4's full-band completeness rule); an absent cell is a
structural failure, not an implicit zero.

## 6. Parameters: the engineering pilot, then the freeze

Numerical parameters are not invented from memory and not left open. Two stages:

**Stage 1, engineering pilot (non-evidentiary).** Restricted to known untwisted `S³`
benchmarks and a PREREGISTERED TUNING SET of non-`2I` synthetic quotient tests, including
the preregistered cyclic cases and `2T`. Its sole purpose is to select the parameter classes listed below: the route
architecture schemes, the eigensolver and rank-algorithm family, the ladders and cutoffs, the
common certification band and matched-accuracy rule, every tolerance and statistic, the
subspace-comparison criterion, the trade-off rubric, the repeat-run count, the resource cap,
and any route-specific diagnostic fields. No `2I`-specific target or coefficient-sector
result may be examined during the pilot.

**The sealed adjudication case is excluded from the pilot (frozen).** The adjudication case
(distinct from every member of the tuning set) does not participate in architecture selection,
tolerance selection, or any other pilot choice. BOTH packets are sealed and hashed before
landing (§ 4.1); the § 4.1 sequence governs when each opens; and both routes must pass rungs
3a and 3b before either touches `S³/2I`. A case that has informed any pilot choice is void as
an adjudication case, and its rungs with it.

Specifically off-limits to the pilot: both packets (Packet I opens at § 4.1 step 2, Packet II
at step 4, and neither during the pilot) and any evaluation of the reference evaluator on the
case-specific input. The GENERIC evaluator may be built and checked at `Γ = 1` during the
pilot, since that check uses only the rung-2 unit-`S³` tower, which is in neither packet.

Because the routes are parameterized by a frozen finite `Γ` (§ 2), the pilot selects
architecture and tolerances that must work for BOTH the tuning-set groups and, later, `2I`,
without ever seeing the adjudication case.

Pilot bookkeeping is author-side and is not published; the pilot's existence and its
selected values are reported in the method note, and its selected values are § 6 below.

**Stage 2, the pre-freeze parameter lock (COMPLETE, 2026-08-01).** The pilot ran on the tuning set below and
its measurements select every value in the tables that follow. Each value is backed by a
recorded measurement; the pilot log and its JSON outputs are author-side evidence and are
not published. Two parameters the drafted table did not anticipate are
frozen here as well, because the pilot found they govern accuracy and would otherwise be
free dials at run time: the adequacy-test accuracy tolerance, and route (a)'s polynomial
augmentation degree with its stencil size. Only with these frozen may either route run on
`S³/2I`.

### 6.1 Tuning set and architecture

| Parameter | Frozen value |
| --- | --- |
| pilot quotient tuning set (non-`2I`, preregistered) | `L(2,1)`, `L(3,1)`, `L(4,1)` (cyclic, homogeneous), `L(7,2)` (cyclic, INHOMOGENEOUS, exercising a genuinely two-sided action), and `2T` (nonabelian, order 24, left). The three homogeneous cyclic cases reproduce published scalar multiplicities exactly through `k = 9`; `L(7,2)` supplies the inhomogeneous two-sided cross-route test; `2T` carries the nonabelian measurements. The sealed adjudication case is excluded by construction and was never run |
| sealed adjudication case: opaque identifier + Packet I and Packet II SHA-256 hashes (identity and values stay sealed, § 4.1) | `M85B-ADJ-01` / `cc8c38f0c6819daa1125b5cfea955bfc6b31d160859e0d226f47aad9399fbecc` / `b89c4a1446fd88df156d2f7505dc3449dcb3f7d4eb9961364d3dc757862ebdf7` |
| source-to-protocol convention maps (generic; no case information) | FILLED: § 11.1-11.5, 11.8, 11.9 |
| route (a) discretization + identification scheme (embedding grid vs intrinsic charts; constraint construction) | `Γ`-orbit point cloud: seeds Riesz-relaxed BEFORE orbits are generated, then RBF-FD with polyharmonic spline and polynomial augmentation, with the orbit reduction accumulated during assembly. Orbits are generated by the ACTION PAIRS as `4 × 4` rotations, never by quaternion left multiplication. The scalar reduction is plain orbit summation. The one-form reduction is `Ad`-EQUIVARIANT, not componentwise: with `f(γ · x) = Ad(v_γ⁻¹) f(x)` (§ 2), the accumulation carries the `3 × 3` rotation, `R[aM+m, cM+l] = Σ_{j ∈ orbit l} Σ_b A_ab[i_m, j] · [Ad(v_γⱼ⁻¹)]_bc`. The operators themselves are unchanged, since the RBF-FD weights, the left-invariant frame and the Hodge assembly are properties of `S³` and not of `Γ`. A structured tensor grid is not `Γ`-closed for a nonabelian group (90 of 240 sampled images landed off the levels), so the orbit cloud is forced rather than preferred; the relaxation converts a flat error into a converging one |
| route (a) RBF-FD operator: polyharmonic exponent, polynomial augmentation degree, stencil | polyharmonic spline exponent `m = 7`, polynomial degree 4, stencil 110 (roughly twice the polynomial-space dimension) is the SOLE certification configuration under `m8_5b-v1`. The EXPONENT is frozen alongside the degree and stencil: it governs the operator as directly as they do, and `m = 5` versus `m = 7` measurably changes every gate quantity. Measured at 11520 nodes, `m = 7` is better at both banded levels on the level error (6.33e-5 and 3.41e-4 against 9.70e-5 and 3.59e-4), the worst-member error and the spread, and it resolves `λ = 80` two ladder rungs earlier. Degree 3 with stencil 60 is retained as a non-certifying engineering fallback: activating it requires a dated § 12 addendum and a new `configuration_id`, and the original `m8_5b-v1` result remains filed. No pair may be substituted at run time |
| route (b) basis + projector assembly scheme | harmonic expansion on `S³` with multiplicities predicted by character averaging over the ACTION PAIRS, the group closed independently from the raw pairs. **Frozen decomposition and index map**, with the left factor carrying `u` and the right factor `v`, which is forced by the § 2 pullback law since the frame is carried by `Ad(v)` alone: scalars `L²(S³) = ⊕_n V_n(u) ⊗ V_n(v)`, invariant dimension `(1/\|Γ\|) Σ_γ χ_n(u_γ) χ_n(v_γ)`; one-forms `Ω¹ = ⊕_n V_n(u) ⊗ [V_n ⊗ V_2](v)` with `V_n ⊗ V_2 = V_{n+2} + V_n + V_{n−2}` and invariant dimension `(1/\|Γ\|) Σ_γ χ_n(u_γ) χ_m(v_γ)` for each `m`. **Sector, eigenvalue and RANGE, frozen explicitly** so the towers cannot be interchanged while still totalling correctly. `n` labels the SCALAR factor `V_n(u)` throughout, never the resulting representation and never the eigenvalue index: `m = n`, exact, `λ = n(n+2)`, for `n ≥ 1`; `m = n+2`, coexact, `λ = (n+2)²`, for `n ≥ 0`; `m = n−2`, coexact, `λ = n²`, for `n ≥ 2`. A summand outside its range is a COMPUTED ZERO CELL wherever the schema requires that sector, never a silent omission |
| eigensolver / rank-algorithm family | route (a): dense non-symmetric eigensolve of the ORBIT-REDUCED operator, never of the full cloud, and never symmetrized. Route (a) rank measurements: nullity from singular values. Route (b): closed form, no solver. Measured basis: the RBF-FD operator has asymmetry ratio 1.17 and symmetrizing it produces hundreds of spurious modes; the reduced solve is about 200 times faster than the full one |

**Three structural rules, frozen as rules rather than left to the implementer.** Each turns a
correct operator into a wrong answer silently, and each is exactly the kind of step a later
tidy-up reintroduces.

1. **No symmetrization** of the RBF-FD operator, for the reason measured above.
2. **No `.real` on eigenvectors.** The reduced operators are non-symmetric, so a degenerate
   real eigenvalue can return eigenvectors with substantial imaginary content (measured up to
   0.13). Discarding it destroys the subspace: run that way, the first subspace comparison
   reported the gradient image and the exact sector as very nearly ORTHOGONAL, sine of the
   maximum principal angle 1.00, where they in fact coincide. The rule: form the real span of
   `[Re V, Im V]` and orthonormalize with an SVD RANK CUT, or a column-pivoted QR. **An
   unpivoted `qr` is not rank-revealing.** For a complex-conjugate pair `[Re V, Im V]` carries
   up to twice the true rank, and an unpivoted factorization returns the excess as spurious
   orthonormal directions lying outside the subspace, which drives any principal-angle
   measurement to 0.9999 where the correct value is 7.2e-3. The factorization is named because
   the failure is silent and the wrong routine satisfies a generic instruction.
3. **Order is quoted in `h`, never per node doubling.** In three dimensions a node doubling
   shrinks `h` only by `2^(-1/3)`, so a raw per-doubling ratio understates the order
   threefold and makes a converging method read as a plateau.

### 6.1b Two-sided effective action: frozen conventions and mandatory gates

**Interface.** A deck element is the action pair `[u_γ, v_γ]` acting as `x ↦ u_γ x v_γ`
(§ 2). `Γ = 2I` is `[γ, 1]`; any left action is `v = 1`. The pair representation is
two-to-one onto the effective action, `[u, v]` and `[−u, −v]` inducing the same isometry,
so the EFFECTIVE order is the gated quantity. An inhomogeneous lens space is not a left
action, so this interface is required and not optional.

**Action conventions.** Orbits are generated by the action pairs as `4 × 4` rotations,
never by quaternion left multiplication. The scalar reduction is plain orbit summation.
The one-form reduction is `Ad`-equivariant.

**Coefficient law (frozen).** For `ω = f_j σ^j` in the left-invariant coframe,

```text
f(gamma . x)  =  Ad(v_gamma^-1) f(x)
```

The inverse is load-bearing. `Ad(v)` and `Ad(v⁻¹) = Ad(v)ᵀ` are orthogonal with the same
`SU(2)` character, since `χ_n(v⁻¹) = χ_n(v)`, so cross-route agreement and every multiplicity
comparison are insensitive to substituting one for another. That is why the manufactured gate
below is mandatory rather than advisory. There are two candidates here, not three: `Ad` is
orthogonal, so `Ad(v)ᵀ` IS `Ad(v⁻¹)`. The distinguishable hazard at assembly is the
accumulation ORDER, which is gated separately.

**Route (b) summand and index map (frozen).** Left factor carries `u`, right factor `v`,
forced by the coefficient law since the frame is carried by `Ad(v)` alone. `n` labels the
SCALAR factor `V_n(u)` throughout, never the resulting representation and never the
eigenvalue index.

| summand | sector | eigenvalue | admissible range |
| --- | --- | --- | --- |
| `m = n` | exact | `λ = n(n+2)` | `n ≥ 1` |
| `m = n+2` | coexact | `λ = (n+2)²` | `n ≥ 0` |
| `m = n−2` | coexact | `λ = n²` | `n ≥ 2` |

Invariant dimension of each piece is `(1/|Γ|) Σ_γ χ_n(u_γ) χ_m(v_γ)`. A summand outside its
range is a COMPUTED ZERO CELL wherever the schema requires that sector, never a silent
omission.

**Output-band completeness (frozen).** The harmonic cutoff is imposed on the requested
OUTPUT band after all Clebsch-Gordan index shifts. Truncating every parent summand at the
same raw harmonic index is forbidden: every summand capable of contributing to a reported
level through `n_max` must be included, **even when its parent harmonic index lies above
`n_max`**. A coexact level `λ = m²` receives contributions from `n = m−2` and from `n = m`,
whose sum is symmetric under `u ↔ v` while neither is alone, so coexact levels are complete
only for `m ≤ n_max` and levels above that are never certified.

**Mandatory gates.** Each is required before either route runs on `S³/2I`:

| Gate | Requirement |
| --- | --- |
| left-action regression | the action-pair path reproduces the quaternion path with identical multiplicities and identical gate outcomes, numerics under G8 tolerances. Bit-for-bit is not required: the clouds differ in node ORDER, so the reduced operators are permutation-similar rather than equal |
| manufactured pullback | a pointwise test of the coefficient law that rejects, at minimum, omission of the rotation, `Ad(v)`, `Ad(u⁻¹)`, and the left/right entry swap. A `v = 1` configuration cannot expose the swap and does not satisfy this gate |
| accumulation order | the assembled reduction is compared against the full operator on a field extended by the coefficient law. A left action cannot discriminate the orders, since `Ad(1) = I`, so this gate requires a two-sided configuration |
| right-only equivalence | inversion `x ↦ x⁻¹` conjugates `F_{u,v}` to `F_{v⁻¹,u⁻¹}`, so `[1, v]` and `[v, 1]` are isometric quotients and must agree per sector, scalar, exact and coexact, with zero cells emitted |
| effective action and freeness | G1b (§ 8), gated independently by both routes on the action each executes |
| cross-route | agreement on an inhomogeneous tuning case through the complete frozen band, including zero cells |

The two-sided generalization and its frozen conventions were validated before freeze by the
preregistered known-answer, left-action regression, manufactured pullback, left/right
equivalence, and inhomogeneous cross-route tests. The complete results, mutations,
configurations, raw outputs and hashes are preserved in the author-side pre-freeze
validation record, SHA-256
`62b95e5a473db9b1436fcbb1ed38b622f3d80ae23b48e6bde2d9d20d3eab8117`. That record is evidence
for the frozen choices; it is not part of the normative execution contract.

### 6.2 The common certification band

| Parameter | Frozen value |
| --- | --- |
| common certification band: canonical harmonic ceiling | the **second nonzero untwisted scalar invariant level** defines the frozen harmonic ceiling `n_max`. Rungs 4 through 6 include every required scalar, exact, coexact and coefficient-sector record through `n_max`, each sector using its OWN frozen eigenvalue map |
| cross-route matched-accuracy rule | both routes must certify the same `n_max` under the same accuracy threshold, with matched configurations defined by that ceiling and never by nominal resolution |

**`n_max` and `λR²_max` are not interchangeable across sectors.** The scalar, exact one-form
and coexact one-form towers carry different eigenvalue maps (`n(n+2)`, `ℓ(ℓ+2)` and `m²`
respectively), so one eigenvalue ceiling does not describe the same content in all three, and
"the first two nonzero invariant levels" does not by itself fix the `Q`, `Q′`, exact and
coexact records. The harmonic ceiling is therefore the canonical coordinate, the § 5 schema
already carries it as `harmonic_level`, and each sector's eigenvalue ceiling is derived from
`n_max` through that sector's own map. The protocol never states the two as equivalent.

**Required band versus achieved band (frozen).** The REQUIRED certification band is the closed
band from the quotient's zero level through `n_max`. **It does not shrink in response to a
numerical failure.** The § 6.3 stopping rule reports the ACHIEVED consecutive resolved band
together with the first unresolved level; any unresolved or failed level inside the required
band fails the affected mandatory rung. Nothing in the stopping rule may redefine the required
band downward to produce a smaller passing one.

The ceiling is stated as a rule rather than a number because its numeric extent depends on the
quotient under test, and the pilot is barred from examining the target group's spectrum. The
rule is evaluated when the run happens.

**Why the second level and not the third, measured.** At 11520 nodes, at the frozen
certification configuration, route (a) delivers the first nonzero scalar invariant level
(`λ = 48`) at level error 6.33e-5 with the correct multiplicity 7, and the second (`λ = 80`)
at level error 3.41e-4 with the correct multiplicity 9. Both multiplicities match route (b)'s
independent character-averaging prediction. At the third level the multiplicity is wrong, 41
modes found against 26 expected, so the third level is not within demonstrated reach and
`n_max` does not extend to it.

**One-form rungs resolve later than scalar rungs.** At the frozen certification configuration
the one-form clusters do not carry correct multiplicities until 2880 nodes, where the scalar
problem is already accurate at 1440. Route (a) implementations must therefore assemble the
one-form operators SPARSELY: the pilot's dense prototype reaches 5760 nodes, while the top
ladder rung at 11520 nodes would need roughly 9.5 GB for a dense `3N × 3N` operator. The
frozen ladder is unchanged; the assembly is required to be sparse.

### 6.3 Tolerances and statistics

Tolerances are separated by the KIND of quantity they govern. A single relative threshold
cannot govern a zero eigenvalue, an integer multiplicity, a residual and a subspace angle, and
G3 is not implementable while it points at one.

| Parameter | Frozen value | Measured margin |
| --- | --- | --- |
| accuracy tolerance, NONZERO eigenvalues (adequacy test) | `\|λ_level − λ_ref\| / λ_ref ≤ 1e-3`, where `λ_level` is the MEAN of the assigned cluster. The certified quantity is the LEVEL VALUE; the scatter of individual members within the cluster is governed separately by the degeneracy-spread gate below, and the two are never conflated | measured at 11520 nodes, frozen configuration: 6.33e-5 at `λ = 48` and 3.41e-4 at `λ = 80`, between threefold and fifteenfold below the threshold. Stated explicitly because the two readings disagree: the worst INDIVIDUAL member at the same rung is 5.58e-4 and 1.46e-3, so a per-member reading of this row would fail its own ladder at `λ = 80` while the level values pass comfortably |
| accuracy tolerance, ZERO modes | `\|λ_num\| ≤ 1e-9` absolute; relative error is undefined here and is never used | measured 3.19e-13, 1.71e-13, 6.68e-14 across the ladder at the certification configuration, over three orders of margin |
| multiplicities | exact integer comparison after the integer-nearness gate below; no relative tolerance applies | |
| degeneracy-clustering threshold (ASSIGNMENT ONLY) | absolute window of 0.35 in `λR²` around a candidate level, deciding cluster MEMBERSHIP only. It is not a spread gate | the window cleanly separated every level used in the cyclic and nonabelian runs |
| degeneracy-spread statistic and threshold | for every NONZERO cluster, `spread := (λ_max − λ_min) / \|λ_mean\|`; threshold 5e-2. No relative spread statistic applies to the zero cluster, where it is undefined | measured on the one-form ladder at 1.37e-2, 3.62e-3, 1.49e-3, and on the scalar band at 11520 nodes at 1.08e-3 (`λ = 48`) and 1.97e-3 (`λ = 80`), between twentyfive- and fortyfold below the threshold at the finest rungs measured |
| exact/coexact orthogonality residual | `residual := ‖UᵀV‖₂ = σ_max(UᵀV)`, where the columns of `U` and `V` are orthonormal bases of the exact and coexact clusters. This is the cosine of the smallest principal angle between them. The largest individual entry of `UᵀV` is basis-dependent and is NOT used. Threshold 1.5e-1 | measured 1.995e-1 and 1.053e-1 at the two rungs where both clusters resolve, converging at order `q ≈ 2.8` in `h`. The threshold is met at 5760 nodes with about a 1.4-fold margin and is NOT met at 2880 (1.995e-1), so this gate, not the eigenvalue gate, is the binding resolution constraint on the one-form rungs |
| Weitzenböck residual | `residual := max_j \|(Δ_H df)_j − λ_f (df)_j\| / max_j \|(df)_j\|` for a scalar eigenfunction `f` of eigenvalue `λ_f`; threshold 1e-7 on this normalized dimensionless residual | measured 6.07e-11, 1.78e-10, 1.51e-10 at harmonic degrees 1 and 2. It GROWS slowly with node count because it is limited by stencil conditioning rather than truncation: this is a round-off ceiling, not a converging quantity, and it is never judged by the convergence statistic |
| eigensolver tolerance | dense direct eigensolve of the reduced operator throughout the frozen ladder; no iterative tolerance is exercised. If a future band pushes the reduced dimension above 4000 an iterative solver becomes necessary, and its tolerance is set by a dated § 12 addendum, never improvised at run time | the reduced dimension stays at most `3 x 480 = 1440` for a deck group of order 24, and the largest such problem solves densely in seconds |
| projector-rank tolerance | 1e-8 relative to the largest singular value | measured kernel separation 1.7e15 |
| multiplicity-rounding (integer-nearness) tolerance | 1e-6 absolute from the nearest nonnegative integer | measured deviations 7e-14 (route b) and comparable for route (a) ranks, eight orders of margin |
| subspace-comparison metric + tolerance | the SINE OF THE MAXIMUM PRINCIPAL ANGLE between the two subspaces, computed as `sqrt(1 − σ_min(UᵀV)²)` where `U` and `V` are orthonormal bases obtained by a RANK-REVEALING factorization (§ 6.1 rule 2). Tolerance 5e-2 | gradient image against the Hodge exact sector at `λ = 48`, both rank 7: 2.1436e-2 and 7.1763e-3 at 2880 and 5760 nodes, converging at order `q ≈ 4.7` in `h`, with 2.3-fold and 7.0-fold margin. At 1440 nodes the scalar cluster has rank 4 against the sector's 6, so the level is UNRESOLVED there and the comparison is not read: the stopping rule's multiplicity condition fires exactly as designed |
| convergence statistic and acceptance rule | monotone decrease of the worst banded error over the final three rungs, order reported in the mesh parameter `h` as `q = 3 log2(ratio)` per node doubling. Applies to discretization-limited quantities only, never to the Weitzenböck residual | `p=3` gives `q ≈ 3.1, 4.0, 4.5`; `p=4` gives `q ≈ 6.8, 4.2, 3.2`; `p=1` is flat at `q ≈ 0.3` to `0.5`, which was the plateau |
| stopping rule for unresolved levels | a NONZERO level is RESOLVED only if all three hold: its cluster is non-empty within the 0.35 assignment window, its multiplicity equals the cross-route predicted multiplicity, and its spread is within the frozen degeneracy-spread threshold. The ZERO level is RESOLVED if its cluster is non-empty, its multiplicity is correct, and every member satisfies `\|λ_num\| ≤ 1e-9`; the relative spread statistic is not applied to it. The ACHIEVED resolved band terminates at the last consecutively resolved level. The first unresolved level and everything above it is reported as UNRESOLVED, never as absent. The achieved band never redefines the § 6.2 required band | multiplicity mismatch is the detector that fired in every failure the pilot saw: scalar `λ = 168` returned 41 modes against 26, and the curl build returned 4 against 5 at one-form `λ = 36`. In both the eigenvalue looked plausible and only the count exposed the failure |
| repeat-run count | 3 | sufficient to demonstrate the reproducibility gate without multiplying the ladder cost |

**The additional residual thresholds are backed by supplemental pilot measurements at the
frozen certification configuration, and every figure quoted beside a threshold traces to a
committed pilot data file.** The Weitzenböck residual is conditioning-limited and is judged
against its absolute threshold rather than the discretization-convergence rule.

**The orthogonality test is not tautological.** The coexact cluster is taken from the Hodge
operator's own spectrum at a perfect-square level, not as the complement of the gradient
image. That is well defined because the two towers never collide: `n(n+2) = (n+1)² − 1` is
never a perfect square. A same-sector control returns 1.000 at every rung, so the statistic
demonstrably registers non-orthogonality rather than being structurally forced small.

### 6.4 Ladders, caps, schema

| Parameter | Frozen value |
| --- | --- |
| resolution ladder (route a) | seed counts 60, 120, 240, 480, giving 1440, 2880, 5760, 11520 nodes for a deck group of order 24, and `\|Γ\| × seeds` in general |
| harmonic cutoff + oversampling (route b) | band limit at least the common certification band, with the run reporting cost at twice that band to demonstrate headroom. No oversampling parameter is frozen because none exists: route (b) evaluates character averages as exact finite sums over the group, with no quadrature to oversample |
| trade-off metrics and route-selection rubric (incl. what "implementation complexity" measures) | as written in the pilot's route-selection rubric, including the four-count definition of implementation complexity and the lexicographic ordering of accuracy, cost, complexity |
| resource cap / timeout | 300 s per ladder rung, fail loud on exceed |
| route-specific diagnostic fields (additions only; the § 5 core schema is already frozen) | route (a): per-rung degeneracy spread, stencil condition estimate, node-separation statistics. Route (b): exact in-band cutoff stability, comparing the certification cutoff against twice that cutoff, plus the integer deviation of every predicted dimension |
| output schema version string | `m8_5b-v1` |

**Matched configurations (frozen definition):** configurations are matched by a COMMON
certification band: the same frozen `n_max`, with each sector evaluated under its own
frozen eigenvalue map and the same applicable accuracy threshold on both routes. Nominal resolution (grid size against harmonic
cutoff) is never the matching variable: comparing those directly is meaningless unless both
routes certify the same spectral content.

**Route-selection rule (frozen):** both routes must first pass every required gate over the
common certification band. Among passing routes, selection follows the preregistered
accuracy/cost/complexity rubric; no failed route may win through lower cost.

**No-widening rule (frozen, ported from the A protocol):** tolerances and convergence gates
may not be relaxed after quotient results are visible. Any change requires preserving the
original run and filing a dated rerun under a new configuration identifier; the original
verdict stands as filed. The deliverable reports stability of every gate under increased
working precision or an equivalent numerical-conditioning check.

## 7. The coexact module and the combined-record state rules

The coexact one-form entry rule (`d` for `d ≥ 2`, `2` for `d = 0`, `3` for `d = 1`; source
convention: level `m`, eigenvalue `m²/R²`) is ASSERTED, with no independently derived target
anywhere in the record. B's operators can produce numerical support or contradiction at the
operator level; the A protocol's optional § 6 module can produce a representation-theoretic
derivation. Neither owns the verdict by fiat. The combined record is classified by the frozen
state rules below, which are consistent with, and do not amend, the A protocol's § 6 standing
rule (that protocol still governs A's own module verdict; these rules govern the combined
record and B's certification consequence).

**Outcome vocabulary (frozen).** A-side outcomes: `structural derivation`,
`numerical agreement`, `contradiction`, `unresolved`, `not run` (A's module is optional).
B-side outcomes: `numerical agreement`, `contradiction`, `unresolved`, `not run`.

**Overlay rule (frozen, applies regardless of the record classification below).** If B's
coexact work was NOT RUN, B's definition of done is not met and no B certification claim is
available: B is certifying one-form infrastructure, so a missing coexact module means B is
incomplete, not that the coexact question is merely open. The record is still classified by
the rules below, so both facts are reported.

**State rules (frozen; exhaustive over every combination, applied in order with the first
matching rule deciding).** No outcome is adjudicated ad hoc.

| # | State | Combined record | Certification consequence |
| --- | --- | --- | --- |
| 1 | both modules contradict | contradicted by both modules | blocked pending rule or convention adjudication |
| 2 | one side decisive-positive (derivation or agreement), the other contradicts | conflict; not resolved | blocked pending adjudication |
| 3 | one side contradicts, the other unresolved or not run | contradicted by one module; combined record unresolved | blocked |
| 4 | A structural derivation, B numerical agreement | structurally derived and numerically reproduced | not blocked by this rule |
| 5 | A structural derivation, B unresolved or not run | A-side structural standing retained (per the A protocol's § 6) | B's operator reproduction unresolved; B certification not earned |
| 6 | both modules agree numerically | ASSERTED, numerically corroborated by both modules; standing unchanged | not blocked by this rule |
| 7 | exactly one module agrees numerically, the other unresolved or not run | ASSERTED, numerically supported by the named module only; standing unchanged | if B is the unresolved side, B certification not earned |
| 8 | neither module decisive | not resolved | not earned |

Rules 4 and 6 are disjoint by construction (rule 4 requires A's derivation, rule 6 requires
A's numerical agreement), so a derivation never reads as mere corroboration. A finite B match
never becomes a derivation, and a valid B contradiction is never dismissed because A produced
an argument: rules 1 through 3 block certification whichever side contradicts.

"Contradiction" requires the disagreement to persist over the final three converged
configurations of the frozen ladder with numerical uncertainty too small to bridge it, under
exact integer comparison after the § 5 integer-nearness gate and a verified convention map
between level indexings. Anything weaker is `unresolved`, not `contradiction`.

## 8. Gates, mutation requirements, and uncertainty

Every line either route prints as PASS is demonstrated to fail under a deliberate defect
before the deliverable ships, via a RUNNABLE mutation harness (a `--mutation-tests` mode or
equivalent) with enforced coverage and a nonzero exit if any mutation fails to redden its
target or any gate is uncovered. Manual attestation does not qualify. Precedent: the
`m7_trivial_ok` defect and its replacement by falsifiable gates.

Minimum gate families (each rung's PASS lines join the requirement; the pilot may add gates,
never remove these):

| # | Gate family |
| --- | --- |
| G1 | group construction, EVERY quotient case: the constructed group is closed and has the frozen case order, and its conjugacy classes partition that order. For the `2I` rungs specifically the order is 120 (the A protocol's G1-G2) |
| G1b | EFFECTIVE ACTION, every quotient case, gated independently by BOTH routes on the action they actually execute: the action pairs reduce to an effective group of exactly the frozen quotient order, with the two-to-one central identification `[u, v] ∼ [−u, −v]` accounted for and not double-counted; every nonidentity effective element acts FREELY, verified as having NO FIXED POINT on `S³`, that is no eigenvalue `+1` of its `4 × 4` rotation, rather than inferred from the declared parameters; and the element-order census matches the frozen group's census. **Canonicalization is applied to the whole generated group, not to the input generators**: closure, uniqueness and census run on the `4 × 4` rotations, which are identical for `[u, v]` and `[−u, −v]`, so two words differing only by the central kernel cannot inflate the apparent order. Verified: adding `[−u, −v]` alongside `[u, v]` leaves the closure order unchanged. A malformed pair representation that silently yields an orbifold, or double-counts the deck group, passes an abstract closure check and fails this one |
| G1b evidence-strength note | Packet I supplies the adjudication case parameters, including the declared group order, because they are required case inputs for the frozen theorem evaluation. Agreement between that declared order and the order generated independently by each route is therefore a falsifiable METADATA-CONSISTENCY check, not an order derivation from withheld information. Closure, the effective-element census, central-kernel deduplication and freeness remain independently derived gates. **The case parameters are consumed by the standalone theorem evaluator; the two primary routes construct the effective action and its order by closing the raw action pairs, and may not use the declared order as an expected-order fixture.** |
| G2 | untwisted eigenvalue gates (rungs 1-2): computed eigenvalues match the pinned references, § 11.1 for the rung-1 scalar spectrum and § 11.2 for the rung-2 one-form towers, within the frozen tolerance, mutation-tested against perturbed references |
| G3 | rung-2 structure gates, each separately mutation-tested: exact one-form multiplicities match the § 11.2 exact tower `λ = ℓ(ℓ+2)`, multiplicity `(ℓ+1)²`; coexact one-form multiplicities match the § 11.2 coexact tower `λ = m²`, multiplicity `2(m²−1)`; the gradient image and the nonconstant exact sector agree under the frozen subspace-comparison criterion (§ 6.3), evaluated only at levels the stopping rule reports as RESOLVED; the harmonic one-form dimension equals zero (§ 11.4); the exact/coexact orthogonality residual is within the frozen 1.5e-1 threshold under the § 6.3 definition (theorem pinned at § 11.3); the Weitzenböck residual is within the frozen 1e-7 absolute threshold under the § 6.3 definition and the declared sign and curvature convention (identity quoted at § 11.5, curvature substitution assembled per § 5) |
| G4 | multiplicity integer-nearness (§ 5) for every measured and predicted multiplicity |
| G5a | precommitted scalar adjudication (rung 3a): BOTH packet hashes verify at their reveal steps; the complete scalar spectrum and multiplicity sequence through the certified band matches Packet II's reference under its recorded indexing map; the comparison harness carries a transcription mutation; the § 4.1 step order is recorded, showing Packet II opened only after both routes' raw outputs were committed |
| G5b | precommitted one-form theorem evaluation (rung 3b): the evaluator's code and its `Γ = 1` recovery of the rung-2 unit-`S³` tower are committed at § 4.1 step 1 (precondition); its case values, computed from Packet I at step 6, match both routes' already-committed outputs through the certified band; its transcription, group-action, and summation gates are each mutation-tested; the class-2 label, the § 0 claim ceiling, and the route-(b) overlap asymmetry are recorded with the result |
| G6 | cross-route agreement gates for rungs 4-6: exact integer equality after integer-nearness, label-free row matching with signature distinctness checked, AND full-band completeness per § 4 (every level-by-sector cell present through the common band including zeros; every route (a) eigenvalue cluster mapped to an allowed level within tolerance; no missing expected cluster, unmatched cluster, omitted zero cell, or spurious sub-ceiling level) |
| G7 | route (a): degeneracy splitting measured at every frozen resolution, with convergence judged only by the frozen statistic; route (b): exact in-band cutoff stability reported per observable, comparing the certification cutoff against twice that cutoff. Route (b) evaluates exact finite group sums, so it has no convergence sequence to report and is never judged by the convergence statistic |
| G8 | reproducibility: repeated runs in the frozen environment reproduce eigenvalues, invariant subspaces, projector ranks, multiplicities, and convergence statistics within their frozen tolerances; bit-for-bit identity is required only for deterministic serialized summaries the implementation explicitly guarantees (the tracked summary artifacts are such summaries, which is what satisfies the platform data policy's regeneration rule); degenerate eigenspaces are compared through projectors or principal angles, never raw eigenvector coordinates |
| G9 | any comparison harness against a transcribed reference carries a transcription mutation (perturb one transcribed cell; the comparison must go red) |

Uncertainty: every route-(a) discretization-limited spectral quantity carries the frozen
convergence statistic; the Weitzenböck gate carries its conditioning-stability check; route
(b) carries the frozen in-band cutoff-stability result; every multiplicity carries its
integer-nearness margin; the trade-off table reports both
routes at matched configurations with the resource measurements that back the cost column.

## 9. Deliverables, provenance, and data policy

- Both prototypes as research scripts (NumPy/SciPy; Taichi only if this graduates to
  production per-frame kernels), plus the summary JSON of every rung, in the repository.
- The measured trade-off table, the route decision with rationale, and the losing prototype
  retained as the M8.4 cross-check tool.
- Heavy arrays are local-only per the data manifest policy (gitignored, never deleted);
  tracked artifacts are the summary JSON/CSV, plots, and the scripts that rebuild
  everything; the manifest table is regenerated by the repository's generator utility, never
  hand-edited. Every run is deterministic from fixed seeds; tracked summary artifacts are
  serialized deterministically (bit-for-bit), and eigenvector-level content reproduces at
  the G8 tolerance level via projectors or principal angles, never raw coordinates.
- Environment records (interpreter, library versions, OS, hardware, seeds) with every
  committed output; SHA-256 hashes of raw outputs.
- The method note follows [`dev_docs/METHOD_NOTE.md`](../../../../../dev_docs/METHOD_NOTE.md)
  in full: equations first, an equation-to-code map with permalinks, results beside their
  gates, at most ~4 inspection artifacts, the recorded adversarial audit before it crosses to
  review, and an explicit list of claims NOT verified. The note distinguishes the § 0
  evidence classes by name and never reports them as equivalent: published-value lookup
  (rungs 1-2, 3a), source-derived theorem evaluation (rung 3b, with its route-(b) overlap
  asymmetry stated), cross-route agreement (rungs 4-6), and post-run author-context reference
  agreement. The not-verified list must include at minimum: the
  truth of the coexact entry rule (per the § 7 matrix state); the M4/M5/M7 family operators
  and every M8.4 claim (§ 0 cap); the M8.3 torsion closed forms; and, unless and until
  M8.5-A's core scalar outcome is `reproduced`, the B certification claim itself (§ 1 gate condition and interim language).

## 10. Pins

| What | Pin |
| --- | --- |
| pre-registration filed against | upstream main `ec877ee0` (2026-07-30; re-pinned at landing if main has moved) |
| M8.2 lock | landed `f18daf27`; close-out `269456b7`; post-lock addendum `0776cc19` (append-only) |
| A protocol | [`m8_5a_reproduction_protocol.md`](m8_5a_reproduction_protocol.md) (landing commit recorded at B's landing) |
| task spec | [`m8_5_task_details.md`](../tasks/m8_5_task_details.md), route table + certification benchmark + DoD |
| author-context reference code (post-run only, § 2) | `scripts/m8_2_first_occurrence.py` at `269456b7`; `scripts/m8_2_indep_reconstruction.py` at `12f4a94a`; `scripts/m8_3_mass_reproducer.py` at its merge commit, recorded at first post-run use |
| procedure record | the PR #350 close-out thread (2026-07-28) |

## 11. Reference pins, verbatim statements, and convention maps

Transferred from the verified source record. Every entry was read from a page actually
fetched; fields that could not be verified are marked UNVERIFIED rather than filled in.
All gates are evaluated at UNIT RADIUS, in the form the sources state them, with the
`R`-dependence carried as the labeled protocol derivation in § 5.

### 11.1 Rung 1: `S³` scalar spectrum and multiplicities

**Polterovich, "Combinatorics of the Heat Trace on Spheres", *Canadian Journal of
Mathematics* 54 (2002), no. 5, 1086-1099, DOI 10.4153/CJM-2002-040-4, § 1.1, p. 1086.**

> Let `S^d` be a sphere with the standard Riemannian metric of curvature +1. The
> Laplace-Beltrami operator `∆` on `S^d` has eigenvalues `λ_{k,d} = k(k + d − 1)`, and each
> `λ_{k,d}` has multiplicity `µ_{k,d}` given by
> `µ_{k,d} = (2k + d − 1)(k + d − 2)! / (k! (d − 1)!)`, `k ≥ 1` and `µ_{0,d} = 1`.

Independent corroboration for `S³` specifically: Lehoucq, Weeks, Uzan, Gausmann and
Luminet, *Class. Quantum Grav.* 19 (2002), no. 18, 4683-4708, DOI
10.1088/0264-9381/19/18/305, § 2.2.

**Convention map.** Source index `k` maps to the protocol's harmonic index `n`; source
dimension label `d` is fixed at 3. Source states a nonnegative Laplacian; the protocol
agrees. Source is unit radius ("curvature +1"); the protocol evaluates at unit radius.
Multiplicities are used as stated. At `d = 3` the statement gives eigenvalue `n(n+2)` with
multiplicity `(n+1)²`.

### 11.2 Rung 2: `S³` one-form exact and coexact towers

**Iwasaki and Katase, "On the Spectra of Laplace Operator on `Λ*(S^n)`", *Proceedings of
the Japan Academy, Series A* 55 (1979), no. 4, 141-145, DOI 10.3792/pjaa.55.141,
Theorem 6, p. 144.** Open access; read from page renders because the OCR layer destroys
sub- and superscripts.

At `n = 3`, `p = 1` the theorem gives the exact tower `λ = ℓ(ℓ+2)` with multiplicity
`(ℓ+1)²`, and the coexact tower `λ = m²` with multiplicity `2(m²−1)`.

`S³`-specific, degeneracy-free alternative: Lauret, "The spectrum on p-forms of a lens
space", *Geom. Dedicata* 197 (2018), 107-122, DOI 10.1007/s10711-018-0322-9, Theorem 2.1,
at `n = 2`, `p = 1`. Sector-labeled alternative: Boucetta, *Publicacions Matemàtiques* 43
(1999), no. 2, 451-483, Théorème 3.1 ii), p. 468, with Prop. 3.15 i) identifying the first
coexact eigenspace with the Killing algebra.

**Convention map.** Source `Δ = dδ + δd`, nonnegative; the protocol agrees. Unit radius.
Sector labels: the source says closed and coclosed where the protocol says exact and
coexact; on `S³` these coincide because `H¹ = 0` (§ 11.4), and the protocol uses
exact/coexact throughout. The coexact tower's starting index differs across sources while
the eigenvalue SET does not; the protocol states its own indexing and does not harmonize
the sources.

### 11.3 Rung 2: Hodge decomposition orthogonality, three named sectors

**Capoferri and Vassiliev, "Beyond the Hodge theorem", *J. London Math. Soc.* 113 (2026),
no. 1, e70431, DOI 10.1112/jlms.70431, eq. (1.1), p. 2.**

> `Ω^k(M) = dΩ^{k−1}(M) ⊕ δΩ^{k+1}(M) ⊕ H^k(M)`, where `dΩ^{k−1}(M)`, `δΩ^{k+1}(M)` and
> `H^k(M)` are the Hilbert subspaces of exact, coexact and harmonic `k`-forms, respectively.

Real-valued forms, orthogonal direct sum, three sectors named exactly as the protocol
names them. Book-level alternative: Warner, *Foundations of Differentiable Manifolds and
Lie Groups*, Springer GTM 94, 1983, DOI 10.1007/978-1-4757-1799-0, Theorem 6.8, p. 223,
with the caveat that its three-line display was recovered only through OCR snippets and
should be re-checked against a physical copy before being quoted.

### 11.4 Rung 2: `H¹(S³) = 0` and the absence of harmonic 1-forms

No single numbered result states this; the verified route is two steps and the protocol
cites it as such rather than implying a single source.

1. **Hatcher, *Algebraic Topology*, Cambridge University Press, 2002, ISBN
   0-521-79540-0, Corollary 2.14, p. 114**: `H̃_n(S^n) ≈ Z` and `H̃_i(S^n) = 0` for
   `i ≠ n`; with the field-coefficient universal-coefficient statement at p. 198.
2. **Warner, Theorem 5.36, p. 206** (de Rham bridge) and **Theorem 6.11, p. 225**
   (harmonic representatives).

Hatcher has no numbered result computing `H^i(S^n)` directly, which is why the citation is
two-step.

### 11.5 Rung 2: Weitzenböck identity on 1-forms, with sign conventions

**Petersen, *Riemannian Geometry*, 3rd edition, Springer, 2016, DOI
10.1007/978-3-319-26654-1, Theorem 9.4.1.**

> The Hodge Laplacian is the Lichnerowicz Laplacian with `c = 1`. Specifically,
> `∆ω = (dδ + δd)(ω) = ∇*∇ω + Ric(ω)`.

Corroborated independently by Petersen and Wink, *SIGMA* 16 (2020), 064, DOI
10.3842/SIGMA.2020.064, Prop. 2.1(a) at `U = 0`; Nicolaescu, *Lectures on the Geometry of
Manifolds*, p. 540; Homma, *Trans. AMS* 358 (2006), no. 1, 87-114, DOI
10.1090/S0002-9947-05-04068-7, p. 109; and Semmelmann and Weingart, *Compositio Math.* 146
(2010), no. 2, 507-540, DOI 10.1112/S0010437X09004333, eq. (3.10).

**Convention map.** All five primary sources agree: `Δ_Hodge = dd* + d*d` nonnegative,
`∇*∇ = −trace(∇²)` nonnegative, curvature term entering with a PLUS. The protocol adopts
that convention.

**UNVERIFIED, deliberately not printed:** the GTM volume number commonly attached to
Petersen could not be confirmed (Crossref volume field null, OpenLibrary series field
null), and the published page of Theorem 9.4.1 is likewise unconfirmed since the
accessible full text is the author's pre-publication manuscript. Only the chapter range,
pp. 333-363, is Crossref-verified. The theorem NUMBER is the citable locator.

**Assembled, not quoted:** `Ric = 2g` on the unit `S³` is derived at protocol level (§ 5)
from constant curvature together with `Ric = (n−1)Kg`, both cited, rather than quoted from
any single page. § 5 labels it as a derivation.

### 11.6 Rungs 3a/3b: opaque case identifier

`M85B-ADJ-01`. The label is opaque: it encodes nothing about the case.

### 11.7 Rungs 3a/3b: packet hashes

Packet I (case-input), SHA-256 `cc8c38f0c6819daa1125b5cfea955bfc6b31d160859e0d226f47aad9399fbecc`.

Packet II (answer), SHA-256 `b89c4a1446fd88df156d2f7505dc3449dcb3f7d4eb9961364d3dc757862ebdf7`.

Listed separately; neither may be inferred from the other. Canonical form for both: keys
sorted, two-space indent, ASCII, LF, single trailing newline.

### 11.8 Rung 3b: the generic one-form multiplicity theorem

**Lauret, "The spectrum on p-forms of a lens space", *Geom. Dedicata* 197 (2018), 107-122,
DOI 10.1007/s10711-018-0322-9, Theorem 2.1**, with Theorem 3.3 supplying the explicit
finite sums.

**Hypotheses, to be checked against the case before evaluation:** the theorem is stated
for `0 ≤ p ≤ n−1` on `Γ\S^{2n−1}` with `Γ` cyclic acting as stated in its § 2, and it
requires the frozen generator convention of that section.

**Use-warning attached to this source:** its § 2 prints the generator with the first
parameter repeated where the last belongs. This was verified in the arXiv version; the
Springer typeset page was NOT independently fetched, so its status there is UNVERIFIED. Use
the corrected form from Lauret, Miatello and Rossetti, *IMRN* 2016, no. 4, 1054-1089, DOI
10.1093/imrn/rnv159, eq. (3.2), which is correct in both.

The case-specific inputs this theorem is evaluated on are sealed in Packet I (§ 4.2); the
theorem itself is generic over a family and carries no case information.

### 11.9 Use-warnings: what must NOT be used, and why

Each source below is excluded from the stated gate because its scope, indexing, bundle
convention or dimensional notation does not match the frozen protocol. **The warning is
limited to that stated use** and is not a general judgement on the source. These exclusions
are part of the contract: each would otherwise produce a plausible and wrong gate.

| Do not use | Reason |
| --- | --- |
| Lauret, Miatello and Rossetti, *IMRN* 2016, Prop. 2.2, for a p-form multiplicity | its printed `p ≥ 1` representation index is shifted relative to the eigenvalue index, and fails the `Γ = 1` sphere check under literal use. Use Lauret 2018 Thm 2.1, or LMR 2015 Thm 1.1 |
| Ikeda and Taniguchi, *Osaka J. Math.* 15 (1978), for a direct `S³` numerical gate | its highest-weight notation reaches a dimension-edge degeneracy in the `S³` specialization: `Λ_j` is defined for `j = 1, ..., m−2`, and `S³` is `m = 2`, so the clause is vacuous and their dominant-weight list contains no plain `Λ₁` there. Not used without a separate convention derivation. Use Lauret Thm 2.1 or Iwasaki-Katase Thm 6 |
| Iwasaki and Katase Theorem 6 for the SCALAR case | Theorem 6 excludes `p = 0`, so it does not support the scalar gate. Use § 11.1 for rung 1 |
| Nash and O'Connor, *J. Math. Phys.* 36 (1995), no. 3, 1462-1505, DOI 10.1063/1.531134, eqs. (3.7), (3.17), (3.19), (3.24), as ordinary multiplicities | these are flat-bundle (TWISTED) sector degeneracies for a nontrivial representation of `π₁`, not an untwisted multiplicity table, and they look identical to ordinary ones. Not used for the untwisted gate without an explicit bundle-character convention map. (Evaluating those equations at `p = 2` returns zero where the untwisted answer is nonzero, which is the protocol's own check that the twist is intended, not a claim made by the source.) |
| Gallot and Meyer as the numerical multiplicity source | Iwasaki and Katase report errors in its multiplicity formulas, so it is not used for the numerical multiplicity gate. Cite it only for the curvature lower bound, via Bär |
| Boucetta's attribution of the `p = 0` multiplicity to Iwasaki-Katase | the cited Iwasaki-Katase theorem excludes `p = 0`, so the attribution cannot support rung 1 |

**Nothing in § 11 may identify the sealed adjudication case.** Slots 11.6 and 11.7 carry
only an opaque label and two hashes; the quotient's identity and action live in Packet I,
and its source, table, row, indexing map and values in Packet II (§ 4.1). A citation naming
the source and row would be an answer-bearing pointer and would unseal the case as surely
as printing the integers, so the pins for rungs 3a and 3b are deliberately not of the same
shape as 11.1 through 11.5. Slot 11.8's theorem is generic over a family and therefore
carries no case information.

## 12. Addenda (post-freeze only)

(none yet)
