# The M8.5-C protocol: target-free qualification of the spectral chassis

> **FILING TEXT, complete (§§ 0 to 16). The redline phase is CLOSED (R1, the seeded R2,
> both redline units' formal sign-off) and § 14 is EXECUTED (eight obligations, two REDs
> raised and resolved at cause). FROZEN at merge; dated addenda only thereafter (§ 16).**
> At filing this document lands as `research/findings/m8_5c_protocol.md`, is FROZEN at merge,
> and admits dated addenda only (§ 16). The freeze boundary, its marker string, the SHA-256 of
> every byte above it, and the exact check command are appended at filing per § 15. The marker
> invariant, everywhere in this protocol: the COMPLETE marker line appears exactly once in the
> file, at the boundary; every check matches the full line only, so inline backticked mentions
> are harmless.
>
> Genesis provenance, none of it operative here: the chassis decision memo (author-side,
> frozen region SHA-256 `44c664d1ceb17949da78e55dcc5fe322cd447375bfe91ab3d33256775f386f4b`,
> three errata below its boundary, all inherited in corrected form; per erratum 1 the
> projector aperture closes by gate 4 alone and no text here cites "gate 11a"), the fourteen
> author rulings D-1 to D-14, and the redline record. Everything operative is IN this document.

## 0. Scope, outcome space, and containment

**What M8.5-C is.** The ONE preregistered qualification attempt for the spectral/Galerkin
chassis: route (b) of the M8.5 charter ("a spectral method in 2I-symmetric harmonics"), built
for the first time as an actual dynamics substrate, TARGET-FREE. The experiment is: determine
whether this chassis is numerically trustworthy, without possessing the target. It is not:
build a chassis capable of answering M8.4. Those are different experiments, and the difference
is the reason this task exists under an engineering charter rather than an M8.4 successor row.

**Outcome space.** Exactly two terminal outcomes, `M8.5-C-QUALIFIED` and `M8.5-C-FAILED`,
defined by the frozen sentences of § 1. ONE ATTEMPT: no repair round after execution begins.
If this protocol is found defective BEFORE execution begins, it is superseded and refiled as a
whole; after execution begins, one attempt means one attempt. Stop conditions and the
disposition of every failure class are § 11's; the containment breach class is § 10's own.

**Cubic only.** The successor law is `cubic_nls` alone, per the memo's decided 11(b).
`saturating` is formally OUT of scope: not deferred, not phase two. The consequence is carried
in § 1's ceiling: any eventual successor experiment is narrower than M8.4 § 9's
two-configuration matrix and every deliverable says so.

**Containment.** A TARGET OBJECT is any of: a nonzero-amplitude nonlinear residual on any
nontrivial `E_ρ`; any M8.4 § 7 scored observable evaluated on any candidate; the M8.4
preregistration document; the McKay distance-to-mass-slot ASSIGNMENT, the mass table, and any
comparison target (mass, `α`, `Λ` values). The build room contains none of these, ever. The
discriminator, stated for the Build Unit reading this section first: the REPRESENTATION-
THEORETIC first-occurrence structure (which harmonic level first carries which sector) is IN
scope, is zero-credit known structure (M8.4 prereg § 8's own category), and is RE-DERIVED
in-room, never loaded; what is excluded is anything assigning those slots to masses or
comparing them to targets. Also in scope: the nine sector bases at ZERO amplitude, whose
operative meaning throughout
this protocol is that NO NONLINEAR TERM IS FORMED on those spaces, the basis fields
themselves being unit-norm objects; everything on `E_R0` including nonlinear runs; and
manufactured operators. Enforcement is executable,
not prose: § 10. The room's enumerated file whitelist is § 12's.

**Standing relations.** M8.4 is CLOSED UNRESOLVED with no target configuration executed
([closeout](m8_4_closeout.md)); this protocol neither reopens it nor inherits its identity.
The successor preregistration `M8.4-R1` is drafted in parallel by the author, held OUT of the
build room, and filed only on the day this protocol adjudicates `M8.5-C-QUALIFIED` with the
frozen cost estimate of § 11; nothing in this room may read it, and nothing in it may read
this room's outputs before adjudication. The grid backend's standing is unchanged: retained
for what `M85B-ADJ-07` certified, excluded from M8.4-lineage dynamics.

**Governance.** The author files and owns this protocol. A fresh-context Build Unit executes
it inside the § 12 room. A context-isolated Adjudication Unit applies the frozen gates to the
pinned output ledger (§ 13) and issues one of the two § 1 sentences; the maintainer reviews
the filed record at merge. Compute ceiling and adjudication mechanics: § 11.

## 1. Claim ceiling and frozen claim language

**The two outcome sentences, frozen byte-for-byte.** The Adjudication Unit issues exactly one
of these, unmodified; every downstream surface quotes, never paraphrases.

```text
M8.5-C-QUALIFIED. The spectral chassis demonstrated its numerical claims on target-free
controls under the frozen gates of this protocol. This is an engineering verdict about an
instrument. It says nothing about MIT dynamics, nothing about M8.4's question, and nothing
about any target observable.
```

```text
M8.5-C-FAILED. The spectral chassis did not demonstrate its numerical claims under the frozen
gates of this protocol. The spectral route closes. The grid backend remains retained for what
M85B-ADJ-07 certified and excluded from M8.4-lineage dynamics, and the M8.4 reopening path
closes with this verdict.
```

**The inherited ceiling, quoted.** The M8.4 preregistration's § 1 governs every artifact this
task produces, quoted verbatim so no M8.5-C surface can drift above the successor's own
ceiling:

> **The limitation that matters most, stated before any run.** OQ1 asks whether a nonlinear
> field equation has solutions "whose energies realize the McKay SLOT STRUCTURE, without
> per-slot tuning." `M4L_Erho` puts one field in each slot by hand: the eight target sectors
> are **installed by construction**, so this family **cannot** answer whether a dynamics
> SELECTS the McKay slots.

> **A second ceiling, on what counts as interesting.** A small-amplitude continuation that
> stays near the free `n = d_ρ` cluster is NOT a striking realization of McKay. Perturbative
> continuation normally preserves ancestry over some interval, so that outcome is expected and
> is reported as **"nonlinear persistence of the installed free structure"**, never as
> "dynamical realization."

The quotation above is reproduced as CEILING LANGUAGE: it carries no target content, and its
presence in this protocol is not an exception to § 0's containment.

**Sentences no M8.5-C artifact may contain, in any wording**: that the chassis is "ready to
answer OQ1" or any equivalent; that any dynamics "selected", "realized", or "produced" the
McKay structure; any numerical value of an M8.4 § 7 scored observable on any candidate; any
comparison against a McKay target. The cubic-only narrowing is stated in every deliverable
that describes scope: this lineage is a cubic-law program, narrower than M8.4 § 9's matrix.
A second inherited limit, stated here so it survives the gap: the chassis's NONLINEAR
identities are tested on `E_R0` and the manufactured fibre-valued extension only; cubic
equivariance on the nontrivial sectors rests on the § 2 derivation and the linear ports,
and on NO nonlinear check anywhere, because § 0 forbids the act. `M8.4-R1` will run
nonlinear dynamics on precisely the sectors this instrument could not test nonlinearly,
and every qualification claim carries that limit.

## 2. The scientific object, the pinned law, and the discrete family

**The object.** The continuum `M4L_Erho` equation is the scientific object. The Galerkin
systems below are its CANONICAL APPROXIMANTS: for cutoff `N` and bundle `ρ`, the truncation

    N_N^ρ(ψ_ρ) = P_N^ρ [ f(⟨ψ_ρ, ψ_ρ⟩) ψ_ρ ]

with `P_N^ρ` the orthogonal Galerkin projector onto the retained section space (complete
levels `n ≤ N`), is uniquely specified with no convention. The cutoff ladder tests CONVERGENCE
toward the continuum; `N` is never a definition and never a model choice.

**The pinned law.** Read from `m4_ewt/wave_engine.py` at repo commit `c9dc3796`, line 29:
"Time-integrated leapfrog vector wave equation: `∂²ψ/∂t² = c²∇²ψ` (the nonlinear restoring
term `−dV(ψ)` is added in P2)". SECOND order in time, leapfrog family. The cubic law is
`v_mode 1`: `V = (c1/4)u²` on `u = ⟨ψ,ψ⟩`, so `dV = c1⟨ψ,ψ⟩ψ`. Standing states `ψ = e^{iωt}φ`
satisfy the finite residual

    R(φ; ω) = c²Δφ − c1⟨φ,φ⟩φ + ω²φ,    ω² → c² d_ρ(d_ρ+2)/R²  at zero amplitude,

in the frozen unit convention `c = R = 1`, so `λ = n(n+2)` exactly and every frozen number in
this protocol reads in those units. Any other unit choice rescales the § 6 ladder and is
prohibited. The coupling is pinned for ALL M8.5-C qualification work: `c1 = +1`, the
restoring sign the § 7 elliptic control requires; `c1 = 0` appears only as the § 7 linear
manufactured operator; `M8.4-R1` is not bound by this pin.

**Structural facts the gates cite.** (i) The fluctuation operator of `R` is REAL-linear on
`R^{2m}`: a complex `m × m` Jacobian drops the `⟨δφ,φ⟩φ` term entirely (relative difference
1.000 against the true real Jacobian; design input `jacobian_check.py`). Every Jacobian in
this protocol is assembled and scored as a real operator on `R^{2m}`. (ii) `iφ` is an exact
kernel vector of that Jacobian at every solution, the `U(1)` gauge mode, and it sits INSIDE
the scored eigenspace. (iii) The residual symmetry of every sector's problem is

    G_ρ = (U(1) × SU(2)_right) / K_ρ,    dim G_ρ = 4,   K_ρ ≅ Z/2,

carried as a DERIVED result: right multiplication commutes with the left `2I` action, so it
descends to the quotient and preserves every `E_ρ`; right translations are isometries, so
they commute with `Δ`; the cubic map is EQUIVARIANT, `N(R_gψ) = R_g N(ψ)`; Peter-Weyl and
Frobenius give level `n` the shape `Hom_2I(V_n, W_ρ) ⊗ V_n^*` with the right action on
`V_n^*`; and multiplicity one AT THE SCORED LEVEL makes the scored space `≅ V_l^*`
irreducible. The multiplicity premise carries TWO named sources: the McKay first-occurrence
tables for the eight sectors at `l = d_ρ` (re-derived, not transcribed, by the § 14 suite),
and the `2I` invariant ring for `R_0`, Molien `(1+t³⁰)/((1−t¹²)(1−t²⁰))`, invariant levels
through 36 at {0, 12, 20, 24, 30, 32, 36} each of multiplicity one, so `H_{R0,12} ≅ V_12^*`,
complex dimension 13, spin 6. Full chain with proofs: the symmetry derivation note filed
beside this protocol. Because complete levels are retained, `P_N^ρ` commutes with the right
action and the SEMI-DISCRETE system retains all of `G_ρ`; its exact invariants are energy,
`U(1)` charge, and the three right-`SU(2)` momenta (gate 5).

**The `K`/`J` witnesses, defined here so this document stands alone.** For an assembled
operator block `A` written in the orthonormal § 3 basis: `K(A) = ‖A − Aᴴ‖₂`, the
self-adjointness defect, and `J(A) = max |Im λ(A)|`, the non-real-spectrum witness. They are
DIFFERENT predicates: a non-Hermitian block can have an entirely real spectrum, so `K` can
fire where `J` does not, and an implementation that conflates them is broken. Gate 1's
machine allowances on the linear diagonal operator, both RELATIVE and of one form:
`K ≤ 100 ε_mach ‖A‖₂` and `J ≤ 100 ε_mach ‖A‖₂`. The relative form for `J` is deliberate:
a general eigensolver returns imaginary parts near `ε_mach ‖A‖₂` on a real block, which at
the top rung (`‖A‖₂ ≈ λ_60 = 3720`) is about `8.3e-13`, so an absolute `1e-12` allowance
would leave a margin of 1.2; the relative allowance is `8.3e-11` there, still nine orders
below the `1.8e-01` that failed `R_2` in P1A, so discrimination is untouched.
The S1b record pinned in § 15 is the instrument's provenance, not its definition; the
definition is this paragraph.

**Implementation status, two layers, never conflated.** The § 14 suite is a PRE-FREEZE
bridge for the REPO SCALAR PRIMITIVES only; that is the entirety of what its green certifies.
The `W_ρ`-valued bases and the production Galerkin system built in-room earn their own green
under the § 3.4 symmetry-realization checks (gates 3 and 5), which reuse the suite's
identities and mutations as NEW in-room implementations with NO inherited credit.

**Fallback, operative as written; TWO named triggers, both pre-freeze only.**
TRIGGER S: the § 14 suite (script 5) not fully green against the scalar primitives at
freeze. TRIGGER L: the § 14 lattice tables (item 4) not armed at freeze while the suite is
green. Under TRIGGER S the FULL degraded contract binds, carried identically by §§ 5 and 8:
(a) branch enumeration runs as the deterministic multi-seed search with no symmetry
classification, branch labels recorded UNAVAILABLE; (b) the zero-count gate runs in
MEASURED-NULLSPACE mode: the inertia count against the same threshold `τ = 1e-8 · ‖J‖₂`,
splitting and
leakage scored on the measured complement, the predicted-versus-measured comparison
SUSPENDED, so a degraded run never dispatches an over-count to the fold detector on symmetry
grounds; (c) gate 5's conserved-set check runs on energy and charge alone. Under TRIGGER L,
two things change and nothing else: enumeration degrades per (a), the multi-seed search with
labels UNAVAILABLE; and the zero-count gate drops its lattice-class cross-check while
KEEPING the predicted-versus-measured form, since `rank_R Z(φ)` needs no lattice. Gate 5 is
untouched.
AFTER execution begins, nothing degrades: a red in a § 3.4 symmetry-realization check is a
gate failure and STOP-QUAL, attributed per § 3.4's port-versus-bases split, never a
contract downgrade.

## 3. Sector bases, fibres, and metrics (the load-bearing build)

Every existing asset is `R_0` machinery; the `W_ρ`-valued section bases for the eight
nontrivial sectors are THE build item this task exists to produce and qualify.

**3.1 Fibre realizations.** For each nontrivial irrep `ρ` of `2I`, an explicit UNITARY matrix
realization on `W_ρ`, built in-room from a verified generating pair of the explicit
120-icosian set. Certification, each item mutation-armed under gate 3: the homomorphism table
exact on all `120 × 120` products; unitarity to rounding on all 120 elements; the trace
character equal to the per-element character table at every element (the table itself armed
by § 14's P1); `ρ(−1) = ±I` per Schur, matching the sector's spinorial parity.

**3.2 Intertwiners and sections.** For each retained level `n` and sector `ρ`, the
intertwiner space `Hom_2I(V_n, W_ρ)` is computed as the nullspace of the stacked equivariance
constraints `ρ(γ) A = A π_n(γ)` over the generating pair; its dimension must equal the
multiplicity from the character route, computed in-room, with the two routes agreeing exactly
(gate 3). Sections are `f_{A_i,j}(x) = A_i π_n(x) e_j`: `i` indexes an orthonormal basis of the
intertwiner space, whose dimension is the multiplicity, and `j` runs over the `n + 1`
multiplet columns. Multiplicity EXCEEDS ONE at levels the ladders retain: first at `n = 15`
for `R8`, 16 for `R7`, 18 for `R5`, 21 for `R6`, 22 for `R3`, 26 for `R4`, 31 for `R1`, 37
for `R2`, and 60 for `R0`; five sectors carry a multiplicity-2 level at the agreement
ladder's bottom rung, and `R0` meets one at Control B's top rung (the pinned 677 modes at
`N = 60` already encode it). The multiplicity-space basis is fixed by the equivariance
constraints only up to unitary mixing, so it is pinned OPERATIONALLY: the deterministic § 3
construction (stacked-constraint SVD in a frozen element order, Lowdin against the level
Gram, sign fixed by the first nonzero component) produces ONE basis object per sector and
rung, hashed into the input manifest; every consumer in this protocol, both § 4.3 routes
included, reads THAT object, so multiplicity-space mixing is common everywhere and cancels
in every comparison. The retained space at cutoff `N` is complete levels only.

**3.3 Grams and orthonormalization.** Per-level Grams are computed by the § 4 production
rule (integrand degree `≤ 2N`, inside its exactness). The basis is orthonormalized per level
by LOWDIN symmetric orthonormalization. Rationale, stated so no reader reasons from a wrong
premise: invariance gives `A(g)† G A(g) = G`, which does NOT make `G` commute with the
represented action; what Lowdin provides is that `G^{1/2} A(g) G^{-1/2}` is unitary, so the
orthonormalized basis carries a UNITARY realization of the same action, and symmetric
orthonormalization is the canonical ordering-independent choice where QR would freeze a
column order. After Lowdin: represented action unitary to rounding (checked), pointwise
equivariance preserved (checked), Casimir labels preserved (gate 2 applies to the
orthonormalized basis).

**3.4 Certification, gate 3 in full.** (a) Reynolds dimension against a method-disjoint
character route, both built in-room, agreeing at every level ANY § 4 rule reads at that
rung, which is `n ≤ 3N`: the § 4.2 monitor projects onto the band `N < n ≤ 3N`, so
certification does not stop at the retained levels; (b) pointwise equivariance of sampled
sections at random points and random group elements, same `n ≤ 3N` range; (c) the
coefficient-to-node map verified round-trip to rounding, same range; (d) the
symmetry-realization checks against the BUILT bases, with the re-run-versus-port split
stated: P0, P1, C5 and C6 of the § 14 suite are group and character theory and re-run
UNCHANGED; C1 through C4 and C7 are PORTED, new implementations reusing the same identities
and mutations with no inherited credit, with their arenas SPLIT by § 0's containment: the
LINEAR ports C1, C2, C4 and C7 (right-action realization, projector commutation,
complete-level retention, spectrum-versus-action) run against the `W_ρ`-valued bases and the
production projector, at zero amplitude, linear operations on fields being outside § 0's
prohibition; the NONLINEAR port C3 (cubic equivariance) runs ONLY on `E_R0` and the
manufactured `E_R0 ⊗ C²` extension, because forming the cubic of a nonzero field on a
nontrivial `E_ρ` is the § 0 forbidden act. Gate 3 qualifies the nontrivial bases and their
right action entirely at zero amplitude, in § 0's operative sense: no nonlinear term is
formed on those spaces, the basis fields themselves being unit-norm; no nonzero-amplitude
nonlinear nontrivial-sector state exists anywhere inside M8.5-C. Port discipline: each
ported check first reproduces
its green parent and
its red mutation on the scalar primitives in-room; a failure THERE is an instrument-port
defect, recorded as such, STOP-QUAL attributed to the instrument, never to the bases and
never a § 2 degradation; port green followed by bases red is a gate 3 failure attributed to
the bases; (e) the in-room first-occurrence RE-DERIVATION over a fixed scan range, labelled per the
column's standing rule as INDEPENDENT-METHOD, NOT BLIND: the expected distances sit in the
room's own scripts and in this very sentence, so nothing here is a discovery; the
character-route scan is method-independent of those transcriptions and must return first
`n` with `mult(ρ, V_n) > 0` equal to `d_ρ` with multiplicity exactly one, all eight
sectors.
The Adjudication Unit, not the Build Unit, compares the in-room tables against the shipped
M8.2 tables (the locked M8.2 reference is the adjudicator's authorized third input, § 11),
preserving the M8.2 lock § 3 independence discipline;
the room never loads the shipped tables (§ 12).

## 4. The projector: production rule, monitor rule, dual implementation

**4.1 The production quadrature, gate 4.** Hopf-coordinate product rule exact to polynomial
degree `4N` on `S³`: `x = (cos η cos ξ₁, cos η sin ξ₁, sin η cos ξ₂, sin η sin ξ₂)`, uniform
`4N + 1` points in each of `ξ₁, ξ₂`, Gauss-Legendre `2N + 1` points in `u = sin²η`, weights
normalized to unit total measure; node count `(2N + 1)(4N + 1)²`. For band-limited `ψ` at
cutoff `N`, every Galerkin coefficient of `⟨ψ,ψ⟩ψ` has degree `≤ 4N`, so the rule returns
`P_N^ρ[⟨ψ,ψ⟩ψ]` to rounding with ZERO aliasing. Synthesis and analysis are SEPARABLE Hopf
transforms, FFTs in `(ξ₁, ξ₂)` with per-mode one-dimensional transforms in `u`: at the top
rungs the node counts (3,613,153 at `N = 48`; the monitor's 23,588,101 at `N = 60`) put a
dense node-by-mode array out of reach, so the fast transform is load-bearing for the § 11
ceiling, not an optimization, and § 14's pre-freeze obligations include a MEASURED wall-clock
at the top rung of each ladder so the ceiling is a live constraint when it freezes.
Gate 4's arms: the node-drop mutation (a rule
exact only to `2N` must err at O(1)), and the dual-implementation agreement of § 4.3. The
design-input record, `1.5e-01` against `2.2e-15`, was taken at `N = 3`, which sits on neither
production ladder: it is layer-1 evidence about the primitives and nothing more. The
production implementation earns gate 4 in-room at every rung of both § 9 ladders.

**4.2 The cascade monitor, gate 6, with its OWN rule.** Aliasing is retired by exactness;
cascade is physics: the cubic generates the band `N < n ≤ 3N` BEFORE projection, and the
monitor measures that ENTIRE band. Its integrands have degree `≤ 6N`, so the monitor carries
its own `6N`-exact rule with its own node-drop mutation (design-input record: `6.0e-15`
against `2.35e-01` under the production rule, whose use here would alias exactly the content
the monitor exists to see). The frozen reading is

    C_N = ‖P_{N<n≤3N} f(ψ)‖₂ / ‖P_{n≤3N} f(ψ)‖₂,

bounded in [0, 1], denominator the TOTAL generated band; the band projection uses the same
§ 3 construction, certified over `n ≤ 3N` per § 3.4. Threshold: `C_N > 0.1` at a rung
flags the run CASCADE-LIMITED in the ledger; all rungs run regardless (§ 9), so the flag
stands in the record and triggers nothing. The flag is a
REPORTED diagnostic: it can neither veto a § 9 CONV label (that would be post-hoc `N`
selection) nor rescue a non-converged one. Arm allocation: the monitor runs on the
production rung spaces and on Control B, which carries the live nonzero reading once
`η > 0`; Control A's manufactured arena (§ 6) sits outside the monitor entirely; the
monitor's mutation is injected high-band content on the production rung spaces at `2×`
threshold; the ledger records which arm carried which.

**4.3 The dual implementation, the same gate.** A second, transform-free route computes the
identical projected cubic by Clebsch-Gordan contraction in coefficient space. The two routes
must agree to rounding on the PREREGISTERED manufactured field set, frozen here: FORTY
pseudorandom band-limited fields at each rung of both § 9 ladders, twenty in the `E_R0`
retained space at that rung and then twenty in the manufactured multi-component extension
(`E_R0 ⊗ C²`) that exercises the fibre-valued code path (forty per rung, twenty per
category; neither count is a rung value); no agreement field lives on a nontrivial `E_ρ`,
per § 10. Generation is SINGULAR, one stream
from `numpy.random.Generator(numpy.random.PCG64(20260901))`: traverse the AGREEMENT ladder
ascending {24, 32, 40, 48}, then the Control-B ladder ascending {36, 44, 52, 60};
within each rung, the twenty scalar fields indexed 0 through 19, then the twenty `⊗ C²`
fields indexed 0 through 19; within a field, one coefficient per basis element in level-major, then
intertwiner-index, then multiplet-index, then component-index order (the `C²` component
innermost, varying fastest; scalar fields have no component index); per coefficient, real
part then imaginary part, each a standard normal draw; the field then normalized to unit L2
norm. Serialization is
canonical: coefficients as little-endian IEEE-754 float64 `(real, imag)` pairs in C-order,
concatenated in the draw order above; the SHA-256 digest of that byte stream is written to
the input manifest BEFORE any gate executes, and the digest, not the recipe, is the
preregistered object of record. Disagreement at any field of any rung is STOP-QUAL. Neither route
materializes the rank-4 coupling tensor (at `M ≈ 1000` retained modes it exceeds `10¹²`
entries); coupling is applied via transforms or via per-level CG contraction only.

**4.4 Exclusion.** The 7200-node cloud is EXCLUDED from the production projector; its moment
qualification is linear-purpose and finite-degree. It may serve as a diagnostic sampling
cross-check only, and no gate reads it.

## 5. Continuation, the primary route

**The solve.** Standing states are zeros of `R(φ; ω)` with `(φ, ω)` explicit unknowns and the
amplitude constraint `‖φ‖₂ = a` as its own row. Newton runs on a SQUARE bordered system with
MATCHED border rows and columns, one pair per predicted continuous zero of the branch under
continuation (the border columns are the measured orbit-tangent vectors, the rows their
adjoints); amplitude `a` is the continuation parameter. Pseudo-arclength is retained SOLELY
as the frozen fallback and activates only on the § 5 fold trigger below, never by judgment;
its algorithm is fully frozen here: Keller pseudo-arclength on the FULL continuation vector
`y = (φ, ω, a)`, all components in the pinned `c = R = 1` units and unweighted, with a
SECANT predictor through the last two accepted points in `y` (on first activation, the
tangent from the bordered system); arclength step equal to the `y`-space distance between
those two points, halved on corrector failure, at most 4 halvings; corrector constraint the
standard arclength hyperplane, orthogonality of the `y`-update to the secant direction,
alongside the unchanged symmetry borders; at most 20 arclength steps per candidate event;
return to the ladder by one
fixed-amplitude solve at the next prescribed rung as soon as the branch's amplitude passes
it. TWO distinct terminations, so a real fold in the amplitude direction is never mislabeled:
if the detector has tripped and `a` then DECREASES monotonically over 5 consecutive arclength
steps, the branch is recorded AMPLITUDE-LIMITED with its fold amplitude (the maximum `a`
reached), a legitimate closure, and its ladder rungs above the fold amplitude are recorded
UNREACHABLE, not failed. Exhausting the 20 steps with neither re-entry nor the
amplitude-limited pattern is the reportable outcome ARCLENGTH-FAIL for that branch, never
silence.
Stopping, frozen, on the FULL AUGMENTED residual, since the bordered system, not `R` alone,
is what Newton solves: the scaled augmented residual is the concatenation of the field rows
`R / (λ_l ‖φ‖₂)` (`λ_l` the control level's free eigenvalue, 168 for Control B), the
amplitude row `(‖φ‖₂ − a) / a`, and each border row `c_k / ‖φ‖₂`, defined:
`c_k(φ) = ⟨ẑ_k, φ − φ_prev⟩`, with `ẑ_k` the unit-normalized orbit-tangent vector measured
at the previous accepted point (at a branch's first point, at the seed), the same vectors
that form the border columns; converged when the ∞-norm of that scaled vector is `≤ 1e-12`.
Control A is
UNBORDERED and retains the `R`-only criterion `‖R‖₂ ≤ 1e-12 · λ_scale,A · ‖φ‖₂` (§ 6). Hard
cap 30 iterations; non-convergence or divergence is the reportable outcome NEWTON-FAIL for
that branch and rung, never silence and never a retry with altered settings.

**The zero-count gate (gate 7's core), non-degraded form.** At each converged state, on the
UNBORDERED real Jacobian, which gate 5 certifies SYMMETRIC:

- MEASURED count: by SYLVESTER INERTIA, `n_neg(J − τI) − n_neg(J + τI)` from two `LDLᵀ`
  factorizations at the shifts `±τ`, `τ = 1e-8 · ‖J‖₂`; `‖J‖₂` is `max|λ|` from the dense
  symmetric eigendecomposition below (the sizes allow it everywhere this gate runs: at most
  1354 real dimensions). The inertia count is AUTHORITATIVE.
- MEASURED-ZERO SUBSPACE, a second frozen operation because inertia yields a count and not
  vectors: the eigenvectors with `|λ| ≤ τ` from one dense symmetric eigendecomposition
  (LAPACK `syevd`, deterministic). Consistency, required: the extracted basis holds EXACTLY
  the inertia count of vectors, each with residual `‖Jv − λv‖₂ ≤ 1e-10 ‖J‖₂`, pairwise
  orthonormal to `1e-12`; any mismatch is an instrument failure, STOP-QUAL attributed to the
  instrument. Splitting and leakage are scored on this basis's orthogonal complement.
- PREDICTED count: `rank_R Z(φ)` with `Z(φ) = span_R{iφ, T_1φ, T_2φ, T_3φ}`, computed by SVD
  of the four-column real matrix at the `1e-8` relative threshold, on the CONVERGED state,
  never the seed; cross-checked against the branch's pinned lattice class. The achievable
  values are 1 at `l = 0`, always 3 at `l = 1`, and 3 or 4 at `l ≥ 2`.
- The gate BINDS AT REGULAR POINTS: measured = predicted passes; measured < predicted is a
  gate FAILURE (the gauge-breaking mutation arms this side); measured > predicted dispatches
  to the FOLD-OR-BIFURCATION-CANDIDATE outcome below, never to STOP-QUAL. PERSISTENT means
  TWO consecutive regular-flagged steps: an over-count on both, with no detector trip at
  either, is the gate failure; a single-point over-count is recorded and continued past.
  Rationale for two, frozen with the rule: the amplitude ladder has four rungs, so a
  three-step definition could never fire on a systematic over-count entering at the tail;
  two consecutive is the minimum separating persistence from a single-point anomaly, and
  the recorded single-point case stays visible in the ledger.

Splitting and leakage are scored on the orthogonal complement of the measured-zero subspace.
Scored quantities per branch, rung and amplitude: cluster position, cluster splitting,
principal angles, leakage, zero count, branch label.

**The observable definitions, frozen, because gate 9 scores nothing it cannot compute.**
Work on the UNBORDERED real Jacobian restricted to `Z(φ)^⊥`, the orthogonal complement of
the measured-zero subspace; let `d_c = D_l − z`, with `z` the measured zero count and `D_l`
the REAL dimension of the scored free block: `2(l+1)` for `H_{R0,12}` (26), and `2(l+1)²`
for a full-S³ scalar level `l` (18 for Control A's level 2).

(a) CLUSTER IDENTIFICATION. At the FIRST amplitude point of a branch, the cluster is the
`d_c` eigenpairs smallest in `|λ|`, and the pick is VALIDATED by a recorded separation
check rather than an appeal to the free limit, which is never evaluated: require
`|λ|_(d_c+1) − |λ|_(d_c) ≥ 0.25 · gap`, with `gap` the control's OWN inter-level distance
(168 for Control B, § 6; `λ₆ − λ₂ = 40` for Control A), so the thresholds are 42 and 10;
the expected bottom-rung value is near `0.9 · gap` for Control B; a violation is recorded
IDENTIFICATION-AMBIGUOUS, the point is unscored, the smallest-`|λ|` pick stands, and the
flag rides the branch label. At each SUBSEQUENT point the identification is
continuity-tracked by a PINNED PROXY, never re-picked by magnitude and never a
combinatorial search: rank the current eigenvectors by `‖P_prev v_i‖₂²`, with `P_prev` the
orthogonal projector onto the previous point's cluster span; take the top `d_c`; break
ties by ascending eigenvalue.

(b) CLUSTER POSITION = the arithmetic mean of the cluster's eigenvalues. (c) CLUSTER
SPLITTING = max minus min of the cluster's eigenvalues.

(d) PRINCIPAL ANGLES. The free reference is `F_⊥`, the orthogonal PROJECTION of the free
level-`l` eigenspace onto `Z(φ)^⊥`, then orthonormalized, its rank RECORDED at every point
(projection, not intersection: the two readings coincide at the ladder's bottom, where
`Z(φ)` sits inside the level-`l` block and the rank degenerates to `2(l+1) − z`, and
diverge above it, where `φ` carries content at levels 20 and beyond; the projection reading
keeps the full free reference and is the pinned one). The angles are the canonical angles
between the cluster's span and `F_⊥`. (e) LEAKAGE = the sine of the largest canonical
angle in (d).

(f) ZERO-COUNT TRANSITIONS, indexed precisely, because scored quantities carry branch,
rung and amplitude while gate 9's `e_r` is a difference ACROSS RUNGS at fixed amplitude.
A `z` that differs between rungs `r` and `r+1` at the same amplitude invalidates that rung
difference: record ZERO-COUNT-TRANSITION with both counts and decline that `e_r`. A `z`
that changes ALONG AMPLITUDE within a rung is a branch event: record it; the rung
differences are untouched. Consequence, pinned here rather than decided inside the run:
four rungs give three differences and two ratios, so declining one difference leaves a
single ratio, and an observable left with a single ratio yields NO CONV label, which is
§ 9's own rule that a non-converged observable yields no label. The transition is a
reportable event, never a failure.

Control A's exact-arithmetic reference values (§ 6, § 14 item 5) are computed under THESE
definitions, so a mismatch there is an implementation error, never a definitional dispute.

**The fold detector, shared, same floor as the count.** One frozen detector serves the
over-count disposition AND the pseudo-arclength trigger: a sign change of the bordered
determinant between ORIENTATION-CONSISTENT consecutive continuation steps, or the bordered
Jacobian's smallest singular value falling below `1e-8` times its largest, the SAME relative
floor as the zero count, so no dead zone exists between "numerically null" and "not yet a
candidate". Orientation is pinned: border tangent vectors are sign-aligned to the previous
step by positive inner product before the determinant is formed, so an SVD or QR sign flip
can never manufacture a spurious sign change. On trigger: the point is recorded
FOLD-OR-BIFURCATION-CANDIDATE with its data, pseudo-arclength continues the branch past it,
and the zero-count comparison is not scored at the flagged point. A detector trip with
measured = predicted at the adjacent regular points is recorded and does not fail anything;
the persistent (two-step) over-count rule above is the only failure path.

**Degraded modes, § 2's TWO triggers, carried here identically.** Under TRIGGER S (the
suite): branch enumeration runs as the deterministic multi-seed search below with no
symmetry classification and branch labels recorded UNAVAILABLE; the zero-count gate runs in
MEASURED-NULLSPACE mode, the inertia count against the same threshold `τ = 1e-8 · ‖J‖₂`,
splitting and
leakage scored on the measured complement, the predicted-versus-measured comparison
SUSPENDED and no over-count ever dispatched to the fold detector on symmetry grounds. Under
TRIGGER L (the lattice, suite green): enumeration alone degrades to the same multi-seed
search with labels UNAVAILABLE; the zero-count gate keeps predicted-versus-measured and
drops only the lattice-class cross-check. The fold detector stays live under both. Nothing
that happens after execution begins activates this paragraph (§ 2).

**Branch enumeration (gate 7 with § 6's Control B).** The frozen procedure: the isotropy
lattice of `G_ρ` acting on `V_l^*`, computed and pinned per § 14 for `l = 0..7` and `l = 12`
with one CANONICAL seed representative per class, fixed in the § 3 basis. For each class in
the pinned lattice: Newton restricted to the class's fixed-point subspace, seeded at the
canonical representative scaled to the bottom amplitude rung, then continued up the § 6
ladder. Equivariant-Branching-Lemma existence is claimed ONLY where the fixed-point subspace
is COMPLEX one-dimensional (one phase orbit, the phase-quotiented variant); every other
class gets the deterministic search with no completeness claim. ALL branches found are
reported, none selected; "no branch found in class X" is a reportable outcome; the
deterministic multi-seed fallback (degraded mode) is itself deterministic AND common across
cutoffs, so § 9 always compares like with like: the `N = 36` rung's twenty SCALAR fields
from the § 4.3 stream are PROJECTED onto `H_{R0,12}`, any with projected norm below `1e-2`
of its pre-projection norm is discarded, the survivors are normalized with order preserved,
and that ONE list, pinned by the manifest, is the seed set reused identically at every
cutoff rung.

## 6. The two continuation controls (gate 7)

**Control A, manufactured, the known-answer arm.** Arena: the TWO-LEVEL truncated Galerkin
system on full-S³ scalar levels {2, 6}, 58 complex modes, a manufactured operator under
gate 10's manufactured category; the manufactured projector `P_{2,6}` truncates the cubic's
out-of-band content BY DEFINITION of the system, so Control A is self-contained, has no
cutoff ladder, and runs ONCE. Pinned construction: path `φ*(s) = s v₂ + (s²/2) v₆` with
path parameter `s ∈ {0.1, 0.3, 0.5}` (`s` is a path label, NOT a § 5 continuation
amplitude; Control A uses none of § 5's border or amplitude machinery), `v₂, v₆` the first
canonical basis vectors of their levels; `ω*² = λ₂ + 1 = 9`, and `9 = n(n+2)` has no
integer solution, so the linear part is nonsingular on EVERY scalar level, not only the two
the path occupies. The forcing `G(s) := −R_{2,6}(φ*(s); ω*)` is added to the manufactured
residual, making `φ*(s)` an EXACT zero at every `s`; the forcing breaks the `U(1)` and
right symmetries, so the reference zero count is 0, Newton runs SQUARE and UNBORDERED at
fixed `ω*`, and a quiescent fold detector is part of the reference. Newton stopping scale,
pinned: `λ_scale,A = λ₂ = 8`. Seed, literal: `φ_seed = φ*(s) + 1e-3 · ‖φ*(s)‖₂ · v_{2,2}`,
with `v_{2,2}` the SECOND canonical basis vector of level 2; Newton must converge back to
`φ*(s)` under the § 5 stopping rule with `λ_l = λ_scale,A`. Reference values of cluster
position, splitting, leakage and zero count at each pinned `s` come from the
implementation-disjoint route: exact-arithmetic assembly of the 116-dimensional real
Jacobian and 50-digit eigenquantities via a disjoint library, pinned in the input manifest
BEFORE any float64 run; splitting and leakage are NONZERO by construction at every `s`
(the G-NULL-c requirement: never an exact-zero or no-change case). Control A arms the
measurement machinery; it cannot arm the symmetry prediction, which is Control B's. Its
correctness rides entirely on the exact-arithmetic reference for its own system; its fields
join no other gate's input set (levels 2 and 6 carry no 2I-invariant content, so they do not
even embed in `E_R0`), and the § 4.2 monitor does not run on its arena.

**Control B, real, the symmetry arm.** Continuation from `H_{R0,12} ≅ V_12^*` (spin 6,
complex dimension 13), permitted because M8.4 prereg § 8 makes `E_R0` nonlinear runs
qualification. Amplitude ladder in the dimensionless coupling `η = |c1| a² / gap`,
`gap = 168` (the λ-distance to the nearest other `R_0` invariant level, § 2's units), closed
rungs `η ∈ {0.05, 0.1, 0.2, 0.4}`; cutoff rungs {36, 44, 52, 60} (§ 9). Control B exercises
the full § 5 machinery where its lattice is richest: finite isotropy classes exist at spin
6, branch enumeration runs the pinned `l = 12` lattice, the zero-count gate runs its
predicted-versus-measured inertia form, and the fold detector is live. Control B is
MACHINERY qualification only: its gates are mechanical consistency (zero counts match at
regular points; enumeration returns the lattice's classes or reportable no-branch outcomes;
§ 9 convergence behaves), and no Control B result carries physics content of any kind.

**Cascade allocation, restated from § 4.2.** The live monitor reading is Control B's once
`η > 0`; the injected-content mutation runs on the production rung spaces; Control A's
arena sits outside the monitor entirely.

## 7. The supportive time arm (gate 8)

Integrator: Stormer-Verlet leapfrog, the pinned law's own family, second order in time.
Both controls evolve on EXACTLY INVARIANT small subspaces and are integrated in their
confined coordinates; one bounded full-width wiring check connects them to the production
system. This keeps the supportive arm far from the § 11 ceiling without weakening it.

**Control (i), linear standing wave.** The manufactured LINEAR operator, `c1 = 0` (gate
10's manufactured category), on `H_{R0,12}`: `ψ(t) = e^{iωt}φ` with `ω = √168` and `φ`
PINNED with every convention explicit, so the pin is checkable by anyone: the canonical
multiplet basis `v_1 .. v_13` carries weights DESCENDING `m = +6 .. −6`;
`φ ∝ Σ_{j=1}^{13} (1 + j/13) e^{i j/3} v_j`, normalized to unit L2. The complex weights are
load-bearing, per the § 14 record: any REAL vector has `M₁ ≡ 0` identically (the generator
is real antisymmetric) and any weight-symmetric modulus profile has `M₃ = 0`, which is
exactly how the original equal-weight pin failed its own assertion. Generators, pinned:
`J₊ v(m) = √((6−m)(6+m+1)) v(m+1)`, `J₋ = J₊ᵀ`, `J₃ = diag(m)`; `T₁ = (J₊ − J₋)/2`,
`T₂ = i(J₊ + J₋)/2`, `T₃ = i J₃`; momenta `M_a(ψ, ψ̇) = Re⟨ψ̇, T_a ψ⟩`. The setup
assertion requires all three `|M_a| ≥ 1e-6 · S_cons` AND two analytic identities of this
specific pin, armed: `|M₁/M₂| = tan(1/3)` to `1e-12` (the pin is a z-rotation by 1/3 of a
real-weighted vector, so the transverse momenta sit at exactly that ratio), and `M₃` equal
to its closed form `+ω Σ_j m_j |1 + j/13|² / Σ_j |1 + j/13|²` (the sign is the content:
`M₃ = Re⟨iωφ, iJ₃φ⟩ = +ω φᴴJ₃φ`, and the descending-weight, ascending-modulus profile
makes the sum negative). A failure of any of these is
an instrument defect, not a pass. The subspace is exactly invariant because the operator is
linear and level-preserving; the confined system is the 13-complex-dimension block. The
nonlinear law does NOT confine a level-12 eigenfunction, which is why `c1 = 0` is pinned
here.

**Control (ii), nonlinear constant section.** `c1 = +1` (§ 2's pin), real initial data
with ZERO initial velocity, `q(0) = q₀ = 1`, `q̇(0) = 0`, on `E_R0`'s constant mode; the
cubic of a constant is a constant, so the one-complex-dimension subspace is exactly
invariant. Exact solution `q(t) = q₀ cn(ωt, k)`, `ω² = c1 q₀² = 1`, `k² = 1/2`; period,
literal: `T_cn = 4K(1/√2) / (√c1 · q₀) = 4K(1/√2) ≈ 7.4163`.

**Steps, ladders, horizon, frozen per control.** `h = T_ctrl / 100` with `T_ctrl` the
control's OWN period (`2π/√168` for (i); `T_cn` for (ii)); closed `dt` ladder
`{h, h/2, h/4, h/8}`; horizon 10 `T_ctrl`. Reference integrator, method-disjoint: adaptive
RK45 at relative tolerance `1e-10`, absolute `1e-13`, on the controls only.

**Full-width wiring check, separate and bounded.** 1,000 steps of the full-width system at
`N = 36`, each control under ITS OWN pinned operator (`c1 = 0` for control (i), a
manufactured operator under gate 10 and the registry's single named law exception, recorded
in the ledger; `c1 = +1` for control (ii)), never under a law that would evict the control
from its subspace; `h_wire = T_min/100`, `T_min = 2π/√λ_36`,
initialized on each control's subspace embedding; the full-width state must match the
confined evolution at every sampled step to `1e-10` relative in L2. This checks the wiring; the
confined runs carry the physics of the arm.

**Error metric, frozen.** The elliptic trajectory crosses zero, so no pointwise relative
error exists; the metric is `err = max_t ‖ψ_h(t) − ψ_ref(t)‖₂ / max_t ‖ψ_ref(t)‖₂`.
Bounds: `err ≤ 1e-3` at the finest rung, contraction `≥ 3×` per rung (second order
predicts about 4×).

**Conservation, split by kind, with a nonzero normalization.** The frozen scale is
`S_cons = max_t (‖ψ_ref‖₂ ‖ψ̇_ref‖₂ + ‖ψ_ref‖₂²)`. Charge and the three right-`SU(2)`
momenta are discrete momentum maps of LINEAR symmetries of the discrete Lagrangian,
conserved by the variational integrator at rounding: the bound is
`|I(t) − I(0)| ≤ 1e-12 · S_cons`, ABSOLUTE and normalized, because control (ii)'s charge
and momenta are identically zero and a relative drift is undefined there. Energy carries
the oscillatory leapfrog envelope, `≤ 5 (ω_ctrl h_r)² / 8` relative to its own mean, plus
the secular test: the linear-fit slope times the horizon must not exceed 10 percent of the
envelope. Allocation, so identical zeros are never read as passes: control (i) carries the LIVE
sharp checks (charge and all three momenta, nonzero by the `φ` pin); control (ii) carries
the live energy check, and its charge and momenta, identically zero by construction (real
data, right-invariant constant mode), are RECORDED as identically zero, never counted as
passes. Conservation is supportive and never sufficient; this arm touches nothing in § 5's
scored quantities.

**Gate 8's mutation, named.** Kick-drift symplectic Euler (first order) on control (ii)
must FAIL the `3×` contraction gate; its first-order error cannot keep up, and an
implementation where it does has a broken error metric or a broken integrator.

## 8. The gate table

| # | Gate | Check | Mutation arm | Green parent | Spec |
| --- | --- | --- | --- | --- | --- |
| 1 | G-LIN wiring | exact-zero `K`/`J` on the linear diagonal operator, per § 2's definitions and machine allowances | injected anti-Hermitian coupling between two retained modes must go red | the diagonal operator itself | § 2 (definitions); S1b provenance via § 15 |
| 2 | G-LABEL | Casimir / round Laplacian through the sampling map on EVERY retained basis element, eigenvalue AND sector, over `n ≤ 3N` | one deliberately mislabeled element must go red | certified § 3 bases | § 3 |
| 3 | sector bases | § 3.4 (a) to (e), including the ported symmetry-realization checks and the port discipline | each § 3.4 item's own arm; ports re-arm on scalar primitives first | § 14 suite greens (layer 1 only) | § 3 |
| 4 | projector exactness | production `4N` rule to rounding; § 4.3 dual-route agreement on the pinned field set | node-drop to `2N` must err O(1); route disagreement is STOP-QUAL | § 4.3 field packets, hash-pinned | § 4 |
| 5 | structural identities | BEFORE any integrator exists, and CONTAINED per § 0: the nonlinear identities run on `E_R0` and the manufactured `E_R0 ⊗ C²` extension ONLY, never on a nontrivial `E_ρ`, because no nonzero-amplitude nonlinear nontrivial-sector state exists anywhere inside M8.5-C. Semi-discrete equivariance under `G_ρ` (ported C3, on those spaces, at production rungs); Jacobian symmetry `‖J − Jᵀ‖ ≤ 1e-12 ‖J‖` (the residual is a gradient); Noether identities on frozen PHASE-SPACE states, `(ψ, ψ̇)` = the rung's § 4.3 scalar fields paired (0,3), (1,4), (2,5), since the pinned law is second order and `F` lives on `(ψ, ψ̇)`: `abs(⟨∇I, F⟩) / (‖∇I‖₂ ‖F‖₂) ≤ 1e-12` for each `I` in {E, Q, M₁, M₂, M₃}, no time integration involved (drift belongs to § 7 alone); zero-denominator convention, frozen WITHOUT any cross-dimensional state scale: the § 14 preflight runs the pairing on the three pinned pairs and confirms every denominator `D_I = ‖∇I‖₂ ‖F‖₂` nonzero (§ 4.3's unit normalization is what guarantees the state norms themselves); a pair with `D_I` exactly zero in float64, unreachable by construction and recorded for completeness, advances deterministically to the next pinned pair with the substitution in the ledger. DEGRADED per § 2 TRIGGER S: the conserved-set check runs on energy and charge alone, the equivariance and Jacobian-symmetry items unchanged. | non-gradient perturbation; symmetry-breaking coupling; intra-level spectrum break (ported C7) | § 14 suite greens (layer 1); § 3 bases | § 2, § 3.4(d) |
| 6 | cascade monitor | `C_N` with its OWN `6N` rule, band certified per § 3.4 | injected high-band content at `2×` threshold; node-drop on the `6N` rule | § 4.2 construction | § 4.2 |
| 7 | continuation controls | Control A against symbolic references; Control B mechanical consistency; zero-count gate; fold detector | gauge-breaking term (measured `1.3e-01` design input); seed-pinned prediction must red | Control A's exact `φ*(s)`, with the perturbed seed RETURNING to it as the parent check | §§ 5, 6 |
| 8 | time arm | both known trajectories within bounds; conservation split sharp/supportive; the full-width wiring check | kick-drift symplectic Euler on control (ii) must fail the `3×` contraction gate; broken-symmetry drift | RK45 reference on the same controls | § 7 |
| 9 | convergence | contraction machinery binds on the Control-B ladder (the converging observables live there); the agreement ladder provides rung coverage for gates 2, 4 and 6, which are pass-fail at rounding on every rung | injected non-contracting error above the floor must red | § 9 scales | § 9 |
| 10 | executable partition | qualification mode cannot load a target configuration; coverage check | deliberately injected violation must be caught | clean-room run record | § 10 |
| 11 | law scope (RECORD, not a check) | recorded DECIDED: option (b), cubic only; `saturating` out of scope | not applicable, a decision record | not applicable | § 0 |

Every arm names its green parent per the standing rule; an arm whose parent is not green is
vacuous and does not count as run. Gates 1 through 10 must be able to go red; row 11 is a
RECORD, not a check, kept in this table so the numbering matches the memo's eleven-gate
floor, and it is exempt by construction: nothing to implement, nothing to fail.

## 9. Ladders, convergence, and CONV (gate 9)

**The closed rung sets.** AGREEMENT ladder `N ∈ {24, 32, 40, 48}`; Control-B ladder
`N ∈ {36, 44, 52, 60}`; time-arm `dt` ladder `{h, h/2, h/4, h/8}`. Four rungs each. The
sets are CLOSED: no rung is ever added, and the ladder is never extended, whatever any
reading says. Binding, stated plainly: the CONTRACTION machinery below binds on the
Control-B ladder alone, where this section's scored observables live; the `dt` ladder is
CLOSED here but its contraction and final-error bounds live entirely in § 7 under gate 8,
with § 7's own `1e-3` threshold, and none of this section's scales apply to it; the
agreement ladder exists for rung coverage of gates 2, 4 and 6, whose checks are pass-fail
at rounding at every rung and contract nothing.

**Scale-normalized errors, with a LOCAL noise floor.** Per scored quantity `q`, the rung
error is `e_r = |q_{r+1} − q_r| / S_q` with frozen nonzero scales: cluster position
`S = λ_l` (168 for Control B); splitting `S = 168` (the § 6 gap); leakage `S = 1` (already
dimensionless); principal angles ABSOLUTE radians, `S = 1`. The floor is `100 ε_mach` in
normalized units, `ε_mach = 2.22e-16` the binary64 MACHINE EPSILON, so the floor is
`2.22e-14`. The contraction rule is LOCAL, per successive ratio, and all rungs always run:
if the preceding error is above the floor, the next must satisfy `e_{r+1} ≤ e_r / 3`; if
the preceding error is at or below the floor, that single ratio is unscored; any later
error that rises above the floor re-enters the ordinary rule; and the final requirement,
`e_last ≤ 1e-6` or `e_last` at or below the floor, ALWAYS binds, so an early floor visit
can never forgive a later excursion. A measured splitting below `1e-6 · S` is reported
BELOW-RESOLUTION, never as a converged value. Mutation: an injected non-contracting error
above the floor must go red.

**CONV, verbatim discipline.** A CONVERGED observable yields its verdict whatever it says. A
NON-converged observable yields NO label. The two outcomes are never merged. The
cascade-limited flag is a reported diagnostic and can neither veto a CONV label nor rescue a
non-converged one. "Stable under the ladder" as prose is not a gate.

## 10. The executable partition (gate 10)

**The line, drawn here and not deferred.** Forming a nonlinear residual on any nontrivial
`E_ρ` at nonzero amplitude IS touching the target. Everything at zero amplitude (§ 0's
operative sense), everything on `E_R0`, and everything on manufactured operators is
qualification. Enforcement is by CODE, not prose; the q3a lesson stands: the contract author
owns any ambiguity, so this section enumerates rather than gestures.

**The arena registry.** The qualification executable carries an ENUMERATED registry of
arenas, each with an identifier and a constructor pinned in the input manifest:

| id | arena | nonlinear-permitted |
| --- | --- | --- |
| `A-R0-N{24,32,40,48,36,44,52,60}` | `E_R0` retained space at the named rung | YES |
| `A-R0C2-N{...}` | the `E_R0 ⊗ C²` manufactured extension at the named rung | YES |
| `A-CTRLA` | Control A's two-level full-S³ arena (§ 6) | YES |
| `A-CTRLI` | control (i)'s `c1 = 0` linear operator on `H_{R0,12}` (§ 7) | YES, trivially: its law has no nonlinear term |
| `A-SECTOR-{R1..R8}-N{...}` | nontrivial `W_ρ`-valued bases and projectors at the named rung | NO: linear operations only |
| `A-MUT-*` | manufactured mutation operators declared per gate | per their declaring gate |

The registry keys on ARENA alone; the law is global `c1 = +1` (§ 2) with ONE named
exception, § 7's control-(i) full-width wiring run under `c1 = 0` on `A-R0-N36`, recorded
in the ledger as the registry's law exception. The NONLINEAR evaluation entry point accepts
a registry identifier and REFUSES any whose column reads NO; there is no code path by
which a nontrivial-sector basis object reaches the
nonlinear evaluator, enforced structurally (the nontrivial-basis module is not in the
nonlinear evaluator's import closure, which § 12's gates verify) and dynamically (the
refusal above, logged).

**The coverage check, two-sided and mutation-armed.** The run FAILS if: (i) any registry
arena was never exercised by its owning gates, so every manufactured category is REACHABLE
and reached, Control A's arena, control (i)'s operator and the `⊗ C²` extension included; or
(ii) the run log contains a nonlinear evaluation tagged with any arena outside the
nonlinear-permitted set. Arms: a deliberately injected out-of-registry nonlinear call must
be caught by (ii); a deliberately skipped registry arena must trip (i).

**Dispositions, all three cases pinned (gate 10 sits outside STOP-QUAL's enumeration).**
(1) A COVERAGE failure, any registry arena never exercised: a qualification failure;
terminate, adjudicate `M8.5-C-FAILED`. (2) A REFUSED prohibited call arising outside an
armed mutation context: the partition HELD, but a production code path attempted a
forbidden arena, which is an instrument defect; terminate, `M8.5-C-FAILED`, attribution
INSTRUMENT in the ledger. (3) An actual forbidden nonlinear evaluation EXECUTED, the
refusal missed, is the CONTAINMENT BREACH: immediate termination; the BREACH record in the
ledger as its own class; the attempt still adjudicates `M8.5-C-FAILED`, preserving § 0's
exactly-two terminal outcomes; and the consequence for `M8.4-R1`'s provenance pre-written
here: a breached room can no longer certify "a chassis that was never in a room with the
target," so the breach record and that voiding statement enter the filed record, and no
successor may inherit the never-in-a-room clause from this room's artifacts.

## 11. Stop conditions, adjudication, and the ceiling

**Execution begins**, defined so "before execution" has a boundary that a preflight cannot
trip: at the first GATE record the commissioned Build Unit writes to the output ledger.
Preflight results (import closure, read closure, room setup) land in the COMMISSIONING
record, never the ledger. A protocol defect found pre-execution supersedes the protocol as
a whole (§ 0); after that first GATE record, one attempt means one attempt.

**STOP-QUAL.** Any failure among gates 1 to 9, or the resource ceiling, terminates the
build before any further work; no discretionary repair round exists. The § 4.3
dual-implementation agreement, the cascade arm, and the manufactured control are
GATE-INTERNAL to gates 4, 6 and 7, named in the memo's own enumeration for emphasis; they
define no additional stop class beyond their gates. Gate 10 carries its own § 10
disposition instead. Instrument-attributed failures (§ 3.4's port discipline, § 5's
basis-consistency mismatch) are STOP-QUAL with the attribution in the ledger.

**TARGET-LOCK and COMMON-LAW**, recorded verbatim for inheritance since no target run
exists inside M8.5-C: once a first target-bearing residual is ever formed (in `M8.4-R1`,
never here), no cutoff, `dt`, integrator, tolerance, normalization, amplitude ladder,
coupling, stopping rule, projector implementation, cascade threshold, branch procedure or
sector treatment may change, and any post-target change to a frozen element voids every
target run made before it, the record keeping both. If any sector requires sector-specific
numerical or physical retuning, the common-law claim FAILS; reporting the sector result is
fine; rescuing the cross-sector claim is not.

**The ceiling.** The single attempt's end-to-end compute is capped at 48 wall-clock hours
on the build machine, all rungs and arms inclusive; crossing the cap is STOP-QUAL. The § 14
pre-freeze wall-clock rehearsal is what makes this a live constraint rather than a hope.

**Adjudication.** A context-isolated Adjudication Unit receives EXACTLY three objects: this
protocol, the pinned output ledger, and the LOCKED M8.2 reference tables, authorized for
the § 3.4 comparison alone; nothing from the room. It applies the frozen gates mechanically
and issues exactly one of § 1's two sentences plus the per-gate verdict table. The
maintainer reviews the filed record at merge. Both the adjudicator design and
the ceiling are inherited verbatim per the standing rulings.

**The reopening criterion, pre-committed.** `M8.5-C-QUALIFIED` plus a frozen cost estimate
under the ceiling files `M8.4-R1` as a fresh preregistration the day of the PASS
adjudication; anything else leaves the chassis as platform infrastructure and files
nothing. `R_0` results may inform the successor's design (M8.4 prereg § 8's pilot
category); nothing from a nontrivial `E_ρ` can, because the room never forms one. The R1
filing must carry two statements verbatim: the § 1 inherited limit, that this instrument
tested no nonlinear nontrivial-sector state and R1 runs nonlinear dynamics on precisely
the sectors the chassis could not test nonlinearly; and, if § 10's breach class ever
fired, the voiding of the never-in-a-room-with-the-target clause.

## 12. The room: commissioning mechanics and the whitelist

**The whitelist, enumerated positively; anything unlisted is out.** This protocol; the
symmetry derivation note filed beside it; the five design-input scripts and their README;
`route_a_nonabelian.py`; the § 4.3 packet generator spec; the isotropy lattice tables and
Control A reference values (§ 14 products, pinned pre-freeze); and the code the Build Unit
writes in-room. Named OUT, so the enumeration cannot be misread as an oversight:
`m8_4_preregistration.md`, `MODELS.md`, `__M8_model_briefing.md`, `m8_2_first_occurrence.py`
and its shipped tables, the M8.5-A artifacts, and `m8_4_closeout.md`. The room builds its
own first-occurrence structure; the ADJUDICATION-side comparison against shipped M8.2
tables (§ 3.4) is the only place the two meet, preserving the M8.2 lock § 3 independence.

**Enforcement, both mutation-armed.** Import closure: every module in the room imports in a
per-module subprocess with a sentinel, so a `sys.exit` in an import cannot false-green.
Runtime FILE-READ closure: the executable cannot read outside the enumerated manifest.
Origin checks probe declared entry points by name, so an absent module cannot silently drop
out of the scan.

**Manifests.** The INPUT manifest (whitelist, arena constructors, § 4.3 digest, § 3
basis-object hashes, Control A references, lattice tables) and the OUTPUT ledger are
separate hashed objects; the input manifest is closed and hashed at commissioning, before
any gate runs. Freeze markers on every frozen in-room document obey one invariant: the
COMPLETE marker
line appears exactly once, at the boundary; every hash check matches the full marker line
only, so inline mentions are harmless.

**Commissioning.** A fresh-context Build Unit, commissioned per the S1b pattern; the
commission is non-governing provenance and is archived after the run. The Adjudication
Unit is commissioned separately and receives only the three § 11 objects.

## 13. Deliverables, provenance, and the output ledger

**The ledger.** Append-only JSON records, one per event, hashed at close. REQUIRED record
types, each with its named fields: GATE (gate id, arena id, rung, parent status, mutation
status, measured values, and a law-exception flag where § 10's single named exception
applies); ARM-ALLOCATION (which control or space carried which live check,
and which identical zeros were RECORDED rather than passed: § 4.2's monitor allocation,
§ 7's conservation allocation); SUBSTITUTION (gate 5's pair advances, § 5's seed
discards); DISPOSITION (NEWTON-FAIL, ARCLENGTH-FAIL, AMPLITUDE-LIMITED with its fold
amplitude and its UNREACHABLE rungs as a named field, FOLD-OR-BIFURCATION-CANDIDATE points
with their data, ZERO-COUNT-TRANSITION with both counts, IDENTIFICATION-AMBIGUOUS with the
separation values, BELOW-RESOLUTION splittings, CASCADE-LIMITED flags); ATTRIBUTION
(instrument-port versus
bases, § 3.4; instrument versus run, § 5; instrument, § 10 case 2); COVERAGE (the § 10
two-sided result);
RESOURCE (cumulative wall-clock against the ceiling, per gate); BREACH (§ 10's class, if
ever). An adjudication cannot proceed on a ledger missing a required type for work that
ran; that is itself an instrument failure.

**What ships.** The Build Unit ships the ledger, the logs, the input manifest, and every
hash. The Adjudication Unit ships the § 1 sentence and the per-gate verdict table, nothing
else. At merge the records land under `research/m8_5c/` beside the protocol's findings
entry, following the column's house layout.

## 14. Pre-freeze obligations

Named work that must land BEFORE this protocol freezes. BLOCKING items have no fallback;
FALLBACK-COVERED items degrade per § 2 if not green.

| # | obligation | status class |
| --- | --- | --- |
| 1 | `right_translation_check.py` green at its pinned hash, re-run record | FALLBACK-COVERED (§ 2 TRIGGER S, the full degraded contract) |
| 2 | re-run records of the other four design-input scripts with the dependency pin | BLOCKING |
| 3 | measured wall-clock at the TOP RUNG of each ladder for the § 4.3 agreement gate, the § 4.2 monitor, and Control B, by author-side SCRATCH code, recorded and discarded, non-inheritable (Control A has no rung cost); this licenses FEASIBILITY of the ceiling only, never a prediction of the Build Unit's in-room cost, and crossing the ceiling is STOP-QUAL regardless of what the rehearsal showed | BLOCKING |
| 4 | the isotropy lattice tables for `l = 0..7` and `l = 12` with one canonical seed representative per class, computed by an ARMED derivation script (design-input script 6, AUTHOR-SIDE: the script never enters the room; its product tables and their hashes do, per § 12's whitelist) | FALLBACK-COVERED (§ 2 TRIGGER L, enumeration-only) |
| 5 | Control A's exact-arithmetic reference values at the three pinned `s`, by the implementation-disjoint route, ARMED, pinned | BLOCKING |
| 6 | control (i)'s nonzero-momenta setup assertion RUN on the pinned 13-vector `φ` | BLOCKING |
| 7 | gate 5's generic-pair denominators `D_I` confirmed nonzero on the three pinned pairs at EVERY rung of both § 9 ladders | BLOCKING |
| 8 | the fallback cross-check: § 2 (both triggers), § 5, gate 5's row, and § 14's own table carry the SAME two degraded contracts, verified by reading all four | BLOCKING |

A BLOCKING item not done means the protocol does not freeze; there is no partially frozen
state. Scratch code from item 3 never enters the room.

## 15. Pins

Two pin classes, kept distinct. FREEZE-TIME pins are values written into this document and
covered by its hash. COMMISSIONING-TIME pins are values fixed at room commissioning and
recorded in the input manifest, which is hashed before any gate runs; they are LISTED here
by name so their absence is checkable, but their values live in the manifest.

**Freeze-time:**

| what | pin |
| --- | --- |
| chassis decision memo | frozen region SHA-256 `44c664d1ceb17949da78e55dcc5fe322cd447375bfe91ab3d33256775f386f4b`, checked with the memo's own boundary command; three errata below its boundary inherited in corrected form; "gate 11a" cited nowhere |
| pinned law | `m4_ewt/wave_engine.py` @ repo commit `c9dc3796`, line 29 quoted in § 2 |
| `route_a_nonabelian.py` | file SHA-256 `e55728534c6b1d611e364551297df7711b95cd87f73d3f97da753a743419931b`; contained in upstream main at the filing pin above |
| `mode_count.py` | `794d7063bea01ad3637a96c33e4e01e7fd71acfc636eb2d443c4f3cecb6959a6` |
| `exact_quad_check.py` | `51fa9b52bc9d3367fbb20baf084235da260c7f15b36fad1b78b31770e7f55ed7` |
| `cascade_quad_check.py` | `7f31bc305150c9d5e37ec312fe8fd00060583a87427feb7a56d6743d2fd66bfa` |
| `jacobian_check.py` | `dce687ea793522be40c4d06d25df19a0dfb99c72309a19361e220077007eb9d8` |
| `right_translation_check.py` | `752baa2de36dbaea9fe4108ec6df9121351b386453d01b06bd8111562de8f547` |
| § 5's observable-definitions block (its own sub-region digest, the object the § 14 records cite: NON-CIRCULAR by construction, since no § 15 edit can move it, unlike a whole-file or frozen-region hash which this table pins the records against) | `823e90662450a0e83297b81d0b217d27c37a7dd8957e8607dcd2c8e316c60544` |
| the symmetry derivation note | filed beside this protocol at filing, named in § 2; SHA-256 `5b231af423066aa9b9902a4ca9b67f19a54158f6ab36339ae711992839820d1a`, recomputed at freeze if the note changes pre-freeze |
| S1b decision rule (`K`/`J` instrument, G-DISCRIM) | frozen region `c44c603a…`, per its own boundary record |
| M8.4 closeout, M8.9 records, M8.5-B preregistration (`M85B-ADJ-07` scope) | the merged findings files, referenced by name; outside the room |
| the rulings record | `M8_5C_OPEN_DECISIONS.md`, author-side provenance, non-operative |
| upstream main at filing | `9ba2a6646630015d41598e1206285748b740aed8` |
| the isotropy lattice tables (§ 14 item 4, with the disposition split MAXIMAL / ABSORBED / GENERIC-STRATUM and the scope note on the generic stratum's reachability) | `efe6c2d64ba9664edaf8b49b9026a683bf93aa78c0c994f09f41bd5a63a7b504` |
| Control A reference values (§ 14 item 5, per-quantity routes recorded) | `dde909f6746994468805c55c1eae192b449b2eb4862e2432f61d55a9f3be4e44` |
| the § 7 amendment's verification record (§ 14 item 6 v2, both analytic arms) | `6615b3966b063ee8c1498f050da559749d105b1cd116c51567f1886dbbcf28a2` |
| design-input script 6, `s14_item4_lattice.py` (AUTHOR-SIDE; the script never enters the room, this hash is provenance) | `f2218ac3a9b2aa4c288463f9ed680aeae5bb457968538974e055fc8e27f0df6e` |
| the § 14 package (provenance, dispositions, supplement, raw records) | ships with the filing under `research/m8_5c/s14/`; `raw/package_hashes.txt` inside it hashes every member |

**Commissioning-time (values in the input manifest):** the § 4.3 packet digest; the § 3
basis-object hash per sector and rung; the arena-constructor registry of § 10; the room
whitelist as enumerated; the Build Unit commission record.

## 16. Addenda (post-freeze only)

Empty at freeze. **Scope, per the platform's § 12.1 standard for freezes from 2026-08-14
on, so a typo never costs a dated addendum:** the append-only, never-in-place rule binds
the COMMITMENTS, which are every gate, threshold, ladder, rung set, definition, frozen
literal, outcome sentence, disposition and stop condition in §§ 0 to 15. Surrounding PROSE
may be corrected in place, and any in-place correction REPUBLISHES the frozen-region digest
in the freeze record below and in every surface that quotes it; the correction and the
superseded digest are recorded in the same commit. A change to a commitment is not an
addendum and not a prose fix but a supersession, and after execution begins (§ 11's first
GATE ledger record) there is no supersession and no in-place change of any kind. The freeze boundary inserted at filing is the full line
`<!-- M85C-FREEZE-BOUNDARY -->`; that string appears in the body only inside this naming
sentence, which is safe because the check command matches the complete marker line alone,
and the frozen region is every byte above the marker, verified by the S1b-pattern `sed`
command recorded beside the digest at filing.
<!-- M85C-FREEZE-BOUNDARY -->

**Freeze record.** The frozen region is every byte ABOVE the boundary marker line. SHA-256:
`e253558b5a767084d4d7777550ac72de5b8a0591ec3d2b847108f04e17c0cc6b`

```bash
sed '/^<!-- M85C-FREEZE-BOUNDARY -->$/,$d' m8_5c_protocol.md | shasum -a 256
```

Filed 2026-08-28; frozen at merge; § 16 governs everything below this line thereafter.

## Addendum 1 (2026-08-29): attempt A1 terminated without adjudication; post-execution supersession granted once

The single commissioned attempt, A1, stopped mid-gate-4 and neither § 1 sentence issues.
The archived record is [`../m8_5c/a1/`](../m8_5c/a1/); the governance thread is
[#501](https://github.com/openwave-labs/openwave/issues/501). This addendum records the
maintainer's ruling and changes no frozen text: the region above the boundary is
byte-identical, digest `e253558b5a767084d4d7777550ac72de5b8a0591ec3d2b847108f04e17c0cc6b`.

**The defect.** Gate 4's frozen mutation arm, "node-drop to `2N` must err O(1)", states a
mathematically false proposition on every arena the gate runs on. Level-`n` content has
Hopf angular parity `n mod 2`; every sector is parity-pure (`ρ(−1) = ±I` is central, § 3.1),
and the `2N` rule's aliasing lattice, multiples of the odd `K = 2N + 1`, sends every alias
into the OPPOSITE parity, which a sector projector never computes; the cheapest same-parity
alias sits at `4N + 2`, outside the degree-`4N` band. Measured in-room: `7.44e-15` where the
clause demands O(1). Independently reproduced twice with separate code (author; maintainer,
including the cross-parity leakage at O(1) where the true value is zero). The design-input
record that certified the arm (`1.5e-01` at `N = 3`) drew mixed-parity fields, levels 0
through 3, a configuration that cannot arise in any legal arena of this protocol. Gate 6's
node-drop arm dies by the same argument, analytically as of this addendum. The unit
recorded a STOP naming the dead arm and no gate-4 verdict of either color; its stop was the
unit's judgment under the commission's non-governing discipline, not a protocol clause.

**The ruling (maintainer, #501, quoted verbatim from its disposition and framing):**

> Read mechanically, § 11's STOP-QUAL on a red arm, or § 13's missing-record clause, both
> terminate in `M8.5-C-FAILED`. The reason that sentence does not issue is not that the
> text is silent; it is that it would state a falsehood about the instrument. The `4N` rule
> was exact to `7.7e-15` at `N = 24`; the chassis did what gate 4 claimed, and the arm that
> was to prove the gate could go red was proven incapable of going red on any legal arena.
> This ruling is therefore a post-execution supersession under a new identity, granted
> once, on a demonstrated and independently reproduced mathematical falsity in a frozen
> operative clause. It is named as such so it cannot be cited as "one attempt was waived".

Ratified on that ground: A1 is TERMINATED and never resumes; no adjudication unit is
commissioned and neither § 1 sentence issues; A1 carries no execution credit anywhere, and
a successor starts at gate 1 in a new room under a NEW protocol identity. The A1 record
stands permanently as the run that discovered the protocol defect, not as evidence for or
against the spectral chassis. Ratification was conditioned on the A1 bytes being filed in
the repository so the two archive hashes are checkable claims, which the filing carrying
this addendum does.

**What the successor must change, recorded here so the requirement survives this document:**
both parity-dead node-drop arms replaced by mutation arms demonstrated live on the actual
gate arenas; every mutation arm, and any design-input record standing in for a gate's
evidence, demonstrated to fire on a field drawn from the arena its gate runs on BEFORE
freeze, with the arena named per arm in the gate table; an administrative PROTOCOL-INVALID
disposition class outside the adjudication unit, requiring a demonstrated internal
contradiction reproduced by the maintainer with independent code before ratification,
terminating the attempt and permitting only a new protocol identity; and the two A1
ledger-schema defects fixed (`cumulative_seconds` truly cumulative; `gates_completed`
distinguishing parent-executed from arm-set-complete). The exact successor arms belong in
the successor's own frozen text, not here.
