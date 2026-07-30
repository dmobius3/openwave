# M8.2: Field-Dynamics Program, Core Contract and Family Modules

> **Status: LOCKED 2026-07-27** ([PR #350](https://github.com/openwave-labs/openwave/pull/350)).
> The core contract (§§ 1-5) and the family modules (§ 6) are frozen as of this date; later
> changes go in a dated addendum at the end of this file, never in-place. The modular
> ARCHITECTURE is complete, and the index-compatibility question is governed by the platform
> standing rule
> [`dev_docs/CROSS_MODEL_TESTING.md`](../../../../../dev_docs/CROSS_MODEL_TESTING.md)
> (2026-07-24), which this document's audit produced (discussion #312). Under its UNIFORM
> DEFAULT no family blocks on an author reply, so nothing here waited on an inbox.
>
> **What the lock does NOT cover.** The native quotient operators for M4 and M5 (§§ 6.1-6.2),
> the native-M7 Hodge decomposition (§ 6.3), and the Zenodo-deposit byte-check of the
> mass-spectrum source (§ 1) are execution-appendix content under § 7: core-internal, not
> author-gated, and gating the first target-bearing run of each family (M8.4) rather than this
> lock. § 8 tracks them per module.
>
> | Component | State |
> | --- | --- |
> | Modular architecture (the split itself, §§ 2-8) | COMPLETE |
> | Core certification logic (§§ 2-5) | conformed to the platform standard this pass |
> | M4 family module (§ 6.1) | AUTHOR-DECLARED geometric displacement (@jeffsyee, 2026-07-24): native/untwisted, T-c N/A SETTLED; internal-triplet table now known inapplicable to native M4 |
> | M5 family module (§ 6.2) | AUTHOR-CLARIFIED (@JarekDuda, 2026-07-25) spacetime/frame tensor: native/untwisted, T-c N/A; "M5 + P" remains an M8-owned object (his latitude is level-choice, not a soldering endorsement); native operator selection remains core-internal |
> | Native M7 module (§ 6.3) | DEFAULT native/untwisted, T-c N/A; column PARKED (author away), so held under the default |
> | `M7_ad` adjoint extension (§ 6.4) | a separate object owned by THIS column; HELD pending @marcf999 / @pwerbos |
>
> **The systematic finding (now the platform standing rule).** M5's `M` (a Lorentz-covariant
> tensor) and M7's `A, J` (spatial one-forms) carry geometric or spacetime indices, not manifest
> internal ones; M4's `ψ` is seeded as a radial spatial displacement and the pinned code does not
> establish its three components as either. Per `CROSS_MODEL_TESTING.md`, no family is treated as
> an internal adjoint multiplet without an AUTHOR declaration (component count and an agent's
> reading of code do NOT count) or a pre-registered soldering/extension. Under the default every
> family is native/untwisted, T-c "not applicable" (a neutral 🚧-weight status, never a failure);
> any soldered result grades "family + prescription P", a different object, not the family; and
> author silence is a valid terminal state, so the default stands and M8.2 locks on schedule.
> Owner: Blake Shatto. Reviewer: maintainer.

## 1. Immutable records (pins)

Certification targets and code are pinned as immutably as the claims they anchor. The family CODE
is pinned at `c9dc3796ba7812a9ddd607647de945abff806057` (abbrev `c9dc3796`): the version actually
read, characterized, and redlined, i.e. the FROZEN family-code version under test, independent of
the branch's rebase onto current upstream. A targeted pin audit (2026-07-27) against current
upstream confirmed the pinned paths: the M4 engine, M5 engine files, M7 functional, and the M8.1
note are byte-identical; the M5 canonical has only documentation additions since (a notation
clarification consistent with our handling, plus one new physics-finding row), the model object
unchanged, so `c9dc3796` is retained as the frozen version.

| Pin | Value |
| --- | --- |
| family-code SHA (verified-at, frozen version under test) | `c9dc3796ba7812a9ddd607647de945abff806057` |
| M4 code of record | `openwave/xperiments/m4_ewt/wave_engine.py` @ `c9dc3796` |
| M5 code of record | `openwave/xperiments/m5_liquid_crystal/research/m5_theory_canonical.md` (spec) + `openwave/xperiments/m5_liquid_crystal/engine2_pde.py` + `openwave/xperiments/m5_liquid_crystal/medium.py` @ `c9dc3796` |
| M7 reference functional | `openwave/xperiments/m7_hydroboros/research/scripts/m7_functional.py` @ `c9dc3796` |
| MIT `2/R²` edge benchmark | in-repo `openwave/xperiments/m8_mit/research/findings/m8_1_method_note.md` @ `c9dc3796` (the 2/R² edge level, independently certified by the M8.1 blind eigensolve; this note is the reproducibility anchor for the 2/R² claim) |
| platform governing contract | `dev_docs/CROSS_MODEL_TESTING.md` @ `1ec2cd97` (added 2026-07-24; its own immutable provenance layer, later than the family-code pin) |

Any later upstream change to a pinned file requires a re-signed module version.

**MIT source papers (content-hashed, not just located).** Record numbers do not freeze content,
so each source is pinned by the SHA-256 of the exact audited PDF. The coexact/galois SSRN PDFs are
byte-identical to their venue-submission copies (cross-checked). The `2/R²` claim anchors to the
SHA-pinned M8.1 note (independently certified), so no separate paper hash is needed for it.

| Source | Identifier | Audited PDF (date) | SHA-256 | size (B) |
| --- | --- | --- | --- | --- |
| mass-spectrum ("The Spectrum") | Zenodo `10.5281/zenodo.18603975` (version DOI) | `SPECTRUMv1.pdf` (2026-02-25) | `31529023af05de3bc79689d1976d7709514c4e222ee572512deb4b7145ca74e5` | 330453 |
| coexact-gap | SSRN `6968698` | `coexact_gap.pdf` (2026-07-16) | `5f4e30691c696eca9efe018a21dcd79bb151a3ba38996c1a2362beeec87b2699` | 504607 |
| galois-pair | SSRN `7129118` | `galois_pair.pdf` (2026-07-16) | `87ea8049bd8eadcc329c6ddbc193f88bfbb44043b9365ea1c21ce4e8fb9c41f8` | 558085 |

Mass-spectrum DOI confirmed via the author's publication registry (`The Spectrum` = Zenodo
`18603975`, 2026 Feb); the hash is of `SPECTRUMv1.pdf`, the audited local copy (2026-02, consistent
with the registry). All three source hashes are final; a byte-check of the mass-spectrum copy
against the immutable Zenodo deposit is available on request.

## 2. Core contract: arena, bundles, and two preconditions

Freeze only what MIT already claims; never count an installed input as an output.

**Governing standard.** This program runs under the platform standing rule
[`dev_docs/CROSS_MODEL_TESTING.md`](../../../../../dev_docs/CROSS_MODEL_TESTING.md). Its operative
clauses, which Precondition A and the family modules implement: (i) **uniform default**, a
borrowed family is native and untwisted unless its author declares an internal representation or
a soldering prescription is supplied and pre-registered; (ii) **what counts as a declaration**,
an author's written statement or the model's canonical spec, NOT a component count and NOT an
agent's reading of the code (so this document's own reading of M4 is not a declaration); (iii)
**soldering clause**, a soldered or extended family is a different object and any result grades
"family + prescription P", owned by this borrowing column, not the family author (a column may
propose such an M8-owned extension for author and maintainer review, but never attribute the
prescription to the native family); (iv) **"not
applicable" is neutral** (🚧-weight, never a negative); (v) **author silence is terminal**, the
default stands and this pre-registration locks on schedule. No part of M8.2 blocks on an inbox.

**Arena.** S³/2I with the Möbius edge; the twisted-Laplacian edge level `2/R²` is verified
(M8.1). The three flat SU(2) connections are `σ ∈ {σ_0, Q, Q'}` (galois-pair § 2.3,
coexact-gap § 2.2): `σ_0` trivial, adjoint `τ_0 = 3R_0`; `Q` (McKay d=1), adjoint
`τ_Q = Sym²Q` (d=2); `Q'` (d=7), adjoint `τ_{Q'} = Sym²Q'` (d=6). Never `Sym²(1)`.

**Precondition A (compatibility, the index question).** A family tests three-connection
structure only if its native field carries a representation the flat connection acts on. A
field whose only index is geometric (spatial / frame / spacetime) does NOT automatically:
its T-c verdict is "not applicable" (native), and it becomes testable only under an explicit
soldering prescription that identifies the geometric index with an internal one, or a separately
defined adjoint-valued extension. A soldering or extension is additional model content belonging
to whoever explicitly defines it: an author-defined prescription may extend the native family,
while an M8-defined prescription is owned and graded by M8 as a separate object. Either way it is
not a neutral transport. On S³, parallelizability makes frame-to-internal soldering
mathematically AVAILABLE but does not select a CANONICAL framing (possible ≠ canonical): the
framing and locking prescription are model content, must be pre-registered, and different
prescriptions define different extended objects.

**Precondition B (background admissibility, the vacuum question).** The KINEMATIC
certification table (§ 3) is a property of the coefficient BUNDLE and needs no background. But
any DYNAMICAL statement (vacuum stability, fluctuation spectrum) is linearized about a
background, and in a nontrivial flat sector a background is admissible only if it is a global
section: either a constant fixed by the holonomy, or an equivariant (generally non-constant)
texture. A pointwise potential minimum that the holonomy does not fix is NOT automatically an
admissible background. This precondition is why the nontrivial-sector `M5 + P` fluctuation tower
is not yet computable (§ 6.2); native untwisted M5 has no such holonomy obstruction (though its
operator is still unresolved, § 6.2).

**Targets (T-a, structural ladder only).** The 8 nontrivial 2I irreps at McKay distances,
`j_first = d/2` (mass-spectrum § 3-4, verified). The 24-entry numeric mass table is OUT
(null `mass-null-v1.0`, p = 0.174).

**T-b1** (`2/R²`, edge) is a verified external benchmark, scored in M8.4 only with a frozen
edge operator + bulk-to-edge correspondence. **T-b2** (`4/R²` standard, `36/R²` Galois) is the
adjoint coexact ONE-FORM benchmark, scored only against an operator sub-block identified with
that Laplacian, and only after M8.1.1; 0-forms use `d(d+2)/R²`. **Timing guard:** if M8.1.1
verifies T-b2 before the first T-b2-eligible target-bearing run, it enters through a newly signed
module version; if verification lands afterward, it is a separate prospective test and cannot
alter the existing verdict.

## 3. The kinematic certification gate (per family, general formula)

At the level of the bare coefficient-bundle harmonic decomposition, the first-occurrence
structure is fixed independently of the potential and before any background is selected.
Reproducing it certifies the bundle implementation only. A nonconstant background may break or
reduce the symmetry of the physical fluctuation operator, so this table is not automatically its
spectral decomposition.

For a block carrying coefficient system `κ_σ` (the block's own internal representation pulled
through the connection), the first harmonic level of irrep ρ is the minimum over the
constituents of `ρ ⊗ κ_σ`:

```text
n^(kappa)_{rho,sigma} = min over alpha in (rho (x) kappa_sigma) of  entry(d_alpha)
```

`κ_σ` is per block, NOT universally adjoint: adjoint / vector block `κ_σ = τ_σ = Sym²σ`
(`3R_0` trivial, `R3` standard, `R4` Galois); trivial scalar block `κ_σ = R_0` (untwisted for
every σ); M5 spin-2 block `κ_σ =` its separately derived 5-dimensional coefficient
representation (σ-dependent, pending V1b). `entry(d)` is set by the operator class: 0-form
(function tower) `entry(d) = d` (eigenvalue `d(d+2)/R²`); coexact 1-form `entry(d) = d`
(`d≥2`), `2` (`d=0`), `3` (`d=1`). Tables are computed and self-checked in
[`scripts/m8_2_first_occurrence.py`](../scripts/m8_2_first_occurrence.py) (orthonormality
1.1e-15; `Sym²Q=R3`, `Sym²Q'=R4`; M4 reproduces mass-spectrum § 4 exactly). No universal
`λ = g(d)` law is frozen.

**This gate certifies the coefficient bundle only.** It presupposes precondition A (which
coefficient bundle) and says nothing about precondition B (which background). A passing table
does not imply an admissible vacuum or a defined fluctuation spectrum exists.

**Independent reproduction (M8.5).** M8.5 must reproduce every certification table through an
INDEPENDENTLY implemented decomposition. It may compare against `m8_2_first_occurrence.py` but
may not call it, import its tables, or share its derived fixtures.

## 4. The success ladder

| Test | What it establishes |
| --- | --- |
| The block's per-connection first-occurrence table reproduced | coefficient-bundle implementation CERTIFIED (kinematic; presupposes A, not B) |
| An admissible global background exists in the sector (const or texture) | precondition B MET (a computed output) |
| Stable vacuum about that background | dynamics ADMISSIBLE (a computed output) |
| One frozen action + coupling set works across the applicable sectors | three-connection COMPATIBILITY (a computed output) |
| The sectors' spectra differ | connection SENSITIVITY (a computed output) |
| Stable nonlinear branches / defects exist | a DYNAMICAL result |
| Their energies obey a separately pre-registered MIT relation | eventual OQ1 SUCCESS |

The first row is kinematic certification. The rest are computed physical outputs; only the
final two establish nonlinear realization and address OQ1 directly.

## 5. Outcome language, four separate axes

Recorded independently, never merged:

| Axis | Values |
| --- | --- |
| Index compatibility | complete / partial / not implementable / not applicable (geometric-only index) / applicable under declared soldering or extension |
| Background admissibility | admissible constant / admissible texture / no admissible background / not evaluated, per sector |
| Vacuum stability | stable / unstable / marginal / not evaluated, per sector |
| Connection sensitivity | sensitive / isospectral / connection-blind / not evaluable |

"Fails the pre-registered survey" = no declared grid/seed point passes; "refuted throughout the
declared range" is reserved for an analytic / certified-interval / exhaustive bound. All-family
failure is scoped to the frozen implementations and domains; OQ1 stays open.

## 6. Family modules

### 6.1 M4 (AUTHOR-DECLARED geometric displacement, @jeffsyee; native/untwisted, T-c N/A)

The pinned code (`wave_engine.py` @ `c9dc3796`): "a single 3-component field `ψ`" evolved by a
plain **vector Laplacian applied componentwise, no div/curl** (lines 29, 118), seeded as a
**radial (longitudinal) displacement `ψ = A·profile(r)·r̂`** from the domain center (lines
51-52, 84-96), with the potential acting on `u = ψ·ψ = ‖ψ‖²` (line 171). The radial
initialization and `r̂` alignment are signs of a physical SPATIAL displacement; the
componentwise Laplacian is index-agnostic. So the index is not established as internal from the
code alone. **RESOLVED by author declaration (@jeffsyee, M4 Q&A, 2026-07-24):** in EWT, ψ is the
geometric spatial displacement of the medium; its three components describe displacement in
physical space, and are NOT internal scalar components nor an internal index soldered to the
spatial one (the geometric-displacement reading, stated just below). Native M4 is therefore geometric and untwisted, T-c
NOT applicable, now by an accepted author declaration rather than the silence default. The
internal-triplet first-occurrence table below is consequently known NOT to describe native M4; it
is the fixture for a separate internally-valued REPLACEMENT model (M4_int, `ψ ∈ Ω⁰(X; E_{τ_σ})`).
A soldered geometric extension (M4 + P) is a different object again: its bundle, operator, and
first-occurrence table must be derived from the actual prescription P and are not automatically the
M4_int table. (Engine credit, per @jeffsyee: @xrodz built the initial M4 model, @lsmolinski
extended it, @jeffsyee gave the theoretical basis.)

The three candidate index readings are now settled by that declaration: **geometric displacement
(chosen)** = native M4, `ψ ∈ Γ(TX)` (a vector field on `X = S³/2I`), untwisted, T-c N/A;
**internal triplet** (rejected as native) would be a REPLACEMENT model M4_int (`ψ ∈ Ω⁰(X; E_{τ_σ})`,
the existing table); **spatial-to-adjoint soldering** would be a distinct M8-owned extension M4 + P
(bundle/operator/table derived from the prescription P, not the M4_int table). Neither is native M4.

| Item | Spec (otherwise complete) |
| --- | --- |
| Bundle | @jeffsyee DECLARED geometric: native M4 is a spatial displacement vector field (`ψ ∈ Γ(TX)`), untwisted; the flat connection does not act on it. The internal-triplet bundle `ψ ∈ Ω⁰(X; E_{τ_σ})` belongs to the internally-valued replacement model M4_int, NOT native M4; a soldered M4 + P has its own bundle derived from P |
| Native quotient operator (core-internal, REMAINING; the only open M4 question, now geometric not ontological) | native M4 is `ψ ∈ Γ(TX)`, flat R³ in the pinned code; transport to S³/2I needs a declared covariant operator (Bochner/rough on the tangent bundle vs Hodge on 1-forms vs frame-scalars vs a parallelization vs an EWT-motivated elastic operator with a longitudinal/transverse split), inequivalent by curvature terms on the curved quotient. Core-internal M8 SPECIFICATION work, frozen in the § 7 execution appendix; the M4 thread is closed on the settled ontology and the operator choice stays inside M8, not reopened as a question to @jeffsyee |
| Operator class | native M4: the geometric quotient operator above (pending, § 7). The twisted FUNCTION Laplacian + potential Hessian (0-form) is the M4_int operator, NOT native M4's |
| Certification table | `n_{ρ,σ}` via § 3; VERIFIED vs mass-spectrum § 4 (all 9 × 3). @jeffsyee declared geometric, so this is NOT native M4's fixture; it is the fixture for the internally-valued model M4_int. A soldered M4 + P needs its own table derived from P. Native M4 fills no T-c cell |
| Vacuum caveat (extensions only, precondition B) | native untwisted M4 has NO Q/Q' sectors, so no holonomy obstruction. For a TWISTED extension (M4_int / M4 + P): a nonzero `v_mode` minimum (`double-well`) is a fixed nonzero vector, generally NOT holonomy-invariant in the Q/Q' sectors, so its appendix must establish a global parallel section or an equivariant texture, not merely "pick the minimum" |
| Coupling taxonomy | action `c1,c2`; discrete `v_mode ∈ {linear, cubic-NLS, saturating, double-well}`; numerics in the § 7 appendix |
| Excluded | `v_mode = density_mod`: a state-dependent-feedback / medium model (potential on a dynamically computed `ρ_local` with clipping), NOT a fixed `V(ψ)`. Out of the Lagrangian survey |

M4 first-occurrence table, the fixture ONLY for the internally-valued model M4_int
(`ψ ∈ Ω⁰(X; E_{τ_σ})`); NOT a prediction or certification target of native M4 (declared geometric
by @jeffsyee), and NOT automatically the table of a soldered M4 + P (which needs its own
derivation from P). Level `n`, eigenvalue `n(n+2)/R²`:

| irrep | `d` | trivial | standard | Galois |
| --- | --- | --- | --- | --- |
| R_1 | 1 | 1 | 1 | 5 |
| R_3 | 2 | 2 | 0 | 4 |
| R_6 | 3 | 3 | 1 | 3 |
| R_7 | 4 | 4 | 2 | 2 |
| R_8 | 5 | 5 | 3 | 1 |
| R_4 | 6 | 6 | 4 | 0 |
| R_5 | 6 | 6 | 4 | 2 |
| R_2 | 7 | 7 | 5 | 3 |
| R_0 | 0 | 0 | 2 | 6 |

### 6.2 M5 (AUTHOR-CLARIFIED spacetime tensor, @JarekDuda; native/untwisted, T-c N/A; native operator selection still core-internal; "M5 + P" an M8-owned object)

Rebased off the current pinned canonical potential of record (superseding the older quartic-LdG).
This is the repository's potential of record, NOT asserted as a uniquely derived or fundamental M5
potential; @JarekDuda describes the deeper anisotropy-generating physics and possible extra
kinetic terms as open (see his follow-up below):

```text
M   real symmetric 4x4,  eta = diag(-1,1,1,1),  M -> O^-T M O^-1  (O in SO(1,3): rotations + boosts)
V(M) = w SUM_{p=1..4} ( Tr_eta(M^p) - C_p )^2 ,   C_p = g^p + 1 + delta^p
M_vac = diag(-g, 1, delta, 0),  V(M_vac)=0
lambda = spectrum(eta M) current ;  Lambda(vac) = spectrum(eta M_vac) = (g,1,delta,0)  (g=8, delta=0.3)
```

**Notation (author @JarekDuda, M5 Q&A 2026-07-25; heed the differences).** λ = current
eigenspectrum, Λ = the vacuum eigenspectrum the potential prefers, O = the transformation (SO(3)
at levels 1-2, SO(1,3) with boosts at level 3). Two cautions: (i) the pinned canonical writes the
transformation as Λ, so this doc uses O to avoid clashing with the vacuum-spectrum Λ; (ii) Jarek
phrased the preferred spectrum as "eigenspectrum of M", while the Lorentz-invariant statement (and
the canonical's convention) is the spectrum of η·M, which this doc uses. M5 has three modeling
levels: (1) EM, S² vacuum, vector `n`; (2) EM+QM, SO(3), 3×3, Λ=(1,δ,0); (3) EM+QM+GEM, SO(1,3),
4×4, Λ=(g,1,δ,0). M8 uses level 3.

**Default and the soldered extension.** Under the governing standard, native M5 is a rank-two
Lorentz-covariant tensor whose transform `M ↦ O⁻ᵀMO⁻¹` (O in SO(1,3)) is a spacetime/frame law, NOT a separate
internal gauge index, so by DEFAULT M5 is untwisted and T-c is not applicable (M5 stays fully in
the program for vacuum stability, spectrum, and defects). M5 becomes three-connection-testable
only through a SOLDERING: identifying the spatial-frame SO(3) with the adjoint holonomy
`Ad∘σ : 2I → SO(3)`. That is economical (no new components, unlike M7_ad) but is a spin-isospin
prescription, not a consequence of Lorentz covariance. Per the soldering clause it defines a
SEPARATE object, "M5 + P", scored under its own name; per the standing rule it may open EITHER
through a @JarekDuda declaration (author-motivated) OR through a separately-approved, M8-owned
pre-registered prescription (attributed to M8, never to native M5 or its author). **@JarekDuda
RESPONDED (M5 Q&A, 2026-07-25), not silent:** M5's field is a spacetime/frame tensor (frame
symmetry SO(1,3), his level 3), NOT an internal multiplet, which CONFIRMS native M5 untwisted /
T-c N/A; and he gave latitude to CHOOSE among the existing M5 levels and use M5 unmodified ("use M5 as it
suits you best, not modifying it"). That latitude is about level-choice, NOT a soldering: he did
not identify the frame SO(3)/SO(1,3) action with MIT's internal connection or prescribe a
soldering. An M8-owned "M5 + P" therefore remains M8-owned and pre-registered by M8 under the
standing rule (an overlay consistent with "not modifying it," not an author endorsement of any
soldering), never a native-M5 claim. Two things remain, at DIFFERENT scopes: native
M5's operator selection (below) is CORE-INTERNAL, required for native M5's own vacuum, spectrum,
and defects that the program tests even under T-c N/A; the background-admissibility and spin-2
items arise ONLY under the soldering (the M5 + P track). Neither is an author gate.

**M5 + P problem (background admissibility, precondition B; the largest math issue; soldered-sector
only).** The spatial part of `M_vac` is `diag(1, δ, 0)`, ANISOTROPIC (three distinct spatial eigenvalues). `Ad∘σ`'s
image is the icosahedral group acting IRREDUCIBLY on ℝ³, so the only holonomy-invariant
symmetric 2-tensors are isotropic (multiples of the identity). Therefore the constant
anisotropic `M_vac` is NOT a parallel background in the Q or Q' sectors. The twisted background
is one of: (i) no parallel pointwise vacuum of this type; (ii) a non-constant equivariant
texture; (iii) a configuration with unavoidable covariant-gradient energy; (iv) a different
invariant background; (v) partial sector incompatibility. V1b must answer "what section is the
fluctuation operator linearized about?" before any tower is computed.

**Native M5 requirement (operator selection; core-internal, NOT soldering-specific).** Required for
native untwisted M5's own vacuum, spectrum, and nonlinear dynamics that the program tests even
under T-c N/A; the same choice then carries into any M5 + P. The canonical specifies a WORKING potential (itself
an open choice, per the canonical and confirmed by @JarekDuda) and NOT a unique fluctuation operator. It records: the kinetic operator `K(M)` is DEGENERATE everywhere
(`Ṁ ∝ η` is an exact null, m5 § 2); `H` is INDEFINITE / unbounded below (m5 § 2); the free-EL
IVP is ILL-POSED (every regularization blows up, diagnostics only); the working formulation is
the free-period least-action BVP, AMENDED to require profile-dynamic (breathing) orbits (rigid
conjugation ruled out on the loop, M5.20.5); and the canonical positive-kinetic stack is
"documented as such, NOT the theory's dynamics." It further rules that a spectral comparison
"under the true L must use" `gen-eig(Hess_V, K10)` (the chirped vacuum ladder), NOT the
unit-inertia ladder. V1b must name the exact operator being certified: (1) the potential Hessian
alone; (2) the `gen-eig(Hess_V, K10)` degenerate-kinetic problem; (3) the BVP second variation;
(4) the regularized diagnostic stack; or several, clearly separated. These carry different null
spaces, signs, and spectra. **@JarekDuda follow-up (2026-07-26):** the author is not certain of the
potential details, holds the potential FORM as a choice "based on agreements" and effective (a
deeper physics presumably sets the anisotropy), and notes Skyrmion-type EXTRA KINETIC TERMS may be
worth including if the current kinetic proves insufficient. So the § 7 execution appendix pins the exact tested
potential and operator against the IMMUTABLE canonical version; any departure (a potential revised
with the author per his "by agreement" framing, or added Skyrmion-type kinetic terms) is an
M8-owned test convention, labeled as such. (The author names the EM/QM/GEM scale hierarchy, `0<δ≪1`, `g≫1`, as the deepest
open problem; that is M5-physics, not an M8.2 gate.)

**Field-space split (potential-level; matches the canonical's audited decomposition,
`m5_21_1e` / `m5_20_4`: "isospectral orbit class" vs "amplitude modes").** `V` depends only on
`spectrum(ηM)`, so the 10-dim fluctuation splits at the POTENTIAL level: 6 orbit-tangent
(V-flat, the O(1,3) conjugation orbit of `M_vac`; discrete stabilizer since the eigenvalues are
distinct) + 4 amplitude (V-stiff, eigenvalue-changing; the four `Tr_η(M^p)=C_p` conditions have
a nondegenerate Vandermonde Jacobian). This is a statement about `V` ALONE. The PHYSICAL
spectrum is the degenerate-indefinite generalized problem of the operator-selection requirement
and does NOT simply inherit this split.

Field-representation content under the (soldered) spatial SO(3), `Sym²(ℝ⊕ℝ³) = 1 ⊕ 3 ⊕ 1 ⊕ 5`:

| Block | SO(3) rep / coefficient | dim | at the diagonal vacuum (LOCAL frame statement) | tower |
| --- | --- | --- | --- | --- |
| time-time | trivial | 1 | amplitude | function, UNtwisted |
| time-space | vector `τ_σ` | 3 | orbit (boosts) | function, twisted → the M4 table |
| space trace | trivial | 1 | amplitude | function, UNtwisted |
| space traceless | spin-2 (irreducible 5) | 5 | 3 orbit (rotations) + 2 amplitude | function, twisted (spin-2) |

The spin-2 row's "3 + 2" is a LOCAL frame statement at the chosen anisotropic vacuum (whose
stabilizer is discrete); it is NOT an SO(3)-invariant subrepresentation, and the flat and stiff
parts do NOT own independent global twisted towers. The spin-2 certification table is for the
FULL 5-dimensional coefficient bundle; the fluctuation operator acts within and mixes it.

**M5 + P track, in order (opens on a @JarekDuda declaration OR an M8-owned prescription; spin-2
table is LAST):** (a) the soldering is declared (author) or approved (M8-owned), else the native
default holds and M8.2 locks without this track; (b) construct/check the admissible twisted
background (the M5 + P background problem); (c) fix the fluctuation operator (the native
operator-selection requirement, needed for native M5 anyway); (d) THEN compute the spin-2 twisted
first-occurrence table (extending the script) and numerically confirm the actual spectral gaps of
the selected operator; (e) freeze the "M5 + P" object and SHA-pin. The vector and scalar blocks'
coefficient content is settled (M4 table + untwisted); their physical spectra still wait on
(b)-(c). None of this blocks the native-M5 default lock, which needs only step (c).

### 6.3 Native M7 (T-c N/A under the default; column PARKED, held)

| Item | Spec |
| --- | --- |
| Native fields | spatial vector fields `A, J` (one-forms, 3 geometric components); div-based charge; massive `J`. Pinned functional uses curls, dot products, divergence charge, NO internal adjoint index |
| Bundle | `A, J ∈ Ω¹(X)`, NO internal index. The flat connection has nothing to act on |
| Linearized space | the full one-form perturbation space, decomposed into exact and coexact Hodge blocks; constraints + Hessian determine the admissible blocks |
| Operator class | coupled 1-form block (curl, `m_J²A·J` mass-mixing, `f(J·J)` Hessian); not the bare Hodge gap |
| T-c | NOT APPLICABLE: native `A, J` carry no internal representation for the flat connection |
| Removed | the twisted coexact first-occurrence table is NOT the native certification target (it certifies an adjoint-coexact sub-block the native field space lacks) |
| Retained | quotient vacuum stability, mode spectrum, nonlinear branches, defects, M4/M5 comparison |
| Vacuum | the native truncation has a known long-wavelength tachyonic band (threshold ≈ 0.786). Whether the frozen quotient realization retains/removes/shifts it is DETERMINED by the § 7 calculation, not pre-declared |

### 6.4 `M7_ad` (a separate object owned by this column; DEFERRED, not author-gated)

A distinct adjoint-valued extension `A, J ∈ Ω¹(X; E_{τ_σ})` (nine-component) could test the three
connections, but per the soldering clause it is a SEPARATE object needing its own action,
contractions, constraints, vacuum analysis, and identifier, and its prescription is M8-owned, not
the M7 authors'. DEFERRED by our own choice (a discretionary hold, not an author gate, since the
prescription is ours): the M7 column is parked (author away), so no M7 outreach goes out now and
native M7 records "not applicable" under the default meanwhile. If pursued later we may consult
@marcf999 / @pwerbos (M7's Ouroboros co-parent) for physical motivation, as a Q&A per
`CROSS_MODEL_TESTING.md` § 6, but the extension is ours to define; a diagonal identification only
if they motivate it.

## 7. Execution appendices (mandatory, pre-run, per family)

The core (§§ 2-5) and each module (§ 6) freeze the arena, bundles, operator classes,
certification tables, admissibility and outcome language, no-search rules, and applicable
benchmarks. They do NOT fix numerics. Before a family's first target-bearing run, a separately
SIGNED execution appendix freezes: the immutable code SHA; the exact action; coupling values +
scan grid; the ADMISSIBLE BACKGROUND (a global parallel section or a declared equivariant
texture, per precondition B, established rather than assumed from a potential minimum); the
selected fluctuation operator (the native M5 operator-selection requirement, and per family where
relevant); the seed / initial-condition
catalogue; continuation + stopping rules; the resolution ladder; convergence tolerances;
stability-test duration; and excluded branches.

No appendix may use target-bearing physical outcomes from its own family or another M8.4 family
to select couplings, seeds, branches, or verdict thresholds. Neutral platform-wide solver
validation, convergence studies, and hardware constraints may be reused if declared and
target-blind.

## 8. Definition of done (per module)

Under the governing standard M8.2 LOCKS on the default without author replies. The remaining
core-internal work is: source pins, plus the native quotient operators for M4, M5, and M7. Author
declarations or M8-owned prescriptions open OPTIONAL tracks that do not block the lock. T-c
classification is settled for all three families; the native dynamics modules are not.

| Module | Remaining to lock (default; core-internal) | Optional track (author declaration or M8-owned prescription) |
| --- | --- | --- |
| Core (§§ 2-5) | maintainer sign-off; the immutable source pins (§ 1, SHA-256 / archived deposit) | — |
| M4 (§ 6.1) | freeze the native quotient transport rule + geometric operator: the pinned code is a FLAT componentwise vector Laplacian, and transport to S³/2I is a choice (Bochner/rough vs Hodge vs frame-scalars vs a parallelization), inequivalent on the curved quotient. T-c stays N/A | @jeffsyee DECLARED geometric (2026-07-24), so the internal-triplet path is closed; the only optional track is an M8-owned soldered "M4 + P", a separate object |
| M5 (§ 6.2) | native operator selection (degenerate/indefinite kinetic; the `gen-eig(Hess_V, K10)` ruling), for native M5's vacuum/spectrum/defects. T-c stays N/A | @JarekDuda clarified spacetime/frame (T-c N/A) and left latitude to use M5 unmodified; an M8-owned "M5 + P" (background problem → spin-2 table → freeze; operator shared with native) opens only if M8 elects it, never attributed to Jarek |
| Native M7 (§ 6.3) | the full-one-form Hodge/constraint decomposition | `M7_ad` (§ 6.4), an M8-owned extension, deferred |
| Author outreach | M4 @jeffsyee ANSWERED 2026-07-24 (geometric, T-c N/A); M5 @JarekDuda ANSWERED 2026-07-25 (spacetime tensor, T-c N/A; an M8-owned M5+P is available under the rule, not endorsed by Jarek); M7 held (parked). Umbrella to Rodrigo posted and answered | — |
| Cross-doc | the OQ1 "energies vs stages" clarification proposed to the maintainer | — |

## Addenda (post-freeze only)

**Post-lock note, 2026-07-29.** The mass-table null-test citation in § 2 refers to the
superseded pre-correction run `mass-null-v1.0` (`p_A = 0.174`). The corrected-table
rerun, `mass-null-v1.1`, gives `p_A = 0.690`, strengthening rather than weakening the
statement that the 24-entry mass comparison is outside M8.2's evidential scope. No M8.2
target, input, computation, or verdict changes.
