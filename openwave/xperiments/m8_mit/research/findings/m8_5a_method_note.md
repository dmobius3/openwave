# M8.5-A method note

Deliverable of the context-isolated independent-method reproduction specified in
`PROTOCOL.md` (frozen 2026-07-30). Implementation: `m8_5a_reproduction.py`. Claim ceiling
(§ 2): context-isolated independent-method reproduction, nothing stronger. Adjudication
against the pinned § 6.1 table is the maintainers' step and is not performed here.

> **Status.** §§ A-H are the implementer's clean-room draft, byte-committed in the § 9
> commitment ([`m8_5a_commitment.md`](m8_5a_commitment.md), commit `dac2b6a1`, draft SHA-256
> `c8ac98a218c0449ff0dfa28eb280b879d9ce44c742c94cf13e7e635ea60ce6a8`); that version is
> retrievable at that commit, and nothing in §§ A-H has been edited since. The maintainer
> layer (§§ I-J: the adjudication outcome and the recorded adversarial audit) was appended
> 2026-08-01, which is what moves this note off DRAFT per the landing map in
> [`../tasks/m8_5_task_details.md`](../tasks/m8_5_task_details.md). The § B line references
> resolve against [`../scripts/m8_5a_reproduction.py`](../scripts/m8_5a_reproduction.py) as
> committed at `dac2b6a1` (the file is frozen; any later edit would be a dated rerun per § 8).

## A. Equations first

### A.1 The group

The packet supplies two unit quaternions `g_1, g_2` with components in `Q(phi)`,
`phi^2 = phi + 1`. The group is the multiplicative closure

```text
G = <g_1, g_2>  under quaternion multiplication in the basis (1, i, j, k).
```

All arithmetic is exact: an element of `Q(phi)` is stored as `a + b*phi` with `a, b`
exact rationals (`fractions.Fraction`); products reduce by `phi^2 = phi + 1`. Conjugacy
classes are orbits under conjugation by the generators (which generate all inner
automorphisms). Inverses use `q^{-1} = conj(q)` for unit quaternions.

### A.2 Inner product and multiplicities

```text
<a, b>_G = (1/|G|) * sum_g a(g) * conj(b(g))     (multiplicities over C)
```

Every character computed in this run is exactly real in `Q(phi)` (traces of the complex
matrices are asserted to have imaginary part exactly zero), so conjugation is the
identity and dual conventions cannot flip any result. Integer-nearness (§ 5.4) is
enforced at tolerance exactly 0: a multiplicity must BE a nonnegative rational integer
(zero `phi`-part, denominator 1). The occurrence test `> 1/2` is applied to the exact
integer only after that gate passes.

### A.3 The scalar tower, exactly

For `g` with `2 cos(alpha) = 2 Re(g) =: t` (an exact element of `Q(phi)`), the power
sums `s_m = e^{i m alpha} + e^{-i m alpha}` satisfy

```text
s_0 = 2,   s_1 = t,   s_m = t*s_{m-1} - s_{m-2}
```

and the restricted `SU(2)` character is the sum of its `n+1` weights:

```text
chi_{V_n}(g) = sum_{k=0..n} e^{i(n-2k) alpha}
             = s_n + s_{n-2} + ... + (s_1 or 1)      (exact in Q(phi))
```

### A.4 Irreducible characters by peeling (the character-extraction method)

Process `n = 0, 1, 2, ...`: subtract from `chi_{V_n}|_G` all known irreducible content
(each coefficient an integer-gated multiplicity). The remainder is a genuine character,
a nonnegative integer combination of not-yet-found irreducibles. If its norm
`<r, r>_G` is exactly 1, `r` IS a new irreducible character. Norm >= 2 remainders are
pooled and re-reduced whenever a new irreducible lands. Extraction completes when

```text
sum_rho (dim rho)^2 = |G|
```

and fails loud if the tower bound `n <= 24` is exhausted first. Correctness needs no
assumption about the group: restrictions of `V_n` span the class functions or the run
fails loud; a norm-1 nonnegative-integer combination of irreducibles is one irreducible.

### A.5 Connections, coefficients, and the frozen embedding

`Q` is identified, never declared (§ 5.3): it is the unique 2-dimensional irreducible
whose character equals `2 cos(theta(g)) = 2 Re(g)` on every class under the packet's
embedding (gate G5). `Q'` is the unique other 2-dimensional irreducible. The explicit
matrices used for the three declared flat-connection classes are:

| sigma | Explicit matrices |
| --- | --- |
| `trivial` | the 1x1 identity on every class |
| `standard` (`Q`) | `su2_matrix(q)`: the packet embedding `q = w+xi+yj+zk -> [[w+xi, y+zi], [-y+zi, w-xi]]` |
| `galois` (`Q'`) | `su2_matrix(galois(q))`: componentwise `phi -> 1-phi`, a quaternion-algebra automorphism, hence a genuine 2-dim representation of the same abstract group; G5 verifies its character is the other 2-dim irreducible |

The coefficient representation is the contract input `tau_sigma = Sym^2(sigma)`,
constructed EXPLICITLY (TASK requirement 2, preferred route): the induced matrix of
`sigma(g)` on the monomial basis `{e_i e_j : i <= j}` is built entry by entry and
`chi_{tau_sigma}(g)` is its trace. The character identity
`(chi^2(g) + chi(g^2))/2` is NOT used in construction; it is checked afterwards by G11
on the same constructed object, with `g^2` evaluated by raw quaternion multiplication.
The two routes genuinely differ, so G11 is a cross-check here, not a restatement.

### A.6 The reproduced quantity

```text
n_first(rho, sigma) = min{ n >= 0 : <chi_{V_n}|_G * chi_{tau_sigma}, chi_rho>_G > 0 }
```

searched for `n <= 24`, NOT-FOUND failing loud. Row identity is label-free: each
irreducible is reported by `(dim, McKay distance)`, where `A_{rho,tau} =
<chi_rho * chi_Q, chi_tau>_G` is derived in-implementation, integer-gated, and the
distance is BFS graph distance from the trivial node over the support `A > 0`.

## B. Equation-to-code map

Line numbers refer to `m8_5a_reproduction.py` at commitment (SHA-256 in
`environment.md`). No repository permalinks exist in this clean room; file plus line
is the stable reference.

| Equation / object | Code |
| --- | --- |
| `Q(phi)` exact arithmetic, Galois map `phi -> 1-phi` | `QPhi`, `m8_5a_reproduction.py:65` |
| Quaternion product, conjugate, norm | `qmul`/`qconj`/`qnorm2`, `m8_5a_reproduction.py:152` |
| Packet parsing plus SHA-256 pin check | `parse_component`/`load_packet`, `m8_5a_reproduction.py:183` |
| A.1 closure `G = <g_1, g_2>` | `generate_group`, `m8_5a_reproduction.py:220`; `check_closure`, `:243` |
| Conjugacy classes as generator-conjugation orbits | `conjugacy_classes`, `m8_5a_reproduction.py:252` |
| A.5 frozen embedding `q -> SU(2)` | `su2_matrix`, `m8_5a_reproduction.py:310` |
| A.5 explicit `Sym^2` matrix on monomial basis | `sym2_matrix`, `m8_5a_reproduction.py:329` |
| A.5 sigma matrices and consumed `tau_sigma` characters | `build_sigma_reps`, `m8_5a_reproduction.py:354`; `build_tau_chars`, `:367` |
| A.3 exact tower via `s_m` recursion over weights | `build_tower`, `m8_5a_reproduction.py:394` |
| A.2 inner product; integer-nearness at tolerance 0 | `make_ip`, `m8_5a_reproduction.py:434`; `multiplicity`, `:449` |
| A.4 peeling extraction | `extract_irreducibles`, `m8_5a_reproduction.py:464` |
| McKay matrix, BFS distance, bipartite witness | `build_mckay`, `m8_5a_reproduction.py:540`; `bfs_distances`, `:568`; `bipartition_ok`, `:585` |
| A.6 first occurrences (scalar) | `scalar_first_occurrences`, `m8_5a_reproduction.py:606` |
| C. coexact module | `coexact_first_occurrences`, `m8_5a_reproduction.py:630`; `coexact_rule_prediction`, `:656` |
| § 7 comparison harness plus doc_typo mutation | `compare_tables`, `m8_5a_reproduction.py:669`; `perturb_rows`, `:697` |
| Pipeline and all gate evaluations G1..G12 | `run_pipeline`, `m8_5a_reproduction.py:714` |
| § 8 mutation harness (12 mutations, coverage enforced) | `MUTATIONS`, `m8_5a_reproduction.py:986`; `run_mutation_tests`, `:1002` |

## C. The coexact module: own harmonic analysis and the general argument

### C.1 The coexact tower on S^3 (derived, not imported)

Peter-Weyl on `S^3 = SU(2)` under left x right translation:

```text
L^2(S^3) = sum_n  V_n (x) V_n
```

Trivializing `Omega^1` by the left-invariant coframe gives
`Omega^1 = L^2 (x) su(2)*`, with right translation acting on `su(2)*` as `V_2` and left
translation acting on functions only, so as a left x right representation

```text
Omega^1 = sum_n V_n (x) (V_n (x) V_2) = sum_n V_n (x) (V_{n-2} + V_n + V_{n+2}).
```

Repeating with the right-invariant coframe forces the symmetric-in-labels form: the
blocks of `Omega^1` are exactly `V_a (x) V_b` with `|a-b| in {0, 2}`, each once.
Exact forms are `d` of nonconstant functions: the `(n, n)` blocks, `n >= 1` (and
`(0,0)` is absent from `Omega^1`). `S^3` has no harmonic 1-forms. What remains is
coexact:

```text
coexact Omega^1 = sum_{m >= 2}  ( V_m (x) V_{m-2}  +  V_{m-2} (x) V_m ).
```

On a 3-manifold, coexact eigenforms are curl eigenmodes and `Delta = (*d)^2`. On the
Maurer-Cartan block (`(0,2)`, level `m = 2`) the curl eigenvalue is `±2/R`, giving
`Delta = 4/R^2`, which is the known Killing-form eigenvalue on `S^3` and fixes the
normalization; the ladder of blocks gives `curl = ±m/R` on the level-`m` pair, so

```text
Delta = m^2 / R^2   on   (V_m (x) V_{m-2}) + (V_{m-2} (x) V_m),   m >= 2,
```

with multiplicity `2(m^2 - 1)`, consistent with the level-2 count of 6. The deck group
acts by LEFT multiplication, so only left factors matter for the quotient; the right
factors contribute nonzero multiplicities, which cannot change a first occurrence
(same waiver as § 5.2). Level `m` therefore contributes `V_m|_G + V_{m-2}|_G`, and the
§ 5.4 pattern gives

```text
m_first(rho, sigma) = min{ m >= 2 :
    <chi_{V_m} chi_{tau_sigma}, chi_rho> > 0  or  <chi_{V_{m-2}} chi_{tau_sigma}, chi_rho> > 0 }
```

searched for `m <= 24`. Convention map: this implementation's level `m` is the source
convention's coexact level `m`; both print eigenvalue `m^2/R^2`. The map is the
identity and is stated explicitly per § 6.

### C.2 The general argument for the entry rule

The asserted rule references only the McKay distance `d`, so it is adjudicated against
the trivial-connection column; the other two columns are reported as data.

**Lemma 1 (first occurrence equals distance).** Let `v_n(rho) = <chi_{V_n}|_G, chi_rho>`.
Clebsch-Gordan on `SU(2)` gives `chi_{V_n} = chi_{V_{n-1}} chi_Q - chi_{V_{n-2}}`
exactly, hence on multiplicity vectors

```text
v_n = A v_{n-1} - v_{n-2},   v_0 = e_trivial,   v_1 = e_Q.
```

By induction `v_n` is supported on distance <= n (applying `A` extends support by one
step; subtraction cannot extend it). At `n = d(rho)`:
`v_d(rho) = (A v_{d-1})(rho) - v_{d-2}(rho)`; the second term vanishes (distance
`d > d-2`), and the first is positive because some geodesic neighbor `tau` of `rho` at
distance `d-1` has `v_{d-1}(tau) > 0` by induction and `A_{rho,tau} > 0`. Hence
`n_first(rho, trivial) = d(rho)` for every irreducible. (Computational witness:
`n_first_trivial_equals_mckay_distance = True` in the module block.)

**Lemma 2 (parity).** The support graph is bipartite (witness: `graph_bipartite =
True`; the BFS parity classes 2-color every edge). By the same recursion,
`v_n(rho) = 0` whenever `n != d(rho) (mod 2)`.

**Rule.** `m_first(rho, trivial) = min{ m >= 2 : v_m(rho) > 0 or v_{m-2}(rho) > 0 }`:

| `d(rho)` | Argument | Result |
| --- | --- | --- |
| `d >= 2` | `v_d(rho) > 0` (Lemma 1) and no `m < d` can reach distance `d` (support bound) | `m_first = d` |
| `d = 0` | at `m = 2`, `v_0(rho) = v_0(trivial) = 1 > 0` | `m_first = 2` |
| `d = 1` | `m = 2` needs `v_2(rho)` or `v_0(rho)`; both vanish (`v_0` by distance, `v_2` by Lemma 2 parity); at `m = 3`, `v_1(rho) > 0` | `m_first = 3` |

This establishes the rule for the full representation family from the operator and
representation structure (the coexact tower's left content plus the restriction
recursion). It is not a pattern fit to the computed range: the computed `m <= 24` table
is the check of the argument, not its source. Both lemma witnesses are computed and
recorded in the module block; the derivation itself is what the § 9 adversarial audit
examines. Verdict claimed: **structurally derived and reproduced** (trivial-column
cells all match; the general argument above). Per the § 6 standing rule, the numerical
match contributes nothing to standing; only the derivation is on offer, and its
weakest link is stated in § H.

## D. Results beside their gates

Full verbatim output in `raw_output.txt`; machine-readable form in `result.json`.

| Result | Gate that binds it | Outcome |
| --- | --- | --- |
| Group order 120, closed | G1 | ✅ PASS |
| 9 classes, sizes [1, 1, 12, 12, 12, 12, 20, 20, 30] | G2 | ✅ PASS |
| 9 irreducibles, exactly orthonormal | G3 (tolerance 0) | ✅ PASS |
| Dimensions [1, 2, 2, 3, 3, 4, 4, 5, 6], squares sum to 120 | G4 | ✅ PASS |
| `Q` identified by `chi = 2 cos theta`; `Q'` the unique other 2-dim | G5 | ✅ PASS |
| McKay matrix: mark condition `A.dims = 2.dims` | G6 | ✅ PASS |
| Distances defined; 9 signatures pairwise distinct | G7 | ✅ PASS |
| All 811 multiplicities exactly nonnegative integers | G8 (tolerance 0) | ✅ PASS |
| `A` symmetric | G9 | ✅ PASS |
| Comparison harness reddens under doc_typo transcription mutation | G10 | ✅ PASS |
| Consumed `tau_sigma` satisfies the Sym^2 character identity classwise | G11 (Addendum 1) | ✅ PASS |
| Consumed tower satisfies `chi_{V_n}(e) = n+1`, `0 <= n <= 24` | G12 (Addendum 1) | ✅ PASS |
| Mutation harness: 12 mutations, every gate covered and reddened, clean run green | § 8 requirement | ✅ PASS (exit 0) |

Scalar first-occurrence table (label-free rows, `(dim, McKay distance)`):

| dim | distance | trivial | standard | galois |
| --- | --- | --- | --- | --- |
| 1 | 0 | 0 | 2 | 6 |
| 2 | 1 | 1 | 1 | 5 |
| 3 | 2 | 2 | 0 | 4 |
| 4 | 3 | 3 | 1 | 3 |
| 5 | 4 | 4 | 2 | 2 |
| 6 | 5 | 5 | 3 | 1 |
| 3 | 6 | 6 | 4 | 0 |
| 4 | 6 | 6 | 4 | 2 |
| 2 | 7 | 7 | 5 | 3 |

Coexact module: all trivial-column cells match the asserted rule; full `m_first` table
for all three columns is in `result.json` and `raw_output.txt`. Verdict:
structurally derived and reproduced (§ C.2), scope stated there.

## E. Method-overlap disclosure (§ 3, required)

| Ingredient | Overlap with object 1 (McKay recursion, literal group-theory inputs) | Overlap with object 2 (explicit quaternions, brute-force classes, Burnside) |
| --- | --- | --- |
| Explicit quaternion arithmetic and brute-force conjugacy classes | none | overlaps: same generic route |
| Character extraction by exact peeling of tower restrictions | none: no group-theory literals imported; characters derived, not assumed | differs: object 2 used Burnside class-sum diagonalization; this run never forms class-sum matrices |
| Clebsch-Gordan recursion used in Lemma 1 (coexact argument only) | the recursion is the same `SU(2)` identity the McKay recursion rests on; here it is a proof device, and the computed tables never use it (the tower is built from weight sums) | none |
| Exact `Q(phi)` arithmetic end to end | differs (object 1 numeric/literal) | differs (object 2 numpy float, seed 11) |

No call, import, copy, read, or execution of either object; no shared derived fixtures.
No third fundamentally different method is claimed; context isolation and
implementation independence are.

## F. Tolerances (§ 8 commitment, fixed before any comparison)

Every tolerance is fixed in the source (`TOLERANCES`, `m8_5a_reproduction.py:44`) at
exactly 0: the entire pipeline is exact over `Q(phi)`; no floating point exists
anywhere in it. Integer-nearness means "is exactly a nonnegative rational integer".
The § 8 stability-under-increased-precision requirement is discharged by exactness:
there is no working precision to increase and no conditioning to degrade; determinism
was additionally confirmed by a byte-identical double run. Tolerances cannot be
widened after unsealing because there are none to widen; any change would be a dated
rerun per § 8.

## G. Inspection artifacts (4)

| # | Artifact | What to inspect |
| --- | --- | --- |
| 1 | `raw_output.txt` | verbatim stdout: gate lines, both tables, mutation-harness lines |
| 2 | `result.json` | fixed § 5 schema, machine-readable |
| 3 | G5/G11 gate lines in `raw_output.txt` | the identified-not-declared `Q` and the consumed-object Sym^2 check |
| 4 | mutation-harness block in `raw_output.txt` | every gate reddened by a deliberate defect; coverage enforced; exit 0 |

## H. Claims NOT verified (§ 9, explicit)

| Claim | Status |
| --- | --- |
| The contract choice `tau_sigma = Sym^2(sigma)` | not verified; its consequences are checked (G11), the choice itself is the lock's § 2 contract input |
| Fitness of any simulation backend | out of scope; that is M8.5-B |
| The M8.3 Reidemeister-torsion closed forms | out of scope (§ 0) |
| The coexact entry rule's truth beyond this derivation | the § 6 verdict claimed is "structurally derived and reproduced", so the rule is not listed as unverified on that ground; the derivation's weakest links are stated below and stand or fall with the adversarial audit |
| Weakest links of § C: the identification of the coexact blocks as `V_m (x) V_{m-2} + V_{m-2} (x) V_m` and the curl normalization `Delta = m^2/R^2` rest on this note's own harmonic analysis (cross-checked against the Killing level `m = 2` count and eigenvalue, and the multiplicity formula `2(m^2-1)`); the rule's scope was read as the trivial-connection column because the rule references only `d` | flagged for the audit |
| That the packet generators are the maintainers' intended input | not verifiable from inside the clean room; G1/G2 bind order and structure, the packet audit (§ 4) binds intent |
| Agreement with the pinned § 6.1 table | not checked here by design; that is the § 7 adjudication, quarantine intact |

## Adversarial-audit note

This draft precedes the recorded adversarial audit required by § 9 before the
deliverable crosses to the author. Items explicitly queued for that audit: the § C.2
derivation (both lemmas and the scope reading), the realness assertion route, and the
peeling algorithm's completion argument.

## I. Adjudication outcome (maintainer layer, 2026-07-31)

Recorded here so the note carries its own result; the § 7 record is canonical.

| Item | Outcome |
| --- | --- |
| scalar table vs pre-registration § 6.1 at `ec877ee0` | **REPRODUCED**: 9 label-free signatures pairwise distinct both sides, 27/27 cells equal, exact integer comparison, no tolerance |
| three-way agreement | holds: this run, § 6.1 (object 1's published table), object 2's reconstruction |
| G10 under the maintainer harness | both transcription mutations redden the comparison |
| harness + record | [`../scripts/m8_5a_adjudication.py`](../scripts/m8_5a_adjudication.py), [`../data/m8_5a_adjudication.json`](../data/m8_5a_adjudication.json), merged separately (#392) after the commitment (#391), ordering auditable |
| claim label | context-isolated independent-method reproduction, the § 2 ceiling; the run does not earn "blind" and does not claim it |

## J. Adversarial audit (§ 9, recorded 2026-08-01)

Run by an independent second agent briefed to REFUTE, with its own scripts and its own
methods: its group construction is an explicit icosian list checked equal to the packet
closure, its classes are trace fibers verified as full conjugation orbits, its character
table comes by Burnside class-algebra splitting over exact `Q(φ)` (the route this
implementation deliberately does not use), and its tower, McKay matrix, and both
first-occurrence tables are its own machinery. Script:
[`../scripts/m8_5a_audit.py`](../scripts/m8_5a_audit.py) (headless, exit 0 only if nothing
is refuted; re-run green by the maintainer); record:
[`../data/m8_5a_audit.json`](../data/m8_5a_audit.json).

| # | Claim attacked | Verdict |
| --- | --- | --- |
| 1 | Lemma 1, `n_first(ρ, trivial) = d(ρ)` | ✅ CONFIRMED: `v_{d-2}(ρ) = 0` exactly and the geodesic-predecessor term of `(A v_{d-1})(ρ)` is strictly positive; nonnegativity of every `v_n` holds because each is a genuine restriction-multiplicity vector, independent of the recursion. All 9 irreducibles verified from the audit's own table |
| 2 | Lemma 2, parity vanishing | ✅ CONFIRMED, and strengthened: bipartiteness is a THEOREM here, not just a computed witness (`−1` is central, `V_n(−1) = (−1)^n`, every irrep carries central sign `(−1)^d`, every McKay edge flips it) |
| 3 | The rule table, all three rows | ✅ CONFIRMED, including the `d = 1` row (exactly one irrep, `Q`: `m = 2` excluded by distance AND parity, `m = 3` fires via `v_1(Q) = 1`). The only graph-specific inputs, connectivity and bipartiteness, are themselves derived, so the protocol § 6 "full family from operator and representation structure" bar is met |
| 4 | The coexact tower, `Δ = m²/R²`, multiplicity `2(m²−1)` | ✅ CONFIRMED by a route the note does not use: Casimir normalization anchored on the function sector, with the level-2 anchor `4/R²` recomputed from raw quaternion algebra; Peter-Weyl block bookkeeping verified from exact Clebsch-Gordan identities |
| 5 | The trivial-column scope reading | ✅ CONFIRMED as faithful: every source statement of the asserted rule is per constituent, the per-constituent entry level IS `m_first(ρ, trivial)`, and every standard/galois cell was verified to be the min-over-constituents with no extra rule content |
| 6 | Realness of all characters | ✅ CONFIRMED as a theorem (all classes inverse-closed) and by three computational routes (inverse closure, class algebra splitting over the real field, Frobenius-Schur indicators all `±1`) |
| 7 | Peeling completion | ✅ CONFIRMED: rank 9 already at `n ≤ 8` (nine distinct trace values, monic polynomials), norm-1 remainders necessarily irreducible; the pool replay on the audit's own data resolved at `n = 8` with no deadlock |
| 8 | Run-record consistency | ✅ CONFIRMED: all 27 scalar cells AND all 27 coexact cells of `result.json` match the audit's independently derived tables, label-free |

**Standing consequence.** The § C.2 derivation held under independent attack, which is what
the `structurally derived and reproduced` verdict was waiting on. Per the § 6 standing rule
the numerical match contributed nothing; the derivation carried it.

**Weaknesses the audit found that this note did not self-flag** (recorded per the audit's
report; none affects a result):

| Weakness | Severity |
| --- | --- |
| § A.5's phrasing invites the false reading that the Galois map fixes the group as a set; it maps it onto the twin icosian copy. The note's actual justification (ring automorphism, hence a genuine representation of the same abstract group) is correct and does not need set-stability, but the subtlety is unstated | informational |
| Lemma 1 as written conflates two inductions (support on `n`, positivity on `d`); the repaired two-induction form is what proves it | minor, exposition |
| Lemma 2's "by the same recursion" leaves base-case parity and connectivity implicit, and cites the computed bipartite witness where the one-line central-character proof exists | minor |
| The `d ≥ 2` row's "no `m < d` can reach distance `d`" silently needs both `v_m` and `v_{m−2}` to vanish; true, but compressed | minor |
| The peeling pool has no a-priori termination guarantee for a general group; only the fail-loud bound protects it. In this run it resolved at `n = 8`; the residual risk is unstated in § A.4 | minor |
| The "`> 1/2` occurrence test applied to an exact integer" is just `> 0` | cosmetic |
