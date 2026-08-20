# M8.8 Method-and-Gate Manifest

## 1. Selected route

**Route class (§6):** Algebraic Reidemeister torsion from the supplied finite based chain
complex, by determinant data. No spectral, zeta-function, heat-kernel, or closed-form
torsion fixture is used.

**Construction within the class:** For each nontrivial irreducible representation ρ of 2I
(dimension d):

1. Evaluate boundary maps d₁, d₂, d₃ under g ↦ ρ(g) to obtain twisted boundary matrices
   D₁(ρ) ∈ ℂ^{2d×d}, D₂(ρ) ∈ ℂ^{2d×2d}, D₃(ρ) ∈ ℂ^{d×2d}.
2. Verify ∂∂ = 0 in the twisted complex.
3. Verify acyclicity (ranks: D₁ has rank d, D₂ has rank d, D₃ has rank d).
4. Compute torsion via the alternating product of determinants of maximal nonsingular
   minors:

   Choose index sets J₁ ⊂ {1,...,2d} (|J₁|=d) such that D₁[J₁,:] is nonsingular.
   Let K₁ = {1,...,2d} \ J₁.
   Choose J₂ ⊂ {1,...,2d} (|J₂|=d) such that D₂[J₂, K₁] is nonsingular.
   Let K₂ = {1,...,2d} \ J₂.
   Verify D₃[:, K₂] is nonsingular.

   Then: τ_ρ = ±Δ₂ / (Δ₁ · Δ₃)

   where Δ₁ = det(D₁[J₁,:]), Δ₂ = det(D₂[J₂, K₁]), Δ₃ = det(D₃[:, K₂]).

5. T²(ρ) = |τ_ρ|² = |Δ₂|² / (|Δ₁|² · |Δ₃|²).

This formula is independent of the choice of J₁, J₂ (standard property of Reidemeister
torsion for a based acyclic complex). The sign ± does not affect |τ|².

**Representation construction:** All 9 irreducible representations of 2I are constructed
from the quaternion-to-SU(2) embedding (the fundamental 2-dimensional representation)
and its Galois conjugate (φ ↦ 1−φ), using the McKay graph (extended Ẽ₈) recursion:

- ρ₀: trivial (dim 1)
- ρ₁: fundamental, from quat→SU(2) map (dim 2)
- ρ₂ = Sym²(ρ₁) (dim 3)
- ρ₃ = V₂⊗ρ₂ − ρ₁ (dim 4)
- ρ₄ = V₂⊗ρ₃ − ρ₂ (dim 5)
- ρ₅ = V₂⊗ρ₄ − ρ₃ (dim 6)
- ρ₇: Galois conjugate of ρ₁ (dim 2)
- ρ₆ = V₂⊗ρ₇ (dim 4)
- ρ₈ = V₂⊗ρ₄ − ρ₃ − ρ₆, or equivalently from §5 of McKay graph (dim 3)

The internal labeling ρ₀–ρ₈ is provisional; final identification uses the §5.5 public
row signature (dimension + characters on s, t, st).

**Exact arithmetic:** All computations in Q(φ,𝐢) = Q((1+√5)/2, √(−1)). Elements of Q(φ)
are pairs (a,b) with a,b ∈ Q, representing a + bφ. Complex values are pairs (x,y) with
x,y ∈ Q(φ), representing x + y𝐢. The squared modulus |z|² = x² + y² ∈ Q(φ) for
z = x + y𝐢.

## 2. Declared native orientation (§5.4)

The native orientation is: T²(ρ) = |τ_ρ|² where τ is the Reidemeister torsion of the
based acyclic complex C_* ⊗_{Z[2I]} V_ρ in the declared bases, computed via the
alternating-product-of-minor-determinants formula above.

This orientation is whatever the formula gives. Under the Cheeger-Müller theorem, it
should equal the analytic T²_RS or its inverse (the global sign of log T² being the
§5.4 bridge). The anchor rule at R7 selects the correct orientation after reveal.

## 3. Conventions consumed

All conventions are read from the construction packet and applied exactly:

| Convention | Value | Source |
| --- | --- | --- |
| module_side | left | construction packet `basing.module_side` |
| vector_convention | row | construction packet `basing.vector_convention` |
| boundary_direction | right action (c · d_n) | construction packet `basing.boundary_direction` |
| evaluation | g ↦ ρ(g), no inverse/transpose/dual | construction packet `basing.evaluation` |
| augmentation | ε: every group element ↦ 1 | construction packet `basing.augmentation` |
| degree_range | [0, 3] | construction packet |
| free_ranks | [1, 2, 2, 1] | construction packet |
| basis_order | degree 1: (e_s, e_t); degree 2: (f_1, f_2) matching relators | construction packet `basing.basis_order` |
| top_closure | single 3-cell, positively oriented, generator of C₃ | construction packet `top_closure` |

**Evaluation procedure:** A boundary map entry ∑ cₖ · gₖ ∈ Z[2I] is evaluated under
representation ρ as ∑ cₖ · ρ(gₖ), a d×d complex matrix. The twisted boundary matrix
D_n(ρ) has block structure r_n × r_{n-1} with each block d×d. The action is: a row
vector c ∈ ℂ^{r_n · d} maps to c · D_n(ρ) ∈ ℂ^{r_{n-1} · d}.

## 4. Instantiated gates

### 4.1 Model gates (§9)

| ID | Gate | What it establishes | Pre-impl status |
| --- | --- | --- | --- |
| M1 | ∂₂∂₁ = 0 and ∂₃∂₂ = 0 over Z[2I] | C_* is a chain complex | PASS (validate_manifest.py) |
| M2 | Free ranks [1,2,2,1], χ = 1−2+2−1 = 0 | Closed 3-manifold, not presentation 2-complex | PASS (validate_manifest.py) |
| M3 | H_*(Z ⊗_{Z[2I]} C_*) ≅ (Z, 0, 0, Z) | Homology of S³/2I | deferred to production |
| M4 | H_*(C_*) ≅ (Z, 0, 0, Z) as Z-modules, with integral saturation certificates for im ∂₁, im ∂₂, im ∂₃ | Universal cover is S³ | deferred to production |
| M5 | Terminal map is augmentation ε | Declared augmentation | PASS (validate_manifest.py) |
| M6 | ∂₁ matches frozen 1-cell correspondence; ε∂₁ = 0 | Generator correspondence | PASS (validate_manifest.py) |
| M7 | Per-irrep acyclicity: R0 non-acyclic (PASS), every nontrivial irrep acyclic (PASS) | Theorem's hypothesis | deferred to production |

**M4 saturation method:** Hermite Normal Form of the Z-expanded boundary matrices.
The Z-expansion of d_n (a matrix over Z[2I]) is an integer matrix of size
(r_n · 120) × (r_{n-1} · 120). The HNF diagonal entries give the elementary divisors;
all-1 diagonals certify saturation. A maximal minor of determinant ±1 is the acceptance
certificate.

### 4.2 Theorem-side gates (§6, §9)

| ID | Gate | What it establishes |
| --- | --- | --- |
| T1 | Invariant positive-definite Hermitian form for each consumed representation, constructed by group averaging | Unitary hypothesis of Cheeger-Müller |
| T2 | Row identity resolved via §5.5 public row signature: dim + characters on s(118), t(80), st(71) as Q(φ) triples | Same flat bundle per label as analytic route |
| T3 | Convention fixture: synthetic exact representation and chain-complex instance, non-unitary where required, processed through the same evaluation code, GREEN under declared conventions, each convention mutated separately with each mutation reddening at least one gate | Character-invisible convention errors |

**T3 fixture design:**

The fixture uses TWO test instances:

**Instance A (evaluation mutation):** Group Z/5 = ⟨g | g⁵ = 1⟩ with a 2-dimensional
representation ρ(g) = P · diag(ω, ω²) · P⁻¹ where ω = e^{2πi/5} and P is a
non-unitary matrix. The character χ(g) = ω + ω² is non-real.

- Gate: character of g under correct evaluation equals ω + ω² (non-real).
- Mutation (g ↦ ρ(g⁻¹)ᵀ): character becomes ω³ + ω⁴ ≠ ω + ω² → GATE FAILS.

**Instance B (boundary direction, vector convention, module side):** The actual 2I chain
complex with a specific nontrivial irrep (e.g., the 2-dimensional fundamental).

- Gate: D₂(ρ) · D₁(ρ) = 0 (twisted ∂∂ = 0).
- Boundary-direction mutation (transpose block structure): the matrix product reverses
  order → ∂∂ ≠ 0 → GATE FAILS.
- Vector-convention mutation (full transpose): equivalent to reversing row/column
  semantics → ∂∂ ≠ 0 → GATE FAILS.
- Module-side mutation (replace each ρ(gₖ) with ρ(gₖ⁻¹)): changes the matrix entries,
  since 2I is non-abelian → ∂∂ ≠ 0 → GATE FAILS.

Each mutation is applied SEPARATELY; only one convention changes at a time.

### 4.3 Derivation-path gates (§7)

| ID | Gate | What it establishes |
| --- | --- | --- |
| D1 | For each irrep: twisted ∂∂ = 0 | Twisted complex is a complex |
| D2 | For each irrep: rank verification for D₁, D₂, D₃ | Acyclicity verified |
| D3 | For each irrep: the determinant factors Δ₁, Δ₂, Δ₃ are recorded and T² is computed from them | Torsion comes from the derivation intermediates |
| D4 | Galois consistency: for each Galois pair, applying σ: φ↦1−φ to T²(ρ) gives T²(σ(ρ)) | Galois action is consistent |
| D5 | Independence of minor choice: torsion computed with two different choices of J₁, J₂ agrees | Well-definedness of torsion |

### 4.4 Mutation tests (§9)

Every gate above is accompanied by a mutation that reddens it:

| Gate | Mutation | Expected result |
| --- | --- | --- |
| M1 (∂∂=0) | Perturb one entry of d₂ by adding an extra term | ∂₃∂₂ ≠ 0 |
| M2 (ranks) | Change free_ranks[3] from 1 to 2 | χ ≠ 0 |
| M3 (augmented homology) | Replace ε with a non-augmentation map | H₀ changes |
| M4 (integral saturation) | Scale d₃ by a non-unit group ring element | Saturation fails at im ∂₃ |
| M5 (augmentation) | Replace ε(g) = 1 with ε(g) = 0 for g ≠ e | ε is not an augmentation |
| M6 (generator corr.) | Swap s and t IDs | d₁ doesn't match |
| M7 (acyclicity) | Use trivial rep (expected non-acyclic) | Non-acyclic → correctly classified |
| T1 (unitarity) | Use a non-unitarizable "representation" | Form is not positive-definite |
| T2 (row signature) | Swap two same-dimension irreps | Characters don't match |
| T3 (fixture eval) | Apply g ↦ ρ(g⁻¹)ᵀ on Z/5 fixture | Character mismatch |
| T3 (fixture dir) | Transpose block structure of D_n | ∂∂ ≠ 0 |
| T3 (fixture vec) | Transpose entire D_n | ∂∂ ≠ 0 |
| T3 (fixture mod) | Replace ρ(g) with ρ(g⁻¹) in evaluation | ∂∂ ≠ 0 |
| D1 (twisted ∂∂) | Perturb one matrix entry | ∂∂ ≠ 0 |
| D2 (rank) | Set one boundary to zero matrix | Rank fails |
| D3 (det factors) | Replace Δ₂ with 1 | T² changes |
| D4 (Galois) | Apply wrong Galois action | Consistency fails |
| D5 (independence) | Use only one choice | Not independently verified (assertion) |

## 5. Overlap disclosure (§6)

**Shared between routes:**
- Both routes use the group 2I and its 9 irreducible representations
- Both routes output values in Q(φ)
- Both routes rely on the Cheeger-Müller theorem for the correspondence T_an = |τ_R|

**Disjoint:**
- The analytic route computes T² from twisted Laplacian spectra via ζ'(0) values
- This route computes T² from boundary map matrices via Reidemeister torsion determinants
- The analytic route uses no chain complex; this route uses no spectral data
- The analytic route processes spectral multiplicities; this route processes matrix ranks
  and determinants

**Libraries:** Python standard library (fractions, hashlib, json), numpy (for integer
matrix operations in saturation certificates), sympy (for Hermite/Smith normal form).
No spectral computation library is used.

## 6. Pre-implementation validation artifacts

| File | Purpose | Result |
| --- | --- | --- |
| `validate_manifest.py` | Verify canonical enumeration, relators, ∂∂=0, augmentation, ranks, SHA-256 hashes | ALL PASS |

### Pre-implementation validation results (from validate_manifest.py)

- Group packet SHA-256 matches protocol pin: PASS
- Construction packet SHA-256 matches protocol pin: PASS
- Group order = 120: PASS
- Canonical enumeration SHA-256 matches protocol: PASS
- Rank 0 = [-2,0,0,0,0,0,0,0] (-1): PASS
- Rank 119 = [2,0,0,0,0,0,0,0] (+1): PASS
- Rank 118 = [1,0,1,0,1,0,1,0] (s): PASS
- Identity at rank 119: PASS
- Relators s³ = t⁵ = (st)² = -1: PASS
- ⟨s,t⟩ generates full group: PASS
- d₂d₁ = 0 over Z[2I]: PASS
- d₃d₂ = 0 over Z[2I]: PASS
- ε(d₁) = 0: PASS
- Free ranks [1,2,2,1], χ = 0: PASS
- d₁ generator correspondence: PASS

## 7. Key element IDs

| Element | ID | 8-int key | Order |
| --- | --- | --- | --- |
| identity (1) | 119 | [2,0,0,0,0,0,0,0] | 1 |
| central element (-1) | 0 | [-2,0,0,0,0,0,0,0] | 2 |
| s | 118 | [1,0,1,0,1,0,1,0] | 6 |
| t | 80 | [0,1,0,0,-1,1,1,0] | 10 |
| st | 71 | [0,0,1,0,-1,1,0,1] | 4 |

MANIFEST STATUS: FINAL; pre-implementation validation complete.
