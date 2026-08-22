# M8.8 Method-and-Gate Manifest

Written before production implementation, per protocol § 8 step 4.

## 1. Route selection

**Route class:** Combinatorial Reidemeister torsion, from the supplied finite based
3-dimensional chain complex, by determinant data.

**Construction:** For each nontrivial irreducible representation ρ of 2I (dimension d),
evaluate the based chain complex C_* ⊗_{Z[2I]} V_ρ by replacing each group-ring element
g with ρ(g). This produces an acyclic complex of C-vector spaces with dimensions
(d, 2d, 2d, d). The torsion τ_ρ is computed as an alternating product of determinants of
square sub-matrices of the evaluated boundary maps. The output is T²(ρ) = |τ_ρ|², the
squared modulus, yielding an element of Q(φ).

**Method disjointness:** This route computes from boundary-map matrices and their
determinants. It does not use spectral multiplicities, twisted spectra, zeta functions,
heat-kernel data, or any analytic torsion fixture.

## 2. Construction registry

| ID | Name | Description |
| --- | --- | --- |
| C-ENUM | Element enumeration | Generate all 120 elements of 2I from group-packet generators; assign canonical IDs by 8-integer sort key; verify SHA-256 |
| C-MULT | Multiplication table | Build the 120×120 multiplication table from the enumerated elements |
| C-SU2 | SU(2) embedding | Map each quaternion q = a+bi+cj+dk to the 2×2 complex matrix [[a+b𝕚, c+d𝕚],[−c+d𝕚, a−b𝕚]] |
| C-SYM | Symmetric powers | Compute Sym^n(ρ_std) for n = 0,...,6 to obtain irreps of dimensions 1,2,3,4,5,6 plus Sym^6 (dim 7) |
| C-GAL | Galois conjugates | Apply σ: φ → 1−φ to matrix entries of Sym^1, Sym^2 for two Galois-conjugate irreps (dims 2,3); Sym^3 is self-conjugate |
| C-PROJ | V₆ extraction | Project V₆ (dim 4) from Sym^6 = V₆ ⊕ V₈ using character projection P = (4/120)Σ_g χ_{V₆}(g)·ρ_{Sym⁶}(g); extract basis from image of P |
| C-CHAR | Character computation | Compute tr(ρ(g)) for each irrep and each group element |
| C-ROWSIG | Row signature | For each irrep, compute (dim, χ(s), χ(t), χ(st)) as the label-free identifier per § 5.5 |
| C-GREVAL | Group-ring evaluation | Replace each group-ring entry Σ c_α g_α with Σ c_α ρ(g_α) to produce d×d matrix blocks |
| C-BDRY | Boundary evaluation | Assemble the full evaluated boundary matrices M_1 (2d×d), M_2 (2d×2d), M_3 (d×2d) in the row-vector right-action convention |
| C-ACYC | Acyclicity check | Verify rank conditions: rank(M_3)=d, rank(M_2)=d, rank(M_1)=d for each nontrivial irrep |
| C-TORSION | Torsion formula | Choose column indices J₃ of M₃ and row indices I₁ of M₁ giving nonsingular d×d minors; the complementary minor M₂[J₃ᶜ,I₁ᶜ] gives τ_ρ = ε · det(M₂[J₃ᶜ,I₁ᶜ]) / (det(M₁[I₁,:]) · det(M₃[:,J₃])), where ε = ±1 cancels in |τ|² |
| C-MODSQ | Squared modulus | Compute T²(ρ) = |τ_ρ|² = τ_ρ · conj(τ_ρ), yielding an element of Q(φ) |
| C-UNITARY | Unitarity verification | Construct the invariant positive-definite Hermitian form by group averaging: H = (1/120) Σ_g ρ(g)* ρ(g); verify H is positive definite and ρ(g)* H ρ(g) = H for all g |
| C-FIXTURE | Convention fixture | A synthetic non-unitary representation and chain-complex instance exercising all four declared conventions; processed through the same parser and evaluation path as the target run |

## 3. Convention registry

| ID | Name | Source | Value |
| --- | --- | --- | --- |
| V-MODSIDE | Module side | packet basing | left: 2I acts on the left of V_ρ |
| V-ROWVEC | Vector convention | packet basing | row: chains are row vectors |
| V-EVAL | Evaluation map | packet basing | g ↦ ρ(g); no inverse, transpose, or dual variant |
| V-BDRYDIR | Boundary direction | packet basing | right: boundary maps act on the right, c ↦ c · d_k |
| V-AUG | Augmentation | packet basing | ε: C_0 → Z sends every group element to 1; ε is the terminal map and is NOT d_1 |
| V-ORIENT | Native orientation | declared here | The computation produces T²_native(ρ) from the supplied based complex with no sign choice. The native orientation is the one determined by the packet's basing and top-closure convention |
| V-QPHI | Q(φ) output encoding | protocol § 5.5 | Normalized triple (a, b, c) for (a + b·φ)/c with c > 0 and gcd(a, b, c) = 1 |

## 4. Pre-reveal gate registry

### 4.1 Model gates

| Gate ID | Gate name | Mutation | What it establishes |
| --- | --- | --- | --- |
| G-M01 | ∂₂∂₃ = 0 | Perturb one entry of d₃; verify product becomes nonzero | Chain-complex condition at degree 2 |
| G-M02 | ∂₁∂₂ = 0 | Perturb one entry of d₂; verify product becomes nonzero | Chain-complex condition at degree 1 |
| G-M03 | Free ranks and χ | Alter one rank; verify χ ≠ 0 | Closed 3-manifold complex, not a presentation 2-complex |
| G-M04 | Augmented homology | Replace one entry of d₂; verify H_* changes | H_*(Z ⊗_{Z[2I]} C_*) ≅ (Z, 0, 0, Z) |
| G-M05 | Universal-cover homology | Scale one row of d₂ by non-unit; verify saturation fails | H_*(C_*) ≅ (Z, 0, 0, Z) integrally with saturation certificates |
| G-M06 | Augmentation is terminal map | Replace ε with a non-augmentation; verify ε∂₁ ≠ 0 | Terminal map is the declared augmentation |
| G-M07 | Generator correspondence | Swap s and t IDs; verify relator check fails | ∂₁ matches abstract_generators; ε∂₁ = 0 |
| G-M08 | Per-irrep acyclicity | For a nontrivial irrep, confirm rank deficiency in mutated complex makes it non-acyclic | R0 non-acyclic (PASS); all nontrivial irreps acyclic (PASS) |

### 4.2 Theorem-side gates

| Gate ID | Gate name | Mutation | What it establishes |
| --- | --- | --- | --- |
| G-T01 | Unitarity | Perturb one representation matrix; verify Hermitian form invariance fails | Unitary hypothesis of Cheeger-Müller |
| G-T02 | Row signature identity | Swap two character values; verify signature match fails | Same flat bundle per label as analytic route |
| G-T03a | Convention: evaluation map | On fixture: replace g ↦ ρ(g) with g ↦ ρ(g⁻¹) (anti-homomorphism); verify ∂∂≠0 or T² changes. All 2I irreps are self-contragredient, so g↦ρ(g⁻¹)ᵀ gives the same T²; the anti-homomorphism mutation is the correct causal test | Evaluation convention exercised causally |
| G-T03b | Convention: boundary direction | On fixture: transpose boundary matrices (left action instead of right); verify gate failure | Boundary direction exercised causally |
| G-T03c | Convention: module side | On fixture: use right-module action instead of left; verify gate failure | Module side exercised causally |
| G-T03d | Convention: vector convention | On fixture: use column vectors instead of row; verify gate failure | Vector convention exercised causally |

### 4.3 Derivation-path gates

| Gate ID | Gate name | Mutation | What it establishes |
| --- | --- | --- | --- |
| G-D01 | ∂∂ = 0 per twisted complex | For one irrep, perturb evaluated boundary; verify M₃·M₂ ≠ 0 | Chain condition on the consumed twisted complex |
| G-D02 | Twisted complex ranks | For one irrep, zero out a row of M₃; verify rank drops | Correct ranks of consumed boundary maps |
| G-D03 | Determinant sub-matrices nonsingular | For one irrep, zero a column of the chosen minor; verify det = 0 | Determinant factors consumed by torsion product |
| G-D04 | Galois consistency | Swap φ ↦ 1−φ on characters of one Galois pair; verify the two T² values are no longer Galois conjugates | Galois action on character field |
| G-D05 | Torsion code-path dependency | Replace torsion formula input with identity matrices; verify output changes | Derivation artifacts are consumed (not bypassed) by the result |

## 5. Native orientation declaration

The native orientation is determined by the construction packet's basing:
- The single 3-cell carries the positive fundamental class (top_closure).
- Boundary maps act on the right in the row-vector convention.
- The evaluation map is g ↦ ρ(g) with no inversion.

T²_native(ρ) := |τ_ρ|² computed in this convention.

The global involution T² ↔ (T²)⁻¹ is the sole admitted bridge (§ 5.4). The
adjudicator selects between the native table and its global inverse at R7 after
the answer packet opens. This manifest does not and cannot know which orientation
matches the analytic convention.

## 6. Overlap disclosure

**Shared with the analytic route:**
- The group 2I and its representation theory (character table, irrep dimensions).
- The golden ratio φ = (1+√5)/2 and Q(φ) arithmetic.
- The concept of torsion as a topological invariant.
- Standard linear algebra (determinants, matrix arithmetic).

**Not shared:**
- The analytic route uses spectral-zeta functions, twisted spectra, and
  quasi-quadratic fits. This route uses boundary-map matrices and their determinants.
- The analytic route works from the spectrum of the Laplacian on twisted forms.
  This route works from the based chain complex.
- No spectral data, heat-kernel data, or zeta values are used here.

## 7. Pre-implementation coverage table

| Registry ID | Checkable? | Validation artifact | Check | Status |
| --- | --- | --- | --- | --- |
| C-ENUM | Yes | validate_enumeration.py | SHA-256 of 120-element enumeration matches protocol § 4.2 hash | PASS |
| C-MULT | Yes | validate_enumeration.py | Identity element acts as identity in table; relators s³=(st)², t⁵=(st)² hold | PASS |
| C-SU2 | Yes | validate_representations.py | SU(2) matrices have det=1, are unitary, and reproduce correct quaternion multiplication | PASS |
| C-SYM | Yes | validate_representations.py | Sym^n matrices satisfy group homomorphism property; characters match Weyl formula | PASS |
| C-GAL | Yes | validate_representations.py | Galois-conjugate characters are related by φ↔1−φ; inner products confirm irreducibility; Sym³ verified self-conjugate | PASS |
| C-PROJ | Yes | validate_representations.py | V₆ projected from Sym⁶; P idempotent, rank 4; homomorphism, character match confirmed | PASS |
| C-CHAR | Yes | validate_representations.py | Character orthogonality relations hold; 9 distinct irreps found with correct dimensions | PASS |
| C-ROWSIG | Yes | validate_representations.py | Row signatures are distinct; dimensions match 1,2,2,3,3,4,4,5,6 | PASS |
| C-GREVAL | Yes | validate_representations.py | Evaluation preserves ring structure: ρ(g₁·g₂) = ρ(g₁)·ρ(g₂) | PASS |
| C-BDRY | Yes | validate_torsion_dry.py | ∂∂=0 for each twisted complex; dimensions match expected | PASS |
| C-ACYC | Yes | validate_torsion_dry.py | Rank conditions verified for each nontrivial irrep; R0 non-acyclic confirmed | PASS |
| C-TORSION | Not before production | — | The torsion formula is algebraically correct by construction; its output is the target result. Cannot be checked target-blind since the result IS the target | — |
| C-MODSQ | Not before production | — | |τ|² = τ·τ̄ is correct by definition; the Q(φ) membership is verified by the Galois-consistency gate. Cannot be pre-checked without computing the actual torsion | — |
| C-UNITARY | Yes | validate_representations.py | Group-averaged form H = (1/120)Σρ(g)*ρ(g) is positive definite and invariant | PASS |
| C-FIXTURE | Yes | validate_fixture.py | Fixture passes all applicable gates under declared conventions; each single-convention mutation reddens at least one gate | PASS |
| V-MODSIDE | Yes | validate_fixture.py | Convention fixture exercises left module action; right-module mutation reddens a gate | PASS |
| V-ROWVEC | Yes | validate_fixture.py | Convention fixture uses row vectors; column-vector mutation reddens a gate | PASS |
| V-EVAL | Yes | validate_fixture.py | Convention fixture uses g↦ρ(g); contragredient mutation reddens a gate on non-unitary fixture | PASS |
| V-BDRYDIR | Yes | validate_fixture.py | Convention fixture uses right action; transposed (left) action mutation reddens a gate | PASS |
| V-AUG | Yes | validate_complex.py | ε∂₁=0 verified; mutation (non-augmentation map) breaks it | PASS |
| V-ORIENT | No | — | Orientation is a declaration; it is not checkable before the answer packet opens | — |
| V-QPHI | Yes | validate_representations.py | Output triples satisfy gcd(a,b,c)=1 and c>0 | PASS |

## 8. Validation artifacts

| File | Purpose |
| --- | --- |
| validate_enumeration.py | Enumerates 2I, verifies SHA-256, relators, generator orders, generation |
| validate_complex.py | Verifies ∂∂=0 over Z[2I], free ranks, χ=0, augmented homology, generator correspondence |
| validate_saturation.py | Computes unimodular saturation certificates for all three boundary maps; establishes H_*(C_*)=(Z,0,0,Z) integrally |
| validate_representations.py | Constructs all 9 irreps, verifies character orthogonality, unitarity, row signatures |
| validate_torsion_dry.py | Evaluates boundary maps for each irrep, verifies ∂∂=0 and acyclicity of twisted complexes |
| validate_fixture.py | Convention fixture with per-convention mutation tests |
| validate_manifest.py | Automated registry-coverage set-equality checker |

## 9. Pre-implementation validation results

### Enumeration (validate_enumeration.py)
- 120 elements generated ✓
- SHA-256 = 27ff780d28d5d854d464ead87e8fc20244fac8334bda9f0600c6ee1b3c89561e ✓
- Rank 0 = [-2,0,0,0,0,0,0,0], Rank 118 = [1,0,1,0,1,0,1,0], Rank 119 = [2,0,0,0,0,0,0,0] ✓
- Relators s³ = (st)² and t⁵ = (st)² ✓
- order(s) = 6, order(t) = 10, order(st) = 4 ✓
- ⟨s,t⟩ generates all 120 elements ✓

### Chain complex (validate_complex.py)
- d₃·d₂ = 0 over Z[2I] ✓
- d₂·d₁ = 0 over Z[2I] ✓
- Free ranks [1,2,2,1], χ = 0 ✓
- ε·d₁ = 0 ✓
- d₁ matches generator correspondence (s→118, t→80) ✓
- Augmented homology: det(d₂_aug) = −1 (unimodular) ✓
- H_*(Z ⊗ C_*) = (Z, 0, 0, Z) ✓

### Saturation (validate_saturation.py)
- rank(d₁_Z) = 119, rank(d₂_Z) = 121, rank(d₃_Z) = 119 ✓
- Unimodular (|det|=1) 119×119 minor of d₃_Z found ✓
- Unimodular 121×121 minor of d₂_Z found ✓
- Unimodular 119×119 minor of d₁_Z found ✓
- All boundary images saturated ⇒ H_*(C_*) = (Z, 0, 0, Z) integrally ✓

### Representations (validate_representations.py)
- SU(2) embedding: homomorphism verified ✓
- Sym^0 through Sym^6: homomorphism verified ✓
- σ(Sym^1), σ(Sym^2): homomorphism verified ✓
- Sym^3 self-conjugate: ⟨Sym³,σ(Sym³)⟩ = 1 ✓
- V₆ character norm = 1, orthogonal to V₀–V₅, V₇, V₈ ✓
- V₆ extracted from Sym⁶ via character projection: P idempotent, rank 4 ✓
- V₆ homomorphism and character match ✓
- Full 9×9 character orthogonality ✓
- Unitarity: invariant Hermitian form H for all 9 irreps ✓
- All 9 row signatures distinct ✓
- Galois pairs: V₁↔V₇, V₂↔V₈; V₀,V₃,V₄,V₅,V₆ self-conjugate ✓
- Dimensions: [1,2,2,3,3,4,4,5,6], sum of squares = 120 ✓

### Twisted complexes (validate_torsion_dry.py)
- V0 (dim 1): ∂∂=0 ✓, ranks (0,2,0), NON-ACYCLIC (trivial rep) ✓
- V1 (dim 2): ∂∂=0 ✓, ranks (2,2,2) ACYCLIC ✓
- V2 (dim 3): ∂∂=0 ✓, ranks (3,3,3) ACYCLIC ✓
- V3 (dim 4): ∂∂=0 ✓, ranks (4,4,4) ACYCLIC ✓
- V4 (dim 5): ∂∂=0 ✓, ranks (5,5,5) ACYCLIC ✓
- V5 (dim 6): ∂∂=0 ✓, ranks (6,6,6) ACYCLIC ✓
- V6 (dim 4): ∂∂=0 ✓, ranks (4,4,4) ACYCLIC ✓
- V7 (dim 2): ∂∂=0 ✓, ranks (2,2,2) ACYCLIC ✓
- V8 (dim 3): ∂∂=0 ✓, ranks (3,3,3) ACYCLIC ✓

### Convention fixture (validate_fixture.py)
- Fixture: non-unitary 2D representation (Sym¹ conjugated by non-unitary P) ✓
- Baseline T² computed and verified real ✓
- G-T03a: evaluation map mutation (g→ρ(g⁻¹)) → dd≠0 REDDENED ✓
- G-T03b: boundary direction mutation (cochain complex) → dd≠0 REDDENED ✓
- G-T03c: module side mutation (transpose = anti-hom) → dd≠0 REDDENED ✓
- G-T03d: vector convention mutation (GR transpose) → dd≠0 REDDENED ✓

### Registry-coverage validator (validate_manifest.py)
- Construction + Convention IDs == Coverage table IDs ✓
- 19 gate IDs registered ✓

---

**MANIFEST STATUS: FINAL; pre-implementation validation complete.**
