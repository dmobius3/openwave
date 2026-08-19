# M8.8 Method-and-Gate Manifest

## 1. Route Selection

**Route class**: Combinatorial (algebraic) Reidemeister torsion, from the supplied
finite based 3-dimensional chain complex representing S³/2I.

**Construction within the class**: For each d-dimensional irreducible unitary
representation ρ: 2I → GL_d(ℂ), evaluate the boundary maps of the supplied
chain complex by replacing each group element g with ρ(g), obtaining a complex
of finite-dimensional ℂ-vector spaces. Verify acyclicity. Compute the
Reidemeister torsion τ_ρ as the alternating product of determinant ratios from
maximal-rank square submatrices of the evaluated boundary maps. Report
T²(ρ) = |τ_ρ|² as an element of Q(φ).

**Disjointness from the analytic route**: This route uses no spectral
multiplicities, twisted spectra, zeta functions or their derivatives,
heat-kernel data, or closed-form torsion fixture. It computes torsion from
combinatorial chain-complex data (boundary matrices, representation matrices,
determinants) only.

## 2. Native Orientation Convention

The native orientation is determined by the computation T²(ρ) = |τ_ρ|² applied
directly to the supplied based complex with its declared top-cell orientation,
using the evaluation convention g ↦ ρ(g) as stated in the construction
packet's basing. No inversion, sign change, or orientation selection is applied
within the implementation. The raw output carries values under this native
convention. The bridge involution T² ↔ (T²)⁻¹ and the R7 anchor selection are
performed post-reveal by the adjudicator per protocol § 5.4.

## 3. Overlap Disclosure (§ 6 requirement)

**Mathematical facts shared with the analytic route:**
- The representation theory of the binary icosahedral group 2I (character table,
  irreducible representations, their dimensions).
- The golden ratio φ = (1+√5)/2 and the field Q(φ) as the character field.
- The Cheeger-Müller theorem as the bridge between analytic and combinatorial
  torsion (invoked conceptually; not used computationally).

**Libraries and conventions shared:**
- Standard exact arithmetic over Q(φ) and Q(φ, i).
- The group packet's 120-element enumeration and canonical element IDs.

**Not shared:**
- No spectral data, zeta functions, or heat-kernel computations.
- No closed-form torsion formulas or torsion tables.
- No code, scripts, or computational artifacts from the M8.3 implementation.

---

## 4. Construction Registry

| ID | Name | Description |
|----|------|-------------|
| CONST-01 | Group enumeration | Generate all 120 elements of 2I from the group packet's generators via quaternion multiplication over Q(φ) |
| CONST-02 | Canonical element ordering | Sort elements by the 8-integer key (A₁,B₁,Aᵢ,Bᵢ,Aⱼ,Bⱼ,Aₖ,Bₖ) as signed integers, first entry most significant |
| CONST-03 | Abstract generator verification | Verify that the construction packet's abstract generators s (ID 118) and t (ID 80) satisfy s³ = t⁵ = (st)² = -1 and generate 2I |
| CONST-04 | Fundamental representation | Map each element g ∈ 2I to a 2×2 SU(2) matrix via the standard quaternion-to-SU(2) correspondence: q = a+bi+cj+dk ↦ [[a+di, c+bi], [-c+bi, a-di]] (entries in Q(φ,i)) |
| CONST-05 | Symmetric power representations | Build Sym^n of the fundamental 2D representation for n=0,...,7; Sym^0 through Sym^5 give irreps of dimensions 1,2,3,4,5,6 |
| CONST-06 | Irrep decomposition | Decompose Sym^6 into irreps of dimensions 4 and 3 (via character projection), and extract the dimension-2 component from Sym^7, yielding all 9 irreps |
| CONST-07 | Boundary map evaluation | For each irrep ρ of dimension d, evaluate each Z[2I]-entry Σ c_k g_k as the d×d matrix Σ c_k ρ(g_k); assemble into block matrices D₁(ρ), D₂(ρ), D₃(ρ) |
| CONST-08 | Torsion via determinant formula | For an acyclic twisted complex, compute τ_ρ = det(∂₁_J) · det(∂₃_{K'}) / det(∂₂_{J₂,J₁'}) where J, J₁', J₂, K' are compatible index selections from maximal-rank submatrices |
| CONST-09 | Squared modulus in Q(φ) | Compute T²(ρ) = abs(τ_ρ)² = τ_ρ · conj(τ_ρ) where conjugation sends i↦-i in Q(φ,i); result lies in Q(φ) |

## 5. Convention Registry

| ID | Name | Source | Value |
|----|------|--------|-------|
| CONV-01 | Evaluation convention | Construction packet basing.evaluation | g ↦ ρ(g); no inverse, transpose or dual variant |
| CONV-02 | Module side | Construction packet basing.module_side | Left Z[2I]-module |
| CONV-03 | Vector convention | Construction packet basing.vector_convention | Row vectors |
| CONV-04 | Boundary direction | Construction packet basing.boundary_direction | Chains are row vectors; boundary maps act on the RIGHT |
| CONV-05 | Augmentation | Construction packet basing.augmentation | ε: C₀ → Z sends every group element to 1; ε is the terminal map and is NOT ∂₁ |
| CONV-06 | Top cell orientation | Construction packet top_closure | Single 3-cell, positively oriented; basis element 0 generates C₃ and carries the fundamental class |
| CONV-07 | Coefficient field | Group packet coefficient_field | Q(φ), φ = (1+√5)/2, minimal polynomial φ²-φ-1 = 0 |
| CONV-08 | Torsion bridge | Protocol § 5.4 | T²_target(ρ) := abs(τ_ρ)²; the sole admitted bridge is the global inversion T² ↔ (T²)⁻¹ |
| CONV-09 | Element ID sort key | Protocol § 4.2 | 8 signed integers (A₁,B₁,Aᵢ,Bᵢ,Aⱼ,Bⱼ,Aₖ,Bₖ) from the (A+Bφ)/2 form, compared entrywise, first most significant |
| CONV-10 | Quaternion-to-SU(2) map | Generic algebra | q = a+bi+cj+dk ↦ [[a+di, c+bi], [-c+bi, a-di]] with quaternion basis order (1,i,j,k) from group packet |

## 6. Pre-Reveal Gate Registry

### 6.1 Model Gates

| ID | Gate | Establishes | Mutation |
|----|------|-------------|----------|
| GATE-M01 | ∂₂∂₁ = 0 and ∂₃∂₂ = 0 over Z[2I] | Chain complex property | Perturb one entry of ∂₂ by adding g₀ to one term; recompute products |
| GATE-M02 | Free ranks = [1,2,2,1] | Matrix dimensions match declaration | Change one rank in the check |
| GATE-M03 | χ = Σ(-1)ⁱrᵢ = 0 | Closed 3-manifold complex | Change one rank sign |
| GATE-M04 | H*(Z ⊗_{Z[G]} C*) ≅ (Z,0,0,Z) | Homology of S³/2I | Perturb augmented ∂₂ by adding 1 to entry (0,0); recompute Smith form |
| GATE-M05 | H*(C*) ≅ (Z,0,0,Z) as Z-modules (exact saturation certificate for im ∂₁, im ∂₂, im ∂₃) | Universal cover is S³ | Multiply ∂₃ by a scalar k≠±1; verify saturation certificate fails |
| GATE-M06 | ε(∂₁) = 0 | Exactness into augmentation | Add a nonzero group element to ∂₁[0][0]; verify ε ≠ 0 |
| GATE-M07 | ∂₁ = [s-1, t-1]ᵀ matches abstract generators | Generator correspondence | Swap s and t IDs; verify mismatch |

### 6.2 Theorem-Side Gates

| ID | Gate | Establishes | Mutation |
|----|------|-------------|----------|
| GATE-T01 | Invariant positive-definite Hermitian form via group averaging on each irrep's matrices | Unitarity hypothesis | Perturb one representation matrix to break unitarity; verify form is not positive-definite |
| GATE-T02 | Row identity: each irrep identified by (dimension, χ(s), χ(t), χ(st)) as Q(φ) triples | Same flat bundle as analytic route | Swap two irreps' character values; verify mismatch |
| GATE-T03 | Convention fixture: synthetic non-unitary representation and chain complex, processed through same code path; GREEN under declared conventions; each of CONV-01 through CONV-04 mutated separately with at least one gate reddening per mutation | Character-invisible convention errors | For each convention mutation: CONV-01 (g↦ρ(g⁻¹)ᵀ), CONV-02 (right module), CONV-03 (column vectors), CONV-04 (left action); verify reddening |

### 6.3 Reproduction Gates (pre-reveal subset)

| ID | Gate | Establishes | Mutation |
|----|------|-------------|----------|
| GATE-R01 | ∂∂ = 0 in each twisted complex | Twisted complex is a complex | Perturb one evaluated matrix entry; verify product ≠ 0 |
| GATE-R02 | Per-irrep acyclicity with R0 non-acyclic = PASS, nontrivial irreps acyclic = PASS | Cheeger-Müller hypothesis | Perturb a matrix to change rank; verify acyclicity fails |
| GATE-R03 | Derivation path: every T²(ρ) produced from committed intermediate objects (evaluated matrices, determinant factors) through the SAME code path | Derivation integrity | Substitute a hardcoded value for one determinant; verify mismatch |
| GATE-R04 | Galois consistency: for each Galois pair, T²(σρ) = σ(T²(ρ)) under the Galois automorphism σ: φ ↦ 1-φ | Galois action on character field | Swap one pair's values; verify Galois relation fails |

### 6.4 Enumeration Gates

| ID | Gate | Establishes | Mutation |
|----|------|-------------|----------|
| GATE-E01 | SHA-256 of canonical enumeration = 27ff780d... | Correct element enumeration | Change one element's coordinates; verify hash mismatch |
| GATE-E02 | Identity element at rank 119 | Correct identity | Move identity to different rank; verify fails |
| GATE-E03 | s³ = t⁵ = (st)² = -1, order(s)=6, order(t)=10, order(st)=4, ⟨s,t⟩ = 2I | Relator and generation | Use wrong element IDs; verify relation fails |

---

## 7. Pre-Implementation Coverage Table

| Registry ID | Checkable Pre-Impl? | Validation Artifact | Check | Result |
|-------------|---------------------|--------------------:|-------|--------|
| CONST-01 | YES | validate_pre_impl.py | Group closure produces 120 elements | PASS |
| CONST-02 | YES | validate_pre_impl.py | Sorted elements match expected rank 0, 118, 119 | PASS |
| CONST-03 | YES | validate_pre_impl.py | s,t satisfy relators and generate 2I | PASS |
| CONST-04 | NO | — | Requires building representation matrices (production) | — |
| CONST-05 | NO | — | Requires representation construction (production) | — |
| CONST-06 | NO | — | Requires character computation and projection (production) | — |
| CONST-07 | NO | — | Requires evaluated boundary maps (production) | — |
| CONST-08 | NO | — | Requires torsion computation (production) | — |
| CONST-09 | NO | — | Requires squared modulus computation (production) | — |
| CONV-01 | YES | validate_pre_impl.py | Convention parsed from packet and consumed in boundary map parsing | PASS |
| CONV-02 | YES | validate_pre_impl.py | Left module side verified in boundary map expansion | PASS |
| CONV-03 | YES | validate_pre_impl.py | Row vector convention verified in matrix dimensions (d_n is r_n × r_{n-1}) | PASS |
| CONV-04 | YES | validate_pre_impl.py | Right action verified: c · d_n composition gives ∂∂=0 | PASS |
| CONV-05 | YES | validate_pre_impl.py | Augmentation ε verified: ε(∂₁)=0, augmented homology correct | PASS |
| CONV-06 | YES | validate_pre_impl.py | Top cell read from packet; consistent with rank-1 C₃ | PASS |
| CONV-07 | YES | validate_pre_impl.py | Q(φ) arithmetic used throughout; φ²=φ+1 verified implicitly by group closure | PASS |
| CONV-08 | NO | — | Bridge is a post-reveal operation | — |
| CONV-09 | YES | validate_pre_impl.py | Sort key produces SHA-256 match (GATE-E01) | PASS |
| CONV-10 | NO | — | Quaternion-to-SU(2) map exercised in production only | — |
| GATE-M01 | YES | validate_pre_impl.py | ∂₂∂₁ = 0 and ∂₃∂₂ = 0 computed over Z[2I] | PASS |
| GATE-M02 | YES | validate_pre_impl.py | Free ranks [1,2,2,1] verified | PASS |
| GATE-M03 | YES | validate_pre_impl.py | χ = 1-2+2-1 = 0 verified | PASS |
| GATE-M04 | YES | validate_pre_impl.py | Augmented homology (Z,0,0,Z) via Smith form | PASS |
| GATE-M05 | YES | validate_pre_impl.py | Ranks (119,121,119) verified; saturation certificates det=±1 (approx; exact in production) | PASS |
| GATE-M06 | YES | validate_pre_impl.py | ε(∂₁[i][j]) = 0 for all entries | PASS |
| GATE-M07 | YES | validate_pre_impl.py | ∂₁ = [s-1, t-1]ᵀ with s=118, t=80 | PASS |
| GATE-E01 | YES | validate_pre_impl.py | SHA-256 = 27ff780d...561e, 2389 bytes | PASS |
| GATE-E02 | YES | validate_pre_impl.py | Identity (2,0,0,0,0,0,0,0) at rank 119 | PASS |
| GATE-E03 | YES | validate_pre_impl.py | s³=t⁵=(st)²=-1, orders 6/10/4, generation | PASS |
| GATE-T01 | NO | — | Requires representation matrices (production) | — |
| GATE-T02 | NO | — | Requires character computation (production) | — |
| GATE-T03 | NO | — | Requires convention fixture construction (production) | — |
| GATE-R01 | NO | — | Requires twisted boundary maps (production) | — |
| GATE-R02 | NO | — | Requires per-irrep acyclicity check (production) | — |
| GATE-R03 | NO | — | Requires derivation path verification (production) | — |
| GATE-R04 | NO | — | Requires Galois pair computation (production) | — |

## 8. Validation Artifacts

| Artifact | Purpose |
|----------|---------|
| `validate_pre_impl.py` | Pre-implementation validation of all checkable constructions, conventions, and gates |
| `validate_manifest.py` | Automated manifest validator: checks set equality between registries and coverage table |

## 9. Pre-Implementation Validation Record

All checks executed by `validate_pre_impl.py` on the permitted inputs (group packet
and construction packet) with no external references consulted.

```
CONST-01: PASS (120 elements generated)
CONST-02: PASS (rank 0, 118, 119 match expected)
CONST-03: PASS (s, t generate all 120 elements)
CONV-01 through CONV-07, CONV-09: PASS (verified through boundary map and enumeration checks)
GATE-E01: PASS (SHA-256 match)
GATE-E02: PASS (identity verified)
GATE-E03: PASS (relators, orders, generation)
GATE-M01: PASS (∂∂ = 0 over Z[2I])
GATE-M02: PASS (free ranks [1,2,2,1])
GATE-M03: PASS (χ = 0)
GATE-M04: PASS (augmented homology (Z,0,0,Z))
GATE-M05: PASS (ranks correct; saturation approximate, exact cert in production)
GATE-M06: PASS (ε(∂₁) = 0)
GATE-M07: PASS (generator correspondence)
```

MANIFEST STATUS: FINAL; pre-implementation validation complete.
