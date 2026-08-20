# Method-and-Gate Manifest (M8.8)

Written BEFORE implementation per protocol § 8 step 4.
Once implementation begins this document is IMMUTABLE.

## 1. Route Selection

**Route class**: Combinatorial Reidemeister torsion from the supplied based chain complex,
computed via determinant data (alternating product of maximal minors of evaluated boundary
matrices). This is within the § 6 permitted class.

**Specific construction**: For each nontrivial irreducible representation ρ of 2I with
dimension d, evaluate the based chain complex

    0 → C_3(ρ) → C_2(ρ) → C_1(ρ) → C_0(ρ) → 0

with dimensions d, 2d, 2d, d obtained by tensoring C_* ⊗_{Z[2I]} V_ρ. After verifying
acyclicity (all ranks = d), compute the Reidemeister torsion τ_ρ via:

    τ_ρ = ε · det(D_3[:, J₃]) · det(D_1[J₁, :]) / det(D_2[J₃^c, J₁^c])

where D_i are the evaluated boundary matrices, J₁ is a set of d row indices of D_1 forming
a nonsingular d×d submatrix, J₃ is a set of d column indices of D_3 forming a nonsingular
d×d submatrix, J₃^c and J₁^c are the complements, and ε is a sign. The target value is

    T²(ρ) = |τ_ρ|²

as the squared modulus, computed exactly in Q(φ).

## 2. Declared Native Orientation

The native orientation is the one induced by the supplied basing without modification:
T²(ρ) = |τ_ρ|² where τ_ρ is computed from the chain complex with the declared bases,
declared boundary direction, and declared top-closure orientation exactly as given in the
construction packet. No sign flip, inversion, or other adjustment is applied to the computed
values.

Whether this native convention agrees with the M8.3 analytic convention or its inverse is
determined post-reveal by the § 5.4 anchor rule at R7.

## 3. Conventions Consumed

From the construction packet (`m8_8_construction_packet.json`):

| Convention | Value consumed |
| --- | --- |
| module_side | left: C_* is a left Z[2I]-module |
| vector_convention | row: chains are row vectors |
| boundary_direction | right action: c ∈ C_n maps to c · ∂_n ∈ C_{n-1} |
| evaluation | g ↦ ρ(g): direct evaluation, no inverse/transpose/dual |
| augmentation | ε: C_0 → Z sends every group element to 1 |
| basis_order | degree 1: (e_s, e_t); degree 2: (f_1, f_2) matching relators |
| top_closure | single 3-cell, positively oriented, fundamental class generator |
| abstract_generators | s = element 118, t = element 80 in canonical enumeration |

From the group packet (`m8_5a_packet.json`):

| Convention | Value consumed |
| --- | --- |
| coefficient_field | Q(φ), φ² = φ + 1 |
| coefficient_form | (a + b·φ)/2 with fixed denominator 2 |
| quaternion_basis | 1, i, j, k |

## 4. Representation Construction

All 9 irreducible representations of 2I are constructed from the natural 2-dimensional
representation ρ₁ via symmetric powers:

| Label | Construction | Dimension |
| --- | --- | --- |
| R0 | trivial | 1 |
| R1 | ρ₁ (natural embedding 2I ⊂ SU(2) ⊂ GL₂(ℂ)) | 2 |
| R2 | σ(ρ₁) (Galois conjugate: φ ↦ 1−φ in matrix entries) | 2 |
| R3 | Sym²(ρ₁) | 3 |
| R4 | Sym²(ρ₂) = σ(Sym²(ρ₁)) | 3 |
| R5 | Sym³(ρ₁) | 4 |
| R6 | Sym³(ρ₂) = σ(Sym³(ρ₁)) | 4 |
| R7 | Sym⁴(ρ₁) | 5 |
| R8 | Sym⁵(ρ₁) | 6 |

The natural representation embeds a unit quaternion q = a+bi+cj+dk as:

    ρ₁(q) = [[a+bi, −c−di], [c−di, a−bi]]

Symmetric powers use the monomial basis {z₁^{k-j} z₂^j}_{j=0}^{k} for Sym^k.
Galois conjugates apply φ ↦ 1−φ to the Q(φ)-coefficients of quaternion coordinates
before computing the matrix.

Identification of representations in the output uses the § 5.5 public row signature:
dimension plus exact Q(φ) characters on s, t, and st.

## 5. Instantiated Gates

### 5.1 Model Gates (on the supplied complex, pre-evaluation)

| ID | Gate | Mutation to redden |
| --- | --- | --- |
| M1 | ∂_n ∂_{n+1} = 0 over Z[2I] for all n | Flip sign of one term in ∂_2 |
| M2 | Free ranks = [1, 2, 2, 1] and χ = 1−2+2−1 = 0 | Alter a declared rank |
| M3 | H_*(Z ⊗_{Z[2I]} C_*) ≅ (Z, 0, 0, Z) via declared augmentation | Replace ε by a non-augmentation map |
| M4 | H_*(C_*) ≅ (Z, 0, 0, Z) as Z-modules with integral saturation certificate for im ∂_1, im ∂_2, im ∂_3 (maximal minor det = ±1) | ∂_3 → 2·∂_3 (changes det of minor, non-saturated image) |
| M5 | Terminal map C_0 → Z is the declared augmentation ε | Replace ε(g) = 1 ∀g by a map sending one generator to 0 |
| M6 | ∂_1 matches the abstract generator correspondence s=118, t=80; and ε(∂_1) = 0 | Swap s and t IDs |
| M7 | Per-irrep acyclicity with expected outcomes: R0 non-acyclic is PASS; every nontrivial irrep acyclic is PASS; a nontrivial irrep non-acyclic is hypothesis failure | Test with R0: correctly reports non-acyclic |

### 5.2 Theorem-Side Gates

| ID | Gate | Mutation to redden |
| --- | --- | --- |
| T1 | Invariant positive-definite Hermitian form by group averaging on each consumed representation; verifies ⟨ρ(g)v, ρ(g)w⟩ = ⟨v, w⟩ | Skip averaging, use a non-invariant form |
| T2 | Row identity via § 5.5 public signature: (dim, χ(s), χ(t), χ(st)) separates all 9 irreps | Swap matrices of two same-dimension irreps |
| T3 | Convention fixture: exact non-unitary 2-dim representation processed through same code path; GREEN under declared conventions; each of the four conventions (evaluation, module side, vector convention, boundary direction) mutated separately, each reddening at least one gate | Each convention mutation is its own sub-test |

### 5.3 Derivation-Path Gates

| ID | Gate | Mutation to redden |
| --- | --- | --- |
| D1 | ∂∂ = 0 per evaluated representation (twisted complex is a complex) | Perturb one evaluated boundary entry |
| D2 | Rank check per evaluated representation (all ranks = dim ρ for acyclic) | Would fail if acyclicity fails |
| D3 | Determinant factors (Δ₁, Δ₂, Δ₃) are nonzero and trace to final T² | Set one factor to 1 instead of computed value |
| D4 | Galois consistency: T²(σ(ρ)) = σ(T²(ρ)) for each Galois pair | Would fail if representation construction or evaluation has a Galois-breaking bug |

## 6. Overlap Disclosure

### Shared with the analytic route:
- The mathematical identity invoked: Cheeger-Müller theorem (the theoretical bridge between routes).
- Standard representation theory of finite groups (character theory, irrep construction via symmetric powers).
- Q(φ) exact arithmetic.
- The supplied based chain complex (the construction packet is shared input).
- The group packet (shared input).

### Not shared with the analytic route:
- No spectral data, zeta functions, or heat-kernel computations.
- No twisted spectra or spectral multiplicities.
- No quasi-quadratic fit or Hurwitz-zeta reduction.
- No closed-form torsion fixture.
- The computational method (alternating product of determinants of boundary-map minors vs. spectral-zeta regularization) is disjoint.
- No code or implementation artifact from M8.3 or the mode-identity-theory artifact.

## 7. Derivation Artifacts (§ 7)

The following intermediate objects are produced BY THE SAME CODE PATH that computes the
final T² values and are committed with the raw output:

1. The canonical element enumeration (120 elements with IDs, SHA-256 verified).
2. Per-irrep evaluated boundary matrices D_1(ρ), D_2(ρ), D_3(ρ).
3. Per-irrep ∂∂ = 0 verification results.
4. Per-irrep rank/acyclicity verification results.
5. Per-irrep determinant factors Δ₁, Δ₂, Δ₃ and the torsion τ_ρ.
6. Per-irrep T² = |τ_ρ|² before normalization to Q(φ) triple.
7. The integral saturation certificates (maximal minor determinants).
8. The convention fixture results (pass under correct conventions, fail under each mutation).

Each artifact is hashed and referenced in RAW_OUTPUT.json.
