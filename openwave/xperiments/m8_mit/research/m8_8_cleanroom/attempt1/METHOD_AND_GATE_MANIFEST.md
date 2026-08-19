# M8.8 Method-and-Gate Manifest

Written BEFORE implementation, per § 8 step 4.

## 1. Selected route

**Route class:** Combinatorial Reidemeister torsion (§ 6), from the supplied finite based
chain complex, by determinant data.

**Construction:** For each irreducible complex representation ρ of 2I (dimension d), form the
based complex of finite-dimensional vector spaces

```
C_*(ρ): 0 → V₃ (dim d) →^{D₃} V₂ (dim 2d) →^{D₂} V₁ (dim 2d) →^{D₁} V₀ (dim d) → 0
```

by evaluating each group-ring entry of the boundary matrices at ρ (per the declared
evaluation convention g ↦ ρ(g)). For each acyclic representation, compute the Reidemeister
torsion τ_ρ via the alternating-determinant formula: choose index sets I' ⊂ {0,...,2d−1}
and J ⊂ {0,...,2d−1}, each of size d, such that the d×d submatrices D₁[I',:], D₂[J',I],
D₃[:,J] are nonsingular (where I = complement of I', J' = complement of J). Then

```
τ_ρ = ε · det(D₂[J', I]) / (det(D₁[I', :]) · det(D₃[:, J]))
```

where ε = ±1 is the permutation sign. The target quantity is

```
T²_target(ρ) = |τ_ρ|²
```

the squared complex modulus, which is independent of ε and of the choice of I', J.

**Arithmetic:** All computation in exact Q(φ, i) arithmetic, where φ = (1+√5)/2 and
i = √(−1). The natural 2-dimensional representation of 2I ⊂ SU(2) has entries in Q(φ, i).
Higher-dimensional irreps are constructed as symmetric powers Sym^n of the natural
representation. Galois-conjugate irreps are obtained by the field automorphism φ ↦ 1−φ. All
T²_target values land in Q(φ) (since |·|² eliminates the imaginary part).

## 2. Method disjointness disclosure (§ 6)

**Shared with the analytic route:**
- The group 2I and its irreducible representations (shared subject matter)
- The coefficient field Q(φ)
- Standard linear algebra (determinants, matrix arithmetic)
- The group packet data

**Not shared:**
- No spectral multiplicities, twisted spectra, or eigenvalue data
- No zeta functions or their derivatives
- No heat-kernel data
- No closed-form torsion fixture or table
- No analytic continuation or regularization
- The chain complex itself: M8.3 never constructed or used this object

The computational paths are joined by the Cheeger–Müller theorem, not by shared machinery.

## 3. Conventions consumed (§ 4.2 basing)

All consumed from the construction packet's `basing` field:

| Convention | Value |
| --- | --- |
| Module side | left Z[2I]-module |
| Vector convention | row vectors |
| Boundary direction | chains act on the RIGHT: c ∈ C_n maps to c · d_n ∈ C_{n−1} |
| Evaluation | g ↦ ρ(g); no inverse, transpose, or dual |
| Augmentation | ε: C₀ → Z sends every group element to 1; terminal map, NOT ∂₁ |
| Element IDs | canonical enumeration from the group packet, 8-integer signed lexicographic key per § 4.2 |

## 4. Declared native orientation (§ 5.4)

**Native convention:** T²_target(ρ) = |τ_ρ|² where τ_ρ is the Reidemeister torsion of the
acyclic based complex C_*(ρ), computed in the construction packet's declared basing, with the
evaluation convention g ↦ ρ(g), and the alternating-determinant formula of § 1 above.

This is the combinatorial orientation inherent in the supplied basing. It may differ from the
analytic convention by the global involution T² ↔ (T²)⁻¹. The § 5.4 anchor rule at R7
resolves this at adjudication; no pre-reveal matching is attempted.

## 5. Instantiated gates (§ 7, § 9)

### 5.1 Model gates

| ID | Gate | Mutation |
| --- | --- | --- |
| M1 | ∂₁∂₂ = 0 and ∂₂∂₃ = 0 over Z[2I] | Perturb one entry of ∂₂; verify product becomes nonzero |
| M2 | Free ranks = [1, 2, 2, 1]; χ = 1−2+2−1 = 0 | Assert ranks and Euler characteristic |
| M3 | H_*(Z ⊗_{Z[2I]} C_*) ≅ (Z, 0, 0, Z): augmented homology via ε | Perturb augmented ∂₁; verify H₀ changes |
| M4 | H_*(C_*) ≅ (Z, 0, 0, Z) as Z-modules (universal cover), with saturation certificate for im ∂₁, im ∂₂, im ∂₃ | Multiply ∂₃ by scalar k≠±1; verify saturation fails. Multiply ∂₂ by central non-unit; verify degree-1 saturation fails |
| M5 | Terminal map C₀ → Z is the declared augmentation ε | Verify ε is not ∂₁ by checking they differ |
| M6 | ε(∂₁) = 0 and ∂₁ matches the frozen 1-cell correspondence: row i of ∂₁ is (s_i − 1) where s_i is the i-th abstract generator | Perturb generator ID; verify ε(∂₁) ≠ 0 or correspondence fails |
| M7 | Per-irrep acyclicity: R0 non-acyclic (PASS), all 8 nontrivial irreps acyclic (PASS) | For each nontrivial irrep, perturb one evaluated boundary entry; verify rank drops |

### 5.2 Theorem-side gates

| ID | Gate | Mutation |
| --- | --- | --- |
| T1 | Unitarity: for each irrep, construct invariant positive-definite Hermitian form H by group averaging H = (1/120) Σ_g ρ(g)*ᴴ ρ(g), verify H is positive definite and ρ(g)*ᴴ H ρ(g) = H for all g | Use a non-group-element matrix; verify invariance fails |
| T2 | Row identity: each irrep identified by its § 5.5 row signature (dimension, χ_ρ(s), χ_ρ(t), χ_ρ(st)) with characters computed via the group packet's generator correspondence | Swap two irrep assignments; verify signature mismatch |
| T3 | Convention fixture: a synthetic non-unitary representation and chain-complex instance, processed through the same parser and evaluator as the target. GREEN under declared conventions. Each of the four conventions (module side, vector convention, evaluation map, boundary direction) mutated separately; each mutation reddens at least one gate | Each mutation documented with the gate it reddens |

### 5.3 Derivation-path gates

| ID | Gate | Mutation |
| --- | --- | --- |
| D1 | ∂∂ = 0 in each evaluated representation: D₁D₂ = 0 and D₂D₃ = 0 | Perturb evaluated boundary entry; verify product nonzero |
| D2 | Evaluated boundary matrices recorded as derivation artifacts | SHA-256 of each matrix in output |
| D3 | Determinant factors: det(D₁[I',:]), det(D₂[J',I]), det(D₃[:,J]) recorded for each irrep | SHA-256 in output |
| D4 | Torsion product: T²_target recomputed from recorded determinant factors matches reported value | Perturb one determinant factor; verify T² changes |
| D5 | Galois consistency: for each Galois pair (ρ, ρ^σ), verify T²(ρ^σ) = σ(T²(ρ)) where σ: φ ↦ 1−φ | Perturb one T² value; verify Galois relation breaks |

### 5.4 Enumeration gate

| ID | Gate | Mutation |
| --- | --- | --- |
| E1 | Canonical enumeration SHA-256 matches the § 4.2 frozen hash: `27ff780d28d5d854d464ead87e8fc20244fac8334bda9f0600c6ee1b3c89561e` | Use wrong sort key (e.g., Q(φ)-normalized triple instead of 8-integer key); verify hash mismatch |

## 6. Representation construction

The 9 irreducible representations of 2I are constructed as follows:

| Label | Dim | Construction |
| --- | --- | --- |
| V₀ | 1 | Trivial: ρ(g) = 1 for all g |
| V₁ | 2 | Natural: quaternion (w,x,y,z) ↦ [[w+zι, x+yι],[−x+yι, w−zι]] (SU(2) embedding, ι = √(−1)) |
| V₂ | 3 | Sym²(V₁) in monomial basis |
| V₃ | 4 | Sym³(V₁) in monomial basis |
| V₄ | 5 | Sym⁴(V₁) in monomial basis |
| V₅ | 6 | Sym⁵(V₁) in monomial basis |
| V₆ | 4 | Galois conjugate of V₃: apply φ ↦ 1−φ to all Q(φ) components |
| V₇ | 2 | Galois conjugate of V₁: apply φ ↦ 1−φ |
| V₈ | 3 | Galois conjugate of V₂: apply φ ↦ 1−φ |

These are matched to the protocol's R0–R8 labels by the § 5.5 row signature after computation.

The Galois conjugates V₆, V₇, V₈ are obtained by applying the Galois automorphism σ: φ ↦ 1−φ
to every Q(φ) component of the representation matrices of V₃, V₁, V₂ respectively. Since
the quaternion coordinates of each group element are in Q(φ), and the SU(2) map is polynomial
in these coordinates, σ produces a distinct representation whenever the original has entries
not fixed by σ.

## 7. Overlap with the analytic route

Disclosed per § 6:

1. **Group structure:** both routes use the same group 2I with the same generators and
   multiplication. This is shared subject matter, not a shared method.
2. **Representation theory:** both routes evaluate the same irreps. The analytic route uses
   them to twist spectra; this route uses them to twist chain groups. The representations
   are the same mathematical objects applied differently.
3. **Q(φ) arithmetic:** both routes express results in Q(φ). The field is determined by the
   group's character ring, not by either route.
4. **No further overlap.** No spectral data, zeta functions, heat kernels, or analytic
   continuation appears in this route.
