# AUDIT_STAGE1: independent answers to worklist items 0–12

All values marked exact were computed in exact arithmetic. I wrote a small class, `mq.py`, for sums
Σ c_s √s with c_s ∈ ℚ(i). Zero tests in it are exact. No exact value in this report was identified from
floating point. Each exact value was also reproduced by an independent float64 quadrature route (item 12 script),
with agreement at 1e-14 or better. The scripts reproduce everything from a clean directory with
`./py run_all.py`, which copies the scripts to `clean_run/` and reruns them: 92 PASS, 0 FAIL, and the exact outputs
are identical. `./py build_results.py` then writes `audit_results.json`.

Every printed PASS line comes with a recorded mutation of the thing it checks, and the check fails on that mutant.
The logs are in `out/checks_*.json` and in the script outputs.

**Conventions I fixed.** These are listed in METHOD.md.
- Quaternion to SU(2): q = w+xi+yj+zk ↦ [[w−iz, −y−ix], [y−ix, w+iz]].
- D^j is Sym^{2j} of the defining representation in the Condon–Shortley orthonormal basis. The phases were checked
  exactly: D^{1/2}(g) = g, and J₊ has positive coefficients.
- Rotations act on sections by (r·ψ)(g) = ψ(r⁻¹g). On fibre vectors at every level this is u ↦ conj(D^j(r)) u.

Every quantity asked for is a function of the block fibre vector u and the sector. A different quaternion
identification conjugates Γ inside SU(2) or complex-conjugates it. A right translation g ↦ gr carries the sections
of one to the other. It leaves every fibre vector unchanged and commutes with Δ, N and all projections. So no
reported value depends on the identification, as the worklist says.

**"P-trick".** I replace Φ = uᵀD(g)η by Φη† = uᵀD(g)P_σ, where P_σ is the exact isotypic projector. Right
multiplication by the isometry η† preserves pointwise norms and commutes with N, DN, Δ, the level projections and
inner products. So all exact values depend only on P_σ. No square roots enter from normalizing η. The explicit η's
are used only for item 1 and for the independent quadrature route.

---

## Item 0

(a) ‖q₁‖² = 1 and ‖q₂‖² = 1, exactly, in ℚ(√5). The closure of {q₁, q₂} under multiplication has **|Γ| = 120**
elements. There are 120 distinct commutators, and they generate a subgroup of order 120. So **Γ equals its derived
subgroup: Γ is perfect.** There are 9 conjugacy classes, with sizes 1, 1, 12, 12, 12, 12, 20, 20, 30.

(b) **⟨3 3; 3 −3 | 6 0⟩ = √231/462 = 1/√924 ≈ 0.0328976.** I computed it from my own Racah-formula
implementation and again from an independent highest-weight/lowering construction; the two agree exactly. It is
positive, as the Condon–Shortley convention requires.

## Item 1

**Declaration: DERIVED.** Everything below was computed from the two generators.
1. Enumerate Γ exactly.
2. Form the class sum Z of q₁'s class (20 elements) in D³. Verify Z² = 5Z exactly, so P₄ = Z/5 and P₃ = I − Z/5
   are the two isotypic projectors, with traces 4 and 3. Verify P² = P, P† = P, and that P commutes with D³(q₁)
   and D³(q₂).
3. Check irreducibility: ⟨χ,χ⟩ = 1 exactly for χ_σ(h) = tr(P_σD³(h)). Both characters are real-valued.
4. For every level, dim Hom(σ, V_j) = (1/120) Σ_h conj(χ_σ(h)) χ_j(h), with χ_j(h) = U_{2j}(Re h) (Chebyshev
   polynomial of the second kind), computed exactly.

From its order and generators I recognize Γ as the binary icosahedral group. I used no character table or other
tabulated data.

| n | 0 | 1..5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3-dim sector | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 1 |
| 4-dim sector | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 2 |

Every odd level is 0, because −1 ∈ Γ acts by (−1)^n on V_{n/2} and trivially on σ. Neither sector has a level below
6, so the block is the lowest eigenspace in both sectors.

**Intertwiners.** Computed in mpmath at 50 significant digits.
- η₃ is an orthonormal basis of range P_σ. This fixes σ(h) := η₃†D³(h)η₃, used at every level.
- At each level, a basis comes from group averaging E = (1/120) Σ_h D^j(h) X D³(h)† applied to η₃. It is then
  Gram–Schmidt orthonormalized in ⟨η, η′⟩ = tr(η†η′)/d.
- Required precision: every residual below 1e-40.
- The check also caught a deliberately wrong group element (D(q₁)η compared with ησ(q₂)).

Largest residuals, as [max |D(h)η − ησ(h)| over h = q₁, q₂; max |η†η − I|; max |η†η′|]:

| sector, n | residuals |
|---|---|
| 3, n=6 | 6.9e-51, 5.9e-51 |
| 3, n=10 | 1.1e-50, 6.7e-51 |
| 3, n=14 | 3.4e-50, 5.7e-51 |
| 3, n=16 | 2.3e-49, 6.2e-50 |
| 3, n=18 | 9.6e-50, 5.8e-50 |
| 4, n=6 | 9.9e-51, 1.5e-50 |
| 4, n=8 | 8.7e-51, 5.4e-51 |
| 4, n=12 | 3.0e-50, 3.2e-50 |
| 4, n=14 | 4.2e-50, 2.1e-50 |
| 4, n=16 | 7.8e-50, 4.7e-50 |
| 4, n=18 (two copies) | 1.2e-49, 1.3e-50; η†η′: 2.7e-50 |

## Item 2

Φ is normalized, so its fibre vector has ‖u‖² = 7/d. The fibre of Π₆N(Φ) is computed exactly as F(u). The test is
whether F(u) − (u†F/u†u)u vanishes exactly. Q is the multiple; it equals ∫|Φ|⁴ (both computed). r̂₆ is computed
exactly. Stationarity uses the exact first derivative of r̂₆ in all 14 real directions.

**At all six points and in both sectors, Π₆N(Φ) = QΦ, and all six are stationary points of r̂₆.**

| U | r̂₆ | stationary | Q (3-dim sector) | Q (4-dim sector) |
|---|---|---|---|---|
| U1 | 1/924 | yes | 1288/1287 | 2289/2288 |
| U2 | 100/231 | yes | 1687/1287 | 168/143 |
| U3 | 24/77 | yes | 175/143 | 161/143 |
| U4 | 463/924 | yes | 1750/1287 | 2751/2288 |
| U5 | 9/35 | yes | 77/65 | 287/260 |
| U6 | 200/903 | yes | 5831/5031 | 609/559 |

**A structural fact, proved exactly (s_verify.py):** on the unit sphere of the block,

∫|Φ|⁴ = 1 + β_σ r̂₆, with β₃ = 28/39 and β₄ = 21/52.

Proof. ∫|Φ_u|⁴ is a Hermitian quartic form in u. It is invariant under u ↦ conj(D(r))u, which has the same
invariants as u ↦ D(r)u. The SU(2)-invariant Hermitian quartics on V₃ form a 4-dimensional space ≅
End_SU(2)(Sym²V₃), because Sym²V₃ = V₆⊕V₄⊕V₂⊕V₀ is multiplicity-free. I evaluated the basis
‖[u⊗u]_K‖² (K = 0, 2, 4, 6) at my 7 points (U1–U6 and item 11's point). The evaluation matrix has rank 4, so
evaluation at these points is injective on that space. The identity holds exactly at all 7 points, hence
everywhere.

Consequence: in each sector, Φ is critical for ∫|Φ|⁴ (that is, Π₆N(Φ) = QΦ) exactly when u is critical for r̂₆, and
the Hessians are proportional. Item 6 confirms this: the form equals (β_σ/4) × the r̂₆ second derivative.

## Item 3

The stabilizer of [u] was computed in two parts:
- The continuous part is the exact kernel of X ↦ P_⊥(i conj(J_X) u).
- The finite part comes from enumerating all rotations that permute the Majorana points (roots of the binary sextic
  of u), then checking D(r)u ∝ u. Residuals are ≤ 4e-15.

"Character" means the λ with r·Φ = λΦ. The last column is the complex dimension of the block subspace on which the
stabilizer acts by that character.

| U | stabilizer in SU(2) | character on Φ | dim |
|---|---|---|---|
| U1 | maximal torus T = {exp(−iθJ_z)} | e^{3iθ} | **1** (ℂv₃) |
| U2 | N(T) = T ∪ T·r_y(π) | T trivial; T·r_y(π) → −1 | **1** (ℂv₀) |
| U3 | binary octahedral, order 48; image: octahedron with vertices ±e_z, (±1,±1,0)/√2 | sign character of O ≅ S₄: +1 on the tetrahedral part, e.g. r_{(√2,0,1)}(2π/3) and r_{(1,1,0)}(π); −1 on r_z(π/2) and r_x(π) | **1** |
| U4 | binary dihedral, order 24; image D₆, axis z; Majorana points form a hexagon on the equator | r_z(π/3) → −1; r_y(π) (through vertices) → +1; r_x(π) (through edge midpoints) → −1 | **1** |
| U5 | cyclic of order 10, {±r_z(2πk/5)} | r_z(2πk/5) → e^{4πik/5} | **2** (span v₂, v₋₃) |
| U6 | binary dihedral, order 12; image D₃, axis z | r_z(2π/3) → +1; r_x(π) → −1 | **2** (span v₀, v₃+v₋₃) |

The element −1 ∈ SU(2) lies in every stabilizer and acts trivially. For U3 to U6 the Lie-algebra stabilizer is 0.
The order census of the SO(3) images identifies the groups: O has 24 elements (9 of order 2, 8 of order 3, 6 of
order 4); D₆ has 12; C₅ has 5; D₃ has 6.

## Item 4

ξ = −g Σ_{n≠6} Π_nN(Φ)/(n(n+2) − 48), so ‖Π_nξ‖²/g² = ‖Π_nN(Φ)‖²/(n(n+2) − 48)². The values are exact at every level
of the sector. Entries in **bold 0** were computed exactly and found to be zero; item 9 proves each of them.

#### 3-dim sector (levels 10, 14, 16, 18)

| U | n=10 | n=14 | n=16 | n=18 |
|---|---|---|---|---|
| U1 | 7/195150384 | 18375/1464486053888 | 7/2749593600 | 245/609749017488 |
| U2 | 1225/48787596 | 8575/2860324324 | **0** | 171500/139734149841 |
| U3 | **0** | 735/220024948 | **0** | 1715/2822912118 |
| U4 | 8575/780601536 | 735/5720648648 | **0** | 1516795/812998689984 |
| U5 | 77/5475600 | 281211/189112352000 | 553/1591200000 | 30233/56458242360 |
| U6 | 319550/88158077163 | 12920075/7517877885232 | 16583/68316230736 | 2355742375/6059914391677302 |

#### 4-dim sector (levels 8, 12, 14, 16, 18; level 18 includes both intertwiner copies)

| U | n=8 | n=12 | n=14 | n=16 | n=18 |
|---|---|---|---|---|---|
| U1 | 63/5360582656 | 7/483225600 | 23625/7209777496064 | 49/28242739200 | 427/1927108005888 |
| U2 | **0** | **0** | 11025/14081596672 | **0** | 74725/110407229504 |
| U3 | **0** | 7/1830400 | 945/1083199744 | **0** | 26901/80296166912 |
| U4 | 1575/31719424 | 189/644300800 | 945/28163193344 | **0** | 2643557/2569477341184 |
| U5 | 3591/346112000 | 63/83200000 | 361557/931014656000 | 34839/147097600000 | 2371131/8029616691200 |
| U6 | 7875/859947712 | 3521/1587595776 | 16611525/37011091127296 | 116081/701717332992 | 4105722425/19152322028017152 |

**ξ transforms by the same character as Φ.**
- Proof: left rotations commute with −Δ, N and Π_n, and preserve the sector (item 9). So r·ξ(Φ) = ξ(r·Φ) = ξ(χΦ).
  Since N(cψ) = |c|²cN(ψ) and |χ| = 1, this equals χξ(Φ).
- Numerical check: all listed stabilizer generators, applied to the exact ξ at every level, U and sector, deviate by
  at most 1.2e-15.

## Item 5

Per g means the value of Π₆DN_Φ[ξ]/g. It does not depend on g. The along-Φ component is ⟨Φ, Π₆DN_Φ[ξ]⟩/g. It is
real in every case, and equals λ₄/g² (item 8).

| U | 3-dim: along Φ | 3-dim: perp norm | 4-dim: along Φ | 4-dim: perp norm |
|---|---|---|---|---|
| U1 | −608931967/36722893315680 | 0 | −368688201/38687492546560 | 0 |
| U2 | −1871763250/229518083223 | 0 | −15820875/15112301776 | 0 |
| U3 | −2203040/944518861 | 0 | −162530109/75561508880 | 0 |
| U4 | −1921938515/459036166446 | 0 | −226441775133/38687492546560 | 0 |
| U5 | −267786421/58544557500 | 7188839√39/81061695000 | −4797453339/2497901120000 | 19565553√39/192146240000 |
| U6 | −336158940460/150812349114141 | 120576113√115/1288994436873 | −11093213145/4965015608696 | 214430223√115/3055394220736 |

Squared perp norms. U5: 51679406167921/168487138365975000000 (3-dim) and 1148432592587427/2840013657395200000000
(4-dim). U6: 1671938888011708435/1661506658289542382018129 (3-dim) and
5287736861620418835/9335433844106948692381696 (4-dim).

Real components, per g:
- **U5:** along e_t, 7188839√39/81061695000 in the 3-dim sector and −19565553√39/192146240000 in the 4-dim sector.
  Along i·e_t: **0** in both.
- **U6:** along τ_x, −120576113√115/1288994436873 in the 3-dim sector and 214430223√115/3055394220736 in the 4-dim
  sector. Along τ_y: **0** in both.

So at U5 and U6 the whole perpendicular part lies along e_t (respectively τ_x). The quadrature route reproduces the
full fibre of Π₆DN_Φ[ξ] to 7e-17.

## Item 6

L is the real-linear operator E ↦ P_⊥Π₆DN_Φ[E] − QE on tangents. It was built exactly as a 14×14 real matrix and
is exactly symmetric. Its quadratic form is the Riemannian Hessian of ∫|Φ|⁴/4 on the sphere. "Cross term" means the
polarized bilinear value B(τ_x, τ_y) = ⟨τ_x, DN_Φ τ_y⟩_ℝ − Q⟨τ_x, τ_y⟩_ℝ. I also computed B(τ_y, τ_x); it is equal.
The r̂₆ cross term is the polarized Hessian, [f″((a+b)/√2) − f″((a−b)/√2)]/2.

| | 3-dim sector | 4-dim sector | r̂₆ second derivative (both sectors) |
|---|---|---|---|
| U5, e_t | −56/165 | −21/110 | −104/55 |
| U5, i·e_t | **0** | **0** | **0** |
| U5, cross (e_t, i·e_t) | **0** | **0** | **0** |
| U6, τ_x | 6440/18447 | 2415/12298 | 920/473 |
| U6, τ_y | 56/429 | 21/286 | 8/11 |
| U6, cross (τ_x, τ_y) | **0** | **0** | **0** |

In every entry, the form equals (β_σ/4) × the r̂₆ value (item 2). The explicit one-variable reductions, exact via
sympy:
- On U5's curve: r̂₆ = (−875S² + 840S + 36)/924 with S = sin²t. Its only critical point in (0, 1) is S = 12/25, a
  maximum along the curve.
- On U6's family: r̂₆(z) = (100|z|⁴ − 20x² + 148y² + 463)/(231(|z|² + 2)²). Its critical real x are 0 and
  ±√(23/10) = ±z₀.

## Item 7

**Free directions.** ker L on the tangent space was computed exactly. In all 12 cases it equals the span of the
infinitesimal rotations, projected perpendicular to Φ. This was checked exactly: each projected rotation direction
lies in the kernel, and the ranks match.
- U1 and U2: real dimension 2. These are the x- and y-rotations; the z-rotation acts on Φ by a phase.
- U3 to U6: real dimension 3. At U5, i·e_t is one of these directions: P_⊥(iJ_zu) = −5i cos t sin t · e_t.

κ is the unique solution that is ⟨·,·⟩_ℝ-orthogonal to ker L. It was solved exactly with a bordered system.

| | ‖κ‖²/g² | components of κ/g | remainder |
|---|---|---|---|
| U1–U4, both sectors | **0** (κ = 0) | – | – |
| U5, 3-dim | 1054681758529/396076284864000000 | e_t: 1026977√39/3930264000; i·e_t: **0** | **0** |
| U5, 4-dim | 2604155538747/234711872512000000 | e_t: −931693√39/1746784000; i·e_t: **0** | **0** |
| U6, 3-dim | 296706102575281/35935889967339860160 | τ_x: 17225159√115/64285514280; τ_y: **0** | **0** |
| U6, 4-dim | 104263765387369/7098447400956021760 | τ_x: −10210963√115/28571339680; τ_y: **0** | **0** |

**The order-a⁵ block equation's component orthogonal to Φ**, as squared norm per g²:
- Without κ: 0 at U1–U4. At U5 and U6 it is the item-5 squared perp norm, which is nonzero.
- With κ: **0** in all 12 cases. The solvability condition (the right-hand side is orthogonal to ker L) holds
  exactly. The quadrature route confirms the κ equation to 8e-18.

So at U5 and U6 a nonzero tilt is required: the expansion without κ is inconsistent at order a⁵.

## Item 8

λ₄ = g⟨Φ, Π₆DN_Φ[ξ + κ]⟩ − gQ⟨Φ, κ⟩. It was computed with κ and without κ, and the two agree exactly in all 12
cases.

| U | λ₄/g², 3-dim sector | λ₄/g², 4-dim sector |
|---|---|---|
| U1 | −608931967/36722893315680 ≈ −1.658e-5 | −368688201/38687492546560 ≈ −9.530e-6 |
| U2 | −1871763250/229518083223 ≈ −8.155e-3 | −15820875/15112301776 ≈ −1.047e-3 |
| U3 | −2203040/944518861 ≈ −2.332e-3 | −162530109/75561508880 ≈ −2.151e-3 |
| U4 | −1921938515/459036166446 ≈ −4.187e-3 | −226441775133/38687492546560 ≈ −5.853e-3 |
| U5 | −267786421/58544557500 ≈ −4.574e-3 | −4797453339/2497901120000 ≈ −1.921e-3 |
| U6 | −336158940460/150812349114141 ≈ −2.229e-3 | −11093213145/4965015608696 ≈ −2.234e-3 |

For reference, λ₂ = gQ with Q from item 2.

**Why the value with κ equals the value without.** Write DN_Φ[κ] = 2Re(Φ̄·κ)Φ + |Φ|²κ. Then:
- Re⟨Φ, DN_Φκ⟩ = 3Re⟨N(Φ), κ⟩ = 3Re⟨Π₆N(Φ), κ⟩ = 3Q Re⟨Φ, κ⟩ = 0.
- Im⟨Φ, DN_Φκ⟩ = Im⟨N(Φ), κ⟩ = Q Im⟨Φ, κ⟩ = 0.

This uses Π₆N(Φ) = QΦ (item 2) and ⟨Φ, κ⟩ = 0. The term −gQ⟨Φ, κ⟩ vanishes for the same reason.

**Sign: λ₄ < 0 for every real g ≠ 0, in all 12 cases.** The argument:
1. By the same computation, ⟨Φ, DN_Φξ⟩ = 3Re⟨N(Φ), ξ⟩ + i Im⟨N(Φ), ξ⟩.
2. Since ξ = −g Σ_{n≠6} Π_nN/μ_n with μ_n = n(n+2) − 48, this gives
   λ₄ = −3g² Σ_{n≠6} ‖Π_nN(Φ)‖²/μ_n.
3. This identity was verified exactly in all 12 cases. Mutant: drop one level from the sum, and the check fails.
4. By item 1 there is no level below 6 in either sector, so every μ_n > 0.
5. Π_nN(Φ) ≠ 0 at some n ≠ 6 in every case; item 4 has nonzero entries at n = 14 and 18 everywhere.
6. Therefore λ₄ < 0, independent of the sign of g.

## Item 9: every vanishing quantity in items 4–8 (and in item 2)

Every zero below was **computed**, in exact arithmetic, and found to be zero. None is a missing value. Each has an
exact argument; none is unresolved.

**Symmetries used, with proofs that they preserve both the equation and the sector.**

**(S1) Left rotations.** (r·ψ)(g) = ψ(r⁻¹g).
- Sector: (r·ψ)(gh) = ψ(r⁻¹g)σ(h), so the sector is preserved.
- Equation: −Δ commutes with r because the round metric is bi-invariant. N commutes because it is pointwise and
  |·| is invariant. λ and g are untouched. So the equation is preserved, and Π_n commutes with r too.
- Consequence: every quantity built equivariantly from Φ (ξ, Π₆DN_Φ[ξ], L, κ) transforms under Stab[Φ] by the
  character χ of item 3. For κ this is because it is the unique minimal-norm solution, and ker L is
  stabilizer-invariant.

**(S2) Antiunitary symmetry S′.**
- Construction. Let η′ := Θ_col(conj η), that is, η′_{k,a} = (−1)^k conj(η_{−k,a}). Then conj(Φ_u) = Φ_{Θu} with η′
  in place of η, and η′ is an intertwiner from conj σ.
- I verified exactly that P_σ Θ_col(conj P_σ) = Θ_col(conj P_σ). So η′ = ηM with M unitary, and
  W := M† ∈ Hom_Γ(σ, conj σ). Hence T(ψ) := conj(ψ)W is defined on the sector and gives T(Φ_u) = Φ_{Θu}.
- T preserves the equation: (T ψ)(gh) = (Tψ)(g)σ(h), so it preserves the sector. Δ is a real operator, and
  |conj(ψ)W| = |ψ|, so N(Tψ) = T N(ψ). λ and g are real.
- Now compose with r₀ = r_y(π), the quaternion j. I verified exactly, for all J ≤ 9, that D^J(r₀)Θx = (−1)^J conj(x).
- So S′ := −r₀∘T is an antilinear symmetry of the equation and of the sector. On block fibres it acts as u ↦ ū.
- All six u are real. So S′Φ = Φ, and ξ, Π₆DN_Φ[ξ], L and κ are all S′-invariant. I checked exactly that their
  block fibres are real.
- e_t and τ_x have real fibres; i·e_t and τ_y have imaginary fibres.

**(S3) Algebra of Sym³.** At U2, U3 and U4, Θu = c·u with c = 1, 1, −1 respectively (exact).
- The level-J fibre of N(Φ) is Σ_L A_L ⊗ Z_L, with A_L = [[Θu⊗u]_L ⊗ u]_J. So A_L = c times an SU(2)-equivariant
  image of u⊗u⊗u ∈ Sym³V₃.
- By weight counting, Sym³V₃ = V₉+V₇+V₆+V₅+V₄+2V₃+V₁. V₈ does not occur, so every A_L with J = 8 vanishes by Schur.
- No symmetry of the equation is used here; this is an identity of the cubic map.

**The zeros.**

| Zero | Where | Reason |
|---|---|---|
| Item 4: ‖Π_nξ‖² = 0 | U2: n = 16 (3-dim); n = 8, 12, 16 (4-dim) | (S1). By item 3, ξ's level-J fibre lies in the χ-subspace of V_J: weight 0 and r_y(π)-eigenvalue −1. r_y(π) acts on v₀ by (−1)^J, so that subspace is 0 for even J. |
| Item 4 | U3: n = 10, 16 (3-dim); n = 8, 16 (4-dim) | (S1). The octahedral χ-subspace of V_J has dimension 0 at n = 8, 10, 16 (item 3). n = 16 also follows from (S3). |
| Item 4 | U4: n = 16, both sectors | (S3). Here the D₆ χ-subspace of V₈ is 1-dimensional, so (S1) alone does not explain this zero. |
| Item 2: Π₆N − QΦ = 0; r̂₆ gradient = 0 | U1–U4 | (S1). The χ-subspace of the block is ℂΦ, and the gradient of an invariant function is χ-covariant. |
| Item 2 | U5, U6 | See below. |
| Item 5: perp part of Π₆DN_Φ[ξ] = 0 | U1–U4 | (S1). It is χ-covariant, and the χ-subspace of the block is ℂΦ. |
| Item 5: components along i·e_t (U5) and τ_y (U6) | both sectors | (S2). The vector is S′-invariant (real fibre), and Re⟨imaginary, real⟩ = 0. |
| Item 6: form on i·e_t; cross term (e_t, i·e_t) | U5 | i·e_t ∝ P_⊥(iJ_zu) is an orbit direction. At a critical point the Hessian of an invariant function vanishes on orbit directions, and L annihilates them (exact kernel). The same holds for r̂₆. |
| Item 6: cross term (τ_x, τ_y) | U6 | (S2). L maps real fibres to real fibres, and the imaginary τ_y is orthogonal to them. For r̂₆: r̂₆(ū) = r̂₆(u), since the CG coefficients are real and Θū = conj(Θu). |
| Item 7: κ = 0 | U1–U4 | (S1). κ is χ-covariant and ⊥ Φ, and the χ-subspace of the block is ℂΦ. |
| Item 7: κ components along i·e_t, τ_y | U5, U6 | (S2). |
| Item 7: κ remainder = 0 | U5, U6 | (S1). κ lies in (χ-subspace) ∩ Φ^⊥ = ℂe_t (resp. ℂτ_x), which is the real span of e_t and i·e_t (resp. τ_x and τ_y). |
| Item 7: perp equation with κ = 0 | all | See below. |
| Item 8: λ₄(with κ) − λ₄(without κ) = 0 | all | The exact identity in item 8. |

**Item 2 at U5 and U6.**
- The residual Π₆N(Φ) − QΦ is complex-orthogonal to Φ and χ-covariant. So it is a multiple of e_t (U5) or τ_x
  (U6).
- Its real e_t (τ_x) component is ¼ of the derivative of ∫|Φ|⁴ along e_t (τ_x). By item 2 this is
  ¼β_σ d r̂₆/dt, which is 0 at S = 12/25 by the explicit curve formula (x-derivative at z₀ for U6).
- Its i·e_t component is a derivative along an orbit direction, so it is 0. Its τ_y component is 0 by (S2), or by
  ∂_y r̂₆ = 0 at y = 0.

**Item 7: why the perpendicular equation with κ is 0 in all cases.**
- ker L is exactly the orbit tangent (items 6/7), and L is invertible on its orthogonal complement.
- The right-hand side P_⊥Π₆DN_Φ[ξ] is orthogonal to the orbit tangent. Reason: the full equation is the gradient of
  a G-invariant functional, so ⟨Xψ, F(ψ, λ)⟩_ℝ ≡ 0 for X ∈ su(2) ⊕ iℝ.
- At order a⁶ this identity gives ⟨XΦ, F₅⟩_ℝ = 0 once F₃ = 0. Since Re⟨XΦ, Φ⟩ = 0, this leaves
  ⟨XΦ, Π₆DN_Φ[ξ]⟩_ℝ = 0.

## Item 10: existence, uniqueness, regularity

**Answer: yes, for all six points U1–U6, in both sectors.** For every sufficiently small a > 0 there is an actual
solution (ψ(a), λ(a)) whose expansion begins as in §2. For U1–U4 the tilt is κ = 0; for U5 and U6 it is the nonzero
κ of item 7. My argument follows.

**Setting.** Let 𝓗 be the real Hilbert space of L² sections of the chosen sector, with ⟨·,·⟩_ℝ, and let
𝓧 = H² sections. The equation is F(ψ, λ) := (−Δ − λ)ψ + gN(ψ) = 0. With respect to ⟨·,·⟩_ℝ this is the gradient of

E_λ(ψ) = ½⟨ψ, (−Δ − λ)ψ⟩ + (g/4)∫|ψ|⁴,

because N is the gradient of ¼∫|ψ|⁴. The compact group G = SU(2)_left × U(1)_phase acts by isometries,
preserves the sector and commutes with F (proved in item 9, (S1), plus the phase invariance of N). So E_λ is
G-invariant.

**Step 1: Lyapunov–Schmidt with blow-up.** Let B be the block. Write ψ = a(Φ̂ + a²w) with Φ̂ = Φ + κ̃, where κ̃ ∈ B is
complex-orthogonal to Φ, and w ⊥ B. Write λ = 48 + a²μ. Then a = ⟨Φ, ψ⟩, the normalization of §2. Since
(−Δ − 48)Φ̂ = 0, dividing F by a³ gives the equation

G(w, κ̃, μ; a²) := (−Δ − 48)w − μ(Φ̂ + a²w) + gN(Φ̂ + a²w) = 0.

This is polynomial in (w, κ̃, μ) and in a² alone.
- Its component in B^⊥ is (−Δ − 48 − a²μ)w + g(1 − Π₆)N(Φ̂ + a²w) = 0. The operator −Δ − 48 is invertible from
  B^⊥ ∩ H² onto B^⊥. The spectrum on B^⊥ is {n(n+2)} over the sector's levels n ≠ 6; the nearest is n = 10 (3-dim
  sector) or n = 8 (4-dim), which gives a gap of 72 or 32.
- N is a real-analytic map of H² into itself, since H² is an algebra in dimension 3.
- So the analytic implicit function theorem gives a unique w = w(κ̃, μ, a²) near w₀ = ξ(Φ̂)/g, analytic in all
  arguments.

This blow-up handles the degeneracy at a = 0. There, ψ = 0 solves the equation for every λ, and the linearization at
λ = 48 has a kernel of real dimension 14. After dividing by a³, the a = 0 problem is regular off B.

**Step 2: the reduced equation.** What remains is the block equation
R(κ̃, μ; a²) := −μΦ̂ + gΠ₆N(Φ̂ + a²w) = 0 in B (real dimension 14).
- Lyapunov–Schmidt preserves the gradient structure. R is the B-gradient of the reduced function
  f(Φ̂, μ; a²) = E(ψ)/a⁴, restricted to the graph of w. This function is invariant under G.
- The real Φ-component of R fixes μ: μ = g⟨Φ, Π₆N(·)⟩ + O(a²), with coefficient −1 in μ, so it is always solvable.
  At a = 0 this gives μ = gQ (items 2 and 8).
- The iΦ-component and the orbit components are identically satisfied: ⟨Xψ, F⟩_ℝ ≡ 0 for all X in the Lie algebra
  of G.
- So by the slice theorem (Palais' symmetric criticality), it is enough to solve the equation projected onto the
  slice S = T ∩ (orbit tangent)^⊥, with κ̃ restricted to S. Here T = Φ^⊥ ∩ B. Near Φ every G-orbit meets S, and a
  critical point of the restriction of an invariant function to S is critical.

**Step 3: IFT on the slice.**
- At a = 0, κ̃ = 0 and μ = gQ, the slice equation holds, because Π₆N(Φ) = QΦ (item 2).
- Its linearization in κ̃ is g·L restricted to S. L is invertible there: ker L on T is exactly the orbit tangent
  (items 6/7, exact kernel computation), and g ≠ 0.
- The analytic IFT gives a unique κ̃(a²) ∈ S near 0, together with μ(a²), for |a| < a₀. At order a² this κ̃ is the
  item-7 κ. That κ is orthogonal to ker L, which is the slice condition.
- So ψ(a) = a(Φ + κ̃(a²) + a²w(…)) is an actual solution. Its expansion is
  ψ = aΦ + a³(κ + ξ) + O(a⁵), λ = 48 + gQa² + λ₄a⁴ + O(a⁶), with λ₄ from item 8.

**Regularity.** a ↦ ψ(a) is real-analytic from (−a₀, a₀) into H², and in fact into C^∞ by elliptic regularity. It
is odd: ψ(−a) = −ψ(a), because G depends on a only through a². a ↦ λ(a) is real-analytic and even, with
λ − 48 = gQa² + O(a⁴) and Q > 0.

**Local uniqueness proved.** There are δ, a₀ > 0 with the following property. Take any solution (ψ, λ) in the sector
with 0 < ‖ψ‖_{H²} < a₀, |λ − 48| < δ, and ψ/‖ψ‖ within δ of the G-orbit of Φ. Then (ψ, λ) = (γ·ψ(a), λ(a)) for some
γ ∈ G and some a ∈ (0, a₀).
- Note that λ − 48 = O(‖ψ‖²) follows from the reduced equation, so it need not be assumed.
- If in addition ψ lies in the slice, that is ⟨Φ, ψ⟩ = a > 0 and the block part minus aΦ is orthogonal to the orbit
  tangent, then ψ = ψ(a) exactly. So uniqueness holds modulo G = SU(2) × U(1), in the stated neighbourhood (cone
  around the orbit, small amplitude).
- The branch is fixed by the stabilizer subgroup Ĥ = {(r, χ(r)⁻¹)}, because the construction is Ĥ-equivariant and
  unique.
- For U1–U4 an alternative proof restricts to Fix(Ĥ), whose block part is ℂΦ (item 3). There the reduced problem
  is the single μ-equation, and existence does not need item 6/7's nondegeneracy. That proof gives uniqueness only
  within Fix(Ĥ).

**Hypotheses, where each is verified, and what fails without it.**

1. **(H1) The kernel of −Δ − 48 on the sector is exactly the level-6 block, of complex dimension 7.**
   - Verified in item 1: dim Hom(σ, V₃) = 1, and n(n+2) = 48 only for n = 6.
   - Without it, the reduction space is wrong, and extra kernel directions would enter the bifurcation equation.
2. **(H2) There is a spectral gap on B^⊥; the nearest levels are 10 and 8.**
   - Verified in item 1.
   - Without it, the range equation of Step 1 is not uniquely solvable.
3. **(H3) Π₆N(Φ) = QΦ.**
   - Verified in item 2, exactly, and by quadrature.
   - Without it, the reduced equation at a = 0 has no solution at Φ, and no branch has leading term aΦ.
4. **(H4) Nondegeneracy modulo symmetry: ker L on T equals the orbit tangent.**
   - Verified in items 6/7: exact kernel, of dimension 2 or 3, equal to the rotation directions.
   - Without it, the IFT of Step 3 fails. Persistence, and in particular uniqueness modulo G, would then need
     higher-order terms such as item 5's tangential component and W₆.
5. **(H5) G acts on the sector and commutes with F.**
   - Verified in item 9 (S1) and by the item-4 character check.
   - Without it, the orbit components are not automatically satisfied. The count of equations (14) against unknowns
     (13 − dim orbit) then fails, and so does the slice reduction.
6. **(H6) Variational structure: F = ∇E_λ with E_λ invariant, and L symmetric.**
   - L's exact symmetry is verified in items 6/7.
   - Without it, the identities ⟨Xψ, F⟩ ≡ 0 and the symmetric-criticality step are not available.
7. **(H7) N is analytic on H².**
   - This is the standard Sobolev algebra property; it is not verified in items 1–9.
   - Without it, the IFT does not apply.
8. **(H8) g is real and nonzero.**
   - This is given in §2.
   - If g = 0, L is multiplied by zero, the problem is linear, and there is no bifurcating branch structure.
9. **Q > 0.**
   - Verified in item 2: Q ≥ (∫|Φ|²)² = 1.
   - It is not needed for existence. It is needed for λ to parametrize the branch: λ − 48 has the sign of g.

**Degeneracies at a = 0 that the argument handles.**
1. The trivial branch ψ ≡ 0 for all λ, and the 14-dimensional kernel at λ = 48: handled by the blow-up and division
   by a³ in Step 1.
2. The phase and rotation orbit directions: they make the linearized reduced equation singular at every a. They are
   handled by the slice and by the identically satisfied orbit components. For U1 and U2 the stabilizer is
   continuous: a torus acts on Φ by a phase. So only 2 rotation directions are free, and the orbit of [Φ] has real
   dimension 3 including the phase.
3. At U5 and U6, the correct next-order direction is Φ + a²κ with κ ≠ 0. Without the tilt the a⁵ equation fails;
   item 7 shows it is nonzero without κ.

## Item 11

On U5's curve at sin²t = 1/4, with e_t = −½v₂ + (√3/2)v₋₃ (normalized):
- **It is not a stationary point of r̂₆.** r̂₆ = 3061/14784. The exact gradient has nonzero components 115/352 and
  −115√3/1056. The curve derivative is dr̂₆/dS = −(1750S − 840)/924 ≠ 0 at S = 1/4. Also, Π₆N(Φ) is not a multiple
  of Φ in either sector.
- **Re⟨Φ, DN_Φ[e_t]⟩ = 805√3/6864 in the 3-dim sector, and 2415√3/36608 in the 4-dim sector.**
- Cross-check: this equals 3Re⟨N(Φ), e_t⟩, computed separately. It also equals (3β_σ/4) dr̂₆/dt, which is
  (3/4)(28/39)(805√3/3696) = 805√3/6864 in the 3-dim sector.

## Item 12

For all 14 cases (U1–U6 and the item-11 point, in both sectors), N(Φ) was evaluated pointwise in float64 on a Haar
quadrature grid on SU(2).
- **Grid.** Euler angles; 24 uniform nodes in α and in γ, and 12 Gauss–Legendre nodes in cos β. This rule is exact
  for every integrand of spin ≤ 23 in α and γ, and of degree ≤ 23 in cos β after the α, γ integration.
- **Projection.** N(Φ) was projected onto every level n ≤ 22, keeping every intertwiner copy (both copies at n = 18
  in the 4-dim sector).
- **−Δξ.** Computed as 4·Casimir, with Casimir = J_x² + J_y² + J_z² built from explicit spin matrices and acting on
  the left index. It was applied to item 4's exact ξ, converted to explicit intertwiner coordinates.

**Result.** max |(−Δ − 48)ξ + g(N − Π₆N)| over all coefficients, levels and cases = **8.5e-15**. I required < 1e-10.
Checks and mutants:
- Mutant: flip the sign of ξ, and the residual is ≥ 0.037. Mutant: replace 4·Casimir by 3·Casimir, and the residual
  is ≥ 0.0065.
- The residual from reconstructing N with every copy is ≤ 2.6e-14. This includes the levels without intertwiners and
  n = 20, 22. Mutant: drop the second copy at n = 18, and the residual is ≥ 0.005.
- Parseval: ∫|N|² − Σ_n ‖Π_nN‖² ≤ 1.8e-14.
- A second grid (30/16/30) agrees to 1.7e-14.
- The quadrature level norms agree with the exact ones to 1.8e-14.

**Why N(Φ) has no component above level 18.** N(Φ)_a = Σ_b conj(Φ_b)Φ_bΦ_a is a sum of products of three spin-3
matrix coefficients; the conjugate is again a spin-3 coefficient. By the Clebsch–Gordan series these span only D^J
with J ≤ 9, that is, levels n ≤ 18. Parseval and the vanishing of the n = 20 and 22 projections confirm this
numerically.

---

## Underdetermined points and the readings I took

1. **Rotations fixing Φ (item 3).** Reported as subgroups of SU(2), for the action (r·ψ)(g) = ψ(r⁻¹g). −1 is always
   included and acts trivially. The "character" is the scalar by which r acts on Φ.
2. **"The cross term" (item 6).** Taken as the polarized bilinear value B(τ_x, τ_y), not 2B; for r̂₆, the polarized
   Hessian. All of these are 0, so the reading does not change any value.
3. **"Every level n of the sector" (item 4).** Taken as n ≤ 18 with Hom ≠ 0, excluding n = 6. ξ has no other
   components.
4. **"Per g" (item 5).** Taken as Π₆DN_Φ[ξ]/g, which does not depend on g. Items 7 and 8 are per g and per g².
5. **Q in item 6.** Taken as item 2's multiple for the normalized Φ, which is ∫|Φ|⁴.
6. **Item 10 uniqueness.** Stated modulo G = SU(2) × U(1) in a cone neighbourhood; see above. The worklist does not
   fix the neighbourhood.
7. **Stationarity (item 2).** Taken on the unit sphere of V₃, in all 13 real tangent directions. I computed all 14
   real directional derivatives; the one along iu is trivially 0.

## Limitations

- The finite stabilizers (item 3) were enumerated in float64, with residuals around 1e-15, and are not certified
  exactly. The exact values of items 2 and 4–8 do not depend on them; they enter only the zero arguments. Those
  arguments are also backed by the exact computations: every zero was computed exactly.
- Hypothesis (H7), the analyticity of N on H², is a standard fact and is not computed.

## Files

- **Method:** `METHOD.md`, with a dated note appended.
- **Library code:** `mq.py` (exact field), `su2.py` (CG in two constructions, D^j, quaternions), `engine.py` (exact
  section calculus), `linalg_mq.py`, `numD.py`, `checks.py`.
- **Scripts:** `s01_group.py` (items 0, 1), `s02_intertwiners.py` (item 1), `s_main.py` (items 2, 4–8, 11),
  `s03_stabilizers.py` (items 3, 4), `s12_quadrature.py` (item 12 and cross-checks), `s_verify.py` (structural
  proofs, item 11).
- **Drivers:** `run_all.py` (clean-directory reproduction), `build_results.py`.
- **Outputs:** `out/`, `audit_results.json`.

## Consulted-material manifest

- **Read:** `BRIEF.md` and `worklist.md` in this room. Nothing else was read. No web, papers, books, or other files
  outside this room.
- **Software:** Python 3.12 standard library (`fractions`, `math`, `itertools`, `json`, `pickle`), numpy (arrays,
  `roots`, `eigh`, `svd`, `leggauss`), mpmath (50-digit arithmetic), and sympy. sympy was used only for simplification
  and differentiation of the explicit r̂₆ formulas and for printing square roots. No library Clebsch–Gordan, Wigner-D
  or group-table routine was used. CG coefficients come from my own Racah formula, cross-checked against my own
  highest-weight construction: 59 coupling triples, 3541 coefficients, 0 mismatches.
- **Item 1 declaration: DERIVED.** The Hom dimensions come from the group enumerated from q₁ and q₂, class-sum
  projectors, and characters computed exactly. I recognized Γ as the binary icosahedral group, but used no tables.
