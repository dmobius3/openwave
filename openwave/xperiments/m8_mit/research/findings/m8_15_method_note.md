# M8.15 method note: the soft-slot branches at finite amplitude

Task: [`tasks/m8_15_task_details.md`](../tasks/m8_15_task_details.md), governing SHA-256 `e0c32d01…`, pre-registered in [#618](https://github.com/openwave-labs/openwave/pull/618), run 2026-10-07 to 2026-10-08 from the frozen solver at MIT `6058f95`. The record is per branch, in [`RECORD.md`](https://github.com/dmobius3/mode-identity-theory/blob/28c8a7648a88c9234d32be753cbbe3563ea2eec3/files/framework/files/working/files/scripts/soft-slot-dynamics/RECORD.md). Adjudication is the maintainer's.

## 1. Equations first

`R = c = g = 1`. A complex vector `Ψ ∈ ℂⁿ` is stored in real form as `(Re Ψ, Im Ψ) ∈ ℝ²ⁿ`, where `J₀` is multiplication by `i`. `⟨a, b⟩_M = aᵀMb` and `‖a‖_M² = aᵀMa`.

### 1.1 The continuum problem

```text
ψ_TT − Δψ + |ψ|²ψ = 0,      ψ : S³ → ℂ^d,   ψ(γx) = ρ(γ) ψ(x)  for γ ∈ 2I acting on the left
ψ = e^{iωT} φ,  σ = ω²:     −Δφ − σφ + |φ|²φ = 0
ε = σ − λ̄₀                    λ̄₀ = the mean of the slot's lowest discrete level at the level run
```

The soft slots are `R4` and `R5`, with lowest level `K = 6` and `λ₀ = 48`, and `R2`, with `K = 7` and `λ₀ = 63`. Each branch is followed on the ladder `ε_k = 0.25 · 2^{k/2}` up to `ε_max = λ₁ − λ₀`, which is 72 in `R4`, 32 in `R5` and 132 in `R2`, with `ε_max` itself appended.

### 1.2 The discretization

The mesh is the canonically refined 600-cell, invariant under left and right 2I. The elements are exact-geometry P2: flat tetrahedra projected radially onto `S³`, with the metric evaluated exactly at the quadrature points. The unknowns are the values at one canonical point per left orbit of nodes, and the equivariance enters through the element blocks:

```text
K = Σ_e ∫_e ∇N_i · ∇N_j dV · ρ(g_i)†ρ(g_j)          M = Σ_e ∫_e N_i N_j dV · ρ(g_i)†ρ(g_j)
V₄(Ψ) = ¼ Σ_q w_q |ψ_h(x_q)|⁴                       f = ∇V₄,  dV₄ = Re⟨f, dΨ⟩
Df(Ψ)[δ] = A δ + C conj(δ)                          A Hermitian, C symmetric
M Ψ'' + K Ψ + f(Ψ) = 0                              H = ½ Π†MΠ + ½ Ψ†KΨ + V₄(Ψ),   N = Im(Ψ†MΠ)
```

Every form omits the common factor 120, which comes from integrating one representative element per 2I orbit. The right-rotation derivatives `X_a`, `a = 1, 2, 3`, are the `L²` projection of the derivative along the left-invariant fields `x ↦ x e_a/2`, assembled with a 5-point rule.

### 1.3 Standing waves: the pinned, bordered Newton

Each representative is oriented so that its pinned subgroup `H` lies in the mesh's right 2I. Newton then runs in `H`'s fixed space:

```text
Fix_H = {Ψ : (conj(a_h) − b_h J) R_h Ψ = Ψ for all h ∈ H}       (b_h = 0 outside R2)
x = x₀ + W_H z,   W_H a basis of Fix_H
W_Hᵀ F = 0,   F = F₀(x, σ) + Σ_k μ_k M p_k + ν M s,   F₀(x, σ) = K x − σ M x + f(x)
xᵀ M x = m,   p_kᵀ M x = 0,   sᵀ M (x − x_start) = 0
stop when   ‖F₀(x, σ) + ν M s‖ / ‖M x‖ + |xᵀ M x − m| / m < 10⁻¹¹,   within 40 iterations
```

- **The gauges `p_k`.** These are the exact fibre generators at the start of the solve: `iΦ`, and in `R2` also `JΦ` and `iJΦ`. Each is kept if it projects onto `Fix_H`.
- **The slice `s`.** It is used for T4 only. It is the discrete right-rotation derivative about the pinned 2-fold axis, applied to the seed, with multiplier `ν`.
- **The stopping test.** It keeps the slice force `ν M s` as part of the solution, since the terms record it with the floors. It leaves out the gauge terms, on the terms' premise that "the gauge multipliers vanish at a solution". § 5.2 records where that premise fails.

### 1.4 Continuation

```text
ε = σ − λ̄₀ is an output; continuation is in m
the seed, at ε_min:  P_H(germ interpolant), scaled to m₀ = ε_min / (120 Q / vol), Q the seed's quartic ratio,
                     vol 120 times the summed element volume
to ε_k:              m from a secant through the last two accepted (m, ε)  (m ∝ ε at the first step),
                     Newton at m, iterated until |ε − ε_k| ≤ 10⁻⁶ ε_k, at most 8 times
Newton failure:      the step in m retried at 1/2, then 1/4; a failure after both ends the branch
fold:                the first accepted ε below its predecessor while m increased; ε(m) maximized by
                     golden section within the last three accepted points, to 10⁻³ relative in m; the branch ends there
```

At the seed there is no previous accepted point, so there is no step to retry. A Newton failure there ends the branch: `Newton failure at the seed`.

### 1.5 The linearization

```text
rotating frame:                 M(ζ'' + 2iω ζ') + L ζ = 0,       L = K − σM + (A + C conj)
first order, x = (ζ, ξ = ζ'):   𝒜 x = λ ℬ x,   𝒜 = [0, I; −L, −2ω M J₀],   ℬ = [I, 0; 0, M]
Ω(x₁, x₂) = ζ₁ᵀMξ₂ − ξ₁ᵀMζ₂ + 2ω ζ₁ᵀ M J₀ ζ₂        B(x₁, x₂) = ζ₁ᵀLζ₂ + ξ₁ᵀMξ₂
Krein value of an eigenvector v:   Re(v† B v);  its sign is the Krein signature
count:   k_r + 2 k_c + 2 k_i⁻ = n(B) − n(B|G),   n(B) = n(L)
```

Shift-invert near 0 computes the `2(j₂ + 1) + 4` eigenvalues nearest 0. The `2(j₂ + 1)` smallest in modulus are the slow set, and the next one gives the slow gap. `n(L)` comes from a complete slice of `L`'s spectrum relative to `M`, upward from `λ_min − σ`. `n(B|G)` is the negative index of `B`'s Gram matrix on the cluster `G`.

### 1.6 The linear verdict (level 3, cross-checked at level 2)

```text
U = span_M{ iΦ (and JΦ, iJΦ in R2), the phase partner's static part, X_aΦ (a = 1, 2, 3);
             in C1 also the rigid family's two tangent directions }
an eigenvector belongs to G if its static part has M-overlap ≥ 0.9 with U; otherwise it is genuine
G also holds the exact fibre vectors and the phase's Jordan pair: (iΦ, 0), (2ω ∂_σΦ, iΦ), and in R2 (JΦ, 0), (iJΦ, 0)
cluster checks:   dim G = k + r_u (frozen per case);  G closed under λ → −λ, λ̄, −λ̄;  Ω-Gram σ_min/σ_max ≥ 10⁻⁶
representatives:  Re λ ≥ −τ_Re and Im λ ≥ −τ_Re
classes:          real (Re λ > τ_Re, |Im λ| ≤ τ_Re), quartet (Re λ > τ_Re), imaginary pair of each Krein sign
floors:           τ₀ = τ_Re = 10 r_f,   r_f⁽³⁾ = 1.040×10⁻⁶,  r_f⁽²⁾ = 3.970×10⁻⁵
lifted modes:     paired |λ| ≤ τ₀;  radical |λ| < ½ max |λ_pred(ε)|;  a genuine |λ| ≤ τ₀ is a zero crossing
kinds:            every lifted mode is radical at T2, T3, T5, C3 and C3s; at T4, a lifted mode whose static part
                  has squared M-overlap ≥ 0.5 with the slice is radical; every other lifted mode is paired
eigenpairs:       ‖𝒜x − λℬx‖ / ((‖𝒜‖₁ + |λ| ‖ℬ‖₁) ‖x‖) ≤ 10⁻¹⁰
cross-level:      classes matched one to one by minimum Σ|λ⁽³⁾ − λ⁽²⁾|,
                  |λ⁽³⁾ − λ⁽²⁾| ≤ 0.01 max(|λ⁽³⁾|, |λ⁽²⁾|) + s₂/ω₀
level stability:  |Re₂ − Re₃| ≤ δ_lev |Re₃|,  Re_h the largest genuine real part at level h,  δ_lev = 10⁻²
leading order:    λ = ε κ τ,  κ = 1 / (8 √λ₀ Q),  τ from the frozen spectra (× 9/16 in R5)
```

- **The types.** A point is **Elliptic** if `k_r = k_c = 0` and the count balances, and **Hyperbolic** if `k_r + k_c > 0` is level-stable. It is **Unresolved** if any check fails, or before a mode's pre-registered first-resolved point.
- **The cross-level check** runs at every third ladder point, at every bisection midpoint, and wherever an instability is a candidate.
- **Changes of type** are bisected to 2% in `ε`, between adjacent ladder points of different decided types, neither Unresolved.
- **Validity at `ε_min`.** Every mode resolved there must match its leading order within 5%, with its type and Krein sign exact.

### 1.7 The persistence test (level 2)

```text
where:        the ladder points nearest ε_max/8 and ε_max/2, each moved down to the last Elliptic point inside the
              slow lineage, above ε_min, with r T ≤ ln(1/(10η)) ≈ 4.61, r the largest real part among the level-2
              lifted modes
lineage:      S = G ⊕ the genuine slow modes, dim S = 2(j₂ + 1), gap |λ_next| ≥ 2 max over S∖G of |λ|,
              every principal cosine with the previous S ≥ 0.9; the first failure ends it; at ε_min it starts
              only if validity and the cross-level check pass
direction:    v = (Re g, Im g; 0),  g ~ N(0, 1)ⁿ from default_rng(20261006): the real parts' standard_normal(n),
              then the imaginary parts'
              slow half (P_S − P_G) v and fast half (1 − P_S) v, by Ω-orthogonal projection,
              each scaled to ‖·‖_E = (η/√2) ‖Φ‖_M,   ‖(ζ, ξ)‖_E² = ζᵀMζ + ξᵀMξ
initial data: Ψ(0) = Φ + δ,  Π(0) = c₀ (Φ + δ) + π,   c₀ = (e^{iθ} − 1 + dt²σ/2)/dt,  θ = 2 arcsin(dt √σ / 2)
reference:    Ψ_ref(n dt) = e^{inθ} Φ, which the scheme carries exactly
scheme:       Störmer–Verlet, consistent mass factorized once,  dt = 0.8 · 2/ω_max,  T = 1000,  η = 10⁻³
measures:     (z_p, z_v) = (Ψ, Π) aligned modulo the phase (Sp(1) in R2) and taken into Φ's rotating frame,
              ζ = z_p − Φ,  ξ = (z_v − c₀Φ) − c₀ζ
              d = ‖ζ‖_M/‖Φ‖_M,   d_⊥ = ‖static part of (1 − P_G)(ζ, ξ)‖_M/‖Φ‖_M,
              d_rot = ‖M-orthogonal projection of ζ onto span X_aΦ‖_M/‖Φ‖_M,   d_ref the same as d for the reference
verdict:      first event: Invalid (drift) if d_rot > √η before d_⊥ > 100η; Fails if d_⊥ > 100η first;
              a tie at one observation is Invalid (drift); otherwise Persists if max d_⊥ ≤ 10η, Unresolved if not
floors:       d_ref ≤ η;  charge drift ≤ 10⁻¹⁰ relative;  E_η = E(0) − E_ref > 0;
              max |E − E(0)| ≤ ¼ E_η;  secular drift ≤ 10⁻² E_η;  all read up to the first event, or T if none
```

`(δ, π)` are the direction's position and rotating-frame velocity parts. `ω_max²` is the largest eigenvalue of `K` relative to `M`. The system is observed every unit of time.

## 2. Equation-to-code map

The links are pinned to MIT `6058f95`, the frozen solver. The physics lives in short modules whose docstrings carry these equations:
- `fem.py` and `nonlinear.py`, the action;
- `pinned.py`, the Newton and its stopping test;
- `linearize.py` and `verdict.py`, the linear verdict;
- `persistence.py`, the dynamics.

`m815_branch.py` is the driver.

| Equation | Code | Line |
| --- | --- | --- |
| § 1.1 `λ̄₀`, the mean of the lowest discrete level | `lowest_level` | [`standing.py:83`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/standing.py#L83-L94) |
| § 1.2 the quadrature | `tet_quadrature` | [`fem.py:15`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/fem.py#L15-L37) |
| § 1.2 elements: exact-geometry P2, metric at the quadrature points | `Mesh._element_matrices` | [`fem.py:62`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/fem.py#L62-L86) |
| § 1.2 `K`, `M` with blocks `ρ(g_i)†ρ(g_j)` | `Mesh.assemble` | [`fem.py:88`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/fem.py#L88-L101) |
| § 1.2 `V₄` | `Quartic.energy` | [`nonlinear.py:20`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/nonlinear.py#L20-L23) |
| § 1.2 `f = ∇V₄` | `Quartic.force` | [`nonlinear.py:32`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/nonlinear.py#L32-L36) |
| § 1.2 `Df[δ] = Aδ + C conj(δ)` | `Quartic.jacobian` (`C` is `B` there) | [`nonlinear.py:38`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/nonlinear.py#L38-L59) |
| § 1.2 `H`, `N`, Störmer–Verlet | `Verlet` | [`integrate.py:9`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/integrate.py#L9-L35) |
| § 1.3 `F₀(x, σ)` | `Standing.residual` | [`standing.py:48`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/standing.py#L48-L49) |
| § 1.3 the gauges `p_k` | `fibre_generators` | [`standing.py:14`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/standing.py#L14-L21) |
| § 1.3 `Fix_H` and its projector | `PinnedCase.__init__` | [`pinned.py:139`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/pinned.py#L139-L154) |
| § 1.3 the slice `s` (T4) | `LevelCase.__init__` | [`m815_common.py:120`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_common.py#L120-L132) |
| § 1.3 the bordered Newton | `solve_pinned` | [`pinned.py:191`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/pinned.py#L191-L241) |
| § 1.3 the stopping test | `solve_pinned` | [`pinned.py:219`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/pinned.py#L219-L220) |
| § 1.4 the ladder | `ladder` | [`continuation.py:24`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/continuation.py#L24-L35) |
| § 1.4 the germ interpolant | `intertwiner`, `germ_seed` | [`standing.py:24`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/standing.py#L24-L37) |
| § 1.4 the seed | `Branch.first` | [`continuation.py:66`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/continuation.py#L66-L77) |
| § 1.4 secant and retries | `Branch.solve_target` | [`continuation.py:79`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/continuation.py#L79-L104) |
| § 1.4 the fold | `Branch.fold_search` | [`continuation.py:106`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/continuation.py#L106-L133) |
| § 1.4 the branch end | `Branch.run` | [`continuation.py:135`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/continuation.py#L135-L159) |
| § 1.5 `L` and the first-order pencil | `Linearization.__init__` | [`linearize.py:11`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/linearize.py#L11-L26) |
| § 1.5 slow eigenvalues | `Linearization.slow_eigs_q` | [`linearize.py:54`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/linearize.py#L54-L67) |
| § 1.5 `Ω` and `B` | `PhaseSpace` | [`verdict.py:29`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/verdict.py#L29-L52) |
| § 1.5 `n(L)` | `negative_index_L` | [`verdict.py:115`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/verdict.py#L115-L139) |
| § 1.5 the count identity | `Verdict.count` | [`verdict.py:292`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/verdict.py#L292-L302) |
| § 1.6 `G`, the 0.9 rule, the cluster checks, Krein values | `Verdict.__init__` | [`verdict.py:142`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/verdict.py#L142-L237) |
| § 1.6 the phase partner `∂_σΦ` | `d_sigma_phi` | [`verdict.py:91`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/verdict.py#L91-L112) |
| § 1.6 representatives | `representatives` | [`verdict.py:83`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/verdict.py#L83-L88) |
| § 1.6 the fast search | `Verdict.fast_search` | [`verdict.py:239`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/verdict.py#L239-L275) |
| § 1.6 classes | `classify` | [`m815_common.py:94`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_common.py#L94-L99) |
| § 1.6 leading order `εκτ` | `kappa`, `predicted` | [`m815_common.py:46`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_common.py#L46-L60) |
| § 1.6 cross-level matching | `match` | [`m815_common.py:75`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_common.py#L75-L91) |
| § 1.6 cross-level agreement | `cross_level` | [`m815_branch.py:123`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_branch.py#L123-L132) |
| § 1.6 the per-point record | `Point` | [`m815_branch.py:55`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_branch.py#L55-L120) |
| § 1.6 validity at `ε_min` | `run_branch` | [`m815_branch.py:210`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_branch.py#L210-L230) |
| § 1.6 the types | `run_branch.decide` | [`m815_branch.py:243`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_branch.py#L243-L276) |
| § 1.6 changes of type | `run_branch` | [`m815_branch.py:284`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_branch.py#L284-L301) |
| § 1.7 the slow lineage | `slow_lineage` | [`m815_branch.py:157`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_branch.py#L157-L179) |
| § 1.7 where the tests run | `run_branch` | [`m815_branch.py:308`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/m815_branch.py#L308-L342) |
| § 1.7 the direction | `direction` | [`persistence.py:40`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/persistence.py#L40-L58) |
| § 1.7 the drift margin | `margin_ok` | [`persistence.py:61`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/persistence.py#L61-L63) |
| § 1.7 alignment modulo the phase | `Aligner` | [`persistence.py:66`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/persistence.py#L66-L89) |
| § 1.7 the run and the measures | `run_test` | [`persistence.py:92`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/persistence.py#L92-L156) |
| § 1.7 the verdict and the floors | `evaluate` | [`persistence.py:159`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/persistence.py#L159-L187) |
| § 1.7 `dt` | `omega_max` | [`integrate.py:38`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/integrate.py#L38-L40) |
| § 1.2 the rotation derivatives `X_a` | `rotation_derivatives` | [`symmetry.py:123`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/symmetry.py#L123-L146) |

## 3. The run

- **Environment.** Python 3.13.13, numpy 2.5.0, scipy 1.18.0, macOS 15.8 on arm64, with one thread and fixed ARPACK starting vectors, as the manifest records. The solver is deterministic under that environment.
- **Controls first.** The closing controls were rerun from the frozen bytes and passed every gate.
- **`FROZEN.json`.** It was written at 2026-10-07 20:10:02 UTC from the committed frozen numbers, before any target computation. The driver refuses a target without it.
- **The targets** ran in waves sized for memory, with T4 first.
- **Deviations,** as `RECORD.md` § 4 lists them. The run log records 1 to 3 and the post-run plan:
  1. **The controls-first stop.** The run's own byte check on `freeze_numbers.json` stopped it before any target, because the regenerated `s₂` moved by at most `1.05×10⁻¹⁰` relative: the frozen spreads predate the one-thread pin. It resumed on the author's go, with the frozen values governing.
  2. **All nine targets launched at once.** An orchestration fault started them together. The seven out-of-wave processes were stopped after about 50 s in setup, with empty logs, and rerun from the start in their waves. The orchestration lies outside the frozen set.
  3. **The driver's IndexError on an empty ladder, in T4/R4 and T4/R5.** It is recorded under the contingency in the author's plan of record, fixed before launch: the branch ends there, by an instrument failure, with no code change. The terms have no crash rule.
  4. **Post-run diagnostics, outside the record:** the exact replays, the post-hoc Newton trace, the gauge-term check, and an independent review's bench checks. None changed a record.
  5. **The controls-first rerun wrote over the frozen control records.** The frozen ones are restored from `6058f95`, and the rerun's are filed in `out/controls_first_2026-10-07/`.
  6. **A pre-flight smoke test** on a control (C2b at bench test levels 2 and 1, to `ε = 2`, `T = 10`), after the go and before launch. Its output was moved out of the repo.

## 4. Results against the frozen terms

The record is per branch. [`RECORD.md`](https://github.com/dmobius3/mode-identity-theory/blob/28c8a7648a88c9234d32be753cbbe3563ea2eec3/files/framework/files/working/files/scripts/soft-slot-dynamics/RECORD.md) gives every number, with its gate.

| Branch | Linear types | Validity at `ε_min` | Persistence, max `d_⊥` (bar `10⁻²`) |
| --- | --- | --- | --- |
| T1/R4, T1/R5 | Unresolved at every point: cluster check, and a paired lifted mode above `τ₀` | **INVALID** (§ 5.1) | Not reached, both |
| T2/R4 | Elliptic 0.25 to 11.31; Unresolved at 16 and 22.63; **Hyperbolic at 32 to 72** | pass, 1.23% | Persists at 5.657 (both tests moved there), `1.34×10⁻³` |
| T2/R5 | Elliptic at all 15 | pass, 0.29% | Persists at 4 and 16, `1.00×10⁻³` and `1.07×10⁻³` |
| T3/R4 | Elliptic at 17; Invalid at 32, an eigenpair residual (§ 5.4) | pass, 0.38% | Persists at 2.828 (both tests moved there), `1.00×10⁻³` |
| T3/R5 | Elliptic 0.25 to 11.31; Unresolved at 16; **Hyperbolic at 22.63 and 32** | pass, 1.28% | Persists at 4, and at 11.31 (moved from 16), `1.07×10⁻³` and `1.45×10⁻³` |
| T5/R2 | Elliptic at 17; Unresolved at 45.25, 128 and 132 | pass, 0.48% | Persists at 16 and 64, `0.98×10⁻³` and `0.94×10⁻³` |
| T4/R4, T4/R5 | none: the branch ended at the seed (§ 5.2, § 5.3) | none | none |

- **No bisection ran.** Every step from Elliptic to Hyperbolic passes through an Unresolved point, so the decided types bracket each change: 11.31 to 32 in T2/R4, and 11.31 to 22.63 in T3/R5.
- **The drift margin also bound at T2/R4 and T3/R4.** The terms' expected-coverage note (lines 255 to 258) named T4 as where real pairs may still appear.
- **Coverage of the slowest slow mode.** Every persistence test that ran spans more than 1.5 periods of it.

## 5. Instrument defects found after the run

None is repaired under M8.15. The terms allow a repair only while no target has run (line 270), and any change after the go is a new pre-registration (line 280).

### 5.1 T1: the overlap rule assigns to `G` the pair nearest the missing resolved mode

- **The missing mode.** At `ε_min`, T1's resolved Krein-negative slow mode is predicted at `|λ| = 1.909×10⁻⁴` (`R4`) and `1.087×10⁻⁴` (`R5`). It is missing from the genuine modes.
- **The slow set at `ε_min`.** It holds 2 fibre modes, 4 lifted modes and 8 genuine ones. The lifted modes are the rotation pair, at `1.03×10⁻⁵` (`R4`) and `8.74×10⁻⁶` (`R5`), and an extra pair, so `G` has dimension 6 against the frozen 4.
- **The evidence points to where the mode went.** The extra pair sits at `1.814×10⁻⁴` (`R4`) and `1.145×10⁻⁴` (`R5`), within 5.0% and 5.3% of the prediction. No genuine mode lies near it, and the nearest other predicted mode is 18 times larger. Its static part has `M`-overlap at least 0.9 with `U`, so the § 1.6 rule puts it in `G`. An independent review's bench check, not declared in advance, puts that mode's static overlap with `U` at 0.968 in the small-amplitude limit, and 0.986 for the lifted quartet it joins at level 2 in `R4`, through the frozen verdict.
- **A second reason in `R4`.** The rotation pair's own lift is `0.986 τ₀` at `ε_min` and exceeds `τ₀` from `ε = 1` to 72, rising to `2.70 τ₀`. So "a paired lifted mode above `τ₀`" holds at 14 of `R4`'s 18 points even without the extra pair. In `R5` it stays between `0.47 τ₀` and `0.86 τ₀`.
- **The verdict.** Validity at `ε_min` fails. The terms make this an instrument defect, never a counterexample to the theorems (line 274).

### 5.2 T4: the stopping test leaves out a nonzero multiplier

**The mechanism.** At a solution of the bordered system of § 1.3:
- `W_Hᵀ F = 0`.
- The phase direction `p = iΦ` lies in `Fix_H`, because `T4`'s fixed-space condition is complex-linear.
- `pᵀ F₀(x, σ) = 0` by phase invariance.
- `F` has no part outside `Fix_H` beyond round-off (measured below), so `F = 0`: `F₀ + ν M s = −μ M p₀`, and the norm term `|xᵀMx − m| / m` vanishes.

So

```text
μ = −ν ⟨p, s⟩_M / ⟨p, p₀⟩_M,      and the stopping test reads   |μ| ‖M p₀‖ / ‖M x‖
```

This is nonzero whenever the slice overlaps the phase direction and carries force. T4's slice has a nonzero phase overlap, and the post-run solve finds a nonzero slice force; the terms permit that force and require it to be recorded (line 98).

**The post-run evidence** ([`out/post_run/`](https://github.com/dmobius3/mode-identity-theory/blob/28c8a7648a88c9234d32be753cbbe3563ea2eec3/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/out/post_run/README.md)), outside the record, from a check declared in the bench log before it ran:
- **The equation is solved.** With every multiplier included, the residual is `3.5×10⁻¹³` (`R4`, level 2), `3.8×10⁻¹³` (`R5`, level 2) and `1.9×10⁻¹²` (`R4`, level 3). The part outside `Fix_H` is at the same level.
- **The multiplier matches.** `μ` equals its prediction to the printed digits.
- **Where the test stalls:**
  - `R4` at level 3: `|μ| ‖M p₀‖ / ‖M x‖ = 1.085×10⁻⁶`, with `M`-cosine 0.1256 between the slice and the phase direction;
  - `R5` at level 2: `1.654×10⁻⁵`;
  - the frozen Newton trace of `R4` holds at `1.08×10⁻⁶` from iteration 3 to 39.
- **The control.** The slice control C3s has a zero cosine (`−5.5×10⁻¹⁶`), and its test passes.

The test is at [`pinned.py:219`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/pinned.py#L219-L220). Under line 44, a defect claim is reproduced by the maintainer with independent code before it is ratified.

### 5.3 The driver on an empty ladder

When the level-3 ladder reaches no point, the driver hands the empty list to the level-2 phase, and `Branch.run` raises IndexError at `targets[0]` before the branch record is written. T4's two records are reconstructed from the run log, with the end reason recovered by an exact replay from the frozen bytes. They carry no ladder point, verdict or persistence entry.

### 5.4 T3/R4 at `ε = 32`: a label

An eigenpair's backward error exceeds the `10⁻¹⁰` bound there. The terms make a point that breaks a residual bound invalid (lines 160 and 273), while the frozen code writes Unresolved. Line 103 also keeps a slow eigenpair only below the bound; the code keeps it and marks the point. Neither label is scored, and the persistence placement is the same under either.

## 6. Provenance

- **Terms:** `tasks/m8_15_task_details.md` at OpenWave `4b1177c8`, SHA-256 `e0c32d01…`.
- **Solver:** MIT `6058f95`. Its manifest pins 23 files, and `SHA256SUMS` 52.
- **Record:** MIT [`28c8a76`](https://github.com/dmobius3/mode-identity-theory/tree/28c8a7648a88c9234d32be753cbbe3563ea2eec3). It holds the branch records, the run logs, `FROZEN.json`, the persistence records, T4's reconstructed records, and the post-run diagnostics in `out/post_run/`.
- **Reproduction:** `REPRODUCE.md` beside the solver.

## 7. What this run does not verify

- **T1 and T4.** There is no finite-amplitude verdict for either, in any slot. A fix to the overlap rule alone would leave T1/R4 Unresolved from `ε = 1` under this floor.
- **The changes of type.** In T2/R4 and T3/R5 they are bracketed, not located.
- **Persistence.**
  - It holds over `T = 1000`, at the points where each test ran, for the one registered direction.
  - No test ran at a Hyperbolic point.
  - Three test runs, covering five scheduled tests, ran below their scheduled points.
- **Levels.** The linear record is at level 3 and the dynamics at level 2. No level-4 check was run.
- **T4/R5's post-run diagnosis** is at level 2 only.

## 8. Audit record

Two independent agents reviewed this note against the frozen terms, the frozen code and the record, before it was filed. On the first draft, each re-walked the equation-to-code map and tried to refute § 5's four claims; neither could, and their corrections are folded in. On the second draft, both returned GO, with four optional wording points, also folded in. One agent's own level-2 run reproduces the numbers of § 5.2 digit for digit.
