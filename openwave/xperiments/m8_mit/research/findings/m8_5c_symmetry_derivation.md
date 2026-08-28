# The residual symmetry of the M4L_Erho problem: derivation note for protocol § 2

> Status: DERIVED. Every claim here is mathematics, closed without numerics, per the M8.5-A § 6
> standing rule that a derivation moves the claim and numbers contribute nothing to it.
> `design_inputs/right_translation_check.py` arms the IMPLEMENTATION bridges only, one labeled
> check per numbered step, so a red tells you which theorem-to-code bridge broke, never whether
> the theorem holds. Terminology: the ACTION leaves the law invariant; the nonlinear map itself
> transforms EQUIVARIANTLY.

## Setting

`S³ = SU(2)` with the round bi-invariant metric; `Γ = 2I ⊂ SU(2)` the binary icosahedral group
(the center `−1 ∈ Γ`). The arena is the left quotient `X = Γ\S³`. For each irrep `ρ` of `Γ` on
`W_ρ`, the flat bundle is the associated bundle

    E_ρ = Γ\(SU(2) × W_ρ),   [x, v] ~ [γx, ρ(γ)v],

whose sections are functions `f: SU(2) → W_ρ` with `f(γx) = ρ(γ)f(x)`.

## The theorem chain

**1. The right action is well defined on the bundle.** `[x, v]·g := [xg, v]` respects the
identification: `[γx, ρ(γ)v]·g = [γxg, ρ(γ)v] = [xg, v]·g`, because left multiplication by `γ`
and right multiplication by `g` commute.

**2. Section equivariance is preserved.** For `f(γx) = ρ(γ)f(x)` and `(R_g f)(x) := f(xg)`:

    R_g f(γx) = f(γxg) = ρ(γ) f(xg) = ρ(γ) R_g f(x),

so `R_g` maps sections of `E_ρ` to sections of `E_ρ`, for EVERY `g ∈ SU(2)`.

**3. Laplacian commutation.** The round metric is bi-invariant, so right translations are
isometries and `[Δ, R_g] = 0`.

**4. Cubic equivariance.** The fibre metric is `Γ`-invariant and `|R_g ψ(x)|² = |ψ(xg)|²`, so

    N(R_g ψ) = |R_g ψ|² R_g ψ = R_g(|ψ|² ψ) = R_g N(ψ):

the map `N` is equivariant; the law `ψ̈ = c²Δψ − c1|ψ|²ψ` is invariant.

**5. Spectral module structure.** Peter-Weyl decomposes `L²(SU(2))` into `⊕_n V_n ⊗ V_n^*`
with the LEFT action on `V_n` and the RIGHT action on `V_n^*`; the two commute. Taking
`Γ`-equivariant parts against `W_ρ` (Frobenius reciprocity) gives, per level `n`,

    H_{ρ,n} ≅ Hom_Γ(V_n, W_ρ) ⊗ V_n^*,

a right-`SU(2)` module whose right action lives entirely on the `V_n^*` factor. The mode-count
factor `(n+1)` in the pinned tables IS this multiplet. (Convention footnote: the equivariant
part can be written against `V_n` or `V_n^*`; 2I has a real character table, so its irreps are
self-dual and the two readings agree.)

**6. Multiplicity one at the scored level, TWO sources.** When
`dim Hom_Γ(V_l, W_ρ) = 1`, step 5 gives `H_{ρ,l} ≅ V_l^*`, IRREDUCIBLE under `SU(2)_right`.
The premise holds at every scored level, from two distinct pinned sources:

- the eight nontrivial sectors at `l = d_ρ`: the McKay first-occurrence tables (first
  occurrence at `n = d_ρ` with multiplicity exactly one, all eight sectors);
- `R_0` at `l = 12` (Control B's space, NOT a first occurrence; `R_0`'s first occurrence is
  `n = 0`, the constants): the invariant ring of `2I`, Molien series
  `(1 + t³⁰)/((1 − t¹²)(1 − t²⁰))`, giving invariant levels through 36 at
  {0, 12, 20, 24, 30, 32, 36}, each with multiplicity one. Hence `H_{R0,12} ≅ V_12^*`,
  complex dimension 13, spin 6.

**7. Galerkin inheritance.** The retained space at cutoff `N` is the sum of COMPLETE levels
`n ≤ N`, each level right-invariant by step 5, so `P_N^ρ` commutes with `SU(2)_right`; with
step 4, the truncated nonlinear system retains the full symmetry

    G_ρ = (U(1) × SU(2)_right)/K_ρ,

where `U(1)` is the phase and `K_ρ` is the finite central kernel: `−1` is central in `2I` and
`ρ` is irreducible, so Schur gives `ρ(−1) = ±I`; then `R_{−1}f(x) = f(x(−1)) = f((−1)x) =
ρ(−1)f(x) = ±f(x)`, the pair `(ρ(−1)⁻¹ phase, −1)` acts trivially, and `K_ρ ≅ Z/2` for every
`ρ`. `dim G_ρ = 4` always.

**8. Corollaries, not new assumptions.**

- **Orbit tangent.** At a state `φ`, `Z(φ) = span_R{iφ, T_1φ, T_2φ, T_3φ}` is the tangent
  space to the `G_ρ`-orbit; `rank_R Z(φ)` is the orbit dimension. This is D-7's predicted
  Jacobian nullity at regular points.
- **Achievable ranks.** The isotropy algebra of `φ ≠ 0` inside `u(1) ⊕ su(2)` has dimension 0
  or 1: `su(2)` has no two-dimensional subalgebra, and no nonzero vector of an irreducible
  `V_l`, `l ≥ 1`, is fixed by all of `su(2)`. Hence rank 1 at `l = 0` (phase only); rank 3
  ALWAYS at `l = 1` (every nonzero vector of `C²` is a weight vector of some torus); ranks
  3 or 4 at `l ≥ 2`, by continuous-isotropy stratum.
- **Branch search space.** D-6's enumeration is the isotropy lattice of `G_ρ` on `V_l^*`,
  pinned for `l = 0..7` and `l = 12`; EBL existence claims only where the fixed-point subspace
  is COMPLEX one-dimensional (one phase orbit), phase-quotiented variant.
- **Noether content for gate 5.** Time translation gives energy; `U(1)` gives charge; the
  three right generators give three momenta. Steps 4 and 7 make all of them EXACT invariants of
  the semi-discrete Galerkin flow. Charge and the three momenta are momentum maps of LINEAR
  symmetries of the discrete Lagrangian, so the variational integrator conserves their discrete
  versions at rounding (discrete Noether); energy carries the standard leapfrog oscillatory
  envelope.

## What the script arms, per step

Two layers, kept distinct: a green below means the REPO SCALAR PRIMITIVES realize the derived
right action. The `W_ρ`-valued sector bases and the real Galerkin system do not exist yet and
receive no credit here; protocol gates 3 and 5 qualify them when they are built.

| check | steps armed | red means |
| --- | --- | --- |
| C1 right-action realization, `l = 1..7, 12` | 1, 2 | the code's `π_l` is not a homomorphism |
| C2 assembled vs analytic coefficient rep | 7 | the quadrature-assembled action disagrees with the derived `I ⊗ P` convention (kron-swap mutation) |
| C3 cubic equivariance | 4 | the projected nonlinearity is not equivariant; the mutation runs under the `6N` rule so a red is the symmetry break, never aliasing |
| C4 complete-level necessity | 7 | complete levels are green; a dropped column leaks under translation (the mutation IS the content) |
| C5 multiplicity, both halves of step 6 | 6 | Molien levels by character + Reynolds; the eight first occurrences DISCOVERED over a fixed scan and compared to the pinned distances; any scored level not multiplicity one |
| C6 2I-commutant census | 6 | `dim = ⟨χ_l, χ_l⟩ = {1,1,1,1,1,2,2,4}` at `l = 1..7, 12`, commutative iff multiplicity-free; the `π₁⊕π₁` mutation fires on noncommutativity |
| C7 level-diagonal spectrum vs assembled right action | 3 | the analytic `n(n+2)` diagonal fails to commute with the ASSEMBLED action (the assembled object is the action, not the Laplacian); intra-level breaking is detected |

Preflights P0 (group forensics) and P1 (the per-element character table, R5 derived in-room
by column orthogonality with a unique sign solution) are armed and
reported like the rest. Each check carries its own real green parent and its own real
mutation; there is no omnibus
PASS, and the script prints the resolved module path, its SHA-256, and library versions so a
green names the implementation that produced it. The first C6 (an SU(2)-commutant of `Sym^n`)
was replaced: Schur makes that 1 for every `n`, so it could not fail.
