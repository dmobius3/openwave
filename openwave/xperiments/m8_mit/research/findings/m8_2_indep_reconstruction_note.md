# M8.2 method note: independent reconstruction of the first-occurrence table

> The review-time verification behind the M8.2 close-out claim "first-occurrence table
> independently reproduced 9/9 by the maintainer", landed as a repository artifact at the
> author's request ([PR #350](https://github.com/openwave-labs/openwave/pull/350) close-out,
> 2026-07-28) so the M8.5 reproduction protocol can quarantine it by name. Object under
> comparison: [`m8_2_preregistration.md § 6.1`](m8_2_preregistration.md). Script:
> [`../scripts/m8_2_indep_reconstruction.py`](../scripts/m8_2_indep_reconstruction.py).
> Claim label: **independent-method reproduction, NOT blind**
> ([roadmap § CONVENTIONS](../m8_roadmap.md#conventions)).

## 0. What this is, and the three objects it sits between

The M8.2 lock (§ 3) requires every certification table to be reproduced through an
independently implemented decomposition. The author added a **context firewall** to that
requirement in the close-out thread, and the M8.5 protocol therefore distinguishes three
objects rather than two:

| # | Object | Method | Status |
| --- | --- | --- | --- |
| 1 | [`../scripts/m8_2_first_occurrence.py`](../scripts/m8_2_first_occurrence.py) | McKay / affine-E8 recursion `V_{n+1} = V_1 ⊗ V_n − V_{n−1}`, hardcoded 9-class list, irrep labels + dims + distances as literals | in repo (M8.2) |
| 2 | [`../scripts/m8_2_indep_reconstruction.py`](../scripts/m8_2_indep_reconstruction.py) | explicit 120 unit quaternions → conjugacy classes by brute-force conjugation → characters by Burnside class-sum diagonalization; dims, adjacency, distances all DERIVED | this note |
| 3 | the M8.5 deliverable-A implementation | the author's frozen protocol, implemented in a fresh context with no M8.2 internals loaded | pending |

Objects 1 and 2 are **forbidden to object 3's implementer** until 3's own source and raw
output are committed; afterwards they serve as adjudication references only. Agreement of
all three is reported as **three-way agreement**, which strengthens provenance without
raising the claim label: object 2 was written by a context that had already seen the target
table, so "blind" is spent and stays spent.

## 1. Equations first

### 1.1 The group, built explicitly

2I ⊂ ℍ, the binary icosahedral group, as 120 unit quaternions (icosians), with
φ = (1+√5)/2:

```text
8    ± 1, ± i, ± j, ± k
16   (± 1 ± i ± j ± k) / 2
96   even permutations of (0, ± 1, ± 1/φ, ± φ) / 2
```

Multiplication is the quaternion product; nothing about the group is assumed beyond this
list, and closure is checked over all 120 × 120 products.

### 1.2 Conjugacy classes and eigen-angles

```text
C(x) = { g x g⁻¹ : g ∈ 2I }                    (orbits, computed by brute force)
q = cos θ + sin θ n̂ ,  |n̂| = 1   ⇒   χ_Q(q) = tr_SU(2) q = 2 cos θ
```

θ is read off the real part of one representative per class, so each class carries an
eigen-angle without any table of angles being supplied.

### 1.3 Characters by Burnside class-sums

Let `Ĉ_i` be the class sums in the center of the group algebra, with structure constants

```text
Ĉ_i Ĉ_j = Σ_k a_ijk Ĉ_k ,        a_ijk = #{(a,b) ∈ C_i × C_j : ab = fixed z ∈ C_k}
```

Central characters `ω_i = |C_i| χ(g_i) / χ(1)` satisfy `ω_i ω_j = Σ_k a_ijk ω_k`, so each
irreducible character is a simultaneous eigenvector of the commuting matrices
`(M_i)_{jk} = a_ijk`. A generic real combination `Σ_i r_i M_i` separates them in one
eigen-decomposition. Normalizing the identity class to `ω_0 = 1`, the degree and the
character follow from column orthogonality:

```text
d = √( |G| / Σ_i ω_i² / |C_i| ) ,        χ_i = d ω_i / |C_i|
```

Nothing here uses the McKay recursion, the E8 diagram, or any irrep labelling.

### 1.4 The three flat connections (a contract input, not a result)

The pre-registration fixes the coefficient bundle as `τ_σ = Sym²(σ)` for
σ ∈ {trivial, Q, Q'}, where Q is the defining 2-dim rep and Q' its Galois partner:

```text
χ_Sym²V (g) = ½ ( χ_V(g)² + χ_V(g²) )
```

Q is *identified*, not declared: among the derived 2-dim irreps, Q is the one whose
character equals `2 cos θ` on every class; Q' is the other.

### 1.5 The McKay graph, derived

```text
V_1 ⊗ R_i = Σ_j A_ij R_j        (A = adjacency; multiplicities from character inner products)
d(R_i)     = graph distance from the trivial node, by BFS on A
```

The affine-E8 mark condition `A · dims = 2 · dims` is then a falsifiable consequence, not
an input.

### 1.6 First occurrence

For the SU(2) tower `V_n` (spin n/2), taken as the explicit weight sum rather than the
closed form:

```text
χ_{V_n}(θ) = Σ_{k=0}^{n} e^{i(n−2k)θ}          ( = sin((n+1)θ) / sin θ )

mult(ρ ; n, σ) = (1/|G|) Σ_i |C_i| χ_{V_n}(θ_i) χ_{τ_σ}(θ_i) χ_ρ(θ_i)

n_first(ρ, σ) = min { n : mult(ρ ; n, σ) > 0 }      eigenvalue n(n+2)/R²,  j = n/2
```

Rows are matched to the published table **label-free**, by the `(dim, distance)`
signature, because the irrep names are an input to object 1 and an output here. That the
9 signatures are distinct is itself checked (C7), since the matching is ill-posed
otherwise.

## 2. Equation-to-code map

| Equation | Function / block | Permalink |
| --- | --- | --- |
| § 1.1 the 120 icosians | `build_2i` | [L146-L172](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L146-L172) |
| § 1.1 quaternion product, conjugate | `qmul`, `qconj` | [L127-L143](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L127-L143) |
| § 1.2 conjugacy classes | C2 block | [L229-L241](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L229-L241) |
| § 1.2 eigen-angles θ_i | `thetas` | [L249-L251](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L249-L251) |
| § 1.3 structure constants a_ijk | `cmat` | [L254-L260](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L254-L260) |
| § 1.3 diagonalization, degree, χ | `mix`, `omega`, `deg` | [L262-L269](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L262-L269) |
| § 1.3 orthonormality (C3) | `gram`, `ortho_err` | [L276-L283](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L276-L283) |
| § 1.4 identify Q vs Q' (C4) | C4 block | [L285-L294](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L285-L294) |
| § 1.4 χ_Sym² | `chi_of_square`, `sym2` | [L297-L308](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L297-L308) |
| § 1.5 adjacency + E8 marks (C6) | C6 block | [L326-L334](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L326-L334) |
| § 1.5 BFS distances, signatures (C7) | C7 block | [L336-L350](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L336-L350) |
| § 1.6 χ_{V_n} weight sum | `chi_su2` | [L197-L202](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L197-L202) |
| § 1.6 multiplicity, n_first | `first_level` | [L353-L358](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L353-L358) |
| § 3 the comparison (C8) | C8 block | [L360-L379](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L360-L379) |
| the transcribed target | `DOC` | [L85-L95](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/scripts/m8_2_indep_reconstruction.py#L85-L95) |

Single file, 492 lines, NumPy only. There is no driver and no import from the M8.2 script:
`python3 m8_2_indep_reconstruction.py` runs the whole reconstruction.

## 3. Result

Raw output of record: [`../data/m8_2_indep_reconstruction_raw.txt`](../data/m8_2_indep_reconstruction_raw.txt)
(reconstruction + mutation suite), structured form
[`../data/m8_2_indep_reconstruction.json`](../data/m8_2_indep_reconstruction.json).

| Item | Value |
| --- | --- |
| Repository commit verified against | `d6790eef13398b809ed1a13486ef179089eba9e7` |
| Environment | python 3.13.9, numpy 2.3.5, macOS-26.5.2-arm64, seed 11 (fixed) |
| Derived class sizes | `[1, 1, 12, 12, 12, 12, 20, 20, 30]` |
| Derived irrep dims | `[1, 2, 2, 3, 3, 4, 4, 5, 6]`, Σd² = 120 |
| Orthonormality error | 2.36e-14 |
| \|χ_Q − 2cos θ\| | 1.78e-15 |
| McKay graph | 8 edges, `A·dims = 2·dims` holds |
| **Rows reproduced** | **9 / 9** |

The reconstructed table, all nine rows identical to § 6.1:

```text
 dim  dist   trivial  standard    galois
   1     0         0         2         6
   2     1         1         1         5
   3     2         2         0         4
   4     3         3         1         3
   5     4         4         2         2
   6     5         5         3         1
   3     6         6         4         0
   4     6         6         4         2
   2     7         7         5         3
```

## 4. Every PASS line is mutation-tested

Standing rule from this task's close-out ([roadmap § CONVENTIONS](../m8_roadmap.md#conventions),
"Self-checks must be able to fail"): a check whose two sides evaluate the same expression
reads exactly like a verified result. `python3 m8_2_indep_reconstruction.py --mutation-tests`
injects one deliberate defect at a time and requires the targeted check to go red. All 8
checks are covered; the suite exits non-zero if any mutation fails to redden its target or
any check has no mutation.

| Check | Mutation that reddens it | Also reddens |
| --- | --- | --- |
| C1 group order + closure | `phi_wrong`: φ → 1.6 | |
| C2 9 classes, sizes divide 120 | `conj_without_inverse`: `g x g` for `g x g⁻¹` | |
| C3 orthonormality | `perturb_char`: one character row × 1.05 | |
| C4 χ_Q = 2 cos θ | `q_is_qprime`: Q' picked as defining | C8 |
| C5 Sym²Q, Sym²Q' irreducible dim 3, distinct | `sym2_as_square`: χ² for Sym²χ | C8 |
| C6 E8 mark condition | `edge_threshold`: edges only above multiplicity 1.5 | C7 |
| C7 distances defined, signatures distinct | `bfs_forget_increment`: parent's distance, not +1 | |
| C8 all 9 rows | `chiv_offbyone`: n weights instead of n+1; `doc_typo`: one target entry altered | |

## 5. Scope honesty: what this does NOT verify

| Not verified | Why it is out of scope |
| --- | --- |
| The **coexact 1-form entry rule** (level 2 at d=0, 3 at d=1, else d) | Nothing here derives it and it has no published target. It stays **ASSERTED**, exactly as the M8.2 script labels it. This artifact reconstructs the scalar (0-form) table only |
| The author's `mass-spectrum.md § 4` table | Object 1 cross-checks against it; this note compares against § 6.1. The two agreeing is object 1's own cross-check, not an independent result here |
| Which model each table is a fixture for (M4_int, the prospective `M7_ad`) | A contract question settled in the pre-registration, not a computation |
| The three flat connections themselves | `τ_σ = Sym²(σ)` is a contract input (§ 1.4); the reconstruction verifies the consequences of that choice, not the choice |
| Anything dynamical | No quotient operator, no spectrum on S³/2I, no Lagrangian. That is M8.4 / M8.5 work |
| The three source PDF hashes in § 1 of the lock | The files live in the gitignored `theory/` folder by design |

## 6. Audit record

| Layer | Status |
| --- | --- |
| Falsifiability of every printed PASS | ✅ mutation suite, § 4, 8/8 covered, exits non-zero on regression |
| Method independence from object 1 | ✅ by construction, § 0 and § 1; no import, no fixture, label-free row matching |
| Context independence | ❌ not claimed. This context had seen the target table, which is why the claim ceiling is independent-method and not blind |
| Independent adversarial audit by a second agent | 🔶 not run on this note. The scheduled adversarial pass is M8.5 deliverable A: a fresh-context implementation of the author's frozen protocol, which is a stronger test of the same claim than a re-read of this page |
