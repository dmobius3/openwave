# AUDIT_STAGE2: comparison with the other worker, reproduction, refutation attempts, and grading of THEOREM.md

My stage-1 files (`METHOD.md`, `AUDIT_STAGE1.md`, `audit_results.json`, the stage-1 scripts and `out/`) were not edited.
Everything new is in new files: the `st2_*.py` scripts, `out2/`, `solver_rerun/`, `mut/`, `solver_rerun_run_all.log`,
`AUDIT_STAGE2.md` and `audit_stage2.json`. Nothing in `solver_work/` was edited.

**One incident to disclose.** While inspecting a data structure, I ran `./py -c "import s01_group"`. That re-executed my
stage-1 script `s01_group.py`, which rewrote `out/group.pkl`, `out/item01.json` and `out/checks_s01.json`. All three
are byte-identical (same SHA-256) to the copies from the stage-1 clean run in `clean_run/out/`. Only their timestamps
changed.

## Summary

| Task | Result |
|---|---|
| 1. Compare | 490 values paired and compared. All 490 agree exactly. Items 0–8 and 10–12 are **CONFIRMED**. Item 9 is **PARTIAL**: the zero sets are identical, but one reason is incomplete and its reason-checker cannot fail. |
| 2. Reproduce | **Yes.** The rerun printed 503 PASS and 0 FAIL. `results.json` and all 11 `results/*.json` files are byte-identical to the originals, and the log matches line by line apart from timings. |
| 3. Refute | No value was refuted. There are four kinds of PASS lines that cannot fail, plus two `results.json` fields that are hard-coded literals. Details below, each with its mutation. |
| 4. Arguments | Items 8, 9 and 10 are correct. Its item-8 sign argument is stronger than mine. Its item-9 "by construction" reason is weaker than it reads. Its item 10 is complete, and its transversality step is more explicit than mine. |
| 5. THEOREM.md | T.a–T.e, E1, E2, the bridge and the S lemma are **ESTABLISHED**. The λ₄ lemma is a **DEFECT**: its sign step omits a necessary hypothesis. |

---

## Task 1: Compare

Script: `st2_compare.py`, output `out2/compare.json` (one record per compared value, with both values and the rule).
- **Rule for exact values:** parse both strings with sympy and require `simplify(radsimp(mine − theirs)) == 0`. Where the
  other worker reports a norm and I report a squared norm, its value is squared first.
- **Integers and booleans:** `==`.
- **Residuals:** each side against its own declared threshold.

For item 3, `st2_item3.py` also applied its closed-form stabilizer generators and characters using **my** D-matrices.

| Item | Mine | Its value | Comparison run | Agree | Label |
|---|---|---|---|---|---|
| 0 | ‖q₁‖²=‖q₂‖²=1; \|Γ\|=120; perfect (derived subgroup 120); ⟨33;3−3\|60⟩=√231/462 | same | 6 values; == and exact equality | 6/6 | **CONFIRMED** |
| 1 | dim Hom table (both sectors, n=0..18); DERIVED; residuals < 1e-40 | same table; DERIVED; residuals ≤ 1.7e-49 | 38 dims ==, declaration ==, residuals against the common 1e-40 | 41/41 | **CONFIRMED** |
| 2 | Q (12 values), r̂₆ (6), all stationary, Π₆N = QΦ everywhere | same | exact equality on Q and r̂₆; == on flags | 36/36 | **CONFIRMED** |
| 3 | stabilizers T, N(T), 2O(48), 2D(24), C₁₀, 2D(12); characters; block χ-subspace dims 1,1,1,1,2,2 | same groups; generators Rz(π/2), R((1,1,0),π/2) / Rz(π/3), Rx(π) / Rz(2π/5) / Rz(2π/3), Rx(π); same characters | dims (6), χ-subspace dims in V_J for J=0..9 (60), orders (6). With my D: its generators fix [u] with its χ (residual ≤ 7.8e-16); they generate groups of order 48, 24, 10, 12; χ is a homomorphism on each whole group (≤ 2.6e-15). Its U2 coset representative Rx(π) and my ry(π) both act by −1 | 72/72, plus the cross-check | **CONFIRMED** (floating-point level on both sides) |
| 4 | ‖Π_nξ‖²/g² at every level, both sectors; same character as Φ | same; it also lists level 6 = 0 | exact equality on 60 ξ-norms (its level-6 entries against 0) and on all ‖Π_nN‖²; == on flags | 144/144 | **CONFIRMED** |
| 5 | along Φ; ⊥ norm; e_t / i·e_t, τ_x / τ_y components | same | exact equality (norms squared) | 48/48 | **CONFIRMED** |
| 6 | forms −56/165, −21/110, 6440/18447, 56/429, 2415/12298, 21/286; zeros; r̂₆″ −104/55, 920/473, 8/11 | same | exact equality | 15/15 | **CONFIRMED** |
| 7 | ‖κ‖², components, remainders, free dims 2/3, ⊥ equation with and without κ | same; its with-κ value is "0 (residual ≤ 1.7e-123)" | exact equality, except the with-κ zero: my exact 0 against its numerical zero | 60/60 | **CONFIRMED** |
| 8 | λ₄/g² with and without κ (12 values, equal); λ₄ < 0 | same | exact equality; both "negative" | 36/36 | **CONFIRMED** |
| 9 | computed zeros in items 4–8, all with reasons, none unresolved | same zero set; also lists Im Q, Im⟨Φ,DNξ⟩ and ‖Π₆ξ‖² | set of zeros mapped to a common vocabulary: both set differences empty; unresolved 0 = 0 | 2/2 | **PARTIAL** (see tasks 3 and 4) |
| 10 | existence yes for U1–U6, both sectors; orbit dims 2,2,3,3,3,3 | same; ranks 10,10,9,9,9,9 | == | 20/20 | **CONFIRMED** |
| 11 | not stationary; r̂₆ = 3061/14784; 805√3/6864 and 2415√3/36608 | same; also \|grad r̂₆\|² = 13225/92928 | exact equality. My gradient components 115/352 and −115√3/1056 have squared sum 13225/92928 | 8/8 | **CONFIRMED** |
| 12 | 8.5e-15 < 1e-10 (float64) | 1.5e-29 < 1e-25 (30 digits) | each below its own requirement; levels 20 and 22 empty on both sides | 2/2 | **CONFIRMED** |

Item 3 depends on floating-point stabilizer searches on both sides (its search ran at 30 digits with tolerance 1e-20).
In task 5 I certified exactly the parts that THEOREM.md needs (`st2_exact_m.py`).

## Task 2: Reproduce

I copied the `.py` files from `solver_work/` and `./py` into `solver_rerun/` (without `results/`, `results.json` or the
log) and ran `./py solver_rerun/run_all.py`. The log is `solver_rerun_run_all.log`.
- The run printed **TOTAL: 503 PASS, 0 FAIL**, matching its RETURN and its `run_all.log`.
- A recursive JSON diff of `solver_rerun/results.json` against `solver_work/results.json` finds **0 differences**.
- All **11** per-script `results/*.json` files are **byte-identical**.
- The log is identical line by line (964 lines) after masking the timings.

**It reproduces its own `results.json`.** Its run_all launches children with `sys.executable -S` instead of nesting
`./py`; the rerun used that as written.

## Task 3: Refute

No value was refuted: every value agrees with my independent exact computation. The weaknesses are in the checks and
in the assembly. Every mutation below ran in a separate copy `mut/<name>/`, made from `solver_rerun/` by `st2_mutations.py`
(log `out2/mutations.json`).

| # | Target | Mutation | Outcome | Verdict |
|---|---|---|---|---|
| m1 | s09 "… explained by ANTIUNITARY" (11 lines) and "no unresolved zero" | set `antiunitary_residual = 1.0` at U5/U6, as if A′ did not fix Φ and ξ | all 12 lines still PASS | **cannot fail**: s09 attaches reason labels without reading any verification |
| m2 | same s09 lines | inject spurious zeros: U1 `form_cross`, U1 `DNxi_perp_comp_ie_t`, U6 `form_ie_t` | PASS, "explained by ANTIUNITARY / ORBIT", at points that have no such tangent | **cannot fail** |
| m5 | s04 antiunitary check, then s09 | drop the minus sign in A′ = −ρ(R_y(π))∘Θ | s04's four A′ lines **FAIL**, which is correct; s09's 11 ANTIUNITARY lines and its "no unresolved zero" line still **PASS** | the s04 check works; confirms m1 in a live run |
| m3 | s10 "48 is an eigenvalue only at level 6" | insert levels 2 and 4 into the 4-dim sector (`item1.json`) | still PASS, along with all H1 and H2 lines | **cannot fail**: n(n+2)=48 only at n=6, and level 6 is always present. No PASS line checks "no level below 6", which the λ₄ sign needs |
| m4 | s03 "every element of the generated group fixes [u] (character values of modulus 1)" | wrong character χ(Rz(π/3)) = +1 at U4 | the generator line FAILs; this line still PASSes | implied by the generator lines; it never tests the characters |
| m9 | s01 "exact dims vanish at every odd n" | add 7 to every even-n dimension | odd-n line PASS (the neighbouring count check FAILs) | low power: it can only fail if χ_σ(−h) ≠ χ_σ(h) |
| m10 | `results.json` field `xi_transforms_by_same_character` | wrong U4 character, run s04 and s99 | s04 covariance lines **FAIL**; `results.json` still says `True` | **literal** in `s99_assemble.py` |
| m11 | `results.json` field `sign_of_lambda4` | make every λ₄ positive in `numeric.json`, run s99 | still "negative for every g != 0" | **literal** in `s99_assemble.py` |
| m6 | s04 "with kappa the equation has no component orthogonal to Φ" | double the eigenvalue in the pseudo-inverse solve | FAILs (5e-4), while the solvability line PASSes | a real check, independent of the solvability line |
| m7 / m8 | s12 item-12 residual line | Casimir factor 4 → 3.9 (U1 only; m8 is the unmutated control) | residual line FAILs (9.9e-4, 5.7e-4); quadrature-vs-algebra line unchanged at 1e-29; control all PASS | real check. Unmutated, it equals the quadrature-vs-algebra line (9.80e-30 vs 9.79e-30), because s12 rebuilds ξ from the algebraic N; its only extra content is the Casimir |

**Hidden conventions and dependencies.**
1. **Choice of σ and η.** σ is made real through Θ-real η, verified numerically (imaginary part ≤ 1.9e-51). Mine uses
   an arbitrary orthonormal η. The values do not depend on this (confirmed).
2. **How the exact values were obtained.**
   - The item-5 and item-7 values (e_t / τ_x components, ‖κ‖²) are *identified* from 80- and 120-digit numerics, not
     derived.
   - The "exact" s02 values are built on Γ-constants (τ_J, ‖r_K‖²). Those constants are themselves identified from
     the same section algebra (the same CG table and intertwiners) that s04 uses.
   - So the s04-vs-s02 agreement lines (< 1e-70) share that machinery and cannot detect an error in it.
3. **What is and is not cross-checked.**
   - Its independent routes are the item-12 quadrature, which checks N but not DN, and the λ₄ identity along Φ.
   - The perpendicular components of Π₆DN_Φ[ξ] and κ have no independent check in its own work.
   - Its Racah CG is cross-checked against sympy only for the spin pairs (3,3), (6,3), (3,2) and (5,4).
   - My exact values agree with all of them, so nothing went wrong. The independence claimed in its RETURN is weaker
     than it reads.

**Readings that differ from mine** (none changes a value):
1. Item 4: it reports level 6 as a computed 0; I left level 6 out (0 by definition).
2. Item 9: it counts Im Q, Im⟨Φ,DNξ⟩ and ‖Π₆ξ‖² as zeros.
3. Item 7: its with-κ zero is numerical; mine is exact.
4. Item 3, U2: it uses the coset representative Rx(π); I used ry(π). Same group, same character.
5. Item 12: 30-digit arithmetic against 1e-25, where I used float64 against 1e-10.
6. Items 6 and 11: the same polarized cross-term reading. It also reports |grad r̂₆|² and dr̂₆/ds.

## Task 4: Its arguments (items 8, 9, 10)

**Item 8: correct and complete. Its sign argument is stronger than mine.**
- *Common part.* Its formula λ₄ = −3g²Σ_{n≠6}‖Π_nN‖²/μ_n and its exact reason why κ does not contribute
  (2X + X̄ with X = Q⟨Φ,κ⟩ = 0) are the same as mine.
- *My sign argument.* It used the computed nonzero level norms (item 4), so it is valid only at the six points.
- *Its sign argument.*
  1. Σ_{n≠6}‖Π_nN‖² = ∫|Φ|⁶ − Q².
  2. By Cauchy–Schwarz this is > 0 unless |Φ|² is constant.
  3. |Φ|² is not constant, because its level-12 part ρ₆(u′)⊗r₆ is nonzero: ‖r₆‖² = d(7−d)/7, and ρ₆(u) ≠ 0 because
     the product of two nonzero binary sextics is nonzero.
- *Scope.* This holds at every critical Φ in either sector, not only the six. It needs "no level below 6", which it
  states as n ≥ 8 (item 1).
- *Exact check.* I verified the quantity exactly at all 12 cases (`out2/theorem_checks.json`). For example, at U6 in
  the 3-dim sector ∫|Φ|⁶ − Q² = 43538303200/351543937329 > 0.

**Item 9: correct reasons, but two are weaker than they read.**
- *What is correct.*
  - The symmetry proofs S1–S3 and A′ are correct. A′ needs σ real; it verifies this numerically. The exact reason is
    that the characters are real and Θ² = +1 in integer spin. In my stage 1 this was the exact identity
    conj(P) = RPR, re-verified today.
  - Its Sym³ argument for U4 at n = 16 is the same as mine.
  - The ORBIT argument (DN_Φ[XΦ] = XN(Φ)) is correct and complete.
  - The item-2 zeros at U5/U6 follow correctly from exact stationarity plus its structural identity.
- *Weakness 1: the with-κ zero.* The reason given for "the equation's ⊥ component with κ = 0" is "by construction".
  - κ = L⁺b gives Lκ = b − P_ker b. So the zero requires P_ker b = 0, and m6 shows this is not automatic.
  - It verifies P_ker b = 0 numerically (≤ 1e-124). The exact reason is the Noether identity, which it gives in item-10
    step 4 but does not cite in item 9. s09 does not record this zero at all.
- *Weakness 2: the machine check.* "117 PASS, none unresolved" certifies only that a label was attached (m1, m2, m5).
  The prose proofs carry the weight.

**Item 10: correct and complete, at the same level as mine.**
- *Architecture.* The same as mine: Lyapunov–Schmidt, the blow-up b = a(Φ+k) with λ = 48 + a²μ, a slice S = W ⊖ T,
  and Noether identities.
- *Transversality.* Its step is more explicit than mine: R ∈ T_Φ ∩ T_v^⊥ = {0}, using only X_x and X_y at U1/U2, where
  the orbit dimension jumps. My stage 1 relied on the slice theorem stated in words.
- *Hypotheses.* H0–H5 are all needed and correctly placed.
  - H2 is exact: the SymPy Hessian rank combined with its structural identity L = κ_d·Hess r̂₆.
  - H0 is true (item 1), but its PASS line cannot fail (m3).
- *Results fields.* `existence_for_all_small_a` is computed from the H1 and H2 flags, not hard-coded.
- *Nothing in it is weaker than it reads.*

## Task 5: Grading THEOREM.md

Supporting scripts: `st2_theorem.py` (log `out2/theorem_checks.log`, data `out2/theorem_checks.json`), `st2_exact_m.py`
and `st2_u5_weights.py`. They run **90 exact or pointwise checks, all PASS**. Each of the 79 checks that has a separate
mutant catches it, and the log marks the rest. All checks are on my stage-1 objects (`out/main_all.pkl`, exact MQ).

**T.a: ESTABLISHED.**
- F is a real-analytic (cubic) map H^{s+2} × ℝ → H^s for s > 3/2, since H^s is an algebra.
- A = −Δ−48 is self-adjoint and elliptic, hence Fredholm of index 0 with kernel K (item 1: dim Hom(σ,V₃) = 1, and
  n(n+2) = 48 only at n = 6).
- Item 1 gives the levels: 6, 10, 14, 16, 18 (3-dim sector) and 6, 8, 12, 14, 16, 18 (4-dim sector). So the gap is
  72 or 32, and no level lies below 6.
- *Remark.* "Nearest other level is 8, gap 32" does not by itself exclude levels 0 and 2 in the 4-dim sector
  (distances 48 and 40). Item 1 does. This matters for the λ₄ lemma.

**T.b: ESTABLISHED.**
- For h ∈ H: h·F(ψ) = F(h·ψ) = F(χψ) = χF(ψ), using rotation invariance and N(χψ) = |χ|²χN(ψ).
- Π_⊥ and A⁻¹ commute with rotations, so ξ ∈ X_H.
- Checked exactly on my ξ at every level:
  - the exact generators of U2, U3, U4, U6 act by χ, and the mutant −χ fails;
  - the torus weights are 3 (U1) and 0 (U2);
  - at U5, every left weight is ≡ 2 (mod 5), which is exactly χ(Rz(2π/5)) = e^{4πi/5}.

**T.c: ESTABLISHED, with one reading made explicit.**
- The LS reduction in X_H is a gradient because F(X_H) ⊂ X_H, so ∇(E|_{X_H}) = F.
- Oddness: F(−ψ) = −F(ψ), and uniqueness of the range solution gives w(−k) = −w(k).
- Radial part: q = (48−λ) + gQ(θ)r + O(r²), with ∂_λq = −1.
- *The reading.* "Tangential part" must be the component in the **moving** tangent space at θ, orthogonal to u(θ).
  This is the θ-gradient of the reduced functional. There, (48−λ)k contributes nothing, and the part is O(a³) and odd,
  so G = Π_{T_θ}N(Φ_θ) + O(r) = ¼∇Q(θ) + O(r).
- With the fixed projection Π_T of the Setting, the term (48−λ)aΠ_TΦ_θ is O(a). Dividing by a³ then works only after
  substituting λ from the radial equation. That gives G = Π_T(N(Φ_θ) − Q(θ)Φ_θ) + O(r), whose θ-derivative at θ₀ is
  again L_T. Either way the Jacobian at (θ₀, 0, 48) is block-triangular, even diagonal, with blocks −1 and ¼Hess_T Q.
- The ¼ comes from dQ(Φ_s)/ds = 4Re⟨N(Φ),E⟩, the first-order form of the bridge.
- The orbit components vanish by invariance: Noether identities plus the transversality of the orbit tangents at
  nearby points.

**T.d: ESTABLISHED.** It is the IFT in (λ, θ) on the blown-up problem, together with uniqueness of the range solution
in a small ball. For solutions with θ near θ₀, μ = (λ−48)/a² is automatically near gQ. Uniqueness holds only within
X_H, as stated. The equivariance under rotations follows because the construction is canonical.

**T.e: ESTABLISHED.**
- An analytic germ's coefficients solve the order-by-order equations. In the stated gauge the formal solution is
  unique: at each order, the range part comes from A⁻¹, the tilt from L_T (invertible by T.c's hypothesis), λ from the
  radial equation, and the orbit components vanish by Noether.
- *Implicit step: reparametrization.* T.c uses a = ‖k‖, with k = a·u(θ). T.e uses a = ⟨Φ,ψ⟩.
  - The change is a_e = a_c·⟨Φ, u(θ(r))⟩ = a_c(1 + O(a_c⁴)). It is analytic and odd, so the conclusions survive.
  - At U6, ⟨Φ,u(z)⟩ is real because z(r) is real: the S lemma forces v_y = 0 at every order.
- *Scope.* T's ansatz allows a block tilt at every odd order. The worklist's §2 ansatz, with ζ ⊥ block, does not at
  a⁵. The two agree through a³, which is what "begins as in §2" requires.

**E1: ESTABLISHED.** `st2_exact_m.py` checks m = 1 exactly for the exact generators (U1, U2 via Rz(π/2) and Rz(π/3),
whose weights mod 12 decide; U3; U4), with u in the χ-subspace. The mutant −χ̄ fails. So K_H = ℂΦ, U(1) leaves one real
amplitude, and T applies with a vacuous nondegeneracy condition.

**E2: ESTABLISHED.** Every stated fact holds:
- **Weights at U5:** v₂ and v₋₃ carry the same phase e^{4πi/5} under Rz(2π/5); the U5 χ-subspace is
  m ≡ 2 (mod 5) = span{v₂, v₋₃}, exactly.
- **U5 orbit direction:** P_⊥(i·conj(J_z)u) = −(2/5)√39 · (i·e_t) exactly, a real multiple.
- **U6 chart:** Rz preserves the chart iff e^{6iθ} = 1, and Rz(π/3) acts as z ↦ −z (exact).
- **No continuous symmetry at U6:** no nonzero element of su(2) preserves K_H (Gram determinant 2916 ≠ 0).
- **U6 dimension:** m = 2 exactly.
- **The tilt equation:** L_T v = −Π_T DN_Φ[ξ] is the tangential part of gΠ_K DN_Φ[v+ξ] − λ₂v − λ₄Φ.
- *Not stated in E2:* it does not say that T.c's nondegeneracy holds. I supply it:
  - U5: Hess_T Q = −224/165 (3-dim) and −42/55 (4-dim).
  - U6: det Hess_T Q = 5770240/7913763 and 202860/879307. All nonzero.

**Bridge: ESTABLISHED.**
- *Derivation.* Take unit E with Re⟨Φ,E⟩ = 0. Then
  |cΦ+sE|⁴ = c⁴|Φ|⁴ + 4c³s|Φ|²Re(Φ†E) + c²s²(4Re(Φ†E)² + 2|Φ|²|E|²) + …, so
  d²/ds²|₀ = −4Q + 4∫(2Re(Φ†E)² + |Φ|²|E|²) = 4(⟨E,DN_ΦE⟩_ℝ − Q).
  Great circles are geodesics, so at a critical point this is the Riemannian Hessian.
- *Check.* In all 12 cases (U5: e_t, i·e_t, cross; U6: τ_x, τ_y, cross; both sectors), Q″ was computed exactly from a
  quartic interpolation of ∫|cΦ+sE|⁴, which uses N only and not DN. It equals 4 × the item-6 form every time; for
  example Q″(τ_x) = 25760/18447 = 4·6440/18447.

**λ₄ lemma: DEFECT.**
- *What is correct.*
  - ⟨Φ,DN_Φ v⟩ = 2X + X̄ holds (checked exactly with a random block v not orthogonal to Φ, all 12 cases; the mutant 3X
    fails).
  - X = ⟨N(Φ),v⟩ = Q⟨Φ,v⟩ = 0 for the tilt (exact, all 12 cases).
  - The formula λ₄ = g⟨Φ,DN_Φξ⟩ = −3g²Σ‖Π_nN‖²/μ_n holds, with the sum over n ≠ 6. As written, the sum includes n = 6,
    where μ₆ = 0.
- *The defect.* The sign claim "λ₄ < 0 since Q ≠ 1" omits its argument and leaves a **necessary hypothesis unstated**:
  every non-block level of the sector must lie above 6, so that every μ_n > 0. Without it, Q ≠ 1 does not fix the sign.
  Nothing in T states it (see T.a).
- *Rescue.*
  1. Item 1 gives μ_n > 0.
  2. Q ≠ 1 implies |Φ|² is not constant (if it were, it would equal 1, and then Q = 1).
  3. By strict Cauchy–Schwarz, Σ_{n≠6}‖Π_nN‖² = ∫|Φ|⁶ − Q² > 0.
  4. Hence λ₄ < 0 for g ≠ 0. Exact values are in `out2/theorem_checks.json`.

**S lemma: ESTABLISHED.** I derived the section-level compatibility rather than assuming it.
- *Section-level form.* Write sections in the P-trick form f = ψη† (f = fP, f(gh) = f(g)D³(h) for h ∈ Γ). Define
  **(Sf)(g) = conj(f(r⁻¹g))·R**, with r = R_z(π/3) and R_{m,m′} = (−1)^m δ_{m,−m′}. Then:
  - Since conj(D) = RDR, (Sf)(gh) = (Sf)(g)D³(h), and Sf·P = Sf. This needs conj(P) = RPR, which holds exactly
    because P is a real-coefficient class-sum polynomial and the characters are real (checked exactly, both sectors).
  - In the η-representation this is ψ ↦ conj(ψ)W with W ∈ Hom_Γ(σ, conj σ).
  - |Sf| = |f| pointwise, so N(Sf) = S N(f). Δ is real and commutes with left translation, so S maps solutions to
    solutions.
- *Action on terms.* On a term xᵀD^J(g)Y: **x ↦ conj(D^J(r))Θx**, and column b of Y ↦ (−1)^b Θ(column −b). On the block,
  Y = P is preserved, so Φ_u ↦ Φ_{conj(D³(r))Θu}.
- *Numerical confirmation.* At random g, on all levels of U6 and U5: formula deviation 4.5e-16; sector preservation
  5.6e-16; commutation with N 1.7e-15. The mutants (r in place of r⁻¹, or R dropped) deviate by 0.77 and 1.17.
- *Chart action (exact).*
  - Θu(z) = −u(−z̄);
  - Rz(π/3): u(z) ↦ −u(−z);
  - S: u(z) ↦ u(z̄), with **no phase**, checked at generic z = 1/3 + 2i/5.
- *Consequences, checked exactly on my objects in both sectors:*
  - SΦ = Φ.
  - Sξ = ξ at every level (10, 14, 16, 18 and 8, 12, 14, 16, 18).
  - The forcing Π₆DN_Φ[ξ] is S-even.
  - Sτ_x = τ_x and Sτ_y = −τ_y.
  - L commutes with S, so B(τ_x,τ_y) = 0.
  - κ is S-even and its τ_y component is exactly 0.
- *Remarks, neither changing the grade.*
  - "Fixes Φ once the phase of Φ is chosen": the freedom must be taken in the lift of S, or in η, which Φ and τ_x share.
    Rephasing u alone would break "τ_x is S-even". With the lift above, no choice is needed.
  - "Hence v_y = 0" also uses L_T(τ_y,τ_y) ≠ 0 (56/429 and 21/286), which is T.c's nondegeneracy hypothesis.

### THEOREM.md against my stage-1 item-10 argument, step by step

| Step | THEOREM.md | My stage 1 | Comparison |
|---|---|---|---|
| Function spaces | H^{s+2} → H^s, s > 3/2 | H² sections, N analytic on H² | equivalent |
| Use of symmetry | restrict first to the twisted fixed space X_H of the stabilizer | full space, G = SU(2)×U(1), slice in the full block; the Fix(Ĥ) route only as an alternative for U1–U4 | T's route is my alternative, extended to m = 2 |
| Blow-up | radial part ÷ a, tangential ÷ a³, via oddness | ψ = a(Φ̂ + a²w), ÷ a³, polynomial in a² directly | equivalent; mine avoids the oddness step and the frame issue noted in T.c |
| Kernel variable and gauge | k = a·u(θ), θ on a slice of ℙ(K_H); gauge ⟨Φ,ψ⟩ = a in T.e | Φ̂ = Φ + κ̃, κ̃ in the slice of Φ^⊥∩B; a = ⟨Φ,ψ⟩ from the start | T needs the reparametrization noted in T.e |
| Nondegeneracy | Hess_T Q invertible on a 0-, 1- or 2-dimensional slice (item 6) | ker L on Φ^⊥∩B (12 real dims) = orbit tangent (items 6/7, exact) | T needs less. Mine implies T's, because L preserves the isotypic splitting and ker L ∩ K_H is the orbit part inside K_H |
| Orbit components | invariance of the reduced functional | Noether identities ⟨Xψ,F⟩ ≡ 0 plus the slice theorem | same idea; neither spells out transversality (the other worker does) |
| Existence | yes, U1–U6 | yes, U1–U6 | same |
| Uniqueness | within X_H, [Π_Kψ] ∈ U, gauge | modulo G in a cone around the orbit, at small amplitude | **mine is stronger**; T's is what its route gives |
| Identification of the coefficients | T.e: formal = Taylor, at all orders | the a² coefficient of κ̃ is the item-7 κ (slice condition = κ ⊥ ker L) | same at the order asked; T is more general |
| Regularity | ψ/a and λ analytic in a² | ψ odd and analytic in a (H², C^∞), λ even | same |
| λ₄ and its sign | lemma with an unstated hypothesis (DEFECT) | complete for the six points, using item 1 (μ_n > 0) and item 4 (nonzero level norms) | mine complete but local; the other worker's Cauchy–Schwarz argument generalizes |
| Degeneracies at a = 0 | trivial branch removed by ÷ a; ∂_λ = −a | trivial branch and the 14-dimensional kernel by blow-up; orbit directions by the slice | same |
| Dependence on item 3 | essential (H, χ, m, K_H) | none on the main route | T's proof inherits item 3's certification level. The facts T uses are now certified exactly (m and ξ ∈ X_H), so its route is complete |

## Files

- **Scripts:**
  - `st2_compare.py` (task 1) and `st2_item3.py` (task 1, item 3);
  - `st2_mutations.py` (task 3);
  - `st2_theorem.py`, `st2_exact_m.py`, `st2_u5_weights.py` (task 5);
  - `st2_build_json.py` (writes `audit_stage2.json`).
- **Outputs:**
  - `out2/compare.json`, `out2/item3_crosscheck.json`, `out2/mutations.json`;
  - `out2/theorem_checks.log` and `.json`, `out2/exact_m.json`, `out2/u5_weights.json`.
- **Copies:** `solver_rerun/` (clean rerun, with its regenerated `results/` and `results.json`), and `mut/<name>/` (one
  mutated copy per mutation).
- **Log:** `solver_rerun_run_all.log`.

## Consulted-material manifest

- **Read in this room:**
  - Briefs, conventions and statements: `BRIEF_stage2.md`, `BRIEF.md`, `worklist.md` (byte-identical to
    `solver_work/worklist.md`), `THEOREM.md`.
  - My stage-1 returns: `METHOD.md`, `AUDIT_STAGE1.md`, `audit_results.json`.
  - The other worker's files: `solver_work/RETURN.md`, `results.json`, `run_all.log`, `run_all.py`, `common.py`,
    `s00`, `s01`, `s02`, `s03`, `s04`, `s09`, `s10`, `s12`, `s13`, `s14`, `s99`.
  - My stage-1 code: `su2.py`, `engine.py`, `s_main.py`, `mq.py`, `numD.py`, `linalg_mq.py`, `s03_stabilizers.py`.
  - My stage-1 outputs: `out/main_all.pkl`, `out/main_all.json`, `out/group.pkl`, `out/item3.json`, `out/item01.json`.
  - The interpreter wrapper `py` and the sandbox profile `.room.sb`.
- **Outside the room:** nothing. No web, papers, books or other files.
- **Software:** Python 3.12 standard library, numpy, mpmath and sympy. The other worker's scripts call
  `sympy.physics.quantum.cg` as their own cross-check; my code does not.
- **Item 1 declaration (stage 1, unchanged): DERIVED.**
