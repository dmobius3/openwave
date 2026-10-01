# M8.14 method note: fixed-edge stability of the conic carrier

Task: [`tasks/m8_14_task_details.md`](../tasks/m8_14_task_details.md), pre-registered in [#608](https://github.com/openwave-labs/openwave/pull/608), run 2026-10-01. Verdict under the frozen rule: **Reproduced.** Two blind rooms, offline, each VALID, each establishing the reduction and the full spectrum and each reproducing the frozen bottom at all five widths; a blind two-stage auditor; both controls pass. Records of the run: [`m8_14/run/`](../m8_14/run/).

## 1. Equations first

`R = 1` in every reported number. `α₀ = πR/(2W)`.

### 1.1 The object

`S³(R) ⊂ ℝ⁴` and `ℝP³(R) = S³(R)/{±I}`. On the great sphere `x₄ = 0`,

```text
X(y, w) = R (cos(y/R) cos(w/R), cos(y/R) sin(w/R), sin(y/R), 0)
```

The band `M(W)` is `X(y, w)` for `0 ≤ y ≤ πR/2` and `X(y, −w)` for `πR/2 < y ≤ πR`, with `|w| ≤ W`, passed to `ℝP³`. The segment `y = πR/2` is the cone point `p`, and `(0, w) ∼ (πR, −w)` is the seam. The metric is `dy² + cos²(y/R) dw²` and the area element `|cos(y/R)| dy dw`. The normal is the constant `ν = e₄`, which the antipodal map reverses, so a normal graph `φν` is defined on `ℝP³` exactly when `φ(πR, −w) = −φ(0, w)`. The admissible class is the completion, in `∫(|∇φ|² + φ²) dA`, of the smooth such `φ` that vanish on `w = ±W` and near `p`.

### 1.2 The second variation (A1)

For the normal-geodesic variation `F_t = cos(tφ/R) X + R sin(tφ/R) e₄`, the induced metric is exactly

```text
g_t = cos²(tφ/R) g + t² dφ ⊗ dφ,     dA_t = cos(tφ/R) √(cos²(tφ/R) + t²|∇φ|²) dA
```

so the first variation vanishes and

```text
A''(0) = q_J[φ] = ∫_M (|∇φ|² − (2/R²) φ²) dA,     J = −Δ − 2/R²
```

with no edge, seam or vertex term, since no integration by parts is taken, and `q_J` extends to the admissible class by continuity.

### 1.3 The reduction (F1a)

On the upper sheet set `(ψ, θ, u) = (y/R, w/R, φ)`, and on the lower sheet `(ψ, θ, u) = (y/R − π, −w/R, −φ)`. Then

```text
ds² = R² (dψ² + cos²ψ dθ²),   φ(πR, −w) = −φ(0, w)  ⇔  u continuous across ψ = 0
```

The seam becomes an interior line, the edges become `θ = ±W/R`, and the two sides of the fiber become the two poles `ψ = ±π/2`, both the cone point. The map is unitary in `L²` and preserves the Dirichlet integral, so the band's problem is the Dirichlet problem on the spherical lune `{|ψ| < π/2, |θ| ≤ W/R}` of opening `2W/R`.

### 1.4 Sectors and the cone (A2)

With `e_k(θ) = sin(ν_k(θ + W/R))`, `ν_k = kπR/(2W) = kα₀`, `k ≥ 1`, sector `k` carries

```text
ℓ_k f = −(cos ψ f′)′ / cos ψ + ν_k² f / cos²ψ     on L²(cos ψ dψ)
```

At a pole the indicial roots are `±ν_k`, and `r^{−ν}` is square-integrable against `r dr` iff `ν < 1`. So a sector is limit-circle iff `ν_k < 1`. For `W ≤ πR/2` every `ν_k ≥ 1` and every sector is limit-point (`ν = 1` at `W = πR/2` is limit-point). At `W/R = 17/10`, `ν₁ = 5π/17 ≈ 0.924`: sector 1 is limit-circle at both poles, and the admissible class selects the Friedrichs realization, the one with no `r^{−ν}` part at either pole.

### 1.5 The spectrum (F1b, A3)

With `x = sin ψ` and `f = (1 − x²)^{ν/2} g`, the sector equation becomes Gegenbauer's, with polynomial solutions `C_n^{(ν+1/2)}` exactly at

```text
λ_{k,n} = (ν_k + n)(ν_k + n + 1) / R²,   k ≥ 1, n ≥ 0
```

and these are complete. The bottom is `α₀(α₀ + 1)/R²` at `(k, n) = (1, 0)`, so `J`'s bottom is `(α₀ − 1)(α₀ + 2)/R²`: positive for `0 < W < πR/2` (index 0, nullity 0), zero at `πR/2`, negative at `17/10`. The ground state is `±|cos(y/R)|^{α₀} cos(πw/2W)`, the sign flipping across the cone.

### 1.6 The frozen cross-check (F2)

```text
a(u, v) = ∫ (|cos y| u_y v_y + u_w v_w / |cos y|) dy dw,    m(u, v) = ∫ |cos y| u v dy dw
nodes  π/2 ∓ (π/2)(k/16)²,  8 cells across,  four uniform refinements  (levels 0..4)
p = log₂((λ₂ − λ₃)/(λ₃ − λ₄)),   λ_ext = λ₄ + (λ₄ − λ₃)/(2^p − 1),   err = |λ_ext − λ₄|
```

P1 elements on `[0, π] × [−W, W]`, the seam anti-periodic, zero data on `w = ±W` and the fiber. Tolerance: a relative `10⁻³` on `α₀(α₀ + 1)`.

## 2. What the rooms derived, none of it supplied

The spec sheet gave the object, the variation, the admissible class and the protocol, and named neither the lune, nor the reduction, nor any formula. Both rooms found the reduction of § 1.3 themselves: Solver A by the latitude-longitude unfolding with the sign `u = −φ` on the lower sheet, Solver B through an explicit unitary `U` onto the lune with the seam as its equator arc. Both then derived § 1.4 and § 1.5 by separation, Frobenius and Gegenbauer, with a completeness argument by Weierstrass density.

| Reading | Taken by | Effect |
| --- | --- | --- |
| T2's core read as smooth on the surface, so all `y`-derivatives anti-periodic at the seam; under values-only matching the operator is not symmetric | both | the form domain, and so T3 to T5, is the same under either reading |
| "uniform refinement" as midpoint insertion, not regraded at `2N` | both | the regraded alternative moves no bottom by more than `3.1 × 10⁻⁷` relative (auditor, `s2_regrade.py`) |
| Richardson with the estimated `p` | both | the `p = 2` reading differs by under `10⁻⁷` relative ([#608 review](https://github.com/openwave-labs/openwave/pull/608)) |

## 3. Equation-to-code map

Links resolve on `main` once the record merges. The rooms' files are in [`m8_14/run/rooms/`](../m8_14/run/rooms/).

| Equation | Solver A | Solver B |
| --- | --- | --- |
| § 1.2 exact density and `A''(0)` | [`t1_second_variation.py:54`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_a/t1_second_variation.py#L54) | [`t1_t3_symbolic.py:43`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_b/t1_t3_symbolic.py#L43) |
| § 1.4 indicial roots and the `ν = 1` borderline | (argued in `RETURN.md` § 3; `ν_k` in [`t3_spectrum.py:98`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_a/t3_spectrum.py#L98)) | [`t1_t3_symbolic.py:58`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_b/t1_t3_symbolic.py#L58), [`:71`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_b/t1_t3_symbolic.py#L71) |
| § 1.5 Gegenbauer eigenfunctions | [`t3_spectrum.py:43`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_a/t3_spectrum.py#L43), [`:65`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_a/t3_spectrum.py#L65) | [`t1_t3_symbolic.py:78`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_b/t1_t3_symbolic.py#L78), [`t3_t5_exact.py:28`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_b/t3_t5_exact.py#L28) |
| § 1.6 weight `\|cos y\|` | [`t4_fem.py:50`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_a/t4_fem.py#L50) | [`t4_fem.py:36`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_b/t4_fem.py#L36) |
| § 1.6 graded nodes and refinement | [`t4_fem.py:54`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_a/t4_fem.py#L54) | [`t4_fem.py:41`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_b/t4_fem.py#L41), [`:48`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_b/t4_fem.py#L48) |
| § 1.1 the seam, sign `−1` with the `w`-reversal | [`t4_fem.py:79`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_a/t4_fem.py#L79) | [`t4_fem.py:69`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_b/t4_fem.py#L69) |
| § 1.6 `p`, `λ_ext`, `err` | [`t4_fem.py:214`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_a/t4_fem.py#L214) | [`t4_fem.py:193`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/solver_b/t4_fem.py#L193) |

The auditor's stage-1 recomputation, by two methods neither of which is T4's: [`fem2d_q2.py:93`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/auditor/fem2d_q2.py#L93), biquadratic elements in `(y, w)` with a cubic grading, and [`fem1d_sectors.py:32`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/auditor/fem1d_sectors.py#L32), exact sines in `w` with degree-8 elements in `y`. Its seam check reads each room's own dof map, [`s2_audit.py:166`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m8_mit/research/m8_14/run/rooms/auditor/s2_audit.py#L166).

The maintainer's route, run before any room opened: [`frozen_check.py`](../m8_14/run/scripts/frozen_check.py) recomputes the frozen table at 30 digits and encodes the frozen outcome rules and control gates, with eight planted cases; [`compare.py`](../m8_14/run/scripts/compare.py) scores a return from its own JSON and re-derives its `p`, `λ_ext` and `err` from its levels, with six planted cases.

## 4. Results against the frozen claims

| ID | Solver A | Solver B | Adjudication |
| --- | --- | --- | --- |
| F1a | the reduction of § 1.3, by unfolding | the same, by the unitary `U` | ✅ agrees in both, ESTABLISHED, SUPPLIED (§ 5.2) |
| F1b | `(ν_k + n)(ν_k + n + 1)/R²`, complete | the same, complete | ✅ agrees in both |
| F2 | Reproduced at all five widths | Reproduced at all five widths | ✅ Reproduced |
| A1 | ESTABLISHED | ESTABLISHED | ✅ audited argument |
| A2 | ESTABLISHED | ESTABLISHED | ✅ audited argument |
| A3 | ESTABLISHED | ESTABLISHED | ✅ audited argument |
| N1 | `2.000000050` | `2.000000050` | ✅ miss `2.5 × 10⁻⁸` relative |
| N2 | `J` bottom `−0.2222301` | `J` bottom `−0.2222301` | ✅ about 1,600 times its error estimate `1.38 × 10⁻⁴` |
| **Verdict** | VALID | VALID | ⭐ **Reproduced** |

The two rooms' finite-element numbers agree at every width and level to `8.2 × 10⁻¹⁴` relative, and with the reviewer's own implementation of the protocol in #608 to `10⁻¹²`.

| `W/R` | extrapolated bottom | `p` | `err` | frozen `α₀(α₀ + 1)` | relative miss | outcome |
| --- | --- | --- | --- | --- | --- | --- |
| 1/4 | 45.761589282 | 1.9967 | 5.60e-3 | 45.761603 | 3.0e-7 | Reproduced |
| 1/2 | 13.011196550 | 1.9994 | 1.18e-3 | 13.011197 | 3.9e-8 | Reproduced |
| 1 | 4.038197481 | 2.0002 | 3.19e-4 | 4.038197 | 1.3e-8 | Reproduced |
| 7/5 | 2.380875543 | 2.0004 | 1.84e-4 | 2.380875 | 2.3e-8 | Reproduced |
| 3/2 | 2.143820313 | 2.0003 | 1.66e-4 | 2.143820 | 2.4e-8 | Reproduced |
| π/2 (N1) | 2.000000050 | 2.0002 | 1.55e-4 | 2 | 2.5e-8 | control passes |
| 17/10 (N2) | 1.777769908 | 1.9997 | 1.38e-4 | 1.777770 | 3.4e-8 | control passes |

Every error estimate sits within the tolerance, read as the frozen text states it, a relative `10⁻³`; the largest relative one is `1.2 × 10⁻⁴`, at `W/R = 1/4`. Under an absolute `10⁻³` reading, which the frozen text does not support, `W/R = 1/4` and `1/2` would be Unresolved. Each outcome is scored on the room's reported `p`, `λ_ext` and `err`, which equal the values recomputed from its own levels to `2.6 × 10⁻¹⁵`. The estimate bounds the extrapolated value's actual error by two to four orders of magnitude, as both rooms state: it measures the finest level, not the extrapolation.

| Diagnostic | Result |
| --- | --- |
| G1, the route | both rooms reduce to the lune before separating; A by coordinates, B by an explicit unitary |
| G2, the order `p` | between 1.9967 and 2.0004 at all seven widths, including `17/10`, where `ν₁ < 1` |
| G3, the ground state | the `(1, 0)` mode in both: Solver B states `±\|cos y\|^{α₀} cos(πw/2W)` with the sign flipping across the cone, the author's form; Solver A gives the same function in lune coordinates, `cos^{α₀}ψ sin(α₀(θ + W/R))`, with the sign `(−1)^{n+k} = −1` on the lower sheet |

**The seam diagnostic** (the dated note in the task's AS RUN section). At the finest level, each room's lowest six levels sit above F1b's formula by `7.7 × 10⁻⁵` to `8.7 × 10⁻⁴` relative at every width, the sign a conforming method must give, and differ from the auditor's values by the same amounts, since the auditor's agree with the formula to `10⁻¹²`. The second level is between 1.32 and 3.00 times the bottom at every width, so neither room returns the repeated bottom that a seam deleted with both ends free would show.

## 5. Audit record

### 5.1 The auditor, blind, two stages

Stage 1 recomputed the seven widths by two methods of its own before seeing any return, with errors on the bottom between `1.8 × 10⁻¹²` and `4.6 × 10⁻¹¹`, and committed them; stage 2, a new launch, found all 17 stage-1 files unchanged.

| Check | Result |
| --- | --- |
| V1, validity | both VALID: every explicit requirement read from the code and confirmed by checks on the rooms' own assembly (mesh, seam sign, eigen-residuals below `4.1 × 10⁻¹²`, Sylvester inertia showing no missed eigenvalue, `p` and the extrapolation recomputed from the levels); reruns reproduce each JSON to `8.4 × 10⁻¹¹` |
| V1, immaterial | Solver B's "exactly zero" free-node gradient on fiber-edge triangles is zero to roundoff (`1.1 × 10⁻¹⁶`); Solver B's planted run switches two defects at once; Solver A's C4.5 is a 5 % gross check |
| V3, comparison | 14 of 14 below a relative `10⁻³`, the largest `3.0 × 10⁻⁷`; no defect to locate and no stage-1 correction |
| G, grading | every step ESTABLISHED in both rooms' T1, T2, T3 and T5, except Solver A's T3 step 3, ESTABLISHED, SUPPLIED: the auditor supplied that the 2D operator is the direct sum of the sector operators, by two integrations by parts in `θ` and Parseval |

Each check the auditor printed as PASS has a planted twin that printed FAIL. One of its own first versions did not fire, and it restated that criterion and reran.

### 5.2 The adversarial audit

An independent agent, given every return, both stages and the maintainer's scripts, tried to refute eight claims with its own code, importing none of the rooms', the auditor's or the maintainer's. Its report is landed as [`m8_14/run/adversarial_audit.md`](../m8_14/run/adversarial_audit.md). It refuted none.

| Claim | Result |
| --- | --- |
| Validity | not refuted: both seams read as `φ(π, −w) = −φ(0, w)` in direction and sign, and confirmed on each room's own map with a test function odd in `w` (error `2.2 × 10⁻¹⁶`); its own P1 assembler reproduces every level value of both rooms to `1.8 × 10⁻¹³` |
| F1a | not refuted, with two supplied parts named in each room: **S-a**, the value-only seam class has the same form closure as the smooth core, asserted in both rooms without the argument (stage 2 graded Solver B's sentence plain ESTABLISHED); **S-b**, the vertex-excluded lune problem is the lune's Dirichlet problem, since points have zero capacity in two dimensions. Under the declared reading of "agrees", both rooms agree |
| F1b, A1 to A3 | not refuted; `A''(0)` from a finite difference of the true area in `ℝ⁴` matches `q_J` to `3.6 × 10⁻⁶` |
| F2, N1, N2 | not refuted, recomputed from each room's JSON without the maintainer's scripts |
| The seam claim | not refuted, on the full protocol at all seven widths and five levels (§ 6) |
| Containment | not refuted |
| The stage-2 supplied grade | correct, and needed only at `W/R = 17/10`: below `πR/2`, Solver A's own essential self-adjointness already gives the direct sum |

| Defect | Effect |
| --- | --- |
| The auditor's seam check (K2) tests the bottom eigenvector, which is even in `w`, so it cannot see a missing `w`-reversal; a no-reversal seam passes it | none on validity, which rests on reading the code and on the odd test function above. K2 is not evidence of the seam's direction |
| The maintainer's seam-variant run was at level 2 and two widths only | the claim stands on the full-protocol run |
| The maintainer's scorer scores the reported triple, with the recompute flags informational | none: every flag is true |
| A line of the maintainer's run notes misread one room's table | not carried into this record |

## 6. The seam

The twist-blindness the task states is sharper than the bottom alone. On the frozen protocol, a flipped seam sign, a seam glued without the `w`-reversal, and a seam cut with zero data on one end only all leave the lowest six levels unchanged; only a seam deleted with both ends free moves them, by repeating the bottom as the second level. The reasons are those the auditor gives at stage 1. The two sheets touch only across the fiber, where the data are zero, so the sign of either sheet is free. `w ↦ −w` is a symmetry of the problem. And the unfolded sector problem is even in `ψ`, so its even and odd modes split between a free end and a fixed one.

The auditor found this at stage 1 (its check C3), and the adversarial audit confirms it on the full protocol at all seven widths and five levels: the flip is identical to `4.7 × 10⁻¹⁵` and the one-sided cut to `8.2 × 10⁻¹⁴`, while the no-reversal glue differs only through the mirrored triangle diagonals of one sheet, by `4.3 × 10⁻⁵` at the finest level and `4.9 × 10⁻⁹` on the bottom, and is identical to `4.6 × 10⁻¹⁴` once those diagonals are mirrored back (the gap shrinks only as `O(h²)` on the split degenerate levels at `W = πR/2`). For those three variants no numerical gate in the protocol can see the error, the seam diagnostic included. The checks that can are reading the code and testing the seam map with a function odd in `w`, and both found each room's seam correct. A check on the bottom eigenvector cannot, since the bottom is even in `w` (§ 5.2).

## 7. Containment record

| Item | Record |
| --- | --- |
| Rooms | four headless one-shot sessions (Solver A, Solver B, the auditor's two stages) under `CLEAN_ROOM_STANDARDS.md` § 3.4: four tools, no MCP server, no instruction file, code only through an interpreter sandboxed to the room, no network |
| Transcripts | every tool call classified against its room: no read outside a room, no network, every command through the room interpreter; two compound commands were refused by the harness and not run |
| Sources | each room's manifest names only its own files and textbook results; no withheld source appears in any transcript |
| Harness context | two rooms' manifests report a session-context line carrying the operator's account e-mail address, loaded without being asked for and not used; no return file or transcript echoes it |
| Isolation of the auditor | stage 2 saw the returns only after its stage-1 files were hashed, and never saw this note, the frozen terms or the frozen values |

## 8. Provenance

| Item | Pin |
| --- | --- |
| The packet | the hashes in the task's AS RUN section, all byte-identical to the files merged in #608 |
| The frozen terms | `43e7165a0ba67eafd722896ce493e21b896bb1f1643aa45d6ed58c39c3612cb9`, the author's span |
| Every landed room file | [`m8_14/run/rooms/SHA256SUMS`](../m8_14/run/rooms/SHA256SUMS) |
| Room model | `claude-opus-5-5`, read from each transcript's init record |

## 9. What this run does not verify

- **Matter.** It certifies the carrier surface's local area stability under fixed-edge normal variations only. It bears on no particle-stability claim, no Derrick row, no field dynamics and no spin-statistics.
- **The twist.** The variations vanish at the cone point, so the result is a lune's Dirichlet stability, which an orientable band of the same shape shares. § 6 shows the seam is invisible to the spectrum in three of four ways it can be wrong.
- **Global or least area.** Nothing about least area, and no width is selected.
- **The Tier 2 bar.** It does not meet it; the author's page counts the result toward that bar only if the conic realization and the fixed edge are adopted on independent grounds.
- **`S³/2I`.** Nothing about the widths at which the band's image embeds there.
- **The arguments as theorems.** A1 to A3, F1a and F1b are audited arguments, graded by AI auditors, never verified theorems. The results the rooms used from memory (Weyl's alternative, deficiency theory, Kato's representation theorem, the Gegenbauer equation and its weight, Weierstrass density, min-max) are named in their manifests and were not re-derived.
- **Records.** No `MODELS.md` cell moves, M8.7's gate is unchanged, and no earlier M8 verdict is reread.
