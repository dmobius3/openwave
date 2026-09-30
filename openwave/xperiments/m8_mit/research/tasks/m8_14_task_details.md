# M8.14: FIXED-EDGE STABILITY OF THE CONIC CARRIER, a blind run of the check the author registered and froze

Pre-registration, for review. It governs nothing until the maintainer's go.

## TASK PLANNING

### The question

The [canonical](../m8_theory_canonical.md)'s arena now carries a projective layer, `ℝP³ = S³/{±I}`, and a Möbius carrier realized through it, whose leading explicit realization is the conic band `M(W)`. The author's carrier page derives that this band is a strictly stable critical point of area when its edge is held fixed, at every embedded width, and on 2026-09-27 it registered and froze an independent check of that result, which has not run.

This task asks one question: **is the second variation of area of `M(W) ⊂ ℝP³(R)`, over normal graphs that vanish on its edge, positive definite for every `0 < W < πR/2`?** It runs the frozen check in two blind rooms, by an analytic rederivation and a finite-element cross-check, and scores it by the author's frozen outcome rules.

### Standing

- **The ruling it files under** is [#512](https://github.com/openwave-labs/openwave/discussions/512#discussioncomment-18415036)'s, as M8.10 to M8.13 did. The task is one pre-registration within the § 12.2 word budget. Nothing further lands on it until it has run. A claim that its frozen text is defective is reproduced by the maintainer with independent code before it is ratified.
- **Run-before-write** (`PR_REVIEW_STANDARDS.md` § 12.2) is met: M8.13's pre-registration ([#583](https://github.com/openwave-labs/openwave/pull/583)) ran and was adjudicated on 2026-09-22, and since then only a living-docs sync ([#606](https://github.com/openwave-labs/openwave/pull/606)) has landed.
- **Not a dynamics programme, and not the closed chassis.** It evolves no field and proposes no Lagrangian. Its solver is a finite-element eigensolve of one surface operator, not the spectral chassis in `2I`-symmetric harmonics that the C2 verdict closed.
- **Why run it.** It runs a registered, frozen check blind, exactly as frozen. Once the reduction holds, the rest is a classical computation, so the run's content is the reduction and the second variation; and N1's classical endpoint makes it a clean calibration of the two-room pipeline.
- **Word budget.** About 3,000 words, inside § 12.2's cap.
- **ID.** M8.14 is the next free ID in creation order, allocated by the row this pull request creates (`ROADMAP_STANDARDS.md` § 6). If a row taking M8.14 lands first, whoever merges second renumbers this one.

### What this task is not

- **Surface, not matter.** It certifies the carrier surface's local area stability only. It makes no particle-stability claim and bears on no Derrick row, no field dynamics and no spin-statistics.
- **Twist-blind.** The variations vanish at the cone point, the one place the twist lives, so the twisted and untwisted Dirichlet problems coincide. What it certifies is a lune's Dirichlet stability, which an orientable band of the same shape would share. It is not a Möbius-specific stability.
- **Local, not global.** It says nothing about least area. It selects no width: the result is one stable band per width, and `W` stays unselected.
- **No bar.** It does not meet the author's Tier 2 bar. The author's page counts the result toward that bar only if the conic realization and the fixed edge are adopted on independent grounds, and both are open there.
- **No `S³/2I` claim.** It says nothing about the widths at which the band's image embeds in `S³/2I`.
- **Records.** No `MODELS.md` cell moves, M8.7's gate is unchanged, and no earlier M8 verdict is reread.

### Ownership and run format

The author freezes the claims below and supplies the frozen terms, which are already public. The run is maintainer-run.

1. **Solvers A and B**, in separate offline rooms, each receive only the spec sheet and commit their returns.
2. **The auditor**, blind like the rooms, works in two stages and holds no expected value at either. At stage 1 it receives its brief and the spec sheet only, recomputes T4 by a different discretization, and commits its scripts, data and table. At stage 2 it receives both returns with their scripts, rules each return's validity, compares, and grades each room's arguments step by step.
3. **Adjudication**, by the maintainer, compares what the rooms established with the claims below, and applies the frozen outcome rules, the control gates and the verdict rule.

### Sources of record

| Source | What it supplies |
| --- | --- |
| MIT [`projective-carrier.md`](https://github.com/dmobius3/mode-identity-theory/blob/main/files/framework/files/working/files/projective-carrier.md) §IX, at commit `7c402c5` (lines 288 to 312 there) | the frozen terms: the rederivation, the cross-check protocol, the tolerance and the outcomes. SHA-256 of the span as the page's header defines it, recomputed from the bytes at that commit: `43e7165a0ba67eafd722896ce493e21b896bb1f1643aa45d6ed58c39c3612cb9` |
| the same page, §I, §V and §VI | the placement `M(W) = (L ∪ −L)/{±I}` in `ℝP³`, the conic band's costs, and `J = −Δ − 2/R²` (Lemma 4) |
| MIT [`first-eigenvalue.md`](https://github.com/dmobius3/mode-identity-theory/blob/main/files/framework/files/bedrock/files/first-eigenvalue.md) | the band `M(W)`, Lemma 3.3 (limit-point at the cone for angular exponent at least 1, limit-circle below) and Remark 1.3 |
| [M8.1](m8_1_task_details.md) (2026-07-21) | the pillar's Neumann spectrum, verified blind; withheld from the rooms (below) |
| This task | F1, F2, A1 to A3, N1, N2 |

## SETTING

`S³(R) ⊂ ℝ⁴` and `ℝP³(R) = S³(R)/{±I}`. On the great sphere `Σ = S³ ∩ {x₄ = 0}`, write `X(y, w) = R(cos(y/R) cos(w/R), cos(y/R) sin(w/R), sin(y/R), 0)`. On `0 ≤ y ≤ πR/2` the band is `X(y, w)`, and on `πR/2 < y ≤ πR` it is `X(y, −w)`, with `|w| ≤ W`. The segment `y = πR/2` is the single point `p_c`, and the quotient by `±I` identifies `(0, w)` with `(πR, −w)`. The identification is the antipodal map. The metric is `dy² + cos²(y/R) dw²`, the area weight is `|cos(y/R)|`, and the smooth locus of `M(W)` is totally geodesic in `ℝP³(R)`, with `Ric(ν, ν) = 2/R²`.

In the strip coordinates the normal is the constant vector `ν = e₄`, the normal of the great sphere `x₄ = 0`. The antipodal map carries `e₄` to `−e₄`, so a normal graph `φν` is well defined on `ℝP³` exactly when `φ(πR, −w) = −φ(0, w)`. Area is varied along such graphs that vanish on the edge `w = ±W`: the completion, in the form norm, of the smooth ones that vanish near `p_c`, with the conic type at `p_c` held fixed. The rooms' [spec sheet](../m8_14/spec_sheet.md) states the same object by formula, with every hypothesis inline.

## THE FIREWALL

- **What the rooms see.** Each solver sees only the [spec sheet](../m8_14/spec_sheet.md). The auditor sees its [brief](../m8_14/auditor_brief.md) and the spec sheet at stage 1, and both returns with their scripts only at stage 2.
- **What is withheld from all three.**
  - this document;
  - the author's pages, the first-eigenvalue paper and the author's scripts;
  - M8.1's task document, method note, scripts and data, which carry the same `α₀(α₀ + 1)/R²` as the pillar's wide Neumann branch;
  - the canonical's §3 rows;
  - the frozen-terms file, which states the target formula in its tolerance clause;
  - every value in the frozen-values table below.
- **Network: offline.** The answer is public in the author's repository, so an online room could find it. A room that reaches any withheld source returns an invalid run (below).

## DISCLOSURE

- **The result is public.** It and its derivation have been public since 2026-09-27, with the frozen terms.
- **Author-side checks.** The author's `identities.py` checks the second-variation identity symbolically (its I9) and the trivialization on the `2/R²` mode (I3). It computes no Dirichlet spectrum, and the author has not run the frozen check.
- **The formula is not new to the platform.** M8.1 verified blind the pillar's Neumann spectrum, whose wide branch is `α₀(α₀ + 1)/R²`: the same function this task's Dirichlet bottom takes, because the Dirichlet transverse sectors share the Neumann problem's nonconstant angular eigenvalues.
- **The author's interest.** The author wants §IX verified, and would cite it only as an independent check of a surface result, under the fences above.

## CANDIDATE PRE-REGISTERED CLAIMS

`R = 1` in every reported number. `α₀ = πR/(2W)`.

### Group F: the frozen check, scored by its own outcome rules

| ID | Claim | Pass | Fail |
| --- | --- | --- | --- |
| F1a | The reduction: the band's Dirichlet problem is unitarily equivalent, as a whole or sector by sector, to the Dirichlet problem on a spherical lune of opening `2W/R` | the room's argument establishes the equivalence, by the author's unitary or by separating variables and unrolling the seam; the maintainer identifies it at adjudication, in whatever form the room states it and the auditor grades it | an established argument reaching a different equivalence |
| F1b | The spectrum: `(ν_m + ℓ)(ν_m + ℓ + 1)/R²`, with `ν_m = mα₀`, `m ≥ 1`, `ℓ ≥ 0` | the room's derivation of the full spectrum is established | an established argument reaching a different spectrum |
| F2 | The finite-element bottom of `−Δ`, by the frozen protocol at `W/R ∈ {1/4, 1/2, 1, 7/5, 3/2}`, matches `α₀(α₀ + 1)/R²` within the frozen tolerance | the frozen "Reproduced" | the frozen "Contradicted"; anything else is "Unresolved" |

**The reading of "agrees", declared before the run.** The frozen text asks the runner to derive "step 2's reduction and the lune's Dirichlet spectrum", and its Contradicted clause fires when "the rederivation disagrees". Both parts stay required. F1a and F1b are graded as arguments under M8.11's rule.
- **Agrees:** both are ESTABLISHED, as written or with a supplied part named. The route may differ from the author's: a room is told neither the reduction nor the lune, and one that separates variables and unrolls the seam establishes the reduction sector by sector.
- **Disagrees:** an established argument reaches a different equivalence or a different spectrum. This alone makes the arm Contradicted.
- **Anything else is Unresolved for the arm.** That includes a correct formula without an established reduction, a partial spectrum such as the bottom alone or one parity family, a GAP, or a DEFECT.
- **Recorded.** The route is recorded as a diagnostic.

### Group A: additions, graded as arguments

| ID | Claim |
| --- | --- |
| A1 | For the normal-geodesic variation with speed `φ`, `φ` smooth in the core, the second variation of area is exactly `q_J[φ] = ∫(|∇φ|² − (2/R²)φ²) dA`, with no boundary or vertex term, and it extends to the whole class by continuity |
| A2 | Every Dirichlet transverse sector is limit-point at `p_c` for `W ≤ πR/2`, so no cone condition is chosen |
| A3 | From F1a, F1b, A1 and A2: `J`'s bottom in the declared class is `(α₀ − 1)(α₀ + 2)/R² > 0` at every embedded width, so index 0 and nullity 0 |

Additions never change the frozen verdict.

### Controls

| ID | Case | Known answer | What it shows |
| --- | --- | --- | --- |
| N1 | the frozen protocol at `W = πR/2` | the lune is a hemisphere, whose Dirichlet bottom is the classical `2/R²`, so `J`'s bottom is 0 | the protocol reaches a classical value at the end of the range |
| N2 | the frozen protocol at `W/R = 17/10` | for the frozen protocol, whose `H¹` mesh selects the Friedrichs (regular-branch) extension, `J`'s bottom is negative (value below); the wide band's cone is limit-circle in its lowest sector, so this is not an extension-independent value | the machinery can return "unstable" |

Neither control leans on the claim under test. N1's value is the classical first Dirichlet level of a hemisphere. N2's sign follows from domain monotonicity: the band at `W/R = 17/10` contains the band at the critical width, so its Dirichlet bottom lies below `2/R²` whatever `α₀(α₀ + 1)` says.

Flipping the seam is not a control: with Dirichlet data at `p_c` the twisted and untwisted problems share a spectrum, so it cannot fail.

### Diagnostics, not claims

| ID | What it records |
| --- | --- |
| G1 | each room's route in F1a and F1b |
| G2 | the estimated order `p` at each width |
| G3 | the ground state's shape, against the author's `±|cos(y/R)|^{α₀} cos(πw/2W)`, sign flipping across the cone |

### The verdict rule

- **Validity first.** The auditor rules each return valid or invalid before any claim is scored. An invalid return never yields Contradicted. A return is invalid if:
  - it breaks an explicit requirement of the spec sheet (the geometry, the seam, the weight, the boundary data, the discretization and protocol, the widths), whatever the numerical effect. A seam error, for instance, can be spectrally invisible here and still invalidates;
  - it used a withheld source;
  - it has an instrument defect that moves a reported bottom by a relative `10⁻³` or more.

  Only incidental defects that break no explicit requirement and move no reported bottom that much may be ruled immaterial.
- **The controls gate validity**, applied by the maintainer, since the auditor holds no expected value. A return whose N1 misses its known value by more than the frozen tolerance, or whose N2 does not give `J` a negative bottom by more than ten times its error estimate, is invalid as an instrument defect.
- **The auditor's recomputation gates validity.** If it disagrees with a return beyond the frozen tolerance, the auditor must locate the defect: in the return, which is then invalid, or in its own computation, which is then corrected and recorded. If it can locate neither, the verdict is Unresolved.
- **One replacement room per room is authorized now.** An invalid return is replaced by a fresh room: a clean context with the spec sheet alone, which never sees the first return, its audit or the auditor's numbers. The replacement is authorized before any run, so it is not a refinement added afterward. If the replacement's return is also invalid, the verdict is Unresolved.
- **Then:**
  - two valid Reproduced returns give Reproduced;
  - any valid Contradicted return gives Contradicted, reported with both returns and with `J`'s sign at the width in question, as the frozen text requires;
  - anything else is Unresolved, and the author's result then stands on its page unverified.

## FROZEN VALUES

| `W/R` | `α₀` | bottom of `−Δ`, `α₀(α₀ + 1)` | bottom of `J`, `(α₀ − 1)(α₀ + 2)` |
| --- | --- | --- | --- |
| 1/4 | 6.283185 | 45.761603 | +43.761603 |
| 1/2 | 3.141593 | 13.011197 | +11.011197 |
| 1 | 1.570796 | 4.038197 | +2.038197 |
| 7/5 | 1.121997 | 2.380875 | +0.380875 |
| 3/2 | 1.047198 | 2.143820 | +0.143820 |
| π/2 (N1) | 1.000000 | 2.000000 | +0.000000 |
| 17/10 (N2) | 0.923998 | 1.777770 | -0.222230 |

The tolerance, the mesh, the refinement count, the order estimate, the extrapolation and the outcome clauses are the frozen span's. They travel as [`frozen_terms_ix.txt`](../m8_14/frozen_terms_ix.txt), byte-identical to the span, whose SHA-256 is the one above. The maintainer adjudicates from it; neither the rooms nor the auditor see it.

## FEASIBILITY, AND THE AUTHOR'S DERIVATION

- **Finite elements.** A graded P1 mesh on a rectangle, seven widths, four refinements each. The finest level has 512 cells along the core and 128 across, seconds per width with a sparse eigensolver.
- **Analytic parts.** Separation of variables, a local limit-point computation at the cone of the kind the first-eigenvalue paper's Lemma 3.3 does, and one second-variation computation.
- **The author's derivation.** It is §IX's four steps: the exact second variation, the reduction by the pillar's unitary, the count on the lune, and the Neumann contrast. The author expects Reproduced.

## TO BE FIXED AT GO

- **Byte pins.** The spec sheet, the auditor brief, and the maintainer's room brief if there is one. The frozen-terms file is already pinned by the author's hash.
- **The exactness rule**, as [#547](https://github.com/openwave-labs/openwave/pull/547)'s: exact means a symbolic derivation, or an identification stating its precision, repeated at a second precision.
- **The argument-grading rule**, M8.11's with M8.13's resolution. The auditor applies it to the rooms' arguments for T1, T2, T3 and T5, and the maintainer maps the graded conclusions to F1a, F1b and A1 to A3. Each step gets one verdict: ESTABLISHED; ESTABLISHED, SUPPLIED, with the supplied part named; GAP; or DEFECT. A passing argument is recorded as an audited argument, never as a verified theorem.
- **The no-search rule.** Every computed number is reported, and nothing is tuned toward a target.

## DEFINITION OF DONE

- **Returns.** Solvers A and B return their derivations and finite-element tables with scripts.
- **Audit.** The auditor commits its own recomputation at stage 1. At stage 2 it returns validity rulings, the comparison, and step-by-step grades of each room's arguments.
- **Adjudication** records F1a, F1b, F2 and the verdict under the rule above, grades A1 to A3, and scores N1 and N2.
- **Records.**
  - A method note: equations first, with an equation-to-code map.
  - A new canonical §3 row, "Fixed-edge stability of the conic carrier", carrying its fences in the row: surface not matter, local not global, twist-blind, not the Tier 2 bar.
  - The roadmap Done row.
  - After the verdict, the author records a dated Result on the carrier page, scored against its frozen terms and quoting this task's declared reading of "agrees".
- **What does not move.** No `MODELS.md` cell, and M8.7's gate is unchanged.
