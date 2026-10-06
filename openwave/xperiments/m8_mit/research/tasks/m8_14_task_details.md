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

## AS RUN (go, 2026-10-01)

The maintainer's go. Everything above this heading is the filed pre-registration, byte for byte: its 17,802 bytes are a prefix of this file and hash to the `72651e2e…` restated in the [#608 review](https://github.com/openwave-labs/openwave/pull/608). This section records what was fixed before the rooms opened, so the posture is readable afterwards rather than reconstructed.

| Item | As run |
| --- | --- |
| Posture | OFFLINE, as the FIREWALL requires. Each room is a headless one-shot session with four tools, no network, no MCP server and no instruction file, running code only through an interpreter sandboxed to the room (`CLEAN_ROOM_STANDARDS.md` § 3.4) |
| Rooms | Solver A, Solver B and the auditor's stage 1, launched independently and in parallel. The auditor's stage 2 is a second launch, opened only after both returns and its own stage-1 commit are hashed |
| Exactness rule | [#547](https://github.com/openwave-labs/openwave/pull/547)'s: exact means a symbolic derivation, or an identification stating its precision, repeated at a second precision |
| Argument-grading rule | M8.11's with M8.13's resolution, as the auditor brief states it: one verdict per step, ESTABLISHED; ESTABLISHED, SUPPLIED with the supplied part named; GAP; or DEFECT. A passing argument is recorded as an audited argument |
| No-search rule | every computed number is reported, and nothing is tuned toward a target. Each brief carries it |
| Replacement rooms | one per room, as authorized in the verdict rule, each a clean context with the spec sheet and the solver brief alone |

### The packet, hashed at the go

| File | Who sees it | SHA-256 |
| --- | --- | --- |
| [`m8_14/spec_sheet.md`](../m8_14/spec_sheet.md) | all rooms | `0ce4e605b0f95f57b6d86cc900583b38696a9a5439cc3a1b194d640cf947f492` |
| [`m8_14/auditor_brief.md`](../m8_14/auditor_brief.md) | the auditor | `a53b66426f7e160cf5c20275141a623f91a852172d2df5fce8ae4361bd890acd` |
| [`m8_14/frozen_terms_ix.txt`](../m8_14/frozen_terms_ix.txt) | no room | `43e7165a0ba67eafd722896ce493e21b896bb1f1643aa45d6ed58c39c3612cb9` |
| `m8_14/run/brief_solver.md` | Solvers A and B | `b69b6fa0a3dea54f97eb60d3f2f014664fee610d47f8eea7eba67aec7607f886` |
| `m8_14/run/brief_auditor_s1.md` | the auditor, stage 1 | `d315153b3aea5acdc685b22e599e5961d9292a951882515c03835e8e4264faea` |
| `m8_14/run/brief_auditor_s2.md` | the auditor, stage 2 | `d90453d3ba22cb2bb5a70fe2e4355d52fdc171c4ab28e89f7826a0d00cc44f1f` |
| `m8_14/run/prompt.txt` | every room, as its opening prompt | `c91894f7d99079f3b0751b1a3b01e1612d3f2f8746bdd7fcd3851b27d778cacd` |

The spec sheet, the auditor brief and the frozen terms are byte-identical to the files merged in #608. The three briefs carry room mechanics only: tools, containment, the exactness and no-search rules, the return format, and the manifest.

**The packet, read again at the go.** The handout gate finds no frozen value and no withheld term in any packet file; its hits are the five widths the spec sheet must carry and the spec sheet's own word "conic", which describes the cone point. Its selftest fires on a planted value and a planted term. Read for structural hints, the two the review named stand as judged there (T2 names the limit-circle case and `17/10`; the auditor brief says a wrong seam sign can leave the numbers unchanged), and the three briefs add none.

### Note, 2026-10-01, at the go: the seam diagnostic

Agreed by the author in the [#608 thread](https://github.com/openwave-labs/openwave/pull/608#issuecomment-5919060713), and added here with the filed text untouched. The adjudication records each return's lowest six levels, at the finest level, against F1b's formula and against the auditor's own values. It is a diagnostic and cannot change the frozen verdict. The bottom is even across the seam, so a return with the seam deleted and both ends left free passes F2, and the second level is what separates it: the bottom repeats there.

⚠️ **Corrected during the run, before anything landed.** As first written at the go, this note followed the author's agreeing comment and named a seam glued without the `w`-reversal as a second case the second level separates. It is not one. The auditor's stage 1 found it invisible, the maintainer's own code confirms it at `W/R = 1` and `7/5` on the third level, and the adversarial audit confirms it on the full protocol (FINDINGS F4). A seam glued without the reversal, a flipped seam sign, and a seam cut with zero data on one end only all leave the lowest six levels unchanged; only the both-ends-free deletion moves the second level. For those three, the validity audit's reading of the code is the only check.

### The maintainer's own route, run before the rooms opened

| Check | Result |
| --- | --- |
| The frozen values table | parsed from this file and recomputed at 30 digits: all 21 computed entries match to the printed six decimals, and `J`'s bottom is the `−Δ` bottom minus 2 at every width. A one-unit change in a sixth decimal fires |
| The outcome rules, encoded | the frozen per-width outcome (Reproduced, Contradicted, Unresolved), F2's rule over the five widths, and the two control gates. Eight planted cases each give the outcome the frozen text assigns them, including a miss inside ten times its error estimate (Unresolved) and a miss of `9 × 10⁻⁴` relative (still Reproduced) |
| The reviewer's implementation of the frozen protocol | re-run from the #608 review at all seven widths, 16 seconds in all, and it reproduces the review's table digit for digit: F2 Reproduced at all five widths, N1 within tolerance, N2's `J` bottom `−0.222230` |
| The seam-diagnostic reference | F1b's lowest six levels at all seven widths, computed from the formula; at `W = πR/2` they are the hemisphere's 2, 6, 6, 12, 12, 12 |

## FINDINGS (run closed 2026-10-01)

Full record: [`../findings/m8_14_method_note.md`](../findings/m8_14_method_note.md). Room returns, hashed: [`../m8_14/run/rooms/`](../m8_14/run/rooms/).

**Reproduced.** Two blind rooms, offline, each ruled VALID by a blind two-stage auditor, each establishing the reduction and the full spectrum and each reproducing the frozen finite-element bottom at all five widths, with both controls passing. By the verdict rule, two valid Reproduced returns give Reproduced.

| F | Finding |
| --- | --- |
| F1 | ⭐ The frozen check reproduces in full. Both rooms reduce the band's Dirichlet problem to the lune of opening `2W/R`, derive `(ν_m + ℓ)(ν_m + ℓ + 1)/R²` with a completeness argument, and run the frozen protocol to a miss of at most `3.0 × 10⁻⁷` relative at the five widths, against `10⁻³` |
| F2 | Neither room was told the lune or the reduction. Solver A found it by unfolding the lower sheet in latitude and longitude with `u = −φ` there; Solver B by an explicit unitary onto the lune, the seam becoming its equator arc |
| F3 | Every graded step is ESTABLISHED in both rooms, one with a supplied part: in Solver A's T3, the auditor supplied that the 2D operator is the direct sum of the sector operators |
| F4 | ⚠️ **The seam is invisible in three of the four ways it can be wrong.** A flipped sign, a glue without the `w`-reversal, and a cut with zero data on one end only all leave the lowest six levels unchanged; only a seam deleted with both ends free shows, as a repeated bottom. Found by the auditor at stage 1 and confirmed by the adversarial audit on the full protocol, all seven widths and five levels. No numerical gate in the protocol sees these three errors, the seam diagnostic included; reading the code and testing the seam map with a function odd in `w` do, and both rooms pass. It corrects the author's agreeing comment and the seam note as first written at the go (AS RUN section) |
| F5 | The two rooms' finite-element numbers agree to `8.2 × 10⁻¹⁴` relative at every width and level, and with the reviewer's own implementation of the protocol in #608 to `10⁻¹²`. Three independent codes reached the same discrete spectrum, so the protocol as frozen determines its numbers |
| F6 | The protocol's error estimate `\|λ_ext − λ_finest\|` measures the finest level, not the extrapolation: it overstates the extrapolated value's actual error by two to four orders of magnitude at every width, as both rooms report |

### The adjudication

| ID | Adjudicated | How |
| --- | --- | --- |
| F1a | ✅ agrees, both rooms, ESTABLISHED, SUPPLIED | the maintainer identifies each room's map with § 1.3 of the method note: an isometry on each sheet, the seam condition becoming continuity across an interior line, the edges `θ = ±W/R`, the two sides of the fiber the two poles. The auditor graded it through the steps that use it; the adversarial audit graded it directly and named two supplied parts in each room, S-a (the value-only seam class has the same form closure as the smooth core) and S-b (the vertex-excluded problem is the lune's Dirichlet problem, points having zero capacity in two dimensions), method note § 5.2 |
| F1b | ✅ agrees, both rooms | the full spectrum with labels, complete, every step ESTABLISHED |
| F2 | ✅ Reproduced, both rooms | the frozen outcome at each width, from each room's own JSON: miss and error estimate within a relative `10⁻³`, the reading the frozen text states, `p` between 1.9967 and 2.0004. Under an absolute reading `W/R = 1/4` and `1/2` would be Unresolved |
| A1 | ✅ audited argument | exact density, no edge, seam or vertex term, extension by continuity |
| A2 | ✅ audited argument | indicial roots `±ν_k`, limit-circle iff `ν_k < 1`, so every sector limit-point for `W ≤ πR/2` |
| A3 | ✅ audited argument | `J`'s bottom `(α₀ − 1)(α₀ + 2)/R² > 0`, index 0, nullity 0, at every embedded width |
| N1 | ✅ passes | `2.000000050` against 2, a miss of `2.5 × 10⁻⁸` relative |
| N2 | ✅ passes | `J`'s bottom `−0.2222301`, about 1,600 times its error estimate |
| Verdict | ⭐ **Reproduced** | two valid Reproduced returns |

### DEVIATIONS LOG

| # | Deviation |
| --- | --- |
| 1 | The seam note as written at the go named a glue without the `w`-reversal as a case the second level separates, taking the author's agreeing comment as given. It is invisible at every level (F4). Corrected in place before anything landed, with the correction marked there |
| 2 | The handout gate's selftest printed FAIL on the spec sheet: the value it plants is a width the sheet must carry, so it was already a hit and the arm could not fire there. Rerun on a packet file with no hits, it passed, and a second arm planting a frozen decimal fired. A maintainer-tool arm dead on its arena, the same shape `PR_REVIEW_STANDARDS.md` D12 names |
| 3 | The maintainer's frozen-table parser stopped with an error on the `π/2` label at its first run, and was fixed before it printed any result. A loud failure, recorded for completeness |
| 4 | Two rooms' manifests report a session-context line carrying the operator's account e-mail address, loaded without being asked for and not used. No return file or transcript echoes it. A containment observation for `CLEAN_ROOM_STANDARDS.md` § 3.4: the restricted headless session still receives that context |
| 5 | The maintainer's seam-variant run was at level 2, not the level 3 its docstring names, and at two widths only. The seam claim is recorded from the adversarial audit's run on the full protocol |
| 6 | The auditor's seam check K2 tests the bottom eigenvector, which is even in `w`, so it cannot see a missing `w`-reversal. Found by the adversarial audit. Validity rests on reading the code and on a test of the seam map with a function odd in `w` |
| 7 | The maintainer's scorer scores each room's reported `p`, `λ_ext` and `err`, and reports their agreement with its own recomputation from the levels without gating on it. Every agreement flag is true here, to `2.6 × 10⁻¹⁵`; the next run should gate on it |
| 8 | A run note of the maintainer's recorded that one room's report omits a width from a table. False: the maintainer's own read had filtered out that row. Found by the adversarial audit and not carried into the record |

## TASK REVIEW (2026-10-01)

Task Duration: 01:11 (from 10:26 to 11:37)
Usage Cap Triggered: NO

| Result | |
| --- | --- |
| ✅ | Reproduced under the frozen verdict rule: two blind rooms, offline, both VALID, each establishing the reduction and the full spectrum and reproducing the frozen bottom at all five widths, the largest miss `3.0 × 10⁻⁷` relative |
| ✅ | N1 and N2 pass: `2.000000050`, and `J`'s bottom `−0.2222301` at about 1,600 times its error estimate |
| ✅ | F1a, F1b and A1 to A3 ESTABLISHED in both rooms, with supplied parts named: one by the auditor in Solver A's T3, and S-a and S-b in the reduction by the adversarial audit |
| ✅ | The adversarial audit refuted none of eight claims, with its own code |
| ✅ | Containment clean: no read outside a room, no network, no withheld source |
| ⚠️ | The seam note as first written at the go named the no-reversal glue as visible in the second level; it is invisible at every level. Corrected in place (F4) |
| ⚠️ | The auditor's seam check K2 cannot see a missing `w`-reversal; validity rests on reading the code and on a test with a function odd in `w` |
| ⚠️ | Three defects in the maintainer's own instruments, none touching the verdict (deviations 5, 7 and 8) |

Issues: none blocking. The filed pre-registration is an exact byte prefix of this file, so its own lint findings stay as filed.

Deviations from plan: eight, logged above. Four were in the maintainer's own instruments or notes, one in the auditor's checks, none in the rooms' mathematics.

Action needed: none outstanding on the platform side. The dated Result on the carrier page is the author's, per the records list.

### Findings

The conic carrier is a strictly stable critical point of area with its edge held fixed at every embedded width, as an audited argument: two blind rooms found the reduction to a lune without being told it, derived its full Dirichlet spectrum, and reproduced the frozen bottom at all five widths. The seam is invisible to every numerical gate of the protocol in three of the four ways it can be wrong, so reading the code is the only check on it.

### Research docs created/updated

- [`tasks/m8_14_task_details.md`](m8_14_task_details.md) (this file): AS RUN, FINDINGS F1-F6, the adjudication, the deviations log
- [`findings/m8_14_method_note.md`](../findings/m8_14_method_note.md): the author-facing record
- [`m8_14/run/`](../m8_14/run/): the room briefs and prompt, the hashed room returns, [`adversarial_audit.md`](../m8_14/run/adversarial_audit.md), the maintainer's scoring scripts
- [`m8_roadmap.md`](../m8_roadmap.md), [`m8_theory_canonical.md`](../m8_theory_canonical.md), [`__M8_model_briefing.md`](../../__M8_model_briefing.md)

## Dated note, 2026-10-01: after the verdict

The author recorded the dated Result that this task's records list asks for, on the carrier page: [The Fixed-Edge Stability Check](https://github.com/dmobius3/mode-identity-theory/blob/df9b3f73a6c37add8573c9b941358369494e9fa5/files/framework/files/working/files/fixed-edge-check.md), at MIT commit `df9b3f7`. It scores the run against the frozen terms, quotes this task's declared reading of "agrees", and reaches Reproduced. It recomputes each width from the rooms' hashed returns, with a scoring script whose planted arms include the absolute reading of the tolerance.

Two notes from [#608](https://github.com/openwave-labs/openwave/pull/608)'s review, which the author agreed to make on the next touch, close here. The filed pre-registration's run-before-write line names [#606](https://github.com/openwave-labs/openwave/pull/606) as the only landing since M8.13's adjudication on 2026-09-22; [#589](https://github.com/openwave-labs/openwave/pull/589) (2026-09-22), M8.13's author package after its verdict, [#595](https://github.com/openwave-labs/openwave/pull/595) (2026-09-24), which completed M8.12's correction record, and [#598](https://github.com/openwave-labs/openwave/pull/598) (2026-09-25), a living-docs sync, also landed in that window. None governs an experiment that had not yet run, so the condition held. Its Sources of record row links the carrier page at `main`; the frozen terms are lines 288 to 312 at commit `7c402c5`, as the row says, [pinned here](https://github.com/dmobius3/mode-identity-theory/blob/7c402c53cfd79a551d914521c0507937c00f4060/files/framework/files/working/files/projective-carrier.md?plain=1#L288-L312). The next row links `first-eigenvalue.md`, for Lemma 3.3 and Remark 1.3, at `main` as well; the page is unchanged since that commit, [pinned here](https://github.com/dmobius3/mode-identity-theory/blob/7c402c53cfd79a551d914521c0507937c00f4060/files/framework/files/bedrock/files/first-eigenvalue.md). The filed bytes above are unchanged.

Since the run, the author's carrier page proves analytically what the method note lists as not verified under global or least area: with the edge fixed as parametrized and the core generating, the conic band has the least area among continuous finite-piecewise-smooth competitors, at every width, and ties only with maps that cover it once ([§IX, Proposition 7](https://github.com/dmobius3/mode-identity-theory/blob/33e9dea0fdcc629748676cf400e34acfb832a2fd/files/framework/files/working/files/projective-carrier.md#ix-fixed-edge-stability), MIT commit `33e9dea`, by Crofton's formula). The Lipschitz class stays open. This is not part of the run, and no verdict, `MODELS.md` cell or gate moves.

## Dated note, 2026-10-05: the Lipschitz class

The 2026-10-01 note's last paragraph left the Lipschitz class open. The author's carrier page has since closed it: [§IX, Proposition 8](https://github.com/dmobius3/mode-identity-theory/blob/fd2962fc8c96bdb2985e9e6eac61aedbc463b646/files/framework/files/working/files/projective-carrier.md#ix-fixed-edge-stability), proved on [The Lipschitz Comparison](https://github.com/dmobius3/mode-identity-theory/blob/fd2962fc8c96bdb2985e9e6eac61aedbc463b646/files/framework/files/working/files/lipschitz-comparison.md) at MIT commit `fd2962f`, extends Proposition 7 to Lipschitz carrier-class competitors at every embedded width, equality case included. With the edge fixed as parametrized and the core generating, the conic band has the least area, and a competitor ties it only by covering it once. This is not part of the run, and no verdict, MODELS.md cell or gate moves.
