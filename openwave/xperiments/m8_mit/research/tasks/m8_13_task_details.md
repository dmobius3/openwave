# M8.13: AUDIT OF THE UNIQUENESS OF THE TOP-MULTIPOLE MAXIMUM AT SPIN 3, the clause of M8.12's G2 that its run left unresolved

Pre-registration, for review. It governs nothing until the maintainer's go.

## TASK PLANNING

### The question

M8.12 reproduced G2 in its value and attainment: the maximum of `r̂₆` on the unit sphere of `V₃` is `463/924`, argued in one room and certified independently by the audit, and the hexagon orbit attains it. It left G2's uniqueness clause unresolved, since neither room nor the audit settled whether any other orbit attains the maximum ([#582](https://github.com/openwave-labs/openwave/pull/582)). That clause rests on the author's argument in `S0_S3_MAXIMUM.md`, landed with the author package at [#581](https://github.com/openwave-labs/openwave/pull/581) and labelled there as not independently audited.

This task asks one question: **is the hexagon orbit the only maximizer of `r̂₆` on the unit sphere of `V₃`?** It answers it twice, by an auditor's own argument committed before it sees the author's, and by grading the author's argument step by step.

### Standing: the #512 ruling

- **Run-before-write** is met on its plain reading: M8.12's pre-registration ran at #578 and closed with its package at #581 and the correction at #582.
- **One program, one pre-registration**, within the § 12.2 budget.
- **A claim that frozen text is defective** is reproduced by the maintainer with independent code before ratification.

M8.13 is allocated here as the next free ID in creation order, under `ROADMAP_STANDARDS.md` § 6. If a row taking M8.13 lands first, whoever merges second renumbers this one, which edits this document and so moves its hash.

### What this task is not

It audits one proposition about one quartic on one seven-dimensional space. It computes no Hessian and no index, reopens none of M8.12's reproduced values, and makes no stability claim. No `MODELS.md` cell moves and M8.7's gate is unchanged.

### Ownership and run format

The author freezes the claims below and supplies the argument under audit, which is already public. One maintainer-run auditor works in two stages, as M8.11's audit did for its theorem. At stage 1 it receives only the worklist and commits its own answer and argument. Stage 2 has two parts, and each return is committed before the next part's files are handed over: at 2a the auditor receives the grading instructions for 2a and N2's variant, and grades the variant's part 3c; at 2b it receives the author's argument with its checker and log, and grades the argument step by step. Then adjudication.

### Sources of record

| Source | What it supplies |
| --- | --- |
| M8.12 ([#578](https://github.com/openwave-labs/openwave/pull/578), [#582](https://github.com/openwave-labs/openwave/pull/582)) | G2 reproduced in value and attainment, its uniqueness clause unresolved; G1, the minimum and its equality case, reproduced in both rooms |
| `research/scripts/m8_12_author/S0_S3_MAXIMUM.md` ([#581](https://github.com/openwave-labs/openwave/pull/581)) | the argument under audit, SHA-256 `7c634a33fdb0ee05d5a934395a414f6346ba26e049b07f2619961792ab725498` |
| `research/scripts/m8_12_author/check_s3_maximum.py` and its log | its checker, SHA-256 `34054633f9602e864c12d017b035e6517ebbc4d8b0ba6c043dc1ef4adc764d3a`, and its log, `37a00d34168b23200a1f07694627dc2a93dd588b9156be5c55661deffab67b1d`: the inequality exactly, the equality set only by a 20001-point sweep |
| This task | U1 to U3, N1, N2 |

## SETTING

`V₃ = ℂ⁷` with the spin-3 operators `f = (fₓ, f_y, f_z)`, time reversal `(Θu)_m = (−1)^{3−m}·conj(u₋ₘ)`, and the top multipole `ρ₆(u) = [u ⊗ Θu]₆` in the Condon-Shortley convention, fixed by `⟨3 3; 3 3 | 6 6⟩ = +1`. Write `r̂₆(u) = ‖ρ₆(u)‖²/‖u‖⁴`, as in M8.12.

The argument under audit runs in four steps. Step 1 writes `r̂₆` on the unit sphere as a fixed combination of three invariant quartics: the magnetization, the singlet-pair weight and the nematic invariant. Step 2 bounds the first two, each with its equality case. Step 3 bounds the third, `TrN̄² ≤ 171/2`, by reducing it to `λ_max(e₁fₓ² + e₂f_y² + e₃f_z²) ≤ 15/√6` over unit traceless `e`, splitting that operator into three 2×2 blocks and bounding each. Step 4 combines the three equality cases.

**Where uniqueness lives.** Not in step 4, which imports step 3's equality result and imposes step 2's two conditions. It lives in the equality cases, and step 3's is the substantive one: the tracing of where each block reaches `15/√6` on the constraint ellipse, and the passage from there to the states attaining `TrN̄² = 171/2`. Step 4 would read as correct even if step 3 had missed an equality point, so step 3's tracing is what this task grades hardest.

**One property of the equality set that the grading has to respect.** Each of its three points is reached by two of the three blocks at once, so no single block is load-bearing for the set. That is a statement about blocks, not about the note's three traced cases: the case `b = 0` with `a > 0` carries the first point for M₁ and, through the symmetry sentence, for M₂, so omitting it would lose that point. Omitting the third block's case instead leaves every point covered while the argument no longer shows that the third block contributes nothing further. N2 uses exactly that omission, and it scores whether the grader notices it.

## THE FIREWALL

The stage-1 room has no access to: this document; `S0_S3_MAXIMUM.md` and its checker and log; the author's check described below; and M8.12's task doc, method note and author package. It receives the worklist only. At 2a it additionally receives `stage2a_grading.md` and N2's variant, `step3_text.md`, and nothing else: the argument itself would expose the variant's omission by comparison, and its checker names the third block's equality point among its exact checks. At 2b it additionally receives `stage2b_grading.md`, and `S0_S3_MAXIMUM.md` at the hash above with its checker and log.

**Network posture: offline, and this is load-bearing for U1.** The author's uniqueness argument has been public in this repository since #581, so an unrestricted online room could locate the very argument under audit. If network access is enabled and that argument is reached at stage 1 from any source, U1 cannot be scored as an independent reproduction, and stage 1 must be rerun offline. In the literature, Kawaguchi and Ueda's review states the bound `TrN̄² ≤ 171/2` without proof, and Romero et al. probe the maximization numerically without settling it; external material located without exposure to the author's argument is recorded as located, and U1 still passes only on the auditor's own completeness argument. The worklist's standing request for anything looked up is the discriminator.

## DISCLOSURE

- **The argument under audit** was written on 2026-09-19 and has been public since #581, with its hash.
- **Its checker** proves the inequalities exactly and corroborates the equality set only by a 20001-point sweep of the circle. A sweep rejects candidates; it cannot show there are no others.
- **An author-side check, written for this task on 2026-09-21**, `m813_equality.py` with its log, determines step 3's equality set exactly: for each block, the equality condition is a conic whose resultant with the ellipse is not zero, so the two conics share no component and Bézout bounds their common points by four; exact factoring over `ℚ(√6)` then enumerates every candidate, and whether `15/√6` is the block's larger eigenvalue there is decided symbolically. It returns the three points below, each reached by two blocks, each a permutation of `(2, −1, −1)/√6` with a two-dimensional top eigenspace: 25 checks, 0 failures. The gate reruns it and requires its output to equal the log byte for byte. It does not grade the note's own tracing. It is withheld from the auditor at both stages, so that the audit's route is its own, and it lands with the author package after the verdict.
- **Review units** read the argument during the M8.12 redlines. One checked the inequalities and both sums of squares, and stated that it did not check the exhaustiveness of the tracing.
- **The author's interest**: the answer settles, at spin 3, a question the author's fourth bedrock paper records as open, and the author intends to cite it only after this audit.

## CANDIDATE PRE-REGISTERED CLAIMS

### Group P: parents, reproduced rather than new

| ID | Claim | Standing | Pass condition | Fail condition |
| --- | --- | --- | --- | --- |
| P1 | The maximum of `r̂₆` on the unit sphere of `V₃` is `463/924`, and the hexagon orbit attains it | M8.12 G2, value and attainment, #582 | the stage-1 argument reaches the same value | a different value |
| P2 | The minimum is `1/924`, attained exactly on coherent states | M8.12 G1, both rooms | the stage-1 argument reaches the same value and set | a different value or set |

### Group U: uniqueness

| ID | Claim | Standing | Pass condition | Fail condition |
| --- | --- | --- | --- | --- |
| U1 | The maximizers of `r̂₆` on the unit sphere of `V₃` are exactly the hexagon orbit, the orbit of `(v₃ + v₋₃)/√2` under rotations and phase | new at run level; the clause M8.12 left unresolved | the auditor's own stage-1 argument, committed before it sees the author's, establishes the complete set, states why it is complete, and says for each case what fails if it is omitted | another maximizer is found, or the stage-1 argument does not establish completeness |
| U2 | The route of `S0_S3_MAXIMUM.md` steps 1 to 4, with any completion the auditor supplies named explicitly, establishes the maximum and its equality set, with step 3 graded in five parts: the block decomposition; each block's bound; the equality set of `λ_max = 15/√6` on the ellipse; the passage from `‖Q‖ = 15/√6` to a top eigenspace of an extremal operator; and that eigenspace as `span{\|3,3⟩ₙ, \|3,−3⟩ₙ}` | the author's argument, public since #581 | every step and every part of step 3 ESTABLISHED, as written or with the supplied part named | any GAP or DEFECT, recorded with its location |
| U3 | G2's uniqueness clause is recorded as an audited argument | derived from U1 and U2 | U1 and U2 both pass, as M8.11 required of its theorem | either fails; a partial verdict then records which one held |

### Controls

| ID | Claim | Standing | Pass condition | Fail condition |
| --- | --- | --- | --- | --- |
| N1 | The stage-1 return proves the minimum and its complete equality set, the coherent states, where the answer is known | known: M8.12 G1 | the known value and set, with a completeness argument | a different set, or no completeness argument, which would mean the stage-1 procedure cannot certify a complete equality set even where the answer is known. It is not a test that the maximum's method transfers: bounding the invariant form term by term reaches the maximum but not the minimum |
| N2 | Given at 2a, before it sees the argument or its checker, a copy of step 3 whose equality tracing omits the third block's case, although its conclusion survives, the grader notices the omission at the equality-set part | new | the 2a return explicitly notes that the tracing sentence omits the third block's equality case, and addresses it. N2 is scored on that content, not on the label: the case may appear in what the grader supplies, in a GAP note, or in a derivation it writes out, and whatever else it supplies, including a branch analysis, neither earns nor costs the pass. The scale-correct grade is ESTABLISHED, SUPPLIED, with the case supplied. A GAP, or an ESTABLISHED that writes out a derivation from the third block's own bound, also passes N2, since each shows the omission was noticed, but each is recorded as a scale misapplication: the case can be completed, so it is not a GAP, and a derivation the grader writes is a supplied part, so it is not ESTABLISHED as written | the 2a return accepts the complete equality set without addressing the third block's equality case, whatever label it gives, which would mean the grading reads conclusions rather than arguments |

### Diagnostics, not claims

| ID | What it records |
| --- | --- |
| D1 | Each equality point is reached by exactly two blocks, so the top eigenvalue `15/√6` has multiplicity two there, one vector from each block, which is the two-dimensional span in U2's last part. Recorded as a consistency property, not adjudicated |

## FROZEN VALUES

| quantity | value | source |
| --- | --- | --- |
| the maximum | `463/924` | `S0_S3_MAXIMUM.md`; M8.12 G2 |
| the nematic bound | `TrN̄² ≤ 171/2`, equivalently `λ_max ≤ 15/√6` over unit traceless `e` | `S0_S3_MAXIMUM.md` step 3 |
| the constraint ellipse | `(2/3)a² + 8b² = 1`, with `e₃ = 2a/3` and `e₁ − e₂ = 4b` | step 3 |
| the equality points on it | `(√6/2, 0)`, `(−√6/4, √6/8)` and `(−√6/4, −√6/8)` | `m813_equality_log.txt` |
| which blocks reach each | the first by M₁ and M₂; the second by M₁ and M₃; the third by M₂ and M₃ | same |
| the same points as `√6·e` | the three permutations of `(2, −1, −1)` | same |
| the top eigenspace at each | two-dimensional; at `√6·e = (−1, −1, 2)` it is `span{v₃, v₋₃}` | same |
| the maximizing set | the orbit of `(v₃ + v₋₃)/√2` under rotations and phase | `S0_S3_MAXIMUM.md` step 4 |
| the minimum and its set | `1/924`, on coherent states | M8.12 G1 |

## FEASIBILITY, AND THE AUTHOR'S DERIVATION

The argument is one page, its checker runs in under a minute, and the author-side check above runs in about a second. The auditor's stage 1 is a genuine derivation, and it may not close: M8.12's solver A argued the maximum's value and stated explicitly that it had not settled uniqueness. U3's rule makes that outcome a partial verdict rather than a failure of the task.

**The author expects part 3c to grade ESTABLISHED, SUPPLIED, not ESTABLISHED as written.** The note states the sign conditions of M₁'s equality cases, `b = 0` with `a > 0` and `a = −2b` with `b > 0`, without deriving them from the squaring branch `a + 6b > 0`, and it does not show that the other branch, `a + 6b ≤ 0`, contains no equality point. Each of those case lines meets the ellipse twice, and only one point of each pair is an equality point. A grader following the scale should supply that branch analysis. That is still a pass, and it is recorded as one the auditor completed. N2's variant inherits the same omission, so its part 3c is expected to need that branch analysis and the third block's equality case supplied: the two texts are expected to earn the same label, and a matching label is not a failure of the control.

## TO BE FIXED AT GO

- Every file the auditor sees, at each stage, byte-pinned: the worklist; the grading instructions for 2a and 2b; N2's variant; `S0_S3_MAXIMUM.md`; `check_s3_maximum.py` and `check_s3_maximum_log.txt`; and the maintainer's room brief, if there is one. The pre-registration's own instruments are pinned in the amendment below, apart from the two logs, which are posted with it.
- The exactness rule, as #547's: exact means a symbolic derivation, or an identification stating its precision, repeated at a second precision.
- The argument-grading rule, M8.11's with its overlap resolved: one verdict per step, ESTABLISHED, ESTABLISHED, SUPPLIED with the supplied part named, GAP or DEFECT, and a passing argument is recorded as an audited argument, never as a verified or proven theorem.

## AMENDMENT (2026-09-22)

The maintainer approved and merged this pre-registration at [#583](https://github.com/openwave-labs/openwave/pull/583#pullrequestreview-5278860214), with two questions and three notes, and pinned its seven instruments by hash in the review. This section records what changed in response. No frozen value moved, and no claim was added or removed.

**Where this sits against the go.** The go was given at 2026-09-22 10:11 EDT and this amendment was written at 10:42 EDT, so it lands after the go rather than before it. The author had no way to know, because the go is given in the maintainer's terminal and is never public. It is accepted, and the record states why it costs nothing rather than leaving a later reader to compare timestamps. The stage-1 packet is byte-identical across the amendment: `worklist.md` is untouched and `step3_text.md` is a pure rename, so the stage-1 room's inputs are the registered ones. No stage-2 file had been handed over. Everything this amendment changes therefore takes effect before it is used, and **the go for stage 2 dates from this merge**.

- **N2 is sequenced, so it cannot be passed by comparing texts.** Handed over together with the argument, the variant differs from it in one clause, so its omission could be found by comparison rather than by reading the variant as an argument. Stage 2 now runs in two parts, each return committed before the next part's files are handed over: at 2a the variant alone, at 2b the argument with its checker and log. The checker waits for 2b as well, since it names the third block's equality point among its exact checks. The variant is renamed `step3_text.md`, byte-identical, so that its file name does not announce it as a variant, and the grading instructions are split to match, with the four verdicts and their rules unchanged.
- **N2 is scored on content, not on the label.** Both texts are expected to earn ESTABLISHED, SUPPLIED, so the label cannot tell them apart. N2's pass condition now scores whether the third block's equality case is addressed, and the pre-registered expectation names what each text is expected to need.
- **The firewall names only what the repository resolves.** An exclusion that named no file in the repository is removed. Stage 1 was already a whitelist, receiving the worklist only.

**The instruments.** Of the seven files the review pinned, two are unchanged and match it byte for byte:

| file | SHA-256 |
| --- | --- |
| `m813_equality.py` | `3577872c112cefe13413eab3a706e834314e635129f7a914f5488d43f3cc090d` |
| `m813_equality_log.txt` | `a30b8198d329d2a038ea3253fce2253b007954db0581fa4b71aeac5d12de0354` |

The gate, its mutation suite and its inventory changed, because this amendment added gates and arms for the sequencing, the scoring rule and the firewall. These hashes supersede the review's for those three files only:

| file | SHA-256 |
| --- | --- |
| `m813_prereg_build.py` | `eb0a440887bb778739464d6fed7bcc350038170f024d7820840a8d1480b66dd6` |
| `m813_prereg_arms.py` | `f068c3d27eac293f424df7069d8559c35e8ee8e224964b9e62938a7206c128b0` |
| `gate_inventory.txt` | `71e9a8887b077c382af0290cfc191f2f455298fb2284326e378619e813b78a6c` |

The gate's log and the mutation suite's log are outputs of runs that read this document, so their hashes cannot be stated here without changing them. They are posted with this amendment, and supersede the review's for those two files only.

## DEFINITION OF DONE

The auditor returns its stage-1 answer and argument, then its 2a grade of N2's variant, then its 2b grades of the argument. Adjudication records U1, U2 and U3, and G2's uniqueness clause moves from unresolved to an audited argument only if U3 passes. If U2 passes with any part graded ESTABLISHED, SUPPLIED, G2's status records what was supplied, for instance "audited argument, with the branch analysis of step 3c supplied by the auditor", so the record says whose argument it was, which is the lesson #582 exists for. No stability claim, no `MODELS.md` cell, and no change to M8.7's gate.

## ADJUDICATION (2026-09-22, maintainer)

The go was given 2026-09-22 10:11 EDT. Stage 1 ran against the registered
`worklist.md`; the [#585](https://github.com/openwave-labs/openwave/pull/585)
amendment merged at `42e65f79` before any stage-2 file was handed over, and
stage 2 ran in its two parts under it.

### Verdicts

| ID | Verdict | Basis |
| --- | --- | --- |
| P1 | ✅ PASS | The stage-1 argument reaches `463/924` by a route sharing nothing with the author's, and the hexagon attains it |
| P2 | ✅ PASS | `1/924` on the coherent states, same value and same set |
| U1 | ✅ PASS, with one asserted sub-step recorded below | The stage-1 argument, committed and hashed before the room saw anything of the author's, establishes the complete maximizer set, states why it is complete, and answers "what fails if omitted" for every case |
| U2 | ✅ PASS | Ten verdicts, no GAP and no DEFECT. Eight of ten ESTABLISHED, SUPPLIED; `M₂`'s bound and step 2 ESTABLISHED as written |
| U3 | ✅ PASS | U1 and U2 both pass, so G2's uniqueness clause moves from unresolved to an **audited argument** |
| N1 | ✅ PASS | The minimum and its complete equality set, with a completeness argument (the strict ordering `n₆ < n₄ < n₂ < n₀` plus `p₆ = 1 ⟺ coherent` by a Cauchy-Schwarz/Vandermonde equality analysis) |
| N2 | ✅ PASS, scale-correct | The 2a return states explicitly that the text never traces `M₃`, supplies its locus, and lands on ESTABLISHED, SUPPLIED, the label the pre-registration names |
| D1 | recorded | Each equality point carries `15/√6` with multiplicity exactly two, one vector from each of the two blocks reaching it. Confirmed by three independent routes |

### What was supplied, per the Definition of Done

G2's uniqueness clause is recorded as **an audited argument, with the branch
analysis of step 3c, `M₂`'s equality locus, the whole of step 3d, and the step
3e eigenspace computation supplied by the auditor.**

The full supplied list, on which the maintainer's independent grading (frozen
at 10:39:42, before stage 2b existed) and the auditor's agree:

| Part | Supplied |
| --- | --- |
| 1 | the `U(1)` type restriction behind the dimension count; invariance of each of the four functions |
| 3a | `‖Q‖ = max_E Tr(QE)`; `Tr(N̄E) = ⟨u, A_E u⟩`; the equivariance licensing "rotate `E` to diagonal"; the max-swap; the block basis |
| 3b `M₁` | `max(a+6b) = √6` on the ellipse, which is what validates the first squaring |
| 3b `M₃` | that `λ_max = −2a + √(30−16a²)` is the on-ellipse form; `max\|a\| = √(3/2)`; that the radicand stays positive |
| 3c | the derivation of the sign conditions; the `a + 6b ≤ 0` branch; `M₂`'s equality locus; that `M₃`'s locus meets the ellipse at all |
| 3d | the `E* = Q/‖Q‖` chain, stated nowhere in the note |
| 3e | the eigenspace computation and its multiplicity two |
| 4 | `\|a₀₀\|² = 4\|α\|²\|β\|²/7` on the span |

⚠️ **One asserted sub-step in the auditor's own stage-1 argument, which it
found and reported itself at 2b.** Stage 1 justified "a Hermitian form on
`Sym²` is determined by its values on the Veronese" with "both sides have `28²`
real dimensions and the map is injective". Equality of dimensions does not give
injectivity. The claim is true and the auditor supplied the one-line proof at
2b. **It was already independently closed on the maintainer side before stage
2b ran**: `poly_identity.py`, hashed into the frozen grading at 10:39:42,
proves the identity it underwrites as an exact polynomial identity in 14 real
variables with no dimension count and no fit, and stage 2b did not finish until
12:05. U1 therefore passes, and the record says whose argument closed which
part, which is the lesson [#582](https://github.com/openwave-labs/openwave/pull/582) exists for.

### The equality set now has four independent derivations

| Route | Method |
| --- | --- |
| author | resultant against the ellipse, Bézout bound, exact factoring over `ℚ(√6)` |
| maintainer | per-block solve of `det(M − (15/√6)I) = 0` against the ellipse with a symbolic larger-eigenvalue test |
| auditor | reduce each conic modulo the ellipse, factor completely over `ℚ` |
| adversarial audit | the same reduction plus Cauchy-Schwarz positivity of every factor, giving `λ*I − M_k ⪰ 0` on the whole ellipse |

All four return `(√6/2, 0)`, `(−√6/4, √6/8)`, `(−√6/4, −√6/8)`, the three
permutations of `(2, −1, −1)/√6`, and no fourth point.

⭐ **The auditor's factorization carries a fact none of the other treatments
states:** `M₃`'s reduced form is `(15/2)(4u+1)²`, a perfect square. `M₃` is
**tangent** to the bound, which is exactly why it contributes no new point.

### The adversarial audit

An independent offline room was given the conventions and eight claims to
refute, with no access to any return and no statement of the maintainer's
conclusions. **It refuted nothing**, and it closed both completeness claims
with exact certificates rather than exhausted searches, which is the only basis
on which a completeness claim should be accepted.

Two findings it raised that the claims did not mention, both verified here:

1. ⭐ **Saturating the nematic bound is necessary for both extrema and
   sufficient for neither.** On the whole plane `span{v₃, v₋₃}`,
   `TrN̄² = 171/2` **identically** (`N̄ = diag(3/2, 3/2, 9)`), while
   `r̂₆ = 2\|α\|²\|β\|² + 1/924` sweeps the entire achievable interval: `1/924`
   at `β = 0`, `463/924` at `\|α\| = \|β\|`, and `7417/23100` strictly between at
   `α = 1, β = 1/2`. So step 3 alone constrains `r̂₆` not at all; `⟨f⟩ = 0` and
   `Θu ∝ u` do the work in step 4. This sharpens what step 4 is for.
2. **"`M₂` = `M₁` at `−b`" is true only up to basis orientation.** In the
   natural parity-adapted basis the block is `[[5a, +√60 b],[+√60 b, −3a−12b]]`,
   which is `diag(1,−1)`-conjugate to `M₁(−b)`, with identical characteristic
   polynomials. Nothing downstream moves, since only spectra are used. **Both
   rooms found this independently**, and both declined to call it a DEFECT for
   the same reason: the note never fixes a basis.

### Instruments, and the equation-to-code map

| Claim | File | Result |
| --- | --- | --- |
| conventions, weight-state spectra, block split, equality set, top eigenspaces, step 4 | [`run/scripts/check_maint.py`](../m8_13/run/scripts/check_maint.py) | 69 passed, 0 failed, 2 deliberately-dead arms |
| step 1 as an exact polynomial identity in 14 real variables | [`run/scripts/poly_identity.py`](../m8_13/run/scripts/poly_identity.py) | residual `0`; arm `1/198 → 1/200` leaves 257 terms |
| the auditor's stage-1 claims, recomputed | [`run/scripts/verify_room.py`](../m8_13/run/scripts/verify_room.py) | 25 passed, 0 failed |
| shared spin-3 algebra | [`run/scripts/spin3.py`](../m8_13/run/scripts/spin3.py) | |

Reproduce: `cd run/scripts && python3 check_maint.py && python3 poly_identity.py
&& python3 verify_room.py` on sympy 1.14.0, numpy 2.5.3, mpmath 1.3.0.

**The rooms' own scripts** are published verbatim under
[`run/rooms/`](../m8_13/run/rooms/) (`stage1/`, `stage2a/`, `stage2b/`, `audit/`),
so every script a return cites by name resolves. Their `.npy` arrays are not
tracked; regenerate them first: `stage1/s03_numeric.py` before `s06` and `s10`,
and `audit/s1_scan.py` before `s5`, `s7` and `s10`. Those scans are seeded
searches that the returns themselves label corroboration only, so a rerun
corroborates and establishes nothing. The room's copy of the author's checker is
not duplicated; it is the pinned file under `scripts/m8_12_author/`.

⚠️ **Two arms in `check_maint.py` are deliberately dead and labelled so.**
Flipping `⟨3 3; 3 3 \| 6 6⟩` is invisible in `‖ρ₆‖²` because `Q = 6` carries a
single term; dropping `Θ`'s `(−1)^{3−m}` phase is invisible **at the hexagon**,
where that phase is `+1` at both `m = ±3`. Each is kept, shown dead on its
arena, and replaced by a live arm. This is D12: an arm has to fire where it is
used.

### What was not verified

The pre-registration's own instruments stay withheld by design until the
package lands, so the gate's 97 checks and the mutation suite's 66 arms remain
the author's report. `m813_equality.py` was not run; its stated conclusion was
reproduced by three other routes. Kawaguchi and Ueda and Romero et al. were not
consulted, so the note's attribution of the bound is recorded as the author's
and unchecked. No stability claim is made, no `MODELS.md` cell moves, and
M8.7's gate is unchanged.

⚠️ **U2's 3c grade is not independent of the 2a grade.** The same room read a
near-identical text twice, so its 2a reading anchors its 2b reading of that
part. N2 itself stays clean, because the 2a return was committed and hashed
before 2b's files existed. The two are not two independent reads and are not
presented as such.

### The author's checker, read as text

27 `gate()` call sites produce the 30 executed checks the note claims. Nothing
gates `TrN̄ = 12`, the variational identity, the equivariance, the max-swap,
**any block's equality locus**, step 3d, step 3e, or the `⊆` direction of step
4. Two gates are weaker than their labels (lines 86 and 87 test constant
identities rather than the maxima they are labelled with), line 61 is inert
because `e1 - e2` is not an atom (line 62 rebuilds the object correctly), and
the strictness sweep excludes exactly **192** of its 20001 points, the ones
nearest the three equality directions. The log's 55 `PASS` lines span ten
headers of which only the 30 under `(L1)`-`(L4)` come from this script.

The note's line 3, "Every step is checked exactly, or armed numerically, by
`check_s3_maximum.py`", overstates that coverage. Its count, "30 checks, 0
failures", is accurate, and the checker reproduces its pinned log byte for byte
on an independent toolchain.

### Record

| Artifact | SHA-256 |
| --- | --- |
| [`run/stage1_return.md`](../m8_13/run/stage1_return.md) | `e9029c840f4e26264e9ef157a00a69558e474a28c891d3a538e4b9fcf991dec8` |
| [`run/stage2a_return.md`](../m8_13/run/stage2a_return.md) | `74fcfdf71ba8f1ad68b686e14147fc3aec6c11adf22eb262682da5e41a2c5ff1` |
| [`run/stage2b_return.md`](../m8_13/run/stage2b_return.md) | `ea6b1f00aaae8c50ecfd7f28dd888e9ddeae5565cb73de1f5e4916c716d83341` |
| [`run/maintainer_grading.md`](../m8_13/run/maintainer_grading.md) | `aa53834e7e2ddf8e8f2be5d1d9c83625da0f7bc934719da63189b8c0d45d4a59` |
| [`run/room_brief_stage1.md`](../m8_13/run/room_brief_stage1.md) | `04adade46974ab68e7dedea641c847fb6bf7e95914752d26a50c44810d7f1d31` |

> The four room returns are **verbatim records** and are published unedited, so
> their bytes match the hashes frozen at the moment each was committed. The
> repository's no-em-dash convention is not applied to them: editing a hashed
> record to match a style rule would break the only thing that makes it a
> record. The surrounding maintainer prose follows the convention.

## TASK REVIEW (2026-09-22)

Task Duration: 02:34 (from 10:11 to 12:45)
Usage Cap Triggered: NO

| Result | |
|---|---|
| ✅ | P1, P2, U1, U2, U3, N1 and N2 pass and D1 is recorded: G2's uniqueness clause is now an audited argument, with the supplied parts named |
| ✅ | The auditor's stage-1 route shares no machinery with the author's and reaches the same maximum, minimum and both extremal sets |
| ✅ | Ten graded parts, no GAP and no DEFECT; nine of ten match the maintainer's grading frozen before stage 2b existed, and on the tenth the maintainer concedes |
| ✅ | Four independent derivations of the equality set agree on three points, and one shows `M₃` is tangent to the bound |
| ✅ | The adversarial audit refuted nothing and closed both completeness claims by certificate, not by search |
| ⚠️ | The author's checker has no exact gate behind the equality set, and its strictness sweep excludes the 192 points nearest the three equality directions |
| ⚠️ | U2's 3c grade is entangled with the 2a grade, since the same room read a near-identical text twice |
| ⚠️ | The auditor's stage-1 argument asserted one load-bearing sub-step; it found and proved it at 2b, and the maintainer's polynomial identity had closed the same ground before 2b ran |

Issues: none blocking.

Deviations from plan: the #585 amendment was written 31 minutes after the go and merged before any
stage-2 file was handed over, so stage 2 ran in its two parts. The maintainer's own checking carried
several defects, none in any room's output, each fixed before use; the two arms that cannot fire on
their arena are kept in `check_maint.py`, labelled.

Action needed: the author package, `m813_equality.py` with its log and the pre-registration's
instruments, is now due; check it against the two pins from the #583 review and the three that #585
superseded.

### Findings

The hexagon orbit is the only maximizer of `r̂₆` on the unit sphere of `V₃`, as an audited argument:
an offline auditor reached it by its own route before seeing the author's, and the author's route
holds at every step once the parts it leaves out are supplied, most substantially step 3d and the
branch analysis of 3c. Saturating the nematic bound is necessary for both extrema and sufficient for
neither, since the whole plane `span{v₃, v₋₃}` saturates it while `r̂₆` sweeps its full range there.

### Research docs created/updated

[Task doc](m8_13_task_details.md), [run record](../m8_13/run/), [rooms' scripts](../m8_13/run/rooms/),
[maintainer instruments](../m8_13/run/scripts/), [roadmap](../m8_roadmap.md),
[briefing](../../__M8_model_briefing.md), [canonical](../m8_theory_canonical.md),
[packet as run](../m8_13/).
