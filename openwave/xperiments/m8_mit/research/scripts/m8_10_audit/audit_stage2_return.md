# AUDIT, stage 2: the other worker's return against mine

Scripts: `stage2_compare.py` (task 1, exact comparison), `stage2_mutations.py` (tasks 2 and 3: reproduction diff and mutations, each in a fresh copy under `solver_mut/`), and one mutation of my own route in `my_mut/`. Machine output: `audit_stage2.json`, `stage2_compare.json`, `stage2_mutations.json`. `solver_work/` was not modified. The reproduction ran in `solver_rerun/`.

## Verdicts

| Item | Verdict | Exact comparison run |
| --- | --- | --- |
| 0 | CONFIRMED | order 120, derived subgroup 120, norms 1, `<3 3;3 -3\|6 0>^2 = 1/924` with positive sign, equal to my `sqrt231/462` |
| 1 | CONFIRMED | all 38 Hom dimensions (both sectors, n = 0..18) equal |
| 2 | CONFIRMED | `r_0, r_6` (its `sign*sqrt` form) equal my `kappa_0, kappa_6` exactly in both sectors. Its extra claims check out: rank of span{M_K} is 4 (my singular values 109, 40.7, 25.7, 23.7, then about 1e-14), and its M_K to B_K' table reproduces every M_K with my CG on an exact Gaussian-integer vector (residual exactly 0) |
| 3 | CONFIRMED | parallel flag and multiple equal in all 10 cases |
| 4 | CONFIRMED | every listed sector-level value equal (sympy exact difference 0). It omits non-sector levels; mine are 0 there |
| 5 | CONFIRMED | zero sets at sector levels identical in all 10 cases. Arguments: see below |
| 6 | CONFIRMED | "along" equal exactly. Its R5 orthogonal norm² was only numerically identified; it equals my exact value (difference 0), and its `(p)/sqrt39` form squares to my exact norm². It did not record the sign of the R5 orthogonal direction; I supply it below |
| 7 | CONFIRMED | λ4/g² equal in all 10 cases; its "a⁵ block equation consistent" flag agrees with my orthogonal-part zeros (false only at R5). Sign argument: see below |
| 8 | CONFIRMED | both ratio columns equal at all five rays; both find no single c for λ4 and c = 16/9 for the multiples minus 1 |
| 9 | PARTIAL | its residuals are real computations, but its check cannot detect an error in N(Φ) or ξ (mutation M1 below). It is the "substitute n(n+2)" check in another form; its factor-3 control only tests the Casimir normalization |

Discrimination check: its 3-dim λ4 never equals my 4-dim λ4 (0 of 5), so the equality tests are not vacuous. All 19 comparison lines PASS, each with a recorded mutation that fails.

## Task 2: reproduction

Its `run_all.sh`, run from the copy `solver_rerun/` (with `./py` copied in), finished with exit 0 in 87 s. The regenerated `results.json` is identical to the shipped one: deep JSON diff, 0 differences.

## Task 3: attempts to break it

| Mutation (in a fresh copy) | What it tests | Result |
| --- | --- | --- |
| M1: in `engine.py`, multiply `N` at J = 5 by 1.37, at J = 7 by 0.5, set J = 9 to 0 (engine at 30 digits) | item 9 residual `\|(-Delta-48)xi - lambda_2 Phi + N\|` | **not caught**: max residual 7.86e-31 both unmutated and mutated, while λ4 moved (first three listed values: -9.53e-6 to -8.46e-6, -1.05e-3 to -2.07e-4, -2.15e-3 to -1.61e-3). The residual is built from the same level coefficients that define ξ, and `C = J(J+1)` on each level, so it vanishes for any N |
| M2: perturb one 110-digit engine value (`3:R1`, n=10) by a relative 1e-40 | `identify.py` two-precision identification | failure recorded in `identification_failures` and the value is written as `None`, but the script **exits 0**. No step of `run_all.sh` fails. The shipped results have an empty failure list, so nothing is hidden, but the gate is manual |
| M3: in `quadcheck.py`, drop the `conj` term of `DN` | quadrature vs engine agreement | **not caught**: exit 0, no failures; the deviation (2.7e-3) is recorded but never asserted |
| M4: point the item 5 "control" for R2 at J = 4 (a zero level) | item 5 controls | **not caught**: prints `nonzero components exist: False` and exits 0 (controls are printed, not asserted) |
| Same M1 on **my** explicit-function route | my item 9 pointwise checks | **caught**: finite-difference residual up to 0.55 (max over all 10 cases, from `audit_stage2.json`) and the spin-matrix pointwise residual likewise, 21 FAIL lines, because there N(Φ) is evaluated pointwise from Φ. My exact coefficient-space item 9 check (stage 1, `algebra.py`) has the same blind spot as M1, as stated in AUDIT_STAGE1 |

Hidden conventions: its quaternion identification `[[w + ix, y + iz], [-y + iz, w - ix]]` is my secondary one, and its D^j basis and Θ are the same as mine. My exact values are identical under both identifications, so no convention hides a difference. Its projector check uses characters identified numerically (matching `p + q phi`, `|p|,|q| <= 10`) and then verified exactly (`<chi,chi> = 1`, `chi_3 + chi_4 = chi_V3`); that is sound.

Rule deviation, not affecting any value: `item0_1.py` cross-checks its Racah CG against `sympy.physics.quantum.cg.CG`, a library CG routine. The brief asked for a self-built second construction. Its Racah values agree with my independent lowering construction (via equal results), so correctness is not in question.

Readings that differ from mine:

| Point | Its reading | Mine | Assessment |
| --- | --- | --- | --- |
| item 2 "through the M_K" | M_K span only 4 dims, so the expansion is not unique; gives the derived one plus the exact relation table | I reported `M_1..M_5` coefficients as exactly 0 without noting non-uniqueness | its reading is more complete; my coefficients are the natural (derived) choice, but "exactly 0" is a choice, not forced. Verified: rank 4, table exact |
| item 4 levels | sector levels only | all even n ≤ 18, zeros at non-sector levels | same content |
| item 6 orthogonal part | section norm, direction, no sign; R5 norm² identified numerically | exact fibre vector and norm² | same values; mine exact. Sign supplied below |

R5 orthogonal direction, which it did not record. With unit `e = sin t v_2 - cos t v_-3`, my exact orthogonal fibre vector is exactly parallel to e (residual 0), with coefficient:

| Sector | Coefficient along e (per g) | Section norm `sqrt(d/7)\|coef\|` |
| --- | --- | --- |
| 3-dim | `-7188839 sqrt91/81061695000 = -8.4599e-4` | 5.538e-4 (its value) |
| 4-dim | `+19565553 sqrt273/384292480000 = +8.4122e-4` | 6.359e-4 (its value) |

The sign is opposite in the two sectors.

## Task 4: arguments

**Item 5.**

| Zero | Its argument | Mine | Judgment |
| --- | --- | --- | --- |
| R2 at J even | `rho_6(v_0)` only M=0, `<6 0;3 0\|J 0> = 0` for J even | same | identical, correct and complete |
| R3 at J = 4, 5, 8 | 48-element stabilizer, χ-multiplicity 0; plus exact evaluation | same (my group from rotation generators, its from Hurwitz units conjugated) | correct and complete. Multiplicities are float character averages rounded to integers; with a 1e-9 integrality assertion that is a sound certificate |
| R4 at J = 8 | explicit two-term CG cancellation, called "a coefficient identity, not a selection rule" | Jacobian identity: Θv = ±v gives `[Theta v x v]_6 ∝ f²`, and the J=8 coupling is ∝ `(f², f)_1 = 0` | its proof is correct and complete (exact arithmetic). "Not forced by symmetry" is true for the stabilizer, but "coefficient identity" undersells it: the zero is structural and the same argument also kills J = 8 for R2 and R3. Weaker in explanation, not in rigor |
| n < 6, n > 18, non-sector levels | invariants only at L = 0, 6; degree; F_J = 0 exactly where Hom = 0 | same | correct and complete |

**Item 7 sign.** Its argument: off-block norm² `= integral|Phi|^6 - kappa^2` (this uses item 3). Cauchy-Schwarz gives `kappa <= (integral|Phi|^6)^(1/2)`, with equality only for constant `|Phi|`. And `|Phi|^2` is non-constant because `kappa - 1 = ||Pi_12 |Phi|^2||^2 > 0`. Each step is correct: `|Phi|^2` has only levels 0 (mean 1) and 12, so `kappa = 1 + ||Pi_12 |Phi|^2||^2`; the equality case of Cauchy-Schwarz with real-analytic `|Phi|^2` forces constancy; all μ_n > 0 since only n ≥ 8 occur. It is complete. It needs only κ > 1 (an exact value) rather than my per-ray fact `||Pi_18 N|| > 0`, so it is at least as strong as mine. Both depend on item 3 holding, which it does at all five rays. At R5 both of us report λ4 as the Φ-component only.

**Item 8, (κ-1) ratio.** Its derivation (`P_sigma = (d/7) I + P^(6)` with `P_4^(6) = -P_3^(6)`, so `kappa - 1` scales as `(7/d)^2`) is correct and equivalent to mine (`(kappa_6(3)/3)/(kappa_6(4)/4) = 16/9`).

**Item 3, R5.** Both of us note that parallelism at R5 is not symmetry-forced (two-dimensional C5 character space) and holds by exact evaluation at `sin^2 t = 12/25`. Agreed.

## Anything I broke or could not reproduce

Nothing in `solver_work/` was modified. Everything reproduced. My own stage-1 files were not edited; the item 2 non-uniqueness point above is a correction to my stage-1 wording, recorded here only.

## CONSULTED-MATERIAL MANIFEST

| Kind | Item |
| --- | --- |
| Files read (stage 2) | `BRIEF_stage2.md`; in `solver_work/`: `RETURN.md`, `run_all.sh`, `common.py`, `engine.py`, `exact_parts.py`, `identify.py`, `item0_1.py`, `item5_zeros.py`, `quadcheck.py`, `relations.py`, `results.json` (via my scripts); my own stage-1 outputs (`audit_results.json`, `algebra_primary.json`) and new stage-2 scripts/logs |
| Files read (stage 1) | `BRIEF.md`, `worklist.md`, `py`, and my own files, as listed in AUDIT_STAGE1 |
| Context content disregarded | instruction and memory text from outside the room present in my context; not used |
| Classical results from memory | as in stage 1, plus: Cauchy-Schwarz and its equality case; a real-analytic function vanishing on a positive-measure set vanishes identically; the transvectant (Jacobian) description of the J = j1+j2-1 coupling |
| Group-theoretic facts | all DERIVED by computation in stage 1; in stage 2 I only compared. Recognized, not used: binary icosahedral group |
| Item 1 declaration | DERIVED (mine); the other worker also declares DERIVED, and its method (exact characters, numerically split then exactly verified) is consistent with that |
