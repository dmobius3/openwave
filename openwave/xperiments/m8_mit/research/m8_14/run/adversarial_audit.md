> Maintainer note: the adversarial audit report, landed verbatim below this line. The scripts and logs it names are kept with the run's local checkpoints, not in the repository; every number marked **[own]** came from them.

# M8.14 adversarial audit of the adjudication

Date 2026-10-01. R = 1 throughout. All scripts and logs are in this folder (`m814_adv/`). Nothing outside it was written, and git was not run. My FE code (`my_fe.py`) imports nothing from the rooms, the auditor or the maintainer. `seam_map_check.py` loads *copies* of the rooms' `t4_fem.py` (in `room_copies/`), but only to read their seam dof maps. It never uses their numbers.

**[own]** marks numbers that came from my own code.

## Verdict table

| # | Claim | Verdict | Evidence (script) |
| --- | --- | --- | --- |
| 1 | VALIDITY: both returns meet every explicit T4 requirement | NOT REFUTED | I read both `t4_fem.py` files line by line. Room A, `dof()` line 79: node `(pi, w_j)` maps to dof `(0, w_{Mw-j}) = (0, -w_j)` with sign -1. Room B, `dof_map` lines 69-71: `(pi, -w_j)` maps to `-(0, w_j)`. Both seams are therefore `phi(pi,-w) = -phi(0,w)`, with the direction and the sign right. I confirmed this dynamically on each room's own map with an asymmetric test function, max error 2.2e-16 at levels 0 and 1. The planted norev, flip and periodic maps all fire (`seam_map_check.py`). Weight, stiffness, Dirichlet rows, grading `(pi/2)(k/16)^2`, 8 cells across, 4 midpoint refinements, p / Richardson / err and all 7 widths are all as specified. **[own]** My independent assembler reproduces every room level value (7 widths x 5 levels x 6 eigenvalues) to max rel 1.8e-13 for A and 1.8e-13 for B (`run_fe.py`, `run_fe.log`). |
| 2 | F1a, the reduction, established in both rooms | NOT REFUTED, with two SUPPLIED parts per room (below) | Both rooms establish a unitary, form-preserving map `U` from the band onto the lune `[-pi/2, pi/2] x [-W, W]`, metric `dpsi^2 + cos^2 psi dtheta^2`. It is a whole-space equivalence, not only sector by sector. **[own]** Symbolic: `-X(y,-w) = X(y-pi,-w)`; the seam condition is equivalent to continuity of `u` at `psi = 0`; the metric is preserved (`sym_checks.py`, S1). The planted "no sign flip" twin fires. |
| 3 | F1b, the spectrum, complete in both rooms | NOT REFUTED | Both completeness arguments (weighted isometry to `L^2((1-x^2)^nu dx)`, Gegenbauer degree exactly n, Weierstrass density, then tensor with a complete `{e_k}`) are correct. **[own]** `(nu+n)(nu+n+1)` holds symbolically in nu for n = 0..5 (the planted `(nu+n)(nu+n+2)` fires). Orthogonality weight `(1-x^2)^nu` and degree n checked at three nu values (`sym_checks.py`, S2). |
| 4 | F2, N1, N2 | NOT REFUTED | Recomputed from each `results.json` without the maintainer's scripts (`score.py`, `score.log`). Both rooms are Reproduced at all 5 widths: rel miss at most 2.98e-7 (W=1/4), rel err at most 1.22e-4 (W=1/4), p in [1.9967, 2.0004], bottoms monotone from above. Reported and recomputed triples agree to 2.6e-15. N1: \|extrap - 2\| = 5.03e-8, tolerance 2e-3. N2: J bottom -0.222230092 against err 1.38e-4, so \|J\|/err = 1609 (> 10). **[own]** My FE gives the same outcome and the same numbers (extrap equal to 3e-12). The maintainer's scripts could not have hidden a failure on these returns. A structural caveat is in D4. |
| 5 | A1, A2, A3 established in both rooms | NOT REFUTED | **[own]** A1: a finite-difference second derivative of the true area in R^4 for an admissible phi at W = 1.4 gives A''(0) = 0.31155050 against q_J = 0.31155162, rel 3.6e-6. The planted twin without `-2 phi^2` fires (`sym_checks.py`, S4). A2: `r cot r = 1 + O(r^2)` and `r^2/sin^2 r = 1 + O(r^2)`, so the indicial roots are +-nu, and `r^{-nu}` is in `L^2(r dr)` iff nu < 1, with nu = 1 log-divergent (S3). nu_k = k alpha0 >= 1 for W <= pi/2. A3 follows: `(alpha0 - 1)(alpha0 + 2) > 0` iff W < pi/2. |
| 6 | Seam claim | NOT REFUTED (one precision added) | **[own]** `run_fe.py` runs all 7 widths x 5 levels on the frozen protocol, lowest six, relative to the spec seam.<br>Flip sign: identical to 4.7e-15 at every level.<br>Cut, D at y=pi, y=0 free: identical to 8.2e-14 (cut D at y=0 / N at y=pi: 7.5e-14).<br>No reversal: differs by 3.1e-2 at L0 and 4.3e-5 at L4. The bottom differs by at most 4.9e-9 at L4. Every gap shrinks under refinement, and the extrap differs by at most 3.2e-6 (W=1/4) (`variant_protocol.log`). Running no-reversal on a mesh whose lower sheet uses the mirrored diagonal reproduces the spec spectrum to 4.6e-14, which shows the residual is purely mesh asymmetry.<br>Both ends free: the second level equals the bottom to 8e-14 at every width.<br>Both ends D: the bottom disappears, and the levels pair up as the odd spectrum doubled. |
| 7 | CONTAINMENT | NOT REFUTED | jq over all 4 transcripts (`results_*.txt` = extracted tool_results). Every Read/Write/Edit path is inside its own room. Every Bash command begins with `./py` and uses relative paths. The 2 non-`./py` compounds (A: a `for` loop; B: `PLANT=1 ./py ...`) were harness-denied and never ran. Tools were Bash/Edit/Read/Write only and `mcp_servers` was empty. The `.room.sb` profiles are identical modulo room name: deny read/write `/Users` and `/private/tmp` except the room, `deny network*`. No withheld-source marker appears in any tool_result or assistant text (grep for author/page/repo tokens: 0 hits). Spec sheet, auditor brief and BRIEF hashes are identical between the repo packet, the stage briefs and the rooms. |
| 8 | Stage-2 SUPPLIED grade, Room A T3 step 3 | NOT REFUTED | The supplied argument is correct: on the core, Parseval in theta gives `a(u) = sum a_k(u_k)`. Truncation in k plus 1D approximation in each sector shows that the form closure is the direct sum of the sector closures, so the operator is the direct sum of the sector operators. Precision: for W <= pi/2 the supply is not even needed, because Room A's own T2 step 1 (`closure T = direct sum of closure T_k`) plus essential self-adjointness already gives A = direct sum of A_k. It is needed only at W = 17/10 (Friedrichs). |

## F1a graded step by step (claim 2)

### Room A, RETURN.md section 1, "the unfolding"

| Step | Grade |
| --- | --- |
| Hypotheses G1-G3; the map `(psi, theta, u) = (y, w, phi)` on the upper sheet and `(y - pi, -w, -phi)` on the lower | ESTABLISHED |
| 1 The metric becomes `dpsi^2 + cos^2 psi dtheta^2` | ESTABLISHED (S1 metric) |
| 2 The seam becomes the interior line `psi = 0`, and the seam condition is continuity of u (under R1, smoothness of u) | ESTABLISHED (S1 seam) |
| 3 The fiber becomes the two poles, and `w = +-W` becomes `theta = +-W` | ESTABLISHED |
| 4 Integrands preserved, so U is unitary on L^2 and preserves the form | ESTABLISHED |
| Conclusion: "every problem below is a problem on the lune L, Dirichlet on the sides, poles excluded; an honest lune for W <= pi/2" | ESTABLISHED, SUPPLIED (S-a, S-b) |

### Room B, RETURN.md section 0, "Reduction used throughout"

| Step | Grade |
| --- | --- |
| Unfolding `-X(y,-w) = X(y',w')`; M(W) is the lune with its two vertices identified to p; the seam becomes the equator arc | ESTABLISHED (S1) |
| Functions `psi = U phi` via `e4 -> -e4`; the seam condition is continuity across the equator | ESTABLISHED |
| U unitary onto `L^2(L_W)` and preserves the Dirichlet integral; edges become sides; "vanish near p" becomes "vanish near both vertices" | ESTABLISHED |
| T2 reading D0: "the form closure is unchanged, because kinked functions lie in the H^1-closure of D0" | ESTABLISHED, SUPPLIED (S-a). Stage 2 graded this plain ESTABLISHED, but the sentence is asserted, not argued |
| Identification with the lune's Dirichlet problem | ESTABLISHED, SUPPLIED (S-b) |

Supplied parts, which are the same for both rooms:

| ID | Supplied part |
| --- | --- |
| S-a | The admissible class (value-only seam, H^1 completion) maps onto the H^1-closure of the smooth vertex-excluded core on the lune. U carries a value-seam function to a Lipschitz function, zero on the sides and near the poles. Its sector coefficients are Lipschitz with compact support in (-pi/2, pi/2), so k-truncation plus 1D mollification approximates it in the form norm by core functions. Room A states it in R1 ("the completion only sees the trace"), Room B in its T2 reading, both without the argument. |
| S-b | The vertex-excluded lune problem is the lune's Dirichlet problem (`H^1_0(L)`). This holds for every W because points have zero H^1-capacity in 2D. For W <= pi/2 it also follows from each room's own established T2: the minimal operator is essentially self-adjoint, and the standard Dirichlet Laplacian of the lune extends it. |

Both rooms establish the equivalence as a whole space, not sector by sector, onto the lune of opening 2W/R. Under the task's declared reading of "agrees" ("ESTABLISHED, as written or with a supplied part named"), F1a agrees in both rooms. The maintainer should record S-a and S-b as named supplied parts.

## Defects demonstrated

| ID | Where | Defect | Demonstration | Effect on verdict |
| --- | --- | --- | --- | --- |
| D1 | Stage-2 auditor, `s2_audit.py` K2 | The seam check rebuilds the bottom eigenvector through the room's own dof map and tests `phi(pi,-w) + phi(0,w) = 0` plus correlation with the closed-form bottom. The bottom mode is even in w, so the test cannot see a missing w-reversal. Its planted twin (the rooms' sign plant) only falsifies the sign. | `k2_blindness.py` **[own]**: the no-reversal map passes the K2-style test (residual 3.6e-15, corr 1.000000). Flip and periodic fail. A w-odd mode (k=2) would catch it (residual 2.0). | None. Seam direction is correct by code reading and by `seam_map_check.py`. K2 should not be cited as evidence of seam direction. |
| D2 | Maintainer `seam_variants.py` | The docstring says "level 3" but the code calls `assemble(W, 2, v)`, so it runs level 2 only, at W in {1, 1.4} only. The run state (G1f(1)) already records that level 2 was mistaken for a decision on the no-reversal glue. | Read the code. | The seam claim holds on the full protocol (my run, claim 6). Record it from a run at all 5 levels and 7 widths, not from level 2. |
| D3 | Maintainer run state, G1f(2) | "Solver B RETURN T5 table omits W=1 (json has it)" is false. Room B RETURN.md line 266 is the T5 row `\| 1 \| (π − 2)(π + 4)/4 \| 2.038197427067 \| + \| 0 \| 0 \|`, and the T3 table has W=1 at line 122. This is probably a parser hit on the level-1 rows `\| 1 \| ...`. | `grep -n "^\| 1 \|" RETURN.md` | Do not record it as an artifact. |
| D4 | Maintainer `compare.py` / `frozen_check.py` | F2 and the controls are scored on the room's REPORTED `p, extrap, err`. The recompute flags (`p_consistent`, `err_consistent`, `reading`) are informational and do not gate F2, so a return whose reported triple disagreed with its own levels would still score on the reported numbers unless someone reads the flags. | Read the code. | None here: reported equals recomputed to 2.6e-15 (`score.py`), and all flags are true in `cmp_a.json` / `cmp_b.json`. |
| D5 | My own instrument (recorded, not tuned) | (a) My first seam check in `run_fe.py` had the same blindness as D1, and its planted norev twin did not fire (`run_fe.log`, the FAIL `[planted]` line). I replaced it with `seam_map_check.py`, an asymmetric test function on the map, whose twins fire. (b) A dense-LAPACK-vs-eigsh check failed its own 1e-10 threshold at W = 1/4, 1/2, 1, level 1 (1.5e-10 to 3.0e-10). The eigsh values match both rooms to 1.8e-13, so the dense path is the less accurate one. Both are reported as they ran. | `run_fe.log` | None. |

## Why the seam variants do what they do (claim 6)

| Variant | Mechanism | Result **[own]** |
| --- | --- | --- |
| Flip sign | Zero data on the fiber means no element couples free dofs of the two sheets. Multiplying one sheet by -1 is an exact unitary, both continuous and discrete. | Exact (4.7e-15) |
| No reversal | `w -> -w` on one sheet is an isometry of that half (metric depends on y only, Dirichlet sides symmetric), so the continuum spectrum is identical. Discretely it mirrors that sheet's triangle diagonals. | O(h^4) gap, and O(h^2) on the split degenerate clusters at W = pi/2. Exact (4.6e-14) once the lower sheet uses the mirrored diagonal |
| Cut, D at y=pi, y=0 free | The lune is symmetric under `psi -> -psi`. Its modes split into even (n even, Neumann at the equator) and odd (n odd, Dirichlet). The upper sheet with y=0 free carries the even part, and the lower sheet with y=pi Dirichlet carries the odd part. The spec-seam mesh, unfolded, is mirror-symmetric about the equator, so this holds discretely too. | Exact (8.2e-14) |
| Both free | Two copies of the even part | Bottom repeated as level 2 |
| Both D | Two copies of the odd part | Bottom (1,0) absent |

So a seam error of the flip, no-reversal or D/N-cut kind is invisible to every numerical gate in the protocol, including the maintainer's F1b seam diagnostic. Only code reading or a dof-map test with a w-odd function catches it. Both rooms pass that test.

## Reading sensitivity (not a defect)

"Error estimate within the tolerance" is scored as `err <= 1e-3 * target`, which matches the frozen text ("to a relative 10^-3"). Under an absolute 1e-3 reading, W = 1/4 (err 5.6e-3) and W = 1/2 (err 1.18e-3) would be Unresolved. The frozen text does not support that reading, but the adjudication should state the relative reading explicitly.

## Suspicions (not demonstrated)

None material. Both rooms take Weyl's alternative, deficiency theory, Kato/Friedrichs and min-max from memory, as their manifests declare. The grades are audited arguments, not verified theorems.

## Not checked

| Item | Why |
| --- | --- |
| Whether a copy of a withheld source exists at a path the sandbox allowed (anything outside `/Users`, `/private/tmp`, and the temp dirs) | The profile is `allow default` with targeted denies. No room attempted any outside read, so no path was exercised |
| Content of thinking blocks, if any were redacted in the stream | I grepped only the plain-text fields |
| The maintainer's `fe608.py` and its seam-variant numbers | Not imported or rerun, by brief. My own code covers the claim |
| Room scripts' T1/T3 symbolic files | Not rerun. Their content is covered by my S1-S4 checks and by stage 2's G-checks |

## Script index (this folder)

| Script | Output | Purpose |
| --- | --- | --- |
| `my_fe.py` | (module) | Own P1 assembler: mpmath tanh-sinh y-moments, exact sympy w-sections, both diagonals, 8 seam variants, eigsh at sigma = -1 |
| `run_fe.py` | `fe_out.json`, `run_fe.log` | Spec protocol at 7 widths x 5 levels compared with both rooms; all seam variants; mirror test |
| `score.py` | `score_out.json`, `score.log` | F2 / N1 / N2 from the rooms' `results.json` and from my FE, with no maintainer import |
| `seam_map_check.py` | stdout | Rooms' seam maps tested dynamically with an asymmetric function; planted twins |
| `k2_blindness.py` | `k2_blindness.log` | D1 demonstration |
| `sym_checks.py` | `sym_checks.log` | S1 unfolding, S2 Gegenbauer spectrum, S3 indicial and limit-point, S4 A1 by finite differences in R^4 |
| `frozen_table_check.py` | `frozen_table_check.log` | Task-doc FROZEN VALUES recomputed at 40 digits (all PASS) |
| `variant_protocol.log` | | Frozen protocol (p, extrap) applied to the flip / no-reversal / cut variants |
