# M8.9 S1b Adjudication of exposed target reading (attempt q3a)

## 1. Room manifest verification

`ADJ_ROOM_MANIFEST.json` SHA-256:
`68b92e1d3fcc22096bb0095c2b55556d31a21df50a65eab3532e04d8be1c470a`
**MATCHES** the out-of-band value given in the commission prompt.

## 2. Frozen contract hash verification

All six frozen-region hashes verified using the boundary-marker extraction specified in each
document's own freeze record:

| document | expected | computed | status |
| --- | --- | --- | --- |
| `S1B_DECISION_RULE.md` | `c44c603a...7d297` | `c44c603a...7d297` | MATCH |
| `S1B_ADDENDUM_1.md` | `6da36a1c...46746` | `6da36a1c...46746` | MATCH |
| `S1B_ADDENDUM_2.md` | `14011c33...1a222` | `14011c33...1a222` | MATCH |
| `S1B_ADDENDUM_3.md` | `e3304fe9...6e28c` | `e3304fe9...6e28c` | MATCH |
| `S1B_ADDENDUM_4.md` | `98484b0e...00634` | `98484b0e...00634` | MATCH |
| `S1B_ADDENDUM_5.md` | `7ca059f0...2b3e5` | `7ca059f0...2b3e5` | MATCH |

Full hashes:
- `S1B_DECISION_RULE.md`: `c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297`
- `S1B_ADDENDUM_1.md`: `6da36a1c672772e8c731d4d66f20a77e81779d19d5cb883e617c028d43f46746`
- `S1B_ADDENDUM_2.md`: `14011c338ce331d9bb3424f0ed5aaff31a1c89a9ceafbc2bae41a91ebb41a222`
- `S1B_ADDENDUM_3.md`: `e3304fe993f7a2ea523cfcbc4110c0f37174caea95383b17e465c8fa85b6e28c`
- `S1B_ADDENDUM_4.md`: `98484b0edd2de97e34093a564674cd8f128a7bca38f96960c1007f7f45f00634`
- `S1B_ADDENDUM_5.md`: `7ca059f074059f925a7f231fcb1ac93932e890121ef9abd1975a4df78542b3e5`

## 3. Evidence manifest verification

### 3.1 OUTPUT_MANIFEST.json root of trust

Per addendum 5 section A5.2, the output manifest hash must be verified before any record it pins
is used:

    evidence/OUTPUT_MANIFEST.json
    expected: f5401c5179d0cde42a8763de175542d3d04aaa316be0329ec30a42cd9d6bc3a4
    computed: f5401c5179d0cde42a8763de175542d3d04aaa316be0329ec30a42cd9d6bc3a4
    status:   MATCH

### 3.2 Target records (addendum 4 section A4.3)

    evidence/target_n12.json
    expected: d96c7fdb3ab192b48f04c49040dbb8f6795379a35fb04330078bff22115a5b91
    computed: d96c7fdb3ab192b48f04c49040dbb8f6795379a35fb04330078bff22115a5b91
    status:   MATCH

    evidence/target_n20.json
    expected: 4e8c0c72fac0a8e2f127fc1629bb4afa4b2d90e9b12b9df97714d56b1057fed7
    computed: 4e8c0c72fac0a8e2f127fc1629bb4afa4b2d90e9b12b9df97714d56b1057fed7
    status:   MATCH

### 3.3 Gate and evidence records verified against authenticated manifest

Every evidence file in the room verified against its entry in the now-authenticated
`OUTPUT_MANIFEST.json`. The manifest uses paths relative to the q3a attempt directory
(e.g. `results/G-ALIGN.json`); the adjudication room maps these to `evidence/G-ALIGN.json`.

| room file | manifest entry | computed hash | status |
| --- | --- | --- | --- |
| `G-REAL.json` | `results/G-REAL.json` | `2fa2c238...0ee10` | MATCH |
| `G-RANK_n12.json` | `results/G-RANK_n12.json` | `827c7233...19786` | MATCH |
| `G-RANK_n20.json` | `results/G-RANK_n20.json` | `ea166655...6e968` | MATCH |
| `G-SUBSPACE_n12.json` | `results/G-SUBSPACE_n12.json` | `85be5bcd...280c6` | MATCH |
| `G-SUBSPACE_n20.json` | `results/G-SUBSPACE_n20.json` | `e5451517...504fd` | MATCH |
| `G-ALIGN.json` | `results/G-ALIGN.json` | `7da3ece9...be9f0` | MATCH |
| `G-SAMPLE_n12_svd.json` | `results/G-SAMPLE_n12_svd.json` | `92f29f08...fd80d` | MATCH |
| `G-SAMPLE_n12_avg.json` | `results/G-SAMPLE_n12_avg.json` | `6ff4cfb7...03438` | MATCH |
| `G-SAMPLE_n20_svd.json` | `results/G-SAMPLE_n20_svd.json` | `f862a267...4d24c` | MATCH |
| `G-SAMPLE_n20_avg.json` | `results/G-SAMPLE_n20_avg.json` | `ac9829b3...e713d` | MATCH |
| `G-BASIS_n12_svd.json` | `results/G-BASIS_n12_svd.json` | `4c8e6fe3...558bc` | MATCH |
| `G-BASIS_n12_avg.json` | `results/G-BASIS_n12_avg.json` | `2cfb9df9...628bd9` | MATCH |
| `G-BASIS_n20_svd.json` | `results/G-BASIS_n20_svd.json` | `a99322de...fbb0d4` | MATCH |
| `G-BASIS_n20_avg.json` | `results/G-BASIS_n20_avg.json` | `a8724d33...766e2f` | MATCH |
| `G-DISCRIM.json` | `results/G-DISCRIM.json` | `881b1794...ce301` | MATCH |
| `G-WIRE.json` | `results/G-WIRE.json` | `7dee538d...cdf1c1e` | MATCH |
| `ladder_controls.json` | `results/ladder_controls.json` | `04aad6be...edde0` | MATCH |
| `adjudicator.json` | `results/adjudicator.json` | `84a101bc...7b0170` | MATCH |
| `target_n12.json` | `results/target_n12.json` | `d96c7fdb...5b91` | MATCH |
| `target_n20.json` | `results/target_n20.json` | `4e8c0c72...fed7` | MATCH |

All 20 evidence files match their manifest entries.

### 3.4 Expected manifest discrepancy: QUALIFICATION_NOTE.md

`QUALIFICATION_NOTE.md` full-file hash `237b64d55d4a31ff480b71495306e3107b17710577a5226127fc6bfd16e1960b`
does NOT match its manifest entry `098631325e9692033597ea3f63b0c5955178c929a9330570573f9767a362175c`.

Truncating immediately before `## Output Manifest` (at byte 4591, removing the preceding blank
line's newline and everything after): hash
`098631325e9692033597ea3f63b0c5955178c929a9330570573f9767a362175c` **MATCHES exactly**.

The append point is precisely located: the original note ended at byte 4591 with
`- All inputs unchanged: YES\n`, and the appended content is the `## Output Manifest` section
(preceded by a blank line). This confirms the discrepancy is entirely from post-manifest
appending, as disclosed in addendum 4 section A4.3 and TASK.md.

**The note diverges NOWHERE else.** The note is non-adjudicative per addendum 5 section A5.3:
no branch rule takes an input from it.

### 3.5 Manifest entries not present in the room

The manifest contains 34 entries. The room holds 21 of these (20 verified above plus
QUALIFICATION_NOTE.md with its explained discrepancy). The following 13 entries are NOT present
in the adjudication room and cannot be checked here:

- `s1b_qual.py`
- `svd_n20.log`
- `run.log` (also known to differ from manifest, non-adjudicative)
- `svd_n20_precompute.py`
- `results/vacuity_audit.json`
- `results/C_svd_n12.npy`
- `results/C_svd_n20.npy`
- `results/bridge_n12.json`
- `results/bridge_n20.json`
- `results/svd_n12.json`
- `results/svd_n20.json`
- `results/contract_contradictions.json`
- `results/infrastructure.json`

None of these absent files is required as input to any branch rule. The gate records, target
records, ladder controls, and adjudicator cases that the branch algorithm needs are all present
and verified.

## 4. Branch algorithm trace

The branch algorithm is defined in the frozen decision rule, first match wins, evaluated in the
order below. Every quantity read is from an evidence file verified against the authenticated
manifest in sections 3.1-3.3 above. `adjudicator.json` contains SYNTHETIC reachability and
precedence cases; it is not used as live evidence and is not consulted below.

### Rule 1: any instrument gate fails -> S1b-DEFECT

All eight gates and the ladder controls are checked below. Every gate record carries a `pass`
field; the adjudication independently confirms each against the frozen criteria.

**G-REAL** (pointwise realization, the load-bearing gate):
- `worst_residual_correct` = 1.790e-15, criterion < 1e-10: **PASS**
- `best_residual_no_transpose` = 9.443e-01, criterion > 1e-10: **PASS**
- Levels swept: [3, 4, 5], all >= 3 as required (n=2 trap avoided)
- Record `pass`: true

**G-RANK at n=12:**
- Returned dimension = 13, expected = 13: **PASS**
- Gap state: "measured", gap = 7.218e+09
- Record `pass`: true

**G-RANK at n=20:**
- Returned dimension = 21, expected = 21: **PASS**
- Gap state: "measured", gap = 2.958e+09
- Record `pass`: true

**G-SUBSPACE at n=12:**
- k_svd = 13, rank_avg = 13: dimensions agree
- sin(theta_C) = 7.289e-11, criterion <= 1e-6: **PASS**
- Rank taken at absolute cutoff 1e-8: rank_abs_1e8 = 13
- Record `pass`: true

**G-SUBSPACE at n=20:**
- k_svd = 21, rank_avg = 21: dimensions agree
- sin(theta_C) = 8.086e-11, criterion <= 1e-6: **PASS**
- Rank taken at absolute cutoff 1e-8: rank_abs_1e8 = 21
- Record `pass`: true

**G-ALIGN:**
- structural_predicate: true
- seed_orbits_sorted: true
- arm_red: true (mutation detected)
- Record `pass`: true. **PASS**

**G-SAMPLE at n=12 (svd route):**
- rank = 13 (expected 13): exact match
- kappa = 37.07, criterion <= 1e6: **PASS**
- Record `pass`: true

**G-SAMPLE at n=12 (avg route):**
- rank = 13 (expected 13): exact match
- kappa = 37.07, criterion <= 1e6: **PASS**
- Record `pass`: true

**G-SAMPLE at n=20 (svd route):**
- rank = 21 (expected 21): exact match
- kappa = 627.8, criterion <= 1e6: **PASS**
- Record `pass`: true

**G-SAMPLE at n=20 (avg route):**
- rank = 21 (expected 21): exact match
- kappa = 627.8, criterion <= 1e6: **PASS**
- Record `pass`: true

**G-BASIS at n=12 (svd route):**
- A_covariance_rel = 9.735e-16, criterion <= 1e-10: **PASS**
- P_diff_rel = 5.093e-16, criterion <= 1e-10: **PASS**
- Record `pass`: true

**G-BASIS at n=12 (avg route):**
- A_covariance_rel = 8.445e-16, criterion <= 1e-10: **PASS**
- P_diff_rel = 7.450e-16, criterion <= 1e-10: **PASS**
- Record `pass`: true

**G-BASIS at n=20 (svd route):**
- A_covariance_rel = 5.704e-16, criterion <= 1e-10: **PASS**
- P_diff_rel = 6.848e-16, criterion <= 1e-10: **PASS**
- Record `pass`: true

**G-BASIS at n=20 (avg route):**
- A_covariance_rel = 5.909e-16, criterion <= 1e-10: **PASS**
- P_diff_rel = 7.107e-16, criterion <= 1e-10: **PASS**
- Record `pass`: true

**G-DISCRIM:**
- Green parent: K = 0.0 (criterion: <= 100 eps ||A||_2), J = 0.0 (criterion: <= 1e-12): **PASS**
- Arm A: |K - 0.5| = 0.0 (criterion: <= 1e-12), J = 0.0 (criterion: <= 1e-12): **PASS**
- Arm B: |K - 2.0| = 0.0 (criterion: <= 1e-12), |J - 2.0| = 1.4e-14 (criterion: <= 1e-12): **PASS**
- Record `pass`: true

**G-WIRE:**
- ||L Q_0||_{M_h,F} = 5.457e-09, criterion <= 1e-8: **PASS**
- arm_red: true (mutation detected, ||L_mut Q_0|| = 1.000e-04, four orders above gate)
- Record `pass`: true

**Ladder controls** (addendum 1 section A1.1):
- k=2 exact-zero control: J = [0, 0, 0] -> COLLAPSES via rule 1. **Correct.**
- k=3 ill-conditioning refusal: J = [6.47e-06, 9.10e-11, 2.19e-17] -> AMBIGUOUS via rule 2. **Correct.**
- Synthetic rule-4 triplet: J = [1e-8, 1e-12, 1e-12] -> COLLAPSES via rule 4. **Correct.**
- Rule-5 fall-through: J = [1e-8, 2e-11, 2e-11] -> AMBIGUOUS via rule 5. **Correct.**
- Edge cases (2.0x and 0.5x tail ratios): both COLLAPSE via rule 4. **Correct.**

**All instrument gates PASS. All ladder controls fire at their required rules.
Rule 1 does NOT fire.**

---

### Rule 2: any target has an AMBIGUOUS J on EITHER construction -> S1b-NO_LABEL

The J ladder verdicts for each target and construction are read from the verified target records.

**n=12, svd route:**
J readings: [2.2010292238211204, 2.2010292238212394, 2.2010292238212394]
Applying the frozen ladder algorithm (first match wins):
1. All zero? No (J(64) = 2.201). Skip.
2. J(50) < 0.5*J(30)? 2.201 < 1.101? No. J(50) > 2*J(30)? 2.201 > 4.402? No. Skip.
3. J(50) >= 0.5*J(64)? 2.201 >= 1.101? Yes. -> **PERSISTS** (rule 3)

**n=12, avg route:**
J readings: [2.201029212772619, 2.2010292127725837, 2.2010292127725837]
Same analysis: J(50)/J(64) ~ 1.000. -> **PERSISTS** (rule 3)

**n=20, svd route:**
J readings: [19.25417414284196, 19.254174142841695, 19.254174142841695]
Same analysis: J(50)/J(64) ~ 1.000. -> **PERSISTS** (rule 3)

**n=20, avg route:**
J readings: [19.254174160953806, 19.25417416095387, 19.25417416095387]
Same analysis: J(50)/J(64) ~ 1.000. -> **PERSISTS** (rule 3)

No AMBIGUOUS verdict on any construction at any target.
**Rule 2 does NOT fire.**

---

### Rule 3: the two constructions disagree on the J reading for any target -> S1b-NO_LABEL

**n=12:** svd says PERSISTS, avg says PERSISTS. **Agree.**
**n=20:** svd says PERSISTS, avg says PERSISTS. **Agree.**

**Rule 3 does NOT fire.**

---

### Rule 4: any target is PERSISTS on BOTH constructions -> S1b-SPECTRAL

**n=12:** svd PERSISTS, avg PERSISTS. Both constructions: **YES.**
**n=20:** svd PERSISTS, avg PERSISTS. Both constructions: **YES.**

Both targets independently satisfy this condition. Rule 4 requires only "any target" to qualify.

**Rule 4 FIRES. Outcome: S1b-SPECTRAL.**

---

### Rules 5-7 (evaluated for completeness, though rule 4 has already fired)

**Rule 5: all J collapsed or zero on both, and constructions disagree on
||K_n||_2 > K_floor(n) for any target -> S1b-NO_LABEL**

Precondition "all J collapsed or zero on both" is FALSE: all four J readings are PERSISTS.
Rule 5 **cannot fire**.

**Rule 6: all J collapsed or zero on both, and some target has ||K_n||_2 > K_floor(n)
on BOTH -> S1b-ADJOINT**

Same precondition is FALSE. Rule 6 **cannot fire**.

**Rule 7: else -> S1b-NULL**

Not reached; rule 4 fires first.

---

### Supporting data from the target records (recorded, not routing)

For transparency, the K and J values from the verified target records:

| quantity | n=12 svd | n=12 avg | n=20 svd | n=20 avg |
| --- | --- | --- | --- | --- |
| ||A_n||_2 | 175.184 | 175.184 | 508.332 | 508.332 |
| ||K_n||_2 | 3.607 | 3.607 | 64.410 | 64.410 |
| K_floor(n) | 6.851e-07 | 6.851e-07 | 9.926e-07 | 9.926e-07 |
| K above floor | true | true | true | true |
| J(64) | 2.201 | 2.201 | 19.254 | 19.254 |
| J(50) | 2.201 | 2.201 | 19.254 | 19.254 |
| J ladder | PERSISTS | PERSISTS | PERSISTS | PERSISTS |
| kappa(V) | 8.390 | 8.390 | 12.899 | 12.899 |
| J_bf | 3.264e-11 | 3.264e-11 | 1.456e-10 | 1.456e-10 |
| K_floor dominant | discrepancy | discrepancy | discrepancy | discrepancy |

The J readings are massive compared to any floor: J = 2.201 at n=12 and J = 19.254 at n=20,
both stable across all three precision rungs and both constructions to twelve significant
figures. The Bauer-Fike bound J_bf is 10 to 11 orders of magnitude below the readings, with
kappa(V) well within the demotability threshold of 1e8.

## 5. Forced outcome

**The frozen branch algorithm forces outcome S1b-SPECTRAL, via rule 4.**

The condition is satisfied at both targets independently: J is PERSISTS on both the SVD and
averaging constructions at n=12 and at n=20. Rule 4 is the first rule in the frozen precedence
order whose condition holds.

## 6. Licensed sentence

From the frozen decision rule, under `S1b-SPECTRAL`:

> the trivial fibre alone can produce a non-real compressed action, so nontrivial fibre
> transport is NOT NECESSARY for that phenomenon, and the base discretization or the scalar
> quotient reduction is strongly implicated.

Nothing beyond this is licensed by S1b-SPECTRAL. Per the decision rule, both S1b-SPECTRAL and
other outcomes leave M8.9 OPEN: a separately frozen interaction or reduction comparison MAY be
commissioned, but S1b-SPECTRAL does not commission one and does not authorize S2 as filed.

## 7. Ambiguity assessment

No ambiguity was encountered in applying the frozen branch algorithm to these particular values.

The readings are so large (J ~ 2.2 and ~19.3, stable to twelve figures across precision rungs
and both constructions) that no edge case, threshold proximity, or rounding question arises.
Every gate passes with wide margins. Every ladder verdict is PERSISTS by rule 3 with ratios
indistinguishable from 1.0. The two constructions agree on every verdict. The frozen rules
force a unique outcome on these bytes without any interpretive choice.

---

*Adjudication issued by the commissioned Adjudication Unit, 2026-08-26.*
*No operator was constructed, no SVD or eigensolve was run, no gate or target quantity was
recomputed. This trace applies the frozen branch algorithm to the pinned bytes and reports
the rule the frozen precedence order forces.*
