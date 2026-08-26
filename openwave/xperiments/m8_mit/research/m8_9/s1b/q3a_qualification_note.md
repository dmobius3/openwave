# S1b Qualification Note — attempt q3a

## Q0: Provenance

**Date:** 2026-08-26
**Attempt ID:** q3a (supplied in handoff prompt, not chosen)

### Environment
- Python 3.13.13
- numpy 2.5.0
- scipy 1.18.0
- mpmath 1.3.0

### Room import gate
`PYTHONPATH=".:m8_5b:m8_5b/pilot:m8_5b/production:m8_5b/gates" python3 room_import_gate.py`
Result: **ROOM READY**
- 10 entry points; 15 required-closure modules; 50 supplied modules
- 35 supplied modules outside the required closure: listed, NOT imported
- 0 required-module failures
- All four frozen contract documents verified OK
- Prior script `prior/s1b_qualification.py` verified OK

### Room manifest
`ROOM_MANIFEST.json` SHA-256: `7d4e0fdb60d430971359cb3f38e7eef44c9259be10d1f377160d76702623da2e`
Matches handoff-prompt value (which appears in no file in the room). Verified independently via `shasum -a 256`.

### Frozen document hashes (independently verified)
| Document | Boundary | Expected SHA-256 | Verified |
|---|---|---|---|
| `contract/S1B_DECISION_RULE.md` | `<!-- FREEZE-BOUNDARY -->` | `c44c603a...7d297` | YES |
| `contract/S1B_ADDENDUM_1.md` | `<!-- ADDENDUM-BOUNDARY -->` | `6da36a1c...6746` | YES |
| `contract/S1B_ADDENDUM_2.md` | `<!-- ADDENDUM2-BOUNDARY -->` | `14011c33...a222` | YES |
| `contract/S1B_ADDENDUM_3.md` | `<!-- ADDENDUM3-BOUNDARY -->` | `e3304fe9...e28c` | YES |

### Prior script (archival only, per A2.3)
`prior/s1b_qualification.py` SHA-256: `5a9e04845375c4d12c3a475607f20a1d5f13cc82829d64caa022c82d5e784802` — matches pinned hash.
NOT executed, NOT relocated, NOT patched. Cannot run in this room by construction (manifest assertion at line 148 references round-1 manifest hash). This is a recorded fact, not an obstacle.

### Import attestation
No OpenWave or project code, and no required-closure module, was imported from outside this room. The standard library (Python 3.13.13) and third-party packages (numpy 2.5.0, scipy 1.18.0, mpmath 1.3.0) are external by necessity and are covered by the environment record above.

**Q0 PASS.**

## G-REAL
- worst_correct: 1.79e-15
- best_no_transpose: 9.44e-01
- **PASS**

## G-DISCRIM
- green: K=0.00e+00, J=0.00e+00
- arm A: K=0.500000, J=0.00e+00
- arm B: K=2.000000, J=2.000000
- **PASS**

## Ladder Controls (A1.1)
- k=2: COLLAPSES via rule 1 (exact-zero control)
- k=3: AMBIGUOUS via rule 2 (ill-conditioning refusal)
- synthetic (1e-8,1e-12,1e-12): COLLAPSES via rule 4
- edges at ratio 2.0 and 0.5: both COLLAPSES via rule 4
- (1e-8,2e-11,2e-11): AMBIGUOUS via rule 5

## Adjudicator (Q7)
- All 7 outcomes reachable: YES
- All 4 precedence collisions correct: YES

## G-ALIGN
- structural predicate: True
- arm (changed r[0]): predicate fails=True
- **PASS**

## G-WIRE
- ||L Q_0||_{M_h,F} = 5.46e-09
- arm: ||L_mut Q_0|| = 1.00e-04, red=True
- **PASS**

## G-RANK n=12
- dim=13 (expected 13)
- gap={'gap': 7218417941.362989, 'state': 'measured'}
- arm: n=11 gives dim=0
- **PASS**

## SVD Bridge n=12
- dim match: True
- basis bitwise match: True
- (diagnostic, routes nothing)

## G-SUBSPACE n=12
- k_svd=13, rank_avg=13
- sin(theta_C)=7.2888e-11 (gate <=1e-6)
- **PASS**

## G-SAMPLE n=12
- svd: kappa=3.7073e+01
- avg: kappa=3.7073e+01
- arm (zeroed column): rank=12, red=True

## G-BASIS n=12
- isometric rotation: covariance and projector invariant on both routes
- arm (non-isometric): ||P'-P||=3.00e+00, red=True

## Target n=12
- sin(theta_C)=7.2888e-11, sin(theta_Q)=2.1197e-10
- K_floor=6.8509e-07 (dominant: discrepancy)
- svd: ||K||=3.6075e+00, J=PERSISTS
- avg: ||K||=3.6075e+00, J=PERSISTS

## G-RANK n=20
- dim=21 (expected 21)
- gap={'gap': 2958163865.946767, 'state': 'measured'}
- arm: n=19 gives dim=0
- **PASS**

## SVD Bridge n=20
- dim match: True
- basis bitwise match: True
- (diagnostic, routes nothing)

## G-SUBSPACE n=20
- k_svd=21, rank_avg=21
- sin(theta_C)=8.0860e-11 (gate <=1e-6)
- **PASS**

## G-SAMPLE n=20
- svd: kappa=6.2781e+02
- avg: kappa=6.2781e+02
- arm (zeroed column): rank=20, red=True

## G-BASIS n=20
- isometric rotation: covariance and projector invariant on both routes
- arm (non-isometric): ||P'-P||=3.00e+00, red=True

## Target n=20
- sin(theta_C)=8.0860e-11, sin(theta_Q)=3.0711e-10
- K_floor=9.9259e-07 (dominant: discrepancy)
- svd: ||K||=6.4410e+01, J=PERSISTS
- avg: ||K||=6.4410e+01, J=PERSISTS

## Adversarial Vacuity Audit (Q8)
- All 8 gates examined
- Key finding: G-BASIS at n=0 has no power (contract addresses via G-WIRE)

## Contract Contradictions (Q9)
- No blocking contradictions found
- 8 items examined and resolved

## End-of-run verification
- All inputs unchanged: YES

## Output Manifest
- 34 files
- total time: 2240.5s
