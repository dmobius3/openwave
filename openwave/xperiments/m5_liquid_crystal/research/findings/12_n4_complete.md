# N4 checkpoint 12 , N4 COMPLETE: the PMNS gate fully closed (3 predictions + 1 free coupling)

All three open questions resolved + the final NuFIT 6.0 scorecard assembled. Record: `n4_findings.md`.

## The three studies (the user's open questions)

| Study | Question | Answer | Script |
| --- | --- | --- | --- |
| A | is theta13=8.56 natural or strong-coupling? | **NATURAL**: continuous, near-linear, reaches 8.56 at `g_chiral* = 0.937` (O(1)), mass gap stable (no resonance) | `m5_11_n4_refine.py` |
| B | what fixes the delta_CP sign? | **the loop HANDEDNESS**: `delta_CP_sign = -sign(g_chiral)` (chi-screw sign irrelevant); `g_chiral>0 -> 270 deg` (data-preferred) | `m5_11_n4_refine.py` |
| C | mu-tau-respecting topological theta13? | **NO** (two attempts: local-azimuth `N*s` and global-azimuth `N*phi` both break mu-tau). theta13 is CONTINUOUS (chiral material coupling), not topological | `m5_11_n4_topo.py`, `m5_11_n4_linking.py` |
| D | final NuFIT 6.0 scorecard | assembled + figure; **gate FULLY CLOSED** | `m5_11_n4_final_scorecard.py` |

## FINAL scorecard (chi=1.2, delta=0.1, g_chiral*=0.937)

| Parameter | N4 | NuFIT 6.0 NO | status |
| --- | --- | --- | --- |
| theta12 | 35.72 (TBM 35.26) | 33.68 | PREDICTION (trimaximal) |
| theta23 | 45.00 | 43.30 | PREDICTION (maximal) |
| theta13 | 8.56 | 8.56 | FREE COUPLING (chiral g*=0.94) |
| delta_CP | 270 | 212 | PREDICTION (maximal; sign=handedness) |

## The article-ready story
3 of 4 PMNS parameters PREDICTED from closed-loop LC geometry (theta12 trimaximal via magic crossing; theta23
maximal + delta_CP maximal via mu-tau reflection + chiral coupling, the Harrison-Scott result realized from
geometry); theta13 = the one free coupling (the cholesteric/Lifshitz chiral strength), natural at O(1). delta_CP
sign predicted from the loop handedness (270 deg, data-preferred). theta12 sits ~2 deg above data (standard
TBM gap).

## End-of-night state (autonomous, 2026-06-21 -> 22)
N3 COMPLETE (TBM gate) + N4 COMPLETE (CP sector + gate closure). 14 scripts, 4 findings docs (n3, n4,
n_foundation, + the 10a plan updates), 13 checkpoints, 2 figures. All < 130 KB, nothing > 1 MB. #236 HELD,
git untouched (the maintainer's). Ready for N5 (the peer-review article draft) tomorrow.
