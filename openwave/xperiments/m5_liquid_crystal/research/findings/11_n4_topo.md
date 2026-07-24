# N4 checkpoint 11 , topological self-linking test: INCONCLUSIVE (honest negative) + the insight

`m5_11_n4_topo.py` (`m5_11_n4_topo_summary.json`). Tested the "article-grade" idea: is the loop handedness a TOPOLOGICAL
integer (self-linking N, `Lk = Tw + Wr`) that quantizes theta13 + delta_CP?

## Result , NOT a clean quantization

Adding N full secondary-frame twists around the loop azimuth (`N*s`), scanned N = -2..2:

| N | theta12 | theta23 | theta13 | delta_CP |
| --- | --- | --- | --- | --- |
| -2 | 32.96 | 44.54 | 0.37 | -18.0 |
| -1 | 32.83 | 44.70 | 0.27 | -18.3 |
| **0** | **35.26** | **45.00** | 0.14 | **-90 (maximal)** |
| 1 | 32.79 | 45.36 | 0.37 | -13.0 |
| 2 | 33.33 | 46.20 | 1.26 | -3.9 |

N != 0 BREAKS mu-tau and degrades TBM (theta12 -> ~33, theta23 wobbles) and delta_CP goes intermediate +
non-antisymmetric. Only N=0 keeps the clean structure. So "theta13 topologically quantized" is NOT supported
by this construction.

## The insight (a real physics lesson, worth keeping)

A mu-tau-SYMMETRIC arrangement forces mu and tau to carry OPPOSITE self-linking (the mu-tau mirror sends
s -> -s, so N -> -N) -> the net handedness tends to CANCEL. So topological self-linking and a mu-tau-symmetric
NONZERO delta_CP are in TENSION. That is exactly WHY the clean maximal delta_CP comes from the CONTINUOUS
chiral (Lifshitz material) coupling (n4_chiral, N=0), not a topological integer.

## Status of N4
- SOLID: delta_CP = +-90 (maximal CP violation), robust, predicted (n4_chiral, N=0). theta13 <-> delta_CP
  unified by the chiral coupling. = the mu-tau-reflection result realized from loop geometry.
- OPEN: a mu-tau-RESPECTING topological framing (resolve the tension above) , next session / Duda.

## End-of-night state (autonomous run, 2026-06-21)
N3 COMPLETE (TBM gate + theta13 crux + scorecard). N4 first pass: delta_CP=maximal SOLID; topological-theta13
INCONCLUSIVE (honest). Clean stopping point. #236 HELD, git untouched (the maintainer's). No raw data > 1 MB.
