# M9.52: \(P_{\mathrm{flat}}\) predicts \(\delta S\) at every \(\alpha\)

> Paper 37: modular \(\mathrm{Tr}(K\Delta C)\) and CHM.
> Paper 45: well-inside enclosed energy. Paper 60: finite
> \(\delta S\) on one enclosing ball is \(2h(\alpha)\).
> This run puts three predictors on 216 balls at three \(\alpha\).

## Equations

\[
P_{\mathrm{flat}}=\sum_B\delta e,\qquad
P_{\mathrm{CHM}}=\sum_B(R^2-r^2)\delta e,\qquad
T_K=\mathrm{Tr}(K_{\mathrm{vac}}\Delta C).
\]

Winner: largest \(\lvert\rho(\delta S,\cdot)\rvert\).

## Verdicts

\(N=12\), \(R=3\), \(216\) balls, \(120\) well-inside.

| \(\alpha\) | \(\rho_{\mathrm{flat}}\) | \(\rho_{\mathrm{CHM}}\) | \(\rho_K\) | winner |
| --- | --- | --- | --- | --- |
| \(0.005\) | \(0.991\) | \(0.963\) | \(0.984\) | \(P_{\mathrm{flat}}\) |
| \(0.02\) | \(0.996\) | \(0.948\) | \(0.974\) | \(P_{\mathrm{flat}}\) |
| \(0.08\) | \(0.999\) | \(0.930\) | \(0.961\) | \(P_{\mathrm{flat}}\) |

Well-inside: same winner at every \(\alpha\). C_track
**PASS**. C_flip **FAIL** (stable). CHM and \(T_K\)
get *worse* as \(\alpha\) grows. Auditor \(N=10\):
\(P_{\mathrm{flat}}\) at \(0.01\) and \(0.06\). Flip
**REFUTED** (stability confirmed).

`WINNER_STABLE`. *computed.* Finite \(\delta S\) on this
cloud tracks enclosed energy, not the modular kernel.
That matches Paper 61: the reusable object is
\(\sum\delta e\). Paper 37's CHM win is a statement
about \(K\), not about finite \(\delta S(\alpha)\).

## Equation-to-code

| Object | Where |
| --- | --- |
| 216 balls, three predictors | `scripts/m9_52_predictor.py` |
| Adversary | `scripts/m9_52_audit_pred.py` |

Paper: [`../latex/62_Predictor.tex`](../latex/62_Predictor.tex).
