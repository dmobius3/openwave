# M9.17--M9.18: first law, bulk vs horizon

> Paper 26 mixed every hop and a flat kernel won. The cut is
> where Jacobson puts Clausius. Split the ball.

## Equations

Same as M9.16. Fit \(K_{\mathrm{CHM}}=aw+b\) and
\(K_{\mathrm{flat}}=\langle K_{\mathrm{NN}}\rangle\) on the
scored set only. Compare \(\delta S\) to \(\mathrm{Tr}(K\Delta C)\).
Bulk: \(r_{\mathrm{mid}}\le 0.5R\). Surface / horizon:
\(r_{\mathrm{mid}}\ge 0.75R\).

C3 primary: \(R_{\mathrm{shape}}(\mathrm{CHM})<R_{\mathrm{shape}}(\mathrm{flat})\).
Tracking floors C0/C1/C2 are the M9.16 numbers, sign-blind on
M9.18.

## M9.17 bulk (locked, \(N=12\), \(R=4\))

| Set | \(n\) | \(\rho_{\mathrm{CHM}}\) | \(R_{\mathrm{CHM}}\) | \(R_{\mathrm{flat}}\) | C3 |
| --- | --- | --- | --- | --- | --- |
| 1d central | 13 | \(-0.951\) | \(0.309\) | \(0.314\) | CHM slightly |
| 3d bulk | 84 | \(-0.428\) | \(0.904\) | \(0.868\) | flat wins |
| 3d surface *diagnostic* | 294 | \(-0.589\) | \(0.808\) | \(0.994\) | CHM, not a pass |

C0 used \(\rho>0.70\) not \(\lvert\rho\rvert\), so the instrument
bit rejected a \(\rho=-0.951\) track. Bulk C3 **FAIL**. Surface
numbers on this sample are **not** a pre-registered pass.

An independent bulk adversary (\(N=11\) and even \(N=8\),
\(R=3\)) **REFUTES** 3d bulk tracking: \(\lvert\rho\rvert\le 0.18\),
\(R_{\mathrm{shape}}\ge 0.98\). C3 is vacuous when both residuals
are noise, and **FAIL** on \(N=8\). Occupancy flickers on odd
cubes; an even \(N=10\) identity control still has
\(\rho_{\mathrm{CHM}}=-0.15\). 1d C0 sign is size-fragile
(\(+0.925\) vs locked \(-0.953\)). This confirms the bulk
negative. It does not audit Paper 27's surface C3.

## M9.18 surface on a new grid (locked)

| Set | \(n\) | \(\rho_{\mathrm{CHM}}\) | \(R_{\mathrm{CHM}}\) | \(R_{\mathrm{flat}}\) | C3 |
| --- | --- | --- | --- | --- | --- |
| 1d ends \(L=24\) | 12 | \(-0.591\) | \(0.807\) | \(0.956\) | CHM |
| 3d \(N=10\), \(R=3\) | 150 | \(-0.393\) | \(0.919\) | \(0.996\) | **PASS** |
| auditor 1d \(L=20\) | 10 | \(+0.646\) | \(0.763\) | \(0.417\) | flat |
| auditor 3d \(N=12\), \(R=3\) | 150 | \(-0.575\) | \(0.818\) | \(0.969\) | **CONFIRMED** |

C0/C1/C2 fail (tracking below \(0.70/0.60/0.70\)).
C3 on the 3d horizon **PASS**, auditor **CONFIRMED**.

## Verdict

`SURFACE_CHM_NOT_SELECTED` as an all-gate pass (tracking floors
fail). C3 only: on the cut, the CHM envelope **beats** a flat
local kernel. In the bulk it does not. That is the horizon
pattern Jacobson needs, at weak correlation.

*computed.* Not \(\eta=1/4G\). Not Einstein. Not de Sitter.

## Equation-to-code

| Object | Where |
| --- | --- |
| Bulk split | `scripts/m9_17_bulk_first_law.py` |
| Surface lock | `scripts/m9_18_surface_first_law.py` |
| Surface audit | `scripts/m9_18_audit_surface.py` |

Paper: [`../latex/27_Horizon_Selects_CHM.tex`](../latex/27_Horizon_Selects_CHM.tex).
