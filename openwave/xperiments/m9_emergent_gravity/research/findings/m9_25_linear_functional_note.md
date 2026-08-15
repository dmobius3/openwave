# M9.25: \(\delta S\) of balls as a linear functional of local energy

> Closest test of ``entanglement \(\Rightarrow\) linearized gravity''.
> Many balls, one weak Gaussian source, fixed hop graph.

## Setup

\(N=12\). \(R=2\) balls at every legal center (\(512\) balls).
\(V=\varepsilon\exp(-r^2/2\sigma^2)\), \(\sigma=2\), \(\varepsilon=0.05\)
and \(0.10\). \(e_i=\sum_j H_{ij}C_{ij}\). Predictors
\(P=\sum_{\mathrm{ball}}w\,\delta e\) with \(w=R^2-r^2\), \(R-r\), \(1\).

## Verdicts

| Gate | Result |
| --- | --- |
| C0 \(\max\|\delta S\|/\langle S\rangle>10^{-4}\) | FAIL \(3.7\times 10^{-5}\) (lock too tight vs area) |
| C1 Pearson\(\delta S(\varepsilon),\delta S(2\varepsilon)\) | PASS \(0.99999999\) |
| C2 CHM beats flat | PASS \(0.168<0.186\) |
| C3 CHM beats linear | FAIL \(0.168>0.165\) |
| C4 \(\lvert\rho(\delta S,P_{\mathrm{CHM}})\rvert>0.60\) | PASS \(-0.986\) |

Auditor \(N=10\), \(216\) balls: \(\rho_{\mathrm{CHM}}=-0.975\),
C4 PASS. C2 **REFUTED** (flat \(0.197\) beats CHM \(0.223\)).

`NOT_CHM_LINEAR_FUNCTIONAL` as an all-gate pass.

What *is* confirmed: \(\delta S\) of balls is a **linear**
functional of local energy (\(\rho_{\varepsilon}=1\),
\(\lvert\rho_{\mathrm{CHM}}\rvert=0.986\)). What is **not**:
the kernel is uniquely CHM / Einstein. A smooth source makes
flat, linear, and CHM almost the same predictor. Auditor: flat
wins.

This is the first-law half of FGHMV as a lattice fact about a
free fermion. It is not \(G_{\mu\nu}=8\pi G T_{\mu\nu}\).

Scripts: `m9_25_linear_functional.py`, `m9_25_audit_linear.py`.
Paper: [`../latex/34_Linear_Functional.tex`](../latex/34_Linear_Functional.tex).
