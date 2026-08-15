# M9.50: \(\kappa\) runs with \(\alpha\). It is \(2h(\alpha)/(\alpha\Delta E)\)

> Papers 46--59 used \(\alpha=0.02\). At that value
> \(S_{\mathrm{global}}=2h(\alpha)\approx\delta S(B_{\mathrm{enc}})\).
> \(h(\alpha)/alpha\) runs as \(-\log\alpha\). A reusable
> first-law constant cannot.

## Equations

\[
h(\alpha)=-\alpha\log\alpha-(1-\alpha)\log(1-\alpha),\qquad
S_{\mathrm{global}}=2h(\alpha).
\]

Enclosing ball \(R=3\):

\[
\kappa(\alpha)=\delta S/P_{\mathrm{flat}},\qquad
r_{sg}=S_{\mathrm{global}}/P_{\mathrm{flat}}.
\]

C_lin: rel range of \(\kappa(\alpha)<0.10\).
C_sg: rel range of \(r_{sg}>0.20\).

## Verdicts

\(N=12\), \(\alpha\in\{0.005,0.01,0.02,0.04,0.08\}\).
\(C\) stays in \([0,1]\).

| \(\alpha\) | \(\kappa\) | \(r_{sg}\) |
| --- | --- | --- |
| \(0.005\) | \(1.343\) | \(1.348\) |
| \(0.01\) | \(1.196\) | \(1.199\) |
| \(0.02\) | \(1.047\) | \(1.050\) |
| \(0.04\) | \(0.897\) | \(0.899\) |
| \(0.08\) | \(0.745\) | \(0.746\) |

Rel range of \(\kappa\): \(0.57\). C_lin **FAIL**.
Rel range of \(r_{sg}\): \(0.57\). C_sg **PASS**.
\(\rho(\kappa,2h(\alpha)/\alpha)=0.999998\).

Auditor \(N=10\): \(\kappa=1.26,1.01,0.85\), rel
\(0.40\). C_lin **REFUTED**. C_sg **CONFIRMED**.

`KAPPA_RUNS_WITH_ALPHA`. *computed.* Paper 46's
``universal \(\kappa\)'' is at fixed \(\alpha\). Across
\(\alpha\), \(\kappa=2h(\alpha)/(\alpha\Delta E)\).
The reusable mass is \(\sum\delta e\), not \(\delta S/\kappa\).
Force-law papers at fixed \(\alpha=0.02\) still used
that same \(M=P_{\mathrm{flat}}\). They stand.
\(\kappa\) is not \(1/T\) independent of the kick.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(\alpha\) scan | `scripts/m9_50_alpha.py` |
| Adversary | `scripts/m9_50_audit_alpha.py` |

Paper: [`../latex/60_Kappa_Runs.tex`](../latex/60_Kappa_Runs.tex).
