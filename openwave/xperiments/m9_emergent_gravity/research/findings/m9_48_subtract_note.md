# M9.48: the first law subtracts the sea

> If \(\rho=e_i(C_0)<0\) went into Gauss, \(a\) would be
> outward. The first law is \(\delta S=\kappa\sum_B\delta e\),
> not \(\kappa\sum_B e\).

## Equations

\[
e_i(C)=\sum_j H_{ij}C_{ij},\qquad
P_{\mathrm{vac}}=\sum_B e(C_0),\qquad
P_{\delta e}=\sum_B\bigl(e(C_1)-e(C_0)\bigr).
\]

\(P_e=\sum_B e(C_1)=P_{\mathrm{vac}}+P_{\delta e}\).
\(M_{\mathrm{FL}}=\delta S/\kappa\). Diagnostic, not a
gate: \(a_{\mathrm{vac}}=-P_{\mathrm{vac}}/R^2\).

## Verdicts

\(N=12\), compact packet. \(P_{\mathrm{vac}}<0\),
\(P_{\delta e}>0\) at every \(R\). At \(R=5\),
\(\lvert P_{\mathrm{vac}}/P_{\delta e}\rvert=2762\).
Raw energy *is* the sea.

\(\rho(\delta S,P_{\delta e})=0.99999\).
\(\rho(\delta S,P_e)=-0.641\). The first law tracks
the packet, not the raw sum. \(M_{\mathrm{FL}}\)
matches \(P_{\delta e}\) on enclosing balls.

\(a_{\mathrm{vac}}\) is outward. That is the diagnostic
of feeding \(E_{\mathrm{vac}}\) to Gauss. It is not a
first-law force.

Auditor \(N=10\): \(\rho_{\delta e}=0.9999\),
\(\rho_e=-0.817\), scale \(973\). C_sign and C_sub
**CONFIRMED**.

`SEA_SUBTRACTED`. *computed.* The Fermi sea is not
repulsive Newton gravity. The instrument subtracts it.
What gravitates is the excess.

## Equation-to-code

| Object | Where |
| --- | --- |
| Split, correlations, diagnostic | `scripts/m9_48_subtract.py` |
| Adversary | `scripts/m9_48_audit_subtract.py` |

Paper: [`../latex/58_Sea_Subtracted.tex`](../latex/58_Sea_Subtracted.tex).
