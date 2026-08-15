# M9.35: both-inside with \(n=69\) --- Gauss wins, CHM loses

> Paper 44's Gauss-versus-boost subset had \(n=11\).
> \(R=3\) gives \(69\) balls that contain both masses.

## Equations

Same two-source \(C\) as M9.34. Balls of radius \(3\).
\(P_{\mathrm{flat}}=\sum_B\delta e\) (enclosed energy).
\(P_{\mathrm{CHM}}=\sum_B(R^2-r^2)\,\delta e\) (boost).

## Verdicts

Instrument holds (\(\rho(K_{\mathrm{vac}})=0.983\)).

| Set | \(n\) | \(\rho_{\mathrm{CHM}}\) | \(\rho_{\mathrm{flat}}\) | \(R_{\mathrm{CHM}}\) | \(R_{\mathrm{flat}}\) |
| --- | --- | --- | --- | --- | --- |
| all | \(216\) | \(0.964\) | \(0.996\) | \(0.266\) | \(0.085\) |
| both-inside | \(69\) | \(0.955\) | \(0.999\) | \(0.298\) | \(0.053\) |

On both-inside, \(\mathrm{std}(P_{\mathrm{CHM}})=0.458\) is
ten times \(\mathrm{std}(P_{\mathrm{flat}})=0.047\). \(\delta S\)
follows the *small* enclosed-energy variation, not the large
boost-weight variation.

C2b PRIMARY **FAIL**. Auditor \(n_{\mathrm{both}}=47\):
C2b **REFUTED** (\(R=0.318>0.053\), \(\rho_{\mathrm{flat}}=0.999\)).

`BOTH_INSIDE_GAUSS`. *computed.* Paper 37's CHM win was
\(R=2\), where ``is it in the ball'' and ``where in the ball''
are mixed. Once both masses sit well inside, the first law
is enclosed energy. That is Gauss, not the CHM boost, not
\(8\pi G\), not Einstein.

## Equation-to-code

| Object | Where |
| --- | --- |
| Pair, \(R=3\), gates | `scripts/m9_35_both_inside.py` |
| Adversary | `scripts/m9_35_audit_both.py` |

Paper: [`../latex/45_Both_Inside_Gauss.tex`](../latex/45_Both_Inside_Gauss.tex).
