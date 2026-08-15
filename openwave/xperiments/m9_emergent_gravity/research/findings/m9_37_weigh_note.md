# M9.37: \(\kappa\) weighs the enclosed mass and finds the site

> Paper 46 gave \(\kappa\). This note uses it.

## Equations

Split well-inside balls (seed 37). Even half:

\[
\kappa=\mathrm{median}(\delta S/P_{\mathrm{flat}}).
\]

Odd half: \(M_{\mathrm{hat}}=\mathrm{median}(\delta S/\kappa)\),
compared to \(M_{\mathrm{enc}}=\mathrm{median}(P_{\mathrm{flat}})\).
Location: intersection of balls with \(\delta S>0.70\max\delta S\);
among survivors, the site in the most high-\(\delta S\) balls.

## Verdicts

\(N=12\), \(R=3\), source \((6,6,6)\). \(120\) well-inside,
\(60/60\) split. \(\kappa=0.968\).

| Gate | Result |
| --- | --- |
| C_vac | PASS \(0.974\) |
| C_mass | PASS, \(\lvert M_{\mathrm{hat}}/M_{\mathrm{enc}}-1\rvert=2\times 10^{-8}\) |
| C_loc | PASS, \(\hat x=(6,6,6)\) |
| C_pred | PASS \(0.999\) |

Global \(\sum\delta e\) is \(28\%\) larger than \(M_{\mathrm{enc}}\):
Gauss does not see the tail.

Auditor \(N=10\), source \((4,5,5)\): C_mass **CONFIRMED**
(\(1.9\%\)). C_loc **REFUTED** (\(\hat x=(4,4,4)\), Chebyshev
\(1\)).

`WEIGHED_AND_LOCATED` on the primary lattice. Location is
exact here, off by one site on the auditor. Not \(1/4G\).
Not Einstein. \(\kappa\) is an instrument that reads enclosed
energy from \(\delta S\) and points at the source.

## Equation-to-code

| Object | Where |
| --- | --- |
| Split, \(\kappa\), locate, gates | `scripts/m9_37_weigh.py` |
| Adversary | `scripts/m9_37_audit_weigh.py` |

Paper: [`../latex/47_Weigh.tex`](../latex/47_Weigh.tex).
