# M9.51: Gauss from \(\sum\delta e\). No \(\kappa\). Shapes hold.

> Paper 60: \(\kappa(\alpha)\) is not a coupling. The mass
> is \(P_{\mathrm{flat}}=\sum_B\delta e\). Rebuild
> \(a(R)=-P_{\mathrm{flat}}(R)/R^2\) and scan \(\alpha\).

## Equations

\[
a(R)=-\frac{1}{R^2}\sum_{i\in B(R)}\delta e_i.
\]

No \(\delta S\). No \(\kappa\). STAR: open compact packet,
slope on enclosing \(R\). SEA: periodic band edges, all \(R\).

## Verdicts

\(N=12\), \(\alpha\in\{0.01,0.02,0.04\}\). \(P_{\mathrm{flat}}\propto\alpha\),
so \(\lvert a\rvert\) scales and the log-log slope cannot
run. Measured spread \(10^{-13}\).

| State | slope (every \(\alpha\)) |
| --- | --- |
| star \(R=3,4,5\) | \(-1.9979\) |
| sea \(R=2,3,4,5\) | \(+0.9673\) |

C_star, C_sea, C_hold **PASS**. Auditor \(N=10\):
star \(-1.949\), sea \(+0.980\), all **CONFIRMED**.

Paper 55's sea slope \(+1.266\) used \(\delta S/\kappa\)
with a single-ball \(\kappa\). That mixed Paper 60's
running \(\kappa\) into the shape. Energy-only Gauss
is closer to \(+1\).

`ENERGY_GAUSS_HOLDS`. *computed.* After \(\kappa\) falls,
the two force laws remain, and they do not depend on
the size of the kick. Still inherited Gauss.

## Equation-to-code

| Object | Where |
| --- | --- |
| Star, sea, \(\alpha\) scan | `scripts/m9_51_energy_gauss.py` |
| Adversary | `scripts/m9_51_audit_energy.py` |

Paper: [`../latex/61_Energy_Gauss.tex`](../latex/61_Energy_Gauss.tex).
